"""Haploblock cluster files from data.haploblocks.org: download, cache and parse.

Each haploblock has one file ``<chrom>_<start>-<end>_cluster.tsv`` with two
columns: cluster ID and haplotype ID ``<sample>_<chrom>_region_<start>-<end>_hap<0|1>``.
"""

from __future__ import annotations

import hashlib
import logging
import re
import subprocess
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from sv_table import canonical_sample_id


log = logging.getLogger("haploblock_clusters")

_USE_SYSTEM_CURL = False
_DOWNLOAD_CLIENT_LOCK = threading.Lock()
_DOWNLOAD_SESSION = threading.local()
DEFAULT_BASE_URL = "https://data.haploblocks.org/haploblock_hashes/1000G"
CLUSTER_FILENAME_RE = re.compile(r"^(chr[^_]+)_(\d+)-(\d+)_cluster\.tsv$")
HAPLOTYPE_RE = re.compile(
    r"^(?P<sample>.+)_(?P<chrom>chr[^_]+)_region_(?P<start>\d+)-(?P<end>\d+)_hap(?P<haplotype>[01])$"
)


@dataclass
class BlockMembership:
    haploblock_id: str
    chrom: str
    start: int
    end: int
    sample_clusters: dict[str, tuple[str, str]]
    total_cluster_haplotypes: Counter[str]


def download_with_system_curl(url: str, timeout: int) -> bytes:
    """Download with the OS HTTPS client, which trusts the macOS keychain."""
    result = subprocess.run(
        ["curl", "--fail", "--silent", "--show-error", "--location", "--max-time", str(timeout), url],
        check=True,
        capture_output=True,
    )
    return result.stdout


def request_with_retries(url: str, retries: int, timeout: int = 60) -> bytes:
    import requests

    global _USE_SYSTEM_CURL
    error: Exception | None = None
    for attempt in range(retries):
        try:
            with _DOWNLOAD_CLIENT_LOCK:
                use_system_curl = _USE_SYSTEM_CURL
            if use_system_curl:
                return download_with_system_curl(url, timeout)
            session = getattr(_DOWNLOAD_SESSION, "session", None)
            if session is None:
                session = requests.Session()
                _DOWNLOAD_SESSION.session = session
            response = session.get(url, timeout=timeout)
            response.raise_for_status()
            return response.content
        except requests.RequestException:
            with _DOWNLOAD_CLIENT_LOCK:
                if not _USE_SYSTEM_CURL:
                    log.warning(
                        "Python's HTTPS client could not reach the server; using the system HTTPS client instead"
                    )
                    _USE_SYSTEM_CURL = True
            try:
                return download_with_system_curl(url, timeout)
            except subprocess.CalledProcessError as exc:
                error = exc
        except subprocess.CalledProcessError as exc:
            error = exc
        if attempt + 1 < retries:
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Failed to download {url} after {retries} attempts: {error}")


def discover_cluster_filenames(base_url: str, chrom: str, retries: int) -> list[str]:
    directory_url = f"{base_url.rstrip('/')}/{chrom}/clusters/"
    log.info("Discovering cluster files for %s", chrom)
    html = request_with_retries(directory_url, retries).decode()
    pattern = re.compile(rf"{re.escape(chrom)}_\d+-\d+_cluster\.tsv")
    names = sorted(set(pattern.findall(html)), key=parse_cluster_filename)
    if not names:
        raise RuntimeError(f"No cluster files found at {directory_url}")
    return names


def parse_cluster_filename(path_or_name: str | Path) -> tuple[str, int, int]:
    name = Path(path_or_name).name
    match = CLUSTER_FILENAME_RE.fullmatch(name)
    if not match:
        raise ValueError(f"Unexpected cluster filename: {name}")
    return match.group(1), int(match.group(2)), int(match.group(3))


def cluster_cache_directory(out_dir: Path, base_url: str, chrom: str) -> Path:
    normalized_url = base_url.rstrip("/")
    source_id = hashlib.sha256(normalized_url.encode()).hexdigest()[:16]
    return out_dir / "_intermediate" / "clusters" / source_id / chrom


def prepare_cluster_files(
    chrom: str,
    base_url: str,
    cache_dir: Path,
    cluster_dir: Path | None,
    workers: int,
    retries: int,
) -> tuple[list[Path], dict]:
    if cluster_dir is not None:
        paths = sorted(cluster_dir.glob(f"{chrom}_*-*_cluster.tsv"), key=parse_cluster_filename)
        if not paths:
            raise ValueError(f"No {chrom} cluster files found in {cluster_dir}")
        return paths, {
            "source": "local",
            "directory": str(cluster_dir.resolve()),
            "discovered": len(paths),
            "downloaded": 0,
            "reused": len(paths),
        }

    names = discover_cluster_filenames(base_url, chrom, retries)
    cache_dir.mkdir(parents=True, exist_ok=True)
    reused = [cache_dir / name for name in names if (cache_dir / name).exists()]
    pending = [name for name in names if not (cache_dir / name).exists()]
    failures: list[str] = []
    log.info("Preparing %d cluster files (%d cached, %d to download)", len(names), len(reused), len(pending))

    def download(name: str) -> Path:
        url = f"{base_url.rstrip('/')}/{chrom}/clusters/{name}"
        destination = cache_dir / name
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        temporary.write_bytes(request_with_retries(url, retries))
        temporary.replace(destination)
        return destination

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(download, name): name for name in pending}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as exc:  # individual failures are collected and reported together
                name = futures[future]
                log.error("Failed to download %s: %s", name, exc)
                failures.append(name)
    if failures:
        raise RuntimeError(
            f"Failed to download {len(failures)} cluster file(s); rerun to resume: {sorted(failures)}"
        )
    paths = sorted((cache_dir / name for name in names), key=parse_cluster_filename)
    return paths, {
        "source": "download",
        "base_url": base_url.rstrip("/"),
        "cache_dir": str(cache_dir.resolve()),
        "discovered": len(names),
        "downloaded": len(pending),
        "reused": len(reused),
    }


def parse_haplotype_id(value: str) -> tuple[str, str, int, int, int]:
    match = HAPLOTYPE_RE.fullmatch(value)
    if not match:
        raise ValueError(f"Unexpected haplotype identifier: {value}")
    return (
        canonical_sample_id(match.group("sample")),
        match.group("chrom"),
        int(match.group("start")),
        int(match.group("end")),
        int(match.group("haplotype")),
    )


def load_block_membership(
    path: Path,
    vcf_samples: set[str],
) -> tuple[BlockMembership, int, int]:
    chrom, start, end = parse_cluster_filename(path)
    haploblock_id = f"{chrom}_{start}_{end}"
    total_counts: Counter[str] = Counter()
    by_sample: dict[str, list[str | None]] = {}
    total_rows = 0
    used_rows = 0

    with path.open() as handle:
        for line in handle:
            cluster_id, member_id = line.rstrip("\n").split("\t")
            sample, member_chrom, member_start, member_end, haplotype = parse_haplotype_id(member_id)
            if (member_chrom, member_start, member_end) != (chrom, start, end):
                raise ValueError(f"Membership coordinates do not match filename in {path}: {member_id}")
            total_rows += 1
            total_counts[cluster_id] += 1
            if sample not in vcf_samples:
                continue
            used_rows += 1
            clusters = by_sample.setdefault(sample, [None, None])
            clusters[haplotype] = cluster_id

    complete_samples = {
        sample: (clusters[0], clusters[1])
        for sample, clusters in by_sample.items()
        if clusters[0] is not None and clusters[1] is not None
    }
    return (
        BlockMembership(
            haploblock_id=haploblock_id,
            chrom=chrom,
            start=start,
            end=end,
            sample_clusters=complete_samples,
            total_cluster_haplotypes=total_counts,
        ),
        total_rows,
        used_rows,
    )
