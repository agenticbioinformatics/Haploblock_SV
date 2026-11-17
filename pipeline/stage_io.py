"""Shared input validation, stage-config handling and CLI error reporting.

Every stage after Stage 1 reads the ``config.yaml`` written by an upstream
stage, adds its own outputs under ``paths`` and writes the extended config to
its own output directory. These helpers keep that contract in one place.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Callable, Iterable, Optional, Sequence

import pandas as pd
import yaml


log = logging.getLogger("haploblock_sv")


def resolve_path(path: str | Path, base_dir: Path) -> Path:
    path = Path(path)
    return path if path.is_absolute() else base_dir / path


def _resolve_registered_paths(paths: dict, base_dir: Path) -> dict:
    resolved = {}
    for name, value in paths.items():
        if isinstance(value, dict):
            resolved[name] = {key: str(resolve_path(item, base_dir)) for key, item in value.items()}
        elif isinstance(value, str) and not value.startswith(("http://", "https://")):
            resolved[name] = str(resolve_path(value, base_dir))
        else:
            resolved[name] = value
    return resolved


def load_config(path: Path, required_paths: Iterable[str] = ()) -> dict:
    """Load a stage config, resolving relative paths against its directory."""
    if not path.is_file():
        raise FileNotFoundError(f"Stage config not found: {path}")
    config = yaml.safe_load(path.read_text()) or {}
    if not isinstance(config.get("paths"), dict):
        raise ValueError(f"{path} has no 'paths' section; pass a config.yaml written by a pipeline stage")
    missing = [name for name in required_paths if name not in config["paths"]]
    if missing:
        raise ValueError(
            f"{path} does not register {missing}; run the upstream stage that produces "
            "them and pass its config.yaml"
        )
    config["paths"] = _resolve_registered_paths(config["paths"], path.parent)
    return config


def read_tsv(path: str | Path, required_columns: Sequence[str] = (), **kwargs) -> pd.DataFrame:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Input table not found: {path}")
    table = pd.read_csv(path, sep="\t", **kwargs)
    missing = sorted(set(required_columns) - set(table.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")
    return table


def read_sample_metadata(config: dict) -> pd.DataFrame:
    metadata = read_tsv(config["paths"]["sample_metadata"], ["sample_id", "population"])
    if metadata.empty:
        raise ValueError(
            f"{config['paths']['sample_metadata']} has no samples; check that the sample "
            "metadata IDs match the VCF sample names"
        )
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    return metadata


def write_stage_config(
    config: dict,
    out_dir: Path,
    paths: Optional[dict] = None,
    thresholds: Optional[dict] = None,
    settings: Optional[dict] = None,
) -> Path:
    """Carry an upstream config forward with this stage's additions."""
    output = dict(config)
    for section, updates in (("paths", paths), ("thresholds", thresholds), ("settings", settings)):
        if updates:
            output[section] = {**config.get(section, {}), **updates}
    config_path = out_dir / "config.yaml"
    config_path.write_text(yaml.safe_dump(output, sort_keys=False))
    return config_path


def run_cli(main: Callable[[Optional[list]], object]) -> None:
    """Run a stage entry point, turning expected input errors into a clear message."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        result = main(sys.argv[1:])
    except (FileNotFoundError, ValueError, KeyError, RuntimeError) as error:
        log.error("%s", error)
        sys.exit(1)
    sys.exit(result if isinstance(result, int) else 0)
