"""Tests for shared stage-config handling."""

from pathlib import Path

import pytest
import yaml

from stage_io import load_config


def test_load_config_resolves_local_paths_and_keeps_urls(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump({"paths": {
        "gtf": "https://example.org/genes.gtf",
        "samples": "samples.tsv",
        "sv_genotypes": {"chr1": "sv_genotypes.chr1.tsv", "chr2": "https://example.org/chr2.tsv"},
    }}))

    paths = load_config(config_path)["paths"]

    assert paths["gtf"] == "https://example.org/genes.gtf"
    assert paths["samples"] == str(tmp_path / "samples.tsv")
    assert paths["sv_genotypes"]["chr1"] == str(tmp_path / "sv_genotypes.chr1.tsv")
    assert paths["sv_genotypes"]["chr2"] == "https://example.org/chr2.tsv"


def test_load_config_rejects_null_paths_with_clear_message(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump({"paths": {"sv_genotypes": {"chr1": None}}}))

    with pytest.raises(ValueError, match=r"paths\.sv_genotypes\.chr1 is empty"):
        load_config(config_path)
