"""Genotype handling shared by the analysis stages.

Converts raw VCF ``GT`` strings to allele counts and alternate-allele dosage,
and provides the multiple-testing correction used by Stages 5 and 6.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from sv_contract import METADATA_COLUMNS


COMMON_GT_DOSAGES = {
    "0|0": 0.0, "0/0": 0.0, "0": 0.0,
    "0|1": 1.0, "1|0": 1.0, "0/1": 1.0, "1/0": 1.0, "1": 1.0,
    "1|1": 2.0, "1/1": 2.0,
    ".": np.nan, ".|.": np.nan, "./.": np.nan,
}
MISSING_GENOTYPES = (".", ".|.", "./.")


def genotype_counts(genotype: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    """Return alternate and called allele counts for one sample across all SVs.

    Any genotype containing ``.`` counts as missing. Multi-digit allele indices
    (e.g. ``12/0``) count as one allele each.
    """
    if pd.api.types.is_numeric_dtype(genotype):
        dosage = pd.to_numeric(genotype, errors="coerce")
        has_call = dosage.notna()
        return dosage.where(has_call, 0).astype(int).to_numpy(), (has_call.astype(int) * 2).to_numpy()

    genotype = genotype.astype("string").str.strip()
    missing = genotype.str.contains(".", regex=False, na=True)
    alternate = genotype.str.count(r"(?:^|[|/])[1-9][0-9]*(?=$|[|/])")
    called = genotype.str.count(r"(?:^|[|/])[0-9]+(?=$|[|/])")
    return (
        alternate.mask(missing, 0).fillna(0).astype(int).to_numpy(),
        called.mask(missing, 0).fillna(0).astype(int).to_numpy(),
    )


def dosage_matrix_from_genotypes(sv: pd.DataFrame, samples: list[str]) -> np.ndarray:
    """Return a variant-by-sample dosage matrix; missing calls are NaN."""
    matrix = np.empty((len(sv), len(samples)), dtype=float)
    for sample_index, sample in enumerate(samples):
        genotype = sv[sample]
        if pd.api.types.is_numeric_dtype(genotype):
            matrix[:, sample_index] = pd.to_numeric(genotype, errors="coerce")
            continue
        mapped = genotype.map(COMMON_GT_DOSAGES)
        unknown = mapped.isna() & genotype.notna() & ~genotype.isin(MISSING_GENOTYPES)
        if unknown.any():
            # Slow path for multi-allelic or haploid encodings outside the lookup table.
            alternate, called = genotype_counts(genotype)
            matrix[:, sample_index] = np.where(called > 0, alternate, np.nan)
        else:
            matrix[:, sample_index] = mapped.to_numpy(dtype=float, na_value=np.nan)
    return matrix


def dosage_table(sv: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    """Retain record metadata beside the numeric sample dosage matrix."""
    return pd.concat(
        [
            sv[METADATA_COLUMNS].reset_index(drop=True),
            pd.DataFrame(dosage_matrix_from_genotypes(sv, samples), columns=samples),
        ],
        axis=1,
    )


def genotype_samples(sv: pd.DataFrame, metadata: pd.DataFrame) -> list[str]:
    """Samples with both a genotype column and population metadata, in metadata order."""
    genotype_columns = set(sv.columns[len(METADATA_COLUMNS):])
    samples = [sample for sample in metadata["sample_id"].astype(str) if sample in genotype_columns]
    if not samples:
        raise ValueError(
            "No sample in the sample metadata has a genotype column in the SV table; "
            "check that sample IDs match"
        )
    return samples


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted q-values, in the input order."""
    p_values = np.asarray(p_values, dtype=float)
    if p_values.size == 0:
        return p_values.copy()
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result
