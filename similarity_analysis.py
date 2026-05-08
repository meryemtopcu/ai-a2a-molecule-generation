"""Similarity analysis: generated ligands vs. ChEMBL reference actives.

Computes Morgan fingerprints (radius=2, 2048 bits) for both generated molecules
(Pocket2Mol and DiffSBDD) and 14 ChEMBL reference actives. Calculates pairwise
Tanimoto similarities and generates heatmaps + distributions.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
from tqdm import tqdm


def extract_smiles_from_download_csv(csv_path: Path) -> list[tuple[str, str]]:
    """Extract SMILES and molecule names from ChEMBL download CSV.
    
    Uses regex to parse the complex CSV format from ChEMBL web download.
    Returns list of (smiles, name) tuples.
    """
    smiles_list = []
    with csv_path.open("r", encoding="utf-8") as f:
        content = f.read()
    
    # Pattern: CHEMBL_ID;""Name"";Type;...;""SMILES""
    pattern = r'CHEMBL(\d+);""([^"]+)"";[^;]+;[^;]+;[^;]+;[^;]+;""([^"]+)""'
    matches = re.findall(pattern, content)
    
    for mol_id, name, smiles in matches:
        if smiles and any(c in smiles for c in "CNOFBrClI"):
            smiles_list.append((smiles, name))
    
    return smiles_list


def load_generated_molecules(csv_path: Path) -> list[tuple[str, str]]:
    """Load SMILES from generated molecule CSV (Pocket2Mol or DiffSBDD)."""
    smiles_list = []
    with csv_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader, start=1):
            if row and "smiles" in row and row["smiles"]:
                smiles = row["smiles"].strip()
                smiles_list.append((smiles, f"Gen_{i}"))
    return smiles_list


def load_sdf_molecules(source_path: Path) -> list[tuple[str, Chem.Mol]]:
    """Load generated molecules directly from one SDF file or a directory of SDF files."""
    path = Path(source_path)
    if not path.exists():
        raise ValueError(f"SDF source not found: {source_path}")

    files = [path] if path.is_file() else sorted(
        f for f in path.iterdir() if f.is_file() and f.suffix.lower() == ".sdf"
    )

    molecules: list[tuple[str, Chem.Mol]] = []
    for sdf_file in files:
        suppl = Chem.SDMolSupplier(str(sdf_file), sanitize=False, removeHs=False)
        for idx, mol in enumerate(suppl):
            if mol is not None:
                label = sdf_file.name if path.is_dir() else f"{sdf_file.stem}#{idx}"
                molecules.append((label, mol))
    return molecules


def load_reference_actives(reference_file: Path) -> list[tuple[str, Chem.Mol]]:
    """Load the 14 ChEMBL reference actives from the download CSV."""
    if not reference_file.is_file():
        raise ValueError(f"Reference file not found: {reference_file}")

    content = reference_file.read_text(encoding="utf-8")
    pattern = r'CHEMBL(\d+);""([^"]+)"";[^;]+;[^;]+;[^;]+;[^;]+;""([^"]+)""'

    molecules: list[tuple[str, Chem.Mol]] = []
    for mol_id, name, smiles in re.findall(pattern, content):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            molecules.append((f"{name} ({mol_id})", mol))

    return molecules


def get_morgan_fingerprint(mol: Chem.Mol):
    """Return a radius-2, 2048-bit Morgan fingerprint for a molecule."""
    try:
        mol_copy = Chem.Mol(mol)
        Chem.SanitizeMol(mol_copy)
        return AllChem.GetMorganFingerprintAsBitVect(mol_copy, radius=2, nBits=2048)
    except Exception:
        return None


def compute_morgan_fingerprints(
    smiles_list: list[tuple[str, str]],
    radius: int = 2,
    nbits: int = 2048,
) -> dict[str, list[Any]]:
    """Compute Morgan fingerprints for all SMILES.
    
    Returns dict with 'fingerprints' (list of fp objects) and 'valid_indices'.
    """
    fingerprints = []
    valid_indices = []
    
    for idx, (smiles, name) in enumerate(
        tqdm(smiles_list, desc="Computing fingerprints", unit="mol")
    ):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
                fingerprints.append(fp)
                valid_indices.append(idx)
        except Exception:
            pass
    
    return {
        "fingerprints": fingerprints,
        "valid_indices": valid_indices,
    }


def compute_similarity_matrix(
    generated_molecules: list[tuple[str, Chem.Mol]],
    reference_fps: list,
) -> np.ndarray:
    """Compute Tanimoto similarity matrix.
    
    Shape: (n_generated, n_reference)
    """
    n_ref = len(reference_fps)
    rows: list[list[float]] = []
    labels: list[str] = []

    for _, (label, mol) in enumerate(
        tqdm(generated_molecules, desc="Computing similarities", unit="generated_mol")
    ):
        gen_fp = get_morgan_fingerprint(mol)
        if gen_fp is None:
            continue
        rows.append([TanimotoSimilarity(gen_fp, ref_fp) for ref_fp in reference_fps])
        labels.append(label)

    if not rows:
        return np.zeros((0, n_ref)), []

    return np.array(rows, dtype=float), labels


def plot_combined_similarity_figure(
    pocket2mol_matrix: np.ndarray,
    diffsbdd_matrix: np.ndarray,
    pocket2mol_labels: list[str],
    diffsbdd_labels: list[str],
    reference_labels: list[str],
    output_path: Path,
) -> None:
    """Create one combined figure with both heatmaps and the max-Tanimoto distribution."""
    pocket2mol_max = pocket2mol_matrix.max(axis=1)
    diffsbdd_max = diffsbdd_matrix.max(axis=1)

    fig = plt.figure(figsize=(18, 14))
    grid = fig.add_gridspec(2, 2, height_ratios=[2.2, 1.0], hspace=0.25, wspace=0.15)

    pocket_ax = fig.add_subplot(grid[0, 0])
    diffsbdd_ax = fig.add_subplot(grid[0, 1])
    dist_ax = fig.add_subplot(grid[1, :])

    pocket_step = max(1, len(pocket2mol_labels) // 20)
    diffsbdd_step = max(1, len(diffsbdd_labels) // 20)

    sns.heatmap(
        pocket2mol_matrix,
        ax=pocket_ax,
        cmap="YlOrRd",
        vmin=0.0,
        vmax=1.0,
        cbar_kws={"label": "Tanimoto similarity"},
        xticklabels=reference_labels,
        yticklabels=pocket_step,
    )
    pocket_ax.set_title("Pocket2Mol vs. A2A reference actives")
    pocket_ax.set_xlabel("A2A reference actives")
    pocket_ax.set_ylabel("Pocket2Mol molecules")
    pocket_ax.tick_params(axis="x", labelrotation=90, labelsize=7)
    pocket_ax.tick_params(axis="y", labelsize=6)

    sns.heatmap(
        diffsbdd_matrix,
        ax=diffsbdd_ax,
        cmap="YlOrRd",
        vmin=0.0,
        vmax=1.0,
        cbar_kws={"label": "Tanimoto similarity"},
        xticklabels=reference_labels,
        yticklabels=diffsbdd_step,
    )
    diffsbdd_ax.set_title("DiffSBDD vs. A2A reference actives")
    diffsbdd_ax.set_xlabel("A2A reference actives")
    diffsbdd_ax.set_ylabel("DiffSBDD molecules")
    diffsbdd_ax.tick_params(axis="x", labelrotation=90, labelsize=7)
    diffsbdd_ax.tick_params(axis="y", labelsize=6)

    bins = np.linspace(0.0, 1.0, 21)
    dist_ax.hist(
        pocket2mol_max,
        bins=bins,
        alpha=0.65,
        color="steelblue",
        label=f"Pocket2Mol (n={len(pocket2mol_max)})",
    )
    dist_ax.hist(
        diffsbdd_max,
        bins=bins,
        alpha=0.65,
        color="coral",
        label=f"DiffSBDD (n={len(diffsbdd_max)})",
    )
    dist_ax.set_title("Max-Tanimoto distribution vs. A2A reference actives")
    dist_ax.set_xlabel("Maximum Tanimoto similarity")
    dist_ax.set_ylabel("Frequency")
    dist_ax.grid(axis="y", alpha=0.3)
    dist_ax.legend()

    fig.suptitle("Similarity Analysis vs. A2A", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved combined similarity figure: {output_path}")


def plot_heatmap(
    matrix: np.ndarray,
    title: str,
    output_path: Path,
    gen_labels: list[str],
    ref_labels: list[str],
) -> None:
    """Create and save heatmap of similarity matrix."""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    sns.heatmap(
        matrix,
        xticklabels=ref_labels,
        yticklabels=gen_labels,
        cmap="YlOrRd",
        cbar_kws={"label": "Tanimoto Similarity"},
        ax=ax,
        vmin=0,
        vmax=1,
    )
    
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("A2A Reference Actives", fontsize=12)
    ax.set_ylabel("Generated Molecules", fontsize=12)
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved heatmap: {output_path}")


def plot_max_tanimoto_distribution(
    pocket2mol_sims: np.ndarray,
    diffsbdd_sims: np.ndarray,
    output_path: Path,
) -> None:
    """Create histogram comparing max-Tanimoto distributions."""
    pocket2mol_max = pocket2mol_sims.max(axis=1)
    diffsbdd_max = diffsbdd_sims.max(axis=1)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    bins = np.linspace(0, 1, 21)
    ax.hist(
        pocket2mol_max,
        bins=bins,
        alpha=0.6,
        label=f"Pocket2Mol (n={len(pocket2mol_max)})",
        color="steelblue",
    )
    ax.hist(
        diffsbdd_max,
        bins=bins,
        alpha=0.6,
        label=f"DiffSBDD (n={len(diffsbdd_max)})",
        color="coral",
    )
    
    ax.set_xlabel("Maximum Tanimoto Similarity", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title("Distribution of Max-Tanimoto vs. ChEMBL Reference Actives", fontsize=14, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(axis="y", alpha=0.3)
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved distribution plot: {output_path}")


def print_statistics(pocket2mol_sims: np.ndarray, diffsbdd_sims: np.ndarray) -> None:
    """Print summary statistics for both tools."""
    print("\n" + "=" * 70)
    print("SIMILARITY ANALYSIS SUMMARY")
    print("=" * 70)
    
    p2m_max = pocket2mol_sims.max(axis=1)
    dff_max = diffsbdd_sims.max(axis=1)
    
    print(f"\nPocket2Mol (n={len(p2m_max)} molecules):")
    print(f"  Max-Tanimoto mean:   {p2m_max.mean():.4f}")
    print(f"  Max-Tanimoto median: {np.median(p2m_max):.4f}")
    print(f"  Max-Tanimoto std:    {p2m_max.std():.4f}")
    print(f"  Min:                 {p2m_max.min():.4f}")
    print(f"  Max:                 {p2m_max.max():.4f}")
    
    print(f"\nDiffSBDD (n={len(dff_max)} molecules):")
    print(f"  Max-Tanimoto mean:   {dff_max.mean():.4f}")
    print(f"  Max-Tanimoto median: {np.median(dff_max):.4f}")
    print(f"  Max-Tanimoto std:    {dff_max.std():.4f}")
    print(f"  Min:                 {dff_max.min():.4f}")
    print(f"  Max:                 {dff_max.max():.4f}")
    
    print("\n" + "=" * 70)


def main() -> None:
    """Run full similarity analysis from direct SDF inputs."""
    base_dir = Path.cwd()
    out_dir = base_dir / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    pocket2mol_source = (
        base_dir
        / "outputs2"
        / "sample_for_pdb_4EIY_clean.pdb_2026_05_02__11_02_39"
        / "SDF"
    )
    diffsbdd_source = base_dir / "outputs2" / "combined_diffsbdd.sdf"
    reference_file = (
        base_dir
        / "ChEMBL reference molecules"
        / "DOWNLOAD-muQF1iIRyjhMqiuLdRuitIuq-3YXqGwME0QTLPrZB4E_eq_.csv"
    )

    print("Loading 14 ChEMBL reference actives...")
    reference_molecules = load_reference_actives(reference_file)
    print(f"Loaded {len(reference_molecules)} reference molecules")

    print("\nLoading generated molecules from SDF files...")
    pocket2mol_molecules = load_sdf_molecules(pocket2mol_source)
    diffsbdd_molecules = load_sdf_molecules(diffsbdd_source)
    print(f"Loaded {len(pocket2mol_molecules)} Pocket2Mol molecules")
    print(f"Loaded {len(diffsbdd_molecules)} DiffSBDD molecules")

    print("\nComputing Morgan fingerprints (radius=2, 2048 bits)...")
    reference_fps = []
    reference_labels = []
    for label, mol in tqdm(reference_molecules, desc="Reference fingerprints", unit="mol"):
        fp = get_morgan_fingerprint(mol)
        if fp is not None:
            reference_fps.append(fp)
            reference_labels.append(label)

    reference_data = {"fingerprints": reference_fps, "names": reference_labels}
    print(f"Valid fingerprints: {len(reference_fps)} ref")

    p2m_matrix, p2m_labels = compute_similarity_matrix(pocket2mol_molecules, reference_fps)
    dff_matrix, dff_labels = compute_similarity_matrix(diffsbdd_molecules, reference_fps)
    print(f"Valid fingerprints: {len(p2m_labels)} Pocket2Mol, {len(dff_labels)} DiffSBDD")

    if p2m_matrix is None or dff_matrix is None:
        print("Error: could not compute similarity matrices.")
        return

    plot_combined_similarity_figure(
        pocket2mol_matrix=p2m_matrix,
        diffsbdd_matrix=dff_matrix,
        pocket2mol_labels=p2m_labels,
        diffsbdd_labels=dff_labels,
        reference_labels=reference_labels,
        output_path=out_dir / "similarity_vs_known_actives.png",
    )

    print_statistics(p2m_matrix, dff_matrix)
    print(f"\nAnalysis complete! Result saved to {out_dir / 'similarity_vs_known_actives.png'}")


if __name__ == "__main__":
    main()
