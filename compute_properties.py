"""Compute properties for molecules in PATH1 and PATH2, save them to a file, and compare them.

This script analyzes molecular structures from two sources (Pocket2Mol and DiffSBDD) and computes
pharmacological properties plus structural similarity to known A2A adenosine receptor actives.

Workflow:
1. Load reference A2A active molecules and compute Morgan fingerprints (radius=2, fpSize=2048)
2. Load molecules from PATH1 (Pocket2Mol, optional) and PATH2 (DiffSBDD, required)
3. Validate molecules: skip any that raise RDKit exceptions during property computation
4. Compute for valid molecules:
   - SMILES representation
   - Molecular Weight (MW): 150-500 Da ideal
   - LogP (partition coefficient): -0.4 to 5.6 ideal
   - HBA (H-bond acceptors): ≤10 ideal
   - HBD (H-bond donors): ≤5 ideal
   - TPSA (topological polar surface area): ≤140 Å² ideal
   - RotBonds (rotatable bonds)
   - QED (quantitative estimate of drug-likeness): >0.5 ideal
   - SA Score (synthetic accessibility): ≤4 ideal (lower = easier to synthesize)
   - Lipinski's Rule of Five (binary sum: 0-5 violations)
    - Tanimoto similarity vs. A2A reference set (max of all pairwise comparisons)
5. Save results to CSV files (separate for PATH1 and PATH2)
6. Compare properties between paired molecules if both paths present

Input Requirements:
- PATH1: outputs2/sample_for_pdb_4EIY_clean.pdb_2026_05_02__11_02_39/SDF
    (Pocket2Mol SDF files, default)
- PATH2: outputs2/4eiy_inpaint.sdf (100 DiffSBDD molecules, default)
- A2A ref: outputs2/a2a_drugs/ADORA2A-world.sdf
    (29 known actives, default)

Output:
- pocket2mol_generated.csv (always created)
- diffsbdd_generated.csv (always created)
CSV columns: smiles, valid, MW, LogP, HBA, HBD, TPSA, RotBonds, QED,
             SA_score, Lipinski, max_tanimoto
"""

#  Standard library imports
import argparse
import csv
import sys
import os
from pathlib import Path

import numpy as np

# pylint: disable=broad-exception-caught

#  RDKit imports for molecular processing
from rdkit import Chem, DataStructs  # Core RDKit: molecular I/O, fingerprint comparison
from rdkit.Chem import (
    Crippen,  # LogP computation
    Descriptors,  # Molecular weight, H-count descriptors
    Lipinski,  # H-donor/acceptor counting (Lipinski rule)
    QED,  # Drug-likeness scoring
    rdFingerprintGenerator,  # Morgan fingerprint generation
    rdMolDescriptors  # TPSA, rotatable bonds counting
)
from combine_sdfs import combine_sdfs

#  Import SA Score calculator (Synthetic Accessibility Score)
sa_score_path = Path(__file__).resolve().parent / 'DiffSBDD-main' / 'analysis' / 'SA_Score'
sys.path.insert(0, str(sa_score_path))
try:
    from sascorer import calculateScore as sa_calculateScore  # pylint: disable=import-error  # type: ignore[import-not-found]
except ImportError:
    print("Warning: sascorer module not found, SA Score will be set to None")
    sa_calculateScore = None

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
except ImportError:
    plt = None
    sns = None


#  ============================================================================
#  Global Configuration
#  ============================================================================
#  Morgan fingerprint generator: radius=2 (includes 2-hop neighbors), fpSize=2048 bits
#  Used for Tanimoto similarity computation (requires consistent fingerprinting)
MORGAN_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


#  ============================================================================
#  Core Property Computation
#  ============================================================================
def _compute_all_properties(mol):
    """Compute all molecular properties in one function to eliminate duplication.
    
    This is a private helper (leading underscore) called by both is_valid_molecule()
    and compute_properties(). Centralizing the logic reduces maintenance burden and
    ensures consistent property calculations.
    
    Args:
        mol: RDKit Mol object (assumed to be valid/non-None)
        
    Returns:
        dict: Property dictionary with keys: MolecularWeight, LogP, NumHAcceptors,
              NumHDonors, TPSA, RotBonds, QED, SA_score, LipinskiRuleOfFive
              
    Raises:
        Exception: If any RDKit descriptor computation fails (e.g., invalid valence,
                   sanitization issues). Caller handles via try/except.
    """
    properties = {}
    #  Molecular weight (g/mol)
    properties['MolecularWeight'] = Descriptors.MolWt(mol)
    #  Lipophilicity (log P, oil-water partition coefficient)
    properties['LogP'] = Crippen.MolLogP(mol)
    #  Hydrogen bond acceptors (e.g., N, O)
    properties['NumHAcceptors'] = Descriptors.NumHAcceptors(mol)
    #  Hydrogen bond donors (N-H, O-H)
    properties['NumHDonors'] = Descriptors.NumHDonors(mol)
    #  Topological Polar Surface Area (Ångström²)
    properties['TPSA'] = rdMolDescriptors.CalcTPSA(mol)
    #  Rotatable bonds (flexible bonds that affect conformation)
    properties['RotBonds'] = rdMolDescriptors.CalcNumRotatableBonds(mol)
    #  Quantitative Estimate of Drug-Likeness (0-1 scale, higher is better)
    properties['QED'] = QED.qed(mol)
    #  Synthetic Accessibility Score (1-10, lower = easier to synthesize, ideal ≤4)
    properties['SA_score'] = compute_sa_score(mol)
    #  Lipinski's Rule of Five: count satisfied rules (higher = more satisfied)
    properties['LipinskiRuleOfFive'] = int(
        Descriptors.MolWt(mol) < 500,
    ) + int(Lipinski.NumHDonors(mol) <= 5) + int(Lipinski.NumHAcceptors(mol) <= 10) + int(
        Crippen.MolLogP(mol) <= 5
    ) + int(rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10)
    return properties


def is_valid_molecule(mol):
    """Check if molecule can be processed without RDKit errors.
    
    Used as a pre-filter before property computation. Some molecules from generative
    models have invalid valence, missing implicit hydrogens, or other RDKit issues.
    Attempting to compute properties (especially QED) will raise exceptions on these.
    This function catches those exceptions and returns False.
    
    Args:
        mol: RDKit Mol object or None
        
    Returns:
        bool: True if all properties can be computed, False otherwise
    """
    try:
        _compute_all_properties(mol)
        return True
    except Exception:  # Broad except is intentional: catch all RDKit errors
        return False


def compute_sa_score(mol):
    """Compute Synthetic Accessibility Score (1-10, lower=easier to synthesize).
    
    SA Score ranges from 1 (easy to synthesize) to 10 (hard to synthesize).
    Drug-like compounds typically have SA Score ≤ 4.
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        float: SA Score (1-10), or None if sascorer unavailable
    """
    if sa_calculateScore is None:
        return None
    try:
        return sa_calculateScore(mol)
    except Exception:
        return None


def get_smiles(mol):
    """Generate SMILES string from molecule.
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        str: SMILES string, or None on failure
    """
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def compute_properties(mol):
    """Compute all molecular properties for a valid molecule.
    
    Assumes mol has already been validated via is_valid_molecule(). If computation
    fails, returns None so caller can skip that molecule without crashing.
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        dict: Property dictionary, or None if any RDKit operation fails
    """
    try:
        return _compute_all_properties(mol)
    except Exception:  # Fallback for edge cases not caught by is_valid_molecule()
        return None


#  ============================================================================
#  Fingerprinting & Similarity Computation
#  ============================================================================
def compute_tanimoto_similarity(mol, reference_data):
    """Compute maximum Tanimoto similarity vs. known A2A actives.
    
    Tanimoto similarity (0-1) measures structural overlap: 1.0 = identical fingerprint,
    0.0 = completely different. We return the max to identify the closest known active.
    
    Args:
        mol: RDKit Mol object
        reference_data: dict with fingerprints and names from
            load_a2a_reference_fingerprints()
        
    Returns:
        float: maximum similarity, or None on failure
    """
    try:
        if not reference_data or not reference_data.get('fingerprints'):
            return None
        #  Generate fingerprint for query molecule
        fp = get_morgan_fingerprint(mol)
        if fp is None:
            return None
        #  Compute similarity to each reference fingerprint
        similarities = []
        for ref_fp in reference_data['fingerprints']:
            sim = DataStructs.TanimotoSimilarity(fp, ref_fp)
            similarities.append(sim)

        if not similarities:
            return None

        return max(similarities)
    except Exception:
        return None


def get_morgan_fingerprint(mol):
    """Generate a Morgan fingerprint from a molecule.
    
    Creates a copy of the molecule and sanitizes it before fingerprinting.
    Sanitization ensures all valence/aromaticity assumptions are consistent.
    We operate on a copy to avoid modifying the original molecule object.
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        ExplicitBitVect: 2048-bit Morgan fingerprint (radius=2), or None
        on error
    """
    try:
        #  Create a copy to avoid side effects
        mol_copy = Chem.Mol(mol)
        #  Sanitize: compute aromaticity, valence, implicit H's, etc.
        Chem.SanitizeMol(mol_copy)
        #  Generate fingerprint using the pre-configured MORGAN_GENERATOR
        return MORGAN_GENERATOR.GetFingerprint(mol_copy)
    except Exception:  # Sanitization may fail on partially valid molecules
        return None


def load_a2a_reference_fingerprints(a2a_sdf_file):
    """Load A2A reference molecules and compute Morgan fingerprints.
    
    The A2A reference set (ADORA2A-world.sdf from ZINC15) contains ~29 known
    adenosine receptor agonists/modulators. We use these to compute structural
    similarity for generated molecules.
    
    Args:
        a2a_sdf_file: path to ADORA2A-world.sdf
        
    Returns:
        dict with fingerprints (list of ExplicitBitVect) and names (list of str)
        for valid molecules
        
    Raises:
        ValueError: if SDF file not found
    """
    sdf_path = Path(a2a_sdf_file)
    if not sdf_path.is_file():
        raise ValueError(f'A2A SDF file not found: {a2a_sdf_file}')

    fingerprints = []
    names = []
    #  Load SDF without automatic sanitization to allow partial molecules
    suppl = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
    for idx, mol in enumerate(suppl):
        #  Skip None (parsing errors) and invalid molecules
        if mol is None or not is_valid_molecule(mol):
            continue
        fp = get_morgan_fingerprint(mol)
        #  Only store fingerprints that were successfully computed
        if fp is not None:
            fingerprints.append(fp)
            # Prefer a stable identifier: use 'zinc_id' if available, then _Name, then SMILES, then index
            mol_name = None
            try:
                if mol.HasProp('zinc_id'):
                    mol_name = mol.GetProp('zinc_id')
                elif mol.HasProp('_Name'):
                    mol_name = mol.GetProp('_Name')
                else:
                    # Fallback to canonical SMILES where possible
                    mol_name = Chem.MolToSmiles(mol) if Chem.MolToSmiles(mol) else None
            except Exception:
                mol_name = None
            if not mol_name:
                mol_name = f'A2A_{idx}'
            names.append(mol_name)

    return {'fingerprints': fingerprints, 'names': names}


#  ============================================================================
#  Molecule Loading (Input Handlers)
#  ============================================================================
def load_molecules_from_path1(path1):
    """Load Pocket2Mol output: SDF files from a directory.

    Pocket2Mol generates one SDF file per sampled molecule. This function loads
    all *.sdf files in the directory, extracting the first non-None molecule from
    each file (expected: one molecule per file).
    
    Args:
        path1: Path to directory containing Pocket2Mol SDF files
        
    Returns:
        list: A list of (name, Mol) tuples containing the filename and RDKit molecule object.

    Raises:
        ValueError: if path1 is not a directory
    """
    path1 = Path(path1)
    if not path1.is_dir():
        raise ValueError(f'PATH1 must be a directory containing SDF files: {path1}')

    molecules = []
    #  Find any file containing '.sdf' (case-insensitive)
    files = sorted([f for f in path1.iterdir() if f.is_file() and '.sdf' in f.name.lower()])
    for sdf_file in files:
        #  Load without sanitization to allow pre-validation
        suppl = Chem.SDMolSupplier(str(sdf_file), sanitize=False, removeHs=False)
        #  Extract first valid molecule
        mol = next((m for m in suppl if m is not None), None)
        if mol is not None:
            molecules.append((sdf_file.name, mol))
    return molecules


def load_molecules_from_path2(path2):
    """Load DiffSBDD output from a single SDF file or a directory of SDF files.

    DiffSBDD usually generates one SDF file containing multiple molecules, but
    if the output is split, this function can load from all .sdf files in a directory.

    Args:
        path2: Path to a single SDF file or a directory containing SDF files.

    Returns:
        list of (name, Mol) tuples where name is 'filename#index'.

    Raises:
        ValueError: if path2 is not a file or directory.
    """
    path = Path(path2)
    files = []
    if path.is_file():
        files.append(path)
    elif path.is_dir():
        # Broaden search to find any file containing '.sdf' (e.g. .sdf_copy)
        files = sorted([f for f in path.iterdir() if f.is_file() and '.sdf' in f.name.lower()])
    else:
        raise ValueError(f'PATH2 must be an SDF file or a directory: {path2}')

    molecules = []
    for sdf_file in files:
        #  Load multi-molecule SDF without sanitization
        suppl = Chem.SDMolSupplier(str(sdf_file), sanitize=False, removeHs=False)
        for idx, mol in enumerate(suppl):
            if mol is not None:
                #  Label: filename#index
                molecules.append((f'{sdf_file.name}#{idx}', mol))
    return molecules


def _fragment_count(mol):
    """Return the number of disconnected fragments in a molecule."""
    try:
        return len(Chem.GetMolFrags(mol, asMols=False))
    except Exception:
        return None


def report_fragmentation(molecules, source_name):
    """Print how many loaded molecules are fragmented for a source."""
    per_file = {}
    total_fragmented = 0
    total_loaded = len(molecules)

    for name, mol in molecules:
        file_name = name.split('#', 1)[0]
        fragment_count = _fragment_count(mol)
        if fragment_count is None:
            per_file.setdefault(file_name, {'fragmented': 0, 'loaded': 0, 'invalid': 0})
            per_file[file_name]['invalid'] += 1
            continue

        per_file.setdefault(file_name, {'fragmented': 0, 'loaded': 0, 'invalid': 0})
        per_file[file_name]['loaded'] += 1
        if fragment_count > 1:
            per_file[file_name]['fragmented'] += 1
            total_fragmented += 1

    if not molecules:
        print(f'{source_name}: no molecules loaded')
        return {'total_loaded': 0, 'total_fragmented': 0, 'per_file': per_file}

    print(f'\nFragmentation check for {source_name}:')
    for file_name in sorted(per_file):
        counts = per_file[file_name]
        if counts['invalid'] and not counts['loaded']:
            print(f'  {file_name}: {counts["invalid"]} invalid molecules, no fragment count available')
        else:
            print(
                f'  {file_name}: {counts["fragmented"]} fragmented molecule(s) '
                f'out of {counts["loaded"]} loaded'
            )
    print(f'  Total fragmented in {source_name}: {total_fragmented} / {total_loaded}')
    return {'total_loaded': total_loaded, 'total_fragmented': total_fragmented, 'per_file': per_file}


def filter_fragmented_molecules(molecules):
    """Keep only molecules that are a single disconnected fragment."""
    filtered = []
    removed = 0
    for name, mol in molecules:
        if _fragment_count(mol) == 1:
            filtered.append((name, mol))
        else:
            removed += 1
    return filtered, removed


def prompt_fragmentation_action(path1_molecules, path2_molecules):
    """Ask whether to continue with fragmented molecules or filter them out."""
    report_fragmentation(path1_molecules, 'PATH1')
    report_fragmentation(path2_molecules, 'PATH2')
    print('\nFragmented molecules detected.')
    print('1. Run anyway')
    print('2. Run without fragmented molecules')
    print('3. QUIT')

    while True:
        choice = input('Select an option [1/2/3]: ').strip().lower()
        if choice in {'1', 'run anyway', 'run', 'y', 'yes'}:
            return 'run_anyway'
        if choice in {'2', 'run without fragmented', 'filter', 'filtered'}:
            return 'filter_fragmented'
        if choice in {'3', 'quit', 'q', 'exit'}:
            return 'quit'
        print('Please enter 1, 2, or 3.')


#  ============================================================================
#  I/O and Comparison Utilities
#  ============================================================================
def save_properties_to_csv(molecules_data, output_file, append=False):
    """Save computed properties to CSV file.
    
    Args:
        molecules_data: list of dicts with keys: smiles, valid, MW, LogP, HBA, HBD, TPSA,
                   RotBonds, QED, SA_score, Lipinski, max_tanimoto
        output_file: path to output CSV file
    """
    fieldnames = ['smiles', 'valid', 'MW', 'LogP', 'HBA', 'HBD', 'TPSA', 'RotBonds',
                  'QED', 'SA_score', 'Lipinski', 'max_tanimoto']

    #  Default behavior: overwrite file. If append=True, append rows and skip
    #  header when file exists.
    def _write_rows(mode, write_header):
        with open(output_file, mode, newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            if write_header:
                writer.writeheader()
            for row in molecules_data:
                writer.writerow(row)

    file_exists = os.path.exists(output_file)
    if append:
        #  Append: if file exists, skip header; otherwise write header then rows
        mode = 'a'
        write_header = not file_exists
    else:
        mode = 'w'
        write_header = True

    _write_rows(mode, write_header)
    if molecules_data:
        action = 'Appended' if append and file_exists else 'Saved'
        print(f"{action} {len(molecules_data)} molecules to {output_file}")
    else:
        if write_header:
            print(f"Created empty CSV header in {output_file}")
        else:
            print(f"No rows to append to {output_file}")


def compare_properties(properties1, properties2):
    """Compare two property dicts and return all differences.
    
    Checks all keys present in either dict (symmetric comparison).
    Returns differences for keys with mismatched values or missing from one dict.
    
    Args:
        properties1, properties2: dicts with property keys/values
        
    Returns:
        dict: {key: (value1, value2)} for differing properties or missing keys
    """
    differences = {}
    #  Get all unique keys from both dicts
    all_keys = set(properties1.keys()) | set(properties2.keys())

    for key in all_keys:
        val1 = properties1.get(key)
        val2 = properties2.get(key)
        #  Record difference if values differ (handles missing keys too)
        if val1 != val2:
            differences[key] = (val1, val2)

    return differences


def compare_property_lists(properties1, properties2):
    """Compare molecules from PATH1 and PATH2.
    
    Assumes molecules are indexed in the same order. Compares properties
    for each pair and collects differences.
    
    Args:
        properties1, properties2: lists of {'name': str, 'properties': dict}
        
    Returns:
        list of dictionaries with index, names, and differences
    """
    differences = []
    #  Pair up molecules (stop at shorter list)
    n_pairs = min(len(properties1), len(properties2))
    for idx in range(n_pairs):
        diff = compare_properties(properties1[idx]['properties'], properties2[idx]['properties'])
        if diff:
            differences.append({
                'index': idx,
                'path1_name': properties1[idx]['name'],
                'path2_name': properties2[idx]['name'],
                'differences': diff,
            })
    return differences


def load_similarity_rows(csv_file):
    """Load max-Tanimoto values from a CSV file.
    
    Args:
        csv_file: Path to the CSV file containing molecular properties.
        
    Returns:
        list: A list of float values representing max_tanimoto similarity.
    """
    rows = []
    if not os.path.exists(csv_file):
        return []
    with open(csv_file, newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            value = row.get('max_tanimoto')
            if value in (None, ''):
                continue
            try:
                rows.append(float(value))
            except ValueError:
                continue
    return rows


def compute_similarity_matrix(molecules, reference_data):
    """Compute a generated-vs-reference Tanimoto matrix for heatmap plotting.
    
    Args:
        molecules: List of (name, Mol) tuples.
        reference_data: Dict containing 'fingerprints' and 'names' of reference actives.
        
    Returns:
        tuple: (numpy.ndarray matrix, list labels, list reference_labels)
    """
    if not molecules or not reference_data or not reference_data.get('fingerprints'):
        return None, [], []

    reference_fingerprints = reference_data['fingerprints']
    reference_names = reference_data.get('names', [])
    matrix = []
    labels = []

    for name, mol in molecules:
        fp = get_morgan_fingerprint(mol)
        if fp is None:
            continue
        row = [DataStructs.TanimotoSimilarity(fp, ref_fp) for ref_fp in reference_fingerprints]
        matrix.append(row)
        labels.append(name)

    if not matrix:
        return None, [], reference_names

    return np.array(matrix, dtype=float), labels, reference_names


def save_similarity_matrix_csv(matrix, row_labels, col_labels, out_file):
    """Save a numeric similarity matrix to CSV.

    The CSV will have the reference labels as columns and the generated molecule
    labels as the first column.
    """
    out_path = Path(out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w', newline='', encoding='utf-8') as fout:
        writer = csv.writer(fout)
        # Prepare column labels: validate provided col_labels or fall back to ref_i
        n_cols = matrix.shape[1]
        if col_labels and isinstance(col_labels, (list, tuple)) and len(col_labels) == n_cols:
            header_cols = [lbl if lbl else f'ref_{i}' for i, lbl in enumerate(col_labels)]
        else:
            header_cols = [f'ref_{i}' for i in range(n_cols)]
        header = ['generated'] + header_cols
        writer.writerow(header)
        for i, row in enumerate(matrix):
            label = row_labels[i] if row_labels and i < len(row_labels) and row_labels[i] else f'gen_{i}'
            writer.writerow([label] + [f'{v:.6f}' for v in row])


def plot_similarity_comparison(
    pocket2mol_csv,
    diffsbdd_csv,
    pocket2mol_matrix=None,
    diffsbdd_matrix=None,
    pocket2mol_labels=None,
    diffsbdd_labels=None,
    reference_labels=None,
    output_dir=None,
):
    """Create the tool-comparison plots and full molecule-vs-reference heatmaps.
    
    Args:
        pocket2mol_csv: Path to Pocket2Mol results CSV.
        diffsbdd_csv: Path to DiffSBDD results CSV.
        pocket2mol_matrix: Tanimoto similarity matrix for Pocket2Mol.
        diffsbdd_matrix: Tanimoto similarity matrix for DiffSBDD.
        pocket2mol_labels: Labels for Pocket2Mol molecules.
        diffsbdd_labels: Labels for DiffSBDD molecules.
        reference_labels: Names of reference actives.
        output_dir: Directory to save generated PNG plots.
        
    Returns:
        dict: Paths to the generated plot files.
    """
    if plt is None or sns is None:
        print('Skipping similarity plots because matplotlib/seaborn is unavailable.')
        return None

    pocket2mol_values = load_similarity_rows(pocket2mol_csv) if pocket2mol_csv else []
    diffsbdd_values = load_similarity_rows(diffsbdd_csv) if diffsbdd_csv else []
    if not pocket2mol_values and not diffsbdd_values:
        print('Skipping similarity plots because no max_tanimoto values were found.')
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    if pocket2mol_values:
        sns.histplot(pocket2mol_values, bins=20, ax=axes[0], color='#1f77b4', kde=True)
        axes[0].set_title('Pocket2Mol max Tanimoto')
        axes[0].set_xlabel('max_tanimoto')
    else:
        axes[0].axis('off')
        axes[0].set_title('Pocket2Mol max Tanimoto (no data)')

    if diffsbdd_values:
        sns.histplot(diffsbdd_values, bins=20, ax=axes[1], color='#d62728', kde=True)
        axes[1].set_title('DiffSBDD max Tanimoto')
        axes[1].set_xlabel('max_tanimoto')
    else:
        axes[1].axis('off')
        axes[1].set_title('DiffSBDD max Tanimoto (no data)')

    for axis in axes:
        axis.set_ylabel('Count')

    fig.suptitle('Max-Tanimoto distribution vs. known A2A actives')
    fig.tight_layout()
    distribution_file = output_dir / 'max_tanimoto_distribution.png'
    fig.savefig(distribution_file, dpi=200, bbox_inches='tight')
    plt.close(fig)

    heatmap_bins = np.linspace(0.0, 1.0, 21)
    pocket_hist, _ = np.histogram(pocket2mol_values or [np.nan], bins=heatmap_bins)
    diffsbdd_hist, _ = np.histogram(diffsbdd_values or [np.nan], bins=heatmap_bins)
    matrix = np.vstack([pocket_hist, diffsbdd_hist])
    fig, ax = plt.subplots(figsize=(10, 2.8))
    sns.heatmap(
        matrix,
        ax=ax,
        cmap='viridis',
        annot=True,
        fmt='d',
        cbar_kws={'label': 'count'},
        yticklabels=['Pocket2Mol', 'DiffSBDD'],
        xticklabels=[f'{edge:.2f}' for edge in heatmap_bins[:-1]],
    )
    ax.set_title('Tool comparison heatmap of max-Tanimoto counts')
    ax.set_xlabel('max_tanimoto bin start')
    heatmap_file = output_dir / 'max_tanimoto_heatmap.png'
    fig.tight_layout()
    fig.savefig(heatmap_file, dpi=200, bbox_inches='tight')
    plt.close(fig)

    reference_labels = reference_labels or []
    matrix_plots = []
    for tool_name, matrix, labels in (
        ('Pocket2Mol', pocket2mol_matrix, pocket2mol_labels),
        ('DiffSBDD', diffsbdd_matrix, diffsbdd_labels),
    ):
        if matrix is None or len(matrix) == 0:
            continue

        fig_width = max(10, len(reference_labels) * 0.45)
        fig_height = max(4, len(labels) * 0.28)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        sns.heatmap(
            matrix,
            ax=ax,
            cmap='mako',
            vmin=0.0,
            vmax=1.0,
            cbar_kws={'label': 'Tanimoto similarity'},
            xticklabels=reference_labels or True,
            yticklabels=labels,
        )
        ax.set_title(f'{tool_name} vs reference actives')
        ax.set_xlabel('Reference actives')
        ax.set_ylabel('Generated molecules')
        fig.tight_layout()
        matrix_file = output_dir / f'{tool_name.lower()}_vs_reference_heatmap.png'
        fig.savefig(matrix_file, dpi=200, bbox_inches='tight')
        plt.close(fig)
        matrix_plots.append(str(matrix_file))

    print(f'Saved similarity plots to {distribution_file} and {heatmap_file}')
    if matrix_plots:
        print('Saved full similarity heatmaps to ' + ', '.join(matrix_plots))
    return {
        'distribution': str(distribution_file),
        'heatmap': str(heatmap_file),
        'matrix_plots': matrix_plots,
    }

#  Duplicate filtering

def filter_duplicates(molecules):
    """Remove duplicate molecules based on canonical SMILES.

    Keeps the first occurrence of each unique SMILES. Returns the filtered
    list and the count of removed duplicates.
    """
    seen = set()
    unique = []
    removed = 0
    for name, mol in molecules:
        try:
            smi = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            #  If SMILES generation fails, treat molecule as unique (validation will remove later)
            unique.append((name, mol))
            continue
        if smi in seen:
            removed += 1
            continue
        seen.add(smi)
        unique.append((name, mol))
    return unique, removed

#  ============================================================================
#  Main Workflow
#  ============================================================================
def main():
    """Main entry point: parse CLI args, load molecules, compute properties.
    
    Workflow:
    1. Parse argparse with defaults pointing to outputs2 directory
    2. Load A2A reference fingerprints
    3. Load PATH1 (optional) and PATH2 (required) molecules
    4. Validate molecules, report invalid counts
    5. Compute properties for valid molecules
    6. Append Tanimoto similarity to each molecule
    7. Save results to CSV files
    8. Output separate CSVs for PATH1 and PATH2
    """
    #  Construct default paths relative to script location for portability
    script_dir = Path(__file__).resolve().parent
    default_out_src = script_dir / 'results'  # Where outputs (CSVs, plots) go
    default_in_src = script_dir / 'outputs2'  # Where input data is located
    parser = argparse.ArgumentParser(
        description='Compute drug-like properties and A2A similarity for generated molecules'
    )
    #  Optional Pocket2Mol input
    parser.add_argument(
        '--path1',
        default=None,
        help=(
            'Directory of single-molecule SDF files from Pocket2Mol '
            '(optional). Provide this explicitly to process PATH1.'
        )
    )
    #  DiffSBDD input (opt-in). If neither --path1 nor --path2 are provided,
    #  both defaults under `outputs2` will be used for backward compatibility.
    parser.add_argument(
        '--path2',
        default=None,
        help='Multi-molecule SDF file from DiffSBDD. Provide explicitly to process PATH2.'
    )
    #  A2A reference fingerprints
    parser.add_argument(
        '--a2a',
        dest='a2a_sdf_file',
        default=str(default_in_src / 'a2a_drugs' / 'ADORA2A-world.sdf'),
        help='A2A actives reference SDF file for similarity computation.'
    )
    parser.add_argument(
        '--append',
        action='store_true',
        help='Append PATH2 results to diffsbdd_generated.csv instead of overwriting'
    )
    parser.add_argument(
        '--plots-only',
        action='store_true',
        help='Generate similarity plots without writing the CSV outputs again'
    )
    parser.add_argument(
        '--outdir',
        default=None,
        help='Directory to save output CSV files and plots (default: results/)'
    )
    parser.add_argument(
        '--csv1',
        default=None,
        help='Path to existing Pocket2Mol CSV file for plotting (replaces default in outdir).'
    )
    parser.add_argument(
        '--csv2',
        default=None,
        help='Path to existing DiffSBDD CSV file for plotting (replaces default in outdir).'
    )
    #  (No JSON output) CSV files are always written to results/
    args = parser.parse_args()

    output_dir = Path(args.outdir) if args.outdir else default_out_src
    output_dir.mkdir(parents=True, exist_ok=True)

    path1 = args.path1
    path2 = args.path2
    a2a_sdf_file = args.a2a_sdf_file

    #  Backward-compatible defaults: if the user did not provide either path,
    #  fall back to the original outputs2 defaults for both PATH1 and PATH2.
    if path1 is None and path2 is None and not (args.csv1 or args.csv2):
        path1 = str(default_in_src / 'sample_for_pdb_4EIY_clean.pdb_2026_05_02__11_02_39' / 'SDF')
        path2 = str(default_in_src / '4eiy_inpaint.sdf')

    #  === PRE-STEP: Combine split SDFs if PATH2 is a directory ===
    if path2 and Path(path2).is_dir():
        combined_sdf_path = output_dir / 'combined_diffsbdd.sdf'
        print(f"Directory detected for PATH2. Combining split SDFs into {combined_sdf_path}...")
        if combine_sdfs(path2, combined_sdf_path):
            path2 = str(combined_sdf_path)
        else:
            print(f"Warning: No valid molecules found to combine in {path2}")

    #  === STEP 1: Load A2A Reference Fingerprints ===
    reference_data = load_a2a_reference_fingerprints(a2a_sdf_file)
    if not reference_data or not reference_data.get('fingerprints'):
        print('Error: Could not load any valid known A2A actives from ADORA2A-world.sdf.')
        return
    print(f'Loaded {len(reference_data["fingerprints"])} known A2A active fingerprints')

    #  === STEP 2: Load Molecules ===
    path1_molecules = []
    if path1:
        if Path(path1).is_dir():
            path1_molecules = load_molecules_from_path1(path1)
            print(f"Loaded {len(path1_molecules)} molecules from PATH1")
        else:
            print(f"PATH1 not found or not a directory, skipping: {path1}")
    else:
        print("PATH1 not provided, skipping Pocket2Mol analysis")

    path2_molecules = []
    if path2:
        path2_path = Path(path2)
        if path2_path.exists():
            path2_molecules = load_molecules_from_path2(path2)
            print(f"Loaded {len(path2_molecules)} molecules from PATH2")
        else:
            print(f"PATH2 not found, skipping: {path2}")
    else:
        print("PATH2 not provided, skipping DiffSBDD analysis")

    if not path1_molecules and not path2_molecules:
        print('Error: Could not load molecules from PATH1 or PATH2.')
        return

    #  === STEP 2.25: Check for fragmented molecules ===
    fragmentation_choice = None
    if path1_molecules or path2_molecules:
        path1_fragmented = sum(1 for _, mol in path1_molecules if _fragment_count(mol) not in (None, 1))
        path2_fragmented = sum(1 for _, mol in path2_molecules if _fragment_count(mol) not in (None, 1))
        if path1_fragmented or path2_fragmented:
            fragmentation_choice = prompt_fragmentation_action(path1_molecules, path2_molecules)
            if fragmentation_choice == 'quit':
                print('Quitting without running property calculations.')
                return
            if fragmentation_choice == 'filter_fragmented':
                if path1_molecules:
                    path1_molecules, removed1 = filter_fragmented_molecules(path1_molecules)
                    print(f'Removed {removed1} fragmented molecule(s) from PATH1')
                if path2_molecules:
                    path2_molecules, removed2 = filter_fragmented_molecules(path2_molecules)
                    print(f'Removed {removed2} fragmented molecule(s) from PATH2')

    #  === STEP 2.5: Remove Duplicates ===
    #  Filter duplicates based on SMILES strings before validation
    path1_duplicates = 0
    path2_duplicates = 0
    if path1_molecules:
        path1_molecules, path1_duplicates = filter_duplicates(path1_molecules)
        if path1_duplicates > 0:
            print(f"Removed {path1_duplicates} duplicate molecules from PATH1")

    path2_molecules, path2_duplicates = filter_duplicates(path2_molecules)
    if path2_duplicates > 0:
        print(f"Removed {path2_duplicates} duplicate molecules from PATH2")

    #  === STEP 3: Validate Molecules ===
    #  Skip molecules that raise RDKit errors during property computation
    print("\nValidating molecules...")
    valid1 = (
        [(name, mol) for name, mol in path1_molecules if is_valid_molecule(mol)]
        if path1_molecules
        else []
    )
    invalid1_count = len(path1_molecules) - len(valid1) if path1_molecules else 0

    valid2 = [(name, mol) for name, mol in path2_molecules if is_valid_molecule(mol)]
    invalid2_count = len(path2_molecules) - len(valid2)

    if path1_molecules:
        print(f"PATH1: {len(valid1)} valid molecules, {invalid1_count} invalid molecules")
    print(f"PATH2: {len(valid2)} valid molecules, {invalid2_count} invalid molecules")

    #  Require at least one valid molecule overall
    if not valid1 and not valid2:
        print('\nError: No valid molecules in PATH1 or PATH2.')
        return

    if valid1 and valid2:
        print(
            f"\nProcessing {len(valid1)} molecules from PATH1 and "
            f"{len(valid2)} molecules from PATH2..."
        )
    elif valid1:
        print(f"\nProcessing {len(valid1)} molecules from PATH1...")
    else:
        print(f"\nProcessing {len(valid2)} molecules from PATH2...")

    #  === STEP 4: Compute Properties and Build CSV Data ===
    csv_data1 = []
    if valid1:
        for name, mol in valid1:
            props = compute_properties(mol)
            if props is not None:
                #  Compute Tanimoto similarity vs. the reference set
                max_tanimoto = compute_tanimoto_similarity(mol, reference_data)
                #  Build CSV row
                csv_row = {
                    'smiles': get_smiles(mol),
                    'valid': 1,
                    'MW': round(props['MolecularWeight'], 2),
                    'LogP': round(props['LogP'], 2),
                    'HBA': props['NumHAcceptors'],
                    'HBD': props['NumHDonors'],
                    'TPSA': round(props['TPSA'], 2),
                    'RotBonds': props['RotBonds'],
                    'QED': round(props['QED'], 4),
                    'SA_score': (
                        round(props['SA_score'], 2)
                        if props['SA_score'] is not None
                        else None
                    ),
                    'Lipinski': props['LipinskiRuleOfFive'],
                    'max_tanimoto': round(max_tanimoto, 4) if max_tanimoto is not None else None,
                }
                csv_data1.append(csv_row)
        print(f"Computed properties for {len(csv_data1)} molecules in PATH1")

    csv_data2 = []
    for name, mol in valid2:
        props = compute_properties(mol)
        if props is not None:
            #  Compute Tanimoto similarity vs. the reference set
            max_tanimoto = compute_tanimoto_similarity(mol, reference_data)
            #  Build CSV row
            csv_row = {
                'smiles': get_smiles(mol),
                'valid': 1,
                'MW': round(props['MolecularWeight'], 2),
                'LogP': round(props['LogP'], 2),
                'HBA': props['NumHAcceptors'],
                'HBD': props['NumHDonors'],
                'TPSA': round(props['TPSA'], 2),
                'RotBonds': props['RotBonds'],
                'QED': round(props['QED'], 4),
                'SA_score': round(props['SA_score'], 2) if props['SA_score'] is not None else None,
                'Lipinski': props['LipinskiRuleOfFive'],
                'max_tanimoto': round(max_tanimoto, 4) if max_tanimoto is not None else None,
            }
            csv_data2.append(csv_row)
    print(f"Computed properties for {len(csv_data2)} molecules in PATH2")

    #  === STEP 5: Save Results to CSV ===
    csv_file1 = Path(args.csv1) if args.csv1 else output_dir / 'pocket2mol_generated.csv'
    csv_file2 = Path(args.csv2) if args.csv2 else output_dir / 'diffsbdd_generated.csv'
    if not args.plots_only:
        #  Only save CSV if the corresponding input path was provided or used as default
        if path1:
            save_properties_to_csv(csv_data1, str(csv_file1))
        if path2:
            save_properties_to_csv(csv_data2, str(csv_file2), append=args.append)
    else:
        print('Plots-only mode: skipping CSV writes')

    pocket2mol_matrix = None
    pocket2mol_labels = []
    if valid1:
        pocket2mol_matrix, pocket2mol_labels, _ = compute_similarity_matrix(
            valid1, reference_data
        )
        # Save full numeric matrix for Pocket2Mol
        if pocket2mol_matrix is not None:
            save_similarity_matrix_csv(
                pocket2mol_matrix,
                pocket2mol_labels,
                reference_data.get('names', []),
                output_dir / 'pocket2mol_similarity_matrix.csv',
            )

    diffsbdd_matrix = None
    diffsbdd_labels = []
    if valid2:
        diffsbdd_matrix, diffsbdd_labels, _ = compute_similarity_matrix(
            valid2, reference_data
        )
        # Save full numeric matrix for DiffSBDD
        if diffsbdd_matrix is not None:
            save_similarity_matrix_csv(
                diffsbdd_matrix,
                diffsbdd_labels,
                reference_data.get('names', []),
                output_dir / 'diffsbdd_similarity_matrix.csv',
            )

    plot_similarity_comparison(
        csv_file1,
        csv_file2,
        pocket2mol_matrix=pocket2mol_matrix,
        diffsbdd_matrix=diffsbdd_matrix,
        pocket2mol_labels=pocket2mol_labels,
        diffsbdd_labels=diffsbdd_labels,
        reference_labels=reference_data.get('names', []),
        output_dir=output_dir,
    )

    #  === STEP 6: Summary ===
    print("\nProcessing complete.")

    def print_summary_stats(name, data_list):
        if not data_list:
            return
        # Extract max_tanimoto values. Handles both dicts (from csv_data) and floats (from load_similarity_rows)
        if isinstance(data_list[0], dict):
            values = [row['max_tanimoto'] for row in data_list if row.get('max_tanimoto') is not None]
        else:
            values = data_list

        if not values:
            return

        print(f"{name} Similarity Summary:")
        print(f"  Mean max Tanimoto: {np.mean(values):.4f}")
        top_k = sorted(values, reverse=True)[:5]
        print(f"  Top 5 max Tanimoto: {', '.join([f'{v:.4f}' for v in top_k])}")

    if args.plots_only:
        print(f'Plots-only mode: reading similarity data from {csv_file1} and {csv_file2}')
        print_summary_stats("Pocket2Mol", load_similarity_rows(csv_file1))
        print_summary_stats("DiffSBDD", load_similarity_rows(csv_file2))
    else:
        if csv_data1:
            print(f"PATH1: {len(csv_data1)} molecules saved to {csv_file1}")
            print_summary_stats("Pocket2Mol", csv_data1)
        if csv_data2:
            print(f"PATH2: {len(csv_data2)} molecules saved to {csv_file2}")
            print_summary_stats("DiffSBDD", csv_data2)


if __name__ == "__main__":
    main()
