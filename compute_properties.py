
"""Compute properties for molecules in PATH1 and PATH2, save them to a file, and compare them.

The following properties are computed:
- Molecular Weight
- LogP
- HBA
- HBD
- TPSA
- QED
- LIPINSKI'S RULE of FIVE
- Tanimoto similarity for PATH1 AND PATH2 vs. known A2A drugs in ADORA2A-world.sdf
"""

import json
from pathlib import Path
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdFingerprintGenerator, rdMolDescriptors


MORGAN_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def is_valid_molecule(mol):
    """Check if molecule can be processed without errors."""
    try:
        properties = {}
        properties['MolecularWeight'] = Descriptors.MolWt(mol)
        properties['LogP'] = Crippen.MolLogP(mol)
        properties['NumHAcceptors'] = Descriptors.NumHAcceptors(mol)
        properties['NumHDonors'] = Descriptors.NumHDonors(mol)
        properties['TPSA'] = rdMolDescriptors.CalcTPSA(mol)
        properties['QED'] = QED.qed(mol)
        properties['LipinskiRuleOfFive'] = int(
            Descriptors.MolWt(mol) < 500,
        ) + int(Lipinski.NumHDonors(mol) <= 5) + int(Lipinski.NumHAcceptors(mol) <= 10) + int(
            Crippen.MolLogP(mol) <= 5
        ) + int(rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10)
        return True
    except Exception:
        return False


def compute_properties(mol):
    """Compute molecular properties. Returns None if computation fails."""
    try:
        properties = {}
        properties['MolecularWeight'] = Descriptors.MolWt(mol)
        properties['LogP'] = Crippen.MolLogP(mol)
        properties['NumHAcceptors'] = Descriptors.NumHAcceptors(mol)
        properties['NumHDonors'] = Descriptors.NumHDonors(mol)
        properties['TPSA'] = rdMolDescriptors.CalcTPSA(mol)
        properties['QED'] = QED.qed(mol)
        properties['LipinskiRuleOfFive'] = int(
            Descriptors.MolWt(mol) < 500,
        ) + int(Lipinski.NumHDonors(mol) <= 5) + int(Lipinski.NumHAcceptors(mol) <= 10) + int(
            Crippen.MolLogP(mol) <= 5
        ) + int(rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10)
        return properties
    except Exception as e:
        return None


def compute_tanimoto_similarity(mol, reference_fingerprints):
    """Return maximum Tanimoto similarity to the known A2A active set."""
    try:
        if not reference_fingerprints:
            return None
        fp = get_morgan_fingerprint(mol)
        if fp is None:
            return None
        similarities = [DataStructs.TanimotoSimilarity(fp, ref_fp)
                        for ref_fp in reference_fingerprints]
        return max(similarities) if similarities else None
    except Exception:
        return None


def get_morgan_fingerprint(mol):
    """Build a Morgan fingerprint from a sanitized molecule."""
    try:
        mol_copy = Chem.Mol(mol)
        Chem.SanitizeMol(mol_copy)
        return MORGAN_GENERATOR.GetFingerprint(mol_copy)
    except Exception:
        return None


def load_a2a_reference_fingerprints(a2a_sdf_file):
    """Load known A2A actives from ADORA2A-world.sdf and build fingerprints."""
    sdf_path = Path(a2a_sdf_file)
    if not sdf_path.is_file():
        raise ValueError(f'A2A SDF file not found: {a2a_sdf_file}')

    fingerprints = []
    suppl = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
    for mol in suppl:
        if mol is None or not is_valid_molecule(mol):
            continue
        fp = get_morgan_fingerprint(mol)
        if fp is not None:
            fingerprints.append(fp)
    return fingerprints


def load_molecules_from_path1(path1):
    """Load all SDF files from a directory for PATH1.

    Each file is treated as a single-molecule SDF and only the first valid
    molecule in each file is used.
    """
    path1 = Path(path1)
    if not path1.is_dir():
        raise ValueError(f'PATH1 must be a directory containing SDF files: {path1}')

    molecules = []
    for sdf_file in sorted(path1.glob('*.sdf')):
        suppl = Chem.SDMolSupplier(str(sdf_file), sanitize=False, removeHs=False)
        mol = next((m for m in suppl if m is not None), None)
        if mol is not None:
            molecules.append((sdf_file.name, mol))
    return molecules


def load_molecules_from_path2(path2):
    """Load all molecules from one multi-molecule SDF for PATH2."""
    path = Path(path2)
    if not path.is_file():
        raise ValueError(f'PATH2 must be a single SDF file: {path2}')

    suppl = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
    molecules = []
    for idx, mol in enumerate(suppl):
        if mol is not None:
            molecules.append((f'{path.name}#{idx}', mol))
    return molecules


def save_properties(properties, output_file):
    with open(output_file, 'w') as f:
        json.dump(properties, f, indent=4)


def load_properties(input_file):
    with open(input_file, 'r') as f:
        properties = json.load(f)
    return properties


def compare_properties(properties1, properties2):
    differences = {}
    for key in properties1:
        if key in properties2:
            if properties1[key] != properties2[key]:
                differences[key] = (properties1[key], properties2[key])
    return differences


def compare_property_lists(properties1, properties2):
    differences = []
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


def main():
    # Hardcoded paths
    path1 = '/Users/ccc/Library/Mobile Documents/com~apple~CloudDocs/_FHNW/4. Semester/AI in Drug Discovery/project/ai-in-drug-discovery-group2/outputs2/sample_for_pdb_4EIY_clean.pdb_2026_05_02__11_02_39/SDF'
    path2 = '/Users/ccc/Library/Mobile Documents/com~apple~CloudDocs/_FHNW/4. Semester/AI in Drug Discovery/project/ai-in-drug-discovery-group2/outputs2/4eiy_inpaint.sdf'
    a2a_sdf_file = '/Users/ccc/Library/Mobile Documents/com~apple~CloudDocs/_FHNW/4. Semester/AI in Drug Discovery/project/ai-in-drug-discovery-group2/outputs2/a2a_drugs/ADORA2A-world.sdf'
    output_file = '/Users/ccc/Library/Mobile Documents/com~apple~CloudDocs/_FHNW/4. Semester/AI in Drug Discovery/project/ai-in-drug-discovery-group2/outputs2/compute_properties.json'

    # Load known A2A active references
    reference_fingerprints = load_a2a_reference_fingerprints(a2a_sdf_file)
    if not reference_fingerprints:
        print('Error: Could not load any valid known A2A actives from ADORA2A-world.sdf.')
        return
    print(f'Loaded {len(reference_fingerprints)} known A2A active fingerprints')

    # Load molecules
    path1_molecules = load_molecules_from_path1(path1)
    path2_molecules = load_molecules_from_path2(path2)

    if not path1_molecules or not path2_molecules:
        print('Error: Could not load one or more molecules from PATH1 or PATH2.')
        return
    
    print(f"Loaded {len(path1_molecules)} molecules from PATH1")
    print(f"Loaded {len(path2_molecules)} molecules from PATH2")
    
    # Validate molecules before processing
    print("\nValidating molecules...")
    valid1 = [(name, mol) for name, mol in path1_molecules if is_valid_molecule(mol)]
    invalid1_count = len(path1_molecules) - len(valid1)
    
    valid2 = [(name, mol) for name, mol in path2_molecules if is_valid_molecule(mol)]
    invalid2_count = len(path2_molecules) - len(valid2)
    
    print(f"PATH1: {len(valid1)} valid molecules, {invalid1_count} invalid molecules")
    print(f"PATH2: {len(valid2)} valid molecules, {invalid2_count} invalid molecules")
    
    if not valid1 or not valid2:
        print('\nError: No valid molecules to compare in one or both paths.')
        return
    
    print(f"\nProcessing {len(valid1)} molecules from PATH1 and {len(valid2)} molecules from PATH2...")

    # Compute properties
    properties1 = []
    for name, mol in valid1:
        props = compute_properties(mol)
        if props is not None:
            props['TanimotoSimilarityA2A'] = compute_tanimoto_similarity(
                mol, reference_fingerprints)
            properties1.append({'name': name, 'properties': props})
    
    properties2 = []
    for name, mol in valid2:
        props = compute_properties(mol)
        if props is not None:
            props['TanimotoSimilarityA2A'] = compute_tanimoto_similarity(
                mol, reference_fingerprints)
            properties2.append({'name': name, 'properties': props})

    print(f"\nComputed properties for {len(properties1)} molecules in PATH1")
    print(f"Computed properties for {len(properties2)} molecules in PATH2")
    
    # Save properties to file
    save_properties({'Path1': properties1, 'Path2': properties2}, output_file)
    print(f"\nProperties saved to {output_file}")

    # Compare properties
    differences = compare_property_lists(properties1, properties2)
    if differences:
        print("Differences found in properties:")
        for entry in differences:
            print(f"Molecule pair {entry['index']}: {entry['path1_name']} vs {entry['path2_name']}")
            for key, (val1, val2) in entry['differences'].items():
                print(f"  {key}: Path1={val1}, Path2={val2}")
    else:
        print("No differences found in properties.")


if __name__ == "__main__":
    main()