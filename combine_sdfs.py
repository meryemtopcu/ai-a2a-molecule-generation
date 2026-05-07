import sys
from pathlib import Path
from rdkit import Chem

def combine_sdfs(input_path, output_file):
    """
    Combines multiple SDF files from a directory into a single multi-molecule SDF file.
    
    Args:
        input_path: Path to a directory containing SDF files or a single SDF file.
        output_file: Path where the combined SDF file will be saved.
    
    Returns:
        bool: True if successful (at least one molecule written), False otherwise.
    """
    path = Path(input_path)
    if path.is_file():
        files = [path]
    elif path.is_dir():
        # Find any file containing '.sdf' (case-insensitive) to catch .sdf_copy etc.
        files = sorted([f for f in path.iterdir() if f.is_file() and '.sdf' in f.name.lower()])
    else:
        return False

    if not files:
        return False

    writer = Chem.SDWriter(str(output_file))
    count = 0
    for sdf_file in files:
        suppl = Chem.SDMolSupplier(str(sdf_file), sanitize=False, removeHs=False)
        for mol in suppl:
            if mol is not None:
                writer.write(mol)
                count += 1
    writer.close()
    
    return count > 0

if __name__ == "__main__":
    if len(sys.argv) == 3:
        combine_sdfs(sys.argv[1], sys.argv[2])