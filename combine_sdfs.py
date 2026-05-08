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
        input_path = sys.argv[1]
        output_file = sys.argv[2]
    elif len(sys.argv) == 2:
        input_path = sys.argv[1]
        p = Path(input_path)
        if p.is_dir():
            output_file = p / "combined.sdf"
        else:
            output_file = p.with_name(p.stem + "_combined.sdf")
    else:
        print("Usage: python combine_sdfs.py <input_dir_or_sdf> [output_file]")
        sys.exit(1)

    success = combine_sdfs(input_path, output_file)
    if success:
        print(f"Combined SDF written to: {output_file}")
        sys.exit(0)
    else:
        print("No molecules written; check input path and SDF files.")
        sys.exit(2)