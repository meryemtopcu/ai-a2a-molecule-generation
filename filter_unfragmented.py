#!/usr/bin/env python3
"""Filter SDF(s) to keep only single-fragment molecules.

Usage:
  python filter_unfragmented.py <input_sdf_or_dir> <output_sdf> [--min-heavy-atoms N]

If <input_sdf_or_dir> is a directory, all files with '.sdf' in the name
(case-insensitive) will be processed.

The script writes a multi-molecule SDF containing only molecules whose
RDKit fragment count is 1. Optionally filters by minimum heavy atom count
or salvages the largest fragment if --salvage is provided.
"""
import sys
from pathlib import Path
from rdkit import Chem


def process_files(inputs, out_path, min_heavy_atoms=0, salvage=False):
    writer = Chem.SDWriter(str(out_path))
    total_in = 0
    written = 0
    skipped_fragment = 0
    skipped_none = 0
    skipped_small = 0

    for inp in inputs:
        suppl = Chem.SDMolSupplier(str(inp), sanitize=False, removeHs=False)
        for mol in suppl:
            total_in += 1
            if mol is None:
                skipped_none += 1
                continue
            try:
                frags = Chem.GetMolFrags(mol, asMols=True)
            except Exception:
                # If RDKit can't compute fragments for this mol, skip it
                skipped_none += 1
                continue
            
            if len(frags) > 1:
                if not salvage:
                    skipped_fragment += 1
                    continue
                # Repair by keeping the fragment with the most heavy atoms
                mol = max(frags, key=lambda f: f.GetNumHeavyAtoms())
            else:
                mol = frags[0]

            if min_heavy_atoms and mol.GetNumHeavyAtoms() < min_heavy_atoms:
                skipped_small += 1
                continue
            writer.write(mol)
            written += 1
    writer.close()
    return {
        'total': total_in,
        'written': written,
        'skipped_fragment': skipped_fragment,
        'skipped_none': skipped_none,
        'skipped_small': skipped_small,
    }


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: python filter_unfragmented.py <input_sdf_or_dir> <output_sdf> [--min-heavy-atoms N] [--salvage]')
        sys.exit(1)
    input_arg = sys.argv[1]
    output_sdf = sys.argv[2]
    min_heavy = 0
    if '--min-heavy-atoms' in sys.argv:
        try:
            idx = sys.argv.index('--min-heavy-atoms')
            min_heavy = int(sys.argv[idx + 1])
        except Exception:
            print('Invalid value for --min-heavy-atoms')
            sys.exit(2)

    salvage = '--salvage' in sys.argv

    p = Path(input_arg)
    inputs = []
    if p.is_file():
        inputs = [p]
    elif p.is_dir():
        inputs = sorted([f for f in p.iterdir() if f.is_file() and '.sdf' in f.name.lower()])
    else:
        print(f'Input path not found: {input_arg}')
        sys.exit(3)

    out_path = Path(output_sdf)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    stats = process_files(inputs, out_path, min_heavy_atoms=min_heavy, salvage=salvage)
    print(f"Processed {stats['total']} input molecules")
    print(f"Written {stats['written']} molecules to: {out_path}")
    print(f"Skipped (None/invalid): {stats['skipped_none']}")
    print(f"Skipped (multi-fragment): {stats['skipped_fragment'] if not salvage else 0}")
    if min_heavy:
        print(f"Skipped (too small < {min_heavy} heavy atoms): {stats['skipped_small']}")
