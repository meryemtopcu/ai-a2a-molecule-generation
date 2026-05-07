#!/usr/bin/env python3
"""Small utility to apply trivial pylint-friendly fixes to Python files.

Current actions:
- Trim trailing whitespace on all lines
- Wrap long comment lines (starting with `#`) to `max_line_length`

Usage: python scripts/pylint_fix.py [paths...]
If no paths provided, runs on current directory recursively.
"""
from __future__ import annotations

import argparse
import os
import sys
import textwrap
from pathlib import Path


def fix_file(path: Path, max_line_length: int) -> bool:
    """Return True if file changed."""
    try:
        text = path.read_text(encoding="utf8")
    except Exception:
        return False

    lines = text.splitlines(keepends=True)
    changed = False
    out_lines = []
    in_triple = False
    triple_delims = ("'''", '"""')

    for line in lines:
        original = line
        # Detect entering/exiting triple-quoted blocks conservatively
        stripped = line.lstrip()
        if any(delim in stripped for delim in triple_delims):
            # toggle heuristically: if delim present, flip state
            # This is conservative and avoids touching docstrings.
            in_triple = not in_triple

        if in_triple:
            # Don't modify content inside triple-quoted strings
            out_lines.append(line)
            continue

        # Trim trailing whitespace but preserve newline
        if line.endswith("\n"):
            core = line[:-1].rstrip()
            new_line = core + "\n"
        else:
            new_line = line.rstrip()

        # Wrap long comment lines only (safe)
        if new_line.lstrip().startswith("#"):
            indent = new_line[: len(new_line) - len(new_line.lstrip())]
            comment = new_line.lstrip()
            # Remove leading '#', preserve single space after it
            if comment.startswith("# "):
                prefix = "# "
                content = comment[2:]
            else:
                prefix = "#"
                content = comment[1:]

            wrapper = textwrap.TextWrapper(width=max_line_length,
                                           subsequent_indent=indent + prefix,
                                           break_long_words=False,
                                           replace_whitespace=False)
            wrapped = wrapper.wrap(content)
            if wrapped:
                wrapped_lines = [indent + prefix + ("" if p == "" else " ") + p if not p.startswith(" ") else indent + prefix + p for p in wrapped]
                # ensure each has newline
                wrapped_lines = [ln + "\n" for ln in wrapped_lines]
                new_line = "".join(wrapped_lines)

        if new_line != original:
            changed = True

        out_lines.append(new_line)

    if changed:
        path.write_text("".join(out_lines), encoding="utf8")

    return changed


def iter_py_files(paths: list[str]) -> list[Path]:
    files: list[Path] = []
    for p in paths or ["."]:
        pth = Path(p)
        if pth.is_file() and pth.suffix == ".py":
            files.append(pth)
        elif pth.is_dir():
            for f in pth.rglob("*.py"):
                files.append(f)
    return files


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", help="Files or dirs to process")
    parser.add_argument("--max-line-length", type=int, default=100)
    args = parser.parse_args(argv)

    files = iter_py_files(args.paths)
    if not files:
        print("No Python files found.")
        return 0

    changed_files = []
    for f in files:
        try:
            if fix_file(f, args.max_line_length):
                changed_files.append(str(f))
        except Exception as exc:
            print(f"Failed to process {f}: {exc}", file=sys.stderr)

    if changed_files:
        print("Modified files:")
        for c in changed_files:
            print(" -", c)
    else:
        print("No changes needed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
