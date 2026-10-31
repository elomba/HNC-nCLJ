#!/usr/bin/env python3
"""Check HNC numerical output against references from the current solver.

Run ``python3 tests/regression.py --update-reference`` only after an intentional,
explained numerical change. Historical files in baseline_results are never updated.
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "bin/hncr_nclj"
FILES = ("Thermohncr.dat", "RDFhncr.dat", "skhnc.dat")
CASES = {"3c": ROOT / "datosLJ.hncr"}


def numeric_rows(path):
    """Load numeric rows while ignoring blank and commented header lines."""
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.size == 0:
        raise ValueError(f"no numeric rows in {path}")
    return data


def compare(actual_path, reference_path, *, thermodynamics=False):
    actual = numeric_rows(actual_path)
    reference = numeric_rows(reference_path)
    label = f"{actual_path.name} vs {reference_path}"
    if actual.shape != reference.shape:
        print(f"  FAIL {label}: shape {actual.shape} != {reference.shape}")
        return False

    deviation = np.abs(actual - reference)
    row, column = np.unravel_index(np.argmax(deviation), deviation.shape)
    atol = 2e-4 if thermodynamics else 1e-8
    rtol = 0.0 if thermodynamics else 1e-6
    passed = bool(np.all(np.isfinite(actual)) and np.all(np.isfinite(reference)))
    passed = passed and bool(np.allclose(actual, reference, atol=atol, rtol=rtol))
    status = "PASS" if passed else "FAIL"
    print(
        f"  {status} {label}: max abs deviation {deviation[row, column]:.6g} "
        f"(row {row + 1}, column {column + 1})"
    )
    return passed


def run_case(case, input_path, reference_root, update):
    print(f"{case}:")
    with tempfile.TemporaryDirectory(prefix=f"hnc-{case}-") as workdir:
        result = subprocess.run(
            [str(BINARY), str(input_path)],
            cwd=workdir,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode:
            print(f"  FAIL solver exited with status {result.returncode}")
            if result.stderr:
                print(result.stderr[-2000:])
            return False

        output_dir = Path(workdir)
        reference_dir = reference_root / case
        try:
            if update:
                reference_dir.mkdir(parents=True, exist_ok=True)
                for filename in FILES:
                    destination = reference_dir / filename
                    shutil.copyfile(output_dir / filename, destination)
                    destination.chmod(0o644)
                print(f"  Updated {reference_dir}")
            results = [
                compare(
                    output_dir / filename,
                    reference_dir / filename,
                    thermodynamics=filename == "Thermohncr.dat",
                )
                for filename in FILES
            ]
            if case == "3c":
                results.extend(
                    compare(output_dir / filename, ROOT / "baseline_results" / filename)
                    for filename in FILES[1:]
                )
            passed = all(results)
        except (OSError, ValueError) as exc:
            print(f"  FAIL {exc}")
            return False
    print(f"{case}: {'PASS' if passed else 'FAIL'}")
    return passed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--update-reference",
        action="store_true",
        help="regenerate current references after an explained numerical change",
    )
    parser.add_argument(
        "--reference-root",
        type=Path,
        default=ROOT / "tests/reference",
        help="reference tree (useful for checking a temporary copy)",
    )
    args = parser.parse_args()
    if not BINARY.is_file():
        parser.error(f"missing solver binary {BINARY}; run `make` first")

    passed = all(
        run_case(case, input_path, args.reference_root, args.update_reference)
        for case, input_path in CASES.items()
    )
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
