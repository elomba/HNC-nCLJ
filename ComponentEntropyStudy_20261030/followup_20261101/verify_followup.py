#!/usr/bin/env python3
"""Verify completion and core transport records of one dense replica."""
import argparse
import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def verify(folder):
    manifest = json.loads((folder / "manifest.json").read_text())
    log = (folder / "lammps_run.log").read_text(errors="replace")
    assert "Total wall time:" in log and "ERROR" not in log, folder
    assert (folder / "data.final").is_file() and (folder / "final.restart").is_file(), folder
    msd = np.atleast_2d(np.loadtxt(folder / "msd_common.dat"))
    fractions = np.array([item["x"] for item in manifest["composition"]["species"]])
    assert msd.shape[1] == len(fractions) + 2
    assert np.allclose(msd[:, 1], msd[:, 2:] @ fractions, rtol=5e-5, atol=2e-5)
    with Dataset(folder / "positions.nc") as dataset:
        n_frames = len(dataset.variables["time"])
    expected = manifest["production_steps"] // 500 + 1
    assert abs(n_frames - expected) <= 1, (folder, n_frames, expected)
    result = {"normal_completion": True, "position_frames": n_frames, "expected_frames": expected}
    (folder / "verified.json").write_text(json.dumps(result, indent=2) + "\n")
    print("Verified", folder, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("replica", type=Path)
    verify(parser.parse_args().replica)
