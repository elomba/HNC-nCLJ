#!/usr/bin/env python3
"""Synthetic checks for transport estimators without a LAMMPS run."""
import numpy as np

from analyze_followup import block_diffusion, fit_msd, msd_curve, viscosity, yeh_hummer


def test_msd_and_blocks():
    rng = np.random.default_rng(123)
    diffusion = 0.2
    time = np.arange(2001, dtype=float)
    increments = rng.normal(0, np.sqrt(2 * diffusion), size=(2000, 80, 3))
    position = np.concatenate([np.zeros((1, 80, 3)), np.cumsum(increments, axis=0)]).astype(
        np.float32
    )
    position -= position.mean(axis=1, keepdims=True)
    typ = np.r_[np.ones(40, dtype=int), np.full(40, 2, dtype=int)]
    lag, curve = msd_curve(time, position, typ, 2, origins=150)
    result = fit_msd(lag, curve)
    assert abs(result["20_100"][0] - diffusion * (1 - 1 / 80)) < 0.05
    assert np.allclose(curve[:, 0], curve[:, 1:].mean(axis=1), atol=1e-5)
    blocks = block_diffusion(time, position, typ, ["a", "b"])
    assert len(blocks["blocks"]) == 5
    assert np.isfinite(blocks["summary"]["net"]["standard_error"])


def test_viscosity_and_yeh_hummer():
    rng = np.random.default_rng(52)
    dt = 0.02
    tau = 0.5
    rho = np.exp(-dt / tau)
    stress = np.empty((60000, 5))
    stress[0] = rng.normal(size=5)
    for i in range(1, len(stress)):
        stress[i] = rho * stress[i - 1] + np.sqrt(1 - rho**2) * rng.normal(size=5)
    eta = viscosity(stress, dt, volume=2.0, temperature=2.0)
    assert abs(eta["cutoffs"]["2"] - tau) < 0.12
    assert "plateau" in eta
    corrected = yeh_hummer(0.2, 1.5, 4000)
    expected = 0.2 + 2.837297 * 2 / (6 * np.pi * 1.5 * (4000 / 0.7) ** (1 / 3))
    assert np.isclose(corrected, expected)


if __name__ == "__main__":
    test_msd_and_blocks()
    test_viscosity_and_yeh_hummer()
    print("Synthetic MSD, block, viscosity, and Yeh–Hummer tests passed")
