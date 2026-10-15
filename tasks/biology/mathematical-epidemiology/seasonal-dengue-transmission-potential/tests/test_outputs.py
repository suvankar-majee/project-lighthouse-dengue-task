import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference as R  # noqa: E402

DATA = "/app/data"
OUT = Path("/app/output")
TRAITS = ["a", "EFD", "pEA", "MDR", "lf", "b", "c", "PDR"]
GRID = np.round(np.arange(16.0, 34.0001, 0.5), 1)
SUMMARY_KEYS = ["mean_temperature_c", "mean_females_per_person", "peak_density_day_of_year",
                "r0_at_mean_temperature", "r0_seasonal"]


@pytest.fixture(scope="session")
def ref():
    return R.full_reference(DATA)


@pytest.fixture(scope="session")
def fits():
    df = pd.read_csv(OUT / "trait_fits.csv")
    assert {"trait", "q", "T0", "Tm"} <= set(df.columns)
    df = df.set_index(df["trait"].astype(str).str.strip())
    assert set(TRAITS) <= set(df.index), "trait_fits.csv must contain all eight traits"
    pars = {t: tuple(float(df.loc[t, c]) for c in ("q", "T0", "Tm")) for t in TRAITS}
    for t, p in pars.items():
        assert all(math.isfinite(x) for x in p) and p[0] > 0 and p[1] < p[2], t
    return pars


@pytest.fixture(scope="session")
def curve():
    df = pd.read_csv(OUT / "r0_constant_temperature.csv")
    assert {"temperature_c", "females_per_person", "r0"} <= set(df.columns)
    df = df.sort_values("temperature_c").reset_index(drop=True)
    assert len(df) == len(GRID) and np.allclose(df["temperature_c"].values, GRID, atol=1e-6), \
        "temperature grid must be 16.0, 16.5, ..., 34.0 degC"
    assert np.all(np.isfinite(df[["females_per_person", "r0"]].values))
    assert np.all(df[["females_per_person", "r0"]].values >= 0)
    return df


@pytest.fixture(scope="session")
def summary():
    with open(OUT / "seasonal_summary.json") as f:
        s = json.load(f)
    for k in SUMMARY_KEYS:
        assert k in s, f"missing key {k}"
        assert isinstance(s[k], (int, float)) and math.isfinite(s[k]), k
    return s


@pytest.fixture(scope="session")
def agent_scale(fits, curve):
    """Density scale implied by the agent's own constant-temperature table."""
    m = R.m_eq(fits, GRID)
    sel = m > 0.02 * m.max()
    assert sel.sum() >= 10, "fitted model gives no viable mosquito population"
    ratios = curve["females_per_person"].values[sel] / m[sel]
    return float(np.median(ratios)), ratios


# ------------------------------------------------------------------ thermal traits
@pytest.mark.parametrize("trait", TRAITS)
def test_thermal_curve_matches_independent_fit(trait, fits, ref):
    T = np.arange(18.0, 34.01, 1.0)
    f_agent = R.tpc(trait, T, *fits[trait])
    f_ref = R.tpc(trait, T, *ref["pars"][trait])
    tol = 0.12 * f_ref + 0.03 * f_ref.max()
    worst = np.max(np.abs(f_agent - f_ref) - tol)
    assert worst <= 0, f"{trait}: fitted curve departs from the data-supported response"


# ------------------------------------------------------------------ constant temperature
def test_density_column_follows_fitted_vector_model(agent_scale):
    k, ratios = agent_scale
    assert k > 0
    assert np.max(np.abs(ratios / k - 1)) < 0.01, \
        "females_per_person is not the equilibrium of the specified vector model"


def test_r0_column_follows_fitted_transmission_model(fits, curve, agent_scale):
    k, _ = agent_scale
    expected = R.r0_const(fits, GRID, k)
    got = curve["r0"].values
    assert np.all(np.abs(got - expected) <= 0.015 * expected + 0.01), \
        "r0 column inconsistent with the fitted traits under the specified model"


def test_constant_temperature_r0_matches_reference(curve, ref):
    r_ref = R.r0_const(ref["pars"], GRID, ref["k"])
    sel = r_ref >= 1.0
    rel = np.abs(curve["r0"].values[sel] / r_ref[sel] - 1)
    assert rel.max() <= 0.06


def test_constant_temperature_density_matches_reference(curve, ref):
    d_ref = ref["k"] * R.m_eq(ref["pars"], GRID)
    err = np.abs(curve["females_per_person"].values - d_ref)
    assert np.all(err <= 0.06 * d_ref + 0.03 * d_ref.max())


# ------------------------------------------------------------------ seasonal analysis
def test_mean_temperature_of_homogenised_climate(summary, ref):
    assert abs(summary["mean_temperature_c"] - ref["Tmean"]) <= 0.1


def test_survey_calibrated_mean_density(summary, ref):
    assert abs(summary["mean_females_per_person"] / ref["mean_density"] - 1) <= 0.03


def test_peak_density_timing(summary, ref):
    d = abs(int(round(summary["peak_density_day_of_year"])) - ref["peak_day"]) % 365
    assert min(d, 365 - d) <= 6


def test_r0_at_mean_temperature(summary, ref):
    assert abs(summary["r0_at_mean_temperature"] / ref["r0_mean_T"] - 1) <= 0.04


def test_seasonal_r0_matches_reference(summary, ref):
    assert abs(summary["r0_seasonal"] / ref["r0_seasonal"] - 1) <= 0.03


def test_seasonal_r0_recomputed_from_agent_model(summary, fits, agent_scale, ref):
    """Periodic R0 recomputed with the agent's own traits and density scale."""
    k, _ = agent_scale
    th, m = R.periodic_vector(fits, ref["Tf"])
    r_agent_model = R.periodic_r0(fits, ref["Tf"], th, m, k)
    assert abs(summary["r0_seasonal"] / r_agent_model - 1) <= 0.025
