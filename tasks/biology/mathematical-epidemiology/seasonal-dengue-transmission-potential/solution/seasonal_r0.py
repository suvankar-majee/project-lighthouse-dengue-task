import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize

DATA = sys.argv[1] if len(sys.argv) > 1 else "/app/data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "/app/output"

FORMS = {"a": "briere", "EFD": "briere", "pEA": "quadratic", "MDR": "briere",
         "lf": "quadratic", "b": "briere", "c": "briere", "PDR": "briere"}
PROB = {"pEA", "b", "c"}
SIGMA_H, GAMMA, PERIOD = 1 / 5.9, 1 / 5.0, 365.0
EFFICIENCY = 0.55
STATION_MOVE, SITE_OFFSET = pd.Timestamp("2018-03-15"), 1.1
KH_OFFSET = 1.5


# --------------------------------------------------------------------------- traits
def tpc(form, T, q, T0, Tm, prob=False):
    T = np.asarray(T, float)
    inside = (T > T0) & (T < Tm)
    if form == "briere":
        v = q * T * (T - T0) * np.sqrt(np.clip(Tm - T, 0.0, None))
    else:
        v = q * (T - T0) * (Tm - T)
    v = np.where(inside, v, 0.0)
    return np.minimum(v, 1.0) if prob else v


def load_traits(path):
    df = pd.read_csv(path)
    df = df[df["notes"].fillna("").str.contains("failure|flagged", case=False) == False].copy()
    df["T"] = np.where(df["temp_unit"] == "F", (df["temperature"] - 32) * 5 / 9, df["temperature"])
    df.loc[df["lab_id"] == "KH-ENT", "T"] += KH_OFFSET
    v = df["value"].astype(float)
    v = np.where(df["unit"] == "percent", v / 100.0, v)
    df["v"] = v
    # constant daily mortality -> lifespan; EIP (days) -> PDR (per day)
    inv = df["trait"].isin(["mu", "EIP"])
    df.loc[inv, "v"] = 1.0 / df.loc[inv, "v"]
    df["trait"] = df["trait"].replace({"mu": "lf", "EIP": "PDR"})
    df = df[np.isfinite(df["v"])]
    return df[["lab_id", "trait", "T", "v"]]


def fit_trait(name, T, y):
    form, prob = FORMS[name], name in PROB
    pos = T[y > 0]

    def sse_profile(T0, Tm):
        # q enters linearly (ignoring the cap at 1, which is checked afterwards)
        g = tpc(form, T, 1.0, T0, Tm)
        if not np.any(g > 0):
            return np.inf, 0.0
        q = float(np.dot(g, y) / np.dot(g, g))
        r = y - tpc(form, T, q, T0, Tm, prob)
        return float(np.dot(r, r)), q

    best = (np.inf, None)
    for T0 in np.arange(pos.min() - 12, pos.min() - 0.05, 0.25):
        for Tm in np.arange(pos.max() + 0.05, pos.max() + 15, 0.25):
            s, q = sse_profile(T0, Tm)
            if s < best[0]:
                best = (s, (q, T0, Tm))
    q0, T00, Tm0 = best[1]

    def obj(p):
        q, T0, Tm = p
        if q <= 0 or T0 >= pos.min() or Tm <= pos.max():
            return 1e12
        r = y - tpc(form, T, q, T0, Tm, prob)
        return float(np.dot(r, r))

    res = minimize(obj, [q0, T00, Tm0], method="Nelder-Mead",
                   options={"xatol": 1e-10, "fatol": 1e-14, "maxiter": 20000, "maxfev": 40000})
    return tuple(res.x) if res.fun <= best[0] else (q0, T00, Tm0)


# --------------------------------------------------------------------------- model
class Model:
    def __init__(self, pars):
        self.p = pars  # name -> (q, T0, Tm)

    def tr(self, name, T):
        q, T0, Tm = self.p[name]
        return tpc(FORMS[name], T, q, T0, Tm, name in PROB)

    def m_eq(self, T):
        T = np.asarray(T, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            mu = 1.0 / self.tr("lf", T)
            mdr = self.tr("MDR", T)
            ell = mdr * (0.5 * self.tr("EFD", T) / mu - 1.0 / self.tr("pEA", T))
            m = 0.5 * mdr * ell / mu
        return np.where(np.isfinite(m) & (m > 0), m, 0.0)

    def r0_const(self, T, k):
        T = np.asarray(T, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            a, b, c, pdr = (self.tr(n, T) for n in ("a", "b", "c", "PDR"))
            mu = 1.0 / self.tr("lf", T)
            r2 = (a * b / mu) * (pdr / (pdr + mu)) * (a * c * k * self.m_eq(T) / GAMMA)
        return np.sqrt(np.where(np.isfinite(r2) & (r2 > 0), r2, 0.0))


def harmonic_fit(t, y, nh=3):
    w = 2 * np.pi / PERIOD
    cols = [np.ones_like(t)]
    for h in range(1, nh + 1):
        cols += [np.cos(h * w * t), np.sin(h * w * t)]
    X = np.column_stack(cols)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)

    def T(tt):
        tt = np.asarray(tt, float)
        out = beta[0] + 0 * tt
        for h in range(1, nh + 1):
            out = out + beta[2 * h - 1] * np.cos(h * w * tt) + beta[2 * h] * np.sin(h * w * tt)
        return out
    return T, beta[0]


def load_climate(path):
    d = pd.read_csv(path)
    d["date"] = pd.to_datetime(d["date"])
    d["T"] = pd.to_numeric(d["tmean_c"], errors="coerce")
    d.loc[d["T"] <= -90, "T"] = np.nan
    # stuck probe before the 2016-09-21 replacement: run of identical readings
    run_id = (d["T"] != d["T"].shift()).cumsum()
    run_len = d.groupby(run_id)["T"].transform("size")
    d.loc[(run_len >= 5) & d["T"].notna(), "T"] = np.nan
    d.loc[d["date"] < STATION_MOVE, "T"] += SITE_OFFSET
    d = d.dropna(subset=["T"])
    ndays = np.where(d["date"].dt.is_leap_year, 366.0, 365.0)
    t = (d["date"].dt.dayofyear.values - 0.5) * PERIOD / ndays
    return harmonic_fit(t, d["T"].values)


def periodic_vector(model, Tfun):
    def rhs(t, y):
        T = Tfun(t)
        ell, m = y
        return [model.tr("EFD", T) * m - model.tr("MDR", T) / model.tr("pEA", T) * ell - ell * ell,
                0.5 * model.tr("MDR", T) * ell - m / model.tr("lf", T)]
    y = np.array([1.0, 1.0])
    for _ in range(500):
        s = solve_ivp(rhs, (0, PERIOD), y, method="LSODA", rtol=1e-10, atol=1e-13)
        ynew = s.y[:, -1]
        done = np.max(np.abs(ynew - y) / np.abs(ynew)) < 1e-9
        y = ynew
        if done:
            break
    else:
        raise RuntimeError("periodic vector solution did not converge")
    return solve_ivp(rhs, (0, PERIOD), y, method="LSODA", rtol=1e-11, atol=1e-14,
                     dense_output=True).sol


def calibrate_scale(path, msol):
    s = pd.read_csv(path)
    s["week_start"] = pd.to_datetime(s["week_start"])
    mbar = []
    for w in s["week_start"]:
        t0 = w.dayofyear - 1.0
        tt = np.linspace(t0, t0 + 7.0, 141) % PERIOD
        mbar.append(np.mean(msol(tt)[1]))
    mbar = np.array(mbar)
    # Poisson MLE for C_i ~ Pois(eff * residents_i * k * mbar_i)
    return s["females_collected"].sum() / np.sum(EFFICIENCY * s["residents"].values * mbar)


def periodic_r0(model, Tfun, msol, k):
    def log_rho(lam):
        def rhs(t, x):
            T = Tfun(t)
            a, b, c, pdr = (model.tr(n, T) for n in ("a", "b", "c", "PDR"))
            mu = 1.0 / model.tr("lf", T)
            M = k * msol(t)[1]
            A = np.array([[-SIGMA_H, 0, 0, a * b / lam],
                          [SIGMA_H, -GAMMA, 0, 0],
                          [0, a * c * M / lam, -(pdr + mu), 0],
                          [0, 0, pdr, -mu]])
            return (A @ x.reshape(4, 4)).ravel()
        X = np.eye(4)
        logscale = 0.0
        for seg in range(12):  # renormalise each month to avoid overflow
            s = solve_ivp(rhs, (seg * PERIOD / 12, (seg + 1) * PERIOD / 12), X.ravel(),
                          method="LSODA", rtol=1e-10, atol=1e-12)
            X = s.y[:, -1].reshape(4, 4)
            nrm = np.linalg.norm(X)
            X /= nrm
            logscale += np.log(nrm)
        return np.log(np.max(np.abs(np.linalg.eigvals(X)))) + logscale

    lo, hi = 0.05, 1.0
    while log_rho(hi) > 0:
        hi *= 2
    return brentq(log_rho, lo, hi, xtol=1e-10)


def main():
    os.makedirs(OUT, exist_ok=True)
    tr = load_traits(os.path.join(DATA, "trait_experiments.csv"))
    pars, rows = {}, []
    for name in FORMS:
        sub = tr[tr["trait"] == name]
        q, T0, Tm = fit_trait(name, sub["T"].values, sub["v"].values)
        pars[name] = (q, T0, Tm)
        rows.append({"trait": name, "q": q, "T0": T0, "Tm": Tm})
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "trait_fits.csv"), index=False)
    model = Model(pars)

    Tfun, Tmean = load_climate(os.path.join(DATA, "station_daily_tmean.csv"))
    msol = periodic_vector(model, Tfun)
    k = calibrate_scale(os.path.join(DATA, "aspirator_survey_2021.csv"), msol)

    grid = np.round(np.arange(16.0, 34.0001, 0.5), 1)
    pd.DataFrame({"temperature_c": grid,
                  "females_per_person": k * model.m_eq(grid),
                  "r0": model.r0_const(grid, k)}).to_csv(
        os.path.join(OUT, "r0_constant_temperature.csv"), index=False)

    tt = np.linspace(0, PERIOD, 36501)
    dens = k * msol(tt)[1]
    peak_day = int(np.floor(tt[np.argmax(dens)])) + 1
    res = {
        "mean_temperature_c": float(Tmean),
        "mean_females_per_person": float(np.trapezoid(dens, tt) / PERIOD),
        "peak_density_day_of_year": peak_day,
        "r0_at_mean_temperature": float(model.r0_const(Tmean, k)),
        "r0_seasonal": float(periodic_r0(model, Tfun, msol, k)),
    }
    with open(os.path.join(OUT, "seasonal_summary.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
