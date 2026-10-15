import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

FORMS = {"a": "briere", "EFD": "briere", "pEA": "quadratic", "MDR": "briere",
         "lf": "quadratic", "b": "briere", "c": "briere", "PDR": "briere"}
PROB = {"pEA", "b", "c"}
SIGMA_H, GAMMA, PERIOD = 1 / 5.9, 1 / 5.0, 365.0
STEPS_PER_DAY = 8


def tpc(name, T, q, T0, Tm, cap=True):
    T = np.asarray(T, float)
    ok = (T > T0) & (T < Tm)
    if FORMS[name] == "briere":
        v = q * T * (T - T0) * np.sqrt(np.clip(Tm - T, 0.0, None))
    else:
        v = q * (T - T0) * (Tm - T)
    v = np.where(ok, v, 0.0)
    return np.minimum(v, 1.0) if (cap and name in PROB) else v


def clean_traits(data):
    df = pd.read_csv(f"{data}/trait_experiments.csv")
    df = df[~df["notes"].fillna("").str.lower().str.contains("failure")]
    T = df["temperature"].astype(float).values.copy()
    f = df["temp_unit"].values == "F"
    T[f] = (T[f] - 32.0) / 1.8
    T[df["lab_id"].values == "KH-ENT"] += 1.5          # setpoint -> realised temperature
    v = df["value"].astype(float).values.copy()
    v[df["unit"].values == "percent"] /= 100.0
    name = df["trait"].values.copy()
    for raw, tgt in (("mu", "lf"), ("EIP", "PDR")):
        sel = name == raw
        v[sel] = 1.0 / v[sel]
        name[sel] = tgt
    keep = np.isfinite(v)
    return pd.DataFrame({"trait": name[keep], "T": T[keep], "v": v[keep]})


def fit_curves(data):
    df = clean_traits(data)
    out = {}
    for name in FORMS:
        s = df[df["trait"] == name]
        T, y = s["T"].values, s["v"].values
        lo, hi = T[y > 0].min(), T[y > 0].max()
        best = None
        for T0 in np.linspace(lo - 10, lo - 0.2, 9):
            for Tm in np.linspace(hi + 0.2, hi + 12, 9):
                g = tpc(name, T, 1.0, T0, Tm)
                q = max(np.dot(g, y) / max(np.dot(g, g), 1e-300), 1e-12)
                r = least_squares(lambda p: y - tpc(name, T, *p, cap=False), [q, T0, Tm],
                                  bounds=([0, lo - 25, hi + 1e-3], [np.inf, lo - 1e-3, hi + 40]),
                                  x_scale=[q, 1, 1], xtol=1e-12, ftol=1e-12, gtol=1e-12)
                if best is None or r.cost < best.cost:
                    best = r
        out[name] = tuple(best.x)
    return out


def traits_at(pars, T):
    return {n: tpc(n, T, *pars[n]) for n in FORMS}


def m_eq(pars, T):
    t = traits_at(pars, T)
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = 1 / t["lf"]
        m = 0.5 * t["MDR"] * t["MDR"] * (0.5 * t["EFD"] / mu - 1 / t["pEA"]) / mu
    return np.where(np.isfinite(m) & (m > 0), m, 0.0)


def r0_const(pars, T, k):
    t = traits_at(pars, T)
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = 1 / t["lf"]
        r2 = (t["a"] ** 2 * t["b"] * t["c"] * k * m_eq(pars, T) * t["PDR"]
              / (mu * (t["PDR"] + mu) * GAMMA))
    return np.sqrt(np.where(np.isfinite(r2) & (r2 > 0), r2, 0.0))


def climatology(data):
    d = pd.read_csv(f"{data}/station_daily_tmean.csv", parse_dates=["date"])
    T = pd.to_numeric(d["tmean_c"], errors="coerce").values
    T[T < -50] = np.nan
    # remove the stuck-probe episode: >= 5 identical consecutive readings
    same = np.r_[False, T[1:] == T[:-1]]
    start = 0
    for i in range(1, len(T) + 1):
        if i == len(T) or not same[i]:
            if i - start >= 5:
                T[start:i] = np.nan
            start = i
    T = T + np.where(d["date"] < "2018-03-15", 1.1, 0.0)
    ok = np.isfinite(T)
    yl = np.where(d["date"].dt.is_leap_year, 366.0, 365.0)
    t = ((d["date"].dt.dayofyear.values - 0.5) * PERIOD / yl)[ok]
    w = 2 * np.pi * t / PERIOD
    X = np.column_stack([np.ones_like(t)] + [f(h * w) for h in (1, 2, 3) for f in (np.cos, np.sin)])
    beta = np.linalg.lstsq(X, T[ok], rcond=None)[0]

    def Tf(tt):
        ww = 2 * np.pi * np.asarray(tt, float) / PERIOD
        s = beta[0]
        for j, h in enumerate((1, 2, 3)):
            s = s + beta[1 + 2 * j] * np.cos(h * ww) + beta[2 + 2 * j] * np.sin(h * ww)
        return s
    return Tf, float(beta[0])


def _grid():
    n = int(PERIOD * STEPS_PER_DAY)
    return np.linspace(0.0, PERIOD, n + 1), PERIOD / n


def periodic_vector(pars, Tf):
    """Adult-female trajectory m*(t) (delta = 1) of the disease-free periodic state,
    sampled on the RK4 grid including half steps. Implicit (Radau) integration: the
    immature equation becomes stiff in summer through the -L^2 term."""
    def rhs(t, y):
        tr = traits_at(pars, Tf(t))
        return [tr["EFD"] * y[1] - tr["MDR"] / tr["pEA"] * y[0] - y[0] ** 2,
                0.5 * tr["MDR"] * y[0] - y[1] / tr["lf"]]
    y = np.array([1.0, 1.0])
    for _ in range(600):
        s = solve_ivp(rhs, (0.0, PERIOD), y, method="Radau", rtol=1e-10, atol=1e-12)
        ynew = s.y[:, -1]
        conv = np.max(np.abs(ynew - y) / ynew) < 1e-9
        y = ynew
        if conv:
            break
    tg, _ = _grid()
    th = np.linspace(0.0, PERIOD, 2 * (len(tg) - 1) + 1)
    s = solve_ivp(rhs, (0.0, PERIOD), y, method="Radau", rtol=1e-11, atol=1e-13, t_eval=th)
    return th, s.y[1]


def survey_scale(data, th, m):
    s = pd.read_csv(f"{data}/aspirator_survey_2021.csv", parse_dates=["week_start"])
    expected = []
    for w in s["week_start"]:
        t0 = w.dayofyear - 1.0
        tt = np.linspace(t0, t0 + 7, 57) % PERIOD
        expected.append(np.interp(tt, th, m).mean())
    expected = 0.55 * s["residents"].values * np.array(expected)
    return s["females_collected"].sum() / expected.sum()


def periodic_r0(pars, Tf, th, m, k):
    tg, h = _grid()
    tr = traits_at(pars, Tf(th))
    a, b, c, pdr, mu = tr["a"], tr["b"], tr["c"], tr["PDR"], 1 / tr["lf"]
    Mv = k * m

    def A(i, lam):
        return np.array([[-SIGMA_H, 0, 0, a[i] * b[i] / lam],
                         [SIGMA_H, -GAMMA, 0, 0],
                         [0, a[i] * c[i] * Mv[i] / lam, -(pdr[i] + mu[i]), 0],
                         [0, 0, pdr[i], -mu[i]]])

    def log_rho(lam):
        X, logs = np.eye(4), 0.0
        for s in range(len(tg) - 1):
            i = 2 * s
            A0, A1, A2 = A(i, lam), A(i + 1, lam), A(i + 2, lam)
            k1 = A0 @ X; k2 = A1 @ (X + 0.5 * h * k1)
            k3 = A1 @ (X + 0.5 * h * k2); k4 = A2 @ (X + h * k3)
            X = X + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            if s % 240 == 0:
                n = np.abs(X).max(); X /= n; logs += np.log(n)
        return np.log(np.abs(np.linalg.eigvals(X)).max()) + logs

    lo, hi = 0.01, 1.0
    while log_rho(hi) > 0:
        lo, hi = hi, 2 * hi
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if log_rho(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)


def full_reference(data="/app/data"):
    pars = fit_curves(data)
    Tf, Tmean = climatology(data)
    th, m = periodic_vector(pars, Tf)
    k = survey_scale(data, th, m)
    dens = k * m
    return {
        "pars": pars, "k": k, "Tf": Tf, "Tmean": Tmean,
        "mean_density": float(np.trapezoid(dens, th) / PERIOD),
        "peak_day": int(np.floor(th[np.argmax(dens)])) + 1,
        "r0_mean_T": float(r0_const(pars, Tmean, k)),
        "r0_seasonal": float(periodic_r0(pars, Tf, th, m, k)),
    }
