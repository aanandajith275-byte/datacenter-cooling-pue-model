"""Runs every study and writes results/*.csv|json. Figures: python make_figures.py

    bash data/fetch_weather.sh      # once: downloads the three EPW weather files
    python run_analysis.py          # ~1 min
"""
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from dccool.plant import ASHRAE_REC_MAX, ASHRAE_REC_MIN, Params, hall, servers
from dccool.psychro import relative_humidity
from dccool.scenarios import DESIGNS
from dccool.transient import Hall, ride_through, time_to
from dccool.weather import CLIMATES, load_epw

os.makedirs("results", exist_ok=True)
W = {c: load_epw(c) for c in CLIMATES}


def annual(res):
    E, IT = res["total_kw"].sum(), res["it_kw"].sum()
    return {"pue": E / IT, "total_mwh": E / 1e3, "it_mwh": IT / 1e3, "overhead_mwh": (E - IT) / 1e3,
            "compute_mwh": res["compute_kw"].sum() / 1e3,
            "server_fan_mwh": res["server_fan_kw"].sum() / 1e3, "leakage_mwh": res["leak_kw"].sum() / 1e3,
            "crah_fan_mwh": res["crah_fan_kw"].sum() / 1e3, "chiller_mwh": res["chiller_kw"].sum() / 1e3,
            "heat_rejection_mwh": res["heat_rejection_fan_kw"].sum() / 1e3, "pumps_mwh": res["pumps_kw"].sum() / 1e3,
            "liquid_plant_mwh": res["liquid_plant_kw"].sum() / 1e3, "ups_mwh": res["ups_loss_kw"].sum() / 1e3,
            "lighting_mwh": res["lighting_kw"].sum() / 1e3,
            "wue_l_per_kwh": res["water_kg_s"].sum() * 3600 / IT,
            "free_hours": int((res["mode"] == "free cooling").sum()), "partial_hours": int((res["mode"] == "partial").sum()),
            "peak_kw": float(res["total_kw"].max()), "max_inlet_c": float(np.max(res["t_in"]))}


def run(c, design, t_sa=None, params=None):
    strat, p, tsa, _ = DESIGNS[design]
    w = W[c]
    return hall(w.t_db.values, w.t_wb.values, w.w.values, tsa if t_sa is None else t_sa, strat, params or p)


# ---------------------------------------------------------------- 1. weather summary
rows = []
for c, w in W.items():
    rows.append({"climate": c, "label": CLIMATES[c][1], "mean_db": w.t_db.mean(), "max_db": w.t_db.max(),
                 "max_wb": w.t_wb.max(), "hours_wb_below_12": int((w.t_wb < 12).sum()),
                 "hours_dp_above_15": int((w.t_dp > 15).sum())})
pd.DataFrame(rows).to_csv("results/weather_summary.csv", index=False)

# ---------------------------------------------------------------- 2. annual simulation, 5 designs x 3 climates
rows, monthly = [], []
for c in W:
    for d in DESIGNS:
        r = run(c, d)
        rows.append({"climate": c, "design": d, **annual(r)})
        df = pd.DataFrame({"month": W[c].month, "total": r["total_kw"], "it": r["it_kw"]}).groupby("month").sum()
        for m, x in df.iterrows():
            monthly.append({"climate": c, "design": d, "month": m, "pue": x.total / x.it})
        if c == "tokyo" and d == "wse":
            pd.DataFrame({"month": W[c].month, "day": W[c].day, "hour": W[c].hour, "t_db": W[c].t_db, "t_wb": W[c].t_wb,
                          "mode": r["mode"], "wse_fraction": r["wse_fraction"], "pue": r["pue"]}).to_csv(
                "results/hourly_tokyo_wse.csv", index=False)
ann = pd.DataFrame(rows)
ann.to_csv("results/annual_designs.csv", index=False)
pd.DataFrame(monthly).to_csv("results/monthly_pue.csv", index=False)
print(ann[["climate", "design", "pue", "total_mwh", "wue_l_per_kwh", "free_hours", "partial_hours"]].round(3).to_string(index=False))

# ---------------------------------------------------------------- 3. supply-air setpoint sweep (Tokyo, WSE)
sweep = []
for tsa in np.arange(16, 30.01, 0.5):
    r = run("tokyo", "wse", t_sa=tsa)
    s = servers(np.array([tsa]), DESIGNS["wse"][1])
    sweep.append({"t_supply": tsa, "t_inlet": float(s["t_in"][0]), "t_cpu": float(s["t_cpu"][0]), **annual(r)})
sweep = pd.DataFrame(sweep)
sweep.to_csv("results/setpoint_sweep_tokyo.csv", index=False)
i_e, i_p = sweep.total_mwh.idxmin(), sweep.pue.idxmin()
lo = sweep[sweep.t_supply == 18].iloc[0]
best = sweep.loc[i_e]
paradox = {"t_supply_min_energy": float(best.t_supply), "t_supply_min_pue": float(sweep.loc[i_p].t_supply),
           "overhead_saving_18_to_best_mwh": float(lo.overhead_mwh - best.overhead_mwh),
           "it_increase_18_to_best_mwh": float(best.it_mwh - lo.it_mwh),
           "net_saving_18_to_best_mwh": float(lo.total_mwh - best.total_mwh),
           "pue_18": float(lo.pue), "pue_best": float(best.pue)}

# ---------------------------------------------------------------- 4. hourly setpoint optimisation
grid = np.arange(17.0, 27.01, 0.5)
opt_rows, opt_hourly = [], None
for c in W:
    for d in ["wse", "ase"]:
        strat, p, tsa_fixed, _ = DESIGNS[d]
        w = W[c]
        tot = np.stack([hall(w.t_db.values, w.t_wb.values, w.w.values, g, strat, p)["total_kw"] for g in grid])
        tin = np.stack([servers(np.full(1, g), p)["t_in"][0] * np.ones(len(w)) for g in grid])
        tot = np.where((tin <= ASHRAE_REC_MAX) & (tin >= ASHRAE_REC_MIN), tot, np.inf)   # ASHRAE recommended inlet range
        k = tot.argmin(0)
        t_opt = grid[k]
        r_opt = hall(w.t_db.values, w.t_wb.values, w.w.values, t_opt, strat, p)
        r_fix = run(c, d)
        a_opt, a_fix = annual(r_opt), annual(r_fix)
        opt_rows.append({"climate": c, "design": d, "fixed_setpoint": tsa_fixed, "pue_fixed": a_fix["pue"],
                         "pue_optimised": a_opt["pue"], "total_fixed_mwh": a_fix["total_mwh"],
                         "total_optimised_mwh": a_opt["total_mwh"],
                         "saving_mwh": a_fix["total_mwh"] - a_opt["total_mwh"],
                         "saving_pct": 100 * (1 - a_opt["total_mwh"] / a_fix["total_mwh"]),
                         "mean_setpoint": float(t_opt.mean())})
        if c == "tokyo" and d == "wse":
            opt_hourly = pd.DataFrame({"t_db": w.t_db, "t_wb": w.t_wb, "t_opt": t_opt, "month": w.month})
opt = pd.DataFrame(opt_rows)
opt.to_csv("results/setpoint_optimisation.csv", index=False)
opt_hourly.to_csv("results/optimal_setpoint_tokyo_wse.csv", index=False)
print(opt.round(3).to_string(index=False))

# ---------------------------------------------------------------- 5. ride-through after a chiller-plant trip
ride = []
for label, fans, st in [("CRAH fans and pumps not on UPS", False, 0), ("fans + pumps on UPS, loop water only", True, 0),
                        ("+ 50 m3 chilled-water tank", True, 50), ("+ 200 m3 chilled-water tank", True, 200)]:
    t, tin, ts = ride_through(Hall(storage_m3=st), fans_on_ups=fans, t_end=3 * 3600 if fans else 600)
    ride.append({"case": label, "storage_m3": st, "fans_on_ups": fans, "t_to_27_s": time_to(t, tin, 27),
                 "t_to_32_s": time_to(t, tin, 32), "t_to_35_s": time_to(t, tin, 35)})
    pd.DataFrame({"t_s": t, "t_inlet": tin, "t_server_mass": ts}).to_csv(f"results/ride_through_{st}_{int(fans)}.csv", index=False)
vols = np.array([0, 10, 25, 50, 100, 150, 200, 300])
curve = []
for v in vols:
    t, tin, _ = ride_through(Hall(storage_m3=v), True, t_end=6 * 3600)
    curve.append({"storage_m3": v, "t_to_32_min": time_to(t, tin, 32) / 60})
pd.DataFrame(ride).to_csv("results/ride_through.csv", index=False)
pd.DataFrame(curve).to_csv("results/ride_through_vs_storage.csv", index=False)
print(pd.DataFrame(ride).to_string(index=False))

# ---------------------------------------------------------------- 6. uncertainty: Monte Carlo over the assumptions
RANGES = {"carnot_frac": (0.35, 0.55), "tower_approach": (3.0, 6.0), "hx_approach": (1.0, 3.0), "coil_eff": (0.7, 0.9),
          "crah_dp_pa": (350.0, 650.0), "fan_eff": (0.5, 0.7), "recirculation": (0.02, 0.10),
          "fan_ramp_start": (23.0, 27.0), "leak_beta": (0.02, 0.05), "leak_frac": (0.05, 0.15),
          "server_fan_frac": (0.01, 0.04), "ups_eff": (0.93, 0.97)}
rng = np.random.default_rng(0)
mc = []
for i in range(400):
    draw = {k: float(rng.uniform(*v)) for k, v in RANGES.items()}
    p = Params().with_(**draw)
    a = annual(run("tokyo", "wse", params=p))
    mc.append({**draw, "pue": a["pue"], "total_mwh": a["total_mwh"]})
mc = pd.DataFrame(mc)
mc.to_csv("results/monte_carlo_tokyo_wse.csv", index=False)
sens = pd.DataFrame([{"parameter": k, "spearman_pue": spearmanr(mc[k], mc.pue)[0],
                      "spearman_energy": spearmanr(mc[k], mc.total_mwh)[0]} for k in RANGES])
sens = sens.reindex(sens.spearman_energy.abs().sort_values(ascending=False).index)
sens.to_csv("results/sensitivity_tokyo_wse.csv", index=False)
print(sens.round(3).to_string(index=False))

# ---------------------------------------------------------------- 7. psychrometric: when is outside air usable?
psy = []
for c, w in W.items():
    rh24 = relative_humidity(24.0 + 0.6, w.w.values)
    ok = (w.t_dp.values <= 15) & (rh24 <= 0.60)
    psy.append({"climate": c, "hours_humidity_ok": int(ok.sum()), "hours_db_below_24": int((w.t_db < 24).sum()),
                "hours_both": int((ok & (w.t_db.values < 24)).sum()), "hours_wb_below_wse_full": int((w.t_wb < 12.5).sum())})
pd.DataFrame(psy).to_csv("results/economiser_hours.csv", index=False)

json.dump({"setpoint_paradox": paradox, "monte_carlo": {"pue_p5": float(mc.pue.quantile(0.05)), "pue_p50": float(mc.pue.median()),
                                                         "pue_p95": float(mc.pue.quantile(0.95))}},
          open("results/summary.json", "w"), indent=1)
print(json.dumps(paradox, indent=1))
