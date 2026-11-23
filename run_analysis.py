"""
Runs the data center cooling study and writes:
  results/scenarios.csv     PUE for every supply-temp x server-dT x outdoor-temp combination
  results/containment.csv   the containment comparison
  figures/*.png             four charts used in the README
Run:  python run_analysis.py
"""
import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dc_model import (Assumptions, power_breakdown, rack_inlet_temp, max_supply_temp,
                      fan_power_kw, ASHRAE_MAX_INLET)

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)

# chart style: light surface, quiet axes, 3 fixed series colors
C = ["#2a78d6", "#eb6834", "#1baf7a"]
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb",
    "axes.edgecolor": "#c9c8c3", "axes.labelcolor": "#52514e", "xtick.color": "#52514e",
    "ytick.color": "#52514e", "text.color": "#0b0b0b", "axes.grid": True, "grid.color": "#ebeae6",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2, "font.size": 10,
})
LIMIT_STYLE = dict(color="#52514e", ls="--", lw=1)

base = Assumptions()

# ---------------------------------------------------------------- 1. scenario sweep
rows = []
supply_temps = np.arange(18, 28, 1.0)
server_dts = [8, 10, 12, 15]
outdoor_temps = [25, 30, 35]
for ts, dt, to in itertools.product(supply_temps, server_dts, outdoor_temps):
    a = base.with_(t_supply=ts, server_dt=dt, t_outdoor=to)
    b = power_breakdown(a)
    rows.append({"t_supply_c": ts, "server_dt_k": dt, "t_outdoor_c": to,
                 "cop": round(b["cop"], 2), "fans_kw": round(b["fans_kw"], 1),
                 "chiller_kw": round(b["chiller_kw"], 1), "pue": round(b["pue"], 3),
                 "inlet_no_containment_c": round(rack_inlet_temp(ts, dt, 0.25), 1),
                 "inlet_with_containment_c": round(rack_inlet_temp(ts, dt, 0.05), 1)})
scen = pd.DataFrame(rows)
scen.to_csv("results/scenarios.csv", index=False)

# ---------------------------------------------------------------- 2. containment study
# Containment does not save energy directly. It lowers recirculation, which lets the
# supply air run warmer while servers stay within the ASHRAE limit. Warmer supply air
# means a higher chiller COP, so the savings come from there.
cases = {"No containment (r = 0.25)": 0.25, "Aisle containment (r = 0.05)": 0.05}
comp = []
for name, r in cases.items():
    ts = min(max_supply_temp(base.server_dt, r), ASHRAE_MAX_INLET)
    ts = np.floor(ts * 2) / 2       # round down to a 0.5 °C setpoint
    b = power_breakdown(base.with_(recirculation=r, t_supply=ts))
    comp.append({"case": name, "recirculation": r, "max_supply_c": ts,
                 "rack_inlet_c": round(b["rack_inlet_c"], 1), "cop": round(b["cop"], 2),
                 "pue": round(b["pue"], 3), "total_kw": round(b["total_kw"], 1),
                 **{k: round(v, 1) for k, v in b.items() if k.endswith("_kw")}})
comp = pd.DataFrame(comp)
comp.to_csv("results/containment.csv", index=False)

# For reference: what if the room is left at 18 °C supply (a common conservative setpoint)
legacy = power_breakdown(base.with_(recirculation=0.25, t_supply=18.0))
best = comp.iloc[1]
inlet_drop = rack_inlet_temp(18, base.server_dt, 0.25) - rack_inlet_temp(18, base.server_dt, 0.05)
annual_saving_mwh = (comp.iloc[0].total_kw - best.total_kw) * 8760 / 1000
legacy_saving_mwh = (legacy["total_kw"] - best.total_kw) * 8760 / 1000

# ---------------------------------------------------------------- figures
# Fig 1: rack inlet temperature vs recirculation
fig, ax = plt.subplots(figsize=(6.4, 4))
r = np.linspace(0, 0.4, 100)
for color, dt in zip(C, [8, 12, 15]):
    ax.plot(r, [rack_inlet_temp(18, dt, ri) for ri in r], color=color, label=f"server ΔT = {dt} K")
ax.axhline(ASHRAE_MAX_INLET, **LIMIT_STYLE)
ax.text(0.005, ASHRAE_MAX_INLET + 0.4, "ASHRAE recommended max inlet (27 °C)", color="#52514e", fontsize=9)
ax.set(xlabel="recirculation fraction r", ylabel="rack inlet temperature (°C)",
       title="Hot-air recirculation raises server inlet temperature\n(supply air = 18 °C)")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("figures/1_inlet_vs_recirculation.png", dpi=150); plt.close(fig)

# Fig 2: PUE vs supply temperature for 3 outdoor temperatures
fig, ax = plt.subplots(figsize=(6.4, 4))
sub = scen[scen.server_dt_k == 12]
for color, to in zip(C, outdoor_temps):
    s = sub[sub.t_outdoor_c == to]
    ax.plot(s.t_supply_c, s.pue, color=color, marker="o", ms=4, label=f"outdoor {to} °C")
ax.set(xlabel="supply air temperature (°C)", ylabel="PUE",
       title="Warmer supply air and cooler weather lower PUE\n(server ΔT = 12 K)")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("figures/2_pue_vs_supply_temp.png", dpi=150); plt.close(fig)

# Fig 3: fan power vs server dT (cubic affinity law)
fig, ax = plt.subplots(figsize=(6.4, 4))
dts = np.linspace(8, 16, 50)
ax.plot(dts, [fan_power_kw(base.with_(server_dt=d)) for d in dts], color=C[0])
ax.set(xlabel="air temperature rise across servers, ΔT (K)", ylabel="fan power (kW)",
       title="Fan power falls steeply as ΔT rises\n(less airflow needed; power ∝ flow³)")
fig.tight_layout(); fig.savefig("figures/3_fan_power_vs_dt.png", dpi=150); plt.close(fig)

# Fig 4: overhead breakdown, three cases
labels = ["18 °C supply,\nno containment", f"{comp.iloc[0].max_supply_c:g} °C supply,\nno containment",
          f"{best.max_supply_c:g} °C supply,\naisle containment"]
parts = [("chiller_kw", "chiller"), ("pumps_kw", "pumps"), ("fans_kw", "fans"),
         ("electrical_loss_kw", "UPS / electrical + lighting")]
data = [legacy, comp.iloc[0].to_dict(), best.to_dict()]
fig, ax = plt.subplots(figsize=(6.8, 4.2))
x = np.arange(3)
bottom = np.zeros(3)
part_colors = ["#2a78d6", "#1baf7a", "#eb6834", "#b4b2aa"]
for (key, name), color in zip(parts, part_colors):
    vals = np.array([d[key] + (d["lighting_kw"] if key == "electrical_loss_kw" else 0) for d in data])
    ax.bar(x, vals, 0.55, bottom=bottom, color=color, label=name, edgecolor="#fcfcfb", linewidth=1.5)
    bottom += vals
for xi, d in zip(x, data):
    ax.text(xi, bottom[xi] + 6, f"PUE {d['pue']:.2f}", ha="center", color="#0b0b0b", fontsize=10)
ax.set_axisbelow(True); ax.grid(axis="x", visible=False)
ax.set_xticks(x, labels)
ax.set(ylabel="non-IT power (kW) for 1 MW of IT", title="Where the overhead goes")
ax.set_ylim(0, bottom.max() * 1.15)
ax.legend(frameon=False, ncol=2, fontsize=9, loc="upper right")
fig.tight_layout(); fig.savefig("figures/4_overhead_breakdown.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- summary
print(f"Scenarios evaluated: {len(scen)}")
print(f"PUE range: {scen.pue.min():.2f} – {scen.pue.max():.2f}")
print(f"Chiller COP range: {scen.cop.min():.1f} – {scen.cop.max():.1f}")
print()
print(comp[["case", "max_supply_c", "rack_inlet_c", "cop", "pue", "total_kw"]].to_string(index=False))
print()
print(f"At 18 °C supply, containment lowers rack inlet temperature by {inlet_drop:.1f} °C.")
print(f"Containment + raised setpoint vs no containment at its max setpoint: "
      f"PUE {comp.iloc[0].pue:.3f} -> {best.pue:.3f}, saving ~{annual_saving_mwh:,.0f} MWh/yr per MW of IT.")
print(f"Vs a legacy 18 °C setpoint: PUE {legacy['pue']:.3f} -> {best.pue:.3f}, "
      f"saving ~{legacy_saving_mwh:,.0f} MWh/yr per MW of IT.")
