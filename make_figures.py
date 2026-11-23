"""Draws every figure in figures/ from results/ (run run_analysis.py first)."""
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon

from dccool import style
from dccool.plant import ASHRAE_REC_MAX, Params, servers
from dccool.psychro import humidity_ratio_from_dewpoint, humidity_ratio_sat
from dccool.scenarios import DESIGNS
from dccool.style import AQUA, BLUE, INK, INK2, MAGENTA, MUTED, ORANGE, VIOLET, YELLOW, subtitle
from dccool.weather import CLIMATES, load_epw

style.apply()
os.makedirs("figures", exist_ok=True)
R = lambda f: pd.read_csv(os.path.join("results", f))           # noqa: E731
ANN = R("annual_designs.csv")
DCOL = {"legacy": MUTED, "contained": YELLOW, "wse": BLUE, "ase": ORANGE, "liquid": AQUA}
DLAB = {"legacy": "Legacy (open aisles, 18 °C, 7 °C water, air-cooled chiller)", "contained": "Contained aisles, 22 °C",
        "wse": "+ water-side economiser, 24 °C", "ase": "Air-side economiser (outside air), 24 °C",
        "liquid": "Direct-to-chip liquid cooling + economiser"}
CLAB = {"tokyo": "Tokyo", "tucson": "Tucson (desert)", "stockholm": "Stockholm"}


def save(fig, name):
    fig.savefig(os.path.join("figures", name))
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- 1. headline: annual PUE and water
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.5, 1]})
ax = axes[0]
designs = list(DESIGNS)
x = np.arange(len(CLIMATES))
wbar = 0.16
for i, d in enumerate(designs):
    v = [ANN[(ANN.climate == c) & (ANN.design == d)].pue.iloc[0] for c in CLIMATES]
    ax.bar(x + (i - 2) * wbar, v, wbar * 0.92, color=DCOL[d], label=DLAB[d])
    for xi, vi in zip(x, v):
        ax.text(xi + (i - 2) * wbar, vi + 0.008, f"{vi:.2f}", ha="center", va="bottom", fontsize=7.5, color=INK, rotation=90)
ax.set_xticks(x, [CLAB[c] for c in CLIMATES])
ax.set_ylim(1.0, 1.66)
ax.set_ylabel("annual PUE (facility energy / IT energy)")
ax.grid(axis="x", visible=False)
ax.legend(loc="upper center", bbox_to_anchor=(0.62, -0.09), fontsize=8.3, ncol=2)
ax.set_title("Annual PUE from 8,760 hours of real weather", fontsize=11.5, pad=10)
ax = axes[1]
MARK = {"tokyo": "o", "tucson": "s", "stockholm": "^"}
for d in designs:
    for c in CLIMATES:
        r = ANN[(ANN.design == d) & (ANN.climate == c)].iloc[0]
        ax.scatter(r.pue, r.wue_l_per_kwh, s=75, color=DCOL[d], marker=MARK[c], zorder=3, edgecolor="white", linewidth=1.2)
for c in CLIMATES:
    ax.scatter([], [], marker=MARK[c], color=INK2, label=CLAB[c])
ax.legend(loc="center right", fontsize=8.5)
ax.set_xlabel("annual PUE (lower = less energy)")
ax.set_ylabel("water use, WUE (litres per IT kWh)")
ax.set_title("Energy vs water", fontsize=11.5, pad=10)
leg = ANN[ANN.design == "legacy"].pue
fig.suptitle(f"From a legacy hall (PUE {leg.min():.2f}–{leg.max():.2f}) to liquid cooling "
             f"(PUE {ANN[ANN.design == 'liquid'].pue.min():.2f}–{ANN[ANN.design == 'liquid'].pue.max():.2f}): "
             "every step is a heat-transfer decision", x=0.01, ha="left", fontsize=12.5, fontweight="bold", color=INK, y=1.04)
save(fig, "01_annual_pue_and_water.png")

# ---------------------------------------------------------------- 2. psychrometric chart, Tokyo, coloured by WSE mode
H = R("hourly_tokyo_wse.csv")
w = load_epw("tokyo")
fig, ax = plt.subplots(figsize=(9.2, 5.6))
tt = np.linspace(-10, 40, 300)
ax.plot(tt, humidity_ratio_sat(tt) * 1000, color=INK2, lw=1.2)
for rh in [0.2, 0.4, 0.6, 0.8]:
    ax.plot(tt, rh * humidity_ratio_sat(tt) * 1000, color=MUTED, lw=0.6, ls=":")
    tl = min(38.5, float(np.interp(25.5, rh * humidity_ratio_sat(tt) * 1000, tt)))
    ax.text(tl, rh * humidity_ratio_sat(tl) * 1000, f"{int(rh * 100)}% RH", fontsize=7.5, color=MUTED, va="bottom", ha="right")
mode_col = {"free cooling": AQUA, "partial": YELLOW, "mechanical": ORANGE}
for m, col in mode_col.items():
    s = H["mode"] == m
    ax.scatter(w.t_db[s], w.w[s] * 1000, s=4, color=col, alpha=0.55, label=f"{m} ({s.sum():,} h)", rasterized=True)
# ASHRAE recommended envelope (inlet): 18-27 C, dew point -9..15 C, RH <= 60 %
w_lo, w_hi = humidity_ratio_from_dewpoint(-9.0) * 1000, humidity_ratio_from_dewpoint(15.0) * 1000
t60 = np.linspace(18, 27, 50)
w60 = np.minimum(0.6 * humidity_ratio_sat(t60) * 1000, w_hi)
poly = [(18, w_lo), (27, w_lo)] + [(t, ww) for t, ww in zip(t60[::-1], w60[::-1])]
ax.add_patch(Polygon(poly, closed=True, fill=False, edgecolor=BLUE, lw=2))
ax.text(22.5, w_lo + 0.4, "ASHRAE recommended\nserver-inlet envelope", color=BLUE, fontsize=8.5, fontweight="bold", ha="center")
ax.set_xlim(-10, 40)
ax.set_ylim(0, 26)
ax.set_xlabel("outdoor dry-bulb temperature (°C)")
ax.set_ylabel("humidity ratio (g water / kg dry air)")
ax.legend(loc="upper left", markerscale=4)
ax.set_title(f"Tokyo: the water-side economiser carries the full load {(H['mode'] == 'free cooling').mean():.0%} of the year")
subtitle(ax, "One dot per hour of the typical year. The limit is the WET-bulb temperature: humid summers, not hot days")
save(fig, "02_psychrometric_tokyo.png")

# ---------------------------------------------------------------- 3. supply temperature sweep
S = R("setpoint_sweep_tokyo.csv")
J = json.load(open("results/summary.json"))["setpoint_paradox"]
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
ax = axes[0]
base = S[S.t_supply == 18].iloc[0]
for col, lab, c in [("chiller_mwh", "chiller", ORANGE), ("crah_fan_mwh", "room (CRAH) fans", BLUE),
                    ("server_fan_mwh", "server fans (counted as IT)", VIOLET), ("leakage_mwh", "CPU leakage (counted as IT)", MAGENTA)]:
    ax.plot(S.t_supply, S[col] - base[col], color=c, label=lab)
ax.plot(S.t_supply, S.total_mwh - base.total_mwh, color=INK, lw=2.6, label="total facility energy")
ax.axhline(0, color=MUTED, lw=0.8)
ax.axvline(J["t_supply_min_energy"], color=INK2, ls=":", lw=1)
ax.set_ylim(-400, 1300)
ax.set_xlabel("supply-air temperature setpoint (°C)")
ax.set_ylabel("change in annual energy vs 18 °C (MWh)")
ax.legend(loc="upper left", fontsize=8.5)
ax.set_title("Warm air saves chiller energy until server fans ramp", fontsize=11.5, pad=10)
ax = axes[1]
SS = S[(S.t_supply >= 18) & (S.t_supply <= 26)]
ax.plot(SS.t_supply, base.overhead_mwh - SS.overhead_mwh, color=BLUE, marker="o", ms=3, label="saving judged by overhead (what PUE sees)")
ax.plot(SS.t_supply, base.total_mwh - SS.total_mwh, color=INK, marker="o", ms=3, label="saving at the utility meter")
ax.fill_between(SS.t_supply, base.total_mwh - SS.total_mwh, base.overhead_mwh - SS.overhead_mwh, color=MAGENTA, alpha=0.15)
ax.annotate("gap = extra server-fan + leakage\nenergy, booked as 'IT'", xy=(23.6, 270), xytext=(18.3, 300), fontsize=8.5,
            color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax.axhline(0, color=MUTED, lw=0.8)
ax.set_xlabel("supply-air temperature setpoint (°C)")
ax.set_ylabel("annual saving vs 18 °C (MWh)")
ax.legend(loc="lower left", bbox_to_anchor=(0.17, 0.03), fontsize=8.5)
ax.set_title("PUE overstates the saving", fontsize=11.5, pad=10)
fig.suptitle(f"From 18 °C to {J['t_supply_min_energy']:.1f} °C the facility overhead falls {J['overhead_saving_18_to_best_mwh']:.0f} MWh, "
             f"but server leakage and fans add {J['it_increase_18_to_best_mwh']:.0f} MWh of 'IT' energy",
             x=0.01, ha="left", fontsize=12.5, fontweight="bold", color=INK, y=1.04)
save(fig, "03_setpoint_tradeoff.png")

# ---------------------------------------------------------------- 4. energy breakdown, Tokyo
fig, ax = plt.subplots(figsize=(10, 4.6))
parts = [("chiller_mwh", "chiller", ORANGE), ("heat_rejection_mwh", "tower / condenser fans", YELLOW),
         ("pumps_mwh", "pumps", MAGENTA), ("liquid_plant_mwh", "liquid loop (dry cooler, CDU pumps)", AQUA),
         ("crah_fan_mwh", "room fans", BLUE), ("ups_mwh", "UPS + distribution losses", VIOLET), ("lighting_mwh", "lighting", MUTED)]
T = ANN[ANN.climate == "tokyo"].set_index("design").loc[designs]
y = np.arange(len(designs))[::-1]
left = np.zeros(len(designs))
for col, lab, c in parts:
    ax.barh(y, T[col], 0.6, left=left, color=c, label=lab, edgecolor="white", linewidth=1.2)
    left += T[col].to_numpy()
for yi, (d, r) in zip(y, T.iterrows()):
    ax.text(left[list(designs).index(d)] + 40, yi, f"{r.overhead_mwh:,.0f} MWh  (PUE {r.pue:.2f})", va="center", fontsize=9, color=INK)
ax.set_yticks(y, [DLAB[d].split(" (")[0] for d in designs], fontsize=9)
ax.set_xlabel("annual non-IT energy for 1 MW of IT (MWh)")
ax.set_xlim(0, left.max() * 1.35)
ax.grid(axis="y", visible=False)
ax.legend(loc="lower right", fontsize=8, ncol=2)
ax.set_title("Tokyo: where the overhead goes")
subtitle(ax, "Once the economiser removes most chiller hours, fans and electrical losses are the overhead that remains")
save(fig, "04_overhead_breakdown_tokyo.png")

# ---------------------------------------------------------------- 5. optimal setpoint vs wet bulb
O = R("optimal_setpoint_tokyo_wse.csv")
OP = R("setpoint_optimisation.csv")
fig, ax = plt.subplots(figsize=(8.6, 4.4))
jit = np.random.default_rng(0).normal(0, 0.08, len(O))
ax.scatter(O.t_wb, O.t_opt + jit, s=4, color=BLUE, alpha=0.35, rasterized=True)
ax.set_xlabel("outdoor wet-bulb temperature (°C)")
ax.set_ylabel("energy-optimal supply air (°C)")
ax.axhline(24, color=ORANGE, lw=1.5, ls="--")
ax.text(O.t_wb.min(), 24.25, "fixed setpoint (24 °C)", color=INK2, fontsize=8.5)
tok = OP[(OP.climate == "tokyo") & (OP.design == "wse")].iloc[0]
ax.set_title(f"Hourly optimisation saves {tok.saving_pct:.1f}% more energy, mostly by running COLDER when cooling is free")
subtitle(ax, "Tokyo, water-side economiser. When the tower alone carries the load, cooler air cuts CPU leakage at no plant cost")
save(fig, "05_optimal_setpoint.png")

# ---------------------------------------------------------------- 6. ride-through
RT = R("ride_through.csv")
fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), gridspec_kw={"width_ratios": [1.5, 1]})
ax = axes[0]
for (_, r), c in zip(RT.iterrows(), [ORANGE, YELLOW, BLUE, AQUA]):
    d = R(f"ride_through_{int(r.storage_m3)}_{int(r.fans_on_ups)}.csv")
    ax.plot(d.t_s / 60, d.t_inlet, color=c, label=r.case)
for lim, lab in [(27, "ASHRAE recommended max"), (32, "A1 allowable max"), (35, "")]:
    ax.axhline(lim, color=MUTED, ls=":", lw=1)
    if lab:
        ax.text(118, lim + 0.3, lab, fontsize=8, color=INK2, ha="right")
ax.set_xlim(0, 120)
ax.set_ylim(20, 45)
ax.set_xlabel("minutes after the chiller plant trips")
ax.set_ylabel("server inlet temperature (°C)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=8.3, ncol=2)
ax.set_title("Thermal mass buys time", fontsize=11.5, pad=10)
ax = axes[1]
V = R("ride_through_vs_storage.csv")
ax.plot(V.storage_m3, V.t_to_32_min, color=BLUE, marker="o")
ax.set_xlabel("chilled-water storage volume (m³)")
ax.set_ylabel("minutes until inlet reaches 32 °C")
slope = np.polyfit(V.storage_m3, V.t_to_32_min, 1)[0]
ax.set_title(f"Each 100 m³ of chilled water ≈ {100 * slope:.0f} more minutes", fontsize=11.5, pad=10)
fans_off = RT[~RT.fans_on_ups.astype(bool)].iloc[0]
fig.suptitle(f"If the room fans lose power, inlets pass 32 °C in {fans_off.t_to_32_s:.0f} s; with fans and pumps on UPS the loop water "
             f"gives {RT.iloc[1].t_to_32_s / 60:.0f} min", x=0.01, ha="left", fontsize=12.5, fontweight="bold", color=INK, y=1.04)
save(fig, "06_ride_through.png")

# ---------------------------------------------------------------- 7. uncertainty
MC = R("monte_carlo_tokyo_wse.csv")
SE = R("sensitivity_tokyo_wse.csv")
fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), gridspec_kw={"width_ratios": [1, 1.25]})
ax = axes[0]
ax.hist(MC.pue, bins=30, color=BLUE, edgecolor="white")
for q in [0.05, 0.5, 0.95]:
    ax.axvline(MC.pue.quantile(q), color=INK2, ls=":" if q != 0.5 else "-", lw=1)
ax.set_xlabel("annual PUE (Tokyo, water-side economiser)")
ax.set_ylabel("Monte Carlo draws")
ax.set_title(f"90% range: {MC.pue.quantile(0.05):.3f} – {MC.pue.quantile(0.95):.3f}", fontsize=11.5, pad=10)
ax = axes[1]
SE = SE.iloc[::-1]
nice = {"ups_eff": "UPS efficiency", "crah_dp_pa": "room-fan pressure rise", "coil_eff": "CRAH coil effectiveness",
        "fan_ramp_start": "server fan-ramp inlet temp.", "fan_eff": "fan efficiency", "tower_approach": "cooling-tower approach",
        "carnot_frac": "chiller (fraction of Carnot)", "hx_approach": "economiser HX approach", "leak_frac": "CPU leakage share",
        "server_fan_frac": "server fan power", "recirculation": "recirculation (contained)", "leak_beta": "leakage temp. coefficient"}
ax.barh(np.arange(len(SE)), SE.spearman_energy, color=[ORANGE if v > 0 else BLUE for v in SE.spearman_energy], height=0.65)
ax.set_yticks(np.arange(len(SE)), [nice[p] for p in SE.parameter], fontsize=9)
ax.axvline(0, color=MUTED, lw=0.8)
ax.set_xlabel("rank correlation with annual energy (Spearman)")
ax.grid(axis="y", visible=False)
ax.set_title("Which assumption matters most", fontsize=11.5, pad=10)
fig.suptitle("400 Monte Carlo draws of 12 uncertain parameters: electrical and fan assumptions outweigh the chiller",
             x=0.01, ha="left", fontsize=12.5, fontweight="bold", color=INK, y=1.04)
save(fig, "07_uncertainty.png")

# ---------------------------------------------------------------- 8. recirculation and server inlet (first-principles basics)
fig, ax = plt.subplots(figsize=(8.4, 4.2))
r = np.linspace(0, 0.4, 100)
for c, dt in zip([BLUE, ORANGE, AQUA], [8, 12, 15]):
    ax.plot(r, 18 + r * dt / (1 - r), color=c, label=f"server ΔT = {dt} K")
ax.axhline(ASHRAE_REC_MAX, color=MUTED, ls="--", lw=1)
ax.text(0.005, ASHRAE_REC_MAX + 0.3, "ASHRAE recommended max inlet (27 °C)", color=INK2, fontsize=8.5)
ax.axvspan(0.02, 0.10, color=AQUA, alpha=0.08)
ax.text(0.06, 18.4, "contained", ha="center", fontsize=8.5, color=INK2)
ax.axvspan(0.20, 0.30, color=ORANGE, alpha=0.08)
ax.text(0.25, 18.4, "open aisles", ha="center", fontsize=8.5, color=INK2)
ax.set_xlabel("recirculation fraction r")
ax.set_ylabel("server inlet temperature (°C)")
ax.legend(loc="upper left")
ax.set_title("Recirculation: T_in = T_supply + r·ΔT/(1 − r)")
subtitle(ax, "Supply air at 18 °C. Containment exists to keep r small, which is what lets the supply air run warmer")
save(fig, "08_recirculation.png")
