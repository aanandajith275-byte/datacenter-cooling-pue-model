"""Hour-by-hour model of a 1 MW data hall: servers, air system, chilled-water plant, economisers,
direct liquid cooling. Every function is vectorised over hours (numpy arrays of weather).

Chain of temperatures for the air-cooled hall (supply air T_sa is the control setpoint):

  server inlet      T_in  = T_sa + r * dT_srv / (1 - r)                (hot-air recirculation, mixing balance)
  server dT         dT_srv = min(dT_nom, T_ex,max - T_in)              (server fans ramp to cap exhaust temp)
  server airflow    m_srv = Q_it / (cp * dT_srv)
  CPU temperature   T_cpu = T_in + dT_srv/2 + theta0 * (m_srv/m_nom)^-0.8   (forced convection, h ~ Re^0.8)
  leakage power     P_leak = L_ref * exp(beta * (T_cpu - T_ref))         (sub-threshold leakage)
  CRAH coil         T_sa = T_ra - eps * (T_ra - T_chws)  ->  T_chws = T_ra - (T_ra - T_sa) / eps
  chiller           COP = eta_c * T_evap / (T_cond - T_evap), capped
  cooling tower     T_cw = T_wb + approach;   water-side economiser covers the load fraction
                    f = (T_chwr - (T_cw + a_hx)) / (T_chwr - T_chws), clipped to [0, 1]
"""
from dataclasses import dataclass, replace

import numpy as np

from .psychro import dew_point, humidity_ratio_from_dewpoint, relative_humidity

CP_AIR, RHO_AIR = 1006.0, 1.18            # J/(kg K), kg/m3
CP_W, RHO_W = 4186.0, 997.0
H_FG = 2.45e6                             # J/kg, latent heat of evaporation near 25 C
ASHRAE_REC_MIN, ASHRAE_REC_MAX, ASHRAE_A1_MAX = 18.0, 27.0, 32.0


@dataclass(frozen=True)
class Params:
    # ---- IT and servers
    it_kw: float = 1000.0             # IT power at reference conditions (inlet <= fan-ramp start, CPU at T_ref)
    dt_nom: float = 12.0              # server air temperature rise with fans at minimum speed, K
    fan_ramp_start: float = 25.0      # inlet temperature above which server fans speed up, C
    server_fan_frac: float = 0.02     # server fan power at minimum speed, fraction of IT
    cpu_frac: float = 0.60            # share of IT power drawn by processors
    leak_frac: float = 0.10           # leakage share of CPU power at T_ref
    t_cpu_ref: float = 60.0           # CPU case temperature at reference conditions, C
    leak_beta: float = 0.035          # 1/K: leakage doubles every ln2/beta = 20 K
    theta0: float = 29.0              # K, CPU-to-local-air rise at nominal airflow
    # ---- air system
    recirculation: float = 0.05       # with hot/cold aisle containment (0.25 without)
    oversupply: float = 1.0           # CRAH airflow / server airflow (bypass air; ~1.2 without containment)
    crah_dp_pa: float = 500.0         # CRAH fan pressure rise at design airflow
    fan_eff: float = 0.60
    coil_eff: float = 0.80            # CRAH coil effectiveness (air side is C_min)
    economiser_dp_pa: float = 150.0   # extra filter/damper pressure drop with outside air
    # ---- chilled-water plant
    chw_dt: float = 6.0               # chilled-water temperature rise, K
    evap_approach: float = 3.0        # T_chws - T_evap, K
    cond_approach: float = 5.0        # T_cond - T_cw (water-cooled), K
    aircooled_cond_approach: float = 12.0
    carnot_frac: float = 0.45         # chiller COP as a fraction of Carnot
    min_lift: float = 15.0            # minimum T_cond - T_evap the compressor needs, K (caps COP near 8.5)
    chws_reset: bool = True           # False: chilled water held at min_chws all year (legacy practice)
    tower_approach: float = 4.0       # T_cw - T_wb, K
    hx_approach: float = 1.5          # water-side economiser plate heat exchanger, K
    tower_fan_kw_per_kw: float = 0.012
    condenser_fan_kw_per_kw: float = 0.025   # air-cooled chiller condenser fans
    pump_head_kpa: float = 250.0
    pump_eff: float = 0.70
    min_chws: float = 7.0
    # ---- direct liquid cooling
    liquid_frac: float = 0.75         # share of IT heat taken by cold plates
    liquid_supply: float = 32.0       # facility water supply temperature (ASHRAE W32), C
    liquid_dt: float = 8.0
    drycooler_approach: float = 6.0
    drycooler_fan_kw_per_kw: float = 0.015
    t_cpu_coldplate_rise: float = 22.0   # K, coolant -> CPU case under a cold plate
    # ---- electrical
    ups_eff: float = 0.95
    lighting_frac: float = 0.01
    # ---- water
    cycles_of_concentration: float = 4.0
    evaporated_fraction: float = 0.85    # share of tower heat rejected by evaporation

    def with_(self, **kw):
        return replace(self, **kw)


def servers(t_sa, p: Params, liquid=False):
    """Server inlet/exhaust temperatures, airflow and IT power for a supply-air temperature (array)."""
    t_sa = np.asarray(t_sa, float)
    t_exmax = p.fan_ramp_start + p.dt_nom
    air_share = (1 - p.liquid_frac) if liquid else 1.0
    q_air = p.it_kw * 1e3 * air_share
    m_nom = q_air / (CP_AIR * p.dt_nom)
    dt = np.full_like(t_sa, p.dt_nom)
    it = np.full_like(t_sa, p.it_kw)
    for _ in range(30):                                            # fixed point: inlet <-> dT <-> IT power
        t_in = t_sa + p.recirculation * dt / (1 - p.recirculation)
        dt = np.clip(np.minimum(p.dt_nom, t_exmax - t_in), 1.0, None)
        m = it * 1e3 * air_share / (CP_AIR * dt)
        fan = p.it_kw * p.server_fan_frac * air_share * (m / m_nom) ** 3
        if liquid:
            t_cpu = p.liquid_supply + p.t_cpu_coldplate_rise + 0 * t_in
        else:
            t_cpu = t_in + dt / 2 + p.theta0 * (m / m_nom) ** -0.8
        leak_ref = p.it_kw * p.cpu_frac * p.leak_frac
        leak = leak_ref * np.exp(p.leak_beta * (t_cpu - p.t_cpu_ref))
        dyn = p.it_kw - leak_ref - p.it_kw * p.server_fan_frac       # compute power, independent of temperature
        it = dyn + leak + fan
    return {"t_in": t_in, "dt_server": dt, "m_server": m, "server_fan_kw": fan, "leak_kw": leak, "it_kw": it,
            "t_cpu": t_cpu, "compute_kw": dyn + 0 * t_in}


def chiller_cop(t_chws, t_cond, p: Params):
    t_ev = t_chws - p.evap_approach + 273.15
    lift = t_cond + 273.15 - t_ev
    # smooth maximum of (lift, min_lift): COP keeps improving as the lift falls, but saturates near the
    # compressor's minimum lift instead of diverging like the Carnot COP
    lift_eff = (np.maximum(lift, 0) ** 4 + p.min_lift ** 4) ** 0.25
    return p.carnot_frac * t_ev / lift_eff


def pump_kw(q_w, dt_w, p: Params):
    v = q_w / (RHO_W * CP_W * dt_w)
    return v * p.pump_head_kpa * 1e3 / p.pump_eff / 1e3


def hall(t_db, t_wb, w_out, t_sa, strategy, p: Params):
    """Facility power (kW) for each hour. strategy in
       'legacy'   air-cooled chiller, no economiser
       'chiller'  water-cooled chiller + cooling tower, no economiser
       'wse'      water-cooled chiller + water-side economiser
       'ase'      direct air-side economiser (outside air within ASHRAE humidity limits) + chiller trim
       'liquid'   direct-to-chip liquid cooling (dry cooler, chiller trim) + 'wse' for the residual air load
    """
    t_db, t_wb = np.asarray(t_db, float), np.asarray(t_wb, float)
    t_sa = np.broadcast_to(np.asarray(t_sa, float), t_db.shape).astype(float)
    liquid = strategy == "liquid"
    s = servers(t_sa, p, liquid=liquid)
    out = {k: v for k, v in s.items()}
    it = s["it_kw"]
    air_share = (1 - p.liquid_frac) if liquid else 1.0
    m_crah = s["m_server"] * p.oversupply * (1 - p.recirculation)   # recirculated air does not pass the CRAH
    m_crah = np.maximum(m_crah, 1e-6)
    m_design = p.it_kw * 1e3 / (CP_AIR * p.dt_nom)
    light = p.lighting_frac * p.it_kw
    q_air = (it - (p.liquid_frac * (it - s["server_fan_kw"]) if liquid else 0)) * 1e3   # W into the room air
    # outside-air economiser state (decided before fan power: it adds filter pressure drop)
    # ASHRAE recommended envelope at the SERVER INLET: dew point <= 15 C and RH <= 60 %. Outside air is mixed
    # with return air, so the inlet has the outdoor humidity ratio but the inlet temperature. Air too dry
    # (dew point < -9 C) is humidified adiabatically (water, negligible energy).
    t_dp = dew_point(w_out)
    rh_in = relative_humidity(s["t_in"], w_out)
    ase_ok = (strategy == "ase") & (t_dp <= 15.0) & (rh_in <= 0.60)
    t_ra_guess = t_sa + q_air / (m_crah * CP_AIR)
    ase_use = ase_ok & (t_db < t_ra_guess)
    dp = p.crah_dp_pa + np.where(ase_use, p.economiser_dp_pa, 0.0)
    v = m_crah / RHO_AIR
    v_design = m_design / RHO_AIR
    crah_fan = v * dp * (m_crah / m_design) ** 2 / p.fan_eff / 1e3   # dp ~ V^2 -> P ~ V^3
    q_room = q_air + (crah_fan + light) * 1e3                         # W to remove from the air
    t_ra = t_sa + q_room / (m_crah * CP_AIR)
    # load on the mechanical plant
    if strategy == "ase":
        # 100% outside air pre-cools from t_ra to max(t_db, t_sa); the coil trims the rest
        free_air = np.where(ase_use, np.clip((t_ra - np.maximum(t_db, t_sa)) / (t_ra - t_sa), 0, 1), 0.0)
    else:
        free_air = np.zeros_like(t_db)
    q_coil = q_room * (1 - free_air)
    t_chws = np.maximum(t_ra - (t_ra - t_sa) / p.coil_eff, p.min_chws) if p.chws_reset else np.full_like(t_db, p.min_chws)
    t_chwr = t_chws + p.chw_dt
    q_liquid = (p.liquid_frac * (it - s["server_fan_kw"]) * 1e3) if liquid else np.zeros_like(t_db)
    # heat rejection
    if strategy == "legacy":
        t_cond = t_db + p.aircooled_cond_approach
        cop = chiller_cop(t_chws, t_cond, p.with_(carnot_frac=p.carnot_frac * 0.9))
        chiller = q_coil / cop / 1e3
        reject_fans = p.condenser_fan_kw_per_kw * (q_coil / 1e3 + chiller)
        wse_frac = np.zeros_like(t_db)
        tower_heat = np.zeros_like(t_db)
        liquid_plant = np.zeros_like(t_db)
    else:
        t_cw = t_wb + p.tower_approach
        wse_frac = np.clip((t_chwr - (t_cw + p.hx_approach)) / (t_chwr - t_chws), 0, 1) if strategy in ("wse", "liquid") \
            else np.zeros_like(t_db)
        q_ch = q_coil * (1 - wse_frac)
        # liquid loop: dry cooler if cold enough, else evaporative tower water, else chiller trim
        liquid_plant = np.zeros_like(t_db)
        if liquid:
            t_ret = p.liquid_supply + p.liquid_dt
            dry = np.clip((t_ret - (t_db + p.drycooler_approach)) / p.liquid_dt, 0, 1)
            tower_l = np.clip((t_ret - (t_cw + p.hx_approach)) / p.liquid_dt, 0, 1) * (1 - dry)
            mech_l = 1 - dry - tower_l
            q_ch = q_ch + q_liquid * mech_l
            liquid_plant = p.drycooler_fan_kw_per_kw * q_liquid / 1e3 * dry + pump_kw(q_liquid, p.liquid_dt, p.with_(pump_head_kpa=150))
            liquid_tower = q_liquid * tower_l
        else:
            liquid_tower = 0.0
        cop = chiller_cop(t_chws, t_cw + p.cond_approach, p)
        chiller = q_ch / cop / 1e3
        tower_heat = q_coil * wse_frac + q_ch + chiller * 1e3 + liquid_tower        # W rejected by the tower
        reject_fans = p.tower_fan_kw_per_kw * tower_heat / 1e3
    pumps = pump_kw(q_coil, p.chw_dt, p) + (pump_kw(tower_heat, 5.0, p.with_(pump_head_kpa=200)) if strategy != "legacy" else 0)
    ups = it * (1 / p.ups_eff - 1)
    total = it + crah_fan + chiller + reject_fans + pumps + liquid_plant + light + ups
    makeup = tower_heat * p.evaporated_fraction / H_FG * (1 + 1 / (p.cycles_of_concentration - 1))   # kg/s
    if strategy == "ase":                     # adiabatic humidification of very dry outside air to -9 C dew point
        w_min = humidity_ratio_from_dewpoint(-9.0)
        oa_flow = m_crah * np.where(ase_use, np.maximum(free_air, 0.0), 0.0)
        makeup = makeup + oa_flow * np.clip(w_min - w_out, 0, None)
    out.update({"crah_fan_kw": crah_fan, "chiller_kw": chiller, "heat_rejection_fan_kw": reject_fans,
                "pumps_kw": pumps, "liquid_plant_kw": liquid_plant,
                "lighting_kw": light + 0 * t_db, "ups_loss_kw": ups, "total_kw": total, "pue": total / it,
                "cop": cop, "t_chws": t_chws, "t_ra": t_ra, "wse_fraction": wse_frac, "air_free_fraction": free_air,
                "water_kg_s": makeup, "mode": _mode(strategy, wse_frac, free_air)})
    return out


def _mode(strategy, wse, free_air):
    f = np.maximum(wse, free_air)
    return np.where(f >= 0.999, "free cooling", np.where(f > 0.001, "partial", "mechanical"))
