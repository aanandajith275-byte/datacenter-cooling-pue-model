"""Cooling-failure ride-through: how fast do server inlets heat up when the chiller plant trips?

Lumped-capacitance nodes (each well mixed; Biot numbers are small for thin server sheet metal and
heat sinks, so a single temperature per node is a fair first approximation):

  server mass   C_s dT_s/dt = P_it - Q_s,          Q_s = UA_s * (T_s - T_air,mean)
  hall air      (only when the CRAH fans stop)     C_a dT_a/dt = Q_s - UA_b (T_a - T_b)
  building      C_b dT_b/dt = UA_b (T_a - T_b)     (floor slab, walls)
  chilled water C_w dT_w/dt = Q_coil + P_fans      (pipe + coil inventory + optional storage tank)

While the CRAH fans run, air circulates in seconds, much faster than any mass heats up, so the air
loop is treated as quasi-steady. With coil effectiveness eps and air capacity rate C = m*cp:
  Q_coil = eps C (T_ex - T_w),  T_in = T_w + Q (1/eps - 1)/C,
  and the server convection gives  Q = UA_s (T_s - T_w) / (1 + UA_s (1/eps - 1/2) / C).
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .plant import CP_AIR, CP_W, RHO_AIR, RHO_W


@dataclass
class Hall:
    it_kw: float = 1000.0
    air_volume_m3: float = 3000.0          # ~600 m2 hall, 5 m including raised floor and ceiling plenum
    server_mass_kg: float = 60_000.0       # ~100 racks x 600 kg of chassis, boards, heat sinks, PSUs
    server_cp: float = 600.0               # J/(kg K), mostly steel and aluminium
    ua_server_kw_k: float = 35.0           # chip/heat-sink to air conductance at nominal airflow
    ua_building_kw_k: float = 10.0         # ~2,000 m2 of slab and walls x h ~ 5 W/(m2 K)
    building_heat_cap_mj_k: float = 250.0  # 0.2 m concrete slab
    loop_water_m3: float = 20.0            # chilled water in pipes, coils and headers
    storage_m3: float = 0.0                # chilled-water storage tank on the loop
    dt_server: float = 12.0
    coil_eff: float = 0.8
    crah_fan_kw: float = 50.0
    t_inlet0: float = 24.6                 # normal server inlet temperature


def ride_through(h: Hall, fans_on_ups=True, t_end=3600.0):
    P = h.it_kw * 1e3
    C = P / h.dt_server                                        # air capacity rate m*cp, W/K
    UAs, UAb = h.ua_server_kw_k * 1e3, h.ua_building_kw_k * 1e3
    Cs = h.server_mass_kg * h.server_cp
    Ca = RHO_AIR * CP_AIR * h.air_volume_m3
    Cb = h.building_heat_cap_mj_k * 1e6
    Cw = RHO_W * CP_W * (h.loop_water_m3 + h.storage_m3)
    k = 1 + UAs * (1 / h.coil_eff - 0.5) / C
    # initial steady state with the chiller running: Q = P
    Tw0 = h.t_inlet0 - P * (1 / h.coil_eff - 1) / C
    Ts0 = Tw0 + P * k / UAs
    Ta0 = h.t_inlet0 + h.dt_server / 2                        # mean hall air (fans-off model)

    def rhs(t, y):
        Ts, Ta, Tb, Tw = y
        if fans_on_ups:
            Q = UAs * (Ts - Tw) / k
            return [(P - Q) / Cs, 0.0, 0.0, (Q + h.crah_fan_kw * 1e3) / Cw]
        Q = UAs * (Ts - Ta)
        return [(P - Q) / Cs, (Q - UAb * (Ta - Tb)) / Ca, UAb * (Ta - Tb) / Cb, 0.0]

    sol = solve_ivp(rhs, (0, t_end), [Ts0, Ta0, Ta0 - 2.0, Tw0], max_step=2.0, dense_output=True)
    t = sol.t
    Ts, Ta, Tb, Tw = sol.y
    if fans_on_ups:
        Q = UAs * (Ts - Tw) / k
        t_in = Tw + Q * (1 / h.coil_eff - 1) / C
    else:
        t_in = Ta                                              # recirculating, fully mixed hall
    return t, t_in, Ts


def time_to(t, temp, limit):
    i = np.argmax(temp >= limit)
    return float(t[i]) if temp[i] >= limit else float("inf")
