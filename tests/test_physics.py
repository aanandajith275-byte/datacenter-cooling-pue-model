"""Physics checks. Run:  python -m pytest -q"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from dccool.plant import CP_AIR, CP_W, RHO_W, Params, chiller_cop, hall, servers  # noqa: E402
from dccool.psychro import humidity_ratio_from_dewpoint, p_ws, wet_bulb  # noqa: E402
from dccool.transient import Hall, ride_through  # noqa: E402

P = Params()


def test_saturation_pressure_matches_ashrae_table():
    # ASHRAE Fundamentals 2017, Table 3: 0.01 C 611.65 Pa, 20 C 2339.3 Pa, 40 C 7384.9 Pa, -10 C (ice) 259.9 Pa
    assert p_ws(np.array([0.01, 20.0, 40.0, -10.0])) == pytest.approx([611.65, 2339.3, 7384.9, 259.9], rel=2e-3)


def test_wet_bulb_known_point_and_saturation():
    W = humidity_ratio_from_dewpoint(np.array([18.45]))          # 30 C, ~50 % RH
    assert wet_bulb(np.array([30.0]), W)[0] == pytest.approx(22.0, abs=0.15)
    Ws = humidity_ratio_from_dewpoint(np.array([25.0]))           # saturated air: wet bulb = dry bulb
    assert wet_bulb(np.array([25.0]), Ws)[0] == pytest.approx(25.0, abs=0.02)


def test_recirculation_mixing_balance():
    s = servers(np.array([20.0]), P.with_(recirculation=0.25))
    t_in, dt = s["t_in"][0], s["dt_server"][0]
    assert t_in == pytest.approx(0.75 * 20.0 + 0.25 * (t_in + dt), abs=1e-6)


def test_server_airflow_energy_balance():
    s = servers(np.array([20.0]), P)
    assert s["m_server"][0] * CP_AIR * s["dt_server"][0] == pytest.approx(s["it_kw"][0] * 1e3, rel=1e-6)


def test_server_fans_ramp_only_above_threshold():
    cold, hot = servers(np.array([20.0]), P), servers(np.array([28.0]), P)
    assert cold["dt_server"][0] == pytest.approx(P.dt_nom)
    assert hot["dt_server"][0] < P.dt_nom and hot["server_fan_kw"][0] > 2 * cold["server_fan_kw"][0]


def test_room_air_energy_balance():
    t = np.array([5.0, 20.0, 30.0])
    r = hall(t, t - 3, humidity_ratio_from_dewpoint(t - 6), 22.0, "wse", P)
    m = r["m_server"] * (1 - P.recirculation)
    q_in = (r["it_kw"] + r["crah_fan_kw"] + r["lighting_kw"]) * 1e3
    assert m * CP_AIR * (r["t_ra"] - 22.0) == pytest.approx(q_in, rel=1e-9)


def test_chiller_never_beats_carnot():
    t_chws = np.linspace(7, 20, 30)
    for t_cond in [15.0, 30.0, 45.0]:
        carnot = (t_chws - P.evap_approach + 273.15) / (t_cond + 273.15 - (t_chws - P.evap_approach + 273.15))
        cop = chiller_cop(t_chws, t_cond, P)
        assert np.all(cop < np.where(carnot > 0, carnot, np.inf))


def test_economiser_limits():
    cold = hall(np.array([0.0]), np.array([-2.0]), humidity_ratio_from_dewpoint(np.array([-5.0])), 24.0, "wse", P)
    hot = hall(np.array([33.0]), np.array([27.0]), humidity_ratio_from_dewpoint(np.array([25.0])), 24.0, "wse", P)
    assert cold["wse_fraction"][0] == 1.0 and cold["chiller_kw"][0] == 0.0
    assert hot["wse_fraction"][0] == 0.0 and hot["chiller_kw"][0] > 0


def test_pue_components_add_up():
    t = np.array([10.0, 25.0])
    r = hall(t, t - 4, humidity_ratio_from_dewpoint(t - 8), 24.0, "wse", P)
    parts = ["it_kw", "crah_fan_kw", "chiller_kw", "heat_rejection_fan_kw", "pumps_kw", "liquid_plant_kw",
             "lighting_kw", "ups_loss_kw"]
    assert sum(r[k] for k in parts) == pytest.approx(r["total_kw"])
    assert np.all(r["pue"] > 1)


def test_ride_through_conserves_energy():
    h = Hall(storage_m3=50)
    t, t_in, ts = ride_through(h, fans_on_ups=True, t_end=1800)
    # all heat (IT + fans) ends up in the water and the server mass
    Cw = RHO_W * CP_W * (h.loop_water_m3 + h.storage_m3)
    Cs = h.server_mass_kg * h.server_cp
    added = (h.it_kw + h.crah_fan_kw) * 1e3 * t[-1]
    # recover the water temperature from T_in = T_w + Q (1/eps - 1)/C with Q = UA (T_s - T_w)/k
    C = h.it_kw * 1e3 / h.dt_server
    UAs = h.ua_server_kw_k * 1e3
    k = 1 + UAs * (1 / h.coil_eff - 0.5) / C
    a = (1 / h.coil_eff - 1) / C
    tw = lambda i: (t_in[i] - a * UAs * ts[i] / k) / (1 - a * UAs / k)          # noqa: E731
    stored = Cw * (tw(-1) - tw(0)) + Cs * (ts[-1] - ts[0])
    assert stored == pytest.approx(added, rel=0.01)
