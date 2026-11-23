"""Sanity checks on the physics. Run:  python -m pytest -q"""
import pytest
from dc_model import (Assumptions, required_airflow, rack_inlet_temp, max_supply_temp,
                      fan_power_kw, chiller_cop, power_breakdown, CP_AIR)


def test_energy_balance_airflow():
    # Q = m_dot * cp * dT must hold exactly
    m_dot, _ = required_airflow(100.0, 10.0)
    assert m_dot * CP_AIR * 10.0 == pytest.approx(100_000.0)


def test_no_recirculation_means_inlet_equals_supply():
    assert rack_inlet_temp(20.0, 12.0, 0.0) == 20.0


def test_recirculation_mixing_balance():
    # T_in should equal the mix of supply air and exhaust (T_in + dT)
    ts, dt, r = 18.0, 12.0, 0.25
    t_in = rack_inlet_temp(ts, dt, r)
    assert t_in == pytest.approx((1 - r) * ts + r * (t_in + dt))


def test_max_supply_keeps_inlet_at_limit():
    ts = max_supply_temp(12.0, 0.25)
    assert rack_inlet_temp(ts, 12.0, 0.25) == pytest.approx(27.0)


def test_fan_power_cubic_law():
    a = Assumptions()
    # halving dT doubles the airflow, so fan power should rise 2^3 = 8x
    ratio = fan_power_kw(a.with_(server_dt=6.0)) / fan_power_kw(a.with_(server_dt=12.0))
    assert ratio == pytest.approx(8.0)


def test_warmer_supply_improves_cop():
    a = Assumptions()
    assert chiller_cop(a.with_(t_supply=24)) > chiller_cop(a.with_(t_supply=18))


def test_pue_is_above_one_and_sums_correctly():
    b = power_breakdown(Assumptions())
    parts = ["it_kw", "fans_kw", "chiller_kw", "pumps_kw", "electrical_loss_kw", "lighting_kw"]
    assert sum(b[p] for p in parts) == pytest.approx(b["total_kw"])
    assert b["pue"] > 1.0
