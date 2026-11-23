"""
First-order thermal and energy model of an air-cooled data center.

All physics here is steady-state energy balance, the same as in a thermodynamics course:
  1. Airflow needed to carry away server heat:   Q = m_dot * cp * dT
  2. Hot-air recirculation (mixing balance):      T_in = T_supply + r * dT / (1 - r)
  3. Fan power (fan affinity laws):               P_fan = V * dP / eta,  with dP ~ V^2
  4. Chiller power from a Carnot-based COP:       COP = eta_c * T_evap / (T_cond - T_evap)
  5. Power Usage Effectiveness:                   PUE = total facility power / IT power

This is a design-stage estimate, not a CFD simulation. Each assumption is a field
in `Assumptions` so it can be changed and its effect tested.
"""
from dataclasses import dataclass, replace

CP_AIR = 1006.0    # specific heat of air, J/(kg K)
RHO_AIR = 1.2      # density of air near 20 °C, kg/m^3
ASHRAE_MAX_INLET = 27.0   # upper end of ASHRAE's recommended server inlet range, °C


@dataclass(frozen=True)
class Assumptions:
    it_load_kw: float = 1000.0      # total server (IT) power; all of it becomes heat
    server_dt: float = 12.0         # air temperature rise across the servers, K
    recirculation: float = 0.25     # fraction of rack inlet air that is recirculated hot exhaust
    t_supply: float = 18.0          # cold air supplied by the cooling units, °C
    t_outdoor: float = 30.0         # outdoor air temperature, °C

    # fans
    design_dp_pa: float = 500.0     # fan pressure rise at the design airflow, Pa
    design_dt: float = 12.0         # server dT at which design_dp_pa applies, K
    fan_eff: float = 0.60           # fan + motor efficiency

    # chiller
    evap_approach: float = 9.0      # supply air -> refrigerant evaporating temperature, K
    cond_approach: float = 8.0      # outdoor air -> refrigerant condensing temperature, K
    carnot_fraction: float = 0.50   # real COP as a fraction of the Carnot COP
    pump_fraction: float = 0.10     # pumps and cooling tower, as a fraction of chiller power

    # electrical
    ups_eff: float = 0.94           # UPS + power distribution efficiency
    lighting_fraction: float = 0.01 # lighting and other loads, as a fraction of IT power

    def with_(self, **kw):
        return replace(self, **kw)


def required_airflow(it_kw, dt):
    """Mass flow (kg/s) and volume flow (m^3/s) of air needed to remove it_kw with a rise of dt."""
    m_dot = it_kw * 1000.0 / (CP_AIR * dt)
    return m_dot, m_dot / RHO_AIR


def rack_inlet_temp(t_supply, dt, r):
    """Server inlet temperature when a fraction r of inlet air is recirculated exhaust.

    Mixing energy balance: T_in = (1-r)*T_supply + r*T_exhaust, with T_exhaust = T_in + dt
    =>  T_in = T_supply + r*dt / (1 - r)
    """
    if not 0 <= r < 1:
        raise ValueError("recirculation fraction must be in [0, 1)")
    return t_supply + r * dt / (1.0 - r)


def max_supply_temp(dt, r, limit=ASHRAE_MAX_INLET):
    """Warmest supply air that still keeps every server inlet at or below `limit`."""
    return limit - r * dt / (1.0 - r)


def fan_power_kw(a: Assumptions):
    """Fan power from the affinity laws: dP scales with V^2, so power scales with V^3."""
    _, v = required_airflow(a.it_load_kw, a.server_dt)
    _, v_design = required_airflow(a.it_load_kw, a.design_dt)
    dp = a.design_dp_pa * (v / v_design) ** 2
    return v * dp / a.fan_eff / 1000.0


def chiller_cop(a: Assumptions):
    """Carnot-based COP. Warmer supply air or cooler outdoor air means a higher COP."""
    t_evap = a.t_supply - a.evap_approach + 273.15
    t_cond = a.t_outdoor + a.cond_approach + 273.15
    if t_cond <= t_evap:
        raise ValueError("condensing temperature must exceed evaporating temperature")
    return a.carnot_fraction * t_evap / (t_cond - t_evap)


def power_breakdown(a: Assumptions):
    """Return each power component in kW and the resulting PUE."""
    it = a.it_load_kw
    fans = fan_power_kw(a)
    heat_to_remove = it + fans          # fan work also ends up as heat in the air
    cop = chiller_cop(a)
    chiller = heat_to_remove / cop
    pumps = a.pump_fraction * chiller
    electrical_loss = it * (1.0 / a.ups_eff - 1.0)
    lighting = a.lighting_fraction * it
    total = it + fans + chiller + pumps + electrical_loss + lighting
    return {
        "it_kw": it, "fans_kw": fans, "chiller_kw": chiller, "pumps_kw": pumps,
        "electrical_loss_kw": electrical_loss, "lighting_kw": lighting,
        "total_kw": total, "cop": cop, "pue": total / it,
        "rack_inlet_c": rack_inlet_temp(a.t_supply, a.server_dt, a.recirculation),
    }


def pue(a: Assumptions):
    return power_breakdown(a)["pue"]
