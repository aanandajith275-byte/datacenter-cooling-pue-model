"""Moist-air psychrometrics (SI), after ASHRAE Handbook - Fundamentals (2017), Chapter 1.

All functions accept numpy arrays. Temperatures in deg C, pressures in Pa.
"""
import numpy as np

P_ATM = 101_325.0


def p_ws(t):
    """Saturation pressure of water vapour (Pa): Hyland-Wexler over ice below 0 C, over water above."""
    T = np.asarray(t, float) + 273.15
    ice = -5.6745359e3 / T + 6.3925247 - 9.6778430e-3 * T + 6.2215701e-7 * T**2 + 2.0747825e-9 * T**3 \
        - 9.4840240e-13 * T**4 + 4.1635019 * np.log(T)
    water = -5.8002206e3 / T + 1.3914993 - 4.8640239e-2 * T + 4.1764768e-5 * T**2 - 1.4452093e-8 * T**3 \
        + 6.5459673 * np.log(T)
    return np.exp(np.where(T < 273.15, ice, water))


def humidity_ratio_from_dewpoint(t_dp, p=P_ATM):
    pw = p_ws(t_dp)
    return 0.621945 * pw / (p - pw)


def humidity_ratio_sat(t, p=P_ATM):
    return humidity_ratio_from_dewpoint(t, p)


def relative_humidity(t, W, p=P_ATM):
    pw = W * p / (0.621945 + W)
    return pw / p_ws(t)


def enthalpy(t, W):
    """Specific enthalpy of moist air, kJ per kg of dry air."""
    t = np.asarray(t, float)
    return 1.006 * t + W * (2501.0 + 1.86 * t)


def wet_bulb(t, W, p=P_ATM, iters=40):
    """Thermodynamic wet-bulb temperature by bisection on ASHRAE eq. 33/35 (adiabatic saturation)."""
    t = np.asarray(t, float)
    W = np.asarray(W, float) * np.ones_like(t)
    lo = np.full_like(t, -40.0)
    hi = t.copy()
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        Ws = humidity_ratio_sat(mid, p)
        above = mid >= 0
        W_mid = np.where(above,
                         ((2501 - 2.326 * mid) * Ws - 1.006 * (t - mid)) / (2501 + 1.86 * t - 4.186 * mid),
                         ((2830 - 0.24 * mid) * Ws - 1.006 * (t - mid)) / (2830 + 1.86 * t - 2.1 * mid))
        # W_mid increases with mid: if it exceeds the actual W, the wet bulb is lower
        hi = np.where(W_mid > W, mid, hi)
        lo = np.where(W_mid > W, lo, mid)
    return 0.5 * (lo + hi)


def dew_point(W, p=P_ATM, iters=50):
    W = np.asarray(W, float)
    pw = W * p / (0.621945 + W)
    lo, hi = np.full_like(W, -60.0), np.full_like(W, 60.0)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        big = p_ws(mid) > pw
        hi, lo = np.where(big, mid, hi), np.where(big, lo, mid)
    return 0.5 * (lo + hi)
