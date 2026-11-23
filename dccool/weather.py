"""Read hourly weather from EnergyPlus EPW files (8760 rows)."""
import os

import numpy as np
import pandas as pd

from .psychro import humidity_ratio_from_dewpoint, wet_bulb

WEATHER_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "weather")
CLIMATES = {  # short name -> (file, label)
    "tokyo": ("JPN_Tokyo.Hyakuri.477150_IWEC.epw", "Tokyo (humid subtropical)"),
    "tucson": ("USA_AZ_Davis-Monthan.AFB.722745_TMY3.epw", "Tucson, Arizona (hot desert)"),
    "stockholm": ("SWE_Stockholm.Arlanda.024600_IWEC.epw", "Stockholm (cold)"),
}


def load_epw(name):
    fn = CLIMATES[name][0] if name in CLIMATES else name
    path = fn if os.path.exists(fn) else os.path.join(WEATHER_DIR, fn)
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} missing: run `bash data/fetch_weather.sh` first")
    cols = ["year", "month", "day", "hour", "minute", "flags", "t_db", "t_dp", "rh", "p_atm"]
    df = pd.read_csv(path, skiprows=8, header=None, usecols=range(10), names=cols)
    p = df["p_atm"].to_numpy(float)
    W = humidity_ratio_from_dewpoint(df["t_dp"].to_numpy(float), p)
    df["w"] = W
    df["t_wb"] = wet_bulb(df["t_db"].to_numpy(float), W, p)
    return df[["month", "day", "hour", "t_db", "t_dp", "rh", "p_atm", "w", "t_wb"]]
