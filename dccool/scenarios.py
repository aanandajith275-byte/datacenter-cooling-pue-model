"""The cooling designs compared in the annual study. Each: plant strategy, parameter overrides, setpoint."""
from .plant import Params

LEGACY = Params(recirculation=0.25, oversupply=1.5, crah_dp_pa=600.0, ups_eff=0.90, chws_reset=False, min_chws=7.0)

DESIGNS = {
    # name: (strategy, params, supply-air setpoint C, label)
    "legacy": ("legacy", LEGACY, 18.0, "Legacy: open aisles, 18 C, 7 C chilled water, air-cooled chiller"),
    "contained": ("chiller", Params(), 22.0, "Contained aisles, 22 C, water-cooled chiller"),
    "wse": ("wse", Params(), 24.0, "+ water-side economiser, 24 C"),
    "ase": ("ase", Params(), 24.0, "Air-side economiser (outside air), 24 C"),
    "liquid": ("liquid", Params(), 24.0, "Direct-to-chip liquid cooling (W32) + economiser"),
}
