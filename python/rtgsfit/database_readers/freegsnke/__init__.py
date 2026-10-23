"""Synthetic MAST-U-like data from a FreeGSNKE equilibrium; the preshot comes from GSFit.

This is the example of a machine other than ST40, and it needs no Tokamak Energy network:

1. FreeGSNKE builds the machine from its pickles in `freegsnke/machine_dir` and solves a static
   free-boundary equilibrium with `ConstrainPaxisIp` profiles and the coil currents in
   `freegsnke/currents_path`;
2. GSFit (settings `gsfit/settings_path`, e.g. `mastu_with_synthetic_data_from_freegsnke`) reads the
   machine and the synthetic measurements from the FreeGSNKE equilibrium (GSFit stores them on
   the time base `[0, 1]`, so they are constant over the replay);
3. the preshot and the signals are then taken from GSFit exactly as for the `gsfit` reader.

Paths may use `~`. `freegsnke`, `gsfit` and `gsfit_rs` are imported lazily.
"""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any

import numpy as np

from rtgsfit.database_readers import PreshotReader
from rtgsfit.database_readers import Signals
from rtgsfit.database_readers import SignalsReader
from rtgsfit.database_readers.gsfit import gsfit_controller
from rtgsfit.database_readers.gsfit import gsfit_signals
from rtgsfit.database_readers.gsfit import preshot_inputs
from rtgsfit.preshot import Preshot
from rtgsfit.preshot import assemble_preshot
from rtgsfit.settings import Settings

MACHINE_FILES = {
    "active_coils_path": "active_coils",
    "passive_coils_path": "passive_coils",
    "limiter_path": "limiter",
    "wall_path": "wall",
    "magnetic_probe_path": "magnetic_probes",
}
"""`build_machine.tokamak` argument -> `<prefix><name>.pickle` in `freegsnke/machine_dir`."""

_EQUILIBRIA: dict[str, Any] = {}


def freegsnke_equilibrium(settings: Settings) -> Any:
    """Solve the FreeGSNKE equilibrium described by `settings["freegsnke"]` (cached)."""
    if settings.xml_text in _EQUILIBRIA:
        return _EQUILIBRIA[settings.xml_text]

    from freegsnke import GSstaticsolver
    from freegsnke import build_machine
    from freegsnke import equilibrium_update
    from freegsnke.jtor_update import ConstrainPaxisIp

    section = settings["freegsnke"]
    machine_dir = Path(section["machine_dir"]).expanduser()
    prefix = section.get("machine_prefix", "")
    tokamak = build_machine.tokamak(
        **{argument: str(machine_dir / f"{prefix}{name}.pickle") for argument, name in MACHINE_FILES.items()}
    )
    grid = section["grid"]
    equilibrium = equilibrium_update.Equilibrium(
        tokamak=tokamak, Rmin=grid["r_min"], Rmax=grid["r_max"], Zmin=grid["z_min"], Zmax=grid["z_max"],
        nx=grid["n_r"], ny=grid["n_z"],
    )  # fmt: skip
    profile = section["profile"]
    profiles = ConstrainPaxisIp(
        eq=equilibrium, paxis=profile["p_axis"], Ip=profile["ip"], fvac=profile["f_vacuum"],
        alpha_m=profile["alpha_m"], alpha_n=profile["alpha_n"],
    )  # fmt: skip
    with open(Path(section["currents_path"]).expanduser(), "rb") as file:
        currents = pickle.load(file)
    for coil, current in currents.items():
        equilibrium.tokamak[coil].current = current

    solver = GSstaticsolver.NKGSsolver(equilibrium)
    solver.solve(eq=equilibrium, profiles=profiles, constrain=None,
                 target_relative_tolerance=float(section.get("tolerance", 1e-9)))  # fmt: skip
    equilibrium._updatePlasmaPsi(equilibrium.plasma_psi)
    tokamak.probes.initialise_setup(equilibrium)

    _EQUILIBRIA.clear()
    _EQUILIBRIA[settings.xml_text] = equilibrium
    return equilibrium


STATIC_TIME = np.array([0.0, 1.0])
"""GSFit's FreeGSNKE reader stores the measurements on this time base, one equilibrium per time."""

_EQUILIBRIUM_LISTS: dict[str, list[Any]] = {}


def _controller(settings: Settings) -> Any:
    """The same equilibrium at both static times, so the measurements are constant over the replay.

    The list is kept so that GSFit's controller cache (keyed on its identity) is hit.
    """
    equilibria = _EQUILIBRIUM_LISTS.get(settings.xml_text)
    if equilibria is None:
        equilibria = [freegsnke_equilibrium(settings)] * len(STATIC_TIME)
        _EQUILIBRIUM_LISTS.clear()
        _EQUILIBRIUM_LISTS[settings.xml_text] = equilibria
    return gsfit_controller(settings, time=STATIC_TIME, freegsnke_eqs=equilibria)


class FreegsnkeReader(PreshotReader, SignalsReader):
    def read_preshot(self, settings: Settings) -> Preshot:
        source = f"freegsnke {settings['freegsnke']['machine_dir']} via gsfit {settings['gsfit']['settings_path']}"
        return assemble_preshot(preshot_inputs(_controller(settings), settings, source))

    def read_signals(self, settings: Settings, preshot: Preshot) -> Signals:
        return gsfit_signals(_controller(settings), settings, preshot)
