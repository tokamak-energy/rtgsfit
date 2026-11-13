"""Build the preshot, and optionally read the signals, from a GSFit controller.

GSFit reads the machine description (coils, passives, sensors, wall) and calculates every Greens
function; this reader rearranges them into an RT-GSFit preshot. It replaces GSFit's old
`rtgsfit_mdsplus` writer, so that everything RT-GSFit needs is defined in this repository.

Settings, all under `<gsfit>`:

* `settings_path`: GSFit settings, e.g. `st40_setup_for_rtgsfit`;
* `pulse` (optional, default `replay/pulse`): the pulse GSFit reads the machine description from;
* `grid/{n_r, n_z, r_min, r_max, z_min, z_max}` (optional): overrides GSFit's grid;
* `sensor_names`: `st40_pcs` (`L007` -> `PSI_FLOOP_007`, `P101` -> `B_BPPROBE_101`,
  `INIVC000` -> `I_ROG_INIVC000`) or `gsfit` (`PSI_<name>`, `B_<name>`, `I_<name>`);
* `rogowski_names`: Rogowski coils in the order they appear in meas (only included ones are used);
* `coils/groups/<COIL>`: the GSFit PF coils connected in series as each RT-GSFit coil;
* `coils/signals`, `coils/matrix`: `PRESHOT.COIL:SIGNALS` and `COIL:MATRIX` (row-major);
  without `coils`, every GSFit PF coil is an RT-GSFit coil driven by a signal of the same name;
* `sensor_replacements/<NAME>/{replacements, coefficients}`;
* `initial/{source, r, z, radius}`: `gsfit_settings` uses GSFit's `rtgsfit_initial_conditions.json`
  (only valid on the grid it was made for), `cosine_bump` a cosine bump centred at (r, z).

Units are those of GSFit: psi in Wb (not Wb/rad), hence the 2 pi in the Poisson right-hand side.
`gsfit`, `gsfit_rs` and `shapely` are imported lazily.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from rtgsfit.constants import MU_0
from rtgsfit.database_readers import PreshotReader
from rtgsfit.database_readers import Signals
from rtgsfit.database_readers import SignalsReader
from rtgsfit.preshot import Preshot
from rtgsfit.preshot import PreshotInputs
from rtgsfit.preshot import SensorReplacement
from rtgsfit.preshot import assemble_preshot
from rtgsfit.preshot import cosine_bump_psi_norm
from rtgsfit.preshot import rectangle_points
from rtgsfit.settings import Settings
from rtgsfit.settings import replay_times

REGULARISATION_SCALE = 0.001 * np.pi
"""Scaling of GSFit's passive regularisation weights, as in GSFit's former `rtgsfit_mdsplus` writer."""


# ----------------------------------------------------------------------------------------------
# Sensor names
# ----------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Sensor:
    """A GSFit sensor: `kind` is `flux_loop`, `bp_probe` or `rogowski`."""

    kind: str
    name: str


_QUANTITY = {"flux_loop": "psi", "bp_probe": "b", "rogowski": "i"}


def sensor_name(sensor: Sensor, style: str) -> str:
    """The name of a constraint in `PRESHOT.CONSTRAINT:NAMES` and `MEAS:NAMES`."""
    if style == "st40_pcs":
        if sensor.kind == "flux_loop":
            return sensor.name.replace("L", "PSI_FLOOP_")
        if sensor.kind == "bp_probe":
            return sensor.name.replace("P", "B_BPPROBE_")
        return f"I_ROG_{sensor.name}"
    if style == "gsfit":
        return f"{_QUANTITY[sensor.kind].upper()}_{sensor.name}"
    raise ValueError(f"unknown gsfit/sensor_names {style!r}; use st40_pcs or gsfit")


def _passive_dof_unit(dof: str) -> str:
    if dof == "constant_current_density":
        return "ampere / metre ** 2"
    if dof.startswith("eig_"):
        return "ampere"
    raise ValueError(f"unsupported passive DOF type {dof!r}")


# ----------------------------------------------------------------------------------------------
# The GSFit controller
# ----------------------------------------------------------------------------------------------
_CONTROLLERS: dict[tuple[Any, ...], Any] = {}


def gsfit_controller(settings: Settings, **gsfit_kwargs: Any) -> Any:
    """A GSFit controller with every object and Greens function set up, but not run.

    `gsfit_kwargs` are passed to GSFit's database reader (e.g. `freegsnke_eqs` and `time`).
    Controllers are cached on the settings, so the preshot and signals readers share one.
    """
    section = settings["gsfit"]
    pulse = int(section.get("pulse", settings["replay"]["pulse"]))
    time = replay_times(settings)
    key = (settings.json_text, id(gsfit_kwargs.get("freegsnke_eqs")))
    if key in _CONTROLLERS:
        return _CONTROLLERS[key]

    from gsfit import Gsfit

    controller = Gsfit(pulseNo=pulse, run_name="RTGSFIT", settings_path=section["settings_path"], write_to_mds=False)
    code_settings = controller.settings["GSFIT_code_settings.json"]
    code_settings["timeslices"]["method"] = "user_defined"
    code_settings["timeslices"]["user_defined"] = [float(t) for t in time]
    for name, value in section.get("grid", {}).items():
        code_settings["grid"][name] = value
    controller.set_environment_variables()
    controller.setup_timeslices()
    controller.setup_objects(**gsfit_kwargs)
    controller.calculate_greens()
    _CONTROLLERS.clear()
    _CONTROLLERS[key] = controller
    return controller


def _included(objects: Any, name: str) -> bool:
    return bool(np.all(objects.get_vec_bool([name, "fit_settings", "include"])))


def constraint_sensors(controller: Any, settings: Settings) -> list[Sensor]:
    """Included flux loops, then Bp probes, then Rogowski coils in the order of `gsfit/rogowski_names`."""
    sensors = [
        Sensor("flux_loop", name) for name in controller.flux_loops.keys() if _included(controller.flux_loops, name)
    ]
    sensors += [
        Sensor("bp_probe", name) for name in controller.bp_probes.keys() if _included(controller.bp_probes, name)
    ]
    rogowski_names = settings["gsfit"].get("rogowski_names", list(controller.rogowski_coils.keys()))
    sensors += [Sensor("rogowski", name) for name in rogowski_names if _included(controller.rogowski_coils, name)]
    return sensors


def _sensor_objects(controller: Any, kind: str) -> Any:
    return {
        "flux_loop": controller.flux_loops,
        "bp_probe": controller.bp_probes,
        "rogowski": controller.rogowski_coils,
    }[kind]


# ----------------------------------------------------------------------------------------------
# Greens functions on the grid edge
# ----------------------------------------------------------------------------------------------
def self_inductance_rectangle(r: float, d_r: float, d_z: float) -> float:
    """Self-inductance [H] of an axisymmetric wire of rectangular cross-section d_r x d_z at radius r."""
    return float(
        MU_0 * r * (
            (1.0 + 2.0 * (d_z / (8.0 * r)) ** 2 + 2.0 / 3.0 * (d_r / (8.0 * r)) ** 2) * np.log(8.0 * r / (d_r + d_z))
            - 0.5 + 0.5 * (d_z / (8.0 * r)) ** 2
        )
    )  # fmt: skip


def edge_greens(grid_r: np.ndarray, grid_z: np.ndarray) -> np.ndarray:
    """[N_EDGE, N_EDGE]: psi [Wb] at each edge point per ampere at each edge point, divided by 2 pi.

    The diagonal is the self-inductance of a surface current along the edge (a strip of width d_z on
    the left and right, d_r on the top and bottom), which is the average over the strip that the
    boundary integral needs. Three corners are not used by RT-GSFit and keep the value of
    `gsfit_rs.greens_py`; the bottom-left corner (index 0) is treated as part of the left side, as in
    GSFit's former writer.
    """
    import gsfit_rs

    n_r, n_z = len(grid_r), len(grid_z)
    d_r, d_z = float(np.mean(np.diff(grid_r))), float(np.mean(np.diff(grid_z)))
    edge_r, edge_z = rectangle_points(grid_r, grid_z)
    n_edge = len(edge_r)
    g_edge = gsfit_rs.greens_py(edge_r, edge_z, edge_r, edge_z, np.full(n_edge, d_r), np.full(n_edge, d_z))
    corners = (0, n_z - 1, n_z + n_r - 2, 2 * n_z + n_r - 3)
    sides = ((0, corners[1], 0.0, d_z), (corners[1] + 1, corners[2], d_r, 0.0),
             (corners[2] + 1, corners[3], 0.0, d_z), (corners[3] + 1, n_edge, d_r, 0.0))  # fmt: skip
    for start, stop, strip_r, strip_z in sides:
        for i in range(start, stop):
            g_edge[i, i] = self_inductance_rectangle(edge_r[i], strip_r, strip_z)
    return np.asarray(g_edge, dtype=np.float64) / (2.0 * np.pi)


# ----------------------------------------------------------------------------------------------
# Plasma mask
# ----------------------------------------------------------------------------------------------
def limiter_mask(
    r_flat: np.ndarray, z_flat: np.ndarray, d_r: float, d_z: float,
    vessel_r: np.ndarray, vessel_z: np.ndarray, limiter_r: np.ndarray, limiter_z: np.ndarray,
) -> np.ndarray:  # fmt: skip
    """Fixed allowed-plasma region, not the instantaneous flux-derived plasma mask.

    Include the vessel interior and points within half a grid diagonal of any limiter sample.
    For ST40 the vessel contour comes from ELMAG via GSFit's first wall unit; other points
    are excluded from carrying plasma current.
    """
    import shapely

    vessel = shapely.geometry.Polygon(np.column_stack((vessel_r, vessel_z)))
    mask = np.asarray(vessel.contains(shapely.points(r_flat, z_flat)), dtype=bool)
    half_diagonal = np.hypot(d_r, d_z) / 2.0
    distance = np.hypot(
        limiter_r[:, np.newaxis] - r_flat[np.newaxis, :], limiter_z[:, np.newaxis] - z_flat[np.newaxis, :]
    )
    return mask | np.any(distance < half_diagonal, axis=0)


# ----------------------------------------------------------------------------------------------
# Reader
# ----------------------------------------------------------------------------------------------
def _sensor_replacements(settings: Settings) -> list[SensorReplacement]:
    section = settings["gsfit"].get("sensor_replacements", {})
    return [
        SensorReplacement(name, tuple(item["replacements"]), tuple(float(c) for c in item["coefficients"]))
        for name, item in section.items()
    ]


@dataclass
class CoilCircuits:
    """RT-GSFit coils (`groups`: name -> GSFit PF coils in series) and the PCS signals that drive them."""

    groups: dict[str, list[str]]
    signals: list[str]
    matrix: np.ndarray  # [N_COIL, N_COIL_SIGNAL]


def coil_circuits(controller: Any, settings: Settings) -> CoilCircuits:
    """From `gsfit/coils`; by default every GSFit PF coil is a coil, driven by a signal of the same name."""
    section = settings["gsfit"].get("coils", {})
    if "groups" in section:
        groups = {name: list(members) for name, members in section["groups"].items()}
    else:
        groups = {name: [name] for name in controller.coils.keys(["pf"])}
    signals = list(section.get("signals", list(groups)))
    if "matrix" in section:
        matrix = np.asarray(section["matrix"], dtype=np.float64).reshape(len(groups), len(signals))
    elif len(signals) == len(groups):
        matrix = np.eye(len(groups))
    else:
        raise ValueError("gsfit/coils/matrix is needed when the number of signals differs from the number of coils")
    return CoilCircuits(groups, signals, matrix)


def preshot_inputs(controller: Any, settings: Settings, source: str) -> PreshotInputs:
    """Everything `assemble_preshot` needs, taken from a set-up GSFit controller."""
    from gsfit_rs.imas import equilibrium_paths as ep
    from gsfit_rs.imas import wall_paths as wp

    plasma, passives, wall = controller.plasma, controller.passives, controller.wall
    gsfit_settings = settings["gsfit"]

    grid_r = np.asarray(plasma.get(ep.time_slice[0].profiles_2d[0].grid.dim1), dtype=np.float64)
    grid_z = np.asarray(plasma.get(ep.time_slice[0].profiles_2d[0].grid.dim2), dtype=np.float64)
    n_r, n_z = len(grid_r), len(grid_z)
    n_grid = n_r * n_z
    d_r, d_z = float(np.mean(np.diff(grid_r))), float(np.mean(np.diff(grid_z)))
    r_2d, z_2d = np.meshgrid(grid_r, grid_z)
    r_flat, z_flat = r_2d.ravel(), z_2d.ravel()

    # PF coils: each RT-GSFit coil is a group of GSFit coils in series
    circuits = coil_circuits(controller, settings)
    groups = circuits.groups
    pf_active_names = list(plasma.get(ep.greens.pf_active[:].name))
    g_grid_coil = np.zeros((n_grid, len(groups)))
    for i_coil, members in enumerate(groups.values()):
        for member in members:
            g_grid_coil[:, i_coil] += np.ravel(plasma.get(ep.greens.pf_active[pf_active_names.index(member)].psi))

    # Passive degrees of freedom, in the order of GSFit's `greens.pf_passive`
    passive_names: list[str] = []
    passive_dof_units: list[str] = []
    g_grid_passive_columns: list[np.ndarray] = []
    passive_dofs: list[tuple[str, str]] = []
    for i_passive, passive in enumerate(passives.keys()):
        for i_dof, dof in enumerate(passives.keys([passive, "dof"])):
            passive_names.append(passive if dof == "constant_current_density" else f"{passive}_{dof.upper()}")
            passive_dof_units.append(_passive_dof_unit(dof))
            passive_dofs.append((passive, dof))
            g_grid_passive_columns.append(np.ravel(plasma.get(ep.greens.pf_passive[i_passive].dof[i_dof].psi)))
    n_passive = len(passive_names)
    g_grid_passive = np.column_stack(g_grid_passive_columns) if n_passive else np.zeros((n_grid, 0))

    # Constraints
    style = gsfit_settings.get("sensor_names", "st40_pcs")
    sensors = constraint_sensors(controller, settings)
    n_constraint = len(sensors)
    g_constraint_grid = np.zeros((n_constraint, n_grid))
    g_constraint_coil = np.zeros((n_constraint, len(groups)))
    g_constraint_passive = np.zeros((n_constraint, n_passive))
    constraint_weight = np.zeros(n_constraint)
    for i, sensor in enumerate(sensors):
        objects = _sensor_objects(controller, sensor.kind)
        g_constraint_grid[i] = objects.get_array1([sensor.name, "greens", "plasma"])
        for i_coil, members in enumerate(groups.values()):
            g_constraint_coil[i, i_coil] = sum(objects.get_f64([sensor.name, "greens", "pf", m]) for m in members)
        for i_dof, (passive, dof) in enumerate(passive_dofs):
            g_constraint_passive[i, i_dof] = objects.get_f64([sensor.name, "greens", "passives", passive, dof])
        weight = objects.get_f64([sensor.name, "fit_settings", "weight"])
        expected = objects.get_f64([sensor.name, "fit_settings", "expected_value"])
        # GSFit stores flux loop weights divided by 2 pi
        constraint_weight[i] = (2.0 * np.pi if sensor.kind == "flux_loop" else 1.0) * weight / expected

    # Regularisation of the passive degrees of freedom
    regularisation_rows: list[np.ndarray] = []
    regularisation_weight: list[float] = []
    regularisation_names: list[str] = []
    i_dof_start = 0
    for passive in passives.keys():
        n_dof = len(passives.keys([passive, "dof"]))
        rows = np.atleast_2d(passives.get_array2([passive, "regularisations"]))
        if rows.size:
            weights = passives.get_array1([passive, "regularisations_weight"])
            for i_reg, row in enumerate(rows):
                full_row = np.zeros(n_passive)
                full_row[i_dof_start : i_dof_start + n_dof] = row
                regularisation_rows.append(full_row)
                regularisation_weight.append(REGULARISATION_SCALE * float(weights[i_reg]))
                regularisation_names.append(f"{passive}_REG_{i_reg:02d}")
        i_dof_start += n_dof
    regularisation = np.vstack(regularisation_rows) if regularisation_rows else np.zeros((0, n_passive))

    # Wall: unit 0 is the vessel, every unit contributes limiter points
    description = wp.description_2d[0].limiter
    n_unit = len(wall.get(description.unit[:].name))
    limiter_r = np.concatenate([wall.get(description.unit[i].outline.r) for i in range(n_unit)]).astype(np.float64)
    limiter_z = np.concatenate([wall.get(description.unit[i].outline.z) for i in range(n_unit)]).astype(np.float64)
    mask_limiter = limiter_mask(
        r_flat, z_flat, d_r, d_z,
        np.asarray(wall.get(description.unit[0].outline.r)), np.asarray(wall.get(description.unit[0].outline.z)),
        limiter_r, limiter_z,
    )  # fmt: skip

    # Initial psi_norm and mask
    initial = gsfit_settings.get("initial", {"source": "gsfit_settings"})
    if initial["source"] == "gsfit_settings":
        conditions = controller.settings["rtgsfit_initial_conditions.json"]
        initial_psi_norm = np.asarray(conditions["flux_norm"], dtype=np.float64)
        initial_mask = np.asarray(conditions["mask"], dtype=np.int32)
        if initial_psi_norm.size != n_grid:
            raise ValueError(
                f"rtgsfit_initial_conditions.json has {initial_psi_norm.size} points but the grid has {n_grid}; "
                "use gsfit/initial/source = cosine_bump"
            )
    elif initial["source"] == "cosine_bump":
        initial_psi_norm = cosine_bump_psi_norm(grid_r, grid_z, initial["r"], initial["z"], initial["radius"])
        initial_mask = (initial_psi_norm < 1.0).astype(np.int32)
    else:
        raise ValueError(f"unknown gsfit/initial/source {initial['source']!r}")

    boundary = settings["preshot"]["boundary"]
    return PreshotInputs(
        mu_0=MU_0,
        grid_r=grid_r,
        grid_z=grid_z,
        poisson_operator="grad_shafranov",
        poisson_rhs_scale=2.0 * np.pi * r_flat * MU_0 * d_z**2,
        mask_limiter=mask_limiter,
        limiter_r=limiter_r,
        limiter_z=limiter_z,
        g_edge_edge=edge_greens(grid_r, grid_z),
        coil_names=list(groups),
        g_grid_coil=g_grid_coil,
        passive_names=passive_names,
        passive_dof_units=passive_dof_units,
        g_grid_passive=g_grid_passive,
        constraint_names=[sensor_name(sensor, style) for sensor in sensors],
        n_flux_loop=sum(sensor.kind == "flux_loop" for sensor in sensors),
        n_bp_probe=sum(sensor.kind == "bp_probe" for sensor in sensors),
        n_rogowski=sum(sensor.kind == "rogowski" for sensor in sensors),
        g_constraint_grid=g_constraint_grid,
        g_constraint_coil=g_constraint_coil,
        g_constraint_passive=g_constraint_passive,
        constraint_weight=constraint_weight,
        regularisation_names=regularisation_names,
        regularisation=regularisation,
        regularisation_weight=np.asarray(regularisation_weight, dtype=np.float64),
        sensor_replacements=_sensor_replacements(settings),
        coil_signals=circuits.signals,
        coil_matrix=circuits.matrix,
        initial_psi_norm=initial_psi_norm,
        initial_mask=initial_mask,
        xpt_fraction=boundary["xpt_fraction"],
        tolerance=boundary["tolerance"],
        ip_min=boundary["ip_min"],
        n_xpt_max=boundary["n_xpt_max"],
        settings_json=settings.json_text,
        source=source,
        psi_to_weber=1.0,
    )


def _experimental(objects: Any, path: list[str], time: np.ndarray) -> np.ndarray:
    value = np.asarray(objects.get_array1([*path, "experimental", "value"]), dtype=np.float64)
    value_time = np.asarray(objects.get_array1([*path, "experimental", "time"]), dtype=np.float64)
    return np.interp(time, value_time, value)


def gsfit_coil_current(controller: Any, circuits: Any, time: np.ndarray) -> np.ndarray:
    """[N_TIME, N_COIL]: GSFit's current in the first GSFit coil of each RT-GSFit coil."""
    return np.column_stack(
        [_experimental(controller.coils, ["pf", members[0], "i"], time) for members in circuits.groups.values()]
    )


def gsfit_signals(controller: Any, settings: Settings, preshot: Preshot) -> Signals:
    """The measurements GSFit read, interpolated onto the replay times.

    Coil currents are those of the first GSFit coil in each group. GSFit has coil currents, not the
    PCS coil signals, so `coil_signal` is the least-squares solution of
    `COIL:MATRIX coil_signal = coil_current`; coils whose row of the matrix is zero get no current,
    as in the real-time code.
    """
    time = replay_times(settings)
    style = settings["gsfit"].get("sensor_names", "st40_pcs")
    by_name = {sensor_name(sensor, style): sensor for sensor in constraint_sensors(controller, settings)}
    for kind, objects in (("flux_loop", controller.flux_loops), ("bp_probe", controller.bp_probes),
                          ("rogowski", controller.rogowski_coils)):  # fmt: skip
        for name in objects.keys():
            by_name.setdefault(sensor_name(Sensor(kind, name), style), Sensor(kind, name))

    meas = np.zeros((len(time), preshot.sizes["N_MEAS"]))
    for i, name in enumerate(preshot["MEAS:NAMES"]):
        if name not in by_name:
            raise ValueError(f"GSFit has no measurement {name!r}; sensor replacements need the st40_mdsplus reader")
        sensor = by_name[name]
        meas[:, i] = _experimental(
            _sensor_objects(controller, sensor.kind), [sensor.name, _QUANTITY[sensor.kind]], time
        )

    coil_current = gsfit_coil_current(controller, coil_circuits(controller, settings), time)
    coil_signal = coil_current @ np.linalg.pinv(preshot["COIL:MATRIX"]).T
    return Signals(time=time, meas=meas, coil_signal=coil_signal)


class GsfitReader(PreshotReader, SignalsReader):
    def read_preshot(self, settings: Settings) -> Preshot:
        controller = gsfit_controller(settings)
        source = f"gsfit {settings['gsfit']['settings_path']} pulse {controller.pulseNo}"
        return assemble_preshot(preshot_inputs(controller, settings, source))

    def read_signals(self, settings: Settings, preshot: Preshot) -> Signals:
        return gsfit_signals(gsfit_controller(settings), settings, preshot)
