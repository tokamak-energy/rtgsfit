"""
Export the analytic-test inputs (raw measurements, coil currents, initial
normalised flux and mask) as text files for tests/rtgsfit_bench.c.

Run after generate_constants_c (which writes data/constants.c, data/flux_loops.txt
and data/bp_probes.txt), e.g. after `pytest`:

    python -m rtgsfit_verify_analytic.export_bench_inputs [out_dir]

The default output directory is data/bench/. A second set with zero plasma
current is written to data/bench_zero/.
"""
import os
import sys

import numpy as np

from rtgsfit_verify_analytic import cnst, measurements, read_constants_c, replay_rtgsfit


def export_bench_inputs(out_dir=None):
    if out_dir is None:
        out_dir = os.path.join(cnst.DATA_DIR, "bench")
    constants_c_path = os.path.join(cnst.DATA_DIR, "constants.c")
    c_dict = read_constants_c.constants_c_dict(constants_c_path)
    data = np.loadtxt(os.path.join(cnst.DATA_DIR, "flux_loops.txt"), dtype=np.float64, skiprows=1)
    fl_coords = data[:, :2].copy()
    data = np.loadtxt(os.path.join(cnst.DATA_DIR, "bp_probes.txt"), dtype=np.float64, skiprows=1)
    bp_probe_coords = data[:, :3].copy()

    cases = [(out_dir, cnst.ANALYTIC_PLASMA_CURRENT), (out_dir + "_zero", 0.0)]
    for case_dir, plasma_current in cases:
        os.makedirs(case_dir, exist_ok=True)
        meas = measurements.generate_measurements(
            fl_coords, bp_probe_coords, cnst.ANALYTIC_RO, cnst.ANALYTIC_ZO, plasma_current
        )
        np.savetxt(os.path.join(case_dir, "meas_pcs.txt"), meas, fmt="%.17g")
        np.savetxt(os.path.join(case_dir, "coil_curr.txt"), np.zeros(c_dict["n_coil"]), fmt="%.17g")
        np.savetxt(
            os.path.join(case_dir, "flux_norm.txt"),
            replay_rtgsfit.initial_flux_norm(c_dict["r_vec"], c_dict["z_vec"]),
            fmt="%.17g",
        )
        np.savetxt(
            os.path.join(case_dir, "mask.txt"),
            np.ones(c_dict["n_r"] * c_dict["n_z"], dtype=np.int32),
            fmt="%d",
        )
        print(f"wrote {case_dir}: n_meas={len(meas)}, n_grid={c_dict['n_r'] * c_dict['n_z']}, "
              f"plasma_current={plasma_current}")


if __name__ == "__main__":
    export_bench_inputs(sys.argv[1] if len(sys.argv) > 1 else None)
