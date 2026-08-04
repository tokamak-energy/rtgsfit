import json
import os
import re

import mdsthin
import numpy as np


def next_test_run_name(pulse_num_write: int) -> str:
    """
    Allocate the next TEST#### run name by inspecting existing RTGSFIT runs.
    """
    existing_run_names = []
    with mdsthin.Connection('smaug') as conn:
        conn.openTree("RTGSFIT", pulse_num_write)
        existing_run_names = np.atleast_1d(
            conn.get('getnci("\\\\RTGSFIT::TOP.*", "node_name")').data()
        )

    latest_index = 0
    for raw_name in existing_run_names:
        name = str(raw_name).strip().upper()
        match = re.fullmatch(r"TEST(\d{4})", name)
        if match is not None:
            latest_index = max(latest_index, int(match.group(1)))
    return f"TEST{latest_index + 1:04d}"

def load_and_prepare_config(run_name: str = None,
                            pulse_num: int = None) -> dict:
    """
    Load the base configuration from a JSON file and augment it with derived values.

    This function reads the default configuration from the `data/default_config.json` file
    in the repository, then computes and adds additional path-related and runtime fields
    based on the repository location and configuration values.

    The resulting configuration dictionary includes:
      - Original values from the JSON file.
      - Paths to key directories and source code.
      - Automatically selected run names based on `pulse_num`.
      - Flags for RTGSFIT and GSFIT runtime state.

    Returns
    -------
    dict
        A configuration dictionary containing both loaded and computed values.
    """

    this_file_path = os.path.abspath(__file__)
    repo_path = os.path.dirname(os.path.dirname(os.path.dirname(this_file_path)))
    rtgsfit_path = os.path.dirname(os.path.dirname(repo_path))

    config_path = os.path.join(repo_path, 'data', 'default_config.json')
    with open(config_path, 'r') as f:
        cfg = json.load(f)

    if pulse_num is not None:
        cfg['pulse_num'] = pulse_num
        cfg["pulse_num_write"] = 52_000_000 + pulse_num
    if run_name is not None:
        cfg['run_name'] = run_name

    cfg['repo_path'] = repo_path
    cfg['data_dir'] = os.path.join(repo_path, 'data')
    cfg['plots_dir'] = os.path.join(repo_path, 'plots')
    cfg["plots_this_run_dir"] = os.path.join(cfg["plots_dir"], f'{cfg["pulse_num"]}_{cfg["run_name"]}')
    cfg['rtgsfit_path'] = rtgsfit_path
    cfg['rtgsfit_src_path'] = os.path.join(rtgsfit_path, 'src')
    cfg["rtgsfit_output_dict_path"] = \
        os.path.join(cfg["data_dir"],
                     f'rtgsfit_output_dict_{cfg["pulse_num"]}_{cfg["run_name"]}.npy')
    cfg["ivc_dict_path"] = os.path.join(cfg["data_dir"],
                                        f'ivc_dict_{cfg["pulse_num"]}_{cfg["run_name"]}.npy')
    cfg["coef_names_path"] = os.path.join(cfg["data_dir"],
                                          f'coef_names_{cfg["pulse_num"]}_{cfg["run_name"]}.txt')
    cfg["meas_names_path"] = os.path.join(cfg["data_dir"],
                                          f'meas_names_{cfg["pulse_num"]}_{cfg["run_name"]}.txt')
    cfg["gsfit_pf_coils_path"] = os.path.join(cfg["data_dir"],
                                              f'gsfit_pf_coil_currents_{cfg["pulse_num"]}_{cfg["run_name"]}.json')

    if cfg['pulse_num'] > 13000:
        cfg['psu2coil_run_name'] = "run05"
    else:
        cfg['psu2coil_run_name'] = "run01"

    cfg["gsfit_replayed"] = False
    cfg["rtgsfit_node_initialised"] = False
    cfg["rtgsfit_compiled"] = False
    cfg["rtgsfit_replayed"] = False

    cfg["j_vrange"] = None

    return cfg
