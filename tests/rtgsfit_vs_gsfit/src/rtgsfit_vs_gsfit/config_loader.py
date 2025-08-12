import os
import json

def load_and_prepare_config() -> dict:
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
