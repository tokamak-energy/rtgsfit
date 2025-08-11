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
    cfg['rtgsfit_path'] = rtgsfit_path
    cfg['rtgsfit_src_path'] = os.path.join(rtgsfit_path, 'src')

    if cfg['pulse_num'] > 13000:
        cfg['psu2coil_run_name'] = "run05"
    else:
        cfg['psu2coil_run_name'] = "run01"

    cfg["gsfit_replayed"] = False
    cfg["rtgsfit_node_initialised"] = False
    cfg["rtgsfit_compiled"] = False
    cfg["rtgsfit_replayed"] = False

    return cfg
