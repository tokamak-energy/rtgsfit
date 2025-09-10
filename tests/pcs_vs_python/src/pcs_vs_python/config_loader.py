import os

import yaml

def load_and_prepare_config(run_name: str = None,
                            pulse_num: int = None) -> dict:
    """
    Load the base configuration from a YAML file and augment it with derived values.

    This function reads the default configuration from the `data/default_config.yaml` file
    in the repository, then computes and adds additional path-related and runtime fields
    based on the repository location and configuration values.

    The resulting configuration dictionary includes:
      - Original values from the YAML file.
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

    config_path = os.path.join(repo_path, 'data', 'default_config.yaml')
    with open(config_path, 'r') as f:
        cfg = yaml.safe_load(f)

    cfg["n_t"] = int((cfg["t_max_approx"] - cfg["t_min"]) / cfg["d_t"]) + 1
    cfg["t_max"] = cfg["t_min"] + (cfg["n_t"] - 1) * cfg["d_t"]

    cfg['repo_path'] = repo_path
    cfg['rtgsfit_src_path'] = os.path.join(cfg['rtgsfit_path'], 'src')
    cfg["plot_dir"] = os.path.join(repo_path, 'plots')

    return cfg

if __name__ == "__main__":
    cfg = load_and_prepare_config()
    for key in cfg.keys():
        print(f"{key}: {cfg[key]}")
