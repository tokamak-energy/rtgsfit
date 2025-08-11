import os
import json

def load_config():
    """
    Load configuration from a JSON file to create a dictionary.

    The dictionary contains information that will be used to set and and run RTGSFIT and GSFIT.

    The JSON file is expected to be in the data directory of the repository.
    """

    this_file_path = os.path.abspath(__file__)
    repo_path = os.path.dirname(os.path.dirname(os.path.dirname(this_file_path)))

    config_path = os.path.join(repo_path, 'data', 'default_config.json')
    with open(config_path, 'r') as f:
        cfg = json.load(f)

    cfg['repo_path'] = repo_path
    cfg['data_dir'] = os.path.join(repo_path, 'data')
    cfg['plots_dir'] = os.path.join(repo_path, 'plots')

    rtgsfit_path = os.path.dirname(os.path.dirname(repo_path))
    cfg['rtgsfit_path'] = rtgsfit_path
    cfg['rtgsfit_src_path'] = os.path.join(rtgsfit_path, 'src')

    if cfg['pulse_num'] > 13000:
        cfg['psu2coil_run_name'] = "run05"
    else:
        cfg['psu2coil_run_name'] = "run01"

    return cfg
