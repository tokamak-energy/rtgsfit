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

    cfg['REPO_PATH'] = repo_path
    cfg['DATA_DIR'] = os.path.join(repo_path, 'data')
    cfg['PLOTS_DIR'] = os.path.join(repo_path, 'plots')

    rtgsfit_path = os.path.dirname(os.path.dirname(repo_path))
    cfg['RTGSFIT_PATH'] = rtgsfit_path
    cfg['RTGSFIT_SRC_PATH'] = os.path.join(rtgsfit_path, 'src')

    if cfg['PULSE_NUM'] > 13000:
        cfg['PSU2COIL_RUN_NAME'] = "RUN05"
    else:
        cfg['PSU2COIL_RUN_NAME'] = "RUN01"

    return cfg
