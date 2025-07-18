"""
Generate a dictionary by reading the values from the
constants.c file.
"""

import re

import numpy as np

def constants_c_dict(constants_c_path,
                     read_g_grid_meas_weight = False):
    """
    Read the constants.c file to extract some of the key constants used in RTGSFIT.
    """

    with open(constants_c_path, 'r') as file:
        content = file.read()

    c_dict = {}

    # Extract N_R and N_Z
    c_dict["n_r"] = int(re.search(r'const int N_R\s*=\s*(\d+);', content).group(1))
    c_dict["n_z"] = int(re.search(r'const int N_Z\s*=\s*(\d+);', content).group(1))

    # Extract N_MEAS
    c_dict["n_meas"] = int(re.search(r'const int N_MEAS\s*=\s*(\d+);', content).group(1))

    # Extact N_COIL
    c_dict["n_coil"] = int(re.search(r'const int N_COIL\s*=\s*(\d+);', content).group(1))

    # Extact N_LCFS_MAX
    c_dict["n_lcfs_max"] = int(re.search(r'const int N_LCFS_MAX\s*=\s*(\d+);', content).group(1))

    # Extract N_COEF
    c_dict["n_coef"] = int(re.search(r'const int N_COEF\s*=\s*(\d+);', content).group(1))

    # Extract R_VEC array
    r_vec_match = re.search(r'const double R_VEC\[\]\s*=\s*{([^}]*)};', content, re.DOTALL)
    r_vec_str = r_vec_match.group(1)
    r_vec_list = list(map(float, r_vec_str.replace('\n', '').split(',')))
    c_dict["r_vec"] = np.array(r_vec_list)

    # Extract Z_VEC array
    z_vec_match = re.search(r'const double Z_VEC\[\]\s*=\s*{([^}]*)};', content, re.DOTALL)
    z_vec_str = z_vec_match.group(1)
    z_vec_list = list(map(float, z_vec_str.replace('\n', '').split(',')))
    c_dict["z_vec"] = np.array(z_vec_list)

    # Extract N_LIMIT
    n_limit_match = re.search(r'const int N_LIMIT\s*=\s*(\d+);', content)
    c_dict["n_limit"] = int(n_limit_match.group(1))

    # Extract LIMIT_IDX array
    limit_idx_match = re.search(r'const int LIMIT_IDX\[\]\s*=\s*{([^}]*)};', content, re.DOTALL)
    limit_idx_str = limit_idx_match.group(1)
    limit_idx_list = list(map(int, limit_idx_str.replace('\n', '').split(',')))
    c_dict["limit_idx"] = np.array(limit_idx_list) 

    # Extract LIMIT_WEIGHT array
    limit_weight_match = re.search(r'const double LIMIT_WEIGHT\[\]\s*=\s*{([^}]*)};', content, re.DOTALL)
    limit_weight_str = limit_weight_match.group(1)
    limit_weight_list = list(map(float, limit_weight_str.replace('\n', '').split(',')))
    c_dict["limit_weight"] = np.array(limit_weight_list)

    if read_g_grid_meas_weight:
        # Extract G_GRID_MEAS_WEIGHT array
        g_grid_meas_weight_match = re.search(r'const double G_GRID_MEAS_WEIGHT\[\]\s*=\s*{([^}]*)};', content, re.DOTALL)
        g_grid_meas_weight_str = g_grid_meas_weight_match.group(1)
        g_grid_meas_weight_list = list(map(float, g_grid_meas_weight_str.replace('\n', '').split(',')))
        c_dict["g_grid_meas_weight"] = np.array(g_grid_meas_weight_list)

    return c_dict