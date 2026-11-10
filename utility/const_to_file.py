import os
import re
from scipy.io import loadmat
import numpy as np
import sys


def const_to_file(const_dict, file_read='constants.h', file_write='constants.c', as_enum=False):

    with open(file_read, 'r') as file:
        file_contents = file.read()

    # Define a function to perform variable replacements
    def replace_variables(match):

        assert((match.lastindex == 4) or (match.lastindex == 5))
        
        var_name = match.group(4).lower()
        var_name_in_file = " ".join([match.group(xx) for xx in range(2, 5)])
        if match.lastindex == 5:
            var_name_in_file = var_name_in_file + match.group(match.lastindex)
        
        if var_name in const_dict:
            replacement = const_dict[var_name]
            if as_enum:
                # Array sizes become enum constants so they are integer constant expressions
                return f'enum {{ {match.group(4)} = {int(replacement)} }};'
            # Check the data type of the replacement value
            if isinstance(replacement, bool):
                return f'{var_name_in_file} = {int(replacement)};'
            elif isinstance(replacement, (int, float)):
                return f'{var_name_in_file} = {replacement};'
            elif isinstance(replacement, np.ndarray) and (match.lastindex == 5):
                data_list = map(str, replacement.flatten().tolist())
                data_list = [x if ((ii+1) % 20) else x+'\n' for ii, x in enumerate(data_list)]
                array_values = ', '.join(data_list)
                return f'{var_name_in_file} = {{{array_values}}};'
            else:
                raise Exception(f'{type(replacement)} is not catered for')
        else:
            raise Exception(f"{var_name} is not in dictionary")
            
        return match.group(0)

    pattern = r'\b(\w+)\s*\b(\w+)\s*\b(\w+)\s*\b(\w+)(\[\s*\d*\s*\])?\s*;'
    
    # Replace variables in the C code
    new_c_code = re.sub(pattern, replace_variables, file_contents)

    with open(file_write, 'w') as file:
        file.write(new_c_code)
    
if __name__ == '__main__':

    if len(sys.argv) != 2:
        print("ERROR: python const_to_file.py <datafile>")
        sys.exit(1)
    dd = loadmat(sys.argv[1], squeeze_me=True)
    # Write both outputs to temporary files and only replace the real ones once both
    # succeed, so a failure never leaves a new constants.c paired with old sizes
    const_to_file(dd, file_write='constants.c.tmp')
    const_to_file(dd, file_read='constants_sizes.h.in', file_write='constants_sizes.h.tmp', as_enum=True)
    os.replace('constants_sizes.h.tmp', 'constants_sizes.h')
    os.replace('constants.c.tmp', 'constants.c')
    

