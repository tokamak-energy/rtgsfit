import pandas as pd

sens_file = '../data/sensors_to_use.csv'
sens = pd.read_csv(sens_file)
sens_name = sens['name'].tolist()
sens_inc = sens['include'].tolist()

with open('../data/included_sensors.txt', 'w') as f:
    sens_use = [];
    for i, si in enumerate(sens_inc):
        if si == 1:
    	    f.write(sens_name[i] + '\n')
f.close()
