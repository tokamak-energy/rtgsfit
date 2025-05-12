#!/bin/bash

# Define the pulse prefix
pulse_prefix=52000000
run_name="SCAN0"

# Define the range of pulse numbers and exclude 13344
pulse_nums=()
for i in {13342..13349}; do
    if [ "$i" -ne 13344 ]; then
        pulse_nums+=($((pulse_prefix + i)))
    fi
done

# Iterate over the pulse numbers and run the Python script
for pulse_num in "${pulse_nums[@]}"; do
    python -c "
import standard_utility as util
print('Creating MDSPlus node for pulse number:', $pulse_num)
util.create_script_nodes(
    script_name='RTGSFIT',
    pulseNo_write=$pulse_num,
    run_name='$run_name',
)
"
done