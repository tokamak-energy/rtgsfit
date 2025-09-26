#!/bin/bash

### Initialisation ###

# Define the pulse prefix
pulse_prefix=52000000
run_name="SCAN3"
run_description="Part of a scan of pulses to test the RT-GSFit code. \
Starting at a later time, namely, 60e-3 s. \
Also changed MC current from equal to MCT to average of MCB and MCT. \
Now using pulse number 99000230 instead of 12001000 in the RT-GSFit Matlab setup code \
also using elmag_run=RUN18 instead of RUN16. \
Where the aim is to fix the issue of assuming the MCB and MCT currents are equal \
when they are actually of opposite sign in the latter part of the pulse."
# Define the range of pulse numbers and exclude 13344
pulse_nums_write=()
pulse_nums=()
for i in {13342..13349}; do
    if [ "$i" -ne 13344 ]; then
        pulse_nums_write+=($((pulse_prefix + i)))
        pulse_nums+=($i)
    fi
done
# Define the path to the Matlab code
matlab_code_path="/home/alex.prokopyszyn/GitLab/Forks/rtgsfit_setup"

### Exceution ###

# Iterate over the pulse numbers and setup the MDSPlus nodes
for pulse_num_write in "${pulse_nums_write[@]}"; do
    python -c "
import standard_utility as util
print('Creating MDSPlus node for pulse number:', $pulse_num_write )
util.create_script_nodes(
    script_name='RTGSFIT',
    pulseNo_write=$pulse_num_write,
    run_name='$run_name',
)
"
done

# Run RT-GSFit Matlab code for each pulse number
echo "Changing directory to $matlab_code_path"
cd "$matlab_code_path" || { echo "Failed to change directory to $matlab_code_path"; exit 1; }
echo "Now in directory: $(pwd)"
# Modify tests/setup_main.m for each pulse and run the RT-GSFit Matlab code for each pulse number
for pulse_num_write in "${pulse_nums_write[@]}"; do
    # Search for the line containing $search_tag and replace it with the new line
    search_tag="% sed_tag:pulse_write_replay_range_of_pulses"
    line_replacement="pulse_write = $pulse_num_write ; $search_tag"
    sed -i "/$search_tag/c\\$line_replacement" tests/setup_main.m
    search_tag="% sed_tag:run_name_replay_range_of_pulses"
    line_replacement="path_to_node = '\\\RTGSFIT::TOP.$run_name.PRESHOT'; $search_tag"
    sed -i "/$search_tag/c\\$line_replacement" tests/setup_main.m
    matlab -nodisplay -nosplash -nodesktop -r "run('tests/setup_main.m');exit;"
done

# Run the replay python script for each pulse number
# Modify make_and_run.sh for each pulse number
# Change the directory to the location the location of the shell script
dirname="$(dirname "$(realpath "$0")")"
cd "$dirname" || { echo "Failed to change directory to $dirname"; exit 1; }
echo "Changed directory to $(pwd)"
# Go up a directory
cd .. || { echo "Failed to change directory to the parent directory"; exit 1; }
echo "Changed directory to $(pwd)"
for pulse_num in "${pulse_nums[@]}"; do
    pulse_num_write=$((pulse_prefix + pulse_num))
    # Make for the current pulse number
    search_tag="# sed_tag:pulse_num_run_name_replay_range_of_pulses"
    line_replacement="make -C src/ SHOT=$pulse_num_write RUN_NAME=$run_name $search_tag"
    sed -i "/$search_tag/c\\$line_replacement" make_and_run.sh
    # Run the modified make_and_run.sh script
    bash make_and_run.sh

    # Modify python script for the current pulse number
    search_tag="# sed_tag:pulse_num_replay_range_of_pulses"
    line_replacement="pulseNo = $pulse_num $search_tag"
    sed -i "/$search_tag/c\\$line_replacement" py-files/replay_rtgsfit_from_python.py
    search_tag="# sed_tag:pulse_num_write_replay_range_of_pulses"
    line_replacement="pulseNo_write = $pulse_num_write $search_tag"
    sed -i "/$search_tag/c\\$line_replacement" py-files/replay_rtgsfit_from_python.py
    search_tag="# sed_tag:run_name_replay_range_of_pulses"
    line_replacement="run_name = '$run_name' $search_tag"
    sed -i "/$search_tag/c\\$line_replacement" py-files/replay_rtgsfit_from_python.py
    search_tag="# sed_tag:run_description_replay_range_of_pulses"
    line_replacement="run_description = '$run_description' $search_tag"
    sed -i "/$search_tag/c\\$line_replacement" py-files/replay_rtgsfit_from_python.py

    # Run the modified python script
    python py-files/replay_rtgsfit_from_python.py

done

echo "Done"