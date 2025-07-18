#!/bin/bash
this_file_path="$(dirname "$(realpath "$0")")"
repo_path="$(realpath "$this_file_path/..")"
constants_c_path="$repo_path/data/constants.c"

# Generate constants.c
python -m rtgsfit_verify_analytic.generate_constants_c

# Change to RTGSFIT directory
cd $repo_path/../..

# Clean src directory
cd src
make clean

# Copy constants.c to src directory
cp $constants_c_path .

make SHOT=0 RUN_NAME="no_mds"

# Run replay_rtgsfit.py
python -m rtgsfit_verify_analytic.replay_rtgsfit