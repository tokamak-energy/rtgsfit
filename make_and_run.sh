# Set-up environment
# scl enable devtoolset-11 bash
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/
# set_global_conda
# conda activate py311

# Copy indexing files
cp files_for_rtgsfit/rtgsfit_idx.txt data/rtgsfit_idx.txt
cp files_for_rtgsfit/pfcoil_index.txt data/pfcoil_index.txt
cp files_for_rtgsfit/meas.txt data/meas.txt
cp files_for_rtgsfit/coil_curr.txt data/coil_curr.txt
cp files_for_rtgsfit/flux_norm.txt data/flux_norm.txt
cp files_for_rtgsfit/mask.txt data/mask.txt
cp files_for_rtgsfit/psi_total.txt data/psi_total.txt

# Make rtgsfit source
rm src/*.o
rm src/constants.c
rm lib/*.so
# make -C src/ DATAFILE=/home/peter.buxton/0_Version_Controlled/rtgsfit/data/12001000_RUN05_for_c.mat
make -C src/ SHOT=52012050 RUN_NAME=GJ_EIG12345

# Make the test program
rm tests/*.o
rm tests/replay_rtgsfit
make -C tests/ -f makefile_test PCS_PATH=/home/alex.prokopyszyn/GitLab/pcs

# Run test program
# (cd tests/ && ./replay_rtgsfit 12050 0.015 0.001 0.19)

