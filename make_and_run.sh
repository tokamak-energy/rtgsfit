# Set-up environment
# scl enable devtoolset-11 bash
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/
# set_global_conda
# conda activate py311

# Copy indexing files
cp /home/peter.buxton/tmp/rtgsfit_github_v2/data/rtgsfit_idx.txt data/rtgsfit_idx.txt
cp /home/peter.buxton/tmp/rtgsfit_github_v2/data/pfcoil_index.txt data/pfcoil_index.txt
cp /home/peter.buxton/tmp/rtgsfit_github_v2/data/meas.txt data/meas.txt
cp /home/peter.buxton/tmp/rtgsfit_github_v2/data/coil_curr.txt data/coil_curr.txt
cp /home/peter.buxton/tmp/rtgsfit_github_v2/data/flux_norm.txt data/flux_norm.txt
cp /home/peter.buxton/tmp/rtgsfit_github_v2/data/mask.txt data/mask.txt
cp /home/peter.buxton/tmp/rtgsfit_github_v2/data/psi_total.txt data/psi_total.txt

# Make rtgsfit source
rm src/*.o
rm src/constants.c
make -C src/ DATAFILE=/home/peter.buxton/0_Version_Controlled/rtgsfit/data/12001000_RUN04_for_c.mat

# Make the test program
rm tests/*.o
rm tests/replay_rtgsfit
make -C tests/ -f makefile_test PCS_PATH=/home/peter.buxton/0_Version_Controlled/pcs

# Run test program
(cd tests/ && ./replay_rtgsfit 12050 0.015 0.004 0.19)

