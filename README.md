# rtgsfit
Real-Time equilibrium reconstruction code

## Compilation
A .mat datafile will be required with the same variable names of that of the 
global constants specified in constants.h. This should 
contain all matrices in row major order, with indexing also in row major order 
and starting from 0.  The conda environment is only required for running the 
python tests

```bash
conda create -n rtgsfit python=3.9
pip install -r requirements.txt
cd src
make DATAFILE=<PATH/TO/DATAFILE.mat>
cd ../tests
make
```




## Program Structure
* shared libraries
* dependancy graph

rtgsfit

lapacke.h
cblas.h
math.h 
float.h
stdio.h
time.h
string.h

## Formatting convention
* Allman bracket style
* 

## Naming convention
* Global constants are fully capitilised with underscores e.g. R_GRID
* All global variables are constant and defined at compile time
* variables and function names should be written in snake case 
i.e. lower_case_with_underscores
* variable names should be a noun and should refer to the quantity rather than 
a Greek letter. e.g. flux instead or psi.
* function names should begin with a verb e.g. 'get_flux_norm' 
* boolean variables are simply integers in C. Maybe use short or other inbuilt types
* Have not used size_t or any other non-inbuilt types
* Boolean names should begin with is_ prefix to differentiate them from other variables
* Boolean names should be always be in the "positive" affirmation, with 1 being
true and 0 being false
* order of variables should go from largest to smallest e.g. GRID_R_MIN
as if it was sucessive objects in the OOP paradigme.  
* Green's functions are prefixed with g_

Names that could be converted to structures



## Function conventions
* docstring format
* order of inputs outputs
* example docstring below ...


## Grid convention
R_GRID[0, 0], Z_GRID[0, 0] - bottom left of grid i.e. (R_MIN, Z_MIN)
R_GRID[0, N_R], Z_GRID[0, N_R] - bottom right of grid i.e. (R_MAX, Z_MIN)
R_GRID[N_Z, 0], Z_GRID[N_Z, 0] - top left of grid i.e. (R_MIN, Z_MAX)
R_GRID[N_Z, N_R], Z_GRID[N_Z, N_R] - top right of grid i.e. (R_MAX, Z_MAX)

[[ (N_Z - 1, 0), (N_Z - 1, 1)  , ... , (N_Z- 1, N_R -1)  ],
 [ (N_Z - 2, 0), (N_Z - 2, 1)  , ... , (N_Z - 2, N_R -1)],
                    ...
 [ (0, 0)      , (0, 1)       , ... , (0, N_R -1)       ]]
 
## Boundary Convention 
* LTRB
* (R_MIN, Z_MIN) -> (R_MIN, Z_MAX) -> (R_MAX, Z_MAX) -> (R_MAX, Z_MIN) -> (R_MIN, Z_MIN)


 
## Testing
* pytest
* seting up python
* freegs

## To Do
* deglss vs dgelsd
* python flux testing
* interpolate hess_rr hess_det
* boundary index convention
* timing
* improve comments code
* write matlab wrapper for code 
* create matlab function that takes pulse and produces the .mat file
* vessel filaments


# To compile
(Peter's notes)
```bash
git checkout replay_rtgsfit
scl enable devtoolset-11 bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/
cd src/
make
```

To check out PCS
/home/peter.buxton/0_Version_Controlled/pcs/model/ST40PCS
gcc -c -o bin/mds_tools.o src/mds_tools.c -Iinclude -I/usr/local/mdsplus/include
gcc -c -o bin/utils.o src/utils.c -Iinclude

To compile tests
```bash
cd ../tests/
make -f makefile_test PCS_PATH=/home/peter.buxton/0_Version_Controlled/pcs
```

./replay_rtgsfit 12050 0.01 0.002 0.2


# Important
"const_to_file.py" is what adds the constants into the c



make DATAFILE=/home/peter.buxton/0_Version_Controlled/rtgsfit/data/12001000_RUN04_for_c.mat


One line running:
```bash
cd src/; make DATAFILE=/home/peter.buxton/0_Version_Controlled/rtgsfit/data/12001000_RUN04_for_c.mat; cd ../tests/; rm replay_rtgsfit; make -f makefile_test PCS_PATH=/home/peter.buxton/0_Version_Controlled/pcs; ./replay_rtgsfit 12050 0.01 0.0004 0.2; cd ../
```

plotting
```bash
python3 ../py-files/plot_timed_data.py
````


backtrace
gdb --args ./replay_rtgsfit 12050 0.02 0.0005 0.2;
run
bt

# Instructions for aleksei

source /opt/intel/oneapi/setvars.sh