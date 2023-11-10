# FortranVortexIdentify
This Package aims to identify discrete vortexs from the velocity field data, which is the result of DaVis PIV process.

## Quick Start
Run these in terminal. Set INPUTDIR, OUTPUTDIR and CORE\_NUM for your data.

```
#!/bin/bash
#Set INPUTDIR and will create new directory in OUTPUTDIR to save result
INPUTDIR=~/Data/beta_2023_pivdata/cone25_4K/230527_0.4Hz_4K_4s
OUTPUTDIR=.
#Paralleled Thread Number CORE_NUM
CORE_NUM=90
make clean
make DATADIR=$INPUTDIR
mpirun -np $CORE_NUM ./vortex_mpi_cmd $INPUTDIR $OUTPUTDIR
```

## Usage
To compute all vortexs, there are three steps. First, you have to set parameters
in *parameters.F9*0. Second, you have to compile the program. And the third is to
processing the program.

Since for different original PIV data, different NX\_DEF, NY\_DEF (
        point number in x, y direction) should be set before compiling. The batch processing
method is also provided in order to set NX\_DEF and NY\_DEF automatically, using shell
script *setnxny.s*h .

### Set Parameters
- X\_MIN\_DEF, X\_MAX\_DEF, Y\_MIN\_DEF, Y\_MAX\_DEF

These 4 parameters determine the range of field, a square area. The Unit is meter. For example, X\_MIN\_DEF -0.1
- INTERPOLATIONPOINTS\_DEF

The points number between X\_MIN\_DEF and X\_MAX\_DEF, and also between Y\_MIN\_DEF and Y\_MAX\_DEF. The larger INTERPOLATIONPOINTS\_DEF, the more accurate for vortexs position. For example, INTERPOLATIONPOINTS\_DEF 2000
- NX\_DEF, NY\_DEF

The points number in x and y direction of the original PIV data. The shell script **setnxny.sh* can read NX\_DEF and NY\_DEF from the first line of data, and write into *parameters.F90*.
You can also set by yourself.
- EFF\_RADIUS\_DEF

The vortexs out of these range will be deleted. e.g. EFF\_RADIUS\_DEF 0.1
- FLITER\_BETA\_DEF

Delete a ratio of vortexs. Larger to have fewer vortewx. e.g. FLITER\_BETA\_DEF 1

### Compile Program
The program uses Makefile to compile.

1.
To compile using parameters in *parameters.F90*, run 
```
make clean
make
```
in terminal.

2.
To set NX\_DEF and NY\_DEF for **specified data**, run
```
make clean
make DATADIR=/directory
```
in terminal, which can set NX\_DEF and NY\_DEF in *parameters.F90* first and then compiling.

3.
To run **single processor version program** using *main_findvortex.f90*, run
```
make clean
make single_processor
```
in terminal.

4.
To **test program** using *main_test.f90*, run
```
make clean
make test
```
in terminal.

### Batch Processing
The example is shown in *./batch_processing_example.sh*.

```
#!/bin/bash
#Set INPUTDIR and will create new directory in OUTPUTDIR to save result
INPUTDIR=~/Data/beta_2023_pivdata/cone25_4K/230527_0.4Hz_4K_4s
OUTPUTDIR=.
#Paralleled Thread Number CORE_NUM
CORE_NUM=90
make clean
make DATADIR=$INPUTDIR
mpirun -np $CORE_NUM ./vortex_mpi_cmd $INPUTDIR $OUTPUTDIR
```

run in terminal.
