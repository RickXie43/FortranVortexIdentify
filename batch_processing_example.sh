#!/bin/bash
#Set INPUTDIR and will create new directory in OUTPUTDIR to save result
INPUTDIR=~/Data/beta_2023_pivdata/cone25_4K/230527_0.4Hz_4K_4s
OUTPUTDIR=.
#Paralleled Thread Number CORE_NUM
CORE_NUM=90
make clean
make DATADIR=$INPUTDIR
mpirun -np $CORE_NUM ./vortex_mpi_cmd $INPUTDIR $OUTPUTDIR
