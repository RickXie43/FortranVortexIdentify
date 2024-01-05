#!/bin/bash
#Set INPUTDIR and will create new directory in OUTPUTDIR to save result
INPUTDIR=~/Data/beta_2023_pivdata/cone25_4K/230527_0.4Hz_4K_4s
OUTPUTDIR=.
#Paralleled Thread Number CORE_NUM
CORE_NUM=90
make clean
make DATADIR=$INPUTDIR
mpirun -np $CORE_NUM ./vortex_mpi_cmd $INPUTDIR $OUTPUTDIR


#if you want to batch processing a series of INPUTDIR, use this:
#Paralleled Thread Number CORE_NUM
CORE_NUM=100
#Set INPUTDIR and will create new directory in OUTPUTDIR to save result
#dir processing
INPUTDIR=~/DataDIRS
OUTPUTDIR=~/ResultsDIR
SUBDIRS=("$INPUTDIR"/*)
for dir in "${SUBDIRS[@]}";
do
	if [[ -d "$dir" ]]; then
		echo "Processing $dir"
		make clean
		make DATADIR="$dir"
		mpirun -np $CORE_NUM ./vortex_mpi_cmd "$dir" $OUTPUTDIR
		echo "Finished processing $dir."
	fi
done

