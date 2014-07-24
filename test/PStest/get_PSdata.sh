#!/bin/bash
#PBS -N get_PSdata
#PBS -q route
#PBS -e log/
#PBS -o log/
##PBS -t 0-12
#PBS -t 0-0
#PBS -l walltime=24:00:00,nodes=1:ppn=1


#
# get_PSdata.sh
#
#  *) creates test: e.g. W2
#	- e.g. W2
#	- see list1 below
#  *) input:
#	- data-dir: mkdir test/PStest/W2 
#	- test/PStest/W2/Run.dat file
#	- sherpa path: see SHERPAdir variable	
#  *) output: 
#	- c++ program: test/W2_PS.cpp
#  	- process-data: test/PStest/W2/collectPS.log 
#  	- PS points: test/PStest/W2/PSpoints.dat
#



##############################
#
# need to provide Run.dat files
# e.g.: test/PStest/WW2j/
# 	test/PStest/process_name/Run.dat
#
# sherpa path
SHERPAdir=/home/hita/workspace/SHERPA-MPI_MAR29/
# BH path
BHdir=/home/hita/workspace/Harald_merge/
#
##############################



#############################################
# need to change also array request: (PBS -t 0-max)
#############################################
# generate all targets
# need array request above (PBS -t 0-12)
#list1=( 2j 3j Wm Wm1 Wm2 Wm3 Z Z1 Z2 Z3 Y1 Y2 Y3 )
# needed estimated timed (sec) for grid construction
#list2=( 100 300 100 100 200 600 100 100 200 1000 100 100 600 )
#############################################
# generate subset of targets
# need array request above (PBS -t 0-9)
#list1=( 2j Z Wm Z1 Wm1 Z2 Wm2 3j Y1 Y2 )
# needed estimated timed (sec) for grid construction
#list2=( 100 100 100 100 100 200 200 300 )
#############################################
# generate single target
# need array request above (PBS -t 0-0)
list1=( Wm3 )
list2=( 600 )
#############################################


cd ${BHdir}/test/PStest/

	label=${list1[$PBS_ARRAYID]}
	prog=${label}_PS
	src=${prog}.cpp
	grid_const_time=${list2[$PBS_ARRAYID]}

	#prepare check program
	cp data/check_program.cpp ../${src}
	eval "sed -i s/STUDYNAME/${label}/g ../${src} "

	#generate PS points	
	cd ${label}/
	# generate_PS.sh STUDYNAME BLACKHATDIR SHERPADIR
	../data/generate_PS.sh ${label} ${BHdir} ${SHERPAdir} ${grid_const_time}


	echo " ======================================="
	echo ${list1[$PBS_ARRAYID]}
	echo "process: "$PBS_ARRAYID
	echo " ======================================="
	
	cd ..

