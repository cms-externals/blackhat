#!/bin/bash
#PBS -N get_targets
#PBS -q route
#PBS -e log/
#PBS -o log/
##PBS -t 0-12
##PBS -t 0-9
#PBS -t 0-0
#PBS -l walltime=24:00:00,nodes=1:ppn=1

testdir=/home/hita/workspace/Harald_merge/test/

#
# get_targets.sh:
#
#  *) creates targets for test: 
#	- e.g. W2
#	- see list1 below
#  *) input: run get_targets first to generate
#	- c++ program: test/W2_PS.cpp
#	- PS points: test/PStest/W2/PSpoints.dat 
#	- process-data: test/PStest/W2/collectPS.log
#  *) output:
#	- target matrix elements stored in:
#	  	e.g. test/PStest/W2/target.dat
#



#############################################
# need to change also array request: (PBS -t 0-max)
#############################################
# generate all targets
# need array request above (PBS -t 0-12)
#list1=( 2j Z Wm Z1 Wm1 Z2 Wm2 3j Y1 Y2 Wm3 Y3 Z3 )
#############################################
# generate subset of targets
# need array request above (PBS -t 0-9)
#list1=( 2j Z Wm Z1 Wm1 Z2 Wm2 3j Y1 Y2 )
#############################################
# generate single target
# need array request above (PBS -t 0-0)
list1=( Wm3 )
#############################################

	echo " ======================================="
	label=${list1[$PBS_ARRAYID]}
	prog=${label}_PS
	src=${prog}.cpp

	sed -i 's/WRITE_TARGETS 1/WRITE_TARGETS 0/g' ${testdir}/${src}

	cd ${testdir}/../optimized/test/
	#pwd
	make ${prog} 2> /dev/null > /dev/null;
	make ${prog} 2> /dev/null > /dev/null;
	#remove taret file
	touch ${testdir}/PStest/${label}/target.dat 
	rm ${testdir}/PStest/${label}/target.dat 
	#write targets
	./${prog}  

	#cleanup
	sed -i 's/WRITE_TARGETS 0/WRITE_TARGETS 1/g' ${testdir}/${src}
	make ${prog} 2> /dev/null > /dev/null;

	echo " ======================================="
	echo ${list1[$PBS_ARRAYID]}
	echo "process: "$PBS_ARRAYID
	echo " ======================================="

