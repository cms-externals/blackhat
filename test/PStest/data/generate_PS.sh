#!/bin/bash

####################
# only edit from here
#testname=STUDYNAME
#BHdir=BLACKHATDIR
#sherparundir=SHERPARUNDIR
testname=$1
BHdir=$2
sherparundir=$3
gridtime=$4
# only edit until here
####################

echo "*******************"
echo "using data"
echo $testname
echo $BHdir
echo $sherparundir
echo $gridtime
echo "*******************"


locdir=${BHdir}/test/PStest/${testname}
cd $locdir

#clean beofre we start
touch {BHsettings*,collectPS.log,PSpoints.dat,Process,r0}
rm -fr {BHsettings*,collectPS.log,PSpoints.dat,Process,r0}


#lib and grid generation
cp ../data/BHsettings_grid BHsettings
$sherparundir/bin/Sherpa
./makelibs
$sherparundir/bin/Sherpa &

#generate grid
pid=$! 
sleep $gridtime;
if ps -p $pid >/dev/null; then kill -9 $pid; fi

#integration
cp ../data/BHsettings_collect BHsettings
$sherparundir/bin/Sherpa

rm -fr {BHsettings*,Process,r0,Analysis,makelibs,Sherpa_References.tex}



