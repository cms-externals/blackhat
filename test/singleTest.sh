t=$1
wd=$2

cd $wd
mkdir -p $wd/status
mkdir -p $wd/logs
touch status/$t.started

date=$(date +%d%m%y_%H%M%S)
make $t
if test -f $wd/.libs/$t ; then
	executable=$wd/.libs/$t
else 
	executable=$wd/$t
fi
$executable 2>&1 $wd/logs/$t.$date
ES=$?

echo $ES > status/$t.exit


 
