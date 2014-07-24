#!/bin/env python

import time
import os
import re
import stat
import shlex, subprocess
import sys


#############################
# start input section

#shorttests=['Y2','Y3']
#shorttests=['W2','W3','Y2','Y3']
#longtests=['W4','W5','Y4']
shorttests=['W2','W3']
longtests=['W4','W5']

#  vector of entries [ [ processes, part, number of PS-points,  precision, \
#                        mass^2=s(c,c+1) c||c+1, tolerance=number of matching digits, \
#                        arithmetics to use =[R,RHP,RVHP,RGMP] ] , ... ]
shortsettings = [ shorttests,\
    [ 'cut', 1000,1000,'1e-88','1e-35','RGMP'],\
    [ 'full', 10,1000,'1e-40','1e-15','RGMP'] ]
    #[ 'cut', 5,1000,'1e-15','1e-7','RHP'] ]

longsettings = [ longtests, \
    [ 'cut', 100,1000,'1e-88','1e-30','RGMP'],\
    [ 'full',3,1000,'1e-40','1e-15','RGMP'] ]
    #[ 'cut', 3,1000,'1e-15','1e-7','RHP'] ]

# end input section
#############################

arraysize=500 # max size of array job
BHDIR='/mt/home/daniel/workspace/BHlib/'
BHRUNDIR=os.getcwd()

log = open('tests.log','w')
#
def ReadProcesses( setting , label ):
    filenames=setting.pop(0)
    nampl=0
    ntot=0
    for file in filenames:
        ffile=open(BHDIR+'test/CollTest/'+file)
        data=ffile.readlines()
        ffile.close()
        list1='list1=( '
        list2='list2=( '
        for line in data:
            if '#' not in line:
                line=line.replace('\t',' ')
                line=line.replace('  ',' ')
                line=line.replace('  ',' ')
                line=line.replace('  ',' ')
                line=line.split(' ')
                pro=line.pop(0)
                pro=pro.replace('\t','')
                for cs in line:
                    nampl+=1
                    csnew=cs.replace('\t','') 
                    csnew=csnew.replace('\n','') 
                    csnew=csnew.replace(' ','') 
                    list1+=pro+' '
                    list2+=csnew+' '
        list1+=')\n'
        list2+=')\n'
    submitfiles=[]
    for set in setting:
        [ part , PS , prec , mass , tol , arithm ] = set
        array=nampl*PS
        ntot+=array
        if array > arraysize :
            arraysizeloc=arraysize
            loopinterval=array/arraysize
            if (array/arraysize)*arraysize!=array:
                loopinterval+=1
        else:
            arraysizeloc=array
            loopinterval=0
        #splitup=range(0,array,array/arraysize)
        #splitup.append(array)
        #list0='list0=( '
        #for ival in splitup:
        #    list0+=str(ival)+' '
        #list0+=')\n'
        name='CollTest_'+label+'_'+part+'.sh'
        submitfiles.append(name)
        test=open(name,'w')
        test.write('#!/bin/bash\n')
        test.write('#PBS -N CollTest_%s_%s\n' % (label,part))
        test.write('#PBS -q route\n')
        test.write('#PBS -e log/CollTest/\n')
        test.write('#PBS -o log/CollTest/\n')
        test.write('##PBS -m ea  # switch for output\n')
        test.write('#PBS -k n\n')
        test.write('#PBS -m n\n')
        test.write('#PBS -t %s-%s\n' % (0,(arraysizeloc-1)) ) 
        test.write('#PBS -V\n')
        test.write('#PBS -l walltime=23:59:00,nodes=1:ppn=1\n')
        test.write('\n')
        test.write('cd %s\n' % BHRUNDIR )
        test.write('\n')
        test.write('# \n')
        test.write('# total number of amplitudes: %s \n' % nampl)
        test.write('# total number of PS points: %s \n' % PS)
        test.write('# \n')
        test.write('\n')
        test.write('\n')
        test.write('export BH_GMP_PRECISION=%s \n' % prec)
        test.write('MASS=\"%s\"\n' % mass)
        test.write('TOL=\"%s\"\n' % tol)
        test.write('ARITHM=%s\n' % arithm)
        test.write('PART=%s\n' % part)
        test.write('\n')
        test.write('%s\n' % list1)
        test.write('\n')
        test.write('%s\n' % list2)
        test.write('\n')
	#        test.write('%s\n' % list0)
        if loopinterval > 0:
            test.write('let \" START = $PBS_ARRAYID *  %s \"\n' % (loopinterval) )
            test.write('let \" END = $START + %s \"\n' % (loopinterval-1) )
            test.write('\n')
            test.write('MAX=%s \n' % (array-1))
            test.write('if [ $END -gt $MAX ] ; then\n')
            test.write('    END=$MAX  \n')
            test.write('fi\n')
            #
            test.write('\n')
            test.write('\n')
            test.write('for IVAL in `seq $START $END`;\n')
            test.write('do \n')
        else:
            test.write(' let \" IVAL = $PBS_ARRAYID \"\n')
            # 
        test.write(' let \" PRO = $IVAL / %s \"\n' % PS )
        test.write(' let \" SEED = 1000 + ( $IVAL %%  %s )\" \n' % PS ) 
        test.write('\n')
        test.write(' ./CollinearLimit -pro ${list1[$PRO]} -cs ${list2[$PRO]} -mass $MASS -seed $SEED -arithm $ARITHM -part $PART -tol $TOL \n')
        test.write(' if [ "${?}" -ne "0" ] ; then\n')
        test.write('     echo "${list1[$PRO]}	 ${list2[$PRO]}" >> CollTest.log \n')
        #test.write('      echo  "./CollinearLimit -pro ${list1[$PRO]} -cs ${list2[$PRO]} -mass $MASS -seed $SEED -arithm $ARITHM -part $PART -tol $TOL \" >> CollTest.log  \n')
        test.write(' else \n')
        test.write('      echo "pass">> CollTest.log \n')
        test.write(' fi \n')
        test.write('\n')
        if loopinterval > 0:
            test.write('done \n')
            test.write('\n')
        test.close()
    nprevtests=0
    if os.path.isfile('CollTest.log'):
        log=open('CollTest.log')
        nprevtests=int(log.readline())
        log.close()
    log=open('CollTest.log','w') 
    log.write('%s\n' % (ntot+nprevtests) )
    log.close()    
    for file in submitfiles:
        cmd='qsub '+file   
        print '  ',cmd
        os.system(cmd)

def Summarize():
   log=open('CollTest.log')
   data=log.readlines()
   log.close()
   ntests=int(data.pop(0))
   ndone=len(data)
   reddata=[]
   for entry in data:
     if not 'pass' in entry:
        reddata.append(entry)
   err=len(reddata)
   log=open('CollTest.stat','w')
   log.write('#############################\n')
   log.write('#  nbr ampls %s \n' % ntests)
   log.write('#  nbr err %s \n' % err)
   log.write('#  nbr pass %s \n' % (ndone - err) )
   log.write('#\n')
   log.write('#  nbr left %s \n' % (ntests-ndone) )
   log.write('#############################\n')
   log.write('\n')
   log.write('\n')
   #
   for entry in reddata:
      log.write('%s' % entry ) 
   log.write('\n')
   log.write('#############################\n')
   log.write('#  nbr ampls %s \n' % ntests)
   log.write('#  nbr err %s \n' % err)
   log.write('#  nbr pass %s \n' % (ndone - err) )
   log.write('#\n')
   log.write('#  nbr left %s \n' % (ntests-ndone) )
   log.write('#############################\n')
   log.close()
   os.system('cat CollTest.stat')


if os.path.isfile('CollTest.log'):
    Summarize()
else:
    os.system('make CollinearLimit')
    os.system('touch BH_debug.dat')
    print 'submitting:'
    ReadProcesses(shortsettings, 'short')
    ReadProcesses(longsettings, 'long')




