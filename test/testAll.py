import time
import os
import re
import stat
import shlex, subprocess
import sys


allTests=[
    'W_ME2_check','W1j_ME2_check','W2j_ME2_check','W3j_partial_and_ME2_check','W4j_ME2_check','W5j_ME2_check','W4j_ME2_check_auto',
    'gkm_check','sub_leading_test',
    'Z4q_test','Z2q2G_test','Z2qNg_test','2q2GNg_test',
    'photon_test',
    'ME_5parton_2lept_check','ME_6parton_2lepton_lc_check',
    'Y1j_ME2_check','Y2j_ME2_check','Y3j_ME2_check','Y4j_ME2_check','Y5j_ME2_check',
    'Z1j_ME2_check','Z2j_ME2_check','Z3j_ME2_check','Z4j_ME2_check','Z5j_ME2_check',
    '2j_ME2_check','G3_3j_ME2_check','3j_ME2_check','4j_ME2_check','4j_loop_ME2_check','5j_ME2_check','Vnj_bornME2_check','higgs_check',
    'YY_ME2_check','YY1j_ME2_check','YY2j_ME2_check',
    'Wm_PS','Wm1_PS','Wm2_PS',
    'Z_PS','Z1_PS','Z2_PS',
    '2j_PS','3j_PS',
    'Y1_PS','Y2_PS',
    'Z3_PS','Wm3_PS','Y3_PS']

allTestsPublic=['W_ME2_check',
          'W1j_ME2_check',
          'W2j_ME2_check',
          'W3j_partial_and_ME2_check',
          'Y1j_ME2_check',
          'Y2j_ME2_check',
          'Y3j_ME2_check',
          'Z1j_ME2_check',
          'Z2j_ME2_check',
          'Z3j_ME2_check',
          '2j_ME2_check',
          '3j_ME2_check',
          '4j_ME2_check',
          '4j_loop_ME2_check',
          'Wm_PS','Wm1_PS','Wm2_PS',
          'Wm3_PS',
          'Z_PS',
          'Z1_PS',
          'Z2_PS',
          'Z3_PS',
          'Y1_PS',
          'Y2_PS',
          'Y3_PS',
          '2j_PS',
          '3j_PS']

allTests=allTestsPublic         
          
log = open('tests.log','w')

def getStatus(done,toSubmit,submitted,running,failed):
    mode_ref=os.stat('Test_Start')
    toTest=list(toSubmit)
    for test in toTest:
        try:
            logFile='status/%s.submitted' % test
            mode=os.stat(logFile)
            if mode[stat.ST_MTIME] > mode_ref[stat.ST_MTIME] :
                log.write('%s submitted\n' % test)
                submitted.append(test)
                toSubmit.remove(test)
        except:
                pass
                #print '%s not there' % test
    toTest=list(submitted)
    for test in toTest :
        try:
            logFile='status/%s.started' % test
            mode=os.stat(logFile)
            if mode[stat.ST_MTIME] > mode_ref[stat.ST_MTIME] :
                running.append(test)
                submitted.remove(test)
        except:
                pass
                #print '%s not there' % test
    toTest=list(running)
    for test in toTest:
        try:
            logFile='status/%s.exit' % test
            mode=os.stat(logFile)
            if mode[stat.ST_MTIME] > mode_ref[stat.ST_MTIME] :
                exitStatusFile=open(logFile)
                line=exitStatusFile.readlines()[0]
                exitStatus=int(line)
                done[test]=exitStatus
                running.remove(test)
        except:
                pass
                #print '%s not there' % test
    for f in [j for j in done.keys() if done[j] != 0 ]:
        if f not in failed:
        	failed.append(f)
        


singlescript='/mt/home/daniel/workspace/BHlib/test/singleTest.sh'

def submitBsub(toSubmit,submitted):
    submittedNow=list()
    copyOfToSubmit=list(toSubmit)
    for test in copyOfToSubmit:
        ARGS="%s $PWD" % (test)
        cmd='bsub -q 1nd -J %s %s %s ' % (test,singlescript,ARGS)
        toSubmit.remove(test)
        os.system('%s >& /dev/null' % cmd)
        submitted.append(test)
        submittedNow.append(test)
    return submittedNow


walltimes={'qmedium':"-l walltime=23:59:59",
           'qbatch':"1:59:59",
           'qlong':"-l walltime=47:59:59"
}

localOptions={
'medium.q':'#PBS -l nodes=1:u12'
}
def submitQsub(toSubmit,submitted,queue='qlong'):
    submittedNow=list()    
    copyOfToSubmit=list(toSubmit)
    for test in copyOfToSubmit:
       	pwd=os.getcwd()
        ARGS="%s %s" % (test,pwd)
        jobName=test
        script=open('script_%s.sh' % jobName ,'w')
        script.write("#PBS -N %s \n" % jobName)
        script.write("#PBS -q %s\n"  % queue)
        local=localOptions.get(queue,None)
        if local:
            script.write(local+'\n')
        script.write("%s %s " % (singlescript,ARGS) )

        script.close()

        cmd='qsub %s script_%s.sh' % (walltimes.get(queue,''),jobName)

        toSubmit.remove(test)
        os.system('%s >> tt.log 2>&1 ' % cmd)
        submitted.append(test)
        submittedNow.append(test)
    return submittedNow




def submitBash(toSubmit,submitted,nbrThreads):
    submittedNow=list()
    while len(submitted) < nbrThreads:
        import os
        
        if not toSubmit:
            return
        test=toSubmit.pop()
        os.system("%s %s $PWD &" % (singlescript,test) )
        submitted.append(test)
        submittedNow.append(test)
    return submittedNow

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def submit(toSubmit,submitted):
    #    return submitBash(toSubmit,submitted,1)
    if which('bsub'):
        #print 'Using bsub'
        submittedNow=submitBsub(toSubmit,submitted)
    elif which('qsub'):
        #print 'Using qsub'
        submittedNow=submitQsub(toSubmit,submitted)
    #print 'Using normal jobs'
    else:
        submittedNow=submitBash(toSubmit,submitted,2)
    for s in submittedNow:
        os.system("touch status/%s.submitted" % (s) )

useNcurses=True

if __name__ == '__main__':
    import sys
    resubmitFailed=False
    if 'resubmitFailed' in sys.argv:
        resubmitFailed=True
        print 'Resubmitting failed jobs...'
        time.sleep(1)

    try:
        if useNcurses:
            import curses
            stdscr = curses.initscr()
            curses.resizeterm(40,100)
            curses.noecho()
            curses.cbreak()
            stdscr.keypad(1)
            curses.start_color()
            curses.init_pair(1, curses.COLOR_RED, curses.COLOR_BLACK)
            curses.init_pair(2, curses.COLOR_GREEN, curses.COLOR_BLACK)
            curses.init_pair(3, curses.COLOR_BLUE, curses.COLOR_BLACK)
            curses.init_pair(4, curses.COLOR_YELLOW, curses.COLOR_BLACK)
        done={}
        toSubmit=list(allTests)
        running=[]
        submitted=[]
	failed=[]
        offset=max([len(x) for x in allTests])+5
        if useNcurses:
            stdscr.addstr(0,5,'--- BH tests ---')
        else :
            print '--- BH tests ---'
        firstTime=True
        while toSubmit or submitted or running or firstTime:
            firstTime=False
            getStatus(done,toSubmit,submitted,running,failed)
            if resubmitFailed:
                toSubmit=toSubmit+failed
            submit(toSubmit,submitted)
            #print done
            for row,t in enumerate(allTests):
                if useNcurses:
                    voffset=0
                    if row%2 == 1:
                        voffset=45
                    stdscr.addstr(row/2+1,2+voffset,'%s :' % t)
                else:
                    print ("  %s :" % t),
                if t in done:
                    es=done[t]
                    if es == 0:
                        if useNcurses:
                            stdscr.addstr(row/2+1,offset+voffset,'PASSED    ',curses.color_pair(2) )
                        else:
                            print "  PASSED"
                    else:
                        stdscr.addstr(row/2+1,offset+voffset,'FAILED (%d)  ' % es,curses.color_pair(1) )
                if t in submitted:
                        if useNcurses:
                            stdscr.addstr(row/2+1,offset+voffset,'SUBMITTED  ',curses.color_pair(3) )
                        else:
                            print "  SUBMITTED"

                if t in running:
                        if useNcurses:
                            stdscr.addstr(row/2+1,offset+voffset,'RUNNING  ',curses.color_pair(4) )
                        else:
                            print "  RUNNING"

                if t in toSubmit:
                        if useNcurses:
                            stdscr.addstr(row/2+1,offset+voffset,'WAITING' )
                        else:
                            print "  WAITING"

            if useNcurses:
                stdscr.refresh()
            if len(done) != len(allTests):
                time.sleep(1)
        n=len(allTests)
        if useNcurses:
            stdscr.addstr(n/2+2,1,'Press Crtl+C to quit')
            stdscr.refresh()
        else :
            print 'Press enter to quit'
        #curses.noecho()
        #curses.cbreak()
        #stdscr.timeout(3000000)
        #stdscr.refresh()
        #c=stdscr.getch()
        raw_input()
        if useNcurses:
            curses.nocbreak(); stdscr.keypad(0); curses.echo()        
            curses.endwin()
        #print "here"
        sys.exit(0)
    finally:
        print 'exception thrown....', sys.exc_info()[0]
        if useNcurses:
            curses.nocbreak(); stdscr.keypad(0); curses.echo()        
            curses.endwin()
            
