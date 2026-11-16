#!/usr/bin/env python
"""\
Parallel version of the Quest test suite:

parallel_suite.py [-r]      Run the test suite
parallel_suite.py -s        Just setup the suite, don't run anything
parallel_suite.py -c        Check the results
"""


import getopt
import glob
import os
import re
import shutil
import sys

# qsub header
qsub_header = """\
#!/bin/bash

#PBS -l nodes=2:ppn=2
#PBS -l walltime=96:00:00
#PBS -q nw
#PBS -N quest_test_suite
#PBS -A FY081517
#PBS -j oe
#PBS -o quest_test_suite.log
"""

# Define a few directories we'll use later
QUEST = '/home/rmuller/Programs/Quest/seqquest'
EXE = os.path.join(QUEST,'obj_tp')
TEST = os.path.join(QUEST,'test_suite')

results = {
    'h2.in' : -1.786568,
    'h2_geopt.in' : -1.782201,
    'h2_gga.in' : -1.8399273,
    'h2_geopt_gga.in' : -1.8355158,
    'si_111.in' : -62.258670,
    'si_111_cell.in' : -62.308584,
    'si_111_complex.in' : -62.269961,
    'gan.in' : -58.798245,
    # Long jobs
    'si_112.in' : -125.257303,
    'c6h6.in' : -75.420703,
    'ru_slab.in' : -232.249206,
    'as_test.in' : -28.3216914
    }

def grep(fname,pat,**kwargs):
    """
    Return the last matching line of *pat* in the file contained in fname.
    If called with return_all=True, will return all matching lines.
    """
    if type(pat) == type(''): pat = re.compile(pat)
    result = []
    for line in open(fname):
        if pat.search(line): result.append(line)
    if kwargs.get('return_all',False): return result
    if len(result) > 0: return result[-1]
    return ""

def parseline(line,format):
    """
    Return the line split and with particular formats applied.

    Format characters:
    x     Skip this record
    s     Convert the record to a string
    f     Convert the record to a float
    d     Convert the record to an int
    i     Convert the record to an int

    Examples:
    >>> parseline('H  0.0 0.0 0.0','sfff')
    'H', 0.0, 0.0, 0.0
    >>> parseline('H  0.0 0.0 0.0','xfff')
    0.0, 0.0, 0.0
    """
    xlat = {'x':None,'s':str,'f':float,'d':int,'i':int}
    result = []
    words = line.split()
    if len(words) < len(format): return None
    for i in range(len(format)):
        f = format[i]
        trans = xlat.get(f)
        if trans: result.append(trans(words[i]))
    if len(result) == 0: return None
    if len(result) == 1: return result[0]
    return result

def run_suite(jobs,**kwargs):
    dosubmit = kwargs.get('dosubmit',False)
    
    qsub_lines = [qsub_header]

    preclean()
    # Copy everything we need
    shutil.copy(os.path.join(EXE,'lcao_tp.x'),'quest_tp.e')
    for job in jobs:
        shutil.copy(os.path.join(TEST,job),job)
        qsub_lines.append('# Running %s' % job)
        qsub_lines.append('cp %s lcao.in' % job)
        qsub_lines.append('mpirun -v -np 4 quest_tp.e')
        qsub_lines.append('mv lcao.out %s' % job.replace('.in','.out'))
        qsub_lines.append('rm lcao*')
        qsub_lines.append('')

    # Write/run the script file:
    open('parallel_suite.qsub','w').write('\n'.join(qsub_lines))
    if dosubmit: 
        os.system('qsub parallel_suite.qsub')
    return

def quest_energy(oname):
    """
    Obtain the energy from the quest output file in *fname*.
    """
    line = grep(oname,'TOTAL ENERGY')
    result = parseline(line,'xxxxf')
    if not result: result = 0
    return result

def check_results(jobs,**kwargs):
    tolerance = kwargs.get('tolerance',1e-4)
    njobs = len(jobs)
    npassed = 0
    nfailed = 0
    failed = []
    for job in jobs:
        oname = job.replace('.in','.out')
        print "---Checking the results for %s---" % job
        if os.path.exists(oname):
            print "%s found!" % oname
            print "Energy should be", results[job]
            energy = quest_energy(oname)
            print "Energy is ", energy
            error = abs(energy-results[job])
            if error < tolerance:
                print "Job passed"
                npassed += 1
            else:
                print "****JOB FAILED****"
                nfailed += 1
                failed.append(job)
            print ""
        else:
            print "****JOB %s NOT FOUND****" % job
    if nfailed == 0:
        print "Full test suite passed"
    else:
        print "%d jobs failed, out of a total %d" % (nfailed,njobs)
        print failed
    return

def preclean():
    for fname in glob.glob('lcao*') + glob.glob('*out'):
        os.unlink(fname)
    return
        
def main():
    opts,args = getopt.getopt(sys.argv[1:],'acrs')
    default_jobs = ['h2.in','h2_geopt.in','h2_gga.in','h2_geopt_gga.in',
                    'si_111.in','si_111_cell.in','si_111_complex.in','gan.in']
    long_jobs = ['si_112.in','c6h6.in','ru_slab.in','as_test.in']

    run_type = 'run'
    joblist = default_jobs
    for key,val in opts:
        if key == '-c':
            run_type = 'check'
        if key == '-s':
            run_type = 'setup'
        if key == '-a':
            joblist = default_jobs + long_jobs

    if args:
        jobs = args
    else:
        jobs = joblist
        
    if run_type == 'setup':
        dosubmit = False
    else:
        dosubmit = True

    #print run_type,jobs,dosubmit
    if run_type == 'check':
        check_results(jobs)
    else:
        run_suite(jobs,dosubmit=dosubmit)

if __name__ == '__main__': main()
