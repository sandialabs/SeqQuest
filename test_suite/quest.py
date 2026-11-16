#!/usr/bin/env python
"""
Command-line interface to Quest scripts.

cleanup           Remove Quest temporary files
clean             Equivalent to cleanup
setup    fname    Set up a directory for a Quest run on *fname*
[setroot dname    Set Quest root directory to dname]
run      (fname)  Run the quest job. If *fname* is specified,
                  use that as the input. Otherwise, use lcao.in
all      fname    Setup/run/cleanup Quest job on fname
mdall    fname    Same as all, but don't purge vxyz file before
                  run (useful for MD restart)
mdanal   fname    Analyzes an output file from an MD run
"""

from glob import glob
from os import system,unlink,rename
from os.path import abspath,exists,join
from re import compile
from shutil import copy
from math import sqrt

# Setup directories. The following should work for the standard
#  seqquest directory structure, although not 100% certain it
#  will work with Windows because of Unix specific stuff like
#  ".." and "/"
#
# Until 'quest setroot' function is working, this is the only
#  option

quest_home = abspath("../..")
quest_exe_dir = join(quest_home,'seqquest/obj')
quest_atom_dir = join(quest_home,'atoms/library')
quest_lda_dir = join(quest_atom_dir,'atoms_lda')
quest_pbe_dir = join(quest_atom_dir,'atoms_pbe')
quest_exe = join(quest_exe_dir,'lcao.x')

def isnear(a,b,tol=0.001): return abs(a-b)<tol

def grep(fname,pat,**kwargs):
    """
    Return the last matching line of *pat* in the file contained in fname.
    If called with return_all=True, will return all matching lines.
    """
    if type(pat) == type(''): pat = compile(pat)
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

def get_functional(fname):
    """
    Read the quest input file in *fname* and return the DFT functional,
    defaulting to lda if none is found.
    """
    file = open(fname)
    functional = 'lda'
    while 1:
        line = file.readline()
        if not line: break
        line = line.strip()

        if line.startswith('functional'):
            functional = file.readline().strip().lower()
    return functional

def get_atom_files(fname):
    """
    Return a list of all of the atom file in quest input file *fname*.
    """
    file = open(fname)
    atom_files = []
    while 1:
        line = file.readline()
        if not line: break
        line = line.strip().lower()
        if line.startswith('atom file'):
            atom_files.append(file.readline().strip())
    return atom_files

def quest_setup(fname,**kwargs):
    """
    Set up files in preparation to run a quest calculation.
    """
    verbose = kwargs.get('verbose',False)
    raise_errors = kwargs.get('raise_errors',False)
    remove_vxyz = kwargs.get('remove_vxyz',True)
    parallel_run = kwargs.get('parallel_run',False)

    
    # Copy input file and executable
    if fname != 'lcao.in' :
        if verbose:
            print "Copyring input file from %s to lcao.in" % fname
        copy(fname,'lcao.in')
    if verbose:
        print "Copying quest executable from %s to lcao.x" % quest_exe

    if parallel_run:
        copy(quest_tp_exe,'lcao.x')
    else:
        copy(quest_exe,'lcao.x')

    # Get atom files
    functional = get_functional(fname)
    atom_files = get_atom_files(fname)
    if functional == 'pbe':
        atom_search_dir = quest_pbe_dir
    else:
        atom_search_dir = quest_lda_dir
    for atom_file in atom_files:
        target = join(atom_search_dir,atom_file)
        if verbose:
            print "Searching for %s in %s" % (atom_file,atom_search_dir)
        if exists(target):
            if verbose:
                print "     found!"
            copy(target,atom_file)
        else:
            if verbose:
                print "     not found; please copy by hand"
            if raise_errors:
                raise "Atom file %s not found" % atom_file
    if remove_vxyz and exists('lcao.vxyz'):
        unlink('lcao.vxyz')
    return

def quest_setup_qsub(fname,**kwargs):
    kwargs['parallel_run'] = True
    quest_setup(fname,**kwargs)
    jobname = fname.replace('.in','')
    procs = kwargs.get('procs',1)
    nodes,rem = divmod(procs,2)
    if rem: nodes += 1
    exefile = 'lcao.x'
    qname = fname.replace('.in','.qsub')
    open(qname,'w').write(qsub_template % locals())
    print "The file %s has been created and can now be submitted" % qname
    print "to the batch system."
    return
    
def quest_run(fname,**kwargs):
    """
    Run a quest job.
    """
    do_preclean = kwargs.get('do_preclean',False)
    do_clean = kwargs.get('do_clean',False)
    do_setup = kwargs.get('do_setup',False)

    oname = quest_oname(fname)

    if do_preclean:
        quest_cleanup(fname,**kwargs)
    if do_setup:
        kwargs['raise_errors']=True
        quest_setup(fname,**kwargs)

    system("./lcao.x >& %s " % oname )

    energy = quest_energy(fname)

    if do_clean: quest_cleanup(fname,**kwargs)
    return energy

def quest_filenamer(fname,ending):
    """Given the name of a Quest input file, rename the file
    intelligently"""
    if not ending.startswith("."):
        ending = "." + ending
    if fname.endswith(".in"):
        return fname.replace(".in",ending)
    return fname + ending

def quest_oname(fname): return quest_filenamer(fname,".out")
def quest_hname(fname): return quest_filenamer(fname,".hist")
def quest_pname(fname): return quest_filenamer(fname,".post")
def quest_pname(fname): return quest_filenamer(fname,".vec")

def quest_cleanup(fname,**kwargs):
    """
    Delete temp files after a Quest run.
    """
    if exists("lcao.hist"):
        rename("lcao.hist",quest_hname(fname))
    if exists("lcao.post"):
        rename("lcao.post",quest_pname(fname))
    files = glob('lcao*') + glob('core.*')
    keep_list = ['lcao.x','lcao.in','lcao.vxyz']

    for fname in files:
        if fname not in keep_list: unlink(fname)
    return

def quest_energy(fname):
    """
    Obtain the energy from the quest output file in *fname*.
    """
    oname = quest_oname(fname)
    line = grep(oname,'TOTAL ENERGY')
    result = parseline(line,'xxxxf')
    if not result: result = 0
    return result

def quest_get_md_data(fname):
    file = open(fname)
    mddata = []
    while 1:
        line = file.readline()
        if not line: break

        if line.startswith(' QMDSUBVCOM'):
            vx,vy,vz = parseline(file.readline(),'fff')
        elif line.startswith(' WARNING: likely uninit'):
            # Reprint the warning messages now, since they
            #  likely signify an error
            print line,
        elif line.startswith(' <QMD'):
            nat = parseline(file.readline(),'xd')
            Tinit, Tcurr = parseline(file.readline(),'xxxxff')
            Ekin,Epot,Etot = parseline(file.readline(),'xxfff')
            P = parseline(file.readline(),'xxf')
            istep,dt,ntsteps,nteqsteps = parseline(file.readline(),'xxxxxdfdd')
            vcomx,vcomy,vcomz = parseline(file.readline(),'xxxxfff')
            skip = file.readline()
            skip = file.readline()
            mddatum = []
            for i in range(nat):
                iat,mass,x,y,z,vx,vy,vz = parseline(file.readline(),'dfffffff')
                mddatum.append((mass,x,y,z,vx,vy,vz))
            mddata.append((nat,Tinit,Tcurr,Ekin,Epot,Etot,P,istep,dt,ntsteps,
                           nteqsteps,vcomx,vcomy,vcomz,mddatum))
    return mddata

def quest_extract_trajectory_record(mddata,recordkey):
    index = dict(nat=0,Tinit=1,Tcurr=2,Ekin=3,Epot=4,Etot=5,P=6,istep=7,
                 dt=8,nsteps=9,nteqsteps=10,vcomx=11,vcomy=12,vcomz=13,
                 atomdata=14)
    if recordkey not in index: raise "Invalid key: %s" % recordkey
    i = index[recordkey]
    return [datum[i] for datum in mddata]

def quest_md_total_momentum(vxyz):
    px = py = pz = 0
    for mass,x,y,z,vx,vy,vz in vxyz:
        px += mass*vx
        py += mass*vy
        pz += mass*vz
    return px,py,pz

def quest_time_anal(fname):
    """\
    Summarize the timer data from a Quest output file.
    """
    import re
    # I believe this regular expression will work regardless of whether the wall clock
    #  or the standard timers are used.
    wpat = re.compile("([\w\s\/]+)\+\s*([\d\.]+)\s*=\s*[\d\.]*\s*seconds \w+\s*\(\s*([\d\.]+)\s*usr\s*([\d\.]+)\s*sys\)")
    timers = {}
    for line in open(fname):
        if wpat.match(line):
            match = wpat.match(line)
            words = match.groups()
            title = words[0].strip()
            total = float(words[1].strip())
            user = float(words[2].strip())
            sys = float(words[3].strip())
            #print title,user,sys
            tot_tot,user_tot,sys_tot = timers.get(title,(0,0,0))
            tot_tot += total
            user_tot += user
            sys_tot += sys
            timers[title] = tot_tot,user_tot,sys_tot
    # Sort the values by the total (u+s) time
    templist = [(t,title,u,s) for title,(t,u,s) in timers.items()]
    templist.sort()
    templist.reverse()
    print "%20s %8s %8s %8s %8s" % ("Heading","Total","User","System","Loss")
    print "%20s %8s %8s %8s %8s" % ("-------","-----","----","------","------")
    tott = totu = tots = 0
    for t,title,u,s in templist:
        print "%20s %8.3f %8.3f %8.3f %8.3f" % (title,t,u,s,t-u-s)
        tott += t
        totu += u
        tots += s
    print "%20s %8s %8s %8s %8s" % ("-------","-----","----","------","------")
    print "%20s %8.3f %8.3f %8.3f %8.3f" % ("Totals",tott,totu,tots,tott-totu-tots)
    return
            
def quest_md_isok(fname):
    oname = quest_oname(fname)
    quest_get_md_data(oname) # Just to make sure everything can be parsed
    # Had to clobber these, since Peter changed the final MD output.
    # This routine now just makes sure that the file exists and can
    # be parsed.
    #md_exit_lines = grep(oname,'MD Run Completed',return_all=True)
    #isok = len(md_exit_lines) == 1
    isok = True
    return isok

def command_line():
    """
    Command-line interface to Quest scripts.

    cleanup           Remove Quest temporary files
    clean             Equivalent to cleanup
    setup    fname    Set up a directory for a Quest run on *fname*
    [setroot dname    Set Quest root directory to dname]
    run      (fname)  Run the quest job. If *fname* is specified,
                      use that as the input. Otherwise, use lcao.in
    all      fname    Setup/run/cleanup Quest job on fname
    mdall    fname    Same as all, but don't purge vxyz file before
                      run (useful for MD restart)
    mdanal   fname    Analyzes an output file from an MD run
    """
    import sys

    commands = ['cleanup','clean','setup','run','all','setroot','mdall',
                'qsub','mdanal','tanal']


    if len(sys.argv) <= 1:
        print __doc__
        sys.exit()

    command = sys.argv[1]

    if len(sys.argv) > 2:
        argument = sys.argv[2].lower()
    else:
        argument = None

    assert command in commands

    if command.startswith('clean'):
        quest_cleanup("test.in")
    elif command == 'setup':
        quest_setup(argument)
    elif command == 'run':
        energy = quest_run(argument)
        print "energy = ",energy
    elif command == 'all':
        energy = quest_run(argument,do_preclean=True,
                           do_setup=True,do_clean=True)
        print "energy = ",energy
    elif command == 'mdall':
        # same as all, but doesn't remove the vxyz file, so an MD
        # run can be restarted
        energy = quest_run(argument,do_preclean=True,
                           do_setup=True,do_clean=True,remove_vxyz=False)
        print "energy = ",energy
    elif command == 'qsub':
        quest_setup_qsub(argument)
    elif command == 'mdanal':
        data = quest_get_md_data(argument)
        print "%6s %10s %12s %12s %12s" % ("Step","E total","Px","Py","Pz")
        for datum in data:
            istep = datum[7]
            Etot = datum[5]
            vxyz = datum[14]
            px,py,pz = quest_md_total_momentum(vxyz)
            print "%-6d %10.4f %12.4e %12.4e %12.4e" % (istep,Etot,px,py,pz)
    elif command == 'tanal':
        quest_time_anal(argument)
    return

if __name__ == '__main__': command_line()
