c***********************************************************************
c     QMD - Quest MD Module
c     Written: Rick Muller 2004
c
c     Public subroutines
c      Configuring the MD module ...
c      QMDSETIWR( IWR ) - Set the write unit for QMD output
c      QMDTYPSET( iphase, mdlabel )      Set MD method/thermostat
c         iphase=1: initialization/equilibration
c         iphase=2: temperature ramp
c         iphase=3: equilibration
c         iphase=4: molecular dynamics
c      QMDTEMPSET( ti )          Set MD ionic temperature
c      QMDTSTEPSET( dt )         Set time step interval
c      QMDNEQSTEPSET( neqsteps ) Set number of equilibration steps
c      QMDNTSTEPSET( ntsteps )   Set number of thermalized MD steps
c      QMDSTEPSGET( nsteptot )   Returns total number of steps in MD
c      Launching a configured MD run ...
c      QMDINIT           Initialize the Quest MD run
c      Updating/tracking a MD run in progress ...
c      QMDUPDATE         Update the MD variables at the end of a
c                          geometry optimization step
c      QMDPRINT          Print a snapshot of the MD configuration
c      QMDCONVCHK        Check for whether MD is done
c      QMDEPSET          Set the potential energy from the El Energy
c
c     Private subroutines
c      GET_VCOM          Compute the velocity of the center of mass
c      GET_EKIN          Convert kinetic energy to temperature
c      GET_TEMP_E        Compute the temperature from the kinetic energy
c      GET_TEMP_V        Compute the temperature from the velocities
c      QMDSUBVCOM        Remove center of mass velocity
c      QMDSUBFCOM        Subtract the center of mass forces
c      QMDUPDATE_DAMPED  Minimization MD runs using a constant velocity damping
c      QMDUPDATE_TSCALE  Velocity scaled to temperature update
c      QMDUPDATE_NVE_LF  NVE ensemble, Leapfrog Verlet update
c      QMDUPDATE_LANGVN  NVT ensemble, Langevin thermostat, Leapfrog
c      QMDUPDATE_HOOVER  NVT ensemble, Nose-Hoover thermostat, Leapfrog
c      QMDUPDATE_BEREND  NVT ensemble, Berendsen thermostat, Leapfrog
c      QMDVINIT          Initialize velocities to a temperature (boltzmann)
c      QMDVDAMP          Damp the velocities by a factor
c      QMDVRESET         Reset all of the atomic velocities to a new temp
c
c     Submodules:
c      TMON /qmdtmon/    Used to monitor temperatures over atom types
c
c     QMD Module data:
c      natmd         i   Parameter - number of atoms module permits
c     Private data stored in common block /questmd/
c      tempinit      f   Target temperature
c      tempprnt      f   Current temperature (in synch with position)
c      tempcurr      f   Current temperature (in synch with velocity)
c      timestep      f   Time step
c      Ekinprnt      f   Kinetic energy (in synch with position)
c      Ekincurr      f   Kinetic energy (in synch with velocity)
c      Epot          f   Potential energy
c      shv           f   Nose-Hoover s (sometimes eta) coupling to heat bath
c      vel           f*  Velocities: leapfrog, kept at half steps w.r.t. the coords
c                        coordinates
c      vprint        f*  Velocities: in synch with positions (interpolated)
c      atom_mass     f*  Atomic Mass of each element
c      ieq_type      i   method for equilibration step
c      md_type       i   Type of Quest MD Run (see QMDTYPSET/QMDUPDATE)
c                        1   NVE Dynamics
c                        2   NVT Dynamics, Berendsen Thermostat
c                        3   NVT Dynamics, Hoover Thermostat
c                        4   NVT Dynamics, Scaled Temperatures
c                        5   Damped dynamics
c                        6   NVP Dynamics, Anderson Thermostat (n/a)
c                        7   Hugoniotstat Dynamics (n/a)
c                        8   NVT Dynamics, Langevin Thermostat (n/a)
c      ntsteps       i   Number of time steps
c      neqtsteps     i   Number of scaled temperature time steps for equilibration
c      md_natoms     i   Number of atoms
c      itimestep     i   Current time step
c      md_nmodes     i   3*Nat-6 or 3*Nat-5; actually just 3*Nat for now
c      imd_wr        i   Output unit for MD, IWR from the main routine
c      irestarted    i   0 if the MD case did not restart; 1 otherwise
c      temp_eq       f   target temperature for equil phase
c      temp_md       f   target temperature for MD phase
c      param_eq      f   parameter for thermostat control, equil phase
c      param_md      f   parameter for thermostat control, MD phase
c      parstat       f   parameter for currently active thermostat
c
c     Private constant data (see QMDCONST)
c      fconst        f   Converts Force in Ry/bohr to ???
c      tconst        f   Converts Temperature from Kelvin to Ry
c***********************************************************************
c Revision history:
c   8Feb13-PAS/2.65: use clean, full atomic symbol tables
c  25Apr08-PAS/2.61b: add energy monitor into monitor function
c  10Feb08-PAS/2.61: fix thermostats/initialization, wire thermostat controls
c  27Jun07-PAS/2.60: add data to restart file for clean thermostat continue
c***********************************************************************
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDSETIWR
c
c
      subroutine QMDSETIWR( IWR )
c---------------------------------------------------------------
c Purpose: set thermostat/type of molecular dyanmics
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c >>>> EXECUTABLE CODE:
c
      imd_wr = IWR
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDTYPSET
c
c
      subroutine QMDTYPSET( imdphase, mdlabel )
c---------------------------------------------------------------
c Purpose: set thermostat/type of molecular dyanmics
c
c Written: Richard P. Muller, 8-May-2004, for v2.58 (qmd)
c
c Revision history:
c   8May04-PAS: extract from RPM's routines into this routine
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
C      COMMON /qmddata/ timestep,itimestep,irestarted,
C     $ md_natoms,md_nmodes,atoms_mass, vel,vprint,
C     $ neqtsteps,ieq_type, temp_eq, param_eq,
C     $ ntsteps  , md_type, temp_md, param_md,
C     $ tempinit,tempcurr,Ekincurr,Epot,shv,parstat,
C     $ imd_wr
c
c input declarations:
      CHARACTER*(*)  mdlabel
c local declarations:
      CHARACTER*6  mdmeth
c
c The following are MD method id tags; they must be distinct.
      DATA imd_nve   ,imd_berend,imd_hoover,imd_tscale,imd_damped
     $    /    1     ,    2     ,    3     ,    4     ,    5     /
      DATA imd_anders,imd_hugoni,imd_langvn
     $    /    6     ,    7     ,    8     /
c
      DIMENSION  paramsmd(8)
      DATA  paramsmd / 0.d0, 400.d0, 100.d0, 0.98d0, 0.95d0,
     $       0.0d0,    0.d0, 100.d0 /
c paramsmd defines default control parameter for each thermostat:
c  1/NVE - none
c  2/Berend - t/time(fs), s=(1/2t)*(T/Told-1)
c  3/Hoover - t/time(fs), ds=(1/t^2)*(T/Told-1)
c  4/Tscale - f/scaling factor (simple: Tnew=  f*Told + (1-f)*T )
c  5/damped - f/scaling factor (simple: Tnew = f*Told )
c  6/anders - n/a
c  7/hugoni - n/a
c  8/langvn - t/time(fs)
c
c >>>> EXECUTABLE CODE:
c
      mdmeth = mdlabel(1:6)
c
      imdtype = 0
c
      if(     mdmeth .eq. 'nve   ' .or. mdmeth .eq. 'NVE   ' )then
c       NVE/Microcanonical update, leapfrog-verlet
        imdtype = imd_nve
c
      elseif( mdmeth .eq. 'langev' .or. mdmeth .eq. 'LANGEV' .or.
     $        mdmeth .eq. 'Langev' )then
c       NVT/Langevin thermostat, leapfrog-verlet
        imdtype = imd_langvn
c
      elseif( mdmeth .eq. 'berend' .or. mdmeth .eq. 'BEREND' .or.
     $        mdmeth .eq. 'Berend' )then
c       NVT/Berendsen thermostat, leapfrog-verlet
        imdtype = imd_berend
c
      elseif( mdmeth .eq. 'hoover' .or. mdmeth .eq. 'HOOVER' .or.
     $        mdmeth .eq. 'Hoover' .or.
     $        mdmeth .eq. 'NVT' )then
c       NVT/Nose-Hoover thermostat, leapfrog-verlet ***NVT default***
        imdtype = imd_hoover
c
      elseif( mdmeth .eq. 'tscale' .or. mdmeth .eq. 'TSCALE' .or.
     $        mdmeth .eq. 'Tscale' )then
c       Temperature scaling, leapfrog-verlet
        imdtype = imd_tscale
c
      elseif( mdmeth .eq. 'damped' .or. mdmeth .eq. 'DAMPED' )then
c       Damped dynamics (bleeding out K.E.), leapfrog-verlet
        imdtype = imd_damped
c
      elseif( mdmeth .eq. 'anders' .or. mdmeth .eq. 'ANDERS' .or.
     $        mdmeth .eq. 'Anders' )then
c       Anderson thermostat
        imdtype = imd_anders
c
      elseif( mdmeth .eq. 'hugoni' .or. mdmeth .eq. 'HUGONI' .or.
     $        mdmeth .eq. 'Hugoni' )then
c       Hugoniostat
        imdtype = imd_hugoni
c
      else
        write(imd_wr,*) 'QMDTYPSET: unknown MD method tag=',mdmeth
        call STOPXERR( 'qmdtypset: invalid dynamics method' )
c
      endif
c
c  Enter MD type into Module data:
      if( imdphase .eq. 1 )then
c       This is the initial equilibration phase
        ieq_type = imdtype
        param_eq = paramsmd(imdtype)
        write(imd_wr,*) '# equilibration param set to=',param_eq
      elseif( imdphase .eq. 4 )then
c       This is the thermalized MD phase
        md_type = imdtype
        param_md = paramsmd(imdtype)
        write(imd_wr,*) '# MD thermostat param set to=',param_md
      else
c       Default: thermalized MD phase
        call STOPXERR( 'mdtypset code error: invalid MD phase' )
        md_type = imdtype
        param_md = paramsmd(imdtype)
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDSETGETS
c
c
c---------------------------------------------------------------
c Purpose: public routines that set/get configurable MD params
c Written: Richard P. Muller, 8-May-2004, for v2.58 (qmd)
c Revision history:
c  17May04-PAS: renamed w/ QMD prefixes, moved to "public" section
c---------------------------------------------------------------
c
      subroutine QMDTEMPSET( temp )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c Configure temperature for MD run:
      tempinit = temp
      temp_md  = temp
c For the moment, also set the equilibration target temp to this, too:
      call QMDEQTSET( temp )
c
      RETURN
      END
c
      subroutine QMDEQTSET( temp )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c Configure temperature for equilibration
      temp_eq = temp
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDPARAMSET
c
c
      subroutine QMDPARAMSET( parmd )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c Configure thermostat parameter for MD run:
      param_md = parmd
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDEQPARSET
c
c
      subroutine QMDEQPARSET( parmd )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c Configure thermostat parameter for equilibration
      param_eq = parmd
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDTSTEPSET
c
c
      subroutine QMDTSTEPSET( dt_md )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c Configure timestep for MD run:
      timestep = dt_md
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDNEQSTEPSET
c
c
      subroutine QMDNEQSTEPSET( nstepeq )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c Configure number of equilibration steps:
      neqtsteps = nstepeq
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDTSTEPGET
c
c
      subroutine QMDNTSTEPSET( nstepmd )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c Configure number of thermalized MD steps:
      ntsteps = nstepmd
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDSTEPSGET
c
c
      subroutine QMDSTEPSGET( nsteptot )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c Total number of steps is equilibration + thermalized MD steps
      nsteptot = ntsteps + neqtsteps
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDINIT
c
c
      subroutine QMDINIT( iwr,nat,itypa,atmass,ratm,igstart )
c
c     Initialize variables for Quest MD module. Also do some
c     simple sanity checks. Had to pass in ratm so we can optionally
c     update the atom positions from the vxyz file
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c     Argument arrays
      DIMENSION  itypa(natmd),atmass(*),ratm(3,natmd)
c Local declarations:
      PARAMETER  (ntypchk=10)
      DIMENSION  itypchk(ntypchk)
c
c >>>> EXECUTABLE CODE:
c
c Initialize temperature-monitoring apparatus
      call TMONINIT( nat, itypa )
      call TMONRESET( )
c
      imd_wr = iwr
      irestarted = 0
c Set dynamics flag in geometry engine:
      call GTYPSET( 'DYNAMICS' )
c
      if (tempinit .lt. zero) then
         call STOPXERR('QMDINIT -- Temperature < 0')
      else if (tempinit .gt. 2000d0) then
         write(imd_wr,*) 'WARNING: Quest MD -- Temperature high',
     $        tempinit
         write(imd_wr,*) '         I hope you know what you are doing'
      endif
      tempcurr = tempinit
      shv = zero
c
c
      if (timestep .gt. 10d0) then
         write(imd_wr,*) 'WARNING: Quest MD -- Large time step entered',
     $        timestep
      else if (timestep .lt. 0.01d0) then
         write(imd_wr,*) 'WARNING: Quest MD -- Small time step entered',
     $        timestep
      endif
c
      if( ntsteps .le. 0 )then
         write(imd_wr,*) 'QMDINIT: invalid # time steps=',ntsteps
         call STOPXERR('QMDINIT - invalid # time steps')
      endif
c
      if( neqtsteps .lt. 0 )then
         write(imd_wr,*) 'QMDINIT: invalid # equil steps=',neqtsteps
         call STOPXERR('QMDINIT - invalid # equilibration steps')
      endif
      itimestep = 0
c
      if( nat .gt. natmd )then
         write(imd_wr,*) 'QMDINIT: too many atoms for MD module'
         write(imd_wr,*) 'Input number ot atoms, max=',nat,natmd
         call STOPXERR('QMDINIT - too many QMD atoms')
      elseif( nat .le. 0 )then
         write(imd_wr,*) 'QMDINIT: bad number of atoms=',nat
         call STOPXERR('QMDINIT - invalid # atoms')
      else
c        Install number of atoms into the MD module:
         md_natoms = nat
c
c        Set the number of modes to 3*nat (rather than 3*nat-6),
c        For periodic system it includes the lattice modes:
c        This is not really correct, but it's just a scaling factor
c
         md_nmodes = 3*md_natoms
         write(imd_wr,*) 'MD atoms, modes=',md_natoms,md_nmodes
      endif
c
c     Install the atom masses into the MD module
c
cqmd: treatment of atomic masses needs to be handled better
      do  ityp=1,ntypchk
        itypchk(ityp) = 0
      enddo
      do  i=1,md_natoms
        itypatom = itypa(i)
        ximass = atmass(itypatom)
        if( itypchk(itypatom) .eq. 0 )then
c         Check this mass vs. an unitialized mass
          diff = ABS( ximass - 1.d0 )
          if( diff .lt. 5.d-7 )then
            write(imd_wr,*) '***** WARNING: likely uninitialized mass'
            write(imd_wr,*) '>>>>> atom type #, mass=',itypatom,ximass
          endif
          itypchk(itypatom) = 1
        endif
        atom_mass(i) = ximass
      enddo
c
c     Acquire restart data for a continued run, if lcao.vxyz file exists
      call QMDRESTART( ratm )
c
c     If no restart, then Initialize velocities and kinetic energy
      if( irestarted .eq. 0 ) call QMDVINIT
c
c     Set geometry step to correct place.
c       Remember itimestep is zero-indexed, igstart is one-indexed
c       This might be done otherwise/where, but convention is ambiguous
      igstart = itimestep + 1
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDRESTART
c
c
      subroutine QMDRESTART( ratm )
c
c     Check for the existance of a restart file
c     If found, use it for the x,y,z,vx,vy,vz information for all atoms,
c               and use it to smoothly continue the thermostat.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c Output:
      dimension ratm(3,natmd)
c
c Local declarations:
      CHARACTER*3  sym
      CHARACTER*8  wantlbl,readlbl
c
c >>>> EXECUTABLE CODE:
c
      call QMDCONST( fconst, tconst )
      sym = ""
      vx = zero
      vy = zero
      vz = zero
      x = zero
      y = zero
      z = zero
c
      call FLEXIST(iexist, 'lcao.vxyz', 'FOR')
      if (iexist .eq. 1) then
         irestarted = 1
         call FLOPEN(ifl, 'lcao.vxyz', 'OLD', 'FOR',ierror)
         if( ierror.ne.0 )
     $        call STOPXERR( 'QMD-VINIT: error opening restart' )
c
c
         read(ifl,*) nat_restart
         if (nat_restart .ne. md_natoms)
     $        call STOPXERR( 'QMD-VINIT: wrong nat in restart' )
c
c     Skip title line
         read(ifl,*)
c
c      Input the atomic positions and velocities
         Ekin_dt = zero
         do i =1, nat_restart
            read(ifl,'(A2,6f12.8)',err=113,end=113) sym,x,y,z,vx,vy,vz
            vel(1,i) = vx
            vel(2,i) = vy
            vel(3,i) = vz
            ratm(1,i) = x
            ratm(2,i) = y
            ratm(3,i) = z
            xm = atom_mass(i)
            Ekin_dt = Ekin_dt + 0.5d0*xm*(vx*vx+vy*vy+vz*vz)/fconst
         enddo
c
c     Remove net momentum, if any:
         call QMDSUBVCOM
         call GET_EKIN( md_natoms,atom_mass,vel, Ekin )
         call GET_TEMP_E( Ekin,md_nmodes, temp )
         Ekincurr = Ekin
         tempcurr = temp
c
         write(imd_wr,'(a,f10.4)') 'Raw Ekin from vel(dt) = ',Ekin_dt
         write(imd_wr,'(a,f10.4)') 'Processed Ekin        = ',Ekin
         write(imd_wr,'(a,f10.4)') 'Processed tempcurr    = ',tempcurr
c
c Recover thermostat continuation data if it exists
c
 8108    format(a8)
c
c Find entry label for section, skip thermostat restart if do not find
         wantlbl = 'thermost'
         read(ifl,8108,err=200,end=200)  readlbl
         if( readlbl .ne. wantlbl ) goto 200
c
         write(imd_wr,'(/a)') 'QMD Restart: found thermostat data'
c
c Found thermostat continuation data: read it, kill job if corrupted
         wantlbl = 'itimeste'
         read(ifl,8108,err=113,end=113)  readlbl
         if( readlbl .ne. wantlbl ) goto 113
         read(ifl,*)  itime_restart
         itimestep = itime_restart
c
         wantlbl = 'temp    '
         read(ifl,8108,err=113,end=113)  readlbl
c        For backward compatibility
         if( readlbl .ne. wantlbl ) wantlbl = 'tempcurr'
         if( readlbl .ne. wantlbl ) goto 113
         read(ifl,*)  temp_restart
         tempprnt = temp_restart
c
         wantlbl = 'Ekin    '
         read(ifl,8108,err=113,end=113)  readlbl
         if( readlbl .ne. wantlbl ) goto 113
         read(ifl,*)  Ekin_restart
         Ekinprnt = Ekin_restart
c
         wantlbl = 'shv     '
         read(ifl,8108,err=113,end=113)  readlbl
         if( readlbl .ne. wantlbl ) goto 113
         read(ifl,*)  shv_restart
         shv = shv_restart
c
         wantlbl = 'end ther'
         read(ifl,8108,err=113,end=113)  readlbl
         if( readlbl .ne. wantlbl ) goto 113
c
         write(imd_wr,'(a,3f18.8)')
     $    'QMD Restart data: temp,Ekin,shv=',
     $     tempprnt,Ekinprnt,shv
         goto 200
c
  113    continue
c        We had an error processing the MD restart file
         write(imd_wr,*) 'Error processing MD restart file (lcao.vxyz)'
         call STOPXERR( 'Corrupted MD restart file' )
c
  200    continue
         call FLCLOSE( ifl )
c
c        May want to turn this printing off or flag with verbose
         write(imd_wr,*) "Following data input in QMDRestart:"
         do  iatm=1,md_natoms
            write(imd_wr,'(i4,7f10.5)') iatm,atom_mass(iatm),
     $           (ratm(j,iatm),j=1,3),
     $           (vel(j,iatm),j=1,3)
         enddo
c
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDVINIT0
c
c
      subroutine QMDVINIT0
c
c     Set velocities to reproduce the Boltzmann distribution. This
c     was the old method, which has been superceded by the version
c     in QMDVINIT
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c     Local storage
      dimension vrand(3),xsign(3),iseed(4)
c iseed=Seeds for the random number generator
c   - Probably a more intelligent set exists (1-3.lt.4096, 4 is odd)
      data iseed /4,3,2,1/
c
c     This sets the initial temperature to twice Ti; the reasoning behind
c     this is that if we set the velocities to twice the temp, if the
c     molecule is at equilibrium, we have on the average that temperature
c     in the molecule
c
c >>>> EXECUTABLE CODE:
c
      call QMDCONST( fconst, tconst )
      Ekin_per_mode = 0.5*tempinit*tconst
c
      do i = 1, md_natoms
         xm = atom_mass(i)
         vavg = dsqrt(2*Ekin_per_mode*fconst/xm)
c     Find a more intelligent setting for stdev!
         stdev = 0.1*vavg     
         call DLARNV(3,iseed,3,vrand)
         call DLARNV(2,iseed,3,xsign)
c     This is unnecessarily complex. The vx,vy,vz
c     components should be centered around +/-vavg. The call to
c     dlarnv(3) gives a normal distribution centered around vavg.
c     The calls to xsign and dsign simply toggle the +/- part.
c         vel(1,i) = sign(vavg,xsign(1))
c         vel(2,i) = sign(vavg,xsign(2))
c         vel(3,i) = sign(vavg,xsign(3))
         vel(1,i) = sign(vavg,xsign(1))+stdev*vrand(1)
         vel(2,i) = sign(vavg,xsign(2))+stdev*vrand(2)
         vel(3,i) = sign(vavg,xsign(3))+stdev*vrand(3)
         vprint(1,i) = zero
         vprint(2,i) = zero
         vprint(3,i) = zero
      enddo
c
c Remove net momentum of cell:
      call QMDSUBVCOM
c
      call GET_EKIN( md_natoms,atom_mass,vel, Ekin )
      call GET_TEMP_E( Ekin,md_nmodes, temp )
c
      factor = SQRT( tempinit/temp )
      call QMDVDAMP( factor )
c
      call GET_EKIN( md_natoms,atom_mass,vel, Ekin )
      call GET_TEMP_E( Ekin,md_nmodes, temp )
      Ekincurr = Ekin
      tempcurr = temp
      Ekinprnt = Ekin
      tempprnt = temp
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDVINIT2
c
c
      subroutine QMDVINIT2
c
c     Set velocities to 3/2 kT in random directions
c
c      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      IMPLICIT NONE
c
      DOUBLE PRECISION  zero,one,two
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      INTEGER  natmd
      PARAMETER  ( natmd =  400 )
      DOUBLE PRECISION  tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,   vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat
      INTEGER  md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,ieq_type
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
      DOUBLE PRECISION  fconst,tconst
c
c
c  Local declarations
      DOUBLE PRECISION  xm,vavg,Ekin,temp
      INTEGER  iseed
      INTEGER  iatm
      DIMENSION iseed(4)
      DATA  iseed /4,3,2,1/
c Probably a more intelligent set for iseed exists
      DOUBLE PRECISION  theta,phi
      DOUBLE PRECISION  pi,twopi
      DOUBLE PRECISION  Ekin_per_atom
c
c >>>> EXECUTABLE CODE:
c
      call QMDCONST( fconst, tconst )
      pi = DACOS( -1.0d0 )
      twopi = two*pi
c        
c     Set the initial temperature to 3/2 kT.
c     The factor of 2 accounts for the fact that half of this energy will
c     eventually go into potential modes.
      Ekin_per_atom = 2*1.5*tempinit*tconst 
c     
      do  iatm=1,md_natoms
         vavg = dsqrt(2d0*Ekin_per_atom*fconst/atom_mass(iatm))
c         
c        Does the following scheme bias toward equatorial vectors? -PAS
         call DLARUV( iseed,1, theta )
         theta = 2*pi*theta
         call DLARUV( iseed,1, phi )
         phi = pi*phi
c        
         vel(1,iatm) = vavg*dcos(theta)*dsin(phi)
         vel(2,iatm) = vavg*dsin(theta)*dsin(phi)
         vel(3,iatm) = vavg*dcos(phi)
         vprint(1,iatm) = zero
         vprint(2,iatm) = zero
         vprint(3,iatm) = zero
      enddo
c     
c      call QMDTEMPANAL( md_natoms,atom_mass,vel,fconst*tconst,
c     $     imd_wr)
c  Compute the raw temperature:
      call GET_TEMP_V( md_natoms,atom_mass,vel, Ekin,temp )
c  Do analysis of temperature over atom types:
      call TMONTEMP( -2 )
      write(imd_wr,*)   'VINIT: raw temp = ',temp
c
c  Remove net momentum:
      call QMDSUBVCOM
      call GET_TEMP_V( md_natoms,atom_mass,vel, Ekin,temp )
      Ekincurr = Ekin
      tempcurr = temp
      Ekinprnt = Ekin
      tempprnt = temp
c
      call TMONTEMP( -2 )
      write(imd_wr,*)  'VINIT: post subvcom temp = ',temp
c
c  Could, at this point, scale the velocities to get T, but since
c   this is initialization and we will miss target anyway, skip it.
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDVINIT
c
c
      subroutine QMDVINIT
c---------------------------------------------------------------
c Initialize atom velocities to random from Boltzmann distribution
c---------------------------------------------------------------
c
c      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      IMPLICIT NONE
c
      DOUBLE PRECISION  zero,one,two
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      integer natmd
      PARAMETER  ( natmd =  400 )
      DOUBLE PRECISION  tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,   vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat
      INTEGER  md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,ieq_type
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
      DOUBLE PRECISION  fconst,tconst
c
c  Local declarations
      DOUBLE PRECISION  pi,twopi
      INTEGER  iatm,i
      DOUBLE PRECISION  xm,vavg,v2mode,vmode, Ekin,temp
      DOUBLE PRECISION  Ekin_per_atom,Ekin_per_mode
cpas      DOUBLE PRECISION  DACOS,DSQRT
cpas      EXTERNAL  DACOS,DSQRT
c
c We use the DLARUV Lapack routine to generate a vector of random numbers
c from a uniform (0,1) distribution:
c  DLARUV( iseed, N, X )
c   iseed(4) - input/output: all elements 0:4095, iseed(4) must be odd
c   N - input: # random numbers to be generated in vector X(N), N.le.128
c   X(N) - output: vector of N random numbers
c  DLARNV( idist, iseed, N, Z ) offers different distributions:
c   idist=1/uniform(0:1), idist=2/uniform(-1:1), idist=3/normal(gaussian)
c We will use the A&T recommended Box-Muller method
c     u1,u2 from uniform distribution (0:1), to get normal distribution:
c     gdist = sqrt(-2*ln( u1 ) ) * cos( 2pi*u2 ) [also with sin(2pi*u2)]
c
      DOUBLE PRECISION  urndm, gdist
      INTEGER    iseed
      DIMENSION  urndm(2,3)
      DIMENSION  iseed(4)
c Probably could do a more intelligent set of input seeds than 4,3,2,1:
c At least select not so obviously trivial a set
      DATA  iseed / 1342, 731,3621, 331 /
c
c >>>> EXECUTABLE CODE:
c
      call QMDCONST( fconst, tconst )
      pi = DACOS( -1.0d0 )
      twopi = two*pi
c
c  This would set the initial temperature to 3/2 kT:
      Ekin_per_mode = 0.5d0*tempinit*tconst
      Ekin_per_atom = 1.5d0*tempinit*tconst
c  Add extra factor of 2: half this energy will seep into potential modes
c   (Note: this assumes that system is at/near equilibrium)
      Ekin_per_mode = two*Ekin_per_mode
      Ekin_per_atom = two*Ekin_per_atom
c
      write(imd_wr,*) 'Initialize MD velocities, Boltzmann distribution'
c
      do  iatm=1,md_natoms
c         vavg = DSQRT( 2d0*Ekin_per_atom*fconst/atom_mass(iatm) )
C         v2atom = two*Ekin_per_atom*fconst / atom_mass(iatm)
C         vatom = DSQRT( v2atom )
         v2mode = two*Ekin_per_mode*fconst / atom_mass(iatm)
         vmode = DSQRT( v2mode )
c
         call DLARUV( iseed, 6, urndm )
         do  i=1,3
           gdist = SQRT(-two*LOG(urndm(1,i)))*DCOS(twopi*urndm(2,i))
cccc          Should truncate this thing to prevent big excursions?
ccc           if( ABS(gdist) .gt. two ) gdist = two*gdist/ABS(gdist)
           vel(i,iatm) = vmode*gdist
         enddo
         vprint(1,iatm) = zero
         vprint(2,iatm) = zero
         vprint(3,iatm) = zero
      enddo
c
c      call QMDTEMPANAL( md_natoms,atom_mass,vel,fconst*tconst,
c     $     imd_wr )
c
c  Compute the temperature of the raw system:
      call GET_TEMP_V( md_natoms,atom_mass,vel, Ekin,temp )
c  Do analysis of temperature over atom types:
      call TMONTEMP( -2 )
      write(imd_wr,*)   'VINIT: raw temp = ',temp
c
c  Remove net momentum:
      call QMDSUBVCOM
      call GET_TEMP_V( md_natoms,atom_mass,vel, Ekin,temp )
      Ekincurr = Ekin
      tempcurr = temp
      Ekinprnt = Ekin
      tempprnt = temp
c
      call TMONTEMP( -1 )
      write(imd_wr,*)  'VINIT: post subvcom temp = ',temp
c
c  Could, at this point, scale the momenta to get T, but since
c   this is initialization and we will miss target anyway, skip it.
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDVDAMP
c
c
      subroutine QMDVDAMP( velfac )
c
c     Rescale the velocities by a desired factor
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
      do  iatm=1,md_natoms
         vel(1,iatm) = vel(1,iatm)*velfac
         vel(2,iatm) = vel(2,iatm)*velfac
         vel(3,iatm) = vel(3,iatm)*velfac
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GET_VCOM
c
c
      subroutine GET_VCOM( nat,atom_mass,vel, vcomx,vcomy,vcomz )
c
c     Compute the velocity of the center of mass
c
c     nat        The number of atoms
c     atom_mass  The mass of each atom
c     vel        The velocity of each atom
c     vcomx      The CoM velocity, x-direction
c     vcomy      The CoM velocity, y-direction
c     vcomz      The CoM velocity, z-direction
c
c     Rick Muller, September 2007
c
      IMPLICIT NONE
c
      DOUBLE PRECISION  zero,one,two
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
      INTEGER  nat
      DOUBLE PRECISION  atom_mass(nat),vel(3,nat),vcomx,vcomy,vcomz
c
      INTEGER  iatm
      DOUBLE PRECISION  xm,totmass,totatm
c
      vcomx = zero
      vcomy = zero
      vcomz = zero
      totmass = zero
      do  iatm=1,nat
         xm =atom_mass(iatm)
         totmass = totmass + xm
         vcomx = vcomx + atom_mass(iatm)*vel(1,iatm)
         vcomy = vcomy + atom_mass(iatm)*vel(2,iatm)
         vcomz = vcomz + atom_mass(iatm)*vel(3,iatm)
      enddo
c
      totatm = DBLE( nat )
c  Now divide by total mass to get net velocity of CoM
      vcomx = vcomx / totmass
      vcomy = vcomy / totmass
      vcomz = vcomz / totmass
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDSUBVCOM
c
c
      subroutine QMDSUBVCOM
c---------------------------------------------------------------
c Removes center-f-mass velocity from the cell
c---------------------------------------------------------------
c
c Make sure there is no net momentum in the cell
c Remove center-of-mass velocitu from all atoms to achieve this
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c >>>> EXECUTABLE CODE
c
      call GET_VCOM( md_natoms,atom_mass,vel, vcomx,vcomy,vcomz )
c
      write(imd_wr,'(a,1p3d12.4)') 'QMDSUBVCOM: Removing velocity_CM=',
     $ vcomx,vcomy,vcomz
c
      do  iatm=1,md_natoms
         xm = atom_mass(iatm)
         vel(1,iatm) = vel(1,iatm) - vcomx
         vel(2,iatm) = vel(2,iatm) - vcomy
         vel(3,iatm) = vel(3,iatm) - vcomz
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDPRINT
c
c
      subroutine QMDPRINT( ratm )
c     Print a snapshot of the current MD run
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
      COMMON /savepress/ presssave,prgpasave
c
c     Argument arrays
      dimension ratm(3,natmd),fatm(3,natmd)
c     Local variables
      character*2 sym
c
c >>>> EXECUTABLE CODE:
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  Record current MD status to the standard output file
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (itimestep .le. neqtsteps) then
         write(imd_wr,*) '<QMD0>'
      else
         write(imd_wr,*) '<QMD>'
      endif
c
      write(imd_wr,'(A,i8)') 'Nat:           ',md_natoms
c
      write(imd_wr,'(A,3F14.5)')
     $  'MDtemperature(K): target, T(t), T(t+1/2)=',
     $                     tempinit,tempprnt,tempcurr
c
      write(imd_wr,'(A,3F14.5)')
     $  'MDenergy(Ry;t) kin+pot=tot: ',
     $   Ekinprnt, Epot, Ekinprnt+Epot
c
      write(imd_wr,'(A,F12.5)') 'Pressure (GPa) ',prgpasave
c
      write(imd_wr,'(A,i8,F12.5,i8,i8)')
     $  'MDtime(fs) # dt tot eq: ',itimestep,timestep,ntsteps,neqtsteps
c
c     Consider putting in units for velocities (bohr/fs)??
      call GET_VCOM( md_natoms,atom_mass,vel, vcomx,vcomy,vcomz )
      write(imd_wr,'(A,F12.5,F12.5,F12.5)')
     $  'Vcom x y z:    ',vcomx,vcomy,vcomz
c
      write(imd_wr,*) 'Positions (bohr) and velocities (bohr/fs) at t'
      write(imd_wr,*) '#  mass  x   y   z   vx   vy   vz'
      do i = 1, md_natoms
        call round(atom_mass(i),imass)
        write(imd_wr,'(2i3,6f12.8)') i,imass,
     $        (ratm(j,i),j=1,3),
     $        (vprint(j,i),j=1,3)
      enddo
c
      write(imd_wr,*) 'Restart information'
      write(imd_wr,*) md_natoms
      write(imd_wr,'(2A,i8)') 'Quest pos (bohr), at t ',
     $     ' vel (bohr/fs) at t+1/2, for t = ',
     $     itimestep
      do i = 1, md_natoms
        call round(atom_mass(i),imass)
        call mass2sym(imass,sym)
        write(imd_wr,'(A2,6f12.8)') sym,(ratm(j,i),j=1,3),
     $        (vel(j,i),j=1,3)
      enddo
c
      write(imd_wr,'(a)')  'shv      - current Hoover parameter'
      write(imd_wr,'(2x,f18.8)')   shv 
c
      if (itimestep .le. neqtsteps) then
         write(imd_wr,*) '</QMD0>'
      else
         write(imd_wr,*) '</QMD>'
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     Dump the current snapshot to the lcao.vxyz file
c * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      call FLOPEN( ifl,'lcao.vxyz','UNK','FOR',ierror )
      if (ierror .ne. 0)
     $     call STOPXERR( 'QMDRESTART: error opening vxyz' )
c
      write(ifl,*) md_natoms
      write(ifl,'(2A,i8)') 'Quest pos (bohr), at t ',
     $     ' vel (bohr/fs) at t+1/2, for t = ',
     $     itimestep
      do i = 1, md_natoms
        call round(atom_mass(i),imass)
        call mass2sym(imass,sym)
        write(ifl,'(A2,6f12.8)') sym,(ratm(j,i),j=1,3),
     $        (vel(j,i),j=1,3)
      enddo
c
c Write out some information for clean thermostat continuation:
      write(ifl,'(a)')  'thermostat continuation data'
      write(ifl,'(a)')  'itimestep - index of current step in MD'
      write(ifl,'(2x,i12)')   itimestep
      write(ifl,'(a)')  'temp     - in-synch temp(t), not temp(t+dt)'
      write(ifl,'(2x,f18.8)')   tempprnt
      write(ifl,'(a)')  'Ekin     - in-synch Ekin(t), not Ekin(t+dt)'
      write(ifl,'(2x,f18.8)')   Ekinprnt
      write(ifl,'(a)')  'shv      - current Hoover parameter'
      write(ifl,'(2x,f18.8)')   shv 
      write(ifl,'(a)')  'end thermostat data'
c
      call FLCLOSE( ifl )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDTEMPANAL
c
c
      subroutine QMDTEMPANAL(nat,xm,v,xk,imd_wr)
      implicit none
      integer nat,i,imass,imd_wr
      character*2 sym
      double precision xm(nat),v(3,nat),xk,v2
c
      write(imd_wr,*) 'Temperature analysis'
      do i = 1, nat
        call round(xm(i),imass)
        call mass2sym(imass,sym)
        v2 = v(1,i)**2 + v(2,i)**2 + v(3,i)**2
        write(imd_wr,'(i4,x,a4,x,7f10.4)') i,sym,v2,xm(i)*v2/(3*xk)
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDUPDATE
c
c
      subroutine QMDUPDATE( ratm, fatm )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c     Argument arrays
      dimension ratm(3,natmd),fatm(3,natmd)
c The following are MD method id tags; they must match QMDTYPSET
      DATA imd_nve   ,imd_berend,imd_hoover,imd_tscale,imd_damped
     $    /    1     ,    2     ,    3     ,    4     ,    5     /
      DATA imd_anders,imd_hugoni,imd_langvn
     $    /    6     ,    7     ,    8     /
c
c >>>> EXECUTABLE CODE:
c
c     RPM: temporarily suspended, so see if this causes the 
c     temp difference between O and H. 9/4/2007
      call QMDSUBFCOM( fatm )
c
c Default to production MD settings
      tempinit = temp_md
      imdtype = md_type
      parstat = param_md
      if( itimestep .lt. neqtsteps )then
c       Equilibration phase - use equilibration MD settings
        write(imd_wr,*) 'MD-update: equilibration phase'
        tempinit = temp_eq
        imdtype = ieq_type
        parstat = param_eq
      endif
      if( itimestep .eq. neqtsteps )then
c       Turnover of thermostats, reset temperature monitor
        call TMONRESET( )
        shv = zero
      endif
c
      if( imdtype .eq. imd_nve    )then 
c       NVE/Microcanonical update, leapfrog-verlet
         write(imd_wr,*) 'MD-update: NVE/leapfrog'
         call QMDUPDATE_NVE_LF( ratm, fatm )
c
      elseif( imdtype .eq. imd_langvn )then 
c       NVT/Langevin thermostat, leapfrog-verlet
         write(imd_wr,*) 'MD-update: NVT/Langevin/leapfrog'
         call QMDUPDATE_LANGVN( ratm, fatm )
c
      elseif( imdtype .eq. imd_berend )then
c       NVT/Berendsen thermostat, leapfrog-verlet
         write(imd_wr,*) 'MD-update: NVT/Berendsen/leapfrog'
         call QMDUPDATE_BEREND( ratm, fatm )
c
      elseif( imdtype .eq. imd_hoover )then
c       NVT/Nose-Hoover thermostat, leapfrom-verlet
         write(imd_wr,*) 'MD-update: NVT/Hoover/leapfrog'
         call QMDUPDATE_HOOVER( ratm, fatm )
c
      elseif( imdtype .eq. imd_tscale )then
c       Temperature scaling, leapfrog-verlet
         write(imd_wr,*) 'MD-update: tscale/leapfrog'
         call QMDUPDATE_TSCALE( ratm, fatm )
c
      elseif( imdtype .eq. imd_damped )then
c       Damped dynamics, leapfrog-verlet
         write(imd_wr,*) 'MD-update: damped/leapfrog'
         call QMDUPDATE_DAMPED( ratm, fatm )
c
      else
         write(imd_wr,*) 'MD Type = ',imdtype,' not currently supported'
         call STOPXERR('QMDUPDATE - invalid MD type')
      endif
cdbg:start
c     RPM 9/17/2007
c     Remove the net momentum, if any.
c     This call should technically not be needed, and, if it fixes
c     things, we must investigate why it is necessary before checking
c     this into the main code
      call QMDSUBVCOM
cdbg:end
c
      itimestep = itimestep + 1
      maxtsteps = ntsteps + neqtsteps
c
c     Need to exit more gracefully (will let main code handle exit)
      if( itimestep .ge. maxtsteps )then
c        Do not call QMD a final time, here --- not in proper units!
C         call QMDPRINT( ratm )
         write(imd_wr,*) '>>>>> MD run completed at timestep ',itimestep
c        Do not stop here, pop up, let main handle final exit,
c         and let main handle the update and final vxyz properly
C         call STOPXOK( 'QMD Complete' )
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDSUBFCOM
c
c
      subroutine QMDSUBFCOM( fatm )
c---------------------------------------------------------------
c     Purify the forces to remove net forces on the Unit Cell
c     arising from numerical errors in force evaluation.
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c     Argument arrays
      dimension fatm(3,natmd)
c
      ftotx = zero
      ftoty = zero
      ftotz = zero
c
c     Compute the total net force present in the cell:
      do i = 1, md_natoms
         ftotx = ftotx + fatm(1,i)
         ftoty = ftoty + fatm(2,i)
         ftotz = ftotz + fatm(3,i)
      enddo
c
      totatom = DBLE( md_natoms )
      ftotx = ftotx / totatom
      ftoty = ftoty / totatom
      ftotz = ftotz / totatom
c
      write(imd_wr,'(a,1p3d12.4)')
     $ 'QMDSUBFCOM: Removing net cell force/atom=', ftotx,ftoty,ftotz
c
c     Subtract this from force on each atom
      do i = 1, md_natoms
         fatm(1,i) = fatm(1,i) - ftotx
         fatm(2,i) = fatm(2,i) - ftoty
         fatm(3,i) = fatm(3,i) - ftotz
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDUPDATE_NVE_LF
c
c
      subroutine QMDUPDATE_NVE_LF( ratm, fatm )
c---------------------------------------------------------------
c Purpose: do NVE (microcanonical) update, leapfrog verlet
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c     Argument arrays
      DIMENSION  ratm(3,natmd),fatm(3,natmd)
c
c >>>> EXECUTABLE CODE:
c
      call QMDCONST( fconst, tconst )
c
      call TMONPREP( )
c
      Ekin_t = zero
      Ekin_n = zero
      Ekin_at = zero
      Ekin_an = zero
c
      do  i=1,md_natoms
        xm = atom_mass(i)
        ax = fatm(1,i)*fconst/xm
        ay = fatm(2,i)*fconst/xm
        az = fatm(3,i)*fconst/xm
c
c       Advance the velocities and compute Ekin
c
c       Compute the off-synch velocity v(t+dt) for next step
        vnewx = vel(1,i) + timestep*ax
        vnewy = vel(2,i) + timestep*ay
        vnewz = vel(3,i) + timestep*az
        vn2 = vnewx*vnewx + vnewy*vnewy + vnewz*vnewz
        ekinn_atm = 0.5d0*xm*vn2/fconst
        Ekin_n = Ekin_n + ekinn_atm
c
c       But which velocity to use in KE?  v(t+d), v(t)?
c       Maybe neither.  Time-average v^2 around t where a(t) computed?
c       Could consider doing t-ave around vnew, the vel that will 
c       be used in the next update, but then a(t) less valid there.
c
c       Compute the in-synch velocity v(t)
        vtx = vel(1,i) + timestep*ax/2d0
        vty = vel(2,i) + timestep*ay/2d0
        vtz = vel(3,i) + timestep*az/2d0
        vt2 = vtx*vtx + vty*vty + vtz*vtz
        ekint_atm = 0.5d0*xm*vt2/fconst
        Ekin_t = Ekin_t + ekint_atm
c
c       KE averaged over this segment of the trajectory:
c       INT(-dt/2,+dt/2) (v+au)^2 du / dt = v^2+a^2dt^2/12.0
c
        accel2 = ax*ax + ay*ay + az*az
        dt2 = timestep*timestep
c
c        ... centered around in-synch position:
        vt2ave = (vt2 + accel2*dt2/12.0d0)
        akint_atm = 0.5d0 * xm*vt2ave/fconst
        Ekin_at = Ekin_at + akint_atm
c
c        ... centered around off-synch velocity:
        vn2ave = (vn2 + accel2*dt2/12.0d0)
        akinn_atm = 0.5d0 * xm*vn2ave/fconst
        Ekin_an = Ekin_an + akinn_atm
c
        call TMONATOM( i   ,ekint_atm,ekinn_atm,akint_atm,akinn_atm )
c
        vprint(1,i) = vtx
        vprint(2,i) = vty
        vprint(3,i) = vtz
c
c       Advance the positions
        ratm(1,i) = ratm(1,i) + timestep*vnewx
        ratm(2,i) = ratm(2,i) + timestep*vnewy
        ratm(3,i) = ratm(3,i) + timestep*vnewz
c
c       Save the new velocities
        vel(1,i) = vnewx
        vel(2,i) = vnewy
        vel(3,i) = vnewz
c
      enddo
c
c  Update the current temperature
      Ekincurr = Ekin_n
      call GET_TEMP_E( Ekincurr,md_nmodes, temp )
      tempcurr = temp
c  Update the in-synch temperature
      Ekinprnt = Ekin_t
      call GET_TEMP_E( Ekinprnt,md_nmodes, temp )
      tempprnt = temp
c
      write(imd_wr,'(a,2f12.3)') 'DBG/NVE: tempts(t,t+1/2)=',
     $ tempprnt,tempcurr
c
      call TMONSTEP( 'NVE' )
c
      Ekinmon = Ekinprnt
      Epotmon = Epot
c Emdmon should contain the heat bath energy, which is zero for NVE
      Emdmon = zero
      call TMONENGY( Ekinmon,Epotmon,Emdmon )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDUPDATE_TSCALE
c
c
      subroutine QMDUPDATE_TSCALE( ratm, fatm )
c---------------------------------------------------------------
c Purpose: NVE (microcanonical), leapfrog verlet update
c
c     Version of MD that does a Leapfrog step and then rescales
c     the velocities to match the temperature.
c
c     WARNING: This type of dynamics puts more energy into the 
c     heavy atoms than the light ones. Use at your own risk!
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c     Argument arrays
      DIMENSION  ratm(3,natmd),fatm(3,natmd)
c
c Local declarations:
      DATA  itypscale / 1 /
      DATA  lstoutmd / 2 /
c
c >>>> EXECUTABLE CODE:
c
c First, start with an NVE update ...
      call QMDUPDATE_NVE_LF( ratm, fatm )
c
c  ... and then scale the resulting velocities to the temperature
c
c  The simple scaling rule appears to work best (how so?)
c
      if( itypscale .eq. 1 )then
c       Simple Tscale: Tnew = tscale*Told + (1-tscale)*T
        tscale = parstat
        if( tscale .lt. zero .or. tscale .gt. one )then
          write(imd_wr,*) 'Simple Tscale factor must be (0,1):',tscale
          call STOPXERR( 'out-of-range factor for simple Tscale' )
        endif
C        velfac = SQRT( tempinit/tempcurr )
c       Let us try something with a parameter that makes this softer
        temp = tempinit + tscale*(tempcurr-tempinit)
        tfac = temp / tempcurr
        velfac = SQRT( tfac )
c
      if( lstoutmd .gt. 1 )
     $  write(imd_wr,'(a,f12.8,3f12.3)')
     $  'MD-Tscale: tscale, T,Told,Tnew=',
     $   tscale, tempinit,tempcurr,temp
c
      elseif( itypscale .eq. 2 )then
c       Alternate Tscale
c       (Berendsen, Postma, van Gunsteren, DiNola and Haak, 1984)
        tscale = 200d0
        tscale = parstat
        write(imd_wr,*) 'TSCALE: time scaling parameter=',tscale
        velfac = SQRT( one + (timestep/tscale)*(tempinit/tempcurr-one) )
c
      else
c       Invalid Tscale type
        call STOPXERR( 'invalid Tscale method' )
      endif
c
c  Now actually scale the velocities!
      call QMDVDAMP( velfac )
c
c  And get the output temperature
      call GET_EKIN( md_natoms,atom_mass,vel, Ekin )
      call GET_TEMP_E( Ekin,md_nmodes, temp )
      Ekincurr = Ekin
      tempcurr = temp
      Ekinprnt = Ekin
      tempprnt = temp
c
      call TMONSTEP( 'NVE' )
      Ekinmon = Ekinprnt
      Epotmon = Epot
c Emdmon should contain the heat bath energy, but not done yet
      Emdmon = zero
      call TMONENGY( Ekinmon,Epotmon,Emdmon )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDUPDATE_LANGVN
c
c
      subroutine QMDUPDATE_LANGVN( ratm, fatm )
c---------------------------------------------------------------
c Purpose: Langevin thermostat, leapfrom verlet update
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c     Argument arrays
      DIMENSION  ratm(3,natmd),fatm(3,natmd)
c Local declarations:
      DATA  lstoutmd / 2 /
c
c >>>> EXECUTABLE CODE:
c
      call STOPXERR( 'LANGVN thermostat not yet implemented' )
c
      call TMONPREP( )
      t_langvn = parstat
c
      if( lstoutmd .gt. 1 )
     $write(imd_wr,'(a,f12.3)') 'MD-Langvn: t_langvn(fs)=',t_langvn
c 
      call QMDCONST( fconst, tconst )
c
      call QMDUPDATE_NVE_LF( ratm, fatm )
c
      call GET_EKIN( md_natoms,atom_mass,vel, Ekin )
      call GET_TEMP_E( Ekin,md_nmodes, temp )
      Ekincurr = Ekin
      tempcurr = temp
      Ekinprnt = Ekin
      tempprnt = temp
c
      call TMONSTEP( 'LVN' )
c
      Ekinmon = Ekinprnt
      Epotmon = Epot
c Emdmon should contain the heat bath energy, but not done yet
      Emdmon = zero
      call TMONENGY( Ekinmon,Epotmon,Emdmon )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDUPDATE_DAMPED
c
c
      subroutine QMDUPDATE_DAMPED( ratm, fatm )
c---------------------------------------------------------------
c Purpose: damped dynamics, leapfrog verlet update
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c     Argument arrays
      DIMENSION  ratm(3,natmd),fatm(3,natmd)
c Local declarations:
      DATA  lstoutmd / 2 /
c
c >>>> EXECUTABLE CODE:
c
c     Take damping factor from current thermostat definition
      damping = parstat
c
      if( lstoutmd .gt. 1 )
     $write(imd_wr,'(a,f12.8)') 'MD-Damped: damping factor=',damping
c
      if( damping .gt. one .or. damping .lt. zero )then
        write(imd_wr,*) 'MD damping factor is invalid=',damping
        call STOPXERR( 'QMD damping factor is out-of-range' )
      endif
c
      call QMDUPDATE_NVE_LF( ratm, fatm )
      call QMDVDAMP( damping )
c
      call GET_EKIN( md_natoms,atom_mass,vel, Ekin )
      call GET_TEMP_E( Ekin,md_nmodes, temp )
      Ekincurr = Ekin
      tempcurr = temp
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDUPDATE_BEREND
c
c
      subroutine QMDUPDATE_BEREND( ratm, fatm )
c---------------------------------------------------------------
c Purpose: Berendsen thermostart, leapfrog verlet update
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c     Argument arrays
      DIMENSION  ratm(3,natmd),fatm(3,natmd)
c
c Local declarations:
      DATA  lstoutmd / 3 /
c
c >>>> EXECUTABLE CODE:
c
      call QMDCONST( fconst, tconst )
c Berendsen friction constant: t_berend
c  t_berend=time (fs) constant to scale connection to the friction
c  Allen&Tildesley quote Berendsen: t_b=0.4ps=400fs for water
c  The O-H stretch is ~10fs, hence this is over 40 vib periods.
c  t_berend has been converted into a user-settable input parameter
      t_berend = parstat
c
      if( lstoutmd .gt. 1 )
     $write(imd_wr,*) 'MD-BERENDSEN: t_berend(fs)=',t_berend
c
      s = ( one / (two*t_berend) ) * (tempcurr-tempinit)/tempcurr
c
      if( lstoutmd .gt. 3 ) write(imd_wr,'(a,i8,4f14.6)')
     $ 'BEREND/in:',itimestep,tempinit,tempcurr,Ekin,s
c
      call TMONPREP( )
c
      Ekin_t = zero
      Ekin_n = zero
      Ekin_at = zero
      Ekin_an = zero
c
      do  i=1,md_natoms
        xm = atom_mass(i)
        ax = fatm(1,i)*fconst/xm - s*vel(1,i)
        ay = fatm(2,i)*fconst/xm - s*vel(2,i)
        az = fatm(3,i)*fconst/xm - s*vel(3,i)
c
c       Advance the velocities and compute Ekin
c       
c       Compute the off-synch velocity v(t+dt) for next step
        vnewx = vel(1,i) + timestep*ax 
        vnewy = vel(2,i) + timestep*ay
        vnewz = vel(3,i) + timestep*az
        vn2 = vnewx*vnewx + vnewy*vnewy + vnewz*vnewz
        ekinn_atm = 0.5d0*xm*vn2/fconst
        Ekin_n = Ekin_n + ekinn_atm
c       
c       But which velocity to use in KE?  v(t+d), v(t)?
c       Maybe neither.  Time-average v^2 around t where a(t) computed?
c       Could consider doing t-ave around vnew, the vel that will
c       be used in the next update, but then a(t) less valid there.
c       
c       Compute the in-synch velocity v(t)
        vtx = vel(1,i) + timestep*ax/2d0
        vty = vel(2,i) + timestep*ay/2d0 
        vtz = vel(3,i) + timestep*az/2d0 
        vt2 = vtx*vtx + vty*vty + vtz*vtz
        ekint_atm = 0.5d0*xm*vt2/fconst
        Ekin_t = Ekin_t + ekint_atm
c       
c       KE averaged over this segment of the trajectory:
c       INT(-dt/2,+dt/2) (v+au)^2 du / dt = v^2+a^2dt^2/12.0
c
        accel2 = ax*ax + ay*ay + az*az
        dt2 = timestep*timestep
c     
c        ... centered around in-synch position:
        vt2ave = (vt2 + accel2*dt2/12.0d0)
        akint_atm = 0.5d0 * xm*vt2ave/fconst
        Ekin_at = Ekin_at + akint_atm
c
c        ... centered around off-synch velocity:
        vn2ave = (vn2 + accel2*dt2/12.0d0)
        akinn_atm = 0.5d0 * xm*vn2ave/fconst
        Ekin_an = Ekin_an + akinn_atm
c
        call TMONATOM( i   ,ekint_atm,ekinn_atm,akint_atm,akinn_atm )
c
        vprint(1,i) = vtx
        vprint(2,i) = vty
        vprint(3,i) = vtz
c       
c       Advance the positions
        ratm(1,i) = ratm(1,i) + timestep*vnewx
        ratm(2,i) = ratm(2,i) + timestep*vnewy
        ratm(3,i) = ratm(3,i) + timestep*vnewz
c       
c       Save the new velocities
        vel(1,i) = vnewx
        vel(2,i) = vnewy
        vel(3,i) = vnewz
c
      enddo
c
c  Update the current temperature
      Ekincurr = Ekin_n
      call GET_TEMP_E( Ekincurr,md_nmodes, temp )
      tempcurr = temp
c  Update the in-synch temperature
      Ekinprnt = Ekin_t
      call GET_TEMP_E( Ekinprnt,md_nmodes, temp )
      tempprnt = temp
c
      if( lstoutmd .gt. 4 )
     $write(imd_wr,'(a,2f12.3)') 'DBG/BER: tempts(t,t+1/2)=',
     $ tempprnt,tempcurr
c
      call TMONSTEP( 'BER' )
c
      Ekinmon = Ekinprnt
      Epotmon = Epot
c Emdmon should contain the heat bath energy, but not done yet
      Emdmon = zero
      call TMONENGY( Ekinmon,Epotmon,Emdmon )
c
      if( lstoutmd .gt. 3 ) write(imd_wr,'(a,4f12.6)')
     $ 'BEREND/out:',tempinit,tempcurr,Ekin
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDUPDATE_HOOVER
c
c
      subroutine QMDUPDATE_HOOVER( ratm, fatm )
c---------------------------------------------------------------
c Purpose: Nose/Hoover thermostart, leapfrog verlet update
c---------------------------------------------------------------
c
c This does a single-pass, non-self-consistent update.
c I.e., the temperature/velocities used as input to the Hoover are
c   from the previous step, T(t-1/2), but we actually use it as T(t),
c   to modify velocities/temperatures from t-1/2 to t+1/2 (leapfrog).
c However, with a short time step and a longer friction time constant,
c   where we are integrating over a large number of short steps,
c   this should wash out.
c But, it is something to think about how to do better in the future.
c   - 9Feb08-PAS
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
c     Argument arrays
      DIMENSION  ratm(3,natmd),fatm(3,natmd)
c Local declarations
      DATA  lstoutmd / 3 /
c
c >>>> EXECUTABLE CODE:
c
      call QMDCONST( fconst, tconst )
c Hoover friction constant
c  In A&T, X is generated as dX/dt=(f/Q) (Tcurr-Tinit)
c    f=# modes; Q=thermal inertia, but no guidance on how to set.
c    RPM had used Q=-2.d6, f=md_nmodes
c      Q = -2d6
c      xnu = md_nmodes/Q
c  Can also define (like here) as (1/t_h^2)*(Tcurr-Tinit)/Tcurr
c    where t_h is a time constant (for Troom water 20-100fs seems good)
c  This has been turned into a user-settable input parameter
      t_hoover = parstat
      x_hoover = one / ( tempinit * t_hoover**2 )
      etadot = x_hoover * ( tempcurr - tempinit )
c
      shv = shv + timestep*etadot
c
      if( lstoutmd .gt. 1 )then
        write(imd_wr,9010) 'MD-Hoover: t_hoover(fs)=',t_hoover
      endif
      if( lstoutmd .gt. 3 )then
        write(imd_wr,9010) 'MD-Hoover: xnu   =',x_hoover
        write(imd_wr,9010) 'MD-Hoover: etadot=',etadot
      endif
      if( lstoutmd .gt. 2 )then
        write(imd_wr,9011) 'HoovComp(T,Tcurr,shv):',
     $   itimestep,tempinit,tempcurr,shv
      endif
 9010 format(1x,a,1pd15.8)
 9011 format(1x,a,1x,i7,1x,2f12.3,1x,1pd15.8)
c
      call TMONPREP( )
c
      Ekin_t = zero
      Ekin_n = zero
      Ekin_at = zero
      Ekin_an = zero
c
      do  i=1,md_natoms
        xm = atom_mass(i)
        ax = fatm(1,i)*fconst/xm - shv*vel(1,i)
        ay = fatm(2,i)*fconst/xm - shv*vel(2,i)
        az = fatm(3,i)*fconst/xm - shv*vel(3,i)
c
c       Advance the velocities and compute Ekin
c       
c       Compute the off-synch velocity v(t+dt) for next step
        vnewx = vel(1,i) + timestep*ax 
        vnewy = vel(2,i) + timestep*ay
        vnewz = vel(3,i) + timestep*az
        vn2 = vnewx*vnewx + vnewy*vnewy + vnewz*vnewz
        ekinn_atm = 0.5d0*xm*vn2/fconst
        Ekin_n = Ekin_n + ekinn_atm
c       
c       But which velocity to use in KE?  v(t+d), v(t)?
c       Maybe neither.  Time-average v^2 around t where a(t) computed?
c       Could consider doing t-ave around vnew, the vel that will
c       be used in the next update, but then a(t) less valid there.
c       
c       Compute the in-synch velocity v(t)
        vtx = vel(1,i) + timestep*ax/2d0
        vty = vel(2,i) + timestep*ay/2d0 
        vtz = vel(3,i) + timestep*az/2d0 
        vt2 = vtx*vtx + vty*vty + vtz*vtz
        ekint_atm = 0.5d0*xm*vt2/fconst
        Ekin_t = Ekin_t + ekint_atm
c       
c       KE averaged over this segment of the trajectory:
c       INT(-dt/2,+dt/2) (v+au)^2 du / dt = v^2+a^2dt^2/12.0
c
        accel2 = ax*ax + ay*ay + az*az
        dt2 = timestep*timestep
c     
c        ... centered around in-synch position:
        vt2ave = (vt2 + accel2*dt2/12.0d0)
        akint_atm = 0.5d0 * xm*vt2ave/fconst
        Ekin_at = Ekin_at + akint_atm
c
c        ... centered around off-synch velocity:
        vn2ave = (vn2 + accel2*dt2/12.0d0)
        akinn_atm = 0.5d0 * xm*vn2ave/fconst
        Ekin_an = Ekin_an + akinn_atm
c
        call TMONATOM( i   ,ekint_atm,ekinn_atm,akint_atm,akinn_atm )
c
        vprint(1,i) = vtx
        vprint(2,i) = vty
        vprint(3,i) = vtz
c       
c       Advance the positions
        ratm(1,i) = ratm(1,i) + timestep*vnewx
        ratm(2,i) = ratm(2,i) + timestep*vnewy
        ratm(3,i) = ratm(3,i) + timestep*vnewz
c       
c       Save the new velocities
        vel(1,i) = vnewx
        vel(2,i) = vnewy
        vel(3,i) = vnewz
c
      enddo
c
c  Update the current temperature
      Ekincurr = Ekin_n
      call GET_TEMP_E( Ekincurr,md_nmodes, temp )
      tempcurr = temp
c  Update the in-synch temperature
      Ekinprnt = Ekin_t
      call GET_TEMP_E( Ekinprnt,md_nmodes, temp )
      tempprnt = temp
c
      if( lstoutmd .gt. 4 )
     $write(imd_wr,'(a,2f12.3)') 'DBG/HVR: tempts(t,t+1/2)=',
     $ tempprnt,tempcurr
c
      call TMONSTEP( 'HVR' )
      Ekinmon = Ekinprnt
      Epotmon = Epot
c Emdmon should contain the heat bath energy, but not done yet
      Emdmon = zero
      call TMONENGY( Ekinmon,Epotmon,Emdmon )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDCONVCHK
c
c
      subroutine QMDCONVCHK( igstat )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
      maxtsteps = ntsteps + neqtsteps
c
c itimestep is zero-indexed. Increment for proper stop-check
c
      if( (itimestep+1) .lt. maxtsteps )then
c       Set igstat to 0=continue
        igstat = 0
      else
c       Set to 1=out-of-steps, not 2="converged"
        igstat = 1
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GET_TEMP_V
c
c
      subroutine GET_TEMP_V( nat,atom_mass,vel, Ekin,temp)
c
c     A function that computes the kinetic energy from
c     the number of atoms, the mass of each atom, 
c     and the velocities of each atom
c
c     nat        The number of atoms
c     atom_mass  The mass of each atom in atomic units
c     vel        The velocity of each atom in bohr/fs
c     ekin[out]  The kinetic energy in Ry
c     temp[out]  The temperature in Kelvin
c
c     WARNING: for simplicity, this routine assumes that 
c     md_nmodes = 3*nat.
c
c     Rick Muller, September 2008
c
      IMPLICIT NONE
c
      DOUBLE PRECISION  zero,one,two
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
      integer nat
      double precision atom_mass(nat),vel(3,nat),v2,ekin,
     $     temp
c
      call GET_EKIN( nat,atom_mass,vel, ekin )
      call GET_TEMP_E( ekin,3*nat, temp )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GET_EKIN
c
c
      subroutine GET_EKIN( nat,atom_mass,vel, ekin )
c
c     A function that computes the kinetic energy from
c     the number of atoms, the mass of each atom, 
c     and the velocities of each atom
c
c     nat        The number of atoms
c     atom_mass  The mass of each atom in atomic units
c     vel        The velocity of each atom in bohr/fs
c     ekin[out]  The kinetic energy in Ry
c
c     Rick Muller, September 2008
c
      IMPLICIT NONE
c
      DOUBLE PRECISION  zero,one,two
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
      double precision fconst,tconst
      integer nat,i
      double precision atom_mass(nat),vel(3,nat),v2,ekin
      double precision ddot
      external ddot
c
      call QMDCONST( fconst, tconst )
      ekin = zero
      do i=1,nat
         v2 = ddot(3,vel(1,i),1,vel(1,i),1)
         ekin = ekin + 0.5d0*atom_mass(i)*v2/fconst
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GET_TEMP_E
c
c
      subroutine GET_TEMP_E( ekin,nmodes, temp )
c
c     A function that computes the temperature in K from
c     the kinetic energy per mode in Ry
c
c     Rick Muller, September 2008
c
      IMPLICIT NONE
c
      DOUBLE PRECISION  zero,one,two
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
      double precision fconst,tconst
      integer nmodes
      double precision ekin,temp
c
      call QMDCONST( fconst, tconst )
      temp = 2d0*ekin/(tconst*dble(nmodes))
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDEPSET
c
c
      subroutine QMDEPSET(epotin)
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c
      Epot = epotin
      return
      end
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ROUND
c
c
      subroutine round(fnum,inum)
c     Simple rounding
      IMPLICIT NONE
c
      double precision fnum, rem
      integer inum
c
      inum = int(fnum)
      rem = fnum-inum
      if (rem .gt. 0.5) inum = inum + 1
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDCONST
c
c
      subroutine QMDCONST( fconvert, tconvert )
c---------------------------------------------------------------
c Purpose: energy conversion, and temperature conversion factors
c---------------------------------------------------------------
c
c fconst converts atomic mass units into E[Ry]*t[fs]^2/r[bohr]^2
c    Internal units: fs, Ry, bohr, amu; Output=KE in Rydberg
c    Essentially it converts mv^2/2 into KE, where m is in
c    atomic mass units, v is in bohr/fs, and KE is in Rydberg.
c    fconst = ( Ry*(fs^2)/bohr^2 ) / amu
c    Or, how to best compute it: amu*(bohr/fs)^2 = Ry * fconst
c    Equivalently, think of to convert F/m into a, where
c    F is (Ry/bohr), a is bohr/fs^2, and m is amu.
c
c tconst converts Kelvin to Ry: tconst = Ry/Kelvin
c
      IMPLICIT NONE
c
      DOUBLE PRECISION  fconst, tconst
      DOUBLE PRECISION  fconvert, tconvert
      PARAMETER  ( fconst = 0.468786d0, tconst = 6.3358d-6 )
c
      fconvert = fconst
      tconvert = tconst
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MASS2SYM
c
c
      subroutine mass2sym( imass, symb )
c  Convert the rounded mass to a symbol; symbols only defined
c  through Hg (symbols now defined through Lr, masses only through Hg )
c  8Feb13-PAS/2.65: call new (expanded) symbol table; cosmetics
      IMPLICIT NONE
c
      INTEGER  imass,iatno
      CHARACTER*2 symb
c
      call mass2atno( imass,iatno )
C      call atno2symbol( iatno,symb )
      call ELEMZ2SYMB( iatno, symb )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MASS2ATNO
c
c
      subroutine mass2atno( imass, iatno )
c  Convert rounded mass to atomic number; only defined through Hg
c  8Feb13-PAS/2.65: protect from bad (out-of-range) input; cosmetics
c
      IMPLICIT NONE
      INTEGER imass,iatno,  numtable(201)
      DATA    numtable
     $  / 1, 0, 0, 2, 0,      0, 3, 0, 4, 0,     5, 6, 0, 7, 0,
     $    8, 0, 0, 9, 10,     0, 0, 11, 12, 0,   0, 13, 14, 0, 0,
     $   15, 16, 0, 0, 17,    0, 0, 0, 19, 20,   0, 0, 0, 0, 21, 
     $    0, 0, 22, 0, 0,    23, 24, 0, 0, 25,  26, 0, 0, 28, 0,
     $    0, 0, 0, 29, 30,    0, 0, 0, 0, 31,    0, 0, 32, 0, 33,
     $    0, 0, 0, 34, 35,    0, 0, 0, 36, 37,   0, 0, 38, 39, 0,
     $   40, 0, 41, 0, 0,    42, 0, 43, 0, 0,   44, 0, 45, 0, 0,
     $   46, 0, 47, 0, 0,     0, 48, 0, 0, 49,   0, 0, 0, 50, 0,
     $    0, 51, 0, 0, 0,     0, 53, 52, 0, 0,  54, 0, 55, 0, 0,
     $    0, 56, 0, 57, 58,  59, 0, 0, 60, 61,   0, 0, 0, 0, 62,
     $    0, 63, 0, 0, 0,     0, 64, 0, 65, 0,   0, 0, 66, 0, 67,
     $    0, 68, 0, 69, 0,    0, 0, 70, 0, 71,   0, 0, 72, 0, 0,
     $   73, 0, 0, 74, 0,    75, 0, 0, 0, 76,    0, 77, 0, 0, 78,
     $    0, 79, 0, 0, 0,    80                 /
cxxx: check if Z=53,52 are truly in inverse order of mass, for table
c
      if( imass .gt. 0 .and. imass .le. 120 )then
        iatno = numtable(imass)
      else
c       Assign a Z=0 to an invalid/out-of-range input
        iatno = 0
      endif
c     If we can't determine a good atomic number, just return Carbon's
      if( iatno .eq. 0 ) iatno = 6
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATNO2SYMBOL
c
c
cxxx subroutine atno2symbol(iatno,sym) ! is this safe to be retired now?
cxxx This routine now replaced by new general "element" module routine 8Feb13-PAS
      subroutine xxx_at2symb(iatno,sym)
c     Convert the atomic number to a two character symbol
      IMPLICIT NONE
      integer iatno
      character*2 sym,symtable(80)
      data symtable/'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ',
     $     'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl',
     $     'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe',
     $     'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
     $     'Kr', 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
     $     'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ',
     $     'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',
     $     'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
     $     'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'/
      sym = symtable(iatno)
      return
      end
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TMON
c
c
c---------------------------------------------------------------
c Purpose: monitor temperature over atom types in a MD simulation
c
c Routines:
c  TMONINIT( natm, itypa ) - load type-by-atom information into module
c  TMONRESET( ) - clear running totals to zero/initialize sums
c  TMONPREP( ) - clear ekin sums
c  TMONATOM( i,e1,e2,e3,e4 ) - load ekin for specific atom
c  TMONSTEP( ) - do temperature analysis as of this step
c  TMONENGY( ke,pe,xe,tote ) - do energy analysis as of this step
c  TMONTEMP( ) - evaluate temperature over types for given vel()
c
c Revision history:
c  24Sep08-PAS/2.62: add zero=0 initialize into TMONSTEP/ENGY
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TMONINIT
c
c
      subroutine TMONINIT( natm,itypa )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c QMD/monitor data
      PARAMETER  ( ntypmond=6 )
      COMMON /qmdtmon/ sumttyp(ntypmond,4),sumt100(ntypmond,4),
     $ sumengy(4),sume100(4),
     $ ektyp(ntypmond,4),modesty(ntypmond),itypmon(natmd),ntypmon,ntemp
c
c Input arguments:
      DIMENSION  itypa(natm)
c
c >>>> EXECUTABLE CODE
c
      itypmax = 1
      do  iatm=1,natm
        ityp = itypa(iatm)
        if( ityp .gt. itypmax ) itypmax = ityp
        if( ityp .ge. (ntypmond-1) )then
c         With more types than monitor slots, lump excess together
          itypmon(iatm) = ntypmond - 1
        elseif( ityp .ge. 0 )then
c         Put it into its appropriate slot
          itypmon(iatm) = ityp
        endif
      enddo
      if( itypmax .ge. ntypmond ) itypmax = ntypmond - 1
      ntypmon = itypmax
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TMONRESET
c
      subroutine TMONRESET( )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero=0.d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c QMD/monitor data
      PARAMETER  ( ntypmond=6 )
      COMMON /qmdtmon/ sumttyp(ntypmond,4),sumt100(ntypmond,4),
     $ sumengy(4),sume100(4),
     $ ektyp(ntypmond,4),modesty(ntypmond),itypmon(natmd),ntypmon,ntemp
c
c >>>> EXECUTABLE CODE:
c
      ntemp = 0
c
      do  imon=1,4
        do  ityp=1,ntypmond
          sumttyp(ityp,imon) = zero
          sumt100(ityp,imon) = zero
        enddo
      enddo
c
      do  ityp=1,4
        sumengy(ityp) = zero
        sume100(ityp) = zero
      enddo
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TMONPREP
c
      subroutine TMONPREP( )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero=0.d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c QMD/monitor data
      PARAMETER  ( ntypmond=6 )
      COMMON /qmdtmon/ sumttyp(ntypmond,4),sumt100(ntypmond,4),
     $ sumengy(4),sume100(4),
     $ ektyp(ntypmond,4),modesty(ntypmond),itypmon(natmd),ntypmon,ntemp
c
c >>>> EXECUTABLE CODE:
c
      do  imon=1,4
        do  ityp=1,ntypmond
          ektyp(ityp,imon) = zero
        enddo
      enddo
      do  ityp=1,ntypmond
        modesty(ityp) = 0
      enddo
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TMONATOM
c
      subroutine TMONATOM( iatm, ekint_a,ekinn_a,akint_a,akinn_a )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero=0.d0 )
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c QMD/monitor data
      PARAMETER  ( ntypmond=6 )
      COMMON /qmdtmon/ sumttyp(ntypmond,4),sumt100(ntypmond,4),
     $ sumengy(4),sume100(4),
     $ ektyp(ntypmond,4),modesty(ntypmond),itypmon(natmd),ntypmon,ntemp
c
c >>>> EXECUTABLE CODE:
c
      ityp = itypmon(iatm)
      if( ityp .ge. ntypmond ) ityp = ntypmond - 1
c
      ektyp(ityp,1) = ektyp(ityp,1) + ekint_a
      ektyp(ityp,2) = ektyp(ityp,2) + ekinn_a
      ektyp(ityp,3) = ektyp(ityp,3) + akint_a
      ektyp(ityp,4) = ektyp(ityp,4) + akinn_a
      modesty(ityp) = modesty(ityp) + 3
c
      ityp = ntypmond
      ektyp(ityp,1) = ektyp(ityp,1) + ekint_a
      ektyp(ityp,2) = ektyp(ityp,2) + ekinn_a
      ektyp(ityp,3) = ektyp(ityp,3) + akint_a
      ektyp(ityp,4) = ektyp(ityp,4) + akinn_a
      modesty(ityp) = modesty(ityp) + 3
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TMONSTEP
c
      subroutine TMONSTEP( chtstat )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero = 0.d0 )
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c QMD/monitor data
      PARAMETER  ( ntypmond=6 )
      COMMON /qmdtmon/ sumttyp(ntypmond,4),sumt100(ntypmond,4),
     $ sumengy(4),sume100(4),
     $ ektyp(ntypmond,4),modesty(ntypmond),itypmon(natmd),ntypmon,ntemp
c
c Input declarations:
      CHARACTER*3  chtstat
c Local declarations
      DIMENSION  tempty(ntypmond),tavety(ntypmond)
c
      DATA  lvlprnt / 3 /
c
c >>>> EXECUTABLE CODE:
c
      ntemp = ntemp + 1
      rntemp = DBLE( ntemp)
c
      do 1000 imon=1,4
c
        do  ityp=1,ntypmon
          Ekin = ektyp(ityp,imon)
          modes = modesty(ityp)
          call GET_TEMP_E( Ekin,modes, temp )
          tempty(ityp) = temp
          sumttyp(ityp,imon) = sumttyp(ityp,imon) + temp
          sumt100(ityp,imon) = sumt100(ityp,imon) + temp
        enddo
c
        ityp = ntypmond
        Ekin = ektyp(ityp,imon)
        modes = modesty(ityp)
        call GET_TEMP_E( Ekin,modes, temp )
        tempall = temp
        sumttyp(ityp,imon) = sumttyp(ityp,imon) + temp
        sumt100(ityp,imon) = sumt100(ityp,imon) + temp
c
        do  ityp=1,ntypmon
          tavety(ityp) = sumttyp(ityp,imon) / rntemp
        enddo
        ityp = ntypmond
        taveall = sumttyp(ityp,imon) / rntemp
c       
        if( lvlprnt .gt. 2 )then
          write(imd_wr,9100) chtstat,'TEMP',imon,'(typs;all)=',
     $     itimestep,ntemp,(tempty(ityp),ityp=1,ntypmon),tempall
          write(imd_wr,9100) chtstat,'ave temps',imon,'(typs;all)=',
     $     itimestep,ntemp,(tavety(ityp),ityp=1,ntypmon),taveall
        endif
c
 1000 continue
c
 9100 format(1x,a3,a,i1,a,i8,i7,3x,6f10.3)
c
c Also do averages over even nbunch's of steps:
c
      ibunch = 100
      nbunch = ntemp / ibunch
      ntempx = nbunch*ibunch
      if( ntempx .eq. ntemp )then
c       We have an even nbunch, put out the average over ibunch
c
        rntemp = DBLE( ibunch )
        do  imon=1,4
          do  ityp=1,ntypmon
            tavety(ityp) = sumt100(ityp,imon) / rntemp
            sumt100(ityp,imon) = zero
          enddo
          ityp = ntypmond
          taveall = sumt100(ityp,imon) / rntemp
          sumt100(ityp,imon) = zero
c
          if( lvlprnt .gt. 1 )then
            write(imd_wr,9100) chtstat,'last100T',imon,'(typs;all)=',
     $      itimestep,ntemp,(tavety(ityp),ityp=1,ntypmon),taveall
          endif
        enddo
c
      endif
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TMONENGY
c
      subroutine TMONENGY( Ekinmon,Epotmon,Emdmon )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero = 0.d0 )
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c QMD/monitor data
      PARAMETER  ( ntypmond=6 )
      COMMON /qmdtmon/ sumttyp(ntypmond,4),sumt100(ntypmond,4),
     $ sumengy(4),sume100(4),
     $ ektyp(ntypmond,4),modesty(ntypmond),itypmon(natmd),ntypmon,ntemp
c
c Local declarations
      DIMENSION  engyave(4)
c
      DATA  lvlprnt / 3 /
c
c >>>> EXECUTABLE CODE:
c
      rntemp = DBLE( ntemp)
c
        Etotmon = Ekinmon + Epotmon + Emdmon
        sumengy(1) = sumengy(1) + Ekinmon
        sumengy(2) = sumengy(2) + Epotmon
        sumengy(3) = sumengy(3) + Emdmon 
        sumengy(4) = sumengy(4) + Etotmon
        sume100(1) = sume100(1) + Ekinmon
        sume100(2) = sume100(2) + Epotmon
        sume100(3) = sume100(3) + Emdmon 
        sume100(4) = sume100(4) + Etotmon
c
        do  ityp=1,4
          engyave(ityp) = sumengy(ityp) / rntemp
        enddo
c       
        if( lvlprnt .gt. 4 )then
          write(imd_wr,9100) 'MDENGYAVE (ke,pe,xe,tote)=',
     $     itimestep,ntemp,(engyave(ityp),ityp=1,4)
        endif
c
 1000 continue
c
 9100 format(1x,a,i8,i7,3x,4f15.5)
c
c Also do averages over even nbunch's of steps:
c
      ibunch = 100
      nbunch = ntemp / ibunch
      ntempx = nbunch*ibunch
      if( ntempx .eq. ntemp )then
c       We have an even nbunch, put out the average over ibunch
c
        rntemp = DBLE( ibunch )
c
          do  ityp=1,4
            engyave(ityp) = sume100(ityp) / rntemp
            sume100(ityp) = zero
          enddo
c
          if( lvlprnt .gt. 1 )then
            write(imd_wr,9100) 'MDlast100E (ke,pe,xe,tote)=',
     $      itimestep,ntemp,(engyave(ityp),ityp=1,4)
          endif
c
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TMONTEMP
c
c
      subroutine TMONTEMP( itemp )
c---------------------------------------------------------------
c Purpose: look at temperature by atom type for new velocities
c Written: 17Feb08-PAS for 2.61
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0d0, one=1d0, two=2d0 )
c
c QMD Module declarations:
      PARAMETER  ( natmd =  400 )
      DIMENSION  vel(3,natmd),atom_mass(natmd),vprint(3,natmd)
      COMMON /questmd/ tempinit,tempcurr,tempprnt,timestep,
     $     Ekinprnt,Ekincurr,Epot,shv,vel,
     $     atom_mass,md_type,ntsteps,neqtsteps,md_natoms,itimestep,
     $     md_nmodes,imd_wr,irestarted,vprint,temp_eq,temp_md,
     $     param_eq,param_md,parstat,ieq_type
c QMD/monitor data
      PARAMETER  ( ntypmond=6 )
      COMMON /qmdtmon/ sumttyp(ntypmond,4),sumt100(ntypmond,4),
     $ sumengy(4),sume100(4),
     $ ektyp(ntypmond,4),modesty(ntypmond),itypmon(natmd),ntypmon,ntemp
c
c Local declarations
      DIMENSION  tempty(ntypmond)
c
c >>>> EXECUTABLE CODE:
c
      call QMDCONST( fconst, tconst )
c
      do  ityp=1,ntypmond
        ektyp(ityp,1) = zero
        modesty(ityp) = 0
      enddo
c
      do  iatm=1,md_natoms
        xm = atom_mass(iatm)
c
        vtx = vel(1,iatm)
        vty = vel(2,iatm)
        vtz = vel(3,iatm)
        vt2 = vtx*vtx + vty*vty + vtz*vtz
        ekint_atm = 0.5d0*xm*vt2 / fconst
c
        ityp = itypmon(iatm)
        ektyp(ityp,1) = ektyp(ityp,1) + ekint_atm
        modesty(ityp) = modesty(ityp) + 3
        ityp = ntypmond
        ektyp(ityp,1) = ektyp(ityp,1) + ekint_atm
        modesty(ityp) = modesty(ityp) + 3
      enddo
c
      do  ityp=1,ntypmon
        Ekin = ektyp(ityp,1)
        modes = modesty(ityp)
        call GET_TEMP_E( Ekin,modes,temp )
        tempty(ityp) = temp
      enddo
      ityp = ntypmond
      Ekin = ektyp(ityp,1)
      modes = modesty(ityp)
      call GET_TEMP_E( Ekin,modes,temp )
      tempall = temp
c       
      write(imd_wr,9101) 'TEMPS1: (typs;all)=',itemp,
     $ (tempty(ityp),ityp=1,ntypmon),tempall
 9101 format(1x,a,i6,1x,6f10.3)
c
      RETURN
      END
