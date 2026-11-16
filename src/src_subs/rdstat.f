c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RDSTAT
c
c
      subroutine RDSTAT( IWR, istatfl,
     $ dosetup, doiters, doforce, do_geom, do_cell, opt_spin,
     $ iucstart,igstart,itstart, engytotl,excengy,esnsengy,
     $ istart_sp, elecno, spindata,
     $ rprim, natm, ratm, frctot,strtot )
c---------------------------------------------------------------
c Purpose: acquire status of old calculation and restart info
c
c Written: Peter A. Schultz, 29-November-2001, for v2.51
c
c Revision history:
c  21Jan08-PAS/2.61: spin optimization integrated
c  21Jan08-PAS/2.61: bugfix/FORCE(S) keyword
c---------------------------------------------------------------
c
c The assumption is that the restart file is correct, since it
c is also written by this code, so I will do no error checking.
c This is a bare bones restart file.  E.g., No forces, and no
c velocities or inertia factors are currently included from
c damped dynamics/ASD.  Since this file is for internal use
c only, I will have the freedom to modify this in the future.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      LOGICAL  dosetup,doiters,doforce,do_geom,do_cell,opt_spin
      DIMENSION  rprim(3,3), ratm(3,natm)
      DIMENSION  frctot(3,natm), strtot(3,3)
      DIMENSION  elecno(2), spindata(4,2)
c
c local declarations:
      CHARACTER  label*80, readlbl*6
      CHARACTER  statkey*(*), wantlbl*4
      EQUIVALENCE  ( label, readlbl )
c
c >>>> EXECUTABLE CODE:
c
 8000 format(a80)
 9030 format(1x,3f16.10)
c
      REWIND( unit=istatfl )
c In case status file empty and need to initialize in wrstat:
      wantlbl = 'INPU'
c
      read(istatfl,8000, err=1000, end=1000)  label
c
      write(IWR,*) '>>>>> FOUND RESTART FILE - loading data'
      write(IWR,*   )  label
c
      if( readlbl .eq. 'GLINES' )then
c       Formatted data - number of lines of geom/cell info
        read(istatfl,*)  nline
        write(IWR,*   )  nline
        read(istatfl,8000, err=899, end=899) label
        write(IWR,*   )  label
      endif
c
      if( readlbl .eq. 'IGSTEP' )then
c       Recover atomic coordinates
        read(istatfl,*) igstart
        write(IWR,*   ) 'Restarting from geometry#',igstart
        read(istatfl,8000)  label
        write(IWR,*   )  label
        if( readlbl .ne. 'COORDI' ) call STOPXERR( 'RDSTAT: COORD' )
        do  iatm=1,natm
          read(istatfl,*)  (ratm(i,iatm),i=1,3)
          write(IWR,9030)  (ratm(i,iatm),i=1,3)
        enddo
        read(istatfl,8000, err=899, end=899)  label
        write(IWR,*   )  label
      endif
c
      if( readlbl .eq. 'IUCELL' )then
c       Recover cell parameters
        read(istatfl,*) iucstart
        write(IWR,'(a,i4)') 'Restarting from cell#',iucstart
        read(istatfl,8000)  label
        write(IWR,*   )  label
        if( readlbl .ne. 'primit' ) call STOPXERR( 'RDSTAT: primit' )
        do  id=1,3
          read(istatfl,*)  (rprim(i,id),i=1,3)
          write(IWR,9030)  (rprim(i,id),i=1,3)
        enddo
        read(istatfl,8000, err=899, end=899)  label
        write(IWR,*   )  label
      endif
c
      if( readlbl .eq. 'SPININ' )then
c       This is the latest spin polarization:
        read(istatfl,*)  elecno
        write(IWR,9030)  elecno
        read(istatfl,8000, err=899, end=899)  label
        write(IWR,*   )  label
      endif
c
      if( readlbl .eq. 'SETUP ' )then
c       Do not need to redo setup
        dosetup = .false.
        read(istatfl,8000, err=899, end=899)  label
        write(IWR,*   )  label
      endif
c
      if( readlbl .eq. 'SPINEW' )then
c       This is now the current spin polarization in spin opt:
        read(istatfl,*)  elecno
        write(IWR,9030)  elecno
        read(istatfl,8000, err=899, end=899)  label
        write(IWR,*   )  label
      endif
c
cspdat:start - the following should probably be in separate object
      if( readlbl .eq. 'SPDATA' )then
        read(istatfl,*) istart_sp
        write(IWR,'(1x,i4)') istart_sp
        nspdat = MIN( 2, istart_sp )
        do  ispdat=1,nspdat
          read(istatfl,*)  (spindata(i,ispdat),i=1,4)
          write(IWR,'(1x,4f16.8)')  (spindata(i,ispdat),i=1,4)
        enddo
        read(istatfl,8000, err=899, end=899)  label
        write(IWR,*   )  label
      endif
cspdat:end
c
      if( readlbl .eq. 'ITERAT' )then
c       We have this many iterations done:
        read(istatfl,*)  itstart
        write(IWR,'(1x,i4)')  itstart
        read(istatfl,8000, err=899, end=899)  label
        write(IWR,*   )  label
      endif
c
      if( readlbl .eq. 'SPIN-F' )then
c       This is the final converged spin polarization:
        read(istatfl,*)  elecno
        write(IWR,9030)  elecno
        read(istatfl,8000, err=899, end=899)  label
        write(IWR,*   )  label
      endif
c
      if( readlbl .eq. 'SCF-EN' )then
c       We are converged, and need to acquire some energies
        doiters = .false.
        read(istatfl,*)  engytotl,excengy,esnsengy
        write(IWR,9030)  engytotl,excengy,esnsengy
c
        read(istatfl,8000, err=899, end=899)  label
        write(IWR,*   )  label
      endif
c
      if( readlbl .eq. 'FORCES' )then
c       We have forces and stresses
        doforce = .false.
        do  iatm=1,natm
          read(istatfl,*)  (frctot(i,iatm),i=1,3)
        enddo
        read(istatfl,8000)  label
        if( readlbl .ne. 'STRESS' ) call STOPXERR( 'RDSTAT: stress' )
        do  id=1,3
          read(istatfl,*)  (strtot(i,id),i=1,3)
        enddo
c
        read(istatfl,8000, err=899, end=899)  label
        write(IWR,*   )  label
      endif
c
      if( readlbl .ne. 'END ST' )then
c       If it is not an "end status", there is a bug
        call STOPXERR( 'rdstart: unknown extra data in file' )
      endif
c
  899 write(IWR,*) 'Restart file ended'
      write(IWR,*)
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WRSTAT
c
c
      entry      WRSTAT( IWR, istatfl, statkey,
     $ do_geom, do_cell, opt_spin,
     $ iucstep,igstep,iter, engytotl,excengy,esnsengy,
     $ istep_sp, elecno, spindata,
     $ rprim, natm, ratm )
c
c >>>> EXECUTABLE CODE:
c
      wantlbl = statkey(1:4)
 9100 format(a)
c
c jump here from RDSTAT for initialization of status file:
 1000 continue
c
      REWIND( unit=istatfl )
c
      if( wantlbl .eq. 'INPU' )then
        write(istatfl,9100) 'PRESETUP'
        goto 999
      endif
c
      if( wantlbl .eq. 'CELL' .or. wantlbl .eq. 'GEOM' )then
c       Put out geom/cell info
c
c       First, record how many lines of geom data we will dump:
        nline = 0
        if( do_geom ) nline = nline + 3 + natm
        if( do_cell ) nline = nline + 3 + 3
        if( opt_spin ) nline = nline + 2
        write(istatfl,9100) 'GLINES'
        write(istatfl,*   )  nline
c
        if( do_geom )then
          write(istatfl,9100)  'IGSTEP'
          write(istatfl,*   )   igstep
          write(istatfl,9100)  'COORDINATES'
          do  iatm=1,natm
            write(istatfl,*   ) (ratm(i,iatm),i=1,3)
          enddo
        endif
c
        if( do_cell )then
          write(istatfl,9100) 'IUCELL'
          write(istatfl,*   )  iucstep
          write(istatfl,9100) 'primitive lattive vectors'
          do  id=1,3
            write(istatfl,*   ) (rprim(i,id),i=1,3)
          enddo
        endif
c
        if( opt_spin )then
          write(istatfl,9100) 'SPININ'
          write(istatfl,*   )  elecno
        endif
c
        goto 999
c
      else
c       Need to skip through geom/cell info, if we have it.
c
        read(istatfl,8000) label
        if( readlbl .eq. 'END ST' .or.
     $      readlbl .eq. 'PRESET' .or.
     $      readlbl .eq. 'SETUP ' )then
          REWIND( unit=istatfl )
          goto 1100
        elseif( readlbl .eq. 'GLINES' )then
c         Skip through cell/geom info
          read(istatfl,*)  nline
          do  iline=1,nline
            read(istatfl,*)
          enddo
        else
          write(IWR,*) 'label(1:20)=',label(1:20)
          call STOPXERR( 'WRSTAT - bad format' )
        endif
      endif
c
 1100 continue
c
      if( wantlbl .eq. 'SETU' )then
        write(istatfl,9100) 'SETUP COMPLETED'
        goto 999
      else
        read(istatfl,*)
C        read(istatfl,8000, err=899, end=899)  label
      endif
c
      if( wantlbl .eq. 'SPIT' )then
        write(istatfl,9100) 'SPINEW'
        write(istatfl,*   )  elecno
cspdat:start - following should be done in separate object
        write(istatfl,9100) 'SPDATA'
        write(istatfl,*   )  istep_sp
        nspdat = MIN( 2, istep_sp )
        do  ispdat=1,nspdat
          write(istatfl,*   ) (spindata(i,ispdat),i=1,4)
        enddo
cspdat:end
        goto 999
      endif
c
      if( wantlbl .eq. 'ITER' )then
c       Current iteration count, keeping spin-iteration data if exists
        if( opt_spin .and. ( istep_sp .gt. 1 )  )then
c         Skip spin-pol iteration data - it is present
          ndone_sp = istep_sp - 1
          nspdat = MIN( 2, ndone_sp )
          nskip = 2 + 2 + nspdat
          do  iline=1,nskip
            read(istatfl,*)
          enddo
        endif
        write(istatfl,9100) 'ITERATION'
        write(istatfl,*   )  iter
        goto 999
C      else
C        read(istatfl,*)
C        read(istatfl,*)
      endif
c
      if( wantlbl .eq. 'SP-F' )then
c       Record converged spin-polarization:
        write(istatfl,9100) 'SPIN-FINAL'
        write(istatfl,*   )  elecno
        goto 999
      else
c       Skip a couple lines (i.e., final spin) if spin-opt:
        if( opt_spin )then
          read(istatfl,*)
          read(istatfl,*)
        endif
      endif
c
      if( wantlbl .eq. 'SCF-' )then
        write(istatfl,9100) 'SCF-ENERGIES (tot,exc,del-es)'
        write(istatfl,*   )  engytotl,excengy,esnsengy
        goto 999
      else
        read(istatfl,*)
        read(istatfl,*)
      endif
c
      if( wantlbl .eq. 'FORC' )then
        write(istatfl,9100) 'FORCES'
c       Forces got passed in through atom coordinate array:
        do  iatm=1,natm
          write(istatfl,*   ) (ratm(i,iatm),i=1,3)
        enddo
c
        write(istatfl,9100) 'STRESS TENSOR'
        do  id=1,3
          write(istatfl,*   ) (rprim(i,id),i=1,3)
        enddo
        goto 999
      else
        nline = 1 + natm + 1 + 3
        do  iline=1,nline
          read(istatfl,*)
        enddo
      endif
c
      write(IWR,*) 'This restart key unknown:',statkey
      call STOPXERR( 'WRSTAT: unknown status key' )
c
  999 write(istatfl,9100) 'END STATUS FILE'
cpas: machine dependent flush. - I should replace this (utl)
      call FLUSH(istatfl)
c
c    That's all Folks!
c
      RETURN
      END
