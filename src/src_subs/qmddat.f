c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> QMDDAT
c
c
      subroutine QMDDAT( IWR,IDAT )
c--------------------------------------------------------------
c Purpose: set parameters for Quest-MD run
c Written: Richard P. Muller, for 2.58 (MD)
c Revision history:
c   9Apr15-PAS/2.66 lnwrite lines
c  10Feb08-PAS/2.61: wire in more thermostat controls
c  10May04-PAS/2.58: extracted to Quest code, and polished
c--------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( nline = 128 )
      CHARACTER  label*(nline), readlbl*6, wantlbl*6
      EQUIVALENCE  (label,readlbl)
      CHARACTER  mdlabel*8
      INTEGER    istart, iend
c
c >>>> EXECUTABLE CODE:
c
 8000 format(a128)
 9000 format(1x,a)
 9010 format(1x,3f20.8)
c
 100  continue
      read(IDAT,8000,err=1300,end=1300)  label
      call LNWRIT( IWR, label, nline )
c
cqmd: "end dynamics" OR "end md" in qmdinput
      if( readlbl .eq. 'end dy' .or.
     $    readlbl .eq. 'end md' )then
c       End of MD data input section
        goto 900
c
      elseif( readlbl .eq. 'dynami' .or.
     $        readlbl .eq. 'md dat' )then
c       Start of dynamics data input
        write(IWR,*) '>>> Reading MD data'
c
      elseif( readlbl .eq. 'timest' )then
c       Set time step for dynamics (fs)
        read(IDAT,*)  dt_md
        write(IWR,'(1x,f10.5)')  dt_md
        call QMDTSTEPSET( dt_md )
c
      elseif( readlbl .eq. 'md_met' )then
c       Specify MD method/thermostat
        read(IDAT,8000,err=1300,end=1300)  label
        call LNWRIT( IWR, label, nline )
c
        call STRPARS( label, 1,72, istart, iend)
        if( istart .eq. 0 ) call STOPXERR('no md type - blank label')
        mdlabel = '        '
        iend = MIN( iend, istart+7 )
        nlbl = iend - istart + 1
        mdlabel(1:nlbl) = label(istart:iend)
        call QMDTYPSET( 4, mdlabel )
c
      elseif( readlbl .eq. 'md_par' )then
c       Set thermostat param - equilibrarion phase
c        (damp for Tscale/Langvn, time for Berend,Hoover)
        read(IDAT,*)  factormd
        write(IWR,'(1x,f14.5)')  factormd
        call QMDPARAMSET( factormd )
c
      elseif( readlbl .eq. 'md_ste' )then
c       Set number of (thermalized) MD steps
        read(IDAT,*)  nstep
        write(IWR,'(1x,i8)')  nstep
        call QMDNTSTEPSET( nstep )
c
      elseif( readlbl .eq. 'md_tem' )then
c       Set temperature (Kelvin) for MD simulation
        read(IDAT,*)  temp
        write(IWR,'(1x,f14.5)')  temp
        call QMDTEMPSET( temp )
c
      elseif( readlbl .eq. 'eq_met' )then
c       Specify equilibration method/thermostat
        read(IDAT,8000,err=1300,end=1300)  label
        call LNWRIT( IWR, label, nline )
c
        call STRPARS( label, 1,72, istart, iend )
        if( istart .eq. 0 ) call STOPXERR('no md type - blank label')
        mdlabel = '        '
        iend = MIN( iend, istart+7 )
        nlbl = iend - istart + 1
        mdlabel(1:nlbl) = label(istart:iend)
        call QMDTYPSET( 1, mdlabel )
c
      elseif( readlbl .eq. 'eq_par' )then
c       Set thermostat param - equilibrarion phase
c        (damp for Tscale/Langvn, time for Berend,Hoover)
        read(IDAT,*)  factormd
        write(IWR,'(1x,f14.5)')  factormd
        call QMDEQPARSET( factormd )
c
      elseif( readlbl .eq. 'equili' .or. readlbl .eq. 'eq_ste' )then
c       Set number of equilibration steps
        read(IDAT,*)  nstep
        write(IWR,'(1x,i8)')  nstep
        call QMDNEQSTEPSET( nstep )
c
      elseif( readlbl .eq. 'eq_tem' )then
c       Set temperature (Kelvin) for equilibration
        read(IDAT,*)  temp
        write(IWR,'(1x,f14.5)')  temp
        call QMDEQTSET( temp )
c
      else
c       Did not find a matching keyword
        call STOPXERR('qmddat - unknown keyword')
c
      endif
c
      goto 100
c
c  Error Handling
c
 1300 continue
      write(IWR,*)  'QMD data: error/end reading QMD data'
      write(IWR,*) readlbl
      call STOPXERR( 'QMDDAT: format error in input file' )
c
c  Normal return
c
 900  continue
c
      RETURN
c
c    That's all Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DAT0QMD
c
c
      subroutine DAT0QMD
c--------------------------------------------------------------
c Purpose: set defaults for Quest-MD run
c Written: Richard P. Muller, for 2.58 (MD)
c Revision history:
c  10Feb08-PAS/2.61: wire in more thermostat controls
c  10May04-PAS/2.58: extracted to Quest code, and polished
c--------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
c Set default output to Quest output:
      call FLGETIWR( IWR )
      call QMDSETIWR( IWR )
c
c Set default MD method to microcanonical ensemble - NVE
      call QMDTYPSET( 4, 'NVE   ' )
c
c Set default equilibration method to TSCALE
      call QMDTYPSET( 1, 'TSCALE' )
c
c Set default temperature (Kelvin) for MD simulation
      temp = 300.d0
      call QMDTEMPSET( temp )
c
c Set default temperature (Kelvin) for equilibration
      temp = 300.d0
      call QMDTEMPSET( temp )
c
c Set default time step (fs)
      dt_md = 0.5
      call QMDTSTEPSET( dt_md )
c
c Set default number of equilibration steps
      nstep = 0
      call QMDNEQSTEPSET( nstep )
c
c Set default number of (thermalized) MD steps
      nstep = 100 000
      call QMDNTSTEPSET( nstep )
c
      RETURN
c
c    That's all Folks!
c
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End
