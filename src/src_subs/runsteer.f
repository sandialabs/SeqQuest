c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RUNSTEER
c
c
      subroutine RUNSTEER( scfconv, gconv,  cellconv,
     &                     itstop,  igstop, maxucstep )
c---------------------------------------------------------------
c Purpose: Modify convergence parameters for Quest midstream
c Written: Andrew Pineda and Peter Schultz, for 2.66
c Rebision History
c  none
c---------------------------------------------------------------
c
c     This routine looks for alterable parameters in a file called lcao.steer
c     that is read by the master process. Input is via NAMELIST variables, so that
c     the user can just set the parameters he/she wants to change.
c     This should be called at the end of the SCF, geometry and cell loops.
c
      IMPLICIT NONE
c
      CHARACTER filenm*128, label*6
      INTEGER  IWR, iexist, ISTEER
      LOGICAL  lexist
      INTEGER  nprocs, iproc, master, icomm
      INTEGER  len
c
c Input/Output parameters
      REAL*8   scfconv, gconv, cellconv
      INTEGER  igstop, itstop, maxucstep
c
c Local variables:
c     Since we are doing namelist input we make these variable names match
c     the appropriate QUEST keyword as best we can.
      REAL*8   scf_conv, geo_conv, cell_conv
      INTEGER  scf_max_iters, geo_max_iters, cell_max_iters
c ALTERABLE Parameters - These are convergence parameters that can be changed at run-time
c by namelist input to the master process via lcao.steer. In lcao.steer, you must use the
c variable names below.
      NAMELIST /alterable_parameters/ scf_conv, geo_conv, cell_conv,
     &     scf_max_iters, geo_max_iters, cell_max_iters
c
c >>>> EXECUTABLE CODE
c
c     Set up parallel variables
      call MPNODES( nprocs )
      call MPCOMM( icomm )
      call MPNODE0( master )
      call MPNODE( iproc )
      call FLGETIWR( IWR )
c
c Translate between QUEST variables and NAMELIST names
      scf_conv = scfconv
      geo_conv = gconv
      cell_conv = cellconv
      scf_max_iters = itstop
      geo_max_iters = igstop
      cell_max_iters = maxucstep
c
c STEER convergence criteria with a namelist input file.
      if( iproc .eq. master )then
c       Master process checks for 'lcao.steer'
        call FLNAMESP( 'steer', filenm )
c       Get the standard default output unit:
        INQUIRE( FILE=filenm, EXIST=lexist )
        if( lexist )then
          call FLOPENF( ISTEER, filenm, 'OLD', 'FOR' )
          read( ISTEER, '(a)' ) label
          if( label .eq. 'steeri' )then
 9010       format(a)
            write(IWR,*) 'Found steering file=',filenm
            write(IWR,*) 'Changing job parameters ...'
  100       continue
            read(ISTEER,9010)  label
            if( label(1:3) .eq. 'end' )then
              write(IWR,*) 'End of steering file'
              iexist = 1
              goto 200
            elseif( label .eq. 'iterat' )then
              read(ISTEER,*)  scf_max_iters
              write(IWR,*) 'Changed scf iterations to',scf_max_iters
            elseif( label .eq. 'conver' )then
              read(ISTEER,*)  scf_conv
              write(IWR,'(a,f15.8)') 'changing scfconv to ',scf_conv
            elseif( label .eq. 'gsteps' )then
              read(ISTEER,*)  geo_max_iters
              write(IWR,*) 'Changed geom steops to',geo_max_iters
            elseif( label .eq. 'gconv ' )then
              read(ISTEER,*)  geo_conv
              write(IWR,'(a,f15.8)') 'changing gconv to ',geo_conv  
            elseif( label .eq. 'ucstep' )then
              read(ISTEER,*)  cell_max_iters
              write(IWR,*) 'Changed cell steps to',cell_max_iters
            elseif( label .eq. 'uc_con' )then
              read(ISTEER,*)  cell_conv
              write(IWR,'(a,f15.8)') 'changing uc_conv to ',cell_conv
            else
c             Not known --- this is unsalvageable error
              write(IWR,*) '**** ERROR: unknown runsteer label=',label
              iexist = -1
              goto 200
            endif
            goto 100
          else
            REWIND( ISTEER )
            read( ISTEER, NML=alterable_parameters )
            iexist = 1
            goto 200
          endif
          iexist = -1
  200     continue
          call FLDELETE( ISTEER )
        else
          iexist = 0
        endif
      endif
c
      len = 1
      if( nprocs .gt. 1 ) call MPBCASTI( master, len, iexist, icomm )
c
      if( iexist .eq. 1 )then
c       Synchronize alterable parameters across processes.
        if( nprocs .gt. 1 )then
          call MPBCAST8( master, len, scf_conv,       icomm )
          call MPBCAST8( master, len, geo_conv,       icomm )
          call MPBCAST8( master, len, cell_conv,      icomm )
          call MPBCASTI( master, len, scf_max_iters,  icomm )
          call MPBCASTI( master, len, geo_max_iters,  icomm )
          call MPBCASTI( master, len, cell_max_iters, icomm )
        endif
c
c       Write parameters to main output file
        write(IWR,*) 'Steering file with alterrable parameters found,'//
     &              ' here they are:'
        write( IWR, NML=alterable_parameters )
c
c       Assign alterable parameters back to internal variables
        scfconv   = scf_conv
        gconv     = geo_conv
        cellconv  = cell_conv
        itstop    = scf_max_iters
        igstop    = geo_max_iters
        maxucstep = cell_max_iters
c
      elseif( iexist .eq. -1 )then
c       Now that we told everyone that this went bad ...
c        ... choose to die, rather than ignore.
        call STOPXERR( 'ERROR reading steering file' )
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End
