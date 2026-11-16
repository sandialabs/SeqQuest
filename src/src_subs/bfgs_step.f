c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BFGS_STEP
c
c
      subroutine BFGS_STEP( IWR, sufx, igstep,nhistg,gblend, fcnvalue,
     $ lenvec, rvec,dvec,rvecx, wk_bfgs,mxwk )
c---------------------------------------------------------------
c Purpose: driver for geometry minimization using LBFGS
c
c Written: Peter A. Schultz, 19-Decemeber-2002
c          based on original code by David Raczkowski
c
c Revision History:
c  20Dec02-PAS/2.55: adapted from original Raczkowski code
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Variables:
c   sufx = character tag to identify minimized quantity, name files
c   nhistg = limit on number of steps included in BFGS minimization
c   gblend = blend parameter for initial step (>0)
c   fcnvalue = value of quantity as function of current vector
c   lenvec = length of vector to be optimized
c   rvec() = vector to be optimized
c   dvec() = derivatives of quantity with respect to vector
c
c   rvecx(),wk_bfgs() = work space
c   mxwk = amount of work space available in wk_bfgs()
c
c ***** Important Notes *****
c   (1) This routine does not cleanly terminate a BFGS run.
c   Convergence is done externally based on derivatives, but
c   the BFGS minimizer has an internal switch.  It needs that
c   this internal switch be communicated outwards.  At some
c   point, this needs to be done, but not today. 07Jan03-PAS
c   Further note: iflag_bfgs carries the critical info?
c   On return =0, it is complete without errors
c   On return <0, it indicates an error
c   On return >0, to continue appropriately instructed
c   (2) iprint lies - there is internal output in lbfgs_g that
c   cannot be shut off with this flag.  This should also be fixed.
c   I do not want to modify the internals of the BFGS routine.
c
c Input/Output:
      CHARACTER  sufx*(*)
      DIMENSION  rvec(lenvec), dvec(lenvec)
c Scratch:
      DIMENSION  rvecx(lenvec)
      DIMENSION  wk_bfgs(*)
c Local:
      CHARACTER  gfile*(40),bfgsfl*(40)
      DATA  zero,one / 0.d0,1.d0 /
      LOGICAL  printon
c
c Local/BFGS requirements:
      EXTERNAL LB2_G
      COMMON /LB3_G/ IWR_BFGS,IERR_BFGS,GTOL,STPMIN,STPMAX
      LOGICAL  accept_g, diagco
      DIMENSION  iprint(2)
      DATA  iprint / -1, 0 /
c      (1) Every 'i'th iteration, or 0=first/lst, or <0=none.
c      (2) 0=min, 1-3=more
      DATA  eps_bfgs, xtol_bfgs / 1.d-5, 1.d-16 /
c
c >>>> EXECUTABLE CODE:
c
      printon = .true.
c Length of BFGS history (3-7 is recommended):
      nhist = 5
      if( nhistg .gt. 0 ) nhist = nhistg
      write(IWR,*) '>>>> Performing LBFGS step, history limit=',nhist
      call BFGS_WARN( IWR )
c
      IWR_BFGS = IWR
      IERR_BFGS = IWR
c
c Determine and check space needed in wk_bfgs():
      lenwk = lenvec*(2*nhist+1) + 2*nhist
      if( lenwk .gt. mxwk )then
        write(IWR,*) 'BFGS: wk mem needs/max=',lenwk,mxwk
        call STOPXERR( 'BFGS-mem: wk space insufficient' )
      endif
c
c Reverse the sign on the forces being optimized:
      call DSCAL( lenvec, -one, dvec,1 )
c
c Open binary files for bfgs
      gfile  = 'bfgs1_' // sufx
      bfgsfl = 'bfgs2_' // sufx
      call FLOPENB( igfile , gfile  )
      call FLOPENB( ibfgsfl, bfgsfl )
c
      if( igstep .eq. 1 )then
c       On initial step, do simple blend
        if( printon )
     $  write(IWR,*) ' LBFGS: first step: simple blend'
c
        iflag_bfgs = 0
c
        REWIND( igfile )
        call WRITBIG( igfile, lenvec, rvec )
c
c       Put simple blend param into rvecx:
        diagco = .true.
        do  i=1,lenvec
          rvecx(i) = gblend
        enddo
c
        call LBFGS_G( lenvec,nhist, rvec,fcnvalue,dvec,
     $   diagco,rvecx,
     $   iprint,eps_bfgs,xtol_bfgs, wk_bfgs,
     $   iflag_bfgs, accept_g )
c
        call WRITBIG( igfile, lenvec, rvec )
c
      else
c       After first step, enter BFGS minimization
        accept_g = .false.
c
c       Get old geom into rvecx
        diagco = .false.
        REWIND( igfile )
        call READBIG( igfile, lenvec, rvecx )
c
c       Retrieve bfgs history
        REWIND( ibfgsfl)
        read( ibfgsfl )  iflag_bfgs
        call READBIG( ibfgsfl, lenwk, wk_bfgs )
c
        call LBFGS_G( lenvec, nhist, rvec,fcnvalue,dvec,
     $   diagco,rvecx,
     $   iprint,eps_bfgs,xtol_bfgs, wk_bfgs,
     $   iflag_bfgs, accept_g )
c
        if( accept_g ) then
          if( printon )
     $    write(IWR,*) ' LBFGS: vector accepted, step=',igstep
c
c         Get the recent accepted geom ...
          call READBIG( igfile, lenvec, rvecx )
          REWIND( igfile )
c          ... and put into past slot:
          call WRITBIG( igfile, lenvec, rvecx )
c
        else
          if( printon )
     $    write(IWR,*) ' LBFGS: vector not accepted, step=',igstep
c
c         Get old (good) geom (leaving it in "past" slot) ...
          REWIND( igfile )
          call READBIG( igfile, lenvec, rvecx )
c
        endif
c       Put new geom into "next" slot:
        call WRITBIG( igfile, lenvec, rvec )
c
      endif
c
      write(IWR,*) ' LBFGS: output iflag_bfgs=',iflag_bfgs
c
      REWIND( ibfgsfl)
      write( ibfgsfl )  iflag_bfgs
      call WRITBIG( ibfgsfl, lenwk, wk_bfgs )
c
      call FLCLOSE( igfile )
      call FLCLOSE( ibfgsfl )
c
      if( printon )
     $write(IWR,*) '<<<< LBFGS: step completed'
c
      RETURN
c
c    That's all Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BFGS_WARN
c
c
      subroutine BFGS_WARN( IWR )
c
      INTEGER  IWR
c
      write(IWR,1300)
 1300 format(
     $ 1x,55('*'),
     $/1x,'***** CAUTION: LBFGS optimizer not cleanly integrated',
     $/1x,'***** Completion&error signals from LBFGS_G ignored',
     $/1x,'***** Optimization wth LBFGS may not terminate properly',
     $/1x,55('*') )
c
      RETURN
      END
