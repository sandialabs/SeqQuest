c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DAT0CELL
c
c
      subroutine DAT0CELL(
     $ do_cell, ndim, iuctype, icellscheme,strsbroy,
     $ iucstep,iuciter,maxucstep,nhistuc, cellconv,strxtrnl,
     $ strnmax, ucstepfac,cellstepdecr,cellstepincr, ucblend,
     $ strlast )
c---------------------------------------------------------------
c Purpose: initialize user-adjustable parameters for cell opt
c
c Written: Peter A. Schultz, 1-May-2002, v2p53 (cell opt)
c
c Revision history:
c   6Apr13-PAS/2.65: independdent tensor opts
c  23Sep08-PAS/2.62: isotropic cell opt
c  11Sep03-PAS/2.57: remove itmaxd-dependence
c  10Jan03-PAS/2.55: limited-history option
c   1May02-PAS/2.53: adapted Rick Muller's original cell code
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      LOGICAL    do_cell
      DIMENSION  strxtrnl(3,3)
      DIMENSION  strlast(3*3), ucstepfac(3*3)
c
c >>>> EXECUTABLE CODE:
c
c Initialization:
c  iucstep = Cell optimization step number
c  iuciter = Cell step in iterative method (broyden)
      iucstep = 1
      iuciter = 1
c
c Standard user-settable parameters:
c
c  maxucstep = Maximum number of cell steps
      maxucstep = 50
c  nhistuc = Length of history used in blend optimizations
      nhistuc = 10
c  cellconv = Max cell stress convergence criteria (GPa)
      cellconv = 0.0020d0
c  strxtrnl = External stress applied to crystal (GPa)
      call MKZERO( 9, strxtrnl )
c
c Slightly non-standard user-settable parameters:
c
c iuctype=type of constraints on cell opt; 0=none,1=isotropic,>1 restrictions
      iuctype = 0
c  icellscheme = Cell optimization method
c   1=steepest descent; 2=broyden; 0=steepest then broyden
      icellscheme = 0
c  strsbroy = Stress at which to start doing broyden (GPa)
      strsbroy = 5.0d0
c
c Non-standard user-settable parameters:
c
c  strnmax = Size limit on strain in initial cell step (unitless)
      strnmax = 0.05d0
c  ucstepfac = Initial cell stress (GPa) to strain conversion
c       (NB: negative to reflect sign convention for stress)
c       pas: best(?) ~( 1 / (3*bulk_modulus) ), use B=133 GPa
      do  ij=1,3*3
        ucstepfac(ij) = -0.0025d0
      enddo
c  cellstepdecr/incr = Ratio to increase/decrease the step size
      cellstepdecr = 0.5d0
      cellstepincr = 1.5d0
c  ucblend = Initial broyden blend (best as 1.0)
      ucblend = 1.0d0
c  strlast = Previous cell max stress in cell sequence
      do  ij=1,3*3
        strlast(ij) = 9.d9
      enddo
c
      RETURN
c
c    That's all Folks!
c
      END
