c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DAT0GEOM
c
c
      subroutine DAT0GEOM( do_neb,
     $ natm,nafrc,iforce,ifdefct,
     $ nhistg,gblend,gconv,tstep,igstart,igstop,
     $ nslics,     iatframe )
c---------------------------------------------------------------
c Purpose: initialize parameters for geometry phase
c
c Written: Peter A. Schultz, 23-May-2001, for v2.47
c
c Revision history:
c   1Aug17-PAS/2.68: unified SLICS
c   6Feb13-PAS/2.65: new default settings (gconv)
c   9Sep08-PAS/2.62: option to freeze atoms along vector
c    Oct06-RMV/2.60: initial SLIC options
c  25Apr04-PAS/2.58: purge "igtype"
c  11Sep03-PAS/2.57: remove itmaxd-dependence
c  10Jan03-PAS/2.55: limited-history option
c  28Jan02-PAS/2.54: fixed frame data
c   2Oct01-PAS/2.50: damped dynamics params
c---------------------------------------------------------------
c
c Version 2.65 new default setttings
c   gconv:   old = 0.0020 to new = 0.0008 Ry.bohr ~ 0.02 eV/A
c    - using 4-4 amd even 1-4 for many apps, need at least this
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      LOGICAL    do_neb
      DIMENSION  iforce(natm)
cRMV array for projection flag
      DIMENSION  iatframe(3)
c
c >>>> EXECUTABLE CODE:
c
c For geometry phase ...
c  number of atoms to relax (default to full relax all)
      nafrc = natm
      nvgfix = 0
      do  iatm=1,natm
        iforce(iatm) = 1
      enddo
c  nslics = number of constraints
      nslics = 0
c  iatframe = defines frame, (1)=0 -> no fixed frame
      iatframe(1) = 0
      call SLICINIT( 'INIT2 ' )
c  ifdefct = atom to which to apply force defect (imbalance)
      ifdefct = 0
c  gblend = force blend factor (suggested: [0.5:(1.5):6.0])
      gblend = 1.5d0
c  gconv = max force convergence criterion (Ry/bohr)
      gconv = 0.0008d0
c  tstep = time step for dynamics (suggested: [0.1:(1.0):2.0])
      tstep = 1.0d0
c  set geom blend relaxation type(broyden/asd/damped/dmd-at):
      if( do_neb )then
        call GTYPSET( 'ASD     ' )
      else
        call GTYPSET( 'BROYDEN ' )
      endif
c
c  igstart = initial geometry step for this run
      igstart = 1
c  igstop = maximum number of geometry relaxations steps
      igstop = 50
c  nhistg = number of steps in history to use for geom update
      nhistg = 10
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
c End:
