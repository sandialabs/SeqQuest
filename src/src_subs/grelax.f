c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRELAX
c
c
      subroutine GRELAX( IWR, gblend, tstep, engytotl,
     $ igstep,nhistg, nrlxdim, chgmax,
     $ natmg,watmg,nxatm, ratmg,frcg, vatmg,amassg, wk3,wk4,
     $ mxwkb, wkb )
c---------------------------------------------------------------
c Purpose: blend to get new coordinates for next geometry
c
c Written: Peter A. Schultz, 18-May-2001, for v2.47
c
c Revision history:
c  21Jun07-PAS/2.60: BFGS disabled
c   7May04-RPM/2.58: add MD stuff
c  25Apr04-PAS/2.58: purge "igtype" passing
c  11Sep03-PAS/2.57: remove itmaxd-dependence
c  10Jan03-PAS/2.55: added LBFGS, limited-history Broyden
c  19Feb02-PAS/2.54: put step increment into this routine
c   4Feb02-PAS/2.51c: dmd velocity/mass file installed
c   8Dec01-PAS/2.51: file manager/io-units
c  30Oct01-PAS/2.50: modularize relax scheme choice
c   3Oct01-PAS/2.50: damped dynamics relaxation
c  18May01-PAS/2.47: extract from main routine
c---------------------------------------------------------------
c
c Relaxation schemes:
c  ig_type  = choice of geometry scheme from following ...
c  ig_steep = simple steepest descents
c  ig_broy  = modified Broyden (D.D. Johnson)
c  ig_md0   = simple dynamics
c  ig_dmd   = "damped dynamics" a la Jonsson - component-wise
c  ig_dmdat = "damped dynamics" a la Jonsson - atom-by-atom
c  ig_asd   = accelerated steepest descent (modified damped dyn.)
c  ig_lbfgs = limited memory BFGS
c  ig_cg    = conjugate gradient (not yet coded)
c  ig_qmd   = molecular dynamics (RPM)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /GEOMTYPS/ ig_type, ig_steep,ig_broy,ig_md0,
     $ ig_dmd,ig_dmdat,ig_asd,ig_lbfgs,ig_cg, ig_qmd
c
c Input/Output array (atomic positions, velocities):
      DIMENSION  ratmg(*)
      DIMENSION  vatmg(*), watmg(*)
c Input:
      DIMENSION  frcg(*), amassg(*)
      DIMENSION  nxatm(*)
c Scratch arrays:
      DIMENSION  wk3(*),wk4(*), wkb(*)
c
c Local declarations:
      DATA  one / 1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      igtype = ig_type
      if( igtype .eq. ig_steep )then
c
        write(IWR,*) '>>>>> Simple blend for new geometry'
c
        do  irlx=1,nrlxdim
          ratmg(irlx) = ratmg(irlx) + gblend*frcg(irlx)
        enddo
c
      elseif( igtype .eq. ig_broy )then
c
c       Modified broyden blend to get new vector.
c
        call BROYDDJ( IWR, 'geom', igstep,nhistg, gblend,
     $   chgmax,nrlxdim,
     $   frcg, ratmg,   wk3,   wk4 )
c -->    vv1-i vv2-io   vv3-s  vv4-s
c
      elseif( igtype .eq. ig_dmd .or. igtype .eq. ig_dmdat .or.
     $        igtype .eq. ig_asd .or. igtype .eq. ig_md0   )then
c
c       "Damped dynamics" or PAS accelerated steepest descent relaxation
c
        if( igtype .eq. ig_dmd )then
c         Jonsson damped dynamics - component-wise
          igdmdmod = 1
        elseif( igtype .eq. ig_dmdat )then
c         Jonsson damped dynamics - atom-by-atom
          igdmdmod = 2
        elseif( igtype .eq. ig_asd )then
c         PAS accel. steep. desc. (=modified damped dynamics)
          igdmdmod = 3
        elseif( igtype .eq. ig_md0 )then
c         Simple dynamics
          igdmdmod = 0
        endif
c
        call DMD( IWR, 'geom', igstep, igdmdmod,
     $   nrlxdim, natmg,nxatm, tstep,
     $   amassg, ratmg,vatmg,frcg, watmg )
c
      elseif( igtype .eq. ig_lbfgs )then
c
c       Limited history BFGS method
c
C        call BFGS_STEP( IWR, 'geom', igstep,nhistg,gblend, engytotl,
C     $   nrlxdim, ratmg,  frcg,  wk3  ,   wkb    ,mxwkb )
Cc -->             rvec-io gvec-i rvecx-s  wk_bfgs(lots)-s
         call STOPXERR( 'BFGS disabled/removed' )
c
      elseif( igtype .eq. ig_qmd )then
c
c        Quest molecular dynamics update
c
         call QMDUPDATE( ratmg,frcg )
c
      else
c
        write(IWR,*) '***** ERROR: relaxation scheme unknown:', igtype
        call STOPXERR( 'GRELAX: unknown relaxation type' )
c
      endif
      igstep = igstep + 1
c
c    That's all Folks!
c
      RETURN
      END
