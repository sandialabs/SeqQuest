c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SETGR1C
c
c
      subroutine SETGR1C( i1ctrfl,
     $ natm,ntyp,nshld,nald, n1r,n2r,n3r,nptr,nlat1c,mxwf,
     $ itypa, numshl,lshel,nala,ala,cala, ratm,hh,rlat,
     $ lmx1ctr,almnv1c,  mngr3,mxgr3,ngr,
     $ wk, wk2 )
c---------------------------------------------------------------
c Purpose: manage (non-standard) fast grid atom terms for
c          dynamic 1-center scheme.
c
c Written: Peter A. Schultz, 7-May-2001, for v2.47
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  ratm(3,natm),hh(3), rlat(3,*)
      DIMENSION  lmx1ctr(ntyp),almnv1c(ntyp)
      DIMENSION  mngr3(*),mxgr3(*),ngr(*)
c scratch arrays:
      DIMENSION  wk(nptr,10), wk2(*)
c
c >>>> EXECUTABLE CODE:
c
      REWIND( unit=i1ctrfl )
c
c  Compute products of 1-center orbitals on mesh for rho1ctr:
c
      call RH1ORB( i1ctrfl,
     $ natm,ntyp,nshld,nald, n1r,n2r,n3r,nptr,nlat1c,mxwf,
     $ itypa, numshl,lshel,nala,ala,cala, ratm,hh,rlat,
     $ lmx1ctr,almnv1c,  mngr3,mxgr3,ngr,
     $ wk(1,1), wk(1,4), wk(1,9), wk(1,10), wk2 )
c -->  r-3o     dfn-5s   rsq-s    gsum-s    wf-25s
c
c  Compute short range (compact) parts of basis orbs:
c
      call CMPTORB( i1ctrfl,
     $ natm,ntyp,nshld,nald, n1r,n2r,n3r,nptr,nlat1c,mxwf,
     $ itypa, numshl,lshel,nala,ala,cala, ratm,hh,rlat,
     $ lmx1ctr,almnv1c,  mngr3,mxgr3,ngr,
     $ wk(1,1), wk(1,4), wk(1,7), wk(1,8), wk2 )
c -->  r-3io    rv-3s    rsq-s    gsum-s   wf-5s
c
c  Compute "fictitious" gaussians for counter-density
c  subtraction of multipole moments in rho1ctr:
c
      call CNTRRHO( i1ctrfl,
     $ natm,ntyp,nshld,      n1r,n2r,n3r,nptr,nlat1c,mxwf,
     $ itypa, numshl,lshel,               ratm,hh,rlat,
     $ lmx1ctr,almnv1c,  mngr3,mxgr3,ngr,
     $ wk(1,1), wk(1,4), wk2 )
c -->  r-3i     rv-3s    rho-s
c
c    That's all Folks!
c
      RETURN
      END
