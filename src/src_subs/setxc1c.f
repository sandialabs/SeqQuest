c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SETXC1C
c
c
      subroutine SETXC1C( ilocrho,
     $ natm,ntyp,nshld,nald, noad,nrd,nfitd, nlat1c, nang,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ occsh,nfityp, nrad,radmsh, angpts,angylm, radwf,  ratm,rlat,
     $ wk, wksm )
c---------------------------------------------------------------
c Purpose: setup portion of xc 1-center dynamic corrections.
c          This is development code only.
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
      DIMENSION  nocc(ntyp),locc(nshld,ntyp),
     $ nalsh(nshld,ntyp),alsh(nald,nshld,ntyp),calsh(nald,nshld,ntyp)
      DIMENSION  occsh(nshld,ntyp)
      DIMENSION  nrad(ntyp),radmsh(nrd,ntyp)
      DIMENSION  ratm(3,natm),rlat(3,*)
      DIMENSION  nfityp(ntyp), angpts(3,nang)
c output arrays:
      DIMENSION  radwf(nang,nrd,noad,nfitd)
      DIMENSION  angylm(nang,9)
c scratch arrays:
      DIMENSION  wk(nrd*nang,9), wksm(nrd,2)
c
c >>>> EXECUTABLE CODE:
c
c  Compute and store density at atom from tails of other atoms
c  Used to get decent density for xc-fits in v1ctrij/blkvofn
c
      call ATMRHO( ilocrho,
     $ natm,ntyp, nshld,nald, noad, nrd, nlat1c, nang,
     $ itypa, nocc,locc,nalsh,alsh,calsh, occsh,nfityp, nrad,radmsh,
     $ ratm,rlat, angpts,
     $ wk(1,1), wk(1,2),wk(1,5),wk(1,8),wk(1,9) )
c -->  rhorad-s rsph-3s rpt-3s  rsq-s   rsqtol-s
c
c  Generate orbitals on radial mesh, for purpose of later xc fitting
c
      call XCORB( nang,nrd,noad,nfitd,
     $ natm,ntyp,nshld,nald, itypa, numshl,lshel,nala,ala,cala,
     $ nrad,radmsh, angpts, nfityp, radwf,
     $ wksm(1,1), wksm(1,2) )
c
c  Generate spherical harmonics on angular grid
c
      call STORYLM( nang, angpts,angylm )
c
c    That's all Folks!
c
      RETURN
      END
