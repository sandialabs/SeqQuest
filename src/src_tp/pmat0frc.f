c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PMAT0FRC
c
c
      subroutine PMAT0FRC( IWR,IWRF, lvlout,
     $ ndim, dolocxc, mat,nmat, convs,convsl,convii,
     $ nafrc,iforce, frctot,strtot,fdefct, frc1,frc2,frc3, str1,str2,
     $ norb,nk, natm,ntyp,nshld,nald, noad, nlat,nlat2c,nlat3c,
     $ itypa, numshl,lshel,nala,ala,cala,
     $ nocc,locc,nalsh,alsh,calsh, occsh,occij,
     $ alamin, norba, occ,znuc, lmxnlp1,almnnl,
     $ nrd, nrad,nrps,radmsh,radwt,
     $ vpsrad, vatrad,vesrad,vxcrad,excrad, wkrad,
     $ ratm, rlat, coskr,sinkr, nsyma,rprim,rmatsym,naofsym,
     $ wksm , naldsq, gntcomb, maxwk,mxang,
     $ dmat, emat )
c---------------------------------------------------------------
c Purpose: compute contributions from force due to
c          (semi-)analytic matrix elements.
c
c Written: Peter A. Schultz, 21-May-2001 for v2.47
c
c Revision history:
c  12May07-PAS/2.60: add scratch arroy, for tp-serial consistency
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c   6May03-DMC:  Added scratch array, to pass to NLOCFFRC and VLOCFRC
c  23Jul01-PAS/2.49: compact dmat
c  21May01-PAS/2.47: pushed down from main program
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      LOGICAL    dolocxc
      DIMENSION  convs(*)
      DIMENSION  iforce(natm), frctot(3,natm),strtot(3,3),fdefct(3)
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  nocc(ntyp),locc(nshld,ntyp),
     $ nalsh(nshld,ntyp),alsh(nald,nshld,ntyp),calsh(nald,nshld,ntyp)
      DIMENSION  occsh(nshld,ntyp),occij(nshld,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp),occ(nshld,ntyp),znuc(ntyp)
      DIMENSION  lmxnlp1(ntyp),almnnl(ntyp)
      DIMENSION  nrad(ntyp),nrps(ntyp),radmsh(nrd,ntyp),radwt(nrd,ntyp)
      DIMENSION  vpsrad(nrd,4,ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat), coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  rprim(3,3), rmatsym(*),naofsym(natm,*)
      DIMENSION  vatrad(*),vesrad(*),vxcrad(*),excrad(*)
c scratch arrays:
      DIMENSION  frc1(3,natm),frc2(3,natm),frc3(3,natm)
      DIMENSION  str1(3,3),str2(3,3)
      DIMENSION  wkrad(nrd,*), wksm(naldsq,46)
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
c big scratch arrays:
      DIMENSION  dmat(norb*norb)
      DIMENSION  emat(*)
c local declarations:
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
 9020 format(1x,a,3f20.10)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   Compute force contribution from two-center part of Hamiltonian
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Set up cutoffs for 2-center (overlap&kinetic) force calculation:
      call MKCUTS( convs(5), convsl, convii )
c
      call FRC2CTR( ndim,
     $ nafrc,iforce, frc1,frc2, str1,str2,
     $ norb,nk, natm,ntyp,nshld,nald, nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba,
     $ ratm, rlat,coskr,sinkr,
     $ dmat,   emat,  wksm )
c -->  dmat-io emat-i frctp-s(3a)
c
      call SYMFRC( IWRF, 'overlap', str1, frc1, frc3, fdefct,
     $ ndim, natm, nsyma, rprim, rmatsym, naofsym )
c
      call DAXPY( 3*natm, one, frc1,1, frctot,1 )
      do  jd=1,ndim
        do  id=1,ndim
          strtot(id,jd) = strtot(id,jd) + str1(id,jd)
        enddo
      enddo
c
      call SYMFRC( IWRF, 'kinetic energy', str2, frc2, frc3, fdefct,
     $ ndim, natm, nsyma, rprim, rmatsym, naofsym )
c
      call DAXPY( 3*natm, one, frc2,1, frctot,1 )
      do  jd=1,ndim
        do  id=1,ndim
          strtot(id,jd) = strtot(id,jd) + str2(id,jd)
        enddo
      enddo
c
      if( lvlout.gt.2 ) call TIMER('After frc2ctr    ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c      Compute forces due to non-local part of pseudopotential
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Set up cutoffs for non-local potential force calculation:
      call MKCUTS( convs(7), convsl, convii )
c
      call PNLOCFRC( nafrc,iforce, ndim, frc1,str1,
     $ norb,nk, natm,ntyp,nshld,nald, nlat,  nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba, znuc,
     $ lmxnlp1,almnnl,  ratm, rlat,coskr,sinkr,
     $ nrd,nrps, radmsh,vpsrad, wkrad,
     $ dmat, gntcomb, emat,    wksm )
c -->  dmat-io        angsav-s frctp-s(3a)
c
      call SYMFRC( IWRF, 'nlocfrc', str1, frc1, frc3, fdefct,
     $ ndim, natm, nsyma, rprim, rmatsym, naofsym )
c
      call DAXPY( 3*natm, one, frc1,1, frctot,1 )
      do  jd=1,ndim
        do  id=1,ndim
          strtot(id,jd) = strtot(id,jd) + str1(id,jd)
        enddo
      enddo
c
      if( lvlout.gt.2 ) call TIMER('After nlocfrc    ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c    Compute force from local (spherical) part of atom potentials
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Set up cutoffs for local potential force calculation:
      call MKCUTS( convs(8), convsl, convii )
c
      call PVLOCFRC( nafrc,iforce, ndim, frc1,str1,
     $ dolocxc, norb,nk, natm,ntyp,nshld,nald, nlat,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,occij,
     $ lmxnlp1, ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vpsrad, vatrad,vesrad, wkrad,
     $ naldsq, wksm(1,1),wksm(1,7),wksm(1,13),
     $  wksm(1,19),wksm(1,20),wksm(1,21), wksm(1,22),  wksm(1,25),
     $ dmat,   emat )
c -->  dmat-io frctp-s(3a)
c
      call SYMFRC( IWRF, 'vlocfrc', str1, frc1, frc3, fdefct,
     $ ndim, natm, nsyma, rprim, rmatsym, naofsym )
c
      call DAXPY( 3*natm, one, frc1,1, frctot,1 )
      do  jd=1,ndim
        do  id=1,ndim
          strtot(id,jd) = strtot(id,jd) + str1(id,jd)
        enddo
      enddo
c
      if( lvlout.gt.2 ) call TIMER('After vlocfrc    ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    That's all Folks!
c
      RETURN
      END
