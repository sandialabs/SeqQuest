c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PMAT0SET
c
c
      subroutine PMAT0SET( IWR, lvlout, iovlpfl,iham0fl,ihxc0fl,ihe0fl,
     $ atm0engy, mat,nmat, convs,convsl,convii,
     $ norb,nk, natm,ntyp,nshld,nald, noad, nlat,nlat2c,nlat3c,
     $ itypa, numshl,lshel,nala,ala,cala,
     $ nocc,locc,nalsh,alsh,calsh, occsh,occij,
     $ alamin, norba, occ,znuc, lmxnlp1,almnnl,
     $ nrd, nrad,nrps,radmsh,radwt,
     $ vpsrad, vatrad,vesrad,vxcrad,excrad, wkrad, wspl,
     $ ratm, rlat, coskr,sinkr, s1atom,
     $ wksm, gntcomb, maxwk,mxang,
     $ wk )
c---------------------------------------------------------------
c Purpose: compute iteration-independent, (semi-)analytic matrix
c          elements over reference atoms: standard dense matrices.
c
c Written: Peter A. Schultz, 21-May-2001 for v2.47
c
c Revision history:
c  20Jun07-PAS/2.60: merge tp and serial
c  30Jul05-PAS/2.59-tp: master-only Ham writes
c  15Oct04-PAS/2.59-tp: clean-up task-parallel routines, purge
c    scratch mem (by again reordering calls).
c  15Oct04-PAS/2.59: reorder orbital mat routines
c  16Apr03-DMC: add scratch array to call.  Will be passed on to
c    nlocmat and vlocmat, to allow parallelization.
c   3Oct02-PAS/2.54: repack memory and i/o for big ham matrices
c   2Oct02-PAS/2.54: opt of non-local - purge killed variables
c  21Jun01-PAS/2.48: replace STOPs
c  21May01-PAS/2.47: extracted from main routine
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  convs(*)
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  nocc(ntyp),locc(nshld,ntyp),
     $ nalsh(nshld,ntyp),alsh(nald,nshld,ntyp),calsh(nald,nshld,ntyp)
      DIMENSION  occsh(nshld,ntyp),occij(nshld,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp), occ(nshld,ntyp), znuc(ntyp)
      DIMENSION  lmxnlp1(ntyp),almnnl(ntyp)
      DIMENSION  nrad(ntyp),nrps(ntyp),radmsh(nrd,ntyp),radwt(nrd,ntyp)
      DIMENSION  vpsrad(nrd,4,ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat), coskr(nk,nlat),sinkr(nk,nlat)
c output arrays:
      DIMENSION  s1atom(noad,noad,nk,natm)
      DIMENSION  vatrad(*),vesrad(*),vxcrad(*),excrad(*)
c scratch arrays:
      DIMENSION  wkrad(nrd,*), wspl(nrd,2), wksm(nk,*)
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
c big scratch array:
      DIMENSION  wk(*)
c local declarations:
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
c find out who is master:
      call MPNODE0( master )
      call MPNODE( iproc )
c
c partition memory in big workspace:
      i01 = 1
      i02 = i01 + nmat
      i03 = i02 + nmat
      i04 = i03 + nmat
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      nlat2c = nlat
      nlat3c = nlat
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                   Overlap matrix
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Set up cutoffs for overlap calculation:
      call MKCUTS( convs(1), convsl, convii )
c
      call SIJ( norb,nk, natm,ntyp,nshld,nald, nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba,
     $ ratm, rlat,coskr,sinkr,
     $ wk(i01) )
c -->  s-o
c
      call MKS1ATM( s1atom,
     $ norb,nk, natm,ntyp, noad, itypa, norba,
     $ wk(i01) )
c -->  s-io
c
      if( nlat2c.lt.nlat .and. lvlout.gt.3 )
     $ write(IWR,*) 'for 2c terms, nlat reduced from',nlat,' to',nlat2c
c
      if( lvlout.gt.2 ) call TIMER('After sij        ')
c
      if( iproc .eq. master )then
c       Save overlap matrix to disk file iovlpfl:
        REWIND( unit=iovlpfl )
        call H0WRIT( iovlpfl, mat,nk, wk(i01) )
      endif
c
      if( lvlout.gt.2 ) call TIMER('After sij/write  ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                Initialize total Hamiltonian
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Will accumulate total Hamiltonian in wk(i01), clear it here:
      call MKZERO( nmat, wk(i01) )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       Matrix elements of local atomic reference potentials
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  "vlocmat" computes matrix elements of local 1-ctr potentials, ie.,
c  sum of coulomb interaction of nuclei, electrostatic and xc potentials
c  due to spherical electron density distributions of reference atoms.
c  Initialize "atm0engy", the calculation of the difference between
c  the system total energy and the sum of the one-electron energies
c  The 2-ctr terms (N:II) are done in "vlocmii".
c
      call MKZERO( nmat, wk(i02) )
      call MKZERO( nmat, wk(i03) )
      call MKZERO( nmat, wk(i04) )
c
c  Set up cutoffs for local potential calculation:
      call MKCUTS( convs(4), convsl, convii )
c
ctp:
c      call VLOCMAT( IWR,
      call PVLOCMAT( IWR,
     $ norb,nk, natm,ntyp,nshld,nald, nlat3c,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,
     $ ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vatrad, vesrad,vxcrad,excrad, wkrad,
     $ wk(i01), wk(i02), wk(i03), wk(i04) )
c -->  vij-so   vloc1-so vloc2-so vloc3-s
c
      if( lvlout.gt.2 ) call TIMER('After vlocmat    ')
c
      atm0engy = zero
      call PVLOCMII( IWR, atm0engy,
     $ norb,nk, natm,ntyp,nshld,nald, nlat3c,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,occij,
     $ ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vatrad, vesrad,vxcrad,excrad, wkrad,
     $ wk(i01), wk(i02), wk(i03), wk(i04) )
c -->  vij-io   vxc-io   exc-io   vtmp-s(3*norb*norbat)
c
      if( lvlout.gt.2 )then
        write(IWR,*) ' '
        write(IWR,9020) ' from vlocmii, atm0engy=',atm0engy
 9020   format(1x,a,3f20.10)
      endif
c
      if( lvlout.gt.2 ) call TIMER('After vlocmii    ')
c
      if( iproc .eq. master )then
c       Save vxc0 and exc0 matrices to ihxc0fl:
        REWIND( unit=ihxc0fl )
        call H0WRIT2( ihxc0fl, mat,nk, wk(i02),wk(i03) )
      endif
c
      if( lvlout.gt.2 ) call TIMER('After vloc/write ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c           Matrix elements of non-local atomic potentials
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  In "nlocmat" compute matrix elements of non-local pseudopotentials,
c  i.e. angular-momentum-dependent contributions of nuclear potential
c  at short distances from each nucleus.
c
c  Assemble some needed angular integrals
      call MKGAUNT( gntcomb, wksm )
c
c  Set up cutoffs for non-local potential calculation:
      call MKCUTS( convs(3), convsl, convii )
c
c  Check memory sufficiency:
      maxuse = i04 - 1 + mxang
      if( maxuse .gt. maxwk )then
        call STOPXERR( 'mem-mat0/ memory err in mat0set/nlocmat' )
      endif
c
      call MKZERO( nmat, wk(i03) )
c
ctp:
c      call NLOCMAT( IWR,
      call PNLOCMAT( IWR,
     $ norb,nk, natm,ntyp,nshld,nald, nlat3c,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba, znuc,
     $ lmxnlp1,almnnl, ratm, rlat,coskr,sinkr, gntcomb,
     $ nrd,nrps, radmsh,vpsrad, wkrad,
     $ wksm(1,1), wksm(1,2),
c -->  vnlr(k)-s vnli(k)-s
     $ wk(i03), wk(i04),        wk(i02) )
c -->  vnlij-io angsav(mxang)-s vnltp-s/tp
c
c  Accumulate total Hamiltonian in wk01, add nonloc-pseudo matrix:
      if( iproc .eq. master )then
        call DAXPY( nmat, one, wk(i03),1, wk(i01),1 )
      endif
c
      if( lvlout.gt.2 ) call TIMER('After nlocmat    ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                   Kinetic matrix elements
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Set up cutoffs for kinetic calculation:
      call MKCUTS( convs(2), convsl, convii )
c
      call TIJ( norb,nk, natm,ntyp,nshld,nald, nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba,
     $ ratm, rlat,coskr,sinkr,
     $ wk(i02) )
c -->  t-o
c
c  Accumulate Hamiltonian in wk01, add kinetic matrix:
      if( iproc .eq. master )then
        call DAXPY( nmat, one, wk(i02),1, wk(i01),1 )
      endif
c
      if( lvlout.gt.2 ) call TIMER('After tij        ')
c
      if( iproc .eq. master )then
c
c       Save kinetic and non-local matrices to ihe0fl:
        REWIND( unit=ihe0fl )
        call H0WRIT2( ihe0fl, mat,nk, wk(i02),wk(i03) )
c
c       Save total reference Hamiltonian "h0" to iham0fl:
        REWIND( unit=iham0fl )
        call H0WRIT( iham0fl, mat,nk, wk(i01) )
c
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    That's all Folks!
c
      RETURN
      END
