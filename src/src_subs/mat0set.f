c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MAT0SET
c
c
      subroutine MAT0SET( IWR, lvlout,
     $ iovlpfl,iham0fl,ivxc0fl,iexc0fl,ih0tfl,ih0nlfl,
     $ atm0engy, mat, nmat, convs,convsl,convii,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, noad, nlat,nlat2c,nlat3c,
     $ itypa, numshl,lshel,nala,ala,cala,
     $ nocc,locc,nalsh,alsh,calsh, occsh,occij,
     $ alamin, norba, occ,znuc, lmxnlp1,almnnl,
     $ nrd, nrad,nrps,radmsh,radwt,
     $ vpsrad, vatrad,vesrad,vxcrad,excrad, wkrad, wspl,
     $ ratm, rlat, coskr,sinkr, s1atom,
     $ wksm, gntcomb, maxwk,mxang,
     $ wkmat )
c---------------------------------------------------------------
c Purpose: compute iteration-independent, (semi-)analytic matrix
c          elements over reference atoms: standard dense matrices.
c
c Written: Peter A. Schultz, 21-May-2001 for v2.47
c
c Revision history:
c  24Nov09-ACP/2.63: reduced footprint stripe-distributed parallelism
c  10Oct08-PAS/2.62: stripe-distributed parallelism
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
c  Memory footprint
c   USING STANDARD H0STORE
c     2*nmatpr + nmatstr + max(norb*norb, 3*norb*norbat)
c   USING H0STOREKP
c     2*nmatpr + max(nmatpr,mxang) + matstr*ncplx+max(norb*norb, 3*norb*maxnorba)
c   where
c     mat     = norb*norb
c     matstr  = slightly bigger than mat (atom-block striped storage)
c     nmat    = nk*norb*norb      
c     nmatstr = nk*matstr*ncplx
c     nmatpr  = O(nmatstr / nprocstr)
c     maxnorba = max( norba(1:natm) ), i.e. most bf/atom
c---------------------------------------------------------------
c
c The distributed parallelism is done via atom-stripes along the
c upper triangle of the diagonal blocks, and done in such a way as
c to distribute the memory requirements roughly equally.
c Each processor gets a single contiguous sequence of atoms as its own,
c with early processors getting fewer atoms with longer lengths and
c later atoms getting more atoms with shorter lengths.
c
c Using H0STORE, the distributed striped matrices are all-k merged
c onto the master proc into full striped matrices, using nk*matstr space.
c The master reorders these into a square i-j order matrix (size mat),
c k-by-k, one matrix(k) at a time, and writes it out for later use.
c  
c Using H0STOREKP, the distributed striped matrix pieces are merged
c onto the master proc one full striped mat[k] (size matstr) at a time.
c The master reorders this into a square i-j order matrix (size mat),
c and writes each out for later use.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input arrays:
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
c Output arrays:
      DIMENSION  s1atom(noad,noad,nk,natm)
      DIMENSION  vatrad(*),vesrad(*),vxcrad(*),excrad(*)
c Scratch arrays:
      DIMENSION  wkrad(nrd,*), wspl(nrd,2), wksm(nk,*)
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
c Big scratch array:
      DIMENSION  wkmat(*)
c
c Local declarations:
      DATA  zero,one / 0.d0,1.d0 /
      LOGICAL    USE_KPSTORE
c     USE_KPSTORE = .true. use reduced footprint scheme
c     USE_KPSTORE = .false. use original scheme
c
c >>>> EXECUTABLE CODE:
c
      USE_KPSTORE=.true.
c     For nk=1, old scheme is smaller by nmatpr:
      if( nk .eq. 1 ) USE_KPSTORE = .false.
      write(IWR,*) 'mat0set: stripe-parallel, KPSTORE=', USE_KPSTORE
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                 Determine work/memory allocations
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc ) 
      call MPNODE0( master )
      call MPCOMM( icomm )
c
c Determine size of full striped matrix
      call MATTOTSTR( matstr,
     $ norb, nk, ncplx, natm, ntyp, itypa, norba )
c
      nmatstr = ncplx*nk*matstr
c Need maximum value of norba for memory allocation.
      maxnorba = norba(1)
      do  ityp=2,ntyp
        maxnorba = MAX( maxnorba, norba(ityp) )
      enddo
c           
c Figure out how much belongs on local processor:
c
      call MATLOCSTR( iproc, nprocs, nprocstr,
     $ myatm1, myatm2, myj0, maxmatpr, mymatpr,
     $ matstr,
     $ norb, nk, ncplx, natm, ntyp, itypa, norba )
c
      nmatpr = mymatpr*ncplx*nk
c
c   NB: nmatpr = O( nmatstr / nprocstr )
c
c
c Partition memory in big workspace:
c
      if( .not. USE_KPSTORE )then
c       Simple all-k merge onto master, needs:
c        i01(nmatpr),i02(nmatpr),
c        i03(nmatstr) [full all-k striped matrices]
c        i04(mat)
c
        i01 = 1
        i02 = i01 + nmatpr
        i03 = i02 + nmatpr
        i04 = i03 + nmatstr
        i0end = i04 + mat
      else
c       Reduced memory, k-by-k merge onto master, needs:
c        i01(nmatpr),i02(nmatpr)
c        i03(nmatpr) [distributed all-k striped matrices]
c         (NB: angsav can safely overlap thru i04 to end of full arena)
c        i04(matstr*ncplx,3*norb*maxnorba),
c        i05(mat)
c
C        nlocmat_angsav_dim = 5*5*5*5*3*81
        nlocmat_angsav_dim = mxang
        i01 = 1
        i02 = i01 + nmatpr
        i03 = i02 + nmatpr
        i04 = i03 + nmatpr
        i05 = i04 + MAX( matstr*ncplx, 3*norb*maxnorba )
        i0end = i05 + mat
      endif
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      nlat2c = nlat
      nlat3c = nlat
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                   Overlap matrix
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Set up cutoffs for overlap calculation:
c
      call MKCUTS( convs(1), convsl, convii )
c
      call SIJ( norb,nk,ncplx, natm,ntyp,nshld,nald, nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba,
     $ ratm, rlat,coskr,sinkr, wksm,
     $ wkmat(i01) )
c -->  smat-o
c -->  nmatpr 
      if( lvlout.gt.2 ) call TIMER('After sij        ')
c
      if( nlat2c.lt.nlat .and. lvlout.gt.3 )
     $ write(IWR,*) 'for 2c terms, nlat reduced from', nlat,
     $                                          ' to', nlat2c
c
c  Extract 1-center overlaps from distributed  matrix for grid-analytic corrections
c
      call MKS1STR( s1atom,
     $ norb,nk,ncplx, natm,ntyp, noad, itypa, norba,
     $ myj0,myatm1,myatm2,
     $ nprocstr, master,iproc,icomm,
     $ wkmat(i01) )
c -->  smat-io
c -->  nmatpr 
c
c  Store overlap matrix:
c
      if( .not. USE_KPSTORE )then
        call H0STORE( iovlpfl , nprocstr, master,iproc,icomm,
     $   nmatstr,nmatpr, mat,
     $   norb,nk,ncplx, natm,ntyp,itypa, norba,
     $   wkmat(i01),        wkmat(i04) )
c -->    hmatstr-is         hmatijk-s
c -->    nmatpr-i/nmatstr-s mat
      else
        call H0STOREKP( iovlpfl , nprocstr, master, iproc, icomm,
     $   nmatstr, nmatpr, mymatpr, mat,
     $   norb, nk, ncplx, natm, ntyp, itypa, norba,
     $   wkmat(i01), wkmat(i05), matstr, wkmat(i04) )
c -->    hmatstr-i   hmatijk-s           hmatkbuf-s
c -->    nmatpr      mat                 matstr*ncplx 
      endif
c
      if( lvlout.gt.2 ) call TIMER('After sij/store  ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                Initialize total Hamiltonian
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Will accumulate total Hamiltonian in wkmat(i01), clear it here:
c
      call MKZERO( nmatpr, wkmat(i01) )
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
      call MKZERO( nmatpr, wkmat(i02) )
      call MKZERO( nmatpr, wkmat(i03) )
c
c  Set up cutoffs for local potential calculation:
c
      call MKCUTS( convs(4), convsl, convii )
c
      call VLOCMAT( IWR,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat3c,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,
     $ ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vatrad, vesrad,vxcrad,excrad, wkrad,
     $ wkmat(i01), wkmat(i02), wkmat(i03), wkmat(i04) )
c -->  vloc0-o     vxc0-o      exc0-o      vtmp-s/tp
c -->  nmatpr      nmatpr      nmatpr    0(Lmatreduce=false)/nmat(Lmatreduce=true)
c
      if( lvlout.gt.2 ) call TIMER('After vlocmat    ')
c
      atm0engy = zero
c
      call VLOCMII( IWR, atm0engy,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat3c,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,occij,
     $ ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vatrad, vesrad,vxcrad,excrad, wkrad,
     $ wkmat(i01), wkmat(i02), wkmat(i03), wkmat(i04) )
c -->  vloc0-io    vxc0-io     exc0-io     vtmp-s
c -->  nmatpr      nmatpr      nmatpr     3*norb*maxnorba
c
      if( lvlout.gt.2 )then
        write(IWR,*) ' '
        write(IWR,9020) ' from vlocmii, atm0engy=',atm0engy
 9020   format(1x,a,3f20.10)
      endif
c
      if( lvlout.gt.2 ) call TIMER('After vlocmii    ')
c
c Pack up contributions to Hamiltonians from local nuc potentials
c We are sorting into (i,j,k) from (k,ij) here
c
c Start from rear, i03=Exc0, so that full merge does not wipe other matrices
c
c  Store Exc0 matrix:
c
      if( .not. USE_KPSTORE )then
        call H0STORE( iexc0fl , nprocstr, master,iproc,icomm,
     $   nmatstr,nmatpr, mat,
     $   norb,nk,ncplx, natm,ntyp,itypa, norba,
     $   wkmat(i03),        wkmat(i04) )
c -->    hmatstr-is         hmatijk-s
c -->    nmatpr-i/nmatstr-s mat
      else
        call H0STOREKP( iexc0fl , nprocstr, master, iproc, icomm,
     $   nmatstr, nmatpr, mymatpr, mat,
     $   norb,nk,ncplx, natm,ntyp,itypa, norba,
     $   wkmat(i03), wkmat(i05), matstr, wkmat(i04))
c -->    hmatstr-i   hmatijk-s           hmatkbuf-s
c -->    nmatpr      mat                 matstr*ncplx
      endif
c
      if( lvlout.gt.2 .and. nprocstr.gt.1 )
     $ call TIMER('  vlocmat/store-E')
c
c  Store Vxc0 matrix:
c
      if( .not. USE_KPSTORE )then
        call H0STORE( ivxc0fl , nprocstr, master,iproc,icomm,
     $   nmatstr,nmatpr, mat,
     $   norb,nk,ncplx, natm,ntyp,itypa, norba,
     $   wkmat(i02),        wkmat(i04) )
c -->    hmatstr-is         hmatijk-s
c -->    nmatpr-i/nmatstr-s mat
      else
        call H0STOREKP( ivxc0fl , nprocstr, master,iproc,icomm,
     $   nmatstr,nmatpr, mymatpr, mat,
     $   norb,nk,ncplx, natm,ntyp,itypa, norba,
     $   wkmat(i02), wkmat(i05), matstr, wkmat(i04) )
c -->    hmatstr-i   hmatijk-s           hmatkbuf-s
c -->    nmatpr      mat                 matstr*ncplx
      endif
c
      if( lvlout.gt.2 .and. nprocstr.gt.1 )
     $ call TIMER('  vlocmat/store-V')
c
      if( lvlout.gt.2 ) call TIMER('After vloc stores')
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
      maxuse1 = i0end - 1
      maxuse2 = i03 - 1 + mxang
      maxuse = MAX( maxuse1, maxuse2 )
      if( maxuse .gt. maxwk )then
        call STOPXERR( 'mem-mat0/ memory err in mat0set/nlocmat' )
      endif
c
      call MKZERO( nmatpr, wkmat(i02) )
c
c     NLOCMAT works with the matrix in stripe-distributed fashion. 
c
      call NLOCMAT( IWR,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat3c,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba, znuc,
     $ lmxnlp1,almnnl, ratm, rlat,coskr,sinkr, gntcomb,
     $ nrd,nrps, radmsh,vpsrad, wkrad,
     $ wksm(1,1), wksm(1,2),
c -->  vnlr(k)-s  vnli(k)-s
     $ wkmat(i02), wkmat(i03),     wkmat(i03) )
c -->  vnlij-o     angsav(mxang)-s vnltp-s/tp
c -->  nmatpr      mxang           0(Lmatreduce=false)/nmat(Lmatreduce=true)
c
c  Accumulate total Hamiltonian in i01, add non-local PP:
c
      call DAXPY( nmatpr, one, wkmat(i02),1, wkmat(i01),1 )
c
c  Store non-local pseudopotential matrix:
c
      if( .not. USE_KPSTORE )then
        call H0STORE( ih0nlfl , nprocstr, master,iproc,icomm,
     $   nmatstr,nmatpr, mat,
     $   norb,nk,ncplx, natm,ntyp,itypa, norba,
     $   wkmat(i02),        wkmat(i04) )
c -->    hmatstr-is         hmatijk-s
c -->    nmatpr-i/nmatstr-s mat
      else
        call H0STOREKP( ih0nlfl, nprocstr, master, iproc, icomm,
     $   nmatstr, nmatpr, mymatpr, mat,
     $   norb, nk, ncplx, natm, ntyp, itypa, norba,
     $   wkmat(i02), wkmat(i05), matstr, wkmat(i04) )
c -->    hmatstr-i   hmatijk-s           hmatkbuf-s
c -->    nmatpr      mat                 matstr*ncplx
      endif
c
      if( lvlout.gt.2 ) call TIMER('After nlocmat    ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                   Kinetic matrix elements
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Set up cutoffs for kinetic calculation:
c
      call MKCUTS( convs(2), convsl, convii )
c
      call TIJ( norb,nk,ncplx, natm,ntyp,nshld,nald, nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba,
     $ ratm, rlat,coskr,sinkr, wksm,
     $ wkmat(i02) )
c -->  tmat-o
c -->  nmatpr
      if( lvlout.gt.2 ) call TIMER('  tij compute    ')
c
c  Accumulate total Hamiltonian in i01, add kinetic matrix:
c
      call DAXPY( nmatpr, one, wkmat(i02),1, wkmat(i01),1 )
c
c  Store kinetic energy matrix:
c
      if( .not. USE_KPSTORE )then
        call H0STORE( ih0tfl , nprocstr, master,iproc,icomm,
     $   nmatstr,nmatpr, mat,
     $   norb,nk,ncplx, natm,ntyp,itypa, norba,
     $   wkmat(i02),        wkmat(i04) )
c -->    hmatstr-is         hmatijk-s
c -->    nmatpr-i/nmatstr-s mat
      else
        call H0STOREKP( ih0tfl , nprocstr, master,iproc,icomm,
     $   nmatstr,nmatpr, mymatpr, mat,
     $   norb,nk,ncplx, natm,ntyp,itypa, norba,
     $   wkmat(i02), wkmat(i05), matstr, wkmat(i04) )
c -->    hmatstr-i   hmatijk-s           hmatkbuf-s
c -->    nmatpr      mat                 matstr*ncplx
      endif
c
      if( lvlout.gt.2 ) call TIMER('After tij        ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c              Save total reference Hamiltonian
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( .not. USE_KPSTORE )then
        call H0STORE( iham0fl, nprocstr, master,iproc,icomm,
     $   nmatstr,nmatpr, mat,
     $   norb,nk,ncplx, natm,ntyp,itypa, norba,
     $   wkmat(i01),        wkmat(i04) )
c -->    hmatstr-is         hmatijk-s
c -->    nmatpr-i/nmatstr-s mat
      else
        call H0STOREKP( iham0fl, nprocstr, master, iproc, icomm,
     $   nmatstr,nmatpr, mymatpr, mat,
     $   norb,nk,ncplx, natm,ntyp,itypa, norba,
     $   wkmat(i01), wkmat(i05), matstr, wkmat(i04) )
c -->    hmatstr-i   hmatijk-s           hmatkbuf-s
c -->    nmatpr      mat                 matstr*ncplx
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MAT0SET1
c
c
      subroutine MAT0SET1( IWR, lvlout,
     $ iovlpfl,iham0fl,ivxc0fl,iexc0fl,ih0tfl,ih0nlfl,
     $ atm0engy, mat,nmat, convs,convsl,convii,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, noad, nlat,nlat2c,nlat3c,
     $ itypa, numshl,lshel,nala,ala,cala,
     $ nocc,locc,nalsh,alsh,calsh, occsh,occij,
     $ alamin, norba, occ,znuc, lmxnlp1,almnnl,
     $ nrd, nrad,nrps,radmsh,radwt,
     $ vpsrad, vatrad,vesrad,vxcrad,excrad, wkrad, wspl,
     $ ratm, rlat, coskr,sinkr, s1atom,
     $ wksm, gntcomb, maxwk,mxang,
     $ wkmat )
c---------------------------------------------------------------
c Purpose: compute iteration-independent, (semi-)analytic matrix
c          elements over reference atoms: standard dense matrices.
c
c Written: Peter A. Schultz, 21-May-2001 for v2.47
c
c Revision history:
c  25Sep08-PAS/2.62: pass ncplx through mat0set
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
c Input arrays:
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
c Output arrays:
      DIMENSION  s1atom(noad,noad,nk,natm)
      DIMENSION  vatrad(*),vesrad(*),vxcrad(*),excrad(*)
c Scratch arrays:
      DIMENSION  wkrad(nrd,*), wspl(nrd,2), wksm(nk,*)
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
c Big scratch array:
      DIMENSION  wkmat(*)
c
c Local declarations:
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      write(IWR,*) 'mat0set1: standard ijkmat'
c Find out who is master:
      call MPNODE0( master )
      call MPNODE( iproc )
c
c Partition memory in big workspace:
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
      call SIJ1( norb,nk,ncplx, natm,ntyp,nshld,nald, nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba,
     $ ratm, rlat,coskr,sinkr,
     $ wkmat(i01) )
c -->  smat-o
c -->  nmat
c
      if( nlat2c.lt.nlat .and. lvlout.gt.3 )
     $ write(IWR,*) 'for 2c terms, nlat reduced from', nlat,
     $                                          ' to', nlat2c
c
      call TIMER('  Before mks1atm')
      call MKS1ATM( s1atom,
c -->               noad*noad*nk*natm
     $ norb,nk,ncplx, natm,ntyp, noad, itypa, norba,
     $ wkmat(i01) )
c -->  smat-io
c -->  nmat
c
      if( lvlout.gt.2 ) call TIMER('After sij        ')
c
      if( iproc .eq. master )then
c       Save overlap matrix to disk file iovlpfl:
        REWIND( unit=iovlpfl )
        call H0WRIT( iovlpfl, mat,nk, wkmat(i01) )
c
        if( lvlout.gt.2 ) call TIMER('After sij/write  ')
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                Initialize total Hamiltonian
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Will accumulate total Hamiltonian in wkmat(i01), clear it here:
      call MKZERO( nmat, wkmat(i01) )
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
      call MKZERO( nmat, wkmat(i02) )
      call MKZERO( nmat, wkmat(i03) )
      call MKZERO( nmat, wkmat(i04) )
c
c  Set up cutoffs for local potential calculation:
c
      call MKCUTS( convs(4), convsl, convii )
c
      call VLOCMAT1( IWR,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat3c,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,
     $ ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vatrad, vesrad,vxcrad,excrad, wkrad,
     $ wkmat(i01), wkmat(i02), wkmat(i03), wkmat(i04) )
c -->  vloc0-o     vxc0-o      exc0-o      vtmp-s/tp(nmat)
c -->  nmat        nmat        nmat
c
      if( lvlout.gt.2 ) call TIMER('After vlocmat    ')
c
      atm0engy = zero
      call VLOCMII1( IWR, atm0engy,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat3c,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,occij,
     $ ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vatrad, vesrad,vxcrad,excrad, wkrad,
     $ wkmat(i01), wkmat(i02), wkmat(i03), wkmat(i04) )
c -->  vloc0-io    vxc0-io     exc0-io     vtmp-s(3*norb*norbat)
c -->  nmat        nmat        nmat
c
      if( lvlout.gt.2 )then
        write(IWR,*) ' '
        write(IWR,9020) ' from vlocmii, atm0engy=',atm0engy
 9020   format(1x,a,3f20.10)
      endif
c
      if( lvlout.gt.2 ) call TIMER('After vlocmii    ')
c
c Pack up contributions to Hamiltonians from local nuc potentials
c We are sorting into (i,j,k) from (k,i,j) here
c
      if( iproc .eq. master )then
c       Save Exc0 matrix to iexc0fl:
        call VLWRIT( iexc0fl, norb,nk, wkmat(i03),  wkmat(i04) )
c -->                                  vkij(nmat)-i vij(mat)-s
c       Save Vxc0 matrix to ivxc0fl:
        call VLWRIT( ivxc0fl, norb,nk, wkmat(i02),  wkmat(i04) )
c -->                                  vkij(nmat)-i vij(mat)-s
      endif
c
      call VLSORT( norb,nk, wkmat(i01), wkmat(i02) )
c -->                       vkijfrom-i  vijkto-o
      call DCOPY( nmat, wkmat(i02),1, wkmat(i01),1 )
c
      if( iproc .eq. master )then
        if( lvlout.gt.2 ) call TIMER('After vloc/write ')
      endif
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
      call MKZERO( nmat, wkmat(i02) )
c
      call NLOCMAT1( IWR,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat3c,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba, znuc,
     $ lmxnlp1,almnnl, ratm, rlat,coskr,sinkr, gntcomb,
     $ nrd,nrps, radmsh,vpsrad, wkrad,
     $ wksm(1,1), wksm(1,2),
c -->  vnlr(k)-s  vnli(k)-s
     $ wkmat(i02), wkmat(i03),     wkmat(i03) )
c -->  vnlij-o     angsav(mxang)-s vnltp-s/tp(nmat)
c
c  Accumulate total Hamiltonian in wk01, add non-local PP:
      if( iproc .eq. master )then
        call DAXPY( nmat, one, wkmat(i02),1, wkmat(i01),1 )
      endif
c
      if( iproc .eq. master )then
c       Save non-local pseudopotential matrix to ih0nlfl:
        REWIND( unit=ih0nlfl )
        call H0WRIT( ih0nlfl, mat,nk, wkmat(i02) )
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
      call TIJ1( norb,nk,ncplx, natm,ntyp,nshld,nald, nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba,
     $ ratm, rlat,coskr,sinkr,
     $ wkmat(i02) )
c -->  tmat-o
c
c  Accumulate total Hamiltonian in wk01, add kinetic matrix:
      if( iproc .eq. master )then
        call DAXPY( nmat, one, wkmat(i02),1, wkmat(i01),1 )
      endif
c
      if( iproc .eq. master )then
c       Save kinetic energy matrix to ih0tfl:
        REWIND( unit=ih0tfl )
        call H0WRIT( ih0tfl, mat,nk, wkmat(i02) )
      endif
c
      if( lvlout.gt.2 ) call TIMER('After tij        ')
c
      if( iproc .eq. master )then
c       Save total reference Hamiltonian "h0" to iham0fl:
        REWIND( unit=iham0fl )
        call H0WRIT( iham0fl, mat,nk, wkmat(i01) )
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLSORT
c
c
      subroutine VLSORT( norb,nk, vkijfrom, vijkto )
c---------------------------------------------------------------
c Purpose: copy matrices from k-first to k-last form
c          Not very efficient, but I do not care.
c
c Written: Peter A. Schultz, 15-September-1999, for v2.39
c
c Revision History:
c  27Mar00-PAS/2.43: rename "enjisrt" to "vljisrt"
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Output array:
      DIMENSION  vijkto(norb,norb,nk)
c Input array:
      DIMENSION  vkijfrom(nk,norb,norb)
c
c >>>> EXECUTABLE CODE:
c     
      do  k=1,nk
c
        do  j=1,norb
          do  i=1,norb
            vijkto(i,j,k) = vkijfrom(k,i,j)
          enddo
        enddo
c
      enddo
c
c    That's all Folks!
c     
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLWRIT
c
c
      subroutine VLWRIT( imatfile, norb,nk, vkij, vij )
c---------------------------------------------------------------
c Purpose: copy matrices from k-first to k-last form
c          Not very efficient, but I do not care.
c
c Written: Peter A. Schultz, 15-September-1999, for v2.39
c
c Revision History:
c  27Mar00-PAS/2.43: rename "enjisrt" to "vljisrt"
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input array:
      DIMENSION  vkij(nk,norb,norb)
c Scratch array:
      DIMENSION  vij(norb,norb)
c
c >>>> EXECUTABLE CODE:
c     
      REWIND( unit=imatfile )
c
      ik = 0
  100 continue
      ik = ik + 1
c
      k = ik
      do  j=1,norb
        do  i=1,norb
          vij(i,j) = vkij(k,i,j)
        enddo
      enddo
c
      mat = norb*norb
      call H0WRIT( imatfile, mat, 1, vij )
c
      if( ik .lt. nk ) goto 100
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
