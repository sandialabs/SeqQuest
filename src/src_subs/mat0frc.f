c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MAT0FRC
c
c
      subroutine MAT0FRC( IWR,IWRF, idmatfl,idmatsfl,iematfl,iematsfl,
     $ lvlout, nspinf,ndim, dolocxc,
     $ mat,nmat, convs,convsl,convii,
     $ nafrc,iforce, frctot,strtot,fdefct, frc1,frc2,frc3, str1,str2,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, noad, nlat,nlat2c,nlat3c,
     $ itypa, numshl,lshel,nala,ala,cala,
     $ nocc,locc,nalsh,alsh,calsh, occsh,occij,
     $ alamin, norba, occ,znuc, lmxnlp1,almnnl,
     $ nrd, nrad,nrps,radmsh,radwt,
     $ vpsrad, vatrad,vesrad,vxcrad,excrad, wkrad,
     $ ratm, rlat, coskr,sinkr, nsyma,rprim,rmatsym,naofsym,
     $ wksm , naldsq, gntcomb, maxwk,mxang,
     $ dmat )
c---------------------------------------------------------------
c Purpose: compute contributions from force due to
c          (semi-)analytic matrix elements.
c
c Written: Peter A. Schultz, 21-May-2001 for v2.47
c
c Revision history:
c    Jan12-ACO/2.63: install reduced-memory k-parallel dmat distribtino
c  15Oct08-PAS/2.62: stripe-parallel
c  20Jun07-PAS/2.60: merge serial and tp
c  12May07-PAS/2.60: add scratch arroy, for tp-serial consistency
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c   6May03-DMC:  Added scratch array, to pass to NLOCFFRC and VLOCFRC
c  23Jul01-PAS/2.49: compact dmat
c  21May01-PAS/2.47: pushed down from main program
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0, one = 1.d0 )
c
c Input/Output:
      DIMENSION  frctot(3,natm),strtot(3,3),fdefct(3)
c Input arrays:
      LOGICAL    dolocxc
      DIMENSION  convs(*)
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
c Scratch arrays:
      DIMENSION  frc1(3,natm),frc2(3,natm),frc3(3,natm)
      DIMENSION  str1(3,3),str2(3,3)
      DIMENSION  wkrad(nrd,*), wksm(naldsq,46)
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
c
c Big scratch array:
      DIMENSION  dmat(*)
c
c  dmat() total memory requirement:
c    Using dm0fetch, all-k full striped matrices
c      Spin/k=1: nmatpr + max(mat,ncplx*matstr=nmatstr) + mat
c      else:     nmatpr + nmatstr + nspin*mat
c    Using dm0fetchkp, k-parallel reduced memory fetches
c      nmatpr + nmatpr + mat + max( ncplx*matstr, (nspin-1)*mat )
c    Later:
c      nmatpr + max( mxangf, 3*natm )
c
c    Where mxangf is not recorded, but it looks like (not checked):
c      5x5x7x3x9x9 + 25x25 + 2 x 4(25x9x3x3)
c      42525 + 625 + 2x8100 = 59350
c      Exact usage doesn't matter ... this has entire end of arena free
c
c Local declaration;
      LOGICAL    USE_KPSTORE
c
c >>>> EXECUTABLE CODE:
c
      USE_KPSTORE = .true.
c     For nk=1, old scheme is smaller by nmatpr:
      if( nk .eq. 1 ) USE_KPSTORE = .false.
      write(IWR,*) 'mat0frc: stripe-parallel, KPSTORE=',USE_KPSTORE
c
 9020 format(1x,a,3f20.10)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c      Set up parallelism for the atom-centered force calculations
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
c Determine size of full striped matrix
      call MATTOTSTR( matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
c
      nmatstr = ncplx*nk*matstr
c
c Figure out how much belongs on local processor:
      call MATLOCSTR( iproc,nprocs,nprocstr,
     $ myatm1, myatm2, myj0, maxmatpr, mymatpr,
     $ matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
c
      nmatpr = mymatpr*ncplx*nk
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c         Get the density matrices and put into usable form
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c 
      n1dm = nmatstr
      n2dm = mat
      if( .not. USE_KPSTORE )then
c       Using original all-k striped matrices retrieve
        if( nspinf .eq. 2 )then
          if( nk .eq. 1 )then
c           Use the target as temporary storage for up/dn spin retrieve
            n1dm = MAX( n1dm, mat )
          else
c           Need to use extra mat buffer to read both up/dn dmat:
            n2dm = mat + mat
          endif
        endif
        n3dm = 0
      else
c       Using reduced memory single-k-by-k striped matrix retrieve
        n1dm = nmatpr
        n3dm = ncplx*matstr
        if( nspinf .eq. 2 ) n3dm = MAX( n3dm, mat )
      endif
c
c  Load up dmat/emat fetching pointers:
      i1dm = 1
      i2dm = i1dm + n1dm
      i3dm = i2dm + n2dm
c
      i1em = i1dm + nmatpr
      i2em = i1em + n1dm
      i3em = i2em + n2dm
c
c     Determine total required storage for dmat/emat:
      memdm = i3em + n3dm - 1
      if( memdm .gt. maxwk )then
c       We do not have enough space for this ...
        call FLGETIERR( IERR )
        write(IERR,*) 'FORCES/MAT0FRC: memory inadequate'
        write(IERR,*) 'Needed, available memory=',memdm,maxwk
        call STOPXERR( 'MAT0FRC memory shortage' )
      endif
c
c  Fetch and distribute density matrix dmat from idmatfl
c    and then
c  Fetch and distribute energy-weighted density matrix emat from iematfl
c
      if( .not. USE_KPSTORE )then
c       Use all-k striped intermediate
        call DM0FETCH( idmatfl,idmatsfl, nspinf,
     $   norb,nk,ncplx, natm,ntyp,itypa,norba,
     $   dmat(i1dm),                 dmat(i2dm) )
c -->    dmatstr(nmatpr-o)           dmatij(mat)[2mat if spin.and.k>1]
c -->    nmatstr|mat(spin&k=1)-s     mat-s (2mat if spin&k>1)
        call DM0FETCH( iematfl,iematsfl, nspinf,
     $   norb,nk,ncplx, natm,ntyp,itypa,norba,
     $   dmat(i1em),                 dmat(i2em) )
      else
c       Do reduced-memory single-k striped intermediates
        call DM0FETCHKP( idmatfl,idmatsfl, nspinf,
     $   norb,nk,ncplx, natm,ntyp,itypa,norba,
     $   dmat(i1dm),        dmat(i2dm),   dmat(i3dm) )
c -->    dmatstr(nmatpr)-o  dmatij(mat)-s dmatbuf(nc*matstr|mat[spin])-s
        call DM0FETCHKP( iematfl,iematsfl, nspinf,
     $   norb,nk,ncplx, natm,ntyp,itypa,norba,
     $   dmat(i1em),        dmat(i2em),   dmat(i3em) )
      endif
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
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba,
     $ ratm, rlat,coskr,sinkr,
     $ dmat(i1dm), dmat(i1em),  wksm )
c -->  dmat-io     emat-i       frctp-s(3a)
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
      call NLOCFRC( nafrc,iforce, ndim, frc1,str1,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat,  nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba, znuc,
     $ lmxnlp1,almnnl,  ratm, rlat,coskr,sinkr,
     $ nrd,nrps, radmsh,vpsrad, wkrad,
     $ dmat(i1dm), gntcomb, dmat(i1em),    wksm )
c -->  dmat-io              angsav-s       frctp-s(3a)
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
      call VLOCFRC( nafrc,iforce, ndim, frc1,str1,
     $ dolocxc, norb,nk,ncplx, natm,ntyp,nshld,nald, nlat,nlat2c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ alamin, norba, znuc,occsh,occij,
     $ lmxnlp1, ratm, rlat,coskr,sinkr,
     $ nrd,nrad, radmsh,radwt, vpsrad, vatrad,vesrad, wkrad,
     $ naldsq, wksm(1,1),wksm(1,7),wksm(1,13),
     $  wksm(1,19),wksm(1,20),wksm(1,21), wksm(1,22),  wksm(1,25),
     $ dmat(i1dm), dmat(i1em) )
c -->  dmat-io     frctp-s(3a)
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
c         
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:        
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
