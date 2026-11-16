c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VSLOFRC
c
c
      subroutine VSLOFRC( idmatfl, ndim,ncplx,idwf,doblas3,
     $ nafrc,iforce, frcslo,strslo,
     $ norb,nk, natm,ntyp, nshld,nald, nlat,
     $ itypa,numshl,lshel,nala,ala,cala, alamin,
     $ ratm, rlat, wtk,coskr,sinkr,
     $ weight, hh, ndbox,nbox,nrec,inrec,noboxs,inboxs, boxrad,
     $ rmesh, vbox, wkbox, frctp,
     $ wvfcns,wf,
     $ wkbi, dphi,
     $ dmatrd, dmat,  vslo )
c      R8(mat) R4(ncplx*nkloc*mat), ncplx*nkloc*nsumnrbox
c---------------------------------------------------------------
c Purpose: compute force/stress from grid potentials
c
c Written: Peter A. Schultz,  3-June-1997, for v2.20 (boxes)
c
c Revision history:
c   4Dec08-PAS/2.63: k-parallel dmat matrices
c   4Jan08-PAS/2.61: reduced (dist) memory for grid wvfcns
c  11May07-PAS/2.60: merged tp and serial
c   6Jul06-APT/2.59: merged task parallel and image parallel
c  30Jul05-PAS/2.59: master-only matrix files
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c   6May03-DMC: Add frctp to parameter list, to provide a local
c               copy of the frcslo array, for parallelization.
c  23Jul01-PAS/2.49: compact dmat; spin-polarized dft
c   2May01-PAS/2.47: remove disk storage of grid orbitals
c  27Sep99-PAS/2.40: full three-dimensional stress tensor
c  10Jun99-PAS/2.36: pull into blas
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  17Dec98-PAS/2.30: clean force/stress arrays, remove grid es
c                    contribution its own routine ("eslofrc")
c   9Nov98-PAS/2.28: single precision grid orbitals/wf4
c  10Sep98-PAS/2.26: on-the-fly grid orbitals
c  23Jul98-PAS/2.23: adapted for split basis
c---------------------------------------------------------------
c
c Notes:
c  Assumes potential vslo has been shuffled into boxes on input
c
c  Nearby stuff not implemented; to do and keep efficient SUM(dmij*ph-j)
c  construction, do atom-wise sums
c
c  Can try to move twowt factor up into more useful place
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      REAL  dmat
      REAL  wvfcns,wf, wkbi
      REAL  vbox
      REAL  dphi
c
c Output arrays:
      DIMENSION  frcslo(3,natm), strslo(3,3)
c Input arrays:
      LOGICAL    doblas3
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp),
     $ nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp)
      DIMENSION  ratm(3,natm), rlat(3,nlat)
      DIMENSION  wtk(nk),coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  hh(3), ndbox(2,2,3),noboxs(nbox),inboxs(nrec,nbox)
      DIMENSION  wvfcns(2,*),wf(*)
      DIMENSION  vslo(*)
c Scratch arrays:
      DIMENSION  wkbi(*), wkbox(*)
      DIMENSION  dmatrd(norb,norb)
      DIMENSION  dmat(*)
c                dmat(ncplx*nkloc,norb,norb)
      DIMENSION  rmesh(3,*)
      DIMENSION  vbox(*)
      DIMENSION  dphi(*)
c Scratch/parallel:
      DIMENSION  frctp(3,natm), strtp(3,3)
c
c Local declarations:
      DIMENSION  r0box(3),boxctr(3)
      DATA  zero,two / 0.d0,2.d0 /
      DATA  r4zero / 0.e0 /
c
c >>>> EXECUTABLE CODE:
c
c  Clear out force arrays
      call MKZERO( 3*natm, frcslo )
      call MKZERO( 9, strslo )
c
c  Get Local MP info
      call MPNODES( nprocl )
      call MPNODE( iprocl )
      call MPNODE0( masterl )
      call MPCOMM( icomml )
c
      call MPNODES_K( nprocs )
      call MPNODE_K( iproc )
      call MPNODE0_K( master )
      call MPCOMM_K( icomm )
c  Set up k-parallel
      call KPMINE( nk,nkloc,nk0 )
c
c  Get and sort density matrix into k-first form
c
      call KPDMATRD2( idmatfl, norb, ncplx, nk,nk0,nkloc,
     $ dmatrd,      dmat )                   
c      r8(norb^2)-s r4(ncplx*norb^2*nkloc)-o
      ipdm2 = ncplx*nkloc*norb*norb + 1
c
c  Initialize location of grid orbitals
      ir00 = 0
c
      irpot = 0
c
c Cycle over boxes of grid points
c
      do 1000 ibox=1,nbox
c
c       Get basic box info
        call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $   ira0,irb0,irc0, idra,idrb,idrc )
        nrbox = idra*idrb*idrc
c
c       Update location in global mesh potential array
        iv0 = irpot
        irpot = irpot + nrbox
c
c       Get number of orbitals in this box; skip if none
        nobox = noboxs(ibox)
        if( nobox.le.0 ) goto 1000
c
c       Each processor will do every nprocs'th box in this loop.
        if( iproc .ne. MOD( ibox, nprocs ) ) goto 1000
c
c       Generate mesh point coordinates
        call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
        call MKRBOX( hh,r0box, idra,idrb,idrc, rmesh )
c
c       Get grid orbitals:
c
        if( idwf.eq.0 )then
c         Compute orbitals on mesh in this box
c
          ir00 = 0
          iwk1 = 1
          iwk2 = 1
C          iwk2 = iwk1 + 3*nrbox
          iwk3 = iwk2 + nrbox
          iwk4 = iwk3 + nrbox
          iwk5 = iwk4 + nrbox
          iwk6 = iwk5 + nrbox
          iwk7 = iwk6 + nrbox
          iwk8 = iwk7 + nrbox
c
          call BOXWFN( ibox, ncplx,
     $     nk,nk0,nkloc, natm,ntyp, nshld,nald, nlat,
     $     itypa,numshl,lshel,nala,ala,cala,alamin,
     $     ratm,rlat,coskr,sinkr,
     $     nobox,nrbox,nrec,inrec,inboxs, boxrad,boxctr,
     $     wvfcns(1,ir00+1), wf(ir00+1),
     $     rmesh      ,wkbox(iwk2),wkbox(iwk3),wkbox(iwk4),
     $     wkbox(iwk5),wkbox(iwk6),wkbox(iwk7),wkbox(iwk8) )
c
        endif
c
c       Put r*8 potential into r*4 space for use below:
        do  nr=1,nrbox
          vbox(nr) = vslo(iv0+nr)
        enddo
c
c       Carve up space for derivatives of grid wave functions:
        nrkbox = ncplx*nkloc*nrbox
        id1 = 1
        id2 = id1 + 1*nrkbox
        id3 = id1 + 2*nrkbox
        id4 = id1 + 3*nrkbox
c
        if( ncplx.eq.1 )then
c
c         Real orbitals
c
          if( doblas3 )then
            call VSLOF1B( ibox, ndim, frcslo,strslo,
     $       norb, natm,ntyp, nshld,nald, nlat,
     $       itypa,numshl,lshel,nala,ala,cala, alamin,
     $       ratm, rlat, coskr,sinkr,
     $       hh, nbox,nrec,inrec, nrbox,nobox,inboxs, boxrad,boxctr,
     $       wf(ir00+1), rmesh, vbox, wkbox, wkbi,
     $       dphi(id1),dphi(id2),dphi(id3),  dphi(id4),
c -->        dphix-s   dphiy-s   dphiz-s     dphiuv-s
     $       dmat, dmat(ipdm2) )
          else
            call VSLOF1( ibox, ndim, frcslo,strslo,
     $       norb, natm,ntyp, nshld,nald, nlat,
     $       itypa,numshl,lshel,nala,ala,cala, alamin,
     $       ratm, rlat, coskr,sinkr,
     $       hh, nbox,nrec,inrec, nrbox,nobox,inboxs, boxrad,boxctr,
     $       wf(ir00+1), rmesh, vbox, wkbox, wkbi,
     $       dphi(id1),dphi(id2),dphi(id3),  dphi(id4),
c -->        dphix-s   dphiy-s   dphiz-s     dphiuv-s
     $       dmat )
          endif
c
        else
c
c         Complex orbitals
c
          if( doblas3 )then
            call VSLOF2B( ibox, ndim, frcslo,strslo,
     $       norb,nk,nk0,nkloc, natm,ntyp, nshld,nald, nlat,
     $       itypa,numshl,lshel,nala,ala,cala, alamin,
     $       ratm, rlat, coskr,sinkr,
     $       hh, nbox,nrec,inrec, nrbox,nobox,inboxs, boxrad,boxctr,
     $       wvfcns(1,ir00+1), rmesh, vbox, wkbox, wkbi,
     $       dphi(id1),dphi(id2),dphi(id3),  dphi(id4),
c -->        dphix-s   dphiy-s   dphiz-s     dphiuv-s
     $       dmat, dmat(ipdm2) )
          else
            call VSLOF2( ibox, ndim, frcslo,strslo,
     $       norb,nk,nk0,nkloc, natm,ntyp, nshld,nald, nlat,
     $       itypa,numshl,lshel,nala,ala,cala, alamin,
     $       ratm, rlat, coskr,sinkr,
     $       hh, nbox,nrec,inrec, nrbox,nobox,inboxs, boxrad,boxctr,
     $       wvfcns(1,ir00+1), rmesh, vbox, wkbox, wkbi,
     $       dphi(id1),dphi(id2),dphi(id3),  dphi(id4),
c -->        dphix-s   dphiy-s   dphiz-s     dphiuv-s
     $       dmat )
          endif
c
        endif
c
c       Close box loop:
  900   continue
        ir00 = ir00 + nkloc*nobox*nrbox
 1000 continue
c
c The two*weight factor had originally been installed into dmat, now here.
      twowt = two*weight
      call DSCAL( 3*natm, twowt, frcslo, 1 )
      call DSCAL( 3*3   , twowt, strslo, 1 )
c
      if( nprocl .gt. 1 )then
c       Each processor has forces/stress from its part of grid.
c       Need to sum these into a single array, on master (node 0).
c       For now, also broadcast total back to all processors.
c
        nreduc = 3*natm
        call MPREDUC8( masterl, nreduc, frcslo, frctp , icomml )
C        call MPI_Reduce( frcslo, frctp , nreduc,
C     $   MPI_DOUBLE_PRECISION, MPI_SUM,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvslofrc - Reduce frcslo' )
        call DCOPY( nreduc, frctp,1, frcslo,1 )
c
        call MPBCAST8( masterl, nreduc, frcslo, icomml )
C        call MPI_Bcast( frcslo, nreduc, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
crpm        if( ierr.ne.0 ) call STOPXERR( 'pvslofrc - Bcast frcslo' )
c
        nreduc = 9
        call MPREDUC8( masterl, nreduc, strslo, strtp , icomml )
C        call MPI_Reduce( strslo, strtp , nreduc,
C     $   MPI_DOUBLE_PRECISION, MPI_SUM,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvslofrc - Reduce strslo' )
        call DCOPY( nreduc, strtp,1, strslo,1 )
c
        call MPBCAST8( masterl, nreduc, strslo, icomml )
C        call MPI_Bcast( strslo, nreduc, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvslofrc - Bcast strslo' )
c
c       End tp
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPDMATRD2
c
c
      subroutine KPDMATRD2( idmatfl, norb, ncplx, nk,nk0,nkloc,
     $ dmatrd,      dmat )                   
c---------------------------------------------------------------
c Purpose: pack complex dmat from r8 into full-c8 for vslofrc
c
c Written: Peter A. Schultz,  4-December-2008, for v2.63 (k-parallel)
c
c Revision history:
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c Input:
      INTEGER    idmatfl
      INTEGER    norb, ncplx,nk, nkloc,nk0
c Output:
      REAL       dmat
      DIMENSION  dmat(*)
c Scratch:
      DOUBLE PRECISION  dmatrd
      DIMENSION  dmatrd(*)
c
c Local declarations:
      INTEGER    nprocl ,iprocl ,masterl ,icomml
      INTEGER    nproc_k,iproc_k,master_k,icomm_k
      INTEGER    matrd,lenmsg, ik,ikloc, kmasterpr
c
c >>>> EXECUTABLE CODE:
c
c  Get Local MP info
      call MPNODES( nprocl )
      call MPNODE( iprocl )
      call MPNODE0( masterl )
      call MPCOMM( icomml )
c           
      call MPNODES_K( nproc_k )
      call MPNODE_K( iproc_k )
      call MPNODE0_K( master_k )
      call MPCOMM_K( icomm_k )
c       
      matrd = norb*norb
c
      ikloc = 0
      do 100 ik=1,nk
c
        lenmsg = matrd 
        call KPMASTER( nk,nprocl, ik,kmasterpr )
c        This might have been necessary once ...
C        if( nkloc .eq. nk ) kmasterpr = masterl
c       
        if( iprocl .eq. masterl )then
c         Image/local master reads dmat, and send to correct k-master proc
          call READBIG( idmatfl, matrd, dmatrd )
          if( kmasterpr .ne. masterl )then
c           This k-point goes off to kmasterpr
            call MPSENDR8( kmasterpr, lenmsg, dmatrd, icomml )
            goto 100
          endif 
        elseif( iprocl .eq. kmasterpr )then
c         Proc is master of this k-point, get dmatrd from masterl
          call MPRECVR8( masterl, lenmsg, dmatrd, icomml )
        else
c         This proc is not leading this k-point
          goto 100
        endif 
c  k-point master (kmasterpr) packs the data it got into its array.
        ikloc = ikloc + 1
c
        if( ncplx.eq.2 )then
c         Complex density matrices
          call VSLOFDM2( norb, nkloc, ikloc, dmatrd, dmat )
        else
c         Real density matrix
          call VSLOFDM1( norb, dmatrd, dmat )
        endif
c
  100 continue
c
c       Copy dmat from k-point masters to the other processors
c       processing the same k-points.
c
      if( nproc_k .gt. 1 .and.  nkloc .gt. 0 )then
        lenmsg = ncplx*nkloc*norb*norb
        call MPBCAST4( master_k, lenmsg, dmat, icomm_k )
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VSLOFDM2
c
c
      subroutine VSLOFDM2( norb, nkloc, ik, dmatrd, dmat )
c---------------------------------------------------------------
c Purpose: pack complex dmat from r8 into full-c8 for vslofrc
c
c Written: Peter A. Schultz,  4-December-2008, for v2.63 (k-parallel)
c
c Revision history:
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c Input:
      INTEGER  norb, nkloc,ik
      DOUBLE PRECISION  dmatrd
      DIMENSION  dmatrd(norb,norb)
c Output:
      REAL       dmat
      DIMENSION  dmat(2,nkloc,norb,norb)
c
c Local declarations:
      INTEGER    i,j
      REAL       r4zero
      DATA       r4zero /0.e0/
c
c >>>> EXECUTABLE CODE:
c
      do  j=1,norb
        do  i=1,j-1
          dmat(1,ik,i,j) =  dmatrd(j,i)
          dmat(2,ik,i,j) =  dmatrd(i,j)
        enddo
          dmat(1,ik,j,j) =  dmatrd(j,j)
          dmat(2,ik,j,j) =  r4zero
        do  i=j+1,norb
          dmat(1,ik,i,j) =  dmatrd(i,j)
          dmat(2,ik,i,j) = -dmatrd(j,i)
        enddo
      enddo
c
c    That's all Folks
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VSLOFDM1
c
c
      subroutine VSLOFDM1( norb, dmatrd, dmat )
c---------------------------------------------------------------
c Purpose: pack real dmat from r8 into full-r4 for vslofrc
c
c Written: Peter A. Schultz,  4-December-2008, for v2.63 (k-parallel)
c
c Revision history:
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c Input:
      INTEGER  norb
      DOUBLE PRECISION  dmatrd
      DIMENSION  dmatrd(norb,norb)
c Output:
      REAL       dmat
      DIMENSION  dmat(norb,norb)
c
c Local declarations:
      INTEGER    i,j
c
c >>>> EXECUTABLE CODE:
c
      do  j=1,norb
        do  i=1,j
          dmat(i,j) = dmatrd(j,i)
        enddo
        do  i=j+1,norb
          dmat(i,j) = dmatrd(i,j)
        enddo
      enddo
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
