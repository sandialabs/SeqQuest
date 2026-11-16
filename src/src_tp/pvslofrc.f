c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PVSLOFRC
c
c
      subroutine PVSLOFRC( idmatfl, ndim,ncplx,idwf,doblas3,
     $ nafrc,iforce, frcslo,strslo,
     $ norb,nk, natm,ntyp, nshld,nald, nlat,
     $ itypa,numshl,lshel,nala,ala,cala, alamin,
     $ ratm, rlat, wtk,coskr,sinkr,
     $ weight, hh, ndbox,nbox,nrec,inrec,noboxs,inboxs, boxrad,
     $ rmesh, vbox, wkbox, frctp,
     $ wvfcns,wf,
     $ wkbi, dphi,
     $ dmatrd, dmat,  vslo )
c---------------------------------------------------------------
c Purpose: compute force/stress from grid potentials
c
c Written: Peter A. Schultz,  3-June-1997, for v2.20 (boxes)
c
c Revision history:
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
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      REAL  dmat
      REAL  wvfcns,wf, wkbi
      REAL  vbox
      REAL  dphi
c
c output arrays:
      DIMENSION  frcslo(3,natm), strslo(3,3)
c parallel scratch:
      DIMENSION  frctp(3,natm), strtp(3,3)
c input arrays:
      LOGICAL    doblas3
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp),
     $ nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp)
      DIMENSION  ratm(3,natm), rlat(3,nlat)
      DIMENSION  wtk(nk),coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  hh(3), ndbox(2,2,3),noboxs(nbox),inboxs(nrec,nbox)
      DIMENSION  wvfcns(2,*),wf(*)
      DIMENSION  vslo(*)
c scratch arrays:
      DIMENSION  wkbi(*), wkbox(*)
      DIMENSION  dmatrd(norb,norb)
      DIMENSION  dmat(ncplx*nk,norb,norb)
      DIMENSION  rmesh(3,*)
      DIMENSION  vbox(*)
      DIMENSION  dphi(*)
c local declarations:
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
ctp:
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
c  Get and sort density matrix into k-first form
c
      twowt = two*weight
      do 15 k=1,nk
        if( iproc.eq.master ) call READBIG( idmatfl, norb*norb, dmatrd )
ctp:
        if( nprocs.gt.1 ) call MPBCAST8( master, norb*norb, dmatrd, 
     $       icomm )
c
        if( ncplx.eq.2 )then
c         Complex density matrix
          kr = 2*(k-1) + 1
          ki = kr + 1
          do  j=1,norb
            do  i=1,j-1
              dmat(kr,i,j) =  twowt*dmatrd(j,i)
              dmat(ki,i,j) =  twowt*dmatrd(i,j)
            enddo
              dmat(kr,j,j) =  twowt*dmatrd(j,j)
              dmat(ki,j,j) =  r4zero
            do  i=j+1,norb
              dmat(kr,i,j) =  twowt*dmatrd(i,j)
              dmat(ki,i,j) = -twowt*dmatrd(j,i)
            enddo
          enddo
        else
c         Real density matrix
          do  j=1,norb
            do  i=1,j
              dmat(1,i,j) = twowt*dmatrd(j,i)
            enddo
            do  i=j+1,norb
              dmat(1,i,j) = twowt*dmatrd(i,j)
            enddo
          enddo
        endif
   15 continue
c
c Initialize location of grid orbitals
      ir00 = 0
c
      irpot = 0
c
c  Cycle over boxes of grid points
c
      do 1000 ibox=1,nbox
c
c  Get basic box info
        call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $   ira0,irb0,irc0, idra,idrb,idrc )
        nrbox = idra*idrb*idrc
c
        iv0 = irpot
        irpot = irpot + nrbox
c
        nobox = noboxs(ibox)
        if( nobox.le.0 ) goto 1000
c
c       Each processor will do every nprocs element in this loop.
        if( iproc .ne. MOD( ibox, nprocs ) ) goto 900
c
c  Generate coordinates
c
        call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
        call MKRBOX( hh,r0box, idra,idrb,idrc, rmesh )
c
c  Get grid orbitals:
c
        if( idwf.eq.0 )then
c
c         >>>> Compute grid orbitals in this box
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
     $     nk, natm,ntyp, nshld,nald, nlat,
     $     itypa,numshl,lshel,nala,ala,cala,alamin,
     $     ratm,rlat,coskr,sinkr,
     $     nobox,nrbox,nrec,inrec,inboxs, boxrad,boxctr,
     $     wvfcns(1,ir00+1), wf(ir00+1),
     $     rmesh      ,wkbox(iwk2),wkbox(iwk3),wkbox(iwk4),
     $     wkbox(iwk5),wkbox(iwk6),wkbox(iwk7),wkbox(iwk8) )
c
        endif
c
c Put r*8 potential into r*4 space for use below:
        do  nr=1,nrbox
          vbox(nr) = vslo(iv0+nr)
        enddo
c
c  Carve up space for derivatives of grid wave functions:
        nrkbox = ncplx*nk*nrbox
        id1 = 1
        id2 = id1 + 1*nrkbox
        id3 = id1 + 2*nrkbox
        id4 = id1 + 3*nrkbox
c
        if( ncplx.eq.1 )then
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
     $       dmat, dmat(1,1,norb+1) )
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
          if( doblas3 )then
            call VSLOF2B( ibox, ndim, frcslo,strslo,
     $       norb,nk, natm,ntyp, nshld,nald, nlat,
     $       itypa,numshl,lshel,nala,ala,cala, alamin,
     $       ratm, rlat, coskr,sinkr,
     $       hh, nbox,nrec,inrec, nrbox,nobox,inboxs, boxrad,boxctr,
     $       wvfcns(1,ir00+1), rmesh, vbox, wkbox, wkbi,
     $       dphi(id1),dphi(id2),dphi(id3),  dphi(id4),
c -->        dphix-s   dphiy-s   dphiz-s     dphiuv-s
     $       dmat, dmat(1,1,norb+1) )
          else
            call VSLOF2( ibox, ndim, frcslo,strslo,
     $       norb,nk, natm,ntyp, nshld,nald, nlat,
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
c    Close box loop:
c
  900   continue
c
        ir00 = ir00 + nk*nobox*nrbox
 1000 continue
ctp:start
c
      if( nprocs .gt. 1 )then
c       Each processor has forces/stress from its part of grid.
c       Need to sum these into a single array, on master (node 0).
c       For now, also broadcast total back to all processors.
c
        nreduc = 3*natm
        call MPREDUC8( master, nreduc, frcslo, frctp , icomm )
C        call MPI_Reduce( frcslo, frctp , nreduc,
C     $   MPI_DOUBLE_PRECISION, MPI_SUM,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvslofrc - Reduce frcslo' )
        call DCOPY( nreduc, frctp,1, frcslo,1 )
c
        call MPBCAST8( master, nreduc, frcslo, icomm )
C        call MPI_Bcast( frcslo, nreduc, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
crpm        if( ierr.ne.0 ) call STOPXERR( 'pvslofrc - Bcast frcslo' )
c
        nreduc = 9
        call MPREDUC8( master, nreduc, strslo, strtp , icomm )
C        call MPI_Reduce( strslo, strtp , nreduc,
C     $   MPI_DOUBLE_PRECISION, MPI_SUM,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvslofrc - Reduce strslo' )
        call DCOPY( nreduc, strtp,1, strslo,1 )
c
        call MPBCAST8( master, nreduc, strslo, icomm )
C        call MPI_Bcast( strslo, nreduc, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvslofrc - Bcast strslo' )
c
c       End tp
      endif
ctp:end
c
c    That's all Folks!
c
      RETURN
      END

