c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PVSLOMAT
c
c
ctp:
      subroutine PVSLOMAT( ncplx,idwf,doblas3,
     $ norb,nk, natm,ntyp, nshld,nald, nlat,
     $ itypa,numshl,lshel,nala,ala,cala, alamin, ratm,rlat,coskr,sinkr,
     $ weight, hh, ndbox,nbox,nrec,inrec,noboxs,inboxs, boxrad,
     $ wvfcns,wf, vslobox, wkbox,wkbi,
     $ vijk, vkij, vslo )
c---------------------------------------------------------------
c Purpose: compute matrix elements of grid potentials
c
c Written: Peter A. Schultz,  3-June-1997, for v2.20 (boxes)
c
c Revision history:
c  13Feb06-PAS/2.60: mpi-clean
c   6Jul06-APT/    : merged task parallel and image parallel
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c  28Mar03-DMC: Many changes to make this a parallel version.
c     I need to distribute files which were only written on
c     processor 0 to all processors.  Also, this routine does
c     every n'th element on a processor, and then uses MPI_Reduce
c     and MPI_Bcast to accumulate the results on all processors.
c   2May01-PAS/2.47: remove disk storage of grid orbitals
c  10Jun99-PAS/2.36: blas3 implemention and flag
c  19May99-PAS/2.35: subcalls to separate real/complex
c  13Feb99-PAS/2.31: reduce op ct plus blas of bottom grid code
c   9Nov98-PAS/2.28: single precision grid orbitals/wf4
c  10Sep98-PAS/2.26: on-the-fly grid orbitals
c---------------------------------------------------------------
c
c Notes:
c  Assumes potential vslo has been shuffled into boxes on input
c
c  Note imaginary part of the integral is stored in v(j,i,k).
c  Thus im(v(i,j,k)) is minus what I store for (j,i), and is thus
c  the imaginary part of the integral of the potential with
c  conjg(phi-i)*phi-j, as is true for the rest of the Hamiltonian
c
c  Using subcalls to separate real/complex paths incurs a cost of
c  ~5% as opposed to pulling code up here in combined arithmetic.
c  Will pay this cost for purpose of clarity and compactness, and
c  another time worry about refining the code to get last bit
c  of speed.  17May99-PAS/2.35
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      REAL  wvfcns,wf, wkbi
      REAL  vkij
      REAL  vslobox
c
c input/output arrays:
      DIMENSION  vijk(norb,norb,nk)
c input arrays:
      LOGICAL    doblas3
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp),
     $ nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), ratm(3,natm), rlat(3,nlat)
      DIMENSION  coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  hh(3,3), ndbox(2,2,3),noboxs(nbox),inboxs(nrec,nbox)
      DIMENSION  wvfcns(2,*),wf(*)
      DIMENSION  vslo(*)
c scratch array
      DIMENSION  wkbi(*), wkbox(*), vslobox(*)
      DIMENSION  vkij(nk,norb,2*norb)
c local declarations:
      DIMENSION  r0box(3),boxctr(3)
c
c >>>> EXECUTABLE CODE:
c
      call MKSZERO( nk*norb*norb, vkij )
c
c  Get Local MP info
c
      call MPNODES( nprocs )
      call MPNODE0( master )
      call MPNODE( iproc )
      call MPCOMM( icomm )
c
c Initialize location of grid orbitals
      ir00 = 0
c
      irpot = 0
c
c Cycle over boxes of grid points
c
      do 1000 ibox=1,nbox
c
c Get basic box data
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
        if( MOD( ibox, nprocs ) .ne. iproc ) goto 900
c
c  Get grid orbitals:
c
        if( idwf.eq.0 )then
c
c         >>>> Compute grid orbitals in this box
c
          call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,boxctr )
          call MKRBOX( hh,r0box, idra,idrb,idrc, wkbox )
          ir00 = 0
          iwk1 = 1
          iwk2 = iwk1 + 3*nrbox
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
     $     wkbox      ,wkbox(iwk2),wkbox(iwk3),wkbox(iwk4),
     $     wkbox(iwk5),wkbox(iwk6),wkbox(iwk7),wkbox(iwk8) )
c
        endif
c
c Put r*8 potential into r*4 space for use below:
        do  nr=1,nrbox
          vslobox(nr) = vslo(iv0+nr)
        enddo
c
c  Compute matrix elements in this box:
c
        if( ncplx.eq.2 )then
c
c >>>> Complex orbitals
c
          if( doblas3 )then
            call VSLOM2B( ibox,
     $       norb,nk, natm,ntyp, nshld,  itypa,numshl,lshel,
     $       nbox,nrec,inrec, nrbox, nobox,inboxs,
     $       wvfcns(1,ir00+1), vslobox, wkbi,
     $       vkij, vkij(1,1,norb+1) )
          else
            call VSLOM2( ibox,
     $       norb,nk, natm,ntyp, nshld,  itypa,numshl,lshel,
     $       nbox,nrec,inrec, nrbox, nobox,inboxs,
     $       wvfcns(1,ir00+1), vslobox, wkbi,
     $       vkij )
          endif
c
        elseif( ncplx.eq.1 )then
c
c >>>> Real orbitals
c
          if( doblas3 )then
            call VSLOM1B( ibox,
     $       norb, natm,ntyp, nshld,  itypa,numshl,lshel,
     $       nbox,nrec,inrec, nrbox, nobox,inboxs,
     $       wf(ir00+1), vslobox, wkbi,
     $       vkij, vkij(1,1,norb+1) )
          else
            call VSLOM1( ibox,
     $       norb, natm,ntyp, nshld,  itypa,numshl,lshel,
     $       nbox,nrec,inrec, nrbox, nobox,inboxs,
     $       wf(ir00+1), vslobox, wkbi,
     $       vkij )
          endif
c
        endif
c
  900   continue
c
        ir00 = ir00 + nk*nobox*nrbox
 1000 continue
ctp:start
c
      if( nprocs .gt. 1 )then
c       Each local processor now has dHam matrix for its part of grid.
c       Need to sum these into a single array, on master (node 0).
c       For now, also broadcast total back to all processors.
c
        nmat = nk*norb*norb
        call MPREDUC4( master, nmat, vkij, vkij(1,1,norb+1), icomm )
C        call MPI_Reduce( vkij, vkij(1,1,norb+1), nmat, MPI_REAL,
C     $   MPI_SUM, master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pvslomat - Reduce vkij' )
c
        call SCOPY( nmat, vkij(1,1,norb+1),1, vkij,1 )
c
C        Skip the broadcast to all the other procs ...
C        call MPBCAST4( master, nmat, vkij, icomm )
      endif
ctp:end
c
c  Load grid matrix elements into Hamiltonian
c
      do  k=1,nk
        do  j=1,norb
          do  i=1,norb
            vijk(i,j,k) = vijk(i,j,k) + weight*vkij(k,i,j)
          enddo
        enddo
      enddo
c
c    That's all Folks!
c
      RETURN
      END
