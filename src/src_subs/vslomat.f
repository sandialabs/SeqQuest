c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VSLOMAT
c
c
      subroutine VSLOMAT( ncplx,idwf,doblas3,
     $ norb,nk, natm,ntyp, nshld,nald, nlat,
     $ itypa,numshl,lshel,nala,ala,cala, alamin, ratm,rlat,
     $ coskr,sinkr,
     $ weight, hh, ndbox,nbox,nrec,inrec,noboxs,inboxs, boxrad,
     $ wvfcns, wf, vslobox, wkbox, wkbi,
     $ vijk,         vkij,          vslo )
c      R8(mat_nkloc) R4(2mat*nkloc) R8(nptr)
c      =nmatloc   =nmatloc                = nptr
c---------------------------------------------------------------
c Purpose: compute matrix elements of grid potentials
c
c Written: Peter A. Schultz,  3-June-1997, for v2.20 (boxes)
c
c Revision history:
c  17Nov08-PAS/2.63: k-parallel scf matrices
c   4Jan08-PAS/2.61: reduced (dist) memory for grid wvfcns
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
      REAL  weightr4
      REAL  vslobox
c
c Input/Output arrays:
      DIMENSION  vijk(norb,norb,*)
c                vijk(norb,norb,nkloc)
c Input arrays:
      LOGICAL    doblas3
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp),
     $ nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), ratm(3,natm), rlat(3,nlat)
      DIMENSION  coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  hh(3,3), ndbox(2,2,3),noboxs(nbox),inboxs(nrec,nbox)
      DIMENSION  wvfcns(2,*),wf(*)
      DIMENSION  vslo(*)
c Scratch array
      DIMENSION  wkbi(*), wkbox(*), vslobox(*)
      DIMENSION  vkij(*)
c                vkij(nkloc,norb,2*norb)
c
c Local declarations:
      DIMENSION  r0box(3),boxctr(3)
c
c >>>> EXECUTABLE CODE:
c
c  Get Local MP info - note we are working at the k-level
      call MPNODES_K( nprocs_k )
      call MPNODE0_K( master_k )
      call MPNODE_K( iproc_k )
      call MPCOMM_K( icomm_k )
c
c  Set up kparallel
      call KPMINE( nk, nkloc, nk0 )
c
      call MKSZERO( nkloc*norb*norb, vkij )
      if( nkloc .eq. 0 ) goto 1001

      ipvkij2 = nkloc*norb*norb + 1
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
c       Get basic box data
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
c       Each processor will do every nprocs_k'th box in this loop.
        if( MOD( ibox, nprocs_k ) .ne. iproc_k ) goto 1000
c
c       Get grid orbitals:
c
        if( idwf.eq.0 )then
c         Compute orbitals on mesh in this box
c
c         Generate mesh point coordinates
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
     $     nk,nk0,nkloc,  natm,ntyp, nshld,nald, nlat,
     $     itypa,numshl,lshel,nala,ala,cala,alamin,
     $     ratm,rlat,coskr,sinkr,
     $     nobox,nrbox,nrec,inrec,inboxs, boxrad,boxctr,
     $     wvfcns(1,ir00+1), wf(ir00+1),
     $     wkbox      ,wkbox(iwk2),wkbox(iwk3),wkbox(iwk4),
     $     wkbox(iwk5),wkbox(iwk6),wkbox(iwk7),wkbox(iwk8) )
c
        endif
c
c       Put r*8 potential into r*4 space for use below:
        do  nr=1,nrbox
          vslobox(nr) = weight*vslo(iv0+nr)
cwt          vslobox(nr) = vslo(iv0+nr)
        enddo
c
c       ***************************************************
c       Compute contribution to matrix elements in this box
c       ***************************************************
c
        if( ncplx.eq.2 )then
c
c         Complex orbitals
c
          if( doblas3 )then
            call VSLOM2B( ibox,
     $       norb,nkloc, natm,ntyp, nshld,  itypa,numshl,lshel,
     $       nbox,nrec,inrec, nrbox, nobox,inboxs,
     $       wvfcns(1,ir00+1), vslobox, wkbi,
     $       vkij, vkij(ipvkij2) )
c            R4(nkloc*mat)   R4(2*noboxmax*noboxmax*nk)
c            (nmatloc+1)/2    nk*noboxmax**2
          else
            call VSLOM2( ibox,
     $       norb,nkloc, natm,ntyp, nshld,  itypa,numshl,lshel,
     $       nbox,nrec,inrec, nrbox, nobox,inboxs,
     $       wvfcns(1,ir00+1), vslobox, wkbi,
     $       vkij )
          endif
c
        elseif( ncplx.eq.1 )then
c
c         Real orbitals
c
          if( doblas3 )then
            call VSLOM1B( ibox,
     $       norb, natm,ntyp, nshld,  itypa,numshl,lshel,
     $       nbox,nrec,inrec, nrbox, nobox,inboxs,
     $       wf(ir00+1), vslobox, wkbi,
     $       vkij, vkij(ipvkij2) )
c            R4(mat)  R4(noboxmax**2)
          else
            call VSLOM1( ibox,
     $       norb, natm,ntyp, nshld,  itypa,numshl,lshel,
     $       nbox,nrec,inrec, nrbox, nobox,inboxs,
     $       wf(ir00+1), vslobox, wkbi,
     $       vkij )
c            R4(mat)
          endif
c
        endif
c
c       Close box loop:
  900   continue
        ir00 = ir00 + nkloc*nobox*nrbox
 1000 continue
c     Bypass for nkloc = 0:
 1001 continue
c 
      if( nprocs_k .gt. 1 )then
c       Each local processor now has a partial dHam matrix built from its boxes
c       Need to sum these into a single array, on master_k (node 0).
c
        nmatloc = nkloc*norb*norb
        call MPREDUC4( master_k, nmatloc,
     $       vkij,    vkij(ipvkij2),    icomm_k )
c -----> nmatloc(r*4) nmatloc(r*4)
c
        call SCOPY( nmatloc, vkij(ipvkij2),1, vkij,1 )
c
      endif
c
c  Load grid matrix elements into Hamiltonian only on kmasters
c
       if( iproc_k.eq.master_k )then
c       kmaster adds vslomat r4(kij) matrix into total r8(ijk) matrix:
        call VSLOMLOAD( norb, nkloc, vijk, vkij )
cwt        call VSLOMLOAD( norb, nkloc, weight, vijk, vkij )
       endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VSLOMLOAD
c
c
      subroutine VSLOMLOAD( norb, nkloc, vijk, vkij )
cwt      subroutine VSLOMLOAD( norb, nkloc, weight, vijk, vkij )
c---------------------------------------------------------------
c Purpose: load computed r4 vslomat into total r8 Ham
c
c Written: Peter A. Schultz,  4-December-20-8, for 2.62 (k-parallel)
c
c Revision history:
c   4Dec08-PAS/2.62: pushed down from vslomat
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c Input:
      INTEGER    nkloc,norb
      REAL       vkij
      DIMENSION  vkij(nkloc,norb,norb)
c Input/Output:
      DOUBLE PRECISION  vijk
cwt      DOUBLE PRECISION  weight
      DIMENSION  vijk(norb,norb,nkloc)
c
c Local declarations:
      INTEGER    ikloc,i,j
c
c >>>> EXECUTABLE CODE:
c
      do  ikloc=1,nkloc
        do  j=1,norb
          do  i=1,norb
            vijk(i,j,ikloc) = vijk(i,j,ikloc) + vkij(ikloc,i,j)
cwt            vijk(i,j,ikloc) = vijk(i,j,ikloc) + weight*vkij(ikloc,i,j)
          enddo
        enddo
      enddo
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
