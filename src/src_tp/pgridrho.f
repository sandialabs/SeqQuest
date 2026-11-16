c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDRHO
c
c
ctp:
      subroutine PGRIDRHO( ivecfl,idmatfl,
     $ irhoopt,ncplx,idwf,doblas3,
     $ nstate,npop, norb,nk, natm,ntyp, nshld,nald, nlat,
     $ itypa,numshl,lshel,nala,ala,cala, alamin, ratm,rlat,
     $ wtk, coskr,sinkr, eigpop, rtpop,
     $ weight, hh, ndbox,nbox,nrec,inrec,noboxs,inboxs, boxrad,
     $ wvfcns,wf, rhobox, wkbox,wkbi,
     $ dmat, dmatrd, rho,  nptr,rho_tp  )
c---------------------------------------------------------------
c Purpose: computes electron density in coordinate space
c
c Written: Peter A. Schultz,  3-June-1997, for v2.20 (boxes)
c
c Revision history:
c  11Feb07-PAS/2.60: mpi-clean, and use wrapped reduce
c   6Jul06-APT/    : merged task parallel and image parallel
c  30Jul05-PAS/2.59: master-only matrix files
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c  24Mar03-DMC: Many changes to make this a parallel version.  I need
c     to distribute files which were only written on processor 0
c     to all processors.  Also, this routine does every n'th element on
c     a processor, and then uses MPI_Reduce and MPI_Bcast to accumulate
c     results on all processors.
c  20Jul01-PAS/2.49: spin-polarized dft, compact dmat, eigfcn
c  21Jun01-PAS/2.48: replace STOPs
c   2May01-PAS/2.47: remove disk storage of grid orbitals
c  10Jun99-PAS/2.36: blas3 installation
c  13May99-PAS/2.35: subcalls to real/cmplx; eigenfcn density
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  13Feb99-PAS/2.31: reduce op ct plus blas of bottom grid code
c   9Nov98-PAS/2.28: single precision grid orbitals/wf4
c  10Sep98-PAS/2.26: on-the-fly grid orbitals
c---------------------------------------------------------------
c
c  Dave Raczkowski suggested the eigenfunction (ef) scheme to
c  compute the density rather than the density matrix (dm) method,
c  faster for certain smaller problems.
c
c  Using subcalls to separate real/complex & dm/ef-based density
c  schemes under the box loop incurs ~5% cost over having that
c  code pulled up into this routine with combined real/cmplx
c  arithmetic/memory indexing.  For the moment, will pay this
c  cost for the purpose of clarity, and perhaps later
c  restructure the work under the ibox loop to get last little
c  bit of speed.  17May99-PAS
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      REAL  wvfcns,wf, wkbi
      REAL  dmat
      REAL  rhobox
c
c input arrays:
      LOGICAL    doblas3
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp),
     $ nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), ratm(3,natm), rlat(3,nlat)
      DIMENSION  wtk(nk), coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  eigpop(nstate,nk)
      DIMENSION  hh(3,3), ndbox(2,2,3),noboxs(nbox),inboxs(nrec,nbox)
      DIMENSION  wvfcns(2,*),wf(*)
c output array:
      DIMENSION  rho(*)
ctp:start
c parallel scratch:
      DIMENSION  rho_tp(*)
ctp:end
c scratch arrays:
      DIMENSION  rtpop(norb)
      DIMENSION  wkbi(*), wkbox(*), rhobox(*)
      DIMENSION  dmat(ncplx*nk,norb,norb), dmatrd(2*norb*norb)
c local declarations:
      DIMENSION  r0box(3),boxctr(3)
      DATA  zero,two / 0.d0,2.d0 /
      DATA  r4zero,r4two / 0.e0,2.e0 /
c
c >>>> EXECUTABLE CODE:
c
c
c  Get Local MP info
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
c Get dmat/eigvec:
      if( irhoopt .eq. 1 )then
c       Density-matrix based grid density, set up dmat(k,i,j)
c
ctp:
        call PRHDMAT( idmatfl, ncplx,norb,nk,
CTP        call RHDMAT( idmatfl, ncplx,norb,nk,
     $   dmat, dmatrd )
c
      elseif( irhoopt .eq. 2 )then
c       Eigenfunction-based grid density, set up eigvec(nc,n,k,i)
c
        numpop = npop
ctp:
        call PRHEIGVEC( ivecfl, ncplx,norb,nstate,numpop,nk,
CTP        call RHEIGVEC( ivecfl, ncplx,norb,nstate,numpop,nk,
     $   wtk,eigpop, rtpop,
     $   dmat,    dmatrd )
c -->    eigvec-o eigvrd-s
c
      else
        call STOPXERR( 'irhoopt is corrupted in gridrho' )
      endif
c
c Initialize location of grid orbitals in big storage
      ir00 = 0
c
      irpot = 0
c
c  Cycle over boxes of grid points
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
c Get number of orbitals in this box
        nobox = noboxs(ibox)
c
        if( nobox.le.0 ) goto 1000
ctp:
c
c       Each processor will do every nprocs element in this loop.
        if( iproc .ne. MOD( ibox, nprocs ) ) goto 900
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
c Clear box density:
        call MKSZERO( nrbox, rhobox )
c
c  Compute the density in this box:
c
        if( irhoopt.eq.1 )then
c
c >>>> Use density matrix route to compute density
c
          if( ncplx.eq.2 )then
c
c           Complex orbitals/DM scheme
c
            if( doblas3 )then
              call GRHODM2B( ibox,
     $         norb,nk, natm,ntyp, nshld,  itypa,numshl,lshel,
     $         nbox,nrec,inrec, nrbox, nobox,inboxs,
     $         wvfcns(1,ir00+1), rhobox, wkbi,
     $         dmat, dmat(1,1,norb+1) )
            else
              call GRHODM2( ibox,
     $         norb,nk, natm,ntyp, nshld,  itypa,numshl,lshel,
     $         nbox,nrec,inrec, nrbox, nobox,inboxs,
     $         wvfcns(1,ir00+1), rhobox, wkbi,
     $         dmat )
            endif
c
          elseif( ncplx.eq.1 )then
c
c           Real orbitals/DM scheme
c
            if( doblas3 )then
              call GRHODM1B( ibox,
     $         norb, natm,ntyp, nshld,  itypa,numshl,lshel,
     $         nbox,nrec,inrec, nrbox, nobox,inboxs,
     $         wf(ir00+1), rhobox, wkbi,
     $         dmat, dmat(1,1,norb+1) )
            else
              call GRHODM1( ibox,
     $         norb, natm,ntyp, nshld,  itypa,numshl,lshel,
     $         nbox,nrec,inrec, nrbox, nobox,inboxs,
     $         wf(ir00+1), rhobox, wkbi,
     $         dmat )
            endif
c
          endif
c
        elseif( irhoopt.eq.2 )then
c
c >>>> Use eigenfunction route to compute density
c
          if( ncplx.eq.2 )then
c
c           Complex orbitals/EF scheme
c
            if( doblas3 )then
              call GRHOEF2B( ibox, numpop,
     $         norb,nk, natm,ntyp, nshld,  itypa,numshl,lshel,
     $         nbox,nrec,inrec, nrbox, nobox,inboxs,
     $         wvfcns(1,ir00+1), rhobox, wkbi,
     $         dmat, dmat(1,1,norb+1) )
            else
              call GRHOEF2( ibox, numpop,
     $         norb,nk, natm,ntyp, nshld,  itypa,numshl,lshel,
     $         nbox,nrec,inrec, nrbox, nobox,inboxs,
     $         wvfcns(1,ir00+1), rhobox, wkbi,
     $         dmat )
            endif
c
          elseif( ncplx.eq.1 )then
c
c           Real orbitals/EF scheme
c
            if( doblas3 )then
              call GRHOEF1B( ibox, numpop,
     $         norb, natm,ntyp, nshld,  itypa,numshl,lshel,
     $         nbox,nrec,inrec, nrbox, nobox,inboxs,
     $         wf(ir00+1), rhobox, wkbi,
     $         dmat, dmat(1,1,numpop+1) )
            else
              call GRHOEF1( ibox, numpop,
     $         norb, natm,ntyp, nshld,  itypa,numshl,lshel,
     $         nbox,nrec,inrec, nrbox, nobox,inboxs,
     $         wf(ir00+1), rhobox, wkbi,
     $         dmat )
            endif
c
          endif
c
        endif
c
c Update total density with box density:
        do  nr=1,nrbox
          rho(iv0+nr) = rho(iv0+nr) + rhobox(nr)
        enddo
c
  900   continue
        ir00 = ir00 + nk*nobox*nrbox
 1000 continue
ctp:start
c
      if( nprocs .gt. 1 )then
c       Each local processor now has density for its boxes.
c       Need to sum these into a single array, on master (node 0).
c       For now, also broadcast total back to all processors.
c
        call MPREDUC8( master, nptr, rho, rho_tp, icomm )
C        call MPI_Reduce( rho, rho_tp, nptr, MPI_DOUBLE_PRECISION,
C     $   MPI_SUM, master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pgridrho - Reduce rho' )
        call DCOPY( nptr, rho_tp,1, rho,1 )
c
        call MPBCAST8( master, nptr, rho, icomm )
C        call MPI_Bcast( rho, nptr, MPI_DOUBLE_PRECISION,
C     $   master, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'pgridrho - Bcast rho' )
       endif
ctp:end
c
c    That's all Folks!
c
      RETURN
      END
