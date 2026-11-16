c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PEIGSOLVKP
c
c
      subroutine PEIGSOLVKP( masterlist,
     $ icluster, gapp,
     $ do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $ ikoopman,ioccdft, iflcdft,
     $ idhamfl,iham0fl,iovlpfl,
     $ ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $ icloseocc,dokdefct,doblas3, ikdefct,nstbulk,
     $ etemp,edegen, elecno, efermi,egap, esumsp,
     $ norb,nstate, nk, veck,wtk, eigval,eigpop,npop, wksml,
     $ wmat )
c---------------------------------------------------------------
c Purpose: k-parallel-parallel solves for complex density matrices,
c          hamiltonian and overlap matrices of current iteration
c
c Written: Andrew C. Pineda, 16-May-2012, for 2.64
c          Adapted (forked) from  peigsolv
c
c Revision history:
c   1Dec22-PAS/2.69: koopmans/occdft
c   4Dec13-PAS/2.64: cosmetics, un-drivered, integrated
c---------------------------------------------------------------
c
c  On input, idhamfl contains full input delta-hamiltonian.
c  To obtain full hamiltonian, delta-H needs to be combined with
c  H0 (computed in setup and saved in iham0fl; overlap in iovlpfl),
c  done in routine "schroed", which also solves the Schroedinger
c  eqn. Routine "pop" populates the resulting levels and finds the
c  fermi energy, then "rhoij" takes the eigenfcns and occupations
c  and builds and saves to disk the density matrix for each k-vector.
c     
c  When this routine is done:
c   (1) eigenvectors are saved to "ivecfl" (if do_eigvecs)
c   (2) density matrix is saved to "idmat" (and if do_eigpops)
c   (3) energy-weighted density matrix is saved to "iemat"
c   (4) eigenvalue spectrum is returned in "eigval()" (not written)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
ctp: Pull in the mpi stuff
      INCLUDE  'mpif.h'
c Input:
      LOGICAL    do_eigvecs, do_eigpops
      LOGICAL    dokdefct
      LOGICAL    doblas3
      DIMENSION  veck(3,nk),wtk(nk)
      DIMENSION  elecno(2)
      DIMENSION  nstbulk(2)
c Output:
      DIMENSION  efermi(2),egap(2)
      DIMENSION  eigval(nstate*nk,nspin),eigpop(nstate*nk,nspin)
c Scratch:
      DIMENSION  wksml(norb,*)
      DIMENSION  wmat(*)
c Scratch/parallel:
      INTEGER    icluster(*)
c ACP/space for list of kp-masters used in pschroedkp
      INTEGER    masterlist(*)
      DIMENSION  gapp(*)
      DATA       itimlog / 5 /
c
c Local:
      DIMENSION  occlvl(2)
      DATA       occlvl / 2.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
 9001 format(a,6i8)
      if( IWR.gt.0 )then
        write(IWR,*) 'PEIGSOLVKP: k-parallel-parallel solver'
C        write(IWR,*) 'KPP/DEV: smaller prhoij(06) stripe/transpose'
C        write(IWR,*) 'DEV/PEIGSOLVKP: ikoop,icdft=',ikoopman,ioccdft
C        write(IWR,9001) 'solvkp: norb,nstate,nk=',norb,nstate,nk
      endif
c  Solve Schroedinger equation for the various k-vectors
c
      mat = norb*norb
      nmat = nk*mat
c
c Partition space within wmat() for schroed:
      iw1 = 1
      iww = iw1 + mat
      nwsch = norb*(4*norb + 2*nstate)
      iwx = iww + nwsch
c
      call MPCOMM( icomm )
      call MPNODES( nproc )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      if( iproc .eq. masterp )then
c       Go to start of idhamfl, to get ready for dHam reads
        REWIND( unit=idhamfl )
c       Go to start of ivecfl, to get ready for eigvec writes
        REWIND( unit=ivecfl )
      endif
c
      if( itimlog.gt.0 ) call TIMER( 'eigsolvkp/begin      ' )
c
c      if( nproc.eq.1 .or. nk.gt.1 )then
c
      if( ikoopman.gt.0 )then
c       Skip solvers if koopmans is on
        goto 2000
      endif
      if( nproc.eq.1) then
c$$$      if( nproc.eq.1 .or. nk.gt.4 .or. norb.lt.800 .or.
c$$$     $    nproc.lt.(4*nk) )then
c       Parallel over k-points - single-proc solves
c
        do  ispin=1,nspin
c
          call PSCHROED2( do_eigvecs,
     $     idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $     norb,nstate,nk, eigval(1,ispin), wksml,    wksml(1,9),
c -->                                       rwork-7ns seval(n)-s
     $     scut,
     $     wmat(iw1), wmat(iww), wmat(iwx) )
c -->      vmat-s     wmat-6ks   iwork-6ns
c
          if( itimlog.gt.1 ) call TIMER( 'eigsolvkp/psch2 k-par' )
        enddo
c
      else
c       Parallel solvers
c
        call PSCHROEDKP( do_eigvecs,
     $   idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $   norb,nstate,nk, eigval, wksml,
c -->                            wksml-ns
     $   nspin, icluster, gapp,
     $   wmat, masterlist )
c -->    wmat-s
c
        if( itimlog.gt.1 ) call TIMER( 'eigsolvkp/pschkp scalap' )
      endif
c
      if( .not. do_eigpops )then
c       Skip populations and density matrix construction
        RETURN
      endif
c
c  Evaluate density matrix in diagonal (eigfcn) representation
c
 2000 continue
c
      nstk = nstate*nk
c
      npop = 0
      esumsp = zero
      occfull = occlvl(nspin)
c
cdev/ckoop: get koopmans eigvals
c If Koopmans: skipped solver, need to get data from baseline (eigs/vecs)
c   Get eigenvalues from koopman baseline
c   Get spink of koopmans baseline, to know how convert vecs for this job
c
      nspink = nspin
      if( ikoopman .gt. 0 )then
c       With koopmans on, skipped solvers, therefore retrieve old eigenvalues
c       Not strictly true, it worked anyway, but logic is cleaner AND
c         this gives a verification of which koopmans source was input
        if( iproc .eq. masterp )then
c         Master reads this from a file ...
          call FLOPENA( ieofile, 'eigkoop' )
          call RDEO( ieofile, 'eigval', 2, nspink, nspin,nk,nstate,
     $               eigval )
          call FLCLOSE( ieofile )
        endif
c        ... and tells all the other processors
        call MPI_Bcast( eigval(1,1), nstate*nk*nspin,
     $   MPI_DOUBLE_PRECISION,
     $   masterp, icomm, ierr )
        if( ierr.ne.0 ) call STOPXERR( 'eigsolv - Bcast eigval' )
      endif
c
      if( ioccdft .gt. 0 )then
c       Use occdft constrained (preset) occupations from file
        if( iproc .eq. masterp )then
c         Master reads pre-set occupations from file ...
          REWIND( unit=iflcdft )
          call RDEO( iflcdft, 'occups', ioccdft,
     $     mspin,nspin, nk, nstate,
     $     eigpop )
        endif
c        ... and tells all the other processors:
        call MPI_Bcast( eigpop(1,1), nstate*nk*nspin,
     $   MPI_DOUBLE_PRECISION,
     $   masterp, icomm, ierr )
      endif
c
      do  ispin=1,nspin
        elnumbr = elecno(ispin)
c
        if( ioccdft.gt.0 )then
c         Using pre-set occdft occupations
          call POPOCCDFT(  nstate,nk,ispin,
     $     elnumbr, numpop, occfull,
     $     efermi(ispin), egap(ispin), esum1e,
     $     wtk, eigval(1,ispin), eigpop(1,ispin) )

        elseif( ( icloseocc .gt. 0 ) .or. dokdefct )then
c         Use closed-shell occupations, or defect sampling
          call POPCLOSE( icloseocc,dokdefct,ikdefct,
     $     efermi(ispin),egap(ispin), nstbulk(ispin),
     $     etemp,edegen, elnumbr,esum1e, occfull,
     $     nstate,nk, wtk, eigval(1,ispin),eigpop(1,ispin),numpop,
     $     wmat )
c -->      wlvl-s(nstate)
c
        else
c         Use standard metallic occupations (default)
          mspin = 1
          call POP( mspin, popnumbr,
     $     efermi(ispin),egap(ispin),
     $     etemp,        elnumbr,esum1e, occfull,
     $     nstate,nk, wtk, eigval(1,ispin),eigpop(1,ispin),numpop,
     $     wmat(1), wmat(1+nstk),wmat(1+2*nstk),
     $     wmat(1+3*nstk), wmat(1+4*nstk) )
c -->      elvl wlvl klvl nlvl flvl - all scratch (nk*state)
        endif
c
        npop = MAX( npop,numpop )
        esumsp = esumsp + esum1e
      enddo
c
      if( itimlog.gt.0 ) call TIMER( 'eigsolvkp/pop        ' )
c
      if( .not. do_eigvecs )then
c       Without eigvecs, skip density matrix construction
        RETURN
      endif
c
c  Evaluate density matrix in orbital basis
c
      mat = norb*norb
c
      call MPNODES( nproc )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      if( iproc .eq. masterp )then
        REWIND( unit=ivecfl )
        REWIND( unit=idmatfl )
        REWIND( unit=iematfl )
        if( nspin.eq. 2 )then
          REWIND( unit=idmatsfl )
          REWIND( unit=iematsfl )
        endif
      endif
      iw1 = 1
      iw2 = iw1 + mat
      iw3 = iw2 + mat
c
      idmfile = idmatfl
      iemfile = iematfl
      do  ispin=1,nspin
        if( nspin.eq.2 .and. nspink.eq.1 )then
c         Need to re-use spinless base vecs for spin-dn calc
          if( iproc .eq. masterp ) REWIND( unit=ivecfl )
        endif
c
        if( nk .ge. nproc )then
c         Then distribute full k-points over processors
c
          call PRHOIJ2( ivecfl, idmfile,iemfile,
     $     norb,nstate,nk,
     $     eigval(1,ispin),eigpop(1,ispin),wksml,npop, wtk,
     $     wmat(iw1), wmat(iw2), wmat(iw3) )
c -->      dmat-(2)s  emat-s     eigvec-2s(n,nst)
c ------> 2*norb*nstate  mat     2*norb*nstate
c
        else
c         Distribute every k-point over all processors
c
c         Set up sizes of stripe-optimized vectors
          call RHOSTRIP( jstripe, norb, nproc )
          iw2 = iw1 + jstripe*norb
          iw3 = MAX( iw1 + 2*norb*npop, iw2 + jstripe*norb )
c
          call PRHOIJ( ivecfl,idmfile,iemfile,
     $     norb,nstate,nk,
     $     eigval(1,ispin),eigpop(1,ispin),wksml,npop, wtk,
     $     wmat(iw1), wmat(iw2), wmat(iw3) )
c -->      dmat-(2)s  emat-s     eigvec-2s
c ------>  2*norb*npop       ->  2*norb*npop
c ------>  stripe     stripe <-  2*norb*npop
c ------>  stripe     mat
c
        endif
c       Switch to spin-dn:
        idmfile = idmatsfl
        iemfile = iematsfl
c
      enddo
      if( itimlog.gt.0 ) call TIMER( 'eigsolvkp/rhoij      ' )
c
c    That's all Folks!
c
      RETURN
      END

c
c
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPMASTERLIST
c
c
      subroutine KPMASTERLIST( nprocs, iproc, master, icomm,
     $     nprocs_k, iproc_k, master_k, icomm_k,
     $     masterlist, imasters )
c
c  Build a list of kpmasters and compute the number of kpmasters.
c
c  This is needed for the PSCHROEDKP k-parallel solver.
c
      INTEGER  nprocs, iproc, master, icomm
      INTEGER  nprocs_k, iproc_k, master_k, icomm_k
      INTEGER  masterlist, imasters
      INTEGER  jproc
      INTEGER  istate
      DIMENSION masterlist(nprocs)
c    
      imasters = 0
      lenmsg = 1
      do  jproc=0,nprocs-1
        if( iproc .eq. jproc )then
          if( iproc_k .eq. 0 )then
            istate = 0
          else
            istate = 1
          endif
          if( iproc .ne. master )then
            call MPSENDI( master, lenmsg, istate, icomm )
          endif
        endif
        if( iproc .eq. master ) then
          if( jproc .ne. master )then
            call MPRECVI( jproc,  lenmsg, istate, icomm )
          endif
          if( istate.eq.0 )then
            imasters = imasters + 1
            masterlist(imasters) = jproc
          endif
        endif
      enddo
c
c Distribute the list of kp-masters to all processes.      
c
      lenmsg = 1
      call MPBCASTI( master, lenmsg, imasters, icomm )
      lenmsg = imasters
      call MPBCASTI( master, lenmsg, masterlist, icomm )
      do jproc=imasters+1,nprocs, 1
        masterlist(jproc) = -1
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
