c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EIGSOLV
c
c
      subroutine EIGSOLV(
     $ icluster, gapp,
     $ do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $ ikoopman, ioccdft,iflcdft,
     $ idhamfl,iham0fl,iovlpfl,
     $ ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $ icloseocc,dokdefct,doblas3, ikdefct,nstbulk,
     $ etemp,edegen, elecno, efermi,egap, esumsp,
     $ norb,nstate, nk, veck,wtk, eigval,eigpop,npop, wksml,
     $ wmat )
c---------------------------------------------------------------
c Purpose: calling routine to solve for density matrix, given
c          hamiltonian and overlap matrices of current iteration
c
c Written: Peter A. Schultz
c
c Revision history:
c   3Nov24-JJL/2.70: maximum overlap method into serial version
c  24Jul23-PAS/2.69: koopmans/occdft into serial version
c   3May12-ACP/2.64: Add vector/pop flags (PAS) to support nonscf
c  22Sep08-PAS/2.62: spin-split dmat/emat
c  12Jan08-PAS/2.61: spin optimization
c  25Jun07-PAS/2.60: merge tp and serial
c  11Sep05-PAS/2.59: remove schroed lin.dep.-based dynamic mem
c  14Aug05-PAS/2.59: read multi-k Hams
c  15Feb05-PAS/2.59: options to impose 0K occs, and defect sampling
c   8Aug03-PAS/2.56: cosmetics: clean up k-vector usage
c   5Oct01-PAS/2.50: remove backspace, and compress restart
c  20Jul01-PAS/2.49: spin-polarized dft, fpop made scratch, dmat
c  17May01-PAS/2.47: extract listing output (eigenspectrum)
c  23Aug99-PAS/2.38: install zhegvx/selected eigvec solver
c  10May99-PAS/2.35: eigenfunction grid density, pass npop
c   9Apr99-PAS/2.34: complex lapack eigensolver installed
c  31Mar99-PAS/2.33: clean out green fcn stuff, tighten mem use
c   9Feb99-PAS/2.31: pass write unit into routine
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  14Jul98-PAS/2.23: pass option for un-chol'd d-Ham into schroed
c  10Oct97-PAS/2.21: change convention for fpop
c  13Sep93-PAS: rhoofr extracted out of routine, calling arguments
c               rationalized, dimensions fixed, renamed "eigsolv"
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
      DIMENSION  gapp(*)
c
c Local:
      DATA       itimlog / 0 /
      DIMENSION  occlvl(2)
      DATA       occlvl / 2.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
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
      if( itimlog.gt.0 ) call TIMER( 'eigsolv/begin      ' )
c
c Skip eigsolvers with koopmans on, go to populations
      if( ikoopman .gt. 0 )then
c       With koopmans set, skip eigensolve
        goto 2000
      endif
c
      if( nproc.eq.1 .or. nk.gt.1 )then
c       Parallel over k-points - single-proc solves
c
        do  ispin=1,nspin
c
          call SCHROED( do_eigvecs,
     $     idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $     norb,nstate,nk, eigval(1,ispin), wksml,    wksml(1,9), scut,
c -->                                       rwork-7ns seval(n)-s
     $     wmat(iw1), wmat(iww), wmat(iwx) )
c -->      vmat-s     wmat-6ks   iwork-6ns
c
          if( itimlog.gt.1 ) call TIMER( 'eigsolv/psch2 k-par' )
        enddo
c
      else
c       Parallel solvers
c
c       Not available here:
        call STOPXERR( 'Parallel solver not enabled' )
c
      endif
c
      if( .not. do_eigpops )then
c       Skip state populations and density-matrix construction
        RETURN
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  Evaluate density matrix in diagonal (eigfcn) representation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 2000 continue
c
      nstk = nstate*nk
c
      npop = 0
      esumsp = zero
      occfull = occlvl(nspin)
c
c If Koopmans: skipped solver, get eigen data from reference (eigs/vecs)
c   Get eigenvalues from koopman baseline
c   Get spink of koopmans baseline, to know how convert vecs to this job
c
      nspink = nspin
      if( ikoopman .gt. 0 )then
c       With koopmans on, skipped solvers, therefore retrieve old eigenvalues
c       Not strictly true, it worked anyway, but logic is cleaner
          call FLOPENA( ieofile, 'eigkoop' )
          call RDEO( ieofile, 'eigval', 2, nspink, nspin,nk,nstate,
     $               eigval )
          call FLCLOSE( ieofile )
C         write(IWR,*) 'DEVkoop/eigsolv - read eigvals, nspink=',nspink
      endif
c
      if( ioccdft .gt. 0 )then
c       Use occdft constrained (preset) occupations from file
c         Master reads pre-set occupations from file ...
          REWIND( unit=iflcdft )
          call RDEO( iflcdft, 'occups', ioccdft,
     $     mspin,nspin, nk, nstate,
     $     eigpop )
      endif
c
      do  ispin=1,nspin
        elnumbr = elecno(ispin)
c
        if( ioccdft.gt.0 )then
c         Using pre-set occdft occupations
          call POPOCCDFT(  nstate,nk,ispin,
     $     elnumbr, npopsp, occfull,
     $     efermi(ispin), egap(ispin), esum1e,
     $     wtk, eigval(1,ispin), eigpop(1,ispin) )

        elseif( ( icloseocc .gt. 0 ) .or. dokdefct )then
c         Use closed-shell occupations, or defect sampling
          call POPCLOSE( icloseocc,dokdefct,ikdefct,
     $     efermi(ispin),egap(ispin), nstbulk(ispin),
     $     etemp,edegen, elnumbr,esum1e, occfull,
     $     nstate,nk, wtk, eigval(1,ispin),eigpop(1,ispin),npopsp,
     $     wmat )
c -->      wlvl-s(nstate)
c
        else
c         Use standard metallic occupations, spin
          mspin = 1
          call POP( mspin, popnumbr,
     $     efermi(ispin),egap(ispin),
     $     etemp,        elnumbr,esum1e, occfull,
     $     nstate,nk, wtk, eigval(1,ispin),eigpop(1,ispin),npopsp,
     $     wmat(1), wmat(1+nstk),wmat(1+2*nstk),
     $     wmat(1+3*nstk), wmat(1+4*nstk) )
c -->      elvl wlvl klvl nlvl flvl - all scratch (nk*state)
        endif
c
        npop = MAX( npop,npopsp )
        esumsp = esumsp + esum1e
      enddo
c
      if( itimlog.gt.0 ) call TIMER( 'eigsolv/pop        ' )
c
      if( .not. do_eigvecs )then
c       Do not compute density matrix if do not have vectors ...
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
          REWIND( unit=ivecfl )
        endif
c
        if( nk .ge. nproc )then
c         Then distribute full k-points over processors
c
          call RHOIJ( ivecfl, idmfile,iemfile,
     $     norb,nstate,nk,
     $     eigval(1,ispin),eigpop(1,ispin),wksml,npop, wtk,
     $     wmat(iw1), wmat(iw2), wmat(iw3) )
c -->      dmat-(2)s  emat-s     eigvec-2s(n,nst)
c
        else
c         Distribute every k-point over all processors
c
c         Not available in serial
          call STOPXERR( 'Parallel dmat not enabled' )
c 
        endif
c       Switch to spin-dn:
        idmfile = idmatsfl
        iemfile = iematsfl
c
      enddo
      if( itimlog.gt.0 ) call TIMER( 'eigsolv/rhoij      ' )
c
c    That's all Folks!
c
      RETURN
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000
characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
