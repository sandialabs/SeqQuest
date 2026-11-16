c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PEIGSOLV
c
c
      subroutine PEIGSOLV( masterlist,do_kppsolve, do_psolv,
     $ icluster, gapp,
     $ do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $ ikoopman, ioccdft,iflcdft,
     $ idhamfl,iham0fl,iovlpfl,
     $ ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $ icloseocc, dokdefct,doblas3, ikdefct,nstbulk,
     $ etemp,edegen, elecno, efermi,egap, esumsp,
     $ norb,nstate, nk, veck,wtk,
     $ eigval,eigpop,  npop, wksml,
     $ wmat )
c---------------------------------------------------------------
c Purpose: branch to complex solvers for density matrix
c
c Written: Peter A. Schultz, 4-January-2013, for 2.64 (kpp),
c          to accommodate Andrew C. Pineda kpp-solver codes
c
c Revision history:
c   1Dec22-PAS/2.69: koopmans
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
      DIMENSION  masterlist(*)
      LOGICAL    do_kppsolve, do_psolv
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
C      DIMENSION  wk(*)
c Scratch/parallel:
      INTEGER    icluster(*)
      DIMENSION  gapp(*)
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
c
      if( IWR .ge. 0 )then
        write(IWR,*) 'PEIGSOLV: '//
     $               'do_eigvecs,do_eigpops,do_kppsolve,do_psolv=',
     $                do_eigvecs,do_eigpops,do_kppsolve,do_psolv
        write(IWR,*) 'DEV/PEIGSOLV: smaller PRHOIJ stripe/transpose'
      endif
c
c  Solve Schroedinger equation for the various k-vectors
      call MPNODES( nprocs )
c
      if( .not. do_kppsolve )then
c       K-parallel serial solves OR K-serial all-processor-parallel solves
c
        call PEIGSOLV0( do_psolv,
     $   icluster, gapp,
     $   do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $   ikoopman,ioccdft, iflcdft,
     $   idhamfl,iham0fl,iovlpfl,
     $   ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $   icloseocc, dokdefct, doblas3, ikdefct,nstbulk,
     $   etemp,edegen, elecno, efermi,egap, esumsp,
     $   norb,nstate, nk, veck,wtk,
     $   eigval,eigpop,  npop, wksml,
     $   wmat )
c -->    vmat-6s+
c
      else
c       k-parallel-parallel solves
c
c       Invokes ScaLAPACK solvers within k-point groups.
c       Should enable code to scale to 16 processors/k-point or better
c
        call PEIGSOLVKP( masterlist,
     $   icluster, gapp,
     $   do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $   ikoopman,ioccdft, iflcdft,
     $   idhamfl,iham0fl,iovlpfl,
     $   ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $   icloseocc, dokdefct, doblas3, ikdefct,nstbulk,
     $   etemp,edegen, elecno, efermi,egap, esumsp,
     $   norb,nstate, nk, veck,wtk,
     $   eigval,eigpop, npop, wksml,
     $   wmat )
c -->    vmat-6s+
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PEIGSOLV0
c
c
      subroutine PEIGSOLV0( do_psolv,
     $ icluster, gapp,
     $ do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $ ikoopman,ioccdft, iflcdft,
     $ idhamfl,iham0fl,iovlpfl,
     $ ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $ icloseocc, dokdefct,doblas3, ikdefct,nstbulk,
     $ etemp,edegen, elecno, efermi,egap, esumsp,
     $ norb,nstate, nk, veck,wtk,
     $ eigval,eigpop, npop, wksml,
     $ wmat )
c---------------------------------------------------------------
c Purpose: calling routine to solve for density matrix, given
c          hamiltonian and overlap matrices of current iteration
c
c Written: Peter A. Schultz
c
c Revision history:
c   1Dec22-PAS/2.69: koopmans
c   3Jan13-PAS/2.64: cosmetics, un-drivered, and integrated
c   2May12-ACP/2.64: Added alternate entry point to support non-SCF portions of code
c   9Dec08-PAS/2.62: reduced-mem parallel-complex solve; pschroed w/nk<4
c   1May08-PAS/2.61b: patch bug - purged POP argument (edegen)
c  25Jun07-PAS/2.60: merge tp and serial
c   7Jul06-PAS/2.60: k-parallel dmat generation added
c  14Aug05-PAS/2.59: master-only matrix files, read multi-k Hams
c  15Feb05-PAS/2.59: options to impose 0K occs, and defect sampling
c  25Sep04-PAS/tp0.1: clean-up task parallel routines
c  30Apr03-DMC     : Added a parallel version of RHOIJ
c  22Apr03-DMC     : Added a second version of a parallel SCHROED to take
c                    advantage when nk > 1.
c   7Apr03-DMC     : Made parallel by modifying SCHROED to be parallel.
c                    Does not appear to be sufficient work in POP or RHOIJ
c                    to warrant changing at this time.
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
c   (1) eigenvectors are saved to "ivecfl"
c   (2) density matrix is saved to "idmat"
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
      LOGICAL    do_psolv
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
      CHARACTER*6  labelxxx
      DATA       itimlog / 5 /
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
      call FLGETIWR( IWR )
c  Declare dev-version of peigsolv0:
 9001 format(a,6i8)
      if( IWR .ge. 0 )then
        write(IWR,*) 'PEIGSOLV0: do_psolv,do_eigvecs,do_eigpops=',
     $                           do_psolv,do_eigvecs,do_eigpops
      endif
c
c Partition space within wmat() for schroed:
      iw1 = 1
      iww = iw1 + mat
      nwsch = norb*(4*norb + 2*nstate)
      iwx = iww + nwsch
c
      call MPCOMM( icomm )
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      if( iproc .eq. masterp )then
c       Go to start of idhamfl, to get ready for dHam reads
        REWIND( unit=idhamfl )
c       Go to start of ivecfl, to get ready for eigvec writes
        REWIND( unit=ivecfl )
      endif
c
      if( itimlog.gt.0 ) call TIMER( 'peigsolv/begin     ' )
c
C     if( .not. do_psolv )then
c
c Skip eigsolvers with koopmans on, go to populations
      if( ikoopman .gt. 0 )then
c       Koopmans is on, skip eigensolve for new vectors
        goto 2000
      endif
      if( nprocs.eq.1 .or.
     $    (  .not. do_psolv .and.
     $      ( norb.lt.800 .or. nk.gt.4 .or. nprocs.lt.(4*nk) )  )      
     $                      )then
c
c       Record path through solver: single-proc solves
        if( IWR .gt. 0 ) write(IWR,*) 'PEIGSOLV: single-proc solves'
c
        do  ispin=1,nspin
c
c         Parallel over k-points - single-proc solves
          call PSCHROED2( do_eigvecs,
     $     idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $     norb,nstate,nk, eigval(1,ispin), wksml,    wksml(1,9), scut,
c -->                                       rwork-7ns seval(n)-s
     $     wmat(iw1), wmat(iww), wmat(iwx) )
c -->      vmat-s     wmat-6ks   iwork-6ns
c
          if( itimlog.gt.1 ) call TIMER( 'peigsolv/schroedkp ' )
        enddo
c
      else
c       All-proc parallel solvers
c
c       Record path through solver: parallel
        if( IWR.ge.0 ) write(IWR,*) 'PEIGSOLV: all-parallel solves'
c
        call PSCHROED( do_eigvecs,
     $   idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $   norb,nstate,nk, eigval, wksml,
c -->                            wksml-ns
     $   nspin, icluster, gapp,
     $   wmat )
c -->    wmat-s
c
        if( itimlog.gt.1 ) call TIMER( 'peigsolv/pschroed ' )
c
      endif
c
      if( .not. do_eigpops )then
c       Skip populations and density matrix construction
        RETURN
      endif
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
        if( iproc .eq. masterp )then
c         Master reads this from a file
          call FLOPENA( ieofile, 'eigkoop' )
          call RDEO( ieofile, 'eigval', 2, nspink, nspin,nk,nstate,
     $               eigval )
          call FLCLOSE( ieofile )
        endif
c       Tell all the other processors
        call MPI_Bcast( eigval(1,1), nstate*nk*nspin,
     $   MPI_DOUBLE_PRECISION,
     $   masterp, icomm, ierr )
        if( ierr.ne.0 ) call STOPXERR( 'peigsolv - Bcast eigval' )
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

      do 3000 ispin=1,nspin
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
c         Use standard metallic (across k-points) occupations
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
 3000 continue
c
      if( itimlog.gt.0 ) call TIMER( 'peigsolv/pop       ' )
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
      call MPNODES( nprocs )
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
        if( nk .ge. nprocs )then
c         K-parallel, single-processor dmat(k) construction
c
C          write(IWR,*) "DEVkoop/peigsolv0: before k-par prhoij2"
          call PRHOIJ2( ivecfl, idmfile,iemfile,
     $     norb,nstate,nk,
     $     eigval(1,ispin),eigpop(1,ispin),wksml,npop, wtk,
     $     wmat(iw1), wmat(iw2), wmat(iw3) )
c -->      dmat-(2)s  emat-s     eigvec-2s(n,nst)
c
        else
c         All-processor parallel: distribute each k-point
c
c         Set up sizes of stripe-optimized vectors
          call RHOSTRIP( jstripe, norb, nprocs )
          iw2 = iw1 + jstripe*norb
          iw3 = MAX( iw1 + 2*norb*npop, iw2 + jstripe*norb )
c
C          write(IWR,*) "DEVkoop/peigsolv0: stripe-par prhoij"
          call PRHOIJ( ivecfl,idmfile,iemfile,
     $     norb,nstate,nk,
     $     eigval(1,ispin),eigpop(1,ispin),wksml,npop, wtk,
     $     wmat(iw1), wmat(iw2), wmat(iw3) )
c -->      dmat-(2)s  emat-s     eigvec-2s
c ------>  2*norb*npop       ->  2*norb*npop
c ------>  stripe     stripe <-  2*norb*npop
c ------>  stripe     mat
c         Total prhoij mem = max(A,B)
c            A= norb*( max( 2*npop, 2*jstripe ) + 2*npop )
c            B= jstripe*norb + norb^2
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
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
