c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ERLSOLV
c
c
      subroutine ERLSOLV(
     $ idosolv,scut, nspin, ikoopman,ioccdft,iflcdft,
     $ idhamfl,iham0fl,iovlpfl,
     $ ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $ icloseocc, doblas3,
     $ etemp,edegen, elecno, efermi,egap, esumsp,
     $ norb,nstate, wtk, eigval,eigpop,sefp,npop, iwork,ifail, lwork4,
     $ vmat,  tempio )
c---------------------------------------------------------------
c Purpose: solve real (double precision) Schroedinger equation
c
c Written: Dave Raczkowski and Peter A. Schultz
c
c Revision history:
c  24Jul23-PAS/2.69: koopmans+occdft stuff transferred into serial
c  10Jan13-PAS/2.64: add (here, dummy) tempio arg to simplify seq/mpi
c  22Sep08-PAS/2.62: split smat/spin files; 14Nov-perlsolv compatibility
c  08Jun05-RPM/2.59: added ability to do a canonical orthogonalization
c                    to remove linear dependence in the basis set
c  15Feb05-PAS/2.59: add option to impose closed shell occs
c   8Mar02-PAS/2.52: clean unused/const d0; closed-spin bug
c   5Oct01-PAS/2.50: remove backspace, and compress restart
c  20Jul01-PAS/2.49: spin-polarized dft, compact dmat
c  21Jun01-PAS/2.48: replace STOPs
c  17May01-PAS/2.47: extract listing output (eigenspectrum)
c  20Dec00-PAS/2.46: patch for small basis problems
c  20Aug99-PAS/2.38: bugfix in call to dsygvx; full dsygvx
c                    implemented, doblas3 passed
c  10May99-PAS/2.35: eigenfunction grid density, pass npop
c   2Apr99-PAS/2.33: tighten big memory usage
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  11Feb99-PAS/2.31: adapted for archival code
c---------------------------------------------------------------
c
c Notes:
c  The routine dsygvx, unlike dsygv, is not standard lapack.
c  Dave Raczkowski found the routine dsygvx that is used here
c  for solving the eigenproblem.  It is typically 2X faster
c  than solving the full eigenproblem using dsygv. 23Aug99
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
      DIMENSION  wtk(1), elecno(nspin)
      DIMENSION  vmat(norb,norb,*)
      LOGICAL    doblas3
c Output:
      DIMENSION  efermi(nspin),egap(nspin)
      DIMENSION  eigval(nstate,nspin),eigpop(nstate,nspin)
c Scratch:
c     sefp first used to store eigenvalues of S in algorithm
c     to remove linear dependence from basis, and then as
c     occupation factor to set up density matrix construction
      DIMENSION  sefp(nstate)
      DIMENSION  iwork(5*norb),ifail(norb)
c Scratch/parallel:
      DIMENSION  tempio(*)
c
c Local:
      LOGICAL    closeshl
      DIMENSION  npops(2)
      DIMENSION  occlvl(2)
      DATA       occlvl / 2.d0,1.d0 /
      LOGICAL    pos,neg
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
      DATA  elcerr,elconv / 1.d-6,1.d-11 /
c     "occupied" criterion for constructing number of occupied levels
      DATA  occnil / 1.d-9 /
      DATA  idbglvl / 0 /
c
c >>>> EXECUTABLE CODE:
c
      if( idbglvl .gt. 0 ) call FLGETIWR( IWR )
      if( idbglvl .gt. 0 )then
        write(IWR,'(a,3i8)') 'DEV/ERLSOLV: ikoopman,ioccdft=',
     $     ikoopman, ioccdft, iflcdft
      endif
c
      mat = norb*norb
c
      REWIND( unit=ivecfl )
c
      if( ikoopman .gt. 0 )then
c       Skip eigensolve, will be using the koopmans reference eigvecs
        goto  200
      endif
c
      do 100 ispin=1,nspin
c
        call RSCHROED( idosolv,scut,ispin,nspin,
     $   idhamfl,iham0fl,iovlpfl, ivecfl,
     $   norb,nstate,eigval(1,ispin),
     $   sefp,iwork,ifail,lwork4,vmat)
c
  100 continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   Construct level populations, i.e. dmat in diagonal representation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  200 continue
crpm Misc Notes:
c     vmat(,,1)   Eigenvectors
c     vmat(,,2)   Destroyed H
c     vmat(,,3)   S
c     vmat(,,4)   Scratch
c
      npop = 0
      esumsp = zero
      occfull = occlvl(nspin)
c
c If Koopmans: skipped solver, get eigendata from reference (eigs/vecs)
c   Get eigenvalues from koopman baseline
c   Get spink of koopmans baseline, to know how convert vecs to this job
c
      nk = 1
      nspink = nspin
c
      if( ikoopman .gt. 0 )then
c       With koopmans on, skipped solvers, therefore retrieve old eigenvalues
c       Not strictly true, it worked anyway, but logic is cleaner
          call FLOPENA( ieofile, 'eigkoop' )
          call RDEO( ieofile, 'eigval', 2, nspink, nspin,nk,nstate,
     $               eigval )
          call FLCLOSE( ieofile )
C         write(IWR,*) 'DEV/erlsolv - read koop eigvals, nspink=',nspink
      endif
c
      if( ioccdft .gt. 0 )then
c       Use occdft constrained (preset) occupations from file
c         Read pre-set occupations from file ...
          REWIND( unit=iflcdft )
          call RDEO( iflcdft, 'occups', ioccdft,
     $     mspin,nspin, nk, nstate,
     $     eigpop )
      endif
c
c  ... end koopmans/occdft setup
c
c Set closed shell flag if to populate using 0K
      closeshl = .false.
      if( icloseocc .eq. 2 ) closeshl = .true.
c
      do 300 ispin=1,nspin
        elnumbr = elecno(ispin)
c
        if( ioccdft.gt.0 )then
c         Use dictated occupations
          call POPOCCDFT(  nstate,nk,ispin,
     $     elnumbr, npopsp, occfull,
     $     efermi(ispin), egap(ispin), esum1e,
     $     wtk, eigval(1,ispin), eigpop(1,ispin) )
        else
c         Determine occupations from current solve
          call POP1K( closeshl,
     $     fermilvl,gap,ehomo,elumo,
     $     etemp,edegen, elnumbr,occfull,
     $     nstate,npopsp, wtk(1),
     $     eigval(1,ispin),eigpop(1,ispin), sefp )
c -->                                       fpop
c
          efermi(ispin) = fermilvl
          egap(ispin) = gap
          esum1e = zero
          do  i=1,npopsp
            esum1e = esum1e + eigpop(i,ispin)*eigval(i,ispin)
          enddo
        endif
        npops(ispin) = npopsp
        npop = MAX( npop, npopsp )
c
        esumsp = esumsp + wtk(1)*esum1e
  300 continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       Construct density-matrix, e-matrix over basis functions
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( nspin .eq. 2 .or. ikoopman .gt. 0 )then
c       Then vectors are NOT in mem, must read
        REWIND( unit=ivecfl )
      endif
      REWIND( unit=idmatfl )
      REWIND( unit=iematfl )
      if( nspin .eq. 2 )then
c       Rewind the spin-dn mats files
        REWIND( unit=idmatsfl )
        REWIND( unit=iematsfl )
      endif
      idmfile = idmatfl
      iemfile = iematfl
c
      do 400 ispin=1,nspin
c
        if( nspin.eq.2 .and. nspink.eq.1 )then
c         Need to re-use spinless base vecs for spin-dn calc
          REWIND( unit=ivecfl )
        endif
c
        if( nspin .eq. 2 .or. ikoopman .gt. 0 )then
c         Retrieve spin-dependent eigenstates (in mem if non-spin)
          call READBIG( ivecfl, norb*nstate, vmat(1,1,1) )
        endif
c
        call RRHOIJ( doblas3,
     $       idmfile,iemfile,
     $       norb,nstate,
     $       eigval(1,ispin),eigpop(1,ispin),sefp,
     $       npops(ispin),wtk(1),
     $       vmat(1,1,3),vmat(1,1,4),vmat(1,1,1),vmat(1,1,2))
c       Switch to dn-spin:
        idmfile = idmatsfl
        iemfile = iematsfl
  400 continue
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
