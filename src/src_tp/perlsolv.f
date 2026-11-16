c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PERLSOLV
c
c
      subroutine PERLSOLV( icluster, gapp, nlr,nlc,
     $ idosolv,scut, nspin,ikoopman,ioccdft,iflcdft,
     $ idhamfl,iham0fl,iovlpfl,
     $ ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $ icloseocc, doblas3,
     $ etemp,edegen, elecno, efermi,egap, esumsp,
     $ norb,nstate, wtk, eigval,eigpop,fpop,npop, iwork,ifail, lwork8,
     $ vmat,  tempio )
c---------------------------------------------------------------
c Purpose: solve real (double precision) Schroedinger equation
c
c Written: Dave Raczkowski, Peter A. Schultz, and Dean P. McCullough
c
c Revision history:
c  16Jan23-PAS/2.69: koopmans+occdft stuff
c   3Jan13-PAS/2.64: cosmetics, un-drivered, and integrated
c   2May12-ACP/2.64: Added flags for non-scf, and distribute >2GB matrix
c  20Dec07-PAS/2.60d: close BLACS INIT with EXIT in PSOLVSIZ
c  30Aug07-RPM/2.60: infrostructure to reduce parallel memory
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  30Jul05-PAS/2.59: master-only matrix files, pass scut
c  15Feb05-PAS/2.59: add option to impose closed shell occs
c  27Oct04-PAS/tp0.1: clean-up/debug task-parallel routines
c   5Jul03-DMC: Need eigsolv to be norb*nk.  However in wkmem, was only
c               setting up nstate*nk.
c  10Feb03-DMC: Modified extensively to use ScaLAPACK and PBLAS routines.
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
c  Call parameter idosolv was used to control which version of the
c  LAPACK routine DSYGV ws used.  Here, we will always use PDSYGVX.
c  The parameter will be used only to control whether all, or only
c  nstate eigenvalues/eigenvectors will be requested.
c  DMC 12/Feb/03
c
c  The size of iwork in the calling routine needs to be increased to
c  at least 6*norb (vs. 5*norb).
c
c  The doblas3 flag is ignored in solvers.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c     Upper limit on size of array that PDGEMR2D can distribute.
c     MAXBLOCKR= 2**31/sizeof(DOUBLE PRECISION) = 2**31/8.
      PARAMETER  ( maxblockr = 268 435 456 )
c
ctp: Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Input:
      DIMENSION  wtk(1), elecno(nspin)
      DIMENSION  vmat(nlr,nlc,*)
      LOGICAL    doblas3
c Output:
      DIMENSION  efermi(nspin),egap(nspin)
      DIMENSION  eigval(nstate,nspin),eigpop(nstate,nspin)
c Scratch:
      DIMENSION  fpop(norb)
      DIMENSION  iwork(6*norb),ifail(norb)
c Scratch/parallel:
      INTEGER    icluster(*)
      DOUBLE PRECISION  gapp(*)
      DIMENSION  tempio(norb*norb)
c
ctp: parallel local
c Local/parallel
      PARAMETER  ( NB = 32 )
      INTEGER    ICONTXT
      INTEGER    DESCA(9), DESCWORK(9)
      DATA       itimlog / 3 /
      LOGICAL    do_eigvecs, do_eigpops
      DATA       do_eigvecs, do_eigpops / .true., .true. /
      LOGICAL    DEBUG_SOLV
      DATA       DEBUG_SOLV / .true. /
c Local:
      LOGICAL    closeshl
      DIMENSION  npops(2)
      DIMENSION  occlvl(2)
      DATA       occlvl / 2.d0,1.d0 /
      CHARACTER  jobz*1
      LOGICAL    pos,neg
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
      DATA  elcerr,elconv / 1.d-6,1.d-11 /
c     for constructing number of occupied levels
      DATA  occnil / 1.d-9 /
c
c >>>> EXECUTABLE CODE:
c
      if( itimlog.gt.0 ) call TIMER( 'perlsolv/begin ' )
      call FLGETIWR( IWR )
cxxx: declare dev-version of perlsolv
      if( IWR .ge. 0 )then
        write(IWR,*) 'PERLSOLV: kpp-compatible version'
      endif
c
      if( do_eigvecs )then
c       Compute eigenvalues and eigenvectors
        jobz = 'V'
      else
c       Compute eigenvalues only
        jobz = 'N'
      endif
c
      mat = norb*norb
c
c  Get Local MP info
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
c
c Initialize ScaLAPACK Grid
c
c     Create a close to square grid of processors for solver
      rnprocs = DBLE( nprocs ) + 0.1
      npcol = SQRT( rnprocs )
      nprow = nprocs / npcol
      npsolver = nprow*npcol
c
c  Create BLACS context using system context, 
c  i.e. the local MP communicator
c
      ICONTXT = icomm
      call BLACS_GRIDINIT( ICONTXT, 'R', nprow, npcol )
      call BLACS_GRIDINFO( ICONTXT, nprow, npcol, myrow, mycol )
c
      if( myrow .ne. -1 )then
c       Initialize matrix descriptions, for procs in solver context
c
c       Set up the descriptor for a distributed matrix
        call DESCINIT( DESCA, norb,norb, NB,NB, 0,0,
     $   ICONTXT, nlr, INFO )
        if( info.ne.0 ) call STOPXERR( 'perlsolv: INIT for DESCA' )
c        nlocr = NUMROC(norb,NB,myrow,0,nprow)
c        nlocc = NUMROC(norb,NB,mycol,0,npcol)
c       This call to PSOLVSIZES is not uesed (yet).
        call PSOLVSIZES( ICONTXT, norb,nstate,nlocr,nlocc,mylwork )
c
c       Set up the descriptor for full local matrix
        call DESCINIT( DESCWORK, norb,norb, norb,norb, 0,0,
     $   ICONTXT, norb, INFO )
        if( info.ne.0 ) call STOPXERR( 'perlsolv: INIT for DESCWORK' )
      endif
c
c Set files to correct locations
c
      if( iproc .eq. masterp )then
        REWIND( unit=ivecfl )
        REWIND( unit=idhamfl )
      endif
c
      if( ikoopman.gt.0 )then
c       Skip solvers with koopmans
        goto 1500
      endif
c
      do 1000 ispin=1,nspin
        if( myrow .eq. -1 ) goto 800
c
c       Read in array to work space, and  distribute to all processors.
c       If array is too big to fit on one processor, this can be done in
c       several steps.
c
c       If nspin = 1, this read not done in the original serial code.
c       However, I (DMC) do not  know where it comes from, and the fill
c       might not be distributed across the processors, thus reload it.
c       Found I also need to rewind it.
c
c       Retrieve all of delta-Ham onto master ...
        if( iproc .eq. masterp )then
          call READBIG( idhamfl, mat, tempio )
        endif
c        ... and distribute d-Ham onto all processors.
        if( mat .le. maxblockr )then
c          Matrix small enough to distribute from full local storage
           call PDGEMR2D( norb,norb,
     $          tempio, 1,1, DESCWORK,
     $          vmat(1,1,1), 1,1, DESCA,  ICONTXT )
        else
c         Distribute the matrix in (column) stripes that are smaller than 2GB each.
c         Should be good as long the distributed matrix is no larger than 2GB on any processor.
          nblock = maxblockr / norb
          nleft = norb - nblock
          noffset = 1
          moffset = 1
  151     continue
          if( noffset.le.norb )then
            call PDGEMR2D( norb,nblock,
     $           tempio(moffset), 1,1, DESCWORK,
     $           vmat(1,1,1), 1, noffset, DESCA,  ICONTXT )
            noffset = noffset + nblock
            moffset = moffset + norb*nblock
            if( nleft.ge.nblock ) then
c             Have enough to do another full block.
              nleft = nleft - nblock
            else
c             Have less than a full block to send.
              nblock = nleft
              nleft = 0
            endif
            goto 151
          endif
        endif
c
c       Retrieve iter-independent Hamiltonian ...
        if( iproc .eq. masterp )then
          REWIND( unit=iham0fl )
          call READBIG( iham0fl, mat, tempio )
        endif
c        ... and distribute over processors:
        if( mat .le. maxblockr )then
          call PDGEMR2D( norb,norb,
     $         tempio, 1,1, DESCWORK,
     $         vmat(1,1,2), 1,1, DESCA,  ICONTXT )
        else
c         Distribute the matrix in (column) stripes that are smaller than 2GB each.
c         Should be good as long the distributed matrix is no larger than 2GB on any processor.
          nblock = maxblockr / norb
          nleft = norb - nblock
          noffset = 1
          moffset = 1
  152     continue
          if( noffset.le.norb )then
            call PDGEMR2D( norb,nblock,
     $           tempio(moffset), 1,1, DESCWORK,
     $           vmat(1,1,2), 1, noffset, DESCA,  ICONTXT )
            noffset = noffset + nblock
            moffset = moffset + norb*nblock
            if( nleft.ge.nblock )then
c             Have enough to do another full block.
              nleft = nleft - nblock
            else
c             Have less than a full block to send.
              nblock = nleft
              nleft = 0
            endif
            goto 152
          endif
        endif
c
c       Retrieve overlap matrix ...
        if( iproc .eq. masterp )then
          REWIND( unit=iovlpfl )
          call READBIG( iovlpfl , mat, tempio )
        endif
c        ... and distribute over processors:
        if( mat .le. maxblockr )then
           call PDGEMR2D( norb,norb,
     $          tempio, 1,1, DESCWORK,
     $          vmat(1,1,3), 1,1, DESCA,  ICONTXT )
        else
c         Distribute the matrix in (column) stripes that are smaller than 2GB each.
c         Should be good as long the distributed matrix is no larger than 2GB on any processor.
          nblock = maxblockr / norb
          nleft = norb - nblock
          noffset = 1
          moffset = 1
  153     continue
          if( noffset.le.norb )then
            call PDGEMR2D( norb,nblock,
     $           tempio(moffset), 1,1, DESCWORK,
     $           vmat(1,1,3), 1, noffset, DESCA,  ICONTXT )
            noffset = noffset + nblock
            moffset = moffset + norb*nblock
            if( nleft.ge.nblock )then
c             Have enough to do another full block.
              nleft = nleft - nblock
            else
c             Have less than a full block to send.
              nblock = nleft
              nleft = 0
            endif
            goto 153
          endif
        endif
c
c       Combine iter-independent Ham and del-Ham to get total Ham
c
        do  i=1,norb
c          do  ii=1,norb
c            vmat(ii,i,2) = vmat(ii,i,2) + vmat(ii,i,1)
c          enddo
          call PDAXPY( norb, one, vmat(1,1,1), 1,i, DESCA, 1,
     $                            vmat(1,1,2), 1,i, DESCA, 1 )
        enddo
c
crpm  Do we really need this now that we're using PDAXPY? Guess it can't hurt
c        call BLACS_BARRIER( ICONTXT, 'All' )
c
c  Total hamiltonian is now in vmat(,,2), overlap in vmat(,,3)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                    Solve Schroedinger Equation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c lwork is workspace available within vmat(,,4), and normally you
c would set lwork=mat(dist) and that would be more than ample.
c But for small problems (e.g. H-atom), we crash as the routine wants
c more than the small mat it gets.  AAAArrrrrgggghhhh!
c Hence I pass in the full free space available ...
c
        lwork = lwork8
c
        if( idosolv .eq.1 )then
          jstate = norb
        elseif( idosolv .eq. 2 )then
          jstate = nstate
        else
          call STOPXERR( 'perlsolv: invalid solver option idosolv' )
        endif
c
        rmone = 2*PDLAMCH( ICONTXT, 'S' )
c
        call PDSYGVX( 1,jobz, 'I','L',
     $   norb,  vmat(1,1,2),1,1,DESCA,  vmat(1,1,3),1,1,DESCA,
     $   zero,zero,1,jstate,  rmone, nstout, nvout,
     $   fpop           ,rmone,
     $   vmat(1,1,1),1,1,DESCA,  vmat(1,1,4),
c -->    Z                       WORK
     $   lwork,iwork,6*norb,ifail,icluster,gapp,info )
c
        if( info .ne. 0 )then
          call STOPXERR( 'perlsolv: pdsygvx error' )
        endif
        call DCOPY( nstate, fpop,1, eigval(1,ispin),1 )
c
c       If no eigvecs, we are done with this spin-Ham:
        if( .not. do_eigvecs ) goto 1000
c
c       Save the the eigenvectors in vmat(,,1) to ivecfl
c
c       Bring back eigenvectors to master processor
        if( mat .le. maxblockr )then
          call PDGEMR2D( norb,norb,
     $         vmat(1,1,1), 1,1, DESCA,
     $         tempio, 1,1, DESCWORK,  ICONTXT )
         else
c          Distribute the matrix in (column) stripes that are smaller than 2GB each.
c          Should be good as long the distributed matrix is no larger than 2GB on any processor.
           nblock = maxblockr / norb
           nleft = norb - nblock
           noffset = 1
           moffset = 1
  154      continue
           if( noffset.le.norb )then
             call PDGEMR2D( norb,nblock,
     $            vmat(1,1,1), 1,1, DESCA, 
     $            tempio(moffset), 1, noffset, DESCWORK,  ICONTXT )
             noffset = noffset + nblock
             moffset = moffset + norb*nblock
             if( nleft.ge.nblock )then
c              Have enough to do another full block.
               nleft = nleft - nblock
             else
c              Have less than a full block to send.
               nblock = nleft
               nleft = 0
             endif
             goto 154
           endif
        endif
c
c       Processors not in ScaLAPACK solver context return here:
  800   continue
Cc
Cc       Record eigenvectors and values
Cc
Cc       Broadcast eigenvectors from master to all processors
C        call MPI_Bcast( vmat(1,1,4), norb*nstate, MPI_DOUBLE_PRECISION,
C     $   0, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'perlsolv: Bcast eigvec' )
Cc
c       Save eigenvectors to disk for later
        if( do_eigvecs .and. ( iproc .eq. masterp ) )then
          call WRITBIG( ivecfl, norb*nstate, tempio )
        endif
c
c       Close loop over spins:
 1000 continue
c
c  If some processors not used in the ScaLAPACK context ...
      if( nprocs .ne. npsolver )then
c        ... also need to broadcast eigvals to those not having them
        call MPI_Bcast( eigval(1,1),nstate*nspin, MPI_DOUBLE_PRECISION,
     $   0, icomm, ierr )
        if( ierr.ne.0 ) call STOPXERR( 'perlsolv: Bcast eigval' )
      endif
ctp:
      if( itimlog.gt.1 ) call TIMER( 'erlsolv/schroed' )
c
      if( .not. do_eigpops )then
c       Skip populations and density matrix
        RETURN
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   Construct level populations, i.e. dmat in diagonal representation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
 1500 continue
c
      npop = 0
      esumsp = zero
      occfull = occlvl(nspin)
c
cdev/ckoop: get koopmans eigvals
c If Koopmans: skipped solver, get eigendata from reference (eigs/vecs)
c   Get eigenvalues from koopman baseline
c   Get spink of koopmans baseline, to know how convert vecs to this job
c
      nk = 1
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
c  ... end koopmans/occdft setup
c
c     Set closed shell flag if to populate using 0K
      closeshl = .false.
      if( icloseocc .eq. 2 ) closeshl = .true.
c
      do 2000 ispin=1,nspin
        elnumbr = elecno(ispin)
c
        if( ioccdft.gt.0 )then
          call POPOCCDFT(  nstate,nk,ispin,
     $     elnumbr, npopsp, occfull,
     $     efermi(ispin), egap(ispin), esum1e,
     $     wtk, eigval(1,ispin), eigpop(1,ispin) )
c
        else
          call POP1K( closeshl,
     $     fermilvl,gap,ehomo,elumo,
     $     etemp,edegen, elnumbr,occfull,
     $     nstate,npopsp, wtk(1),
     $     eigval(1,ispin),eigpop(1,ispin), fpop )
c
          efermi(ispin) = fermilvl
          egap(ispin) = gap
          esum1e = zero
          do  i=1,npopsp
            esum1e = esum1e + eigpop(i,ispin)*eigval(i,ispin)
          enddo
        endif
c
        npops(ispin) = npopsp
        npop = MAX( npop, npopsp )
c
        esumsp = esumsp + esum1e
 2000 continue
c
      if( itimlog.gt.1 ) call TIMER( 'erlsolv/pop    ' )
c
      if( .not. do_eigvecs )then
c       Without eigvecs, skip density matrix construction
        RETURN
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       Construct density-matrix, e-matrix over basis functions
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( nspin .eq. 2 )then
c       Need to recover eigvecs ... they have been stepped on
        if( iproc .eq. masterp )then
          REWIND( unit=ivecfl )
        endif
      endif
      if( iproc .eq. masterp )then
        REWIND( unit=idmatfl )
        REWIND( unit=iematfl )
        if( nspin .eq. 2 )then
c         Rewind the spin-dn mats files
          REWIND( unit=idmatsfl )
          REWIND( unit=iematsfl )
        endif
      endif
      idmfile = idmatfl
      iemfile = iematfl
c
      do 4000 ispin=1,nspin
c       Skip if this processor is not involved
        if( myrow .eq. -1 ) goto 2995
        if( nspin.eq.2 .and. nspink.eq.1 )then
c         Need to re-use spinless base vecs for spin-dn calc
          if( iproc .eq. masterp ) REWIND( unit=ivecfl )
        endif
c
        numpop = npops(ispin)
        if( nspin .eq. 2 )then
c         Retrieve spin-dependent eigenvectors (in mem if non-spin)
          if( iproc .eq. masterp )then
            call READBIG( ivecfl, norb*nstate, tempio )
          endif
c          ... and distribute over processors
          if( mat .le. maxblockr )then
c           Matrix is small enough to be distributed monolithically via BLACS ...
            call PDGEMR2D( norb,norb,
     $           tempio, 1,1, DESCWORK,
     $           vmat(1,1,1), 1,1, DESCA,  ICONTXT )
          else
c           Matrix is larger than BLACS can handle in one chunk ...
c           Distribute matrix in (column) stripes smaller than 2GB each.
c           Should be good as long the distributed matrix <2GB on any processor.
            nblock = maxblockr / norb
            nleft = norb - nblock
            noffset = 1
            moffset = 1
 2155       continue
            if( noffset.le.norb )then
              call PDGEMR2D( norb,nblock,
     $             tempio(moffset), 1,1, DESCWORK,
     $             vmat(1,1,1), 1, noffset, DESCA,  ICONTXT )
              noffset = noffset + nblock
              moffset = moffset + norb*nblock
              if( nleft.ge.nblock )then
c               Have enough to do another full block.
                nleft = nleft - nblock
              else
c               Have less than a full block to send.
                nblock = nleft
                nleft = 0
              endif
              goto 2155
            endif
          endif
        endif
c
        if( numpop .lt. 1 )then
c         No occupied levels: dmat=emat=0, and skip construction
crpm      Testing whether we can use PDLASET instead of MKZERO here:
c          call MKZERO( mat, vmat(1,1,3) )
c          call MKZERO( mat, vmat(1,1,4) )
          call PDLASET( 'A',norb,norb,zero,zero,vmat(1,1,3),1,1,DESCA )
          call PDLASET( 'A',norb,norb,zero,zero,vmat(1,1,4),1,1,DESCA )
          goto 2990
        endif
c
c       Scale level occupations with appropriate sample weights
        wt = wtk(1)
        do  n=1,numpop
          fpop(n) = wt*eigpop(n,ispin)
        enddo
c
c >>>>> Closed shell, use Raczkowski blas code:
c
c       Make a copy of the eigenvectors ...
        do  i=1,norb
           call PDCOPY( norb,vmat(1,1,1), 1,i, DESCA, 1,
     $                       vmat(1,1,2), 1,i, DESCA, 1 )
        enddo
c
c        ni = (numpop-1) / (NB*npcol)
c        ... and scale copy vmat(:,i,2) by occupations fpop(i)
c        do  i=0,ni
c          ki = i*NB
c          kii = (i*npcol + mycol)*NB
c          nj = MIN( NB, numpop-kii )
c          if( nj .ge. 0 )then
c            do  j=1,nj
c              do  jj=1,norb
c                vmat(jj,ki+j,2) = vmat(jj,ki+j,2)*fpop(kii+j)
c              enddo
c            enddo
c          endif
c        enddo
crpm    Replace the explicit math with a call to PDSCAL:
        do  i=1,numpop
          call PDSCAL( norb, fpop(i),vmat(1,1,2),1,i,DESCA, 1 )
        enddo
c
c       To make dscal happy (not return a floating point error)
crpm    Testing whether PDLASET will work here:
c        call MKZERO( mat, vmat(1,1,3) )
c        call MKZERO( mat, vmat(1,1,4) )
        call PDLASET( 'A',norb,norb,zero,zero,vmat(1,1,3),1,1,DESCA )
        call PDLASET( 'A',norb,norb,zero,zero,vmat(1,1,4),1,1,DESCA )
c
c       Make density matrix:
        call PDGEMM( 'n','t', norb,norb,numpop,one,
     $   vmat(1,1,2),1,1,DESCA, vmat(1,1,1),1,1,DESCA,
     $   zero, vmat(1,1,3),1,1,DESCA )
c
c       Make e-matrix:
c
c        ... scale pop-wt'd eigenvectors with eigenvalues
c        do  i=0,ni
c          ki = i*NB
c          kii = (i*npcol + mycol)*NB
c          nj = MIN( NB, numpop-kii )
c          if( nj .ge. 0 )then
c            do  j=1,nj
c              do  jj=1,norb
c                vmat(jj,ki+j,2) = vmat(jj,ki+j,2)*eigval(kii+j,ispin)
c              enddo
c            enddo
c          endif
c        enddo
crpm    Replace the explicit math with a call to PDSCAL:
        do  i=1,numpop
          call PDSCAL( norb, eigval(i,ispin),vmat(1,1,2),1,i,DESCA, 1 )
        enddo
c
        call PDGEMM( 'n','t', norb,norb,numpop,one,
     $   vmat(1,1,2),1,1,DESCA, vmat(1,1,1),1,1,DESCA,
     $   zero, vmat(1,1,4),1,1,DESCA )
c
c       Scale dmat and emat by weight function:
        wt = wtk(1) 
crpm        call DSCAL( 2*mat, wt, vmat(1,1,3),1 )
        do  i=1,norb
          call PDSCAL( norb,wt,vmat(1,1,3),1,i,DESCA, 1 )
          call PDSCAL( norb,wt,vmat(1,1,4),1,i,DESCA, 1 )
        enddo
c
c     Zero out the imaginary portion (upper triangle) of the two matrices
crpm  Using the PDLASET routine to do this. This is a slightly non-standard
c     way to call the routine, but it insures that the upper triangle
c     without the diagonal is zeroed. The default call also zeros the diagonal,
c     so for this call, you basically give arguments for the first off-diag
c     element, and then zero N-1 elements in the triangle.
        call PDLASET( 'U',norb-1,norb-1,zero,zero,
     $       vmat(1,1,3),1,2,DESCA )
        call PDLASET( 'U',norb-1,norb-1,zero,zero,
     $       vmat(1,1,4),1,2,DESCA )
c        nic = (norb-1) / (NB*npcol)
c        nir = (norb-1) / (NB*nprow)
c        do  ic=0,nic
c          do  ir=0,nir
c            jc = ic*npcol + mycol
c            jr = ir*nprow + myrow
c            nc = MIN( NB, norb - jc*NB )
c            nr = MIN( NB, norb - jr*NB )
c            kr = ir*NB
c            kc = ic*NB
c            if( jc .gt. jr )then
c              do  i=1,nr
c                do  j=1,nc
c                  vmat(i+kr,j+kc,3) = zero
c                  vmat(i+kr,j+kc,4) = zero
c                enddo
c              enddo
c            elseif( jc. eq. jr )then
c              do  i=1,nr-1
c                do  j=i+1,nc
c                  vmat(i+kr,j+kc,3) = zero
c                  vmat(i+kr,j+kc,4) = zero
c                enddo
c              enddo
c            endif
c          enddo
c        enddo
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c               Write out density-matrix, e-matrix:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c       If we skipped constructing zero-matrices, return here:
 2990   continue
c
c       Non-context processors return here:
 2995   continue
        if( itimlog.gt.2 ) call TIMER( 'erlsolv/dmat compute' )
c
c       Save density matrix dmat to idmatfl:
c
        if( myrow .ne. -1 )then
c         Reassemble full dmat on master from distributed dmat
          if( mat .le. maxblockr )then
            call PDGEMR2D( norb,norb,
     $           vmat(1,1,3), 1,1, DESCA,
     $           tempio, 1,1, DESCWORK,  ICONTXT )
          else
c           Should be good as long distributed matrix <2GB on any processor.
            nblock = maxblockr/norb
            nleft = norb - nblock
            noffset = 1
            moffset = 1
 3156       continue
            if( noffset.le.norb )then
              call PDGEMR2D( norb,nblock,
     $             vmat(1,1,3), 1, noffset, DESCA,
     $             tempio(moffset), 1,1, DESCWORK,  ICONTXT )
              noffset = noffset + nblock
              moffset = moffset + norb*nblock
              if( nleft.ge.nblock )then
c               Have enough to do another full block.
                nleft = nleft - nblock
              else
c               Have less than a full block to send.
                nblock = nleft
                nleft = 0
              endif
              goto 3156
            endif
          endif
        endif
        if( itimlog.gt.3 ) call TIMER( 'erlsolv/dmat collect' )
c
Cc       Broadcast dmat from master to all processors
C        call MPI_Bcast( vmat(1,1,1), mat, MPI_DOUBLE_PRECISION,
C     $   0, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'perlsolv: Bcast dmat' )
C        if( itimlog.gt.3 ) call TIMER( 'erlsolv Bcast dmat' )
c
        if( iproc .eq. masterp )then
c         Only master saves dmat to idmatfl
          call WRITBIG( idmfile, mat, tempio )
          if( itimlog.gt.3 ) call TIMER( 'erlsolv/dmat write' )
        endif
c
c       Save e-matrix emat to iematfl:
c
        if( myrow .ne. -1 )then
c         Reassemble full emat on master from distributed emat
          if( mat .le. maxblockr )then
            call PDGEMR2D( norb,norb,
     $           vmat(1,1,4), 1,1, DESCA,
     $           tempio, 1,1, DESCWORK,  ICONTXT )
          else
c           Should be good as long distributed matrix <2GB on any processor.
            nblock = maxblockr / norb
            nleft = norb - nblock
            noffset = 1
            moffset = 1
 3157       continue
            if( noffset.le.norb )then
              call PDGEMR2D( norb,nblock,
     $             vmat(1,1,4), 1, noffset, DESCA,
     $             tempio(moffset), 1,1, DESCWORK,  ICONTXT )
              noffset = noffset + nblock
              moffset = moffset + norb*nblock
              if( nleft.ge.nblock )then
c               Have enough to do another full block.
                nleft = nleft - nblock
              else
c               Have less than a full block to send.
                nblock = nleft
                nleft = 0
              endif
              goto 3157
            endif
          endif
        endif
        if( itimlog.gt.3 ) call TIMER( 'erlsolv/emat collect' )
c
Cc       Broadcast emat from master to all processors
C        call MPI_Bcast( vmat(1,1,1), mat, MPI_DOUBLE_PRECISION,
C     $   0, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'perlsolv: Bcast emat' )
C         if( itimlog.gt.3 ) call TIMER( 'erlsolv/emat Bcast' )
c
        if( iproc .eq. masterp )then
          call WRITBIG( iemfile, mat, tempio )
          if( itimlog.gt.3 ) call TIMER( 'erlsolv/emat write' )
        endif
c
c       Switch to spin-dn:
        idmfile = idmatsfl
        iemfile = iematsfl
c
c       Close loop over spins:
 4000 continue
c
cxxx: should GRIDEXIT be used to clear context for out-of-context procs, too?
      if( myrow .ne. -1 ) call BLACS_GRIDEXIT( ICONTXT )
c
      if( itimlog.gt.1 ) call TIMER( 'erlsolv/rhoij  ' )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DISTRMAT
c
c
      subroutine DISTRMAT( norb, mat, maxblockr,
     $ DESCWORK,DESCA,ICONTXT,
     $ rmatfull,rmatdist )
c---------------------------------------------------------------
c Purpose: distribute matrix over processes, with >2GB handling
c
c **** this is just draft, not tested or invoked yet *****
c
c Revision history:
c   2Jan12-PAS/2.64: extracted from ACP code in persolv
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c Input:
      INTEGER  norb, mat, maxblockr
      INTEGER  ICONTXT
      INTEGER  DESCA(9), DESCWORK(9)
      REAL*8   rmatfull(*)
c Output:
      REAL*8   rmatdist(*)
c Local declarations:
      INTEGER  nblock, nleft,noffset,moffset
c
c >>>> EXECUTABLE CODE:
c
crpm_tp          call PDGEMR2D( norb,norb, vmat(1,1,4), 1,1,
c
      if( mat .le. maxblockr )then
c       Matrix is small enough to be distributed monolithically via BLACS ...
c
        call PDGEMR2D( norb,norb,
     $                 rmatfull, 1,1, DESCWORK,
     $                 rmatdist, 1,1, DESCA, ICONTXT )
c
      else
c       Distribute the matrix in (column) stripes that are smaller than 2GB each.
c       Should be good as long the distributed matrix is no larger than 2GB on any processor.
c
        nblock = maxblockr / norb
        nleft = norb - nblock
        noffset = 1
        moffset = 1
  151   continue
        if( noffset.le.norb )then
          call PDGEMR2D( norb,nblock,
     $                   rmatfull(moffset), 1,1,       DESCWORK,
     $                   rmatdist,          1,noffset, DESCA, ICONTXT )
          noffset = noffset + nblock
          moffset = moffset + norb*nblock
          if( nleft.ge.nblock ) then
c           Have enough to do another full block.
            nleft = nleft - nblock
          else
c           Have less than a full block to send.
            nblock = nleft
            nleft = 0
          endif
          goto 151
        endif
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PSOLVSIZES
c
c
      subroutine PSOLVSIZES( ICONTXT, norb,neig,nlocrow,nloccol,lwork )
c---------------------------------------------------------------
c Purpose: determine solver needs for serial solver (dsygvx)
c Written: Richard P. Muller
c Revision history:
c  17Jan13-PAS/2.64: tidied up a little during kpp solver work
c---------------------------------------------------------------
c
c On input:
c    ICONTXT:  <0 indicates need to create and probe context for info
c            0,>0 indicates probe existing context (ICONTXT) for info
c Note:
c     The output here is only the _minimum_ required to solve,
c     and not the amount needed to make this optimal - 17Jan13
c
      IMPLICIT NONE

crpm  This is only here for the MPI_ALLREDUCE command. If we don't need
c     the command, the include statement can go:
      INCLUDE  'mpif.h'
c
      INTEGER     NB
      PARAMETER  (NB=32)
c     NB      Blocking size (=32)
c
c Functions (these are scalapack functions)
      INTEGER  NUMROC,ICEIL
c
c Input arguments
      INTEGER  ICONTXT
      INTEGER  norb,neig
c     norb    Number of orbitals
c     neig    Number of eigenvalues to compute
c
c Output arguments
      INTEGER  nlocrow,nloccol,lwork
c     nlocrow Number of rows required for local storage
c     nloccol Number of columns required for local storage
c     lwork   Size needed for the secratch array WORK (minimum)
c
c Local variables
      INTEGER  nn,np0,mq0,nprow,npcol,myrow,mycol,nenb,
     $         nprocs,iproc,masterp,icomm,npsolver,JCONTXT
c
      DOUBLE PRECISION  rnprocs
      LOGICAL  init_blacs
c
c >>>> EXECTUABLE CODE:
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
c
c  Create a close to square grid of processors for solver
      rnprocs = DBLE( nprocs ) + 0.1
      npcol = SQRT( rnprocs )
      nprow = nprocs / npcol
      npsolver = nprow*npcol
c
c  Create the BLACS context if we need to:
      if( ICONTXT .lt. 0 )then
c       We do not have an existing context ... create one
        JCONTXT = icomm
        call BLACS_GRIDINIT( JCONTXT, 'R', nprow, npcol )
        init_blacs = .true.
      else
c       We have an existing context, use it ...
        JCONTXT = ICONTXT
        init_blacs = .false.
      endif
c  Get information about BLACS processor grid:
      call BLACS_GRIDINFO( JCONTXT, nprow, npcol, myrow, mycol )
c
      nn = MAX( norb,NB )
      nenb = MAX( neig,NB )
      np0 = NUMROC( nn,NB,0,0,nprow )
      mq0 = NUMROC( nenb,NB,0,0,npcol )
c
      nlocrow = NUMROC( norb,NB,myrow,0,nprow )
      nloccol = NUMROC( norb,NB,mycol,0,npcol )
c
      lwork = 5*norb + MAX( 5*nn,np0*mq0 + 2*NB*NB)
     $     + ICEIL( neig,nprow*npcol )*nn
c
c     nlocrow/col have to be at least 1, since we're using them
c     as the leading dimension of matrices
      nlocrow = MAX( nlocrow, 1 )
      nloccol = MAX( nloccol, 1 )
c
      if( init_blacs )then
c       Close up this BLACS Context if we created it
        if( myrow .ne. -1 ) call BLACS_GRIDEXIT( JCONTXT )
      endif
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
