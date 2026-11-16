c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PSCHROED2
c
c
      subroutine PSCHROED2( do_eigvecs,
     $ idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $ norb,nstate,nk, eigval, rwork, seval, scut,
     $ vmat, wmat, iwork )
c---------------------------------------------------------------
c Purpose: solve Schroedinger equation given input ham and ovlp,
c          with k-resolved distribution over processors
c          using LAPACK on each processor.
c
c Written: Peter A. Schultz,  8-April-1999, for v2.34
c
c Revision history:
c   2May12-ACP/2.64: Added flag to skip eigvecs to support non-SCF
c  25Jun07-PAS/2.60: merge tp and serial
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  18Aug05-PAS/2.59: master-only matrix files; read dHams
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c  22Apr03-DMC: This version of a parallel schroed routine is designed
c               to take advantage of nk>1, by parallelizing on k.
c  24Jul01-PAS/2.49: spin-polarized dft
c  21Jun01-PAS/2.48: replace STOPs
c  17May01-PAS/2.47: extract writes
c  23Aug99-PAS/2.38: install zhegvx expert driver eigensolve
c  10May99-PAS/2.35: eigvals cut from ivecfl output
c   9Apr99-PAS/2.34: complex lapack eigensolver installed
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  10Dec98-PAS/2.29: promote source to explicit double precision
c   5Sep98-PAS/2.25: print out spectral limits
c  14Jul98-PAS/2.23: option to use unchol'd delta-Ham
c  12Aug94-PAS/2.11: cleanup
c---------------------------------------------------------------
c
c When this routine is done
c   (1) it writes out eigenvectors to ivecfl
c   (2) it returns eigenvalues in "eigval()"
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
ctp: Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c MPI information
      INTEGER status(MPI_STATUS_SIZE)
c
c Input:
      LOGICAL    do_eigvecs
c Output: eigenvalue spectrum
      DIMENSION  eigval(nstate,nk)
c Scratch arrays:
      DIMENSION  rwork(*)
      DIMENSION  seval(*)
      DIMENSION  iwork(norb,*)
      DIMENSION  vmat(norb*norb)
      DIMENSION  wmat(norb*norb,*)
c
c Local declarations:
      CHARACTER  jobz*1,uplo*1
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
cxxx: declare dev-version of pschroedkp:
      if( IWR.ge.0 ) write(IWR,*)
     $ 'PSCHROED2: k-parallel single-processor solver, nk=',nk
      mat = norb*norb
c
c Do eigenvectors with jobz='V' (jobz="N" turns off eigenvectors):
      if( do_eigvecs )then
c       Compute eigenvalues and eigenvectors
        jobz = 'V'
      else
c       Compute eigenvalues only
        jobz = 'N'
      endif
      uplo = 'L'
c
c  Get Local MP info
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
      if( iproc .eq. masterp )then
        REWIND( unit=iham0fl )
        REWIND( unit=iovlpfl )
      endif
c
c    Bloch vector loop:
c
      do 1000 kblk=1,nk, nprocs
        do 1010 kk=1,nprocs
          k = kk + kblk - 1
          if( k .gt. nk ) goto 1011
c
c         This sets it up so that masterp goes *last*
          kproc = MOD( kk, nprocs )
c
          if( iproc .eq. masterp )then
c
c           Retrieve delta-Hamiltonian from idhamfl
            call READBIG( idhamfl, mat, vmat )
c
c           Retrieve reference Hamiltonian from iham0fl
            call READBIG( iham0fl, mat, wmat(1,1) )
c
c           Add del-Ham to ref-Ham to get total Ham
            do  ij=1,mat
              wmat(ij,1) = wmat(ij,1) + vmat(ij)
            enddo
c
c           Retrieve overlap matrix from iovlpfl
            call READBIG( iovlpfl , mat, wmat(1,2) )
c
            if( kproc .ne. masterp )then
c             Send H and S together to node that will do this work
              lenmsg = 2*mat
              call MPSENDR8( kproc, lenmsg, wmat(1,1), icomm )
C              call MPI_SEND( wmat(1,1), 2*mat, MPI_DOUBLE_PRECISION,
C     $         kproc,k, icomm, ierr )
C              if( ierr.ne.0 ) call STOPXERR( 'pschroed2 - Send err' )
            endif
c
          elseif( iproc .eq. kproc )then
c
            lenmsg = 2*mat
            call MPRECVR8( masterp, lenmsg, wmat(1,1), icomm )
C            call MPI_RECV( wmat(1,1), lenmsg,  MPI_DOUBLE_PRECISION,
C     $       masterp,k, icomm, status, ierr )
C              if( ierr.ne.0 ) call STOPXERR( 'pschroed2 - Send err' )
c
          endif
c
c         Only the node(s) with work will continue
          if( iproc .ne. kproc ) goto 1010
c
c         Load overlap into lower triangle full complex matrix wmat(3:4)
          call C16LOW( norb, wmat(1,2), wmat(1,3) )
c
c         Load Ham into lower triangle of full complex matrix wmat(1:2)
          call DCOPY( mat, wmat(1,1),1, vmat,1 )
          call C16LOW( norb, vmat     , wmat(1,1) )
c
c         Need 3*norb-2 for r*8 workspace "rwork"
c         lwork is (complex) length for workspace "l"
          lwork = (mat) / 2
c
c         Ham is in wmat(,1:2), overlap in wmat(,3:4)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                  Solve Schroedinger equation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
C         We will always use the expert driver for parallel.
c         We use idosolv only to control how many eigenvalues/vectors
c         to solve for.
c
          if( idosolv .eq. 1 )then
            jstate = norb
          elseif( idosolv .eq. 2 )then
            jstate = nstate
          else
            call STOPXERR( 'pschroed2 - invalid solver option' )
          endif
c
c          ***** Selected eigenvalue Hermitian eigenproblem *****
c
          call ZHEGVX( 1,jobz, 'I', 'L',
     $     norb, wmat(1,1),norb,         wmat(1,3),norb,
c -->            H(2,n,n)-i              S(2,n,n)-i
     $     -100d0,-20d0, 1,jstate, 2*DLAMCH('S'), nstout,
     $     eigval(1,k),
     $     wmat(1,5),norb,    vmat     ,lwork, rwork,
c -->      eigvec(2,n,nstate) work(2,lwork)    (7n)
     $     iwork(1,2),iwork, INFO )
c -->      int(5n)    int(n)
c
          if( INFO .ne. 0 )then
            call STOPXERR( 'pschroed2 - diagonalization error' )
          endif
c
          if( do_eigvecs )then
c           Put eigenvectors in the right place ...
            call DCOPY( 2*norb*nstate, wmat(1,5),1, wmat(1,1),1 )
          endif
c
 1010   continue
 1011   continue
c
c       Now collect results of work at master, and store from there
c
        do 1020 kk=1,nprocs
          k = kk + kblk - 1
          if( k .gt. nk ) goto 1021
c
          kproc = MOD( kk, nprocs )
c
c         * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c          Send eigenvalues to everyone from kproc that owns this k-pt
c         * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
          call MPBCAST8( kproc, nstate, eigval(1,k), icomm )
c
c         Without eigvecs, skip their retrieval and storage
          if( .not. do_eigvecs ) goto 1020
c
c         * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                          Record eigenvectors
c         * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
          nvec = 2*norb*nstate
c
c         Transfer the eigenvectors to masterp for storage
          if( iproc.eq.kproc .and. iproc.ne.masterp )then
c           Send eigvecs to masterp (safe - Bcast above shields)
            call MPSENDR8( masterp, nvec, wmat(1,1), icomm )
c
          elseif( iproc.eq.masterp .and. kproc.ne.masterp )then
c           Receive eigvecs from kproc (safe - Bcast above shields)
            call MPRECVR8( kproc, nvec, wmat(1,1), icomm )
c
          elseif( iproc.eq.masterp .and. kproc.eq.masterp )then
c           Receives stomped on masterp eigvecs, must refresh
c           (because masterp does last k-pt, not first)
            call DCOPY( nvec, wmat(1,5),1, wmat(1,1),1 )
          endif
c
c         Save eigenfunction coefficients to ivecfl
          if( iproc .eq. masterp )then
            call WRITBIG( ivecfl, 2*norb*nstate, wmat(1,1) )
          endif
c
 1020   continue
 1021   continue
c
 1000 continue
C      call MPBARRIER( icomm )
C      call MPI_ABORT( icomm, ierrors, istatus )
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
