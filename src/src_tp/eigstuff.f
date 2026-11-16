c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EIGSTUFF
c
c
c===============================================================
c Purpose: a set of routines to help handle the parallel solver stuff.
c Written; Peter A. Schultz, 9-December-2008, for 2.62 (k-parallel)
c Revision history:
c  none
c===============================================================
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EIGPROCS
c
c
      subroutine EIGPROCS( nprocs, nprow,npcol, npsolver )
c---------------------------------------------------------------
c Purpose: return nearly square processor grid for solver
c---------------------------------------------------------------
c
c Generate the best (most square) processor grid within nprocs
c
c nprocs = total number of processors available in this communicator
c npcol = number of columns in solver processor grid
c nprow = number of rows in solver processor grid
c npsolver = npcol*nprow = number of processors in solver grid
c
      IMPLICIT NONE
c
c Input:
      INTEGER  nprocs
c Output:
      INTEGER  npcol,nprow, npsolver
c
c Local:
      REAL     rnprocs
c
c >>>> EXECUTABLE CODE:
c
      rnprocs = DBLE( nprocs ) + 0.1
      npcol = SQRT( rnprocs )
      nprow = nprocs/npcol
      npsolver = nprow*npcol
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EIGMDIST
c
c
      subroutine EIGMDIST( matdist, norb, npsolver )
c---------------------------------------------------------------
c Purpose: return default distributed matrix size
c Revision history:
c  13Jan13-PAS/2.63a: fix a dumb mistake, very narrow failure range
c---------------------------------------------------------------
c
c norbdist = rank of matrices for which we will use a reduced
c            distributed memory allocation.
c
c norbdist is set to a somewhat larger value to ensure that
c we do not run into block-size problems, but small enough
c that using a full matrix assignment is sure to give no problems.
c
      IMPLICIT NONE
c
      INTEGER  NB
      DATA     NB / 32 /
c
c Output:
      INTEGER  matdist
c Input:
      INTEGER  norb, npsolver
c
c Local:
      INTEGER  mat, matpercpu
      INTEGER  norbdist
      DATA     norbdist / 800 /
c
c >>>> EXECUTABLE CODE
c
      mat = norb*norb
      if( norb .ge. norbdist .and. npsolver .gt. 1 )then
        matpercpu = (mat+npsolver-1) / npsolver
c       Add 32*norb to give slack for block sizes coarseness
c        and uneven distribution across processors
        matdist = matpercpu + NB*(norb+NB)
      else
c       Just give the full matrix space
        matdist = mat
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPPMEMX
c
c
      subroutine KPPMEMX( per_process, do_eigvecs,
     $ norb,nstate,  memeig )
c---------------------------------------------------------------
c Purpose: compute memory requirements (size of wmat argument) of 
c          PSCHROEDKP for wkmem using Scalapack calls.
c
c Written:  Andrew C. Pineda,  9-February-2010 for v2.62
c
c Revision history:
c  4Jan13-PAS/2.64: adapted from pschroedkpsz by acp
c---------------------------------------------------------------
c
c Return argument:
c    memeig = process value (per_process = .true.)
c           = max value over all processes (per_process = .false.)
c
      INTEGER     NB
      PARAMETER  (NB=32)
      INCLUDE  'mpif.h'
c
      INTRINSIC  MAX
c
c Input arguments
      INTEGER  norb
      LOGICAL  per_process
      LOGICAL  do_eigvecs
c
c Output arguments
      INTEGER  memeig
c
c Local variables
      integer  memeiglocal
      logical  per_process_local
c
      integer  mat, matdist
      integer  nwpreph, nwpreps, nwprep, nwsolv
      integer  nwclose, nwpsch, nwleft
      integer  iwHmat, iwOvlp, iwVecs
      integer  iwZwork, iwRwork, iwIwork, iwiwk1, iwlast, iwEnd
      integer  lwzwork,lzwork, lrwork, liwork
      integer  iwMat1, iwMat2, iwVeco
c
      integer  CONTEXT_K
      integer  nprocs_k, icomm
      integer  nprow_k, npcol_k, npsolver_k
c
      integer  ierr
c
      LOGICAL  USE_SCALAPACK_FORMULAS
      DATA    USE_SCALAPACK_FORMULAS / .true. /
C      LOGICAL  OPTIMAL_EIGENVECTORS
C      DATA    OPTIMAL_EIGENVECTORS   / .true. /
c
      LOGICAL  DEBUG_MEM
      DATA     DEBUG_MEM / .false. /
c
c >>>> EXECUTABLE CODE:
c
      if( DEBUG_MEM )then
        call FLGETIWR( IWR )
        write(IWR,*) 'DEV/KPPMEMX: do_eigvecs=',do_eigvecs
      endif
      mat = norb*norb
c  Get Local MP info
c     Image level info
      call MPCOMM( icomm )
c     K-point level info
      call MPNODES_K( nprocs_k )
c
      memeiglocal = 0
c
      if( nprocs_k .eq. 1 )then
c       Using Lapack serial eigensolver
c 
c       Only one processor assigned to the MPI context this processor
c       belongs to, set up memory to use LAPACK. 
c
        matdist = mat
        lwzwork = MAX( 4*norb, mat )
        lrwork  = 7*norb
        liwork  = 5*norb
c
c       Close serial solver memory needs
      else
c       USING SCALAPACK
c
c       Get the process grid dimensions, given the context size
c
        call EIGPROCS( nprocs_k, nprow_k, npcol_k, npsolver_k )
c
c       Determine pzhegvx workspace requirements
c
        if( USE_SCALAPACK_FORMULAS )then
          if( DEBUG_MEM )then
c           Use default mpi/blacs-avoiding naive work space:
            call PSCHWKSZ0(            npsolver_k, do_eigvecs,
     $       norb, matdist, lwzwork, lzwork, lrwork, liwork )
            lwork = lzwork
          endif
c         Get EXACT mpi/blacs-aware work space using this context:
C          CONTEXT_K = -1
          call MPCOMM_K( icomm_k )
          CONTEXT_K = icomm_k
          call BLACS_GRIDINIT( CONTEXT_K, 'R', nprow_k, npcol_k )
          call PSCHWKSZX( CONTEXT_K, npsolver_k, do_eigvecs,
     $     norb, matdist, lwzwork, lzwork, lrwork, liwork )
          lwork = lzwork
          call BLACS_GRIDEXIT( CONTEXT_K )
        else
c         Use default mpi/blacs-avoiding naive work space:
          call PSCHWKSZ0(            npsolver_k, do_eigvecs,
     $     norb, matdist, lwzwork, lzwork, lrwork, liwork )
          lwork = lzwork
        endif
      endif
c
c     Initialize memory allocations within icomm_k.
c     MEM-Solver: d-C16-H, d-C16-S, d-C16-Vec, d-C16-Zwork, d-R8-Rwork)
c
cxxx264: rednudant (unused here) temp storage for solvers:
      iwHmat  = 1
      iwOvlp  = iwHmat + 2*matdist
      iwVecs  = iwOvlp + 2*matdist
      iwZwork = iwVecs + 2*matdist
      iwRwork = iwZwork + lwzwork
      iwIwork = iwRwork + lrwork
      iwiwk1  = iwIwork + liwork
      iwlast  = iwiwk1  + norb
      iwEnd   = iwlast  - 1
c
      iwMat1  = iwOvlp
      iwMat2  = MAX( iwMat1 + mat, iwVecs )
      iwVeco  = iwZwork
c
c     The solver memory allocations have the following requirements:
      nwpreph = 2*matdist + 3*mat
      nwpreps = 4*matdist + 2*mat
      nwprep = MAX( nwpreph, nwpreps )
      nwsolv = 6*matdist + lwzwork + lrwork + liwork + norb
c     nwclose = 6*matdist + 2*norb*nstate
c     because all vectors are brought back into full norb*norb ...
cxx264: nwclose (eigeenec store) needs to be dealt with (norb*nstate)
      nwclose = 6*matdist + 2*norb*norb
      nwpsch = MAX( nwprep, nwsolv, nwclose )
c
      nwleft = nwpsch - nwsolv
      if( DEBUG_MEM )then 
        write(IWR,*) 'DEV/KPPMEMX: spare solver mem=',nwleft
      endif
c
      memeiglocal = MAX( iwEnd, iwMat2+2*mat, iwVeco+2*mat )
c
c >>>>> Put it all together and report to final result
c
      per_process_local = per_process
      if( per_process_local ) then
c       Report the local memeig value
        memeig = memeiglocal
      else
c       Report the maximum value of memeig for any process in icomm
c       Need barrier to prevent a crash
        call MPBARRIER( icomm )
        icount = 1
        call MPI_ALLREDUCE( memeiglocal, memeig, icount,
     $       MPI_INTEGER, MPI_MAX, icomm,
     $       ierr )
      endif
      if( DEBUG_MEM )then
        write(IWR,'(a,4i12)') 'DEV/KPPMEMX: works,Z,R,I=',
     $                         lwzwork, lrwork, liwork
        write(IWR,'(a,3i12)') 'DEV/KPPMEMX: nwsch,memeig=',
     $                         nwpsch, memeiglocal, memeig
        write(IWR,*) 'DEV/KPPMEMX: icomm,nprocs_k=',icomm,nprocs_k
        call FLUSH( IWR )
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PSOLWKSZX
c
c
      subroutine PSCHWKSZX( icontxt, nprocs, do_eigvecs,
     $ norb, matdist, lwzwork, lzwork, lrwork, liwork )
c---------------------------------------------------------------
c Purpose: compute exact memory requirements for PSCHROEDKP solvers,
c          using Scalapack (p)zhegvx calls.
c      
c Written:  Peter A. Schultz, 10-January-2013, for 2.64
c           Derived from PSCHROEDKPSZ by Andrew C. Pineda (9Feb10)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c  Work array space for
c    (1)  serial solver zhegvx
c    (2)  parallel solver pzhegvx
c
      INTEGER     NB
      PARAMETER  (NB=32)
C      INCLUDE  'mpif.h'
c
      INTRINSIC  MAX
c
c Input arguments
      INTEGER  icontxt
      INTEGER  nprocs
      LOGICAL  do_eigvecs
c
      INTEGER  norb
c
c Output arguments
      INTEGER  matdist
      INTEGER  lwzwork, lzwork, lrwork, liwork
c
c Local variables
      INTEGER  mat
      INTEGER  nprow, npcol, npsolver
c
      INTEGER  CONTEXT
      INTEGER  myrow, mycol
      INTEGER  NN, NP0, NQ0, NEIG, MQ0
      INTEGER  ANB, SQNPC, NPS, NHETRD_LWOPT, NHEGST_LWOPT
c
      LOGICAL  USE_SCALAPACK_FORMULAS
      LOGICAL  OPTIMAL_EIGENVECTORS
      DATA    USE_SCALAPACK_FORMULAS / .true. /
      DATA    OPTIMAL_EIGENVECTORS   / .true. /
c
      LOGICAL  DEBUG_MEM
      DATA     DEBUG_MEM  / .false. /
c
c >>>> EXECUTABLE CODE:
c
      if( DEBUG_MEM )then
        call FLGETIWR( IWR )
        write(IWR,*) 'DEV/WKSZ-X: norb,nprocs=',norb, nprocs
        write(IWR,*) 'DEV/WKSZ-X: do_eigvecs=',do_eigvecs
      endif
c
      mat = norb*norb
c
      if( nprocs .eq. 1 )then
c       Using Lapack serial eigensolver
c 
c       Only one processor assigned to the MPI context this processor
c       belongs to, set up memory to use LAPACK. 
c
        matdist = mat
cxxx264        lwzwork = 2*MAX( 2*norb, mat ) ! I am sure we don't need 2mat -pas
        lwzwork =   MAX( 4*norb, mat )
        lzwork  = lwzwork / 2
        lrwork  = 7*norb
        liwork  = 5*norb
cxxx264 - need to scope a do_eigvecs-aware memory for single-proc
C        memeig = 2*mat + 2*mat + 2*norb*nstate + mat + 7*norb
c
c       Close serial solver zhegvx memory needs
      else
c       Parallel, USING SCALAPACK
c
c       Get the process grid dimensions, given the context size
        call EIGPROCS( nprocs, nprow, npcol, npsolver )
c       NOTE: following is a crude approximation to matdist, NOT exact
        call EIGMDIST( matdist, norb, npsolver )
c
        if( DEBUG_MEM )then
          write(IWR,'(a,4i12)')'WKSZX: mat,matdist=',mat,matdist
          write(IWR,'(a,4i6)')'WKSZX: nprocs,nprow,npcol,npsolver=',
     $                                nprocs,nprow,npcol,npsolver
        endif
c
c   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c    Simplified appxximate work memory requirements
c   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c      Simplified naive values
        lwzwork = 2*matdist
        lzwork = lwzwork / 2
        lrwork = 2*matdist
        liwork = 6*norb
        if( .not. do_eigvecs )then
c         No eigenvectors, this becomes much less:
          lzwork = (NB+1) * (norb+1)
          lrwork =9*norb
        endif
c
c   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c    Determine EXACT memory requirements, from Scalapack docs (per ACP)
c   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c      ***** Get some of the basic components of the formulas
c
        NN = MAX( norb, NB, 2 )
c       NUMROC is a SCALAPACK function to compute the amount of a 
c       distributed array that is resident on a given processor.
c       It is EXACT.
c       NP0 == number of rows assigned to root process
c           == upper bound on rows assigned to any process.
        NP0   = NUMROC( norb, NB, 0, 0, nprow )
c       NQ0 == number of columns assigned to root process
c           == upper bound on columns assigned to any process.
        NQ0   = NUMROC( norb, NB, 0, 0, npcol )
        NEIG  = norb
        MQ0   = NUMROC( MAX( NEIG,NB,2 ), NB, 0, 0, npcol )
c
        if( DEBUG_MEM )
     $  write(IWR,'(a,6i9)') 'WKSZX: norb,nb,nn,np0,nq0.mq0=',
     $                               norb,nb,nn,np0,nq0,mq0
c
c      ***** Minimum size for PZHEGVX argument WORK: lzwork
c            lzwork is the number of complex elements in wmat(iwZwork)
c
        if( do_eigvecs ) then
          lzwork = norb + NB*(NP0+MQ0+NB)
c
          if( DEBUG_MEM ) write(IWR,*) 'do_eigvecs: base lzwork=',lzwork
          if( OPTIMAL_EIGENVECTORS )then
c           For optimal performance, lzwork needs to be bigger
            if( icontxt .ge. 0 )then
c             Query existing context
              CONTEXT = icontxt
            else
c             Create trial context to get data
              call BLACS_GRIDINIT( CONTEXT, 'R', nprow, npcol )
            endif
            call BLACS_GRIDINFO( CONTEXT, nprow, npcol, myrow, mycol )
            if( DEBUG_MEM )
     $      write(IWR,*) 'optvecs: nprow,npcol=',nprow,npcol,myrow,mycol
            if( myrow .ne. -1 )then
c             OK we're in the BLACS grid/context
              ANB   = PJLAENV( CONTEXT, 3, 'PZHETTRD',
     $             'L', 0, 0, 0, 0 )
              SQNPC = INT( SQRT( DBLE( nprow * npcol ) ) )
              NPS   = MAX( NUMROC( norb, 1, 0, 0, SQNPC ), 2*ANB )
              NHETRD_LWOPT = 2*(ANB+1)*(4*NPS+2) + (NPS+4)*NPS
              NHEGST_LWOPT = (2*NP0+NQ0+NB)*NB
              lzwork = MAX( lzwork, norb+NHETRD_LWOPT, NHEGST_LWOPT )
cxxx264:      Moved gridexit outside this if-block (to balance gridinit) - A bug?
cxxx264              call BLACS_GRIDEXIT( CONTEXT )
            endif
c           Exit the BLACS, we'll recreate it on each call to the solver.
            if( icontxt .lt.0 ) call BLACS_GRIDEXIT( CONTEXT )
c
            if( DEBUG_MEM ) write(IWR,*) 'end optvec: lzwork=', lzwork
          endif
        else
c         Size of WORK array when computing eigenvalues only
          lzwork = norb + MAX( NB*(NP0+1), 3 )
          if( DEBUG_MEM ) write(IWR,*) 'eigval only: lzwork=', lzwork
        endif
c        Factor of 2 in next line  because wmat(iwZwork) is C*16, allocating into R*8
        lwzwork = 2*lzwork
c
c      ***** Minimum size of RWORK: lrwork
c
c       ICEIL(X,Y) is a SCALAPACK function that computes CEILING(DBLE(X)/DBLE(Y)).
        if( do_eigvecs )then
c         Size of RWORK array when computing eigenvalues and eigenvectors
          lrwork = 4*norb + MAX( 5*NN, NP0*MQ0 )
     $         + ICEIL( NEIG, nprow*npcol )*NN
        else
c         Size of RWORK array when computing eigenvalues only
          lrwork = 4*norb + 5*NN
        endif
cxxx264: pas: why was this 2x of lrwork in the ACP code?  Seems wrong (thinks it's c16?).
cxxx264        lrwork = 2*lrwork
c
c        ***** Minumum size of liwork.
c         Noting iwIwork is integer array, wmat is R*8, and sizeof(INTEGER) <= sizeof(REAL*8)
c
        liwork = 6*MAX( norb, nprow*npcol+1, 4 )
c       
c       Close pzhegvx parallel solver work array memory needs
      endif
c
      if( DEBUG_MEM )then
        write(IWR,'(a,3i12)') 'DEVVWKSZ-X: works Z,R,I=',
     $                         lwzwork,lrwork,liwork
        call FLUSH( IWR )
      endif
cxx:end
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
