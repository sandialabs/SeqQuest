c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PSCHROED
c
c
      subroutine PSCHROED( do_eigvecs,
     $ idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $ norb,nstate,nk, eigval, wksml,
     $ nspin, icluster, gapp,
     $ wmat )
c---------------------------------------------------------------
c Purpose: solve Schroedinger equation given input ham and ovlp,
c          single solve, distributed over processors
c
c Written: Peter A. Schultz,  8-April-1999, for v2.34
c
c Revision history:
c  22Aug17-PAS/2.68: fix idle-proc (solver context) bug
c  10Jan13-PAS/2.64: cosmetics, un-drivered, and integrated
c   2May12-ACP/2.64: Added alternate entry point to support non-SCF portions of code
c   2May12-ACP/2.64: Workaround for 2GB limit in ScaLAPACK
c   9Dec08-PAS/2.62: reduce mem complex-parallel solve
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  18Aug05-PAS/2.59: master-only matrix files, dHam from disk
c  27Oct04-PAS/tp0.1: clean-up task-parallel routines
c   7Apr03-DMC: Parallelized by moving loop on nspin into routine,
c               and using ScaLAPACK to solve eigen problem.
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
c Notes:
c  The call parameter idosolv was used to control which version of the
c  LAPACK routine ZHEGV was used.  Here, we always use PZHEGVX.
c  parameter will be used only to control whether the number of
c  eigenvectors/values returned is all (norb) or just nstate.
c  DMC 8Apr03
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c     MAXBLOCK= 2**31/sizeof(DOUBLE COMPLEX) = 2**31/16.
      PARAMETER ( maxblockz = 67 108 864 )
c
ctp: Pull in the mpi stuff
      INCLUDE  'mpif.h'
c
      LOGICAL    do_eigvecs
c Output:
      DIMENSION  eigval(nstate,nk,nspin)
c Scratch:
      DIMENSION  wksml(*)
      DIMENSION  wmat(*)
c Scratch/parallel:
      INTEGER    icluster(*)
      DIMENSION  gapp(*)
c
c Local declarations:
      CHARACTER  jobz*1
      DATA       iprint / 1 /
ctp: Local for ScaLAPACK
      PARAMETER  (NB = 32)
      INTEGER    CONTEXT
      INTEGER    DESCA(9), DESCWORK(9)
c
      LOGICAL    USE_SCALAPACK_FORMULAS
      DATA       USE_SCALAPACK_FORMULAS / .true. /
      LOGICAL    DEBUG_MEM
      DATA       DEBUG_MEM / .false. /
c
c >>>> EXECUTABLE CODE:
c
      if( iprint.gt.0 .or. DEBUG_MEM )then
        call FLGETIWR( IWR )
        write(IWR,*) 'PSCHROED: parallel distributed solver'
      endif
c
c     Do eigenvectors with jobz='V' (jobz="N" turns off eigenvectors):
      if( do_eigvecs )then
c       Compute eigenvalues and eigenvectors
        jobz = 'V'
      else
c       Compute eigenvalues only
        jobz = 'N'
      endif
c
c  Get Local MP info
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
c
c  Initialize ScaLAPACK Grid
c
c  Define processor grid dimensions of largest nearly square rectangular array
c  of processors within nprocs. npsolver is the number of processors in the grid
c  (npsolver.le.nprocs).
c
      call EIGPROCS( nprocs, nprow,npcol, npsolver )
      if( iprint.gt.0 ) write(IWR,'(a,4i8)')
     $ ' SOLVER grid:  procs, rowsXcols=solver =',
     $      nprocs,  nprow,npcol, npsolver
c
c  Create BLACS context using system context, 
c  i.e. the local MP communicator
      CONTEXT = icomm
      call BLACS_GRIDINIT( CONTEXT, 'R', nprow, npcol )
      call BLACS_GRIDINFO( CONTEXT, nprow, npcol, myrow, mycol )
c
      if( myrow .ne. -1 )then
c       Initialize matrix descriptions, for procs in solver context
c
c       Set up the descriptor for a distributed matrix
c       First, need to limit lead dimension in local space
c       nblkrow = ceil(norb/NB) = # of blocks in a row
        nblkrow = (norb+NB-1)/NB
c       npblkrow = ceil(nblkrow/nprow)
        npblkrow = (nblkrow+nprow-1)/nprow
        LLDA = npblkrow*NB
        call DESCINIT( DESCA, norb,norb, NB,NB, 0,0,
     $   CONTEXT, LLDA, INFO )
        if( INFO.ne.0 ) call STOPXERR( 'pschroed - DESCINIT/DESCA' )
c
c       Set up the descriptor for full local matrix
c
c       call DESCINIT( DESC, M, N, MB, NB, RSRC, CSRC, ICTXT, MXLLDA, INFO)
c       MB=M, NB=N means only the processor (RSRC,CSRC) has any data.
        call DESCINIT( DESCWORK, norb,norb, norb,norb, 0,0,
     $   CONTEXT, norb, INFO )
        if( INFO.ne.0 ) call STOPXERR( 'pschroed - DESCINIT/DESCWORK' )
      else
c       This processor not in solver context, and leaves
        goto 2000
      endif
      if( iprint.gt.2 )then
        write(IWR,*) 'DESC= dtype ctxt m n mb nb csrc rsrc lld'
        write(IWR,'(a,9i8)') '  DESCA=',DESCA
        write(IWR,'(a,9i8)') '  DWORK=',DESCWORK
      endif
c
c  Initialize memory allocations:
c
      mat = norb*norb
      call EIGMDIST( matdist, norb, npsolver )
c
c Allocate "simple" memory sizes for solver work arrays (using crude assumptions).
c Work (Zwork) is complex, r8 length (in wmat) and c16 length (in solver) differ 2x
      lwzwork = mat
      lzwork = lwzwork / 2
C      lrwork = mat
      lrwork = matdist
      liwork = 6*norb
c
      if( DEBUG_MEM )then
        write(IWR,'(a,4i10)') '  PSCH: old native works Z,R,I=',
     $                         lwzwork,lrwork,liwork
c       Use default mpi/blacs-avoiding naive work space:
        call PSCHWKSZ0(            npsolver  , do_eigvecs,
     $   norb, matdist, lwzwork, lzwork, lrwork, liwork )
        lwork = lzwork
        write(IWR,'(a,4i10)') '  PSCH: wksz0  works Z,R,I=',
     $                         lwzwork,lrwork,liwork
      endif
c
      if( USE_SCALAPACK_FORMULAS )then
        if( DEBUG_MEM .or. iprint.gt.2 )
     $  write(IWR,*) '  PSCH: WKSZ/SCAL FORMULAS, CONTEXT=',CONTEXT
c       Get EXACT mpi/blacs-aware work space using this context:
        call PSCHWKSZX( CONTEXT  , npsolver  , do_eigvecs,
     $   norb, matdist, lwzwork, lzwork, lrwork, liwork )
        lwork = lzwork
        if( DEBUG_MEM .or. iprint.gt.2 )
     $  write(IWR,'(a,4i10)') '  PSCH: wkszx  works Z,R,I=',
     $                         lwzwork,lrwork,liwork
      else
c       Use default mpi/blacs-avoiding naive work space:
        call PSCHWKSZ0(            npsolver_k, do_eigvecs,
     $   norb, matdist, lwzwork, lzwork, lrwork, liwork )
        lwork = lzwork
      endif
c
      if( .not. do_eigvecs )then
c       No eigenvectors, this becomes much less:
        lzwork = (NB+1) * (norb+1)
        lrwork =9*norb
      endif
c
c The solver memory allocations have the following requirements:
      nwpreph = 2*matdist + 3*mat
      nwpreps = 4*matdist + 2*mat
      nwprep = MAX( nwpreph, nwpreps )
c     Using a naive allocation for solver work arrays
      nwsolv = 6*matdist + lwzwork + lrwork + liwork + norb
C      nwclose = 6*matdist + 2*norb*nstate
c      because vectors are brought back into full norb*norb ...
      nwclose = 6*matdist + 2*norb*norb
      nwpsch = MAX( nwprep, nwsolv, nwclose )
c
      nwleft = nwpsch - nwsolv
      if( DEBUG_MEM )then
        write(IWR,*) '  PSCHROED: spare solver mem=',nwleft
      endif
c
c MEM-Solver: d-C16-H, d-C16-S, d-C16-Vec, d-C16-Zwork, d-R8-Rwork)
c
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
      if( nwsolv .ne. iwEnd )then
c       Our memory arithmetic has failed this idiot-check
        write(*,*) '***** ERROR/PSCHROED: nwsolv,iwsolv=',nwsolv,iwEnd
        call STOPXERR( 'Code ERROR: nwsolv mem arith in pschroed' )
      endif
c
c MEM-Prep: (d-C16-H, d-C16-S=full-R8-mat, full-C16-mat)
      iwMat1 = iwOvlp
c     Put full-C16-mat after full-R8-mat AND after dist-S
      iwMat2 = MAX( iwMat1 + mat, iwVecs )
c MEM-After: (d-C16-H, d-C16-S, d-C16-Vec, full-C16-Vec)
c     Put full output vectors after dist-vectors:
      iwVeco = iwZwork
c
c Set files to correct locations
c
      if( iproc .eq. masterp )then
        REWIND( unit=ivecfl )
        REWIND( unit=idhamfl )
      endif
c
      do 100 ispin=1,nspin
c
        if( iproc .eq. masterp )then
          REWIND( unit=iham0fl )
          REWIND( unit=iovlpfl )
        endif
c
c       Bloch vector loop:
c
        do 1000 k=1,nk
c
          if( idosolv .eq. 1 )then
            jstate = norb
          elseif( idosolv .eq. 2 )then
            jstate = nstate
          else
            call STOPXERR( 'pschroed - invalid solver option' )
          endif
c
c         Form H0+dH=full Hamiltonian(spin,k), and distribute
c
          if( iproc .eq. masterp )then
c           Retrieve reference Hamiltonian H0 (packed) from iham0fl
            call READBIG( iham0fl, mat, wmat(iwMat2) )
c           Retrieve delta-Hamiltonian dH (packed) from idhamfl
            call READBIG( idhamfl, mat, wmat(iwMat1) )
c           Combine H0 and dH to get full Hamiltonian (packed)
            do  ij=1,mat
              wmat(iwMat1-1+ij) = wmat(iwMat1-1+ij) + wmat(iwMat2-1+ij)
            enddo
c           Load Ham (packed real) into lower triangle of complex matrix
            call C16LOW( norb, wmat(iwMat1), wmat(iwMat2) )
c-->                      n    real*8 nxn    complex*16 nxn (lower triangular)
          endif
c         Distribute complex Ham over processors
          if( norb*norb .le. maxblockz )then
c           The matrix is small enough that pzgemr2d can swallow it whole.    
            call PZGEMR2D( norb,norb,
c                From full matrix on masterp ...
     $           wmat(iwMat2), 1,1, DESCWORK,
c                 ... into distributed matrix:
     $           wmat(iwHmat), 1,1, DESCA,  CONTEXT )
          else
c           Distribute the matrix in (column) stripes that are smaller than 2GB each.
c           Should be good as long the distributed matrix is no larger than 2GB on any processor.
            nblock = maxblockz / norb
            nleft = norb - nblock
            noffset = 1
            moffset = 0
  151       continue
            if( noffset.le.norb )then
              call PZGEMR2D( norb, nblock,
     $             wmat(iwMat2+moffset), 1,1, DESCWORK,
     $             wmat(iwHmat), 1, noffset, DESCA,  CONTEXT )
              noffset = noffset + nblock
              moffset = moffset + 2*norb*nblock
              if( nleft.ge.nblock )then
c               Have enough to do another full block.
                nleft = nleft - nblock
              else
c               Have less than a full block to send.
                nblock = nleft
                nleft = 0
              endif
              goto 151
            endif
          endif
c
c         Retrieve and distribute overlap matrix
c
          if( iproc .eq. masterp )then
c           Retrieve overlap matrix (packed) from iovlpfl
            call READBIG( iovlpfl , mat, wmat(iwMat1) )
c           Load overlap (packed) into lower triangle complex matrix
            call C16LOW( norb, wmat(iwMat1), wmat(iwMat2) )
          endif
c         Distribute complex overlap over processors
          if( norb*norb .le. maxblockz )then
c           The matrix is small enough that pzgemr2d can swallow it whole.    
            call PZGEMR2D( norb,norb,
     $           wmat(iwMat2), 1,1, DESCWORK,
     $           wmat(iwOvlp), 1,1, DESCA,  CONTEXT )
          else
c           Distribute the matrix in (column) stripes that are smaller than 2GB each.
c           Should be good as long the distributed matrix is no larger than 2GB on any processor.
            nblock = maxblockz / norb
            nleft = norb - nblock
            noffset = 1
            moffset = 0
  152       continue
            if( noffset.le.norb )then
              call PZGEMR2D( norb, nblock,
     $             wmat(iwMat2+moffset), 1,1, DESCWORK,
     $             wmat(iwOvlp), 1, noffset, DESCA,  CONTEXT )
              noffset = noffset + nblock
              moffset = moffset + 2*norb*nblock
              if( nleft.ge.nblock )then
c               Have enough to do another full block.
                nleft = nleft - nblock
              else
c               Have less than a full block to send.
                nblock = nleft
                nleft = 0
              endif
              goto 152
            endif
          endif
c
c         Ham is in wmat(iwHmat), overlap in wmat(iwOvlp)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                  Solve Schroedinger equation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
          rmone = 2*PDLAMCH( CONTEXT, 'S' )
c
          call PZHEGVX( 1,jobz, 'I', 'L',
     $     norb, wmat(iwHmat),1,1,DESCA, wmat(iwOvlp),1,1,DESCA,
c -->            H(2,n,n)-i              S(2,n,n)-i
     $     -100d0,-20d0, 1,jstate, rmone, nstout,nvout,
     $     wksml,
c <--      eigval(1,k,ispin),
     $     rmone, wmat(iwVecs),1,1,DESCA,
     $     wmat(iwZwork),lzwork, wmat(iwRwork),lrwork,
     $     wmat(iwIwork),liwork, wmat(iwiwk1),
     $     icluster,gapp, INFO )
c  -->     int(6n)               int(n)
c
c         Check if something went wrong in parallel solver
          if( INFO .ne. 0 )then
            call FLGETIERR( IERRFL )
            write(IERRFL,*) 'PSCHROED/PZHEGVX error, INFO=',INFO
            call STOPXERR( 'PSCHROED/PZHEGVX solver error' )
          endif
c
c         Put eigenvalues into the right place
          call DCOPY( nstate, wksml,1, eigval(1,k,ispin),1 )
c
c         If no eigvecs, we are done with this Ham:
          if( .not. do_eigvecs ) goto 1000
c
c         Collect (only nstate!) eigenvectors onto master
          if( norb*nstate .le. maxblockz )then
c           The matrix is small enough to retrieve whole.
            call PZGEMR2D( norb, nstate,
     $           wmat(iwVecs), 1,1, DESCA,
     $           wmat(iwVeco), 1,1, DESCWORK,  CONTEXT )
          else
c           Retrieve the matrix in (column) stripes that are smaller than 2GB each.
c           Should be good as long the distributed matrix is no larger than 2GB on any processor.
            nblock = maxblockz / norb
            nleft = nstate - nblock
            noffset = 1
            moffset = 0
 1502       if( noffset.le.nstate )then
              call PZGEMR2D( norb, nblock,
     $             wmat(iwVecs+moffset), 1,1, DESCA,
     $             wmat(iwVeco), 1, noffset, DESCWORK,  CONTEXT )
              noffset = noffset + nblock
              moffset = moffset + 2*norb*nblock
              if( nleft.ge.nblock )then
c               Have enough to do another full block.
                nleft = nleft - nblock
              else
c               Have less than a full block to send.
                nblock = nleft
                nleft = 0
              endif
              goto 1502
            endif
          endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c              Record eigenvalues and eigenvectors
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c         Save eigenvectors to ivecfl for later
          if( do_eigvecs .and. ( iproc .eq. masterp ) )then
            call WRITBIG( ivecfl, 2*norb*nstate, wmat(iwVeco) )
          endif
c
 1000   continue
  100 continue
c
c Non-solver procs return here:
c
 2000 continue
c
c If some processors not used in the ScaLAPACK solver context ...
      if( nprocs .ne. npsolver )then
c        ... also need to broadcast eigvals to those not having them
        call MPI_Bcast( eigval(1,1,1), nstate*nk*nspin,
     $   MPI_DOUBLE_PRECISION,
     $   masterp, icomm, ierr )
        if( ierr.ne.0 ) call STOPXERR( 'eigsolv - Bcast eigval' )
      endif
c
      if( myrow .ne. -1 ) call BLACS_GRIDEXIT( CONTEXT )
      call MPBARRIER( icomm )
c      call MPI_ABORT(icomm, ierrors, istatus)
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
