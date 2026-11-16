c     
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PSCHROEDKP
c
c
      subroutine PSCHROEDKP( do_eigvecs,
     $     idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $     norb,nstate,nk, eigval, wksml,
     $     nspin, icluster, gapp,
     $     wmat, masterlist )
c---------------------------------------------------------------
c Purpose: solve Schroedinger equation given input ham and ovlp,
c     fully k-parallel-parallel solve. 
c     
c Written: Andrew C. Pineda,  9-February-2010 for v2.64,
c          (adapted from "pschroed" routine)
c     
c Revision history:
c  16May12-ACP/2.64: Removed ENTRY point.
c  2May12-ACP/2.64: Added alternate entry point to support non-SCF portions of code
c  17Apr12-ACP/2.64: Work around ScaLAPACK 2GB limit
c---------------------------------------------------------------
c     
c  Notes:
c     Each k-point can have multiple processes associated with it.
c     Is so, SCALAPACK is used instead of LAPACK.
c     The master process, with access to all the data on disk, 
c     distributes ham and ovlp for eachk-point to a k-point master.
c     Each k-master leads a parallel solvewithin its MPI context,
c     BLACS context derived from the k-point communicator). 
c
c     The call parameter idosolv was used to control which version of the
c     LAPACK routine ZHEGV was used.  Here, we always use PZHEGVX.
c     parameter will be used only to control whether the number of
c     eigenvectors/values returned is all (norb) or just nstate.
c     Note: this means that SCF treatment of linear dependence
c       (idsolv=3) is not possible with this routine as it stands.
c     DMC 8Apr03
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c     MAXBLOCK= 2**31/sizeof(DOUBLE COMPLEX) = 2**31/16.
      PARAMETER  ( maxblockz = 67 108 864 )
c
ctp: Pull in the mpi stuff
      INCLUDE  'mpif.h'
c
c Input:
      LOGICAL    do_eigvecs
c Output:
      DIMENSION  eigval(nstate,nk,nspin)
c Scratch:
      DIMENSION  wksml(*)
      DIMENSION  wmat(*)
c Scratch/parallel:
      INTEGER    icluster(*)
      DIMENSION  gapp(*)
      INTEGER    masterlist(*)
c
c Local declarations:
      CHARACTER  jobz*1
ctp: Local for ScaLAPACK
      PARAMETER  (NB = 32)
      INTEGER    ANB, SQNPC
      INTEGER    CONTEXT_K
      INTEGER    DESCA(9), DESCWORK(9)
      INTEGER    procindex
c
      LOGICAL    USE_SCALAPACK_FORMULAS
      LOGICAL    OPTIMAL_EIGENVECTORS
      DATA  USE_SCALAPACK_FORMULAS / .true. /
      DATA  OPTIMAL_EIGENVECTORS   / .true. /
      LOGICAL    DEBUG_MEM
      DATA       DEBUG_MEM / .false. /
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
cxxx: declare dev-version of pschroedkp:
      if( IWR.ge.0 ) write(IWR,*)
     $ 'PSCHROEDKP: k-parallel-parellel solves (opt mem)'
c
      mat = norb*norb
c
c     Get Local MP info
c
c     Image level
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
c
c     K-parallel level
      call MPNODES_K( nprocs_k )
      call MPNODE_K( iproc_k )
      call MPNODE0_K( master_k )
      call MPCOMM_K( icomm_k )
c
c     Build a list of the ranks of the kpmasters within the icomm
c     communicator and get a count of kpmasters. The returned list of
c     ranks is in increasing order.
c
      call KPMASTERLIST( nprocs, iproc, masterp, icomm,
     $     nprocs_k, iproc_k, master_k, icomm_k,
     $     masterlist, num_masters )
c
      if( IWR.ge.0 )then
        write(IWR,9131) 'number of kpp masterss=',num_masters
         if( DEBUG_MEM ) write(IWR,9132)
     $     ' kpp masterlist ranks=',(masterlist(i),i=1, num_masters)
      endif
 9131 format( 4x,a,1x,i4)
 9132 format( 4x,a,(16i5))
c
c     Broadcast rank in icomm communicator of kp-master to processes in
c     icomm_k and store in my_kpmaster
c
      lenmsg = 1
      my_kpmaster = iproc
      call MPBCASTI( master_k, lenmsg, my_kpmaster, icomm_k )
c
      if( nprocs_k.eq.1 )then
c       SINGLE PROCESSOR SOLVE - Use LAPACK
c
        matdist = mat
C        lwzwork  = 2*MAX( 2*norb, mat )
        lwzwork  = MAX( 4*norb, mat )
        lzwork  = lwzwork / 2
        lwork  = lzwork
        lrwork  = 7*norb
        liwork  = 5*norb
c
c       Cloase single-processor solve work sizes
      else
c       PRALLEL SOLVER - Use SCALAPACK
c
c      Initialize ScaLAPACK Grid for k-point to be handled by this processor.
c      Define processor grid dimensions of largest nearly square
c       rectangular array of processors within icomm_k.
c       npsolver_k is the number of processors in the grid (npsolver_k.le.nprocs_k).
c
        call EIGPROCS( nprocs_k, nprow_k, npcol_k, npsolver_k )
        if( IWR.gt.0 ) write(IWR,'(a,4i6)') ' procs, rowsxcols=solver:',
     $                nprocs_k,  nprow_k,npcol_k, npsolver_k
c
c     Create BLACS context using system context from 
c     k-point communicator.
c
        CONTEXT_K = icomm_k
        call BLACS_GRIDINIT( CONTEXT_K, 'R', nprow_k, npcol_k )
        call BLACS_GRIDINFO( CONTEXT_K,
     $       nprow_k, npcol_k,
     $       myrow_k, mycol_k )
c
c     myrow_k = -1 if iproc > npsolver_k?
c
c     Initialize matrix descriptions, for procs in solver contexts
c
        if( myrow_k .ne. -1 )then
c         Process in solver context, set up matrix descriptors
c
c         First, set up descriptor for distributed matrix
c
c         Need to limit lead dimension in local space
c          nblkrow = ceiling(DBLE(norb)/DBLE(NB)) = # of blocks in a row
c          npblkrow = ceiling(DBLE(nblkrow)/DBLE(nprow_k))
c          Lead dimension LLDA is row bloacks times blocking factor
          nblkrow = (norb+NB-1) / NB
          npblkrow = (nblkrow+nprow_k-1) / nprow_k
          LLDA = npblkrow*NB
          call DESCINIT( DESCA,  norb,norb, NB,NB,  0, 0,
     $         CONTEXT_K, LLDA, INFO )
          if( INFO.ne.0 )
     $      call STOPXERR( 'pschroedkp - DESCINIT/DESCA' )
c
c        Set up the descriptor for full local matrix.
c
c        call DESCINIT( DESC, M, N, MB, NB, RSRC, CSRC, ICTXT, MXLLDA, INFO)
c        MB=M, NB=N means only the processor (RSRC,CSRC)=(0,0) has any data.
c
          call DESCINIT( DESCWORK,  norb,norb, norb,norb,  0, 0,
     $         CONTEXT_K, norb, INFO )
          if( INFO.ne.0 )
     $      call STOPXERR(' pschroedkp - DESCINIT/DESCWORK' )
        else
c         This processor is not in solver context, and leaves
          goto 990
        endif
c
        if( USE_SCALAPACK_FORMULAS )then
c         Get EXACT mpi/blacs-aware work space using this context:
          call PSCHWKSZX( CONTEXT_K, npsolver_k, do_eigvecs,
     $     norb, matdist, lwzwork, lzwork, lrwork, liwork )
          lwork = lzwork
        else
c         Use default mpi/blacs-avoiding naive work space:
          call PSCHWKSZ0(            npsolver_k, do_eigvecs,
     $     norb, matdist, lwzwork, lzwork, lrwork, liwork )
          lwork = lzwork
        endif
c
        if( DEBUG_MEM ) write(IWR,'(a,3i12)')'DEV/SCHKPP: works Z,R,I=',
     $                                        lwzwork,lrwork,liwork
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                Do the memory allocations for parallel solves
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     MEM-Solver: d-C16-H, d-C16-S, d-C16-Vec, d-C16-lwork, d-R8-rwork)
c
c     The solver memory allocations have the following requirements:
        nwpreph = 2*matdist + 3*mat
        nwpreps = 4*matdist + 2*mat
        nwprep  = MAX( nwpreph, nwpreps )
        nwsolv  = 6*matdist + lwzwork + lrwork + liwork + norb
c     nwclose = 6*matdist + 2*norb*nstate
c     because vectors are brought back into full norb*norb ...
        nwclose = 6*matdist + 2*norb*norb
        nwpsch  = MAX( nwprep, nwsolv, nwclose )
c
        nwleft = nwpsch - nwsolv
        if( DEBUG_MEM )then
          write(IWR,*) 'PSCHROEDKP: spare solver mem=',nwleft
        endif
c
c     end of  if (nprocs_k.gt.1) then
      endif
c
      iwHmat  = 1
      iwOvlp  = iwHmat + 2*matdist
      iwVecs  = iwOvlp + 2*matdist
      iwZwork = iwVecs + 2*matdist
c
      iwRwork = iwZwork + lwzwork
      iwIwork = iwRwork + lrwork
      iwiwk1  = iwIwork + liwork
      iwlast  = iwwk1   + norb
      iwEnd   = iwlast  - 1
c
c     MEM-Prep: (d-C16-H, d-C16-S=full-R8-mat, full-C16-mat)
c
      iwMat1 = iwOvlp
c
c     Put full-C16-mat after full-R8-mat AND after dist-S
c
      iwMat2 = MAX( iwMat1 + mat, iwVecs )
c
c     MEM-After: (d-C16-H, d-C16-S, d-C16-Vec, full-C16-Vec)
c     Put full output vectors after dist-vectors:
c
      iwVeco = iwZwork
c
c     Set files to correct locations
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
        if( do_eigvecs )then
c         Compute eigenvalues and eigenvectors
          jobz = 'V'
        else
c         Compute eigenvalues only
          jobz = 'N'
        endif
c
c     Bloch vector loop: individual k point problems are distributed to
C     blocks of processes controled by the k-point masters.
c
        do 1000 kblk = 1, nk, num_masters
          my_k = -1
          do 1010 kk = 1, num_masters
            k = kk + kblk - 1  
c     k out of range, everyone skips to end of loop if they are not the image master
c     which must remain to do I/O.
            if( k .gt. nk )then
              if( (my_k.eq.-1) .and. (iproc.eq.masterp) )then
                my_k = 0   
              endif
              goto 1011
            endif
c     Since masterlist is in increasing order, we ensure that masterp does the last k-point
            procindex = 1 + MOD( kk, num_masters )
c     Select the k-point master for this k.
            kmasterpr = masterlist(procindex)
            if( my_kpmaster .eq. kmasterpr )then
c             Set parameters for the group of processes working on this k-point.
              my_k = k
              if( idosolv .eq. 1 )then
                jstate = norb
              elseif( idosolv .eq. 2 )then
                jstate = nstate
              else
                call STOPXERR( 'PSCHROEDKP: invalid solver option' )
              endif
            endif
c
c     Form H0+dH=full Hamiltonian(spin,k), and distribute it to the k-point master
c
            if( iproc .eq. masterp )then
c             Retrieve reference Hamiltonian H0 (packed) from iham0fl
              call READBIG( iham0fl, mat, wmat(iwMat2) )
              if( iproc.ne.kmasterpr )then
                lenmsg = mat
                call MPSENDR8( kmasterpr,lenmsg,wmat(iwMat2),icomm )
              endif
c             Retrieve delta-Hamiltonian dH (packed) from idhamfl
              call READBIG( idhamfl, mat, wmat(iwMat1) )
              if( iproc.ne.kmasterpr )then
                lenmsg = mat
                call MPSENDR8( kmasterpr,lenmsg,wmat(iwMat1),icomm )
              endif
            endif
c
c     Receive full Hamiltonian on the k-point master at wmat(iwMat2).
            
            if( iproc .eq. kmasterpr )then
              if( iproc .ne. masterp )then
                lenmsg = mat
                call MPRECVR8( masterp,lenmsg,wmat(iwMat2),icomm )
                lenmsg = mat
                call MPRECVR8( masterp,lenmsg,wmat(iwMat1),icomm )
              endif
c             Combine H0 and dH to get full Hamiltonian (packed)
              do  ij=1,mat
                wmat(iwMat1-1+ij) =
     $               wmat(iwMat1-1+ij) + wmat(iwMat2-1+ij)
              end do
c             Load Ham (packed real) into lower triangle of complex matrix
              call C16LOW( norb, wmat(iwMat1), wmat(iwMat2) )
c----------------------------> n   real*8 nxn     complex*16 nxn
c---------------------------->                  (lower triangular)
c
c             Close if(iproc.eq.kmasterpr) 
            endif
c
c     Distribute complex Ham over the subset of processes belonging to CONTEXT_K using SCALAPACK.
c     If I have done this right this is the group of processes headed up by kpmasterpr.
c
c     iwMat2 --> full matrix on kmasterpr
c     iwHmat --> distributed matrix within icomm_k
c
            if( my_kpmaster .eq. kmasterpr )then
              if( nprocs_k.eq.1 )then
c               Copy iwMat2 to iwHmat
                call DCOPY( 2*mat, wmat(iwMat2), 1,
     $               wmat(iwHmat), 1)
              else
c     Distribute the matrix if iproc is part of the correct CONTEXT_K/icomm_k.
c     All processes in this CONTEXT_K/icomm_k must participate in this call.
                if( norb*norb.le.maxblockz )then
c                 Matrix small enough that pzgemr2d can swallow it whole.    
                  call PZGEMR2D( norb,norb,
     $                 wmat(iwMat2), 1,1, DESCWORK,
     $                 wmat(iwHmat), 1,1, DESCA,
     $                 CONTEXT_K )
                else
c                 Distribute the matrix in (column) stripes smaller than 2GB each.
c                 Good if distributed matrix is no larger than 2GB on any processor.
                  nblock = maxblockz / norb
                  nleft = norb - nblock
                  noffset = 1
                  moffset = 1
  151             continue
                  if( noffset.le.norb )then
                    call PZGEMR2D( norb,nblock,
     $                   wmat(iwMat2-1+moffset), 1, 1, DESCWORK,
     $                   wmat(iwHmat), 1, noffset, DESCA,
     $                   CONTEXT_K )
                    noffset = noffset + nblock
                    moffset = moffset + 2*norb*nblock
                    if( nleft.ge.nblock )then
c                     Have enough to do another full block.
                      nleft = nleft - nblock
                    else
c                     Have less than a full block to send.
                      nblock = nleft
                      nleft = 0
                    endif
                    goto 151
                  endif
                endif
              endif
            endif
c
c     Retrieve overlap matrix on masterp and send it to k-point master
c
            if( iproc .eq. masterp )then
c             Retrieve overlap matrix (packed) from iovlpfl
              call READBIG( iovlpfl , mat, wmat(iwMat1) )
              if( iproc .ne. kmasterpr )then
                lenmsg = mat
                call MPSENDR8( kmasterpr,lenmsg,wmat(iwMat1),icomm )
              endif
            endif
c
c     Receive complex overlap on k-point master
c
            if( iproc .eq. kmasterpr )then
              if( iproc .ne. masterp )then
                lenmsg = mat
                call MPRECVR8( masterp,lenmsg,wmat(iwMat1),icomm )
              endif
c             Load overlap (packed) into lower triangle complex matrix
              call C16LOW( norb, wmat(iwMat1), wmat(iwMat2) )
c----------------------------n   real*8 nxn    complex*16 nxn
c----------------------------(lower triangular)
            endif
c
c     Distribute complex overlap over processes controlled by k-point master
c
            if( my_kpmaster .eq. kmasterpr )then
              if( nprocs_k.eq.1 )then
c               Copy complex n X n matrix into place.
                call DCOPY( 2*mat, wmat(iwMat2), 1,
     $               wmat(iwOvlp), 1 )
              else
c     Distribute the matrix if iproc is part of the correct CONTEXT_K/icomm_k.
c     All processes in this CONTEXT_K/icomm_k must participate in this call.
                if( norb*norb .le. maxblockz )then
c                 Matrix is small enough that pzgemr2d can swallow it whole.                      
                  call PZGEMR2D( norb,norb,
     $                 wmat(iwMat2), 1,1, DESCWORK,
     $                 wmat(iwOvlp), 1,1, DESCA,
     $                 CONTEXT_K )
                else
c                 Distribute matrix in (column) stripes smaller than 2GB each.
c                 Good if distributed matrix no larger than 2GB on any processor.
                  nblock = maxblockz / norb
                  nleft = norb - nblock
                  noffset = 1
                  moffset = 1
  152             continue
                  if( noffset.le.norb )then
c     We handle the offset to the full matrix here so that ScaLAPACK doesn't have to.
c     NB: This foreshadows what we want to do to distribute the matrix in stripes.
                    call PZGEMR2D( norb,nblock,
     $                   wmat(iwMat2-1+moffset), 1,1, DESCWORK,
     $                   wmat(iwOvlp), 1,noffset, DESCA,
     $                   CONTEXT_K )
                    noffset = noffset + nblock
                    moffset = moffset + 2*norb*nblock
                    if( nleft.ge.nblock )then
c                     Have enough to do another full block.
                      nleft = nleft - nblock
                    else
c                     Have less than a full block to send.
                      nblock = nleft
                      nleft = 0
                    endif
                    goto 152
                  endif
                endif
              endif
            endif
c     
c     End of 1st kk loop (Eigenproblem distribution)
c     
 1010     continue
 1011   continue
c
c     Ham is in wmat(iwHmat), overlap in wmat(iwOvlp)
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                      Solve Schroedinger equation
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Processors in a context do the solve for their assigned k-point
c
c     If I'm part of a CONTEXT_K that is processing a k-point (my_k>=0),
c     call the SCALAPACK solver. Note that there are up to num_masters
c     parallel solves happening here as there are num_masters CONTEXT_Ks. 
c
        if( my_k .ge. 1 )then
c
          if( nprocs_k.eq.1 )then
c           Only one processor assigned to this k-point, use LAPACK
c
            rmone = 2*DLAMCH('S')
c
            call ZHEGVX( 1,   jobz, 'I',  'L', norb,
     $       wmat(iwHmat), norb,  wmat(iwOvlp),norb,
c -->        A             LDA    B            LDB
     $       -100d0, -20d0,    1, jstate,  rmone,
c -->        VL      VU c      IL IU       ABSTOL
     $       nstout,   wksml,
c -->        M         W[N]
     $       wmat(iwVecs), norb,
c -->        Z             LDZ  [Z(LDZ, max(1,M))]
     $       wmat(iwZwork), lwork,
c -->        WORK           LWORK    [WORK(LWORK)]
     $       wmat(iwRwork),  wmat(iwIwork),  wmat(iwiwk1),
c -->        RWORK[7*N]      IWORK[5*N]      IFAIL[N]
     $       INFO )
c
c           If something failed, go do postmortem
            if( INFO .ne. 0 ) goto 1310
c
c           Close serial LAPACK solve
          else
c           Parallel Scalapack solver
c
            if( myrow_k .ne. -1 ) then
c
              rmone = 2*PDLAMCH( CONTEXT_K, 'S' )
c
              call PZHEGVX( 1, jobz, 'I', 'L', norb,
     $             wmat(iwHmat), 1,1,DESCA, 
     $             wmat(iwOvlp), 1,1,DESCA,
     $             -100d0,-20d0, 1,jstate, rmone,
     $             nstout, nvout, wksml, rmone, 
c <--                             eigval(1,k,ispin),
     $             wmat(iwVecs), 1,1,DESCA,
     $             wmat(iwZwork), lwork,
     $             wmat(iwRwork), lrwork,
     $             wmat(iwIwork), liwork,
     $             wmat(iwiwk1),
     $             icluster,  gapp, INFO )
c
c             If something failed, go do postmortem
              if( INFO .ne. 0 ) goto 1320
c
c             Close parallel solver branch
            endif
          endif
        endif
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                      Record eigenvalues and eigenvectors
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Retrieve results from kpmasters
        if( my_k .ge. 0 )then
          do 2010 kk = 1, num_masters
            k = kk + kblk - 1
c           If k out of range, everyone skips to end of loop
            if( k .gt. nk )then
              goto 2011
            endif
            procindex = 1 + MOD( kk,num_masters )
            kmasterpr = masterlist(procindex)
c
c     Put eigenvalues into the right place on masterp.
c
            if( iproc .eq. kmasterpr )then
              if( iproc .eq. masterp )then
                call DCOPY( nstate, wksml, 1,
     $               eigval(1,k,ispin), 1 )
              else
                lenmsg = nstate
                call MPSENDR8( masterp,lenmsg,wksml,icomm )
              endif
            endif
            if( iproc .eq. masterp )then
              if( iproc .ne. kmasterpr )then
                lenmsg = nstate
                call MPRECVR8( kmasterpr,lenmsg,
     $               eigval(1,k,ispin),icomm )
              endif
            endif
c
            if( do_eigvecs )then
c     Optionally retrieve eigenvectors for current k-point from wmat(iwVecs) and
c     put them in wmat(iwVeco) on the k-point master (kmasterpr).
c     We only retrieve the vectors for nstate states (length norb) we need.
c     There should be at least nstate states in wmat(iwVecs), more if
c     jstate=norb.
              if( my_kpmaster .eq. kmasterpr )then
                if( nprocs_k.eq.1 )then
                  call DCOPY( 2*norb*nstate, wmat(iwVecs), 1,
     $                 wmat(iwVeco), 1 )
                else
c                 Matrix small enough that PZGEMR2D won't choke reading it all at once.
                  if( norb*nstate.le.maxblockz )then
                    call PZGEMR2D( norb, nstate,
     $                   wmat(iwVecs), 1, 1, DESCA,
     $                   wmat(iwVeco), 1, 1, DESCWORK,
     $                   CONTEXT_K )
                  else
c                   Retrieve matrix in (column) stripes smaller than 2GB each.
c                   Good if distributed matrix is no larger than 2GB on any processor.
                    nblock = maxblockz / norb
                    nleft = nstate - nblock
                    noffset = 1
                    moffset = 1
  153               continue
                    if( noffset.le.nstate )then
                      call PZGEMR2D( norb, nblock,
     $                     wmat(iwVecs), 1, noffset, DESCA,
     $                     wmat(iwVeco-1+moffset), 1, 1, DESCWORK,
     $                     CONTEXT_K )
                      noffset = noffset + nblock
                      moffset = moffset + 2*norb*nblock
                      if( nleft.ge.nblock )then
c                       Have enough to do another full block.
                        nleft = nleft - nblock
                      else
c                       Have less than a full block to send.
                        nblock = nleft
                        nleft = 0
                      endif
                      goto 153
                    endif
                  endif
                endif
              endif
c
c     Send eigenvectors from kpmaster for this k-point to masterp
c     if kpmaster is not the same as masterp.
c
              if( iproc .eq. kmasterpr )then
                if( iproc .ne. masterp )then
c                 Send only the nstate (complex) vectors we need.
c                 lenmsg is 2*norb*nstate since the wmat data is complex.
                  lenmsg = 2*norb*nstate
                  call MPSENDR8( masterp,lenmsg,wmat(iwVeco),icomm )
                endif
              endif
c
c             On masterp, save the eigenvectors to ivecfl.
c
              if( iproc .eq. masterp )then
c               Receive nstate vectors we need from kmasterpr if needed. 
c               The factor of 2 is because the data is complex.
                if( iproc .ne. kmasterpr )then
                  lenmsg = 2*norb*nstate
                  call MPRECVR8( kmasterpr,lenmsg,wmat(iwVeco),
     $                 icomm )
                endif
c               wmat(iwVeco) on masterp is filled either by PZGEMR2D or MPRECVR8
c               nstate <= jstate so we should always have data.
                isize = 2*norb*nstate
                call WRITBIG( ivecfl, isize, wmat(iwVeco) )
              endif
            endif
c
c           End of 2nd kk loop
 2010     continue
 2011     continue
        endif
c       Close kblk loop
 1000 continue
c     Close ispin loop
  100 continue
c
c Non-solver processors re-enter here:
  990 continue
c
c  Need to broadcast eigvals to everyone
c
      lenmsg = nstate*nk*nspin
      call MPBCAST8( masterp, lenmsg, eigval(1,1,1), icomm )
c
c  Everyone exits their BLACS context
c
      if( nprocs_k.gt.1 )then
        if( myrow_k .ne. -1 ) call BLACS_GRIDEXIT( CONTEXT_K )
      endif
c
      call MPBARRIER( icomm )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                           NORMAL EXIT
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      RETURN
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                     SOLVER ERROR POSTMORTEMS
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c           ***** Lapack/ZHEGVX error postmortem *****
 1310 continue
      call FLGETIERR( IERRFL )
      write(IERRFL,*) '***** ERROR: PSCHROEDKP/ZHEGVX, INFO=',INFO
      if( INFO .lt. 0 )then
        write(IERRFL,*) -INFO,'th argument of ZHEGVX had  illegal value'
        if( INFO .eq. -20 ) write(IERRFL,*) 'lwork = ', lwork
      else
        if( INFO .le. norb )then
          write(IERRFL,*)  INFO,' eigenvectors failed to converge'
        else
          write(IERRFL,*)
     $     'factorization of overlap matrix could not be completed '//
     $     'no eigenvalues or eigenvectors computed'
        endif
      endif
      call FLUSH( IERRFL )
c
      call STOPXERR( 'PSCHROEDKP/ZHEGVX solver error' )
      goto 1399
c
c           ***** Scalapack/PZHEGVX error postmortem *****
 1320 continue
      call FLGETIERR( IERRFL )
      write(IERRFL,*) '***** ERROR: PSCHROEDKP/PZHEGVX, INFO=',INFO
c
      if(MOD( INFO,2 ) .ne. 0 )then
        write(IERRFL,*) 'One or more eigenvectors failed to converge'
        call print_ifail_array( IERRFL, wmat(iwiwk1), norb )
      endif
      if( MOD( INFO/2 ,2 ) .ne. 0 )then
        write(IERRFL,*) 'One or more eigenvectors '//
     $   'could not be reorthogonalized due to insufficient workspace'
        write(IERRFL,*) 'Here is the list of clusters:'
c       The cluster list is a zero terminated list
        ikl = 1
        write(IERRFL,*) 'Cluster #    Starts    Ends   Size    Gap'
 1321   write(IERRFL,*)  ikl, icluster(2*ikl-1),icluster(2*ikl),
     $    icluster(2*ikl)-icluster(2*ikl-1)+1, gapp(ikl)
        if( icluster(2*ikl+1).ne.0 )then
          ikl = ikl + 1
          goto 1321
        endif
        write(IERRFL,*)'Here are the workspace sizes:'
        write(IERRFL,*)'lwork,lrwork,liwork=',lwork,lrwork,liwork
      endif
      if( MOD( INFO/4, 2 ) .ne. 0 )
     $ write(IERRFL,*) 'Space limit prevented PZHEGVX from ' //
     $  'computing all eigenvectors between VL(-100) and VU(-20).'
      if( MOD( INFO/8, 2 ) .ne. 0 )
     $ write(IERRFL,*) 'Scalapack PZSTEBZ failed to compute eigenvalues'
      if( mod( INFO/16, 2 ) .ne. 0 )then
        write(IERRFL,*) 'Overlap matrix is not positive definite'
        write(IERRFL,*) 'norb=',norb
        call print_ifail_array( IERRFL, wmat(iwiwk1), 1 )
      endif
      if( INFO.eq.-28 )
     $  write(IERRFL,*) 'LWORK= ', lwork,', but needs=',wmat(iwZwork)
      if( INFO.eq.-30 )
     $  write(IERRFL,*) 'LRWORK=',lrwork,', but needs=',wmat(iwRwork)
      if( INFO.eq.-32 )
     $  write(IERRFL,*) 'LRWORK=',liwork,', but needs=',wmat(iwIwork)
c
      call STOPXERR( 'PSCHROEDKP/PZHEGVX solver error' )
c
 1399 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> print_ifail_array
c
      subroutine print_ifail_array(iunit, iarray, n)
      INTEGER  iunit
      INTEGER  n
      INTEGER  iarray(n)
      INTEGER  j
      do  j=1,n
        write(iunit,*)'ifail(',j,')=', iarray(j)
      enddo
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
