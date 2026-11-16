c 
c***********************************************************************
c Purpose: Routines to write/read all k-resolved orbital matrices
c***********************************************************************
c
c For single-processor i/o:
c  H0WRIT - write full set of matrixes from processor
c  H0READ - read full set of matrices onto processor
c For k-parallel distrobuted matrices from master i/o:
c  H0WRITKP - write full set of matrices, from k-parallel matrices on k-mastyers
c  H0READKP - read full set of matrices, onto k-parallel matrices on k-mastyers
c  H0WRITKPDIST - write full set of matrices, from k-parallel fully distributed matrices
c  H0READKPDIST - read full set of matrices, onto k-parallel fully distributed matrices
c  KPDMATRD - read full set of matrices, onto k-parallel matrices on k-masters
c
c***********************************************************************
c Notes:
c  KPDMATRD looks redundant with H0READKP, should merge at some point
c  There is also a "KPDMATRD2" embedded in vslofrc, where the matrix needs the extra
c    reorganization into a full complex (rather than compact real/imaginary triangles).
c    The standard routines here cannot be used for vslofrc, but their existence
c    is noted here to prevent perhaps creating another redundancy in the future.
c***********************************************************************
c Revision history:
c   1Mat12-ACP&PAS/2.63: multi-matrix reads purged, k-parallel reads added
c***********************************************************************
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> H0WRIT
c
c
      subroutine H0WRIT( ifile, mat,nk, vij1 )
c---------------------------------------------------------------
c Purpose: write out a single k-resolved matrix to a disk file
c          The k-resolved matrix is resident on a single process
c
c Written: Peter A. Schultz, 3-October-2002
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Input array:
      DIMENSION  vij1(mat,nk)
c
c >>>> EXECUTABLE CODE:
c
      do  k=1,nk
        call WRITBIG( ifile , mat, vij1(1,k) )
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> H0READ
c
c
      subroutine H0READ( ifile, mat,nk, vij1 )
c---------------------------------------------------------------
c Purpose: read out full k-resolved matrix from file
c          The k-resolved matrix is resident on a single process
c
c Written: Peter A. Schultz, 13-August-2005
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Output array:
      DIMENSION  vij1(mat,nk)
c
c >>>> EXECUTABLE CODE:
c
      do  k=1,nk
        call READBIG( ifile , mat, vij1(1,k) )
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> H0WRITKP
c
c
      subroutine H0WRITKP( ifile, mat, nk, vij1 , vbuffer )
c---------------------------------------------------------------
c Purpose: Write k-resolved, k-distributed matrix
c          resident on k-point masters to master disk file.
c
c Written: Andrew C. Pineda, Mov. 2009-Mar-2010, for 2.63
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Input array:
      DIMENSION  vij1(mat,*)
c Scratch array:
      DIMENSION  vbuffer(*)
c -->            vbuffer(mat)
c
c >>>> EXECUTABLE CODE:
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )

      call MPNODES_K( nprocs_k )
      call MPNODE_K( iproc_k )
      call MPNODE0_K( master_k )
c
c  If processor is not image master or k-point master, skip
      if( iproc .ne. master .and. iproc_k .ne. master_k ) goto 999
c
      ikloc = 1
      do  500 k=1,nk
c       Which processor is k-master, and has matrix for this k-point?
        call KPMASTER( nk, nprocs, k, kmasterpr )

c       If running task-parallel/serial
        if( nprocs .eq. nprocs_k ) kmasterpr = master
c
        if( iproc .eq. master )then
          if( iproc .eq. kmasterpr )then
c           The matrix is local to master, write directly
            call WRITBIG( ifile , mat, vij1(1,ikloc) )
            ikloc = ikloc + 1
          else
c           The master receives matrix from k-master that owns it, using tickets
            call MPSENDTICKET( kmasterpr, icomm )
            call MPRECVR8( kmasterpr, mat, vbuffer, icomm ) 
            call WRITBIG( ifile , mat, vbuffer )
          endif
        else
          if( iproc .eq. kmasterpr )then
c           This k-master sends its matrix to master, waiting for ticker
            call MPRECVTICKET( master, icomm )
            call MPSENDR8( master, mat, vij1(1,ikloc), icomm )
            ikloc = ikloc + 1
          endif
        endif
  500 continue
c
  999 continue
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> H0READKP
c
c
      subroutine H0READKP( ifile, mat, nk, vij1, vbuffer )
c---------------------------------------------------------------
c Purpose: Read k-resolved matrix from master disk file,
c          and distribute to k-masters.
c
c Written: Andrew C. Pineda, Nov 2009-Mar2010, for 2.63
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Input array:
      DIMENSION  vij1(mat,*)
c Scratch array:
      DIMENSION  vbuffer(*)
c -->            vbuffer(mat)
c
c >>>> EXECUTABLE CODE:
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
      call MPNODES_K( nprocs_k )
      call MPNODE_K( iproc_k )
      call MPNODE0_K( master_k )
c
c If processor is not image master or k-point master, skip
      if( iproc .ne. master .and. iproc_k .ne. master_k ) goto 999
c
      ikloc = 1
      do 500 k=1,nk
C       Which processor is k-master, and needs matrix for this k-point
        call KPMASTER( nk, nprocs, k, kmasterpr )

c       If running rask-parallel/serial
        if( nprocs .eq. nprocs_k ) kmasterpr = master
c
        if( iproc .eq. master )then
          if( iproc .eq. kmasterpr ) then
c           This matrix belongs on master, read in directly
            call READBIG( ifile, mat, vij1(1,ikloc) )
            ikloc = ikloc + 1
          else
c           The master sends this matrix to k-master that owns it
            call READBIG( ifile, mat, vbuffer )
            call MPSENDTICKET( kmasterpr, icomm )
            call MPSENDR8( kmasterpr, mat, vbuffer, icomm )
          endif
        else
          if( iproc .eq. kmasterpr )then
c           THe k-master receives the matrix it owns from master
            call MPRECVTICKET( master, icomm )
            call MPRECVR8( master, mat, vij1(1,ikloc), icomm )
            ikloc = ikloc + 1
          endif
        endif
  500 continue
c
  999 continue
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> H0WRITKPDIST
c
c
      subroutine H0WRITKPDIST( ifile, mat, nk, vij1 , vbuffer )
c---------------------------------------------------------------
c Purpose: Write k-resolved,  k-distr=ibuted matrix,
c          distributed over each k-parallel group, to master disk file.
c
c Written: Andrew C. Pineda, 13 Nov 2009, for 2.63
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input array:
      DIMENSION  vij1(mat,*)
c Scratch array:
      DIMENSION  vbuffer(*)
c -->            vbuffer(mat)
c
c >>>> EXECUTABLE CODE:
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
      ikloc = 1
      do 500 k=1,nk
c       What range of processes is matrix(k) distributed over?
        call KPWHERE( nk, k, kproc_min, kproc_max )
c
        do 400  kproc=kproc_min,kproc_max, 1
          if( iproc .eq. master )then
            if( iproc .eq. kproc )then
c             This part of matrix(k) is local to master, write directly
              call WRITBIG( ifile , mat, vij1(1,ikloc) )
              ikloc = ikloc + 1
            else
c             Receive part of matrix(k) owned elsewhere, then write
              call MPSENDTICKET( kproc, icomm )
              call MPRECVR8( kproc, mat, vbuffer, icomm ) 
              call WRITBIG( ifile, mat, vbuffer )
            endif
          else
            if( iproc .eq. kproc )then
c             THis processor sends its part of matrix(k) to master to be written
              call MPRECVTICKET( master, icomm )
              call MPSENDR8( master, mat, vij1(1,ikloc), icomm )
              ikloc = ikloc + 1
            endif
          endif
  400   continue
c
  500 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> H0READKPDIST
c
c
      subroutine H0READKPDIST( ifile, mat, nk, vij1, vbuffer )
c---------------------------------------------------------------
c Purpose: Read k-resolved, k-distr=ibuted matrices from master file,
c          and distribute over all processes in each k-parallel group.
c
c Written: Andrew C. Pineda, 13 Nov 2009, for 2.63
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Output array:
      DIMENSION  vij1(mat,*)
c Scratch array:
      DIMENSION  vbuffer(*)
c -->            vbuffer(mat)
c
c >>>> EXECUTABLE CODE:
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )

      ikloc = 1
      do 500  k=1,nk
c       What range of processes is matrix(k) distrobuted over?.
        call KPWHERE( nk, k, kproc_min, kproc_max )
c
        do 400 kproc=kproc_min,kproc_max,1
          if( iproc .eq. master )then
            if( iproc .eq. kproc )then
c             Data is local to master, read in directly
              call READBIG( ifile, mat, vij1(1,ikloc) )
              ikloc = ikloc + 1
            else
c             Send data to processoer that owns it
              call READBIG( ifile, mat, vbuffer )
              call MPSENDTICKET( kproc, icomm )
              call MPSENDR8( kproc, mat, vbuffer, icomm )
            endif
          else
            if( iproc .eq. kproc )then
c             Receive data fr
              call MPRECVTICKET( master, icomm ) 
              call MPRECVR8( master, mat, vij1(1,ikloc), icomm )
              ikloc = ikloc + 1
            endif
          endif
  400   continue
c
  500 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPDMATRD
c
c
      subroutine KPDMATRD( idmatfl, norb, ncplx, nk,nk0,nkloc,
     $ dmatrd,      dmatloc )                   
c ---> R8(mat)      R8(mat,nkloc)
c---------------------------------------------------------------
c Purpose: pack complex dmat from r8 into k-distributed-c8 for atdefct
c
c Written: Peter A. Schultz, 24-January-2009, for v2.62 (k-parallel)
c
c Revision history:
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c Input:
      INTEGER    idmatfl
      INTEGER    norb, ncplx,nk, nkloc,nk0
c Output:
      DOUBLE PRECISION  dmatloc
      DIMENSION  dmatloc(norb*norb,*)
c Scratch:
      DOUBLE PRECISION  dmatrd
      DIMENSION  dmatrd(*)
c
c Local declarations:
      INTEGER    nprocl ,iprocl ,masterl ,icomml
      INTEGER    nproc_k,iproc_k,master_k,icomm_k
      INTEGER    matrd,lenmsg, ik,ikloc, kmasterpr
c
c >>>> EXECUTABLE CODE:
c
c  Get Local MP info
      call MPNODES( nprocl )
      call MPNODE( iprocl )
      call MPNODE0( masterl )
      call MPCOMM( icomml )
c           
      call MPNODES_K( nproc_k )
      call MPNODE_K( iproc_k )
      call MPNODE0_K( master_k )
      call MPCOMM_K( icomm_k )
c       
      matrd = norb*norb
c
      ikloc = 0
      do 100 ik=1,nk
c
        lenmsg = matrd 
        call KPMASTER( nk, nprocl, ik, kmasterpr )
        if( nkloc .eq. nk ) kmasterpr = masterl
c       
        if( iprocl .eq. masterl )then
c         Image/local master reads dmat, and sends it to correct k-master proc
          call READBIG( idmatfl, matrd, dmatrd )
          if( kmasterpr .ne. masterl )then
c           This k-point goes off to kmasterpr
            call MPSENDR8( kmasterpr, lenmsg, dmatrd, icomml )
            goto 100
          endif 
        elseif( iprocl .eq. kmasterpr )then
c         Proc is master of this k-point, get dmatrd from masterl
          call MPRECVR8( masterl, lenmsg, dmatrd, icomml )
        else
c         This proc is not leading this k-point
          goto 100
        endif 
        ikloc = ikloc + 1
c
        if( nk.gt.1 )then
c         Only need to cache this if there is more than 1 k-point
c         (NB: _not_ nkloc.gt.1, master needs this space to store&read
          call DCOPY( matrd, dmatrd,1, dmatloc(1,ikloc),1 )
        endif
c
  100 continue
c
C      if( nproc_k .gt. 1 .and.  nkloc .gt. 0 )then
Cc       Send around each k-group its r8 dmat
C        lenmsg = nkloc*norb*norb
Cc       Actually, this bcast is moot, as only k-masters work in atdefect
C        if( nk.eq.1 )then
C          call MPBCAST8( master_k, lenmsg, dmatrd, icomm_k )
C        else
C          call MPBCAST8( master_k, lenmsg, dmatloc, icomm_k )
C        endif
C      endif
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
