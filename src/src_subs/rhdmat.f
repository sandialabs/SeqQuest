c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHDMAT
c
c
      subroutine RHDMAT( idmatfl, ncplx, norb,
     $ nk,nk0,nkloc,
     $ dmat,             dmatrd )
c      R4(nc*nkloc*mat)  R8(mat)
c---------------------------------------------------------------
c Purpose: read eigenvectors and sort them into an order most
c          convenient for use in the density routines in real.
c
c Written: Peter A. Schultz, 23-July-2001, for v2.49 (spin)
c
c Revision history:
c  26Nov08-PAS/2.63: k-parallel dmat distribution
c   6Jul06-APT/    : merged task parallel and image parallel
c  30Jul05-PAS/2.59: master-only matrix files
c---------------------------------------------------------------
c
c We are only doing triangles of density matrix, hence we need
c need to double the off-diagonals.
c One might think that transposing the dmat's here would lead
c to better memory access in dmat in the box density routines.
c However, tests reversing i-j result in slower execution, at
c least in tests on one DEC Alpha.  Hence, leave it. - Aug01-PAS
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  dmat
      REAL  r4zero, r4two
c
c Output array:
cpas: could pack the dmat here by 2X
      DIMENSION  dmat(ncplx*nkloc,norb,norb)
c Scratch array:
      DIMENSION  dmatrd(norb,norb)
c
c Local declarations
      DATA  r4zero,r4two / 0.e0,2.e0 /
      DATA  two / 2.d0 /
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
        call KPMASTER( nk,nprocl, ik,kmasterpr )
        if( nkloc .eq. nk ) kmasterpr = masterl
c
        if( iprocl .eq. masterl )then
c         Image/local master reads dmat, and send to correct k-master proc
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
        if( ncplx.eq.2 )then
c         Complex density matrix
          kr = 2*(ikloc-1) + 1
          ki = kr + 1
          do  j=1,norb
            do  i=1,j-1
              dmat(kr,i,j) =  two*dmatrd(j,i)
              dmat(ki,i,j) =  two*dmatrd(i,j)
            enddo
            dmat(kr,j,j) =  dmatrd(j,j)
            dmat(ki,j,j) =  r4zero
            do  i=j+1,norb
              dmat(kr,i,j) =  two*dmatrd(i,j)
              dmat(ki,i,j) = -two*dmatrd(j,i)
            enddo
          enddo
        else
c         Real density matrix
          do  j=1,norb
            do  i=1,norb
              dmat(1,i,j) = two*dmatrd(i,j)
            enddo
            dmat(1,j,j) = dmat(1,j,j) / r4two
          enddo
        endif
c
  100 continue
c
C      if( nkloc .ne. nk .and. nkloc .gt. 0 .and. nproc_k .gt. 1 )then
      if( nproc_k .gt. 1 .and. nkloc .gt. 0 )then
c       Send around each k-group its r4 dmat
        lenmsg = ncplx*nkloc*norb*norb
        call MPBCAST4( master_k, lenmsg, dmat, icomm_k )
      endif
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
