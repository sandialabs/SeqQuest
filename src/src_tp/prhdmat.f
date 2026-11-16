c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHDMAT
c
c
      subroutine PRHDMAT( idmatfl, ncplx,norb,nk,
     $ dmat, dmatrd )
c---------------------------------------------------------------
c Purpose: read eigenvectors and sort them into an order most
c          convenient for use in the density routines in real.
c
c Written: Peter A. Schultz, 23-July-2001, for v2.49 (spin)
c
c Revision history:
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
      DIMENSION  dmat(ncplx*nk,norb,norb)
c Scratch array:
      DIMENSION  dmatrd(norb,norb)
c
c Local declarations
      DATA  r4zero,r4two / 0.e0,2.e0 /
      DATA  two / 2.d0 /
c
c >>>> EXECUTABLE CODE:
c
      matrd = norb*norb
c
c  Get Local MP info
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c     
      do 100 k=1,nk
c
        if( iproc.eq.master) call READBIG( idmatfl, matrd, dmatrd )
        if( nprocs.gt.1 ) call MPBCAST8( master, matrd, dmatrd, icomm )
c
        if( ncplx.eq.2 )then
c         Complex density matrix
          kr = 2*(k-1) + 1
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
c    That's all Folks!
c
      RETURN
      END
