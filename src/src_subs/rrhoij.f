c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RRHOIJ
c
c
      subroutine RRHOIJ(doblas3,
     $     idmatfl,iematfl,norb,nstate,
     $     eigval,eigpop,fpop,numpop, wt,
     $     dmat,emat,eigvec,work)
c---------------------------------------------------------------
c Purpose: construct density matrix in orbital basis, real
c
c Written: Rick Muller, June 4, 2007
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
      integer n,i,j,mat,numpop
      integer idmatfl,iematfl,norb,nstate
      double precision eigval(nstate),eigpop(nstate),
     $     fpop(nstate)
      double precision wt,dm,em
      LOGICAL    doblas3
      double precision dmat(norb,norb), emat(norb,norb)
      double precision eigvec(norb,norb), work(norb,norb)
c
      double precision zero,half,one,two
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
c
c >>>> EXECUTABLE CODE:
c
      mat = norb*norb
      if( numpop .lt. 1 )then
c     No occupied levels: dmat=emat=0, and skip construction
         call MKZERO( mat, dmat )
         call MKZERO( mat, emat )
         goto 390
      endif
c
c     Scale level occupations with appropriate sample weights
      do  n=1,numpop
         fpop(n) = wt*eigpop(n)
      enddo
c
      if( doblas3 ) then
c
c >>>>> Closed shell, use Raczkowski blas code:
c
         call DCOPY( norb*numpop, eigvec,1,
     $        work,1 )
         do  i=1,numpop
            call DSCAL( norb, fpop(i), work(1,i),1 )
         enddo
c
c     Make density matrix:
         call DGEMM( 'n','t', norb,norb,numpop,one,
     $        work,norb, eigvec,norb,
     $        zero, dmat,norb )
c
c     Make e-matrix:
         do  i=1,numpop
            call DSCAL( norb, eigval(i), work(1,i),1 )
         enddo
c
         call DGEMM( 'n','t', norb,norb,numpop,one,
     $        work,norb,  eigvec,norb,
     $        zero, emat,norb)
c     
c     Clear out imaginary part of matrices:
         do  j=1,norb
            do  i=1,j-1
               dmat(i,j) = zero
               emat(i,j) = zero
            enddo
         enddo
c
c     Scale dmat and emat by weight function:
         call DSCAL( 2*mat, wt, dmat,1 )
c
        else
c
c >>>>> General occupations code to generate density matrix:
c
c     Transpose the eigenvector matrix:
           do  n=1,numpop
              do  i=1,norb
                 work(n,i) = eigvec(i,n)
              enddo
           enddo
c     
           do  i=1,norb
              do  n=1,numpop
                 eigvec(n,i) = fpop(n)*work(n,i)
              enddo
           enddo
c
           do  j=1,norb
c     Set imaginary triangle to zero
              do  i=1,j-1
                 dmat(i,j) = zero
                 emat(i,j) = zero
              enddo
c     Construct real triangle
              do  i=j,norb
                 dm = zero
                 em = zero
c     Tripling up on this loop (3 i's at a time) gets you another
c     25%, not enough to merit complicating the code.
                 do  n=1,numpop
                    dm = dm + eigvec(n,i)*work(n,j)
                    em = em + eigvec(n,i)*work(n,j)*eigval(n)
                 enddo
                 dmat(i,j) = dm
                 emat(i,j) = em
              enddo
           enddo
c
        endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c               Write out density-matrix, e-matrix:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
 390    continue
c
c       Save density matrix to idmatfl:
        call WRITBIG( idmatfl, mat, dmat )
c
c       Save e-matrix to iematfl:
        call WRITBIG( iematfl, mat, emat )
c
c    That's all, Folks!
c
      RETURN
      END
