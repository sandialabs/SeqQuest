c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PRODHAB
c
c
      subroutine PRODHAB( a,b,nd,n,w )
c---------------------------------------------------------------
c Purpose: multiply h.c.(a) * b, where a is upper triangular
c      output result in b, leaving a unaffected
c
c Revision history:
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  10Dec97-PAS/2.21: Changed cplx to real; flip loop order (2X)
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  a(2,nd,nd),b(2,nd,nd), w(2,nd)
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      do 50 i=1,n
        do  j=1,n
c         w(j) = (0.,0.)
          w(1,j) = zero
          w(2,j) = zero
        enddo
c loop order gave up to 2X improvement for bigger n
        do  j=1,n
          do  k=1,j
c           w(j) = w(j) + CONJG(a(k,j))*b(k,i)
c parens needed for precision, and is a bit faster, too
            w(1,j) = w(1,j) + ( a(1,k,j)*b(1,k,i)
     $                        + a(2,k,j)*b(2,k,i) )
            w(2,j) = w(2,j) + ( a(1,k,j)*b(2,k,i)
     $                        - a(2,k,j)*b(1,k,i) )
          enddo
        enddo
        do  j=1,n
c         b(j,i) = w(j)
          b(1,j,i) = w(1,j)
          b(2,j,i) = w(2,j)
        enddo
   50 continue
c
      RETURN
      END
