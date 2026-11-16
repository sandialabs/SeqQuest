c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ICOPYOFF
c
c
      subroutine ICOPYOFF( n, i1vec,n1off, i2vec,n2off )
c---------------------------------------------------------------
c Purpose: copy integer vector, with offsets
c
c Written: Peter A. Schultz, 20-February-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  i1vec(*), i2vec(*)
c
c >>>> EXECUTABLE CODE:
c
      do  i=1,n
        i2vec(n2off+i) = i1vec(n1off+i)
      enddo
c
c    That's all Folks!
c
      RETURN
      END
