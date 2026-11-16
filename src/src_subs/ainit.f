c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> AINIT
c
c
      subroutine AINIT( a )
c---------------------------------------------------------------
c Purpose: initialize a-matrix for non-local matrix element
c
c Revision history:
c   7Mar02-PAS/2.52: ang constants extracted
c  18Nov99-PAS/2.42: cosmetic cleanup
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  a(9,9)
c
c local declarations:
      DIMENSION  vnorm(5)
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Get angular normalizations:
      call VNORMS( vnorm )
c
      do  j=1,9
        do  i=1,9
          a(i,j) = zero
        enddo
      enddo
c
      a(1,1) =  vnorm(1)
      a(2,2) = -vnorm(2)
      a(3,3) = -vnorm(2)
      a(4,4) =  vnorm(2)
      a(5,5) =  vnorm(3)
      a(6,6) = -vnorm(3)
      a(7,7) = -vnorm(3)
      a(8,8) =  vnorm(4)
      a(9,9) =  vnorm(5)
c
c    That's all Folks!
c
      RETURN
      END
