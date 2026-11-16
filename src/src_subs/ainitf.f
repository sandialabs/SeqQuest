c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> AINITF
c
c
      subroutine AINITF( a,grada )
c---------------------------------------------------------------
c Purpose: initialize a() for non-local force matrix element
c
c Revision history:
c   7Mar02-PAS/2.52: ang constants extracted
c  18Nov99-PAS/2.42: cosmetic cleanup
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  a(9,9),grada(3,9,9)
c
c local declarations:
      DIMENSION  vnorm(5)
      DATA  zero,two,four / 0.d0,2.d0,4.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Get angular normalizations:
      call VNORMS( vnorm )
c
      do  j=1,9
        do  i=1,9
          a(i,j) = zero
          grada(1,i,j) = zero
          grada(2,i,j) = zero
          grada(3,i,j) = zero
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
      grada(1,1,2) = -vnorm(1)
      grada(2,1,3) = -vnorm(1)
      grada(3,1,4) = -vnorm(1)
      grada(2,2,5) =  vnorm(2)
      grada(1,3,5) =  vnorm(2)
      grada(3,3,6) =  vnorm(2)
      grada(2,4,6) = -vnorm(2)
      grada(3,2,7) =  vnorm(2)
      grada(1,4,7) = -vnorm(2)
      grada(1,2,8) =  two*vnorm(2)
      grada(2,3,8) = -two*vnorm(2)
      grada(1,2,9) = -two*vnorm(2)
      grada(2,3,9) = -two*vnorm(2)
      grada(3,4,9) = -four*vnorm(2)
c
c    That's all Folks!
c
      RETURN
      END
