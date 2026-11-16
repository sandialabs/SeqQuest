c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GOVLP
c
c
      subroutine GOVLP( sg, l, n1,a1,c1, n2,a2,c2 )
c---------------------------------------------------------------
c Purpose: compute overlap between pair of contract gaussian orbs
c
c Written: Peter A. Schultz, 27-July-1998, for v2.23 (split basis)
c
c Revision history:
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
      DIMENSION  a1(n1),c1(n1), a2(n2),c2(n2)
c local declarations:
      DATA  zero,two / 0.d0,2.d0  /
c
c >>>> EXECUTABLE CODE:
c
      rtpid2 = rtpi/two
      lp1 = l + 1
c
      s = zero
      do  i2=1,n2
        do  i1=1,n1
          als = a1(i1) + a2(i2)
c 07Dec93-PAS: I choose the following form rather than pure
c   exponentiation of "als" cause this way more likely more accurate
          alsfac = SQRT(als)*als**lp1
          s = s + c1(i1)*c2(i2)/alsfac
        enddo
      enddo
c
      s = rtpid2*s
      do  ll=1,2*l+1,2
        s = s*DBLE(ll)/two
      enddo
c
      sg = s
c
c    That's all Folks!
c
      RETURN
      END
