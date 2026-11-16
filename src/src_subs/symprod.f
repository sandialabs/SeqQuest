c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMPROD
c
c
      subroutine SYMPROD( rmata, rmatb, rmatab, rprim )
c---------------------------------------------------------------
c Purpose: take product of symmetry elements
c
c Written: Peter A. Schultz, 4-January-1999, for v2.30
c
c Revision history:
c  21Jul99-DBR/2.37: bug fix - properly rotate translation
c---------------------------------------------------------------
c
c  Note:
c   We are actually computing the product BA, not AB
c  Variables:
c   rmata(),rmatb() - input symmetry elements
c   rmatab()   - output product of symmetry elements
c   rprim(3,3) - primitive lattice vectors
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output array:
      DIMENSION  rmatab(3,4)
c input arrays:
      DIMENSION  rmata(3,4),rmatb(3,4)
      DIMENSION  rprim(3,3)
c local declarations:
      LOGICAL    outuc, shift
      DIMENSION  t1(3)
c
c >>>> EXECUTABLE CODE:
c
c  Operate on symmetry operation a with operation b
      do  j=1,3
        do  i=1,3
          rmatab(i,j) = rmata(i,1)*rmatb(1,j)
     $                + rmata(i,2)*rmatb(2,j)
     $                + rmata(i,3)*rmatb(3,j)
        enddo
      enddo
c
c  Transform the translation vector as well:
      do  i=1,3
        rmatab(i,4) = rmatb(1,i)*rmata(1,4)
     $              + rmatb(2,i)*rmata(2,4)
     $              + rmatb(3,i)*rmata(3,4)
      enddo
c
c  Complete translation vector:
      t1(1) = rmatab(1,4) + rmatb(1,4)
      t1(2) = rmatab(2,4) + rmatb(2,4)
      t1(3) = rmatab(3,4) + rmatb(3,4)
c
c  Project translation vector into unit cell:
      call PROJUC( outuc,shift, 3, t1,rmatab(1,4), rprim )
c
c    That's all Folks!
c
      RETURN
      END
