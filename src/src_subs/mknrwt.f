c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKNRWT
c
c
      subroutine MKNRWT( alf, wtsum, orig, hh, n1,n2,n3, wtwk )
c---------------------------------------------------------------
c Purpose: construct grid weighting function for "nearby" value
c
c Written: Peter A. Schultz, 3-September-1993
c
c Revision history:
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  orig(3), hh(3,3)
      DIMENSION  wtwk(*)
c local declarations:
      DIMENSION  r1(3),r2(3),r3(3)
c
c >>>> EXECUTABLE CODE:
c
      do  i=1,3
        r3(i) = orig(i) - hh(i,1) - hh(i,2) - hh(i,3)
      enddo
c
      ir = 0
      do 33 j3=1,n3
        do  i=1,3
          r3(i) = r3(i) + hh(i,3)
          r2(i) = r3(i)
        enddo
        do 32 j2=1,n2
          do  i=1,3
            r2(i) = r2(i) + hh(i,2)
            r1(i) = r2(i)
          enddo
          do 31 j1=1,n1
            r1(1) = r1(1) + hh(1,1)
            r1(2) = r1(2) + hh(2,1)
            r1(3) = r1(3) + hh(3,1)
            ir = ir + 1
            wtwk(ir) = r1(1)**2 + r1(2)**2 + r1(3)**2
   31     continue
   32   continue
   33 continue
c
      do  ir=1,n1*n2*n3
        wtwk(ir) = EXP(-alf*wtwk(ir))
        wtsum = wtsum + wtwk(ir)
      enddo
c
      RETURN
      END
