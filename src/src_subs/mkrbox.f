c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKRBOX
c
c
      subroutine MKRBOX( hh,r0, nra,nrb,nrc, r )
c---------------------------------------------------------------
c Purpose: construct coordinates of regular grid from mesh
c          primitive vectors, with offset
c
c Written:  Peter A. Schultz, 1-May-1997, for v2.17
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c Not very optimized as not particularly expensive
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input array: (mesh definition vectors, r-offset)
      DIMENSION  hh(3,3), r0(3)
c output array: (output coordinates of grid points)
      DIMENSION  r(3,*)
c local declarations:
      DIMENSION  r1(3),r2(3),r3(3)
c
c >>>> EXECUTABLE CODE:
c
      do  i=1,3
c       offset initial vector
        r3(i) = r0(i)
      enddo
c
      ir = 0
      do 300 irc=1,nrc
        do  i=1,3
          r2(i) = r3(i)
          r3(i) = r3(i) + hh(i,3)
        enddo
c
        do 200 irb=1,nrb
          do  i=1,3
            r1(i) = r2(i)
            r2(i) = r2(i) + hh(i,2)
          enddo
c
          do 100 ira=1,nra
            ir = ir + 1
c           record coordinates of this grid pt.
            r(1,ir) = r1(1)
            r(2,ir) = r1(2)
            r(3,ir) = r1(3)
            r1(1) = r1(1) + hh(1,1)
            r1(2) = r1(2) + hh(2,1)
            r1(3) = r1(3) + hh(3,1)
  100     continue
  200   continue
  300 continue
c
c    That's all Folks!
c
      RETURN
      END
