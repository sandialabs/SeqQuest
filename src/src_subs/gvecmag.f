c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GVECMAG
c
c
      subroutine GVECMAG( n1r,n2r,n3r,nptr,gprim, gvecsq )
c---------------------------------------------------------------
c Purpose: compute magnitudes of recip. lattice squared for use
c      in iteration phase solution of poisson's equation.
c
c Written by: Peter A. Schultz, for v2.00
c
c Revision history:
c  18Nov99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  gprim(3,3),gvecsq(nptr)
c local declarations:
      DIMENSION  g(3),g2(3),g3(3)
c
c >>>> EXECUTABLE CODE:
c
      f1 = DBLE( (n1r+2)/2 )
      f2 = DBLE( (n2r+2)/2 )
      f3 = DBLE( (n3r+2)/2 )
      do  i=1,3
        g3(i) = - f1*gprim(i,1) - f2*gprim(i,2) - f3*gprim(i,3)
      enddo
c
      nr = 0
      do 30 j3=1,n3r
        do  i=1,3
          g3(i) = g3(i) + gprim(i,3)
          g2(i) = g3(i)
        enddo
        do 20 j2=1,n2r
          do  i=1,3
            g2(i) = g2(i) + gprim(i,2)
            g(i)  = g2(i)
          enddo
          do 10 j1=1,n1r
            g(1) = g(1) + gprim(1,1)
            g(2) = g(2) + gprim(2,1)
            g(3) = g(3) + gprim(3,1)
            nr = nr + 1
            gvecsq(nr) = g(1)*g(1) + g(2)*g(2) + g(3)*g(3)
   10     continue
   20   continue
   30 continue
c
c    That's all Folks!
c
      RETURN
      END
