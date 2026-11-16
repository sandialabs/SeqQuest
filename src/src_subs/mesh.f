c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MESH
c
c
      subroutine MESH( n1r,n2r,n3r,rprim, weight,hh, gprim )
c---------------------------------------------------------------
c Purpose: setup grid intervals, compute recip. lattice vecs
c
c Revision history:
c   8Aug03-PAS/2.56: use KCELL to generate fft cell
c  27Aug02-PAS/2.53g: sign error in rt-hand cell to gprim
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
c Input: real space primitive lattice vectors:
      DIMENSION  rprim(3,3)
c Output: real space integration grid primitive vectors
      DIMENSION  hh(3,3)
c         reciprocal space primitive cell for fft:
      DIMENSION  gprim(3,3)
c
c Local:
      DIMENSION  tvec(3)
c
c >>>> EXECUTABLE CODE:
c
c Generate fft reciprocal space vectors (i.e. full 3-D)
      kdim = 3
      call KCELL( kdim, rprim, gprim )
c
      r1n = DBLE( n1r )
      r2n = DBLE( n2r )
      r3n = DBLE( n3r )
      do  i=1,3
        hh(i,1) = rprim(i,1) / r1n
        hh(i,2) = rprim(i,2) / r2n
        hh(i,3) = rprim(i,3) / r3n
      enddo
c
      call CROSS( hh(1,2), hh(1,3), tvec )
      vol = ABS( hh(1,1)*tvec(1)+hh(2,1)*tvec(2)+hh(3,1)*tvec(3) )
      weight = vol
c
c    That's all Folks!
c
      RETURN
      END
