c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DUPEMAT
c
c
      subroutine DUPEMAT( nd, kv, rv, nv,rvset )
c---------------------------------------------------------------
c Purpose: find a match, if there is one, for a candidate
c          vector in a set of vectors.
c
c Written: Peter A. Schultz, 4-January-1999, for v2.30
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c  Variables:
c   nd      - dimension of vector
c   rv()    - vector to be checked
c   rvset() - set of vectors to be checked against
c   nv      - number of vectors in set to be checked against
c   kv      - index of matched vector (0 if there is none)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  rv(nd), rvset(nd,*)
c
c  criteria for difference:
      DATA  vsmall / 2.d-6 /
c
c >>>> EXECUTABLE CODE
c
      kv = 0
      do 20 iv=1,nv
        do 10 id=1,nd
          dif = ABS( rv(id) - rvset(id,iv) )
          if( dif.gt.vsmall ) goto 20
   10   continue
        kv = iv
        goto 90
   20 continue
c
   90 continue
c
      RETURN
      END
