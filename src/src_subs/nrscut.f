c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NRSCUT
c
c
      subroutine NRSCUT( nrmn,nrmx,nrs,argmx,
     $ exparg,pmag,alsum,radmsh,w1 )
c---------------------------------------------------------------
c Purpose: find nr at which point shtrng exponential underflow
c
c Written: Peter A. Schultz
c
c Revision history:
c  24Sep08-PAS/2.61i: clean tab characters from source
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  radmsh(nrs)
      DIMENSION  w1(nrs)
c
      nrmn = 1
      nrmx = nrs
c
      rn = radmsh(nrs)
      argrn = exparg + alsum*rn*rn - pmag*rn
      if( exparg.lt.argmx .and. argrn.lt.argmx ) RETURN
      argmx2 = argmx - exparg
c
      do  nr=1,nrs
        r = radmsh(nr)
        w1(nr) = alsum*r*r - pmag*r
      enddo
c
      if( argrn.gt.argmx )then
        do 20 nr=nrs,1,-1
          if( w1(nr).lt.argmx2 ) goto 25
          nrmx = nr - 1
   20   continue
        RETURN
   25   continue
      endif
c
      if( exparg.gt.argmx )then
        do 30 nr=1,nrs
          if( w1(nr).lt.argmx2 ) goto 35
          nrmn = nr + 1
   30   continue
        RETURN
   35   continue
      endif
c
      RETURN
      END
