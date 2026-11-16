c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KCELL
c
c
      subroutine KCELL( kdim, rprim, gprim )
c---------------------------------------------------------------
c Purpose: produce BZ cell vectors appropriate to dimensionality
c
c Written: Peter A. Schultz, 5-August-2003, for v2.55a
c
c Revision history:
c   4Sep08-PAS/2.61g: once more fix handedness of gprim vectors
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
c Input: real-space primitive vectors
      DIMENSION  rprim(3,3)
c Output: reciprocal space primitive vectors
      DIMENSION  gprim(3,3)
c Local:
      DIMENSION  rvec(3,3),gvec(3,3)
c
c >>>> EXECUTABLE CODE:
c
c Copy over the primitive vectors:
      do  jd=1,3
        do  i=1,3
          rvec(i,jd) = rprim(i,jd)
        enddo
      enddo
c
c Project out bad stuff from r-vectors:
      if( kdim.lt.3 )then
        do  idk=1,kdim
          dotkk = rvec(1,idk)*rvec(1,idk) + rvec(2,idk)*rvec(2,idk)
     $          + rvec(3,idk)*rvec(3,idk)
          do  idx=kdim+1,3
            fac = ( rvec(1,idk)*rvec(1,idx) + rvec(2,idk)*rvec(2,idx)
     $            + rvec(3,idk)*rvec(3,idx) ) / dotkk
            rvec(1,idx) = rvec(1,idx) - fac*rvec(1,idk)
          enddo
        enddo
      endif
c
c Produce reduced primitive reciprocal vectors:
      call CROSS( rvec(1,2),rvec(1,3), gvec(1,1) )
      call CROSS( rvec(1,3),rvec(1,1), gvec(1,2) )
      call CROSS( rvec(1,1),rvec(1,2), gvec(1,3) )
      vol = rvec(1,1)*gvec(1,1) + rvec(2,1)*gvec(2,1) +
     $      rvec(3,1)*gvec(3,1)
c
c Transfer the results to the output g-prim:
      twopi = pi + pi
      scfac = twopi / vol
      do  jd=1,3
        do  i=1,3
          gprim(i,jd) = scfac*gvec(i,jd)
        enddo
      enddo
      vol = ABS( vol )
c
c    That's all Folks!
c
      RETURN
      END
