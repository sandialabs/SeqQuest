c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMLAT
c
c
      subroutine SYMLAT( nsym, rmatsym, rprim )
c---------------------------------------------------------------
c Purpose: convert translation part of symmetry matrix operators
c          into lattice units.
c
c Written: Peter A. Schultz, 1-May-2002, for v2.53 (cell opt)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c  Variables:
c   nsym  - # (non-e) symmetry elements in group
c   rmatsym() - cartesian representations of symmetry elements
c   rprim(3,3) - primitive lattice vectors
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input/output array:
      DIMENSION  rmatsym(3,4,*)
c input array:
      DIMENSION  rprim(3,3)
c local declarations:
      DIMENSION  off(3), offlat(3)
c
c >>>> EXECUTABLE CODE:
c
      if( nsym.lt.1 ) RETURN
c
      do  isym=1,nsym
c
        do  i=1,3
          off(i) = rmatsym(i,4,isym)
        enddo
c
        call LVPROJ( rprim, off, offlat )
c
        do  i=1,3
          rmatsym(i,4,isym) = offlat(i)
        enddo
c
      enddo
c
      RETURN
c
c    That's all Folks!
c
      END
