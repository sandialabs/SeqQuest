c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PROJUC
c
c
      subroutine PROJUC( outuc,shift, ndim, rv,rvnew, rlat )
c---------------------------------------------------------------
c Purpose: project atoms into unit cell (in periodic dimensions)
c
c Written: Peter A. Schultz, 26-January-1996, for v2.15
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c  Code selects atoms as inside uc shifted by "edge" toward
c  zero sides of primary unit cell.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      LOGICAL  outuc,shift
      DIMENSION  rv(3),rvnew(3),rlat(3,3)
c local declarations:
      DIMENSION  v(3)
      DATA  edge / 0.0001d0 /
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      outuc = .false.
      do  i=1,3
        rvnew(i) = rv(i)
      enddo
c project coordinates onto lattice vectors:
      call PROJLV( rv, rlat(1,1),rlat(1,2),rlat(1,3), v )
c
      shift = .false.
      do 100 idim=1,3
        vchk = v(idim) + edge
        ilat = ABS(vchk)
        if( v(idim)+edge .lt. zero ) ilat = -ilat - 1
c
        if( ilat.ne.0 .and. idim.le.ndim )then
          shift = .true.
          clat = DBLE(ilat)
          do  i=1,3
            rvnew(i) = rvnew(i) - clat*rlat(i,idim)
          enddo
        elseif( ilat.ne.0 )then
          outuc = .true.
        endif
  100 continue
c
      RETURN
      END
