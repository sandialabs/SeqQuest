c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>STR2GPA
c
c
      subroutine STR2GPA( ndim, rprim, stress, strgpa, strlabel )
c---------------------------------------------------------------
c Purpose: convert stress from Ry/bohr^3 to "standard" unit: GPa
c
c Written: Peter A. Schultz, 1-May-2002, v2.53 (cell opt)
c
c Revision history:
c  14Oct02-PAS/2.53g: fix 2D conversion to SI stress
c  19May02-PAS/2.53b: 1D stress units fixed and changed
c   1May02-PAS/2.53: pushed down from "symfrc"
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
      DIMENSION  rprim(3,3), stress(3,3)
c Output:
      DIMENSION  strgpa(3,3)
      CHARACTER*8  strlabel
c Local:
      DIMENSION    v(3)
c
c   Units conversion from (Ry/au**ndim) for output ...
c     (1 GPa =1x10d10 N/m^2)
c     (Rydberg/bohr^3 = 1.471078x10d4 GPa = 147.1078x10d12 N/m^2)
c     (Rydberg/bohr^2 = 77.8461x10+2 N/m)
c     (Rydberg/bohr   = 41.1944x10-8 N)
      DIMENSION  strunits(3)
      DATA  strunits / 41.11944d1 , 77.8461d0 , 1.471078d4 /
      CHARACTER*8  strlbl(3)
c      ... and name of output units for N-dimensional stress
      DATA  strlbl / '10-7 N  ','10+2 N/m','GPa     ' /
c
c >>>> EXECUTABLE CODE:
c
      if( ndim.eq.0 ) RETURN
c
c Output stress result ...
c  ... work out "volume" factor for stress units conversion
      if( ndim.eq.1 )then
        vol = rprim(1,1)
      elseif( ndim.eq.2 )then
        call CROSS( rprim(1,1),rprim(1,2), v(1) )
        vol = SQRT( v(1)**2 + v(2)**2 + v(3)**2 )
      elseif( ndim.eq.3 )then
        call CROSS( rprim(1,2),rprim(1,3), v(1) )
        vol = ABS( rprim(1,1)*v(1) + rprim(2,1)*v(2) +
     $             rprim(3,1)*v(3) )
      else
        call STOPXERR( 'ndim-syf' )
      endif
c
      strlabel = strlbl(ndim)
c  ... convert stress into the appropriate units for dimension:
      strfac = strunits(ndim) / vol
      do  jd=1,ndim
        do  id=1,ndim
          strgpa(id,jd) = strfac*stress(id,jd)
        enddo
      enddo
c
      RETURN
c
c    That's all Folks!
c
      END
