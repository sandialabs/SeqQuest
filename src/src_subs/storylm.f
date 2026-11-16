c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STORYLM
c
c
      subroutine STORYLM( nang, angsamp,sphharm )
c---------------------------------------------------------------
c Purpose: calculates real spherical harmonics on the
c          McLaren 72-points on the unit sphere
c          and stores them for use in the xc-potential fit.
c
c Written: PJF
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  21Apr93-PAS: general angle dimension/pretty up code
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input array - angular sampling
      DIMENSION  angsamp(3,nang)
c output array - real spherical harmonics on sphere
      DIMENSION  sphharm(nang,9)
c local declarations:
      DIMENSION  ylm(2,6)
c
c      Loop over all nang points:
c
      do 10 iang=1,nang
c
        call YLMREAL( ylm,2,angsamp(1,iang) )
c
c  Rearrange order of ylm's and store
c
c > Y(0,0), proportional to 1
        sphharm(iang,1) = ylm(1,1)
c > Real [ Y(1,1) ], proportional to x
        sphharm(iang,2) = ylm(1,4)
c > Imag [ Y(1,1) ], proportional to y
        sphharm(iang,3) = ylm(2,4)
c > Y(1,0), proportional to z
        sphharm(iang,4) = ylm(1,2)
c > Imag [ Y(2,2) ], proportional to xy
        sphharm(iang,5) = ylm(2,6)
c > Imag [ Y(2,1) ], proportional to yz
        sphharm(iang,6) = ylm(2,5)
c > Real [ Y(2,1) ], proportional to zx
        sphharm(iang,7) = ylm(1,5)
c > Real [ Y(2,2) ], proportional to x**2-y**2
        sphharm(iang,8) = ylm(1,6)
c > Y(2,0), proportional to 2.*z**2-x**2-y**2
        sphharm(iang,9) = ylm(1,3)
c
   10 continue
c
      RETURN
      END
