c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VNORMS
c
c
      subroutine VNORMS( vnorm )
c---------------------------------------------------------------
c Purpose: initialize orbital angular normalizations
c
c Written: Peter A. Schultz, for v2.52 (linux)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c This set of constants normalizes the angular functions:
c   1,         (x,y,z),     (xy,yz,zx),   (x**2-y**2),   (3*z**2-r**2)
c They equal, respectively:
c   sqrt(4pi), sqrt(4pi/3), sqrt(4pi/15), sqrt(16pi/15), sqrt(16pi/5)
c
      PARAMETER  ( c_s    = 3.544907701811 d0,
     $             c_p    = 2.046653415893 d0,
     $             c_xy   = 0.9152912328638d0,
     $             c_x2y2 = 1.830582465728 d0,
     $             c_dz2  = 3.170661838085 d0 )
c
      DIMENSION  vnorm(*)
c
c >>>> EXECUTABLE CODE:
c
c s
      vnorm(1) = c_s
c x,y,z
      vnorm(2) = c_p
c xy,yz,zx
      vnorm(3) = c_xy
c x2-y2
      vnorm(4) = c_x2y2
c dz2: 3z2-r2
      vnorm(5) = c_dz2
c
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VNORMS9
c
      entry VNORMS9( vnorm )
c s
      vnorm(1) = c_s
c x,y,z
      vnorm(2) = c_p
      vnorm(3) = c_p
      vnorm(4) = c_p
c xy,yz,zx
      vnorm(5) = c_xy
      vnorm(6) = c_xy
      vnorm(7) = c_xy
c x2-y2
      vnorm(8) = c_x2y2
c dz2: 3z2-r2
      vnorm(9) = c_dz2
c
      RETURN
c
c    That's all Folks!
c
      END
