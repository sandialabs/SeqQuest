c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDGGA
c
c
      subroutine GRIDGGA( xcrhocut, weight,nptr,
     $ rhotot, vxctot,exctot,
     $ rho, g_vxc,g2_exc,gpggp )
c---------------------------------------------------------------
c Purpose: compute spin-unpolarized GGA functionals over grid,
c          with derivative term needed for stress
c
c Written: Peter A. Schultz, 20-July-2001, for 2.49 (spin)
c
c Revision history:
c  11Jun03-PAS/2.55: split off stress-specific call
c   5Apr02-PAS/2.53: gga stress
c---------------------------------------------------------------
c
c Variables:
c  Input arrays:
c    rho()    - total electron density on grid
c    g_vxc()  - | grad[rho] | "grho"
c    g2_exc() - grad**2[rho]  "g2rho"
c    gpggp() - grad[rho] DOT grad | grad[rho] | "gpggp"
c  Output arrays:
c    g_vxc()  - computed xc-potential on grid (in rydberg)
c    g2_exc() - computed xc-energy density on grid
c  Output values:
c    rhotot  - grid integrated total electron density
c    vxctot  - grid integrated Vxc:rho
c    exctot  - grid integrated total xc-energy
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c input arrays:
      DIMENSION  rho(nptr)
      DIMENSION  g_vxc(nptr), g2_exc(nptr), gpggp(nptr)
c
c >>>> EXECUTABLE CODE:
c
      rhotot = zero
      vxctot = zero
      exctot = zero
c
      do  nr=1,nptr
        rhpt = rho(nr)
        rhotot = rhotot + rhpt
c
        if( rhpt .gt. xcrhocut )then
          grh  = g_vxc(nr)
          g2rh = g2_exc(nr)
          ggrh = gpggp(nr)
c
          call DFTGGA( rhpt, grh,g2rh,ggrh,
     $                 vxcpt,excpt, dxcpt,dxc2pt )
c
          vxctot = vxctot + vxcpt*rhpt
          exctot = exctot + excpt*rhpt
          g_vxc(nr) = vxcpt
          g2_exc(nr) = excpt
        else
          g_vxc(nr) = zero
          g2_exc(nr) = zero
        endif
c
      enddo
c
c  Scale totals by mesh weighting:
c
      rhotot = weight*rhotot
      vxctot = weight*vxctot
      exctot = weight*exctot
c
c    That's all Folks!
c
      RETURN
      END
