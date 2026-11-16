c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDLDA
c
c
      subroutine GRIDLDA( xcrhocut,weight,nptr,
     $ rhotot, vxctot,exctot,
     $ rho, vxc,exc )
c---------------------------------------------------------------
c Purpose: compute spin-restricted LDA functionals over grid
c
c Written: Peter A. Schultz, 20-July-2001, for 2.49 (spin)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Variables:
c  Input arrays:
c    rho()   - total electron density on grid
c  Output arrays:
c    vxc()   - computed xc-potential on grid (in rydberg)
c    exc()   - computed xc-energy density on grid
c  Output values:
c    rhotot  - grid integrated total electron density
c    vxctot  - grid integrated Vxc:rho
c    exctot  - grid integrated total xc-energy
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c input/output arrays:
      DIMENSION  rho(nptr),  vxc(nptr),exc(nptr)
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
c
          call DFTLDA( rhpt, vxcpt,excpt )
c
          vxctot = vxctot + vxcpt*rhpt
          exctot = exctot + excpt*rhpt
          vxc(nr) = vxcpt
          exc(nr) = excpt
        else
          vxc(nr) = zero
          exc(nr) = zero
        endif
c
      enddo
c
c  Scale totals by mesh weighting:
c
      rhotot = weight*rhotot
      exctot = weight*exctot
      vxctot = weight*vxctot
c
c    That's all Folks!
c
      RETURN
      END
