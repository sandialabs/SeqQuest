c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDSTR
c
c
      subroutine GRIDSTR( xcrhocut, nptr, weight,
     $ rho, g_dxc,g2_dxc2,gpggp )
c---------------------------------------------------------------
c Purpose: compute extra fcnals for non-spin GGA stress
c
c Written: Peter A. Schultz, 11-JUne-2003, for 2.55
c
c Revision history:
c   9Oct03-PAS/2.57: accomodate spin-dependent stress for GGA
c---------------------------------------------------------------
c
c Variables:
c  Input arrays:
c    rho()     - total electron density on grid
c    g_dxc()   - | grad[rho] | "grho"
c    g2_dxc2() - grad**2[rho]  "g2rho"
c    gpggp()   - grad[rho] DOT grad | grad[rho] | "gpggp"
c  Output arrays:
c    g_dxc()   - fcnal for GGA 1st dervative stress
c    g2_dxc2() - fcnal for GGA Laplacian stress
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
c Input:
      DIMENSION  rho(nptr)
c Input/output:
      DIMENSION  g_dxc(nptr), g2_dxc2(nptr)
c Output:
      DIMENSION  gpggp(nptr)
c
c >>>> EXECUTABLE CODE:
c
      exctot = zero
      do  nr=1,nptr
        rhpt = rho(nr)
c
        if( rhpt .gt. xcrhocut )then
          grh  = g_dxc(nr)
          g2rh = g2_dxc2(nr)
          ggrh = gpggp(nr)
c
          call DFTGGA( rhpt, grh,g2rh,ggrh,
     $                 vxcpt,excpt, dxcpt,dxc2pt )
c
          exctot = exctot + rhpt*excpt
c
          g_dxc(nr) = dxcpt
          g2_dxc2(nr) = dxc2pt
        else
          g_dxc(nr) = zero
          g2_dxc2(nr) = zero
        endif
c
      enddo
      exctot = weight*exctot
C      write(6,'(a,3f20.10)') 'STR: energy exctot=',exctot
c
      RETURN
c
c    That's all Folks!
c
      END
