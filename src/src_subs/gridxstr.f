c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDXSTR
c
c
      subroutine GRIDXSTR( xcrhocut, nptr, weight,
     $   rho,   rop,grho,      g2rho,     ggrho,
     $   dxc )
c---------------------------------------------------------------
c Purpose: compute extra fcnal for exchange part of
c          spin-dependent GGA stress
c
c Written: Peter A. Schultz and Ann. E. Mattsson
c          10-October-2003, for 2.57
c
c Revision history:
c  17Spt08-RJM/2.62: added opposite spin (rop) to be passed to AM05 spin
c---------------------------------------------------------------
c
c Variables:
c  Input arrays:
c    rho()    - total spin xc electron density on grid
c    grho()   - | grad[rho] | "grho"
c    g2rho()  - grad**2[rho]  "g2rho"
c    ggrho()  - grad[rho] DOT grad | grad[rho] | "gpggp"
c  Output arrays:
c    dxc()   - fcnal for GGA 1st derivative stress
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0, two = 2.d0 )
c
c Input:
      DIMENSION  rho(nptr),rop(nptr)
      DIMENSION  grho(nptr)
      DIMENSION  g2rho(nptr)
      DIMENSION  ggrho(nptr)
c Output:
      DIMENSION  dxc(nptr)
c
c >>>> EXECUTABLE CODE:
c
      extot = zero
      do  nr=1,nptr
        rh = rho(nr)
        rhop = rop(nr)
        rhpt = rh * two
c
        if( rhpt .gt. xcrhocut )then
          grh = grho(nr) * two
          g2rh = g2rho(nr) * two
          ggrh = ggrho(nr) * two * two
c
          call DFTXGGA( rhpt, rhop, grh, g2rh, ggrh, vxr,exr, dexsr )
c
          extot = extot + exr*rh
          dxc(nr) = dxc(nr) + dexsr
        endif
c
      enddo
      extot = weight*extot
C      write(6,'(a,3f20.10)') 'XSTR: energy extotsp=',extot
c
      RETURN
c
c    That's all Folks!
c
      END
