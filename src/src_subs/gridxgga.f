
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDXGGA
c
c
      subroutine GRIDXGGA( xcrhocut, weight,nptr,
     $ rhtot, vxtot,extot, rho,rop, grho,g2rho,gpggp, vxc,exc )
c---------------------------------------------------------------
c Purpose: compute spin-pol GGA exchange functionals over grid
c
c Written: Peter A. Schultz, 10-August-2001, for 2.49 (spin)
c
c Revision history:
c  17Spt08-RJM/2.62: added opposite spin passing for AM05 spin 
c  10Oct03-AEM+PAS/2.57: accommodate spin stress for GGA
c---------------------------------------------------------------
c
c ***** NOTE: exchange term only, correlation previously *****
c
c Note on spin-dependent exchange:
c   Ex(up,dn) = ( Ex(2*up) + Ex(2*dn) ) / 2
c   where Ex refer to energy (not energy densities)
c   Hence, we want to use *twice* the input density as input
c   to the routines that generate the potentials and energy
c   densities, and gpggp gets doubled twice because it has
c   two factors of density in it.  Once we get the potential
c   and energy density, then we can use the undoubled density
c   to do the integrals, and this takes care of the final
c   division by two.
c
c Note on status of the grid potential and energy-density:
c   We presume that exc and vxc have already been through
c   the correlation gga calculation and on input contain the
c   correlation part of the gga and what we do here is
c   increment exc and vxc with the exchange part of the gga.
c
c Variables:
c  Input arrays:
c    rho()    - total electron density on grid
c    grho()  - | grad[rho] | "grho"
c    g2rho() - grad**2[rho]  "g2rho"
c    gpggp()  - grad[rho] DOT grad | grad[rho] | "gpggp"
c  Input/output arrays:
c    vxc() - computed xc-potential on grid (in rydberg)
c    exc() - computed xc-energy density on grid
c  Output values:
c    rhtot  - grid integrated total electron density
c    vxtot  - grid integrated Vx:rho
c    extot  - grid integrated total Ex-energy
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0, two = 2.d0 )
c input/output arrays:
      DIMENSION  rho(nptr),rop(nptr)
      DIMENSION  grho(nptr), g2rho(nptr), gpggp(nptr)
      DIMENSION  vxc(nptr), exc(nptr)
c
      REAL*8 rhop
c
c >>>> EXECUTABLE CODE:
c
      rhtot = zero
      extot = zero
      vxtot = zero
c
      do  nr=1,nptr
        rh = rho(nr)
        rhop = rop(nr)
        rhtot = rhtot + rh
        rhpt = rh * two
c
        if( rhpt .gt. xcrhocut )then
          grh = grho(nr) * two
          g2rh = g2rho(nr) * two
          ggrh = gpggp(nr) * two * two
c
          call DFTXGGA( rhpt, rhop,grh,g2rh,ggrh, vxr,exr, dexsr )
c
          vxtot = vxtot + vxr*rh
          extot = extot + exr*rh
          vxc(nr) = vxc(nr) + vxr
          exc(nr) = exc(nr) + exr
        endif
c
      enddo
c
c  Scale totals by mesh weighting:
c
      rhtot = weight*rhtot
      vxtot = weight*vxtot
      extot = weight*extot
c
c    That's all Folks!
c
      RETURN
      END
