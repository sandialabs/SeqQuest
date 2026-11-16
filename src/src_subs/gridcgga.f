c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDCGGA
c
c
      subroutine GRIDCGGA( xcrhocut, weight,nptr,
     $ rhuptot, vxcuptot,excuptot,  rhdntot, vxcdntot,excdntot,
     $ rhoup,  g_vxcup, g2_excup, gpggp,
     $ rhodn,  g2u_vxcdn, gd_excdn, grhoup )
c---------------------------------------------------------------
c Purpose: compute spin-polarized GGA correlation fcnals over grid
c
c Written: Peter A. Schultz, 20-July-2001, for 2.49 (spin)
c
c Revision history:
c  10Oct03-AEM+PAS/2.57: accommodate spin stress for GGA
c  26Jun03-AEM/2.56: added g2rhup arg. in DFTCGGA, needed for blypsp.
c  26Jun03-PAS/2.55: restrict rho(up/dn).gt.0, g2rho[up] input
c---------------------------------------------------------------
c
c ***** NOTE: correlation term only, no exchange *****
c
c Note on status of grid arrays _vxc(),_exc():
c  On input to this routine, these arrays contain the density
c  gradients that are inputs to the gga fcns.  Hence, they are
c  initialized here with the values of the correlation part of
c  the gga, and will be updated later with the exchange part.
c
c  Input arrays:
c    rhoup()    - spin-up electron density on grid
c    g_vxcup()  - | grad[rho] | "grho"
c    g2_excup() - grad**2[rho]  "g2rho"
c    gpggp()    - grad[rho] DOT grad | grad[rho] | "gpggp"
c    rhodn()    - spin-down electron density on grid
c    g2u_vxcdn()- grad**2[rhoup]
c    gd_excdn() - | grad[rho-dn] |
c    grhoup()   - | grad[rho-up] |
c  Output arrays [up/dn]:
c    _vxc()  - computed xc-potential on grid (in rydberg)
c    _exc()  - computed xc-energy density on grid
c  Output values [up/dn]:
c    rhotot  - grid integrated total electron density
c    vxctot  - grid integrated Vxc:rho
c    exctot  - grid integrated total xc-energy
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c input/output arrays:
      DIMENSION  rhoup(nptr), grhoup(nptr)
      DIMENSION  g_vxcup(nptr), g2_excup(nptr), gpggp(nptr)
      DIMENSION  rhodn(nptr)
      DIMENSION  g2u_vxcdn(nptr), gd_excdn(nptr)
c
c >>>> EXECUTABLE CODE:
c
      rhuptot = zero
      excuptot = zero
      vxcuptot = zero
      rhdntot = zero
      excdntot = zero
      vxcdntot = zero
c
      do  nr=1,nptr
        rhup = rhoup(nr)
        if( rhup.lt.zero ) rhup = zero
        rhuptot = rhuptot + rhup
        rhdn = rhodn(nr)
        if( rhdn.lt.zero ) rhdn = zero
        rhdntot = rhdntot + rhdn
        rhpt = rhup + rhdn
c
        if( rhpt .gt. xcrhocut )then
          grho = g_vxcup(nr)
          g2rh = g2_excup(nr)
          ggrh = gpggp(nr)
          g2rhup = g2u_vxcdn(nr)
          grhup = grhoup(nr)
          grhdn = gd_excdn(nr)
c
          call DFTCGGA( rhup, grho,g2rh,ggrh, vxcupr,excupr,
     $                  rhdn, grhup,grhdn,g2rhup, vxcdnr,excdnr,
     $                  dcgrho, dcgrhup, dcgrhdn,
     $                  dcg2rhup, dcg2rhdn )
c
          vxcuptot = vxcuptot + vxcupr*rhup
          excuptot = excuptot + excupr*rhup
          g_vxcup(nr) = vxcupr
          g2_excup(nr) = excupr
          vxcdntot = vxcdntot + vxcdnr*rhdn
          excdntot = excdntot + excdnr*rhdn
          g2u_vxcdn(nr) = vxcdnr
          gd_excdn(nr) = excdnr
        else
          g_vxcup(nr) = zero
          g2_excup(nr) = zero
          g2u_vxcdn(nr) = zero
          gd_excdn(nr) = zero
        endif
c
      enddo
c
c  Scale totals by mesh weighting:
c
      rhuptot = weight*rhuptot
      vxcuptot = weight*vxcuptot
      excuptot = weight*excuptot
      rhdntot = weight*rhdntot
      vxcdntot = weight*vxcdntot
      excdntot = weight*excdntot
c
c    That's all Folks!
c
      RETURN
      END
