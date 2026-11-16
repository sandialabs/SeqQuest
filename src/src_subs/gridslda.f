c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDSLDA
c
c
      subroutine GRIDSLDA( xcrhocut, weight,nptr,
     $ rhouptot, vxcuptot,excuptot,  rhodntot, vxcdntot,excdntot,
     $ rhoup, vxcup,excup,  rhodn, vxcdn,excdn )
c---------------------------------------------------------------
c Purpose: compute spin-polarized LDA functionals over igrid
c
c Written: Peter A. Schultz, 20-July-2001, for 2.49 (spin)
c
c Revision history:
c  19Apr03-PAS/2.55: restrict rho-up/dn to be non-negative
c---------------------------------------------------------------
c
c Variables: (up/dn)
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
      DIMENSION  rhoup(nptr),  vxcup(nptr),excup(nptr)
      DIMENSION  rhodn(nptr),  vxcdn(nptr),excdn(nptr)
c
c >>>> EXECUTABLE CODE:
c
      rhouptot = zero
      excuptot = zero
      vxcuptot = zero
      rhodntot = zero
      excdntot = zero
      vxcdntot = zero
c
      do  nr=1,nptr
        rhup = rhoup(nr)
        if( rhup .lt. zero ) rhup = zero
        rhouptot = rhouptot + rhup
        rhdn = rhodn(nr)
        if( rhdn .lt. zero ) rhdn = zero
        rhodntot = rhodntot + rhdn
        rhpt = rhup + rhdn
c
        if( rhpt .gt. xcrhocut )then
c
          call DFTSLDA( rhup, vxcupr,excupr,
     $                  rhdn, vxcdnr,excdnr )
c
          vxcuptot = vxcuptot + vxcupr*rhup
          excuptot = excuptot + excupr*rhup
          vxcup(nr) = vxcupr
          excup(nr) = excupr
          vxcdntot = vxcdntot + vxcdnr*rhdn
          excdntot = excdntot + excdnr*rhdn
          vxcdn(nr) = vxcdnr
          excdn(nr) = excdnr
        else
          vxcup(nr) = zero
          excup(nr) = zero
          vxcdn(nr) = zero
          excdn(nr) = zero
        endif
c
      enddo
c
c  Scale totals by mesh weighting:
c
      rhouptot = weight*rhouptot
      vxcuptot = weight*vxcuptot
      excuptot = weight*excuptot
      rhodntot = weight*rhodntot
      vxcdntot = weight*vxcdntot
      excdntot = weight*excdntot
c
c    That's all Folks!
c
      RETURN
      END
