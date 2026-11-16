c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRIDCSTR
c
c
      subroutine GRIDCSTR( xcrhocut, nptr, weight,
     $   rhoa,  rhob,  g_dxc,   g2_d2up,   gg_d2dn,
     $   grhoa,     grhob,     g2rhoa,     dxcup,     dxcdn )
c---------------------------------------------------------------
c Purpose: compute extra fcnals for correlation part of
c          spin-dependent GGA stress
c
c Written: Peter A. Schultz and Ann. E. Mattsson
c          10-October-2003, for 2.57
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Variables:
c  Input arrays:
c    rhoa()    - total up-spin xc electron density on grid
c    rhob()    - total dn-spin xc electron density on grid
c    g_dxc()   - | grad[rho] | "grho"
c    g2_d2up() - grad**2[rho]  "g2rho"
c    gg_d2dn() - grad[rho] DOT grad | grad[rho] | "gpggp"
c    grhoa()   - | grad[rho-up] |
c    grhob()   - | grad[rho-dn] |
c    g2rhoa()  - grad**2.[rho-up]
c  Output arrays:
c    g_dxc()   - fcnal for GGA total 1st derivative stress
c    g2_d2up() - fcnal for GGA up-spin Laplacian stress
c    gg_d2dn() - fcnal for GGA dn-spin Laplacian stress
c    dxcup()   - fcnal for GGA up-spin 1st derivative stress
c    dxcdn()   - fcnal for GGA dn-spin 1st derivative stress
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
c Input:
      DIMENSION  rhoa(nptr),rhob(nptr)
      DIMENSION  grhoa(nptr),grhob(nptr)
      DIMENSION  g2rhoa(nptr)
c Input/output:
      DIMENSION  g_dxc(nptr),g2_d2up(nptr),gg_d2dn(nptr)
c Output:
      DIMENSION  dxcup(nptr),dxcdn(nptr)
c
c >>>> EXECUTABLE CODE:
c
      ectot  = zero
      ectota = zero
      ectotb = zero
      do  nr=1,nptr
        rhup = rhoa(nr)
        if( rhup.lt.zero ) rhup = zero
        rhdn = rhob(nr)
        if( rhdn.lt.zero ) rhdn = zero
        rhpt = rhup + rhdn
c
        if( rhpt .gt. xcrhocut )then
          grh = g_dxc(nr)
          g2rh = g2_d2up(nr)
          ggrh = gg_d2dn(nr)
          g2rhup = g2rhoa(nr)
          grhup = grhoa(nr)
          grhdn = grhob(nr)
c
          call DFTCGGA( rhup, grh,g2rh,ggrh, vxcupr,excupr,
     $                  rhdn, grhup,grhdn,g2rhup, vxcdnr,excdnr,
     $                  dcgrho, dcgrhup, dcgrhdn,
     $                  dcg2rhup, dcg2rhdn )
c
          ectota = ectota + excupr*rhup
          ectotb = ectotb + excdnr*rhdn
c
          g_dxc(nr) = dcgrho
          dxcup(nr) = dcgrhup
          dxcdn(nr) = dcgrhdn
          g2_d2up(nr) = dcg2rhup
          gg_d2dn(nr) = dcg2rhdn
        else
          g_dxc(nr) = zero
          dxcup(nr) = zero
          dxcdn(nr) = zero
          g2_d2up(nr) = zero
          gg_d2dn(nr) = zero
        endif
c
      enddo
      ectota = weight*ectota
      ectotb = weight*ectotb
      ectot = ectota + ectotb
C      write(6,'(a,3f20.10)') 'CSTR: energy ectot,a,b=',
C     $ ectot,ectota,ectotb
c
      RETURN
c
c    That's all Folks!
c
      END
