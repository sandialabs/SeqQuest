c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DF_SGGA0
c
c
      subroutine DF_GGA0( rho, grho,g2rho,gpggp,  zet,grhup,grhdn,
     $ g, rs, twokf,twoksg, s,u,v, t,uu,vv, ww )
c---------------------------------------------------------------
c Purpose: set up gga inputs from raw electron density/gradients
c
c Written: Peter A. Schultz, January-2001, for v2.46 (gga)
c
c Revision history:
c   3Apr02-PAS/2.53: gga stress
c   8Aug01-PAS/2.49: spin-polarized dft, entry df_gga0x added
c---------------------------------------------------------------
c
c Input quantities:
c   rho   = total electron number density
c   grho  = | GRAD[rho] |
c   g2rho = GRAD**2[rho]
c   gpggp = GRAD[rho] DOT GRAD[ | GRAD[rho] | ]
c   grhup = | GRAD[rhup] |
c   grhdn = | GRAD[rhdn] |
c   zet = (up-dn)/rho = relative spin polarization
c Output quantities:
c   described individually below
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Set parameters needed to construct various terms:
      PARAMETER  ( zero = 0.d0, one = 1.d0, two = 2.d0 )
      PARAMETER  ( twothrd = 0.6666666666666667 d0 )
      PARAMETER  ( thrd = 0.3333333333333333 d0 )
      PARAMETER  ( C_pi = 3.141592653589793 d0 )
      PARAMETER  ( C_pi2= 9.869604401089358 d0 )
c
c C_Rs = (3/(4pi))**(1/3)
      PARAMETER  ( C_Rs = 0.6203504908994001 d0 )
c C_Kf = (3*pi**2)**(1/3)
      PARAMETER  ( C_Kf = 3.093667726280136 d0 )
c C_Ks = sqrt( 4*C_Kf / pi )
      PARAMETER  ( C_Ks = 1.984686395219856 d0 )
c C_Rf = C_Rs*C_Kf = ( 9*pi/4 )**(1/3) = Burke's "alpha"
      PARAMETER  ( C_Rf = 1.919158292677513 d0 )
c
c >>>> EXECUTABLE CODE:
c
c  g = ( (1+zet)**(2/3) + (1-zet)**(2/3) ) / 2
c     With no spin, zet=0, and g=1.
      if( zet .eq. zero )then
        g = one
      else
        g = ( ( one + zet )**twothrd + ( one - zet )**twothrd ) / two
      endif
c
      rhthrd = rho**thrd
c
c  rs = (3/(4pi*rho))**(1/3) = local Seitz radius=alpha/fk
      rs = C_Rs / rhthrd
c
      s0 = grho / rho
      u0 = gpggp / (rho*rho)
      v0 = g2rho / rho
c
c  fk = Kf = C_Kf * rho**(1/3) = local fermi wavevector
      fk = C_Kf*rhthrd
c  s = |grho| / (rho * (2*Kf))
c  u =  gpggp / (rho**2 * (2*Kf)**3)
c  v =  g2rho / (rho * (2*Kf)**2)
      twokf = two * fk
      s = s0 / (twokf)
      u = u0 / (twokf**3)
      v = v0 / (twokf**2)
c
c  sk = Ks = C_Ks * rho**(1/6) = Thomas-Fermi screening wavevector
      sk = C_Ks*SQRT( rhthrd )
c  t  = |grho| / (rho * (2*Ks*g))
c  uu =  gpggp / (rho**2 * (2*Ks*g)**3)
c  vv =  g2rho / (rho * (2*Ks*g))**2)
      twoksg = two * sk * g
      t  = s0 / (twoksg)
      uu = u0 / (twoksg**3)
      vv = v0 / (twoksg**2)
c
c  ww = zero for non-spin-polarized
c  ww = GRAD[rho]*GRAD[zet] / ( rho * (2*Ks*g)**2 ) [Hamann pspot]
c  ww=(|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2 [Burke]
      if( zet .eq. zero )then
        ww = zero
      else
        ww = ( grhup**2 - grhdn**2 - zet*grho**2 )
     $       / ( rho*twoksg )**2
      endif
c
c Put rho into twokf and twoksg for use in stress gga term:
      twokf = rho*twokf
      twoksg = rho*twoksg
c
c    That's all Folks!
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DF_GGA0X
c
c
      entry DF_GGA0X( rho, grho,g2rho,gpggp, twokf, s,u,v )
c---------------------------------------------------------------
c Purpose: set up gga exchange inputs from electron density/gradients
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Input quantities:
c   rho   = electron number density
c   grho  = | GRAD[rho] |
c   g2rho = GRAD**2[rho]
c
c >>>> EXECUTABLE CODE:
c
      rhthrd = rho**thrd
c
      s0 = grho / rho
      u0 = gpggp / (rho*rho)
      v0 = g2rho / rho
c
c  fk = Kf = C_Kf * rho**(1/3) = local fermi wavevector
      fk = C_Kf*rhthrd
c  s = |grho| / (rho * (2*Kf))
c  u =  gpggp / (rho**2 * (2*Kf)**3)
c  v =  g2rho / (rho * (2*Kf)**2)
      twokf = two * fk
      s = s0 / (twokf)
      u = u0 / (twokf**3)
      v = v0 / (twokf**2)
c
c Put rho into twokf for use in stress gga term:
      twokf = rho*twokf
c
c    That's all Folks!
c
      RETURN
      END
