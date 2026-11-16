c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTGGA
c
c
      subroutine DFTGGA( rho, grho,g2rho,gpggp, vxc,exc, dexc,dexc2 )
c---------------------------------------------------------------
c Purpose: compute gga potentials of selected functional flavor
c          ***** spin-unpolarized version *****
c
c Written: Peter A. Schultz, Dec-2000, for v2.46 (gga)
c
c Revision history:
c  13Nov12-JEM/2.65: HSE added
c  24May06-AEM/2.60: AM05 added
c  14Oct03-AEM/2.57: LYP stress term corrected
c   8Oct03-AEM/2.57: argument deleted in DF_PW91C and DF_PBEC
c   7Oct03-AEM/2.56: Spin lyp used instead of non-spin lyp
c  13Jun03-AEM/2.55: dexc2 stress term due to Laplacians added
c  11Jun03-AEM/2.55: blyp added (with stress)
c   3Apr02-PAS/2.53: gga stress
c   7Aug01-PAS/2.49: spin-polarization
c  21Jun01-PAS/2.48: replace STOPs
c---------------------------------------------------------------
c
c Input quantities
c   rho   = electron number density
c   grho  = | GRAD[rho] |
c   g2rho = GRAD**2[rho]
c   gpggp = GRAD[rho] DOT GRAD[ | GRAD[rho] | ]
c Output quantities:
c   vxc = xc-potential in rydbergs
c   exc = xc energy density
c   dexc = ( d[exc*rho] / d[ d[rho]/d[x_i] ] ) / d[rho]/d[x_i]
c        = term needed for stress calculation in gga.
c   dexc2 = 2 * d[exc*rho] / d[ grad**2[rho] ]
c Note:
c  Must convert functionals to Rydbergs
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Set parameters needed to construct various terms:
      PARAMETER  ( zero = 0.d0, one = 1.d0, two = 2.d0 )
c
      PARAMETER  ( xcrhocut = 1.d-21 )
c
      COMMON  /DFTTYPS/ i_fcnal, i_gga,i_spin,
     $ i_capz,i_pw91,i_pbe,i_blyp,i_am05,i_hse,i_hse2
c
c >>>> EXECUTABLE CODE:
c
      if( rho .lt. xcrhocut )then
        exc = zero
        vxc = zero
        dexc = zero
        dexc2 = zero
        RETURN
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c           Need to assemble inputs to gga routines:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c No spin polarization (zet=0), and up/dn gradients are half of total
      zet = zero
      grhup = grho/two
      grhdn = grho/two
c
c LYP uses grho, grhup, grhdn, and g2rho as inparameters.
c Please never change them!
c
      call DF_GGA0( rho, grho,g2rho,gpggp,  zet, grhup,grhdn,
     $ g, rs, twokf,twoksg, s,u,v, t,uu,vv, ww )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                      Flavors of DFT:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Output quantities:
c  vxc = total xc potential (Rydberg units)
c  exc = total xc energy density (Rydberg units)
c Intermediate quantities:
c  ex    = exchange energy
c  vx    = exchange potential
c  ec    = lsd correlation energy
c  vcup  = lsd up correlation potential
c  vcdn  = lsd down correlation potential
c  h     = gradient correction to correlation energy
c  dvcup = gradient correction to up correlation potential
c  dvcdn = gradient correction to down correlation potential
c  ecrs  = d(ec)/d(rs)
c  eczet = d(ec)/d(zet)
c  zlfc  = "spin stiffness"
c
      if( i_fcnal .eq. i_capz )then
c
c       >>>>>>>>> CAPZ/LDA <<<<<<<<<
c
c       Ceperley/Alder parameterized by Perdew/Zunger
c       D.M. Ceperley and B.J. Alder, PRL 45, 566 (1980).
c       J. Perdew and A. Zunger,  PRB 23, 5048 (1981).
c
        call DFTLDA( rho, vxc,exc )
c
      elseif( i_fcnal .eq. i_pw91 )then
c
c       >>>>>>>>> PW91/GGA <<<<<<<<<
c
        call DF_PW91X( rho, s,u,v, ex,vx, dexs )
        call DF_LSDC( rs,zet,ec,vcup,vcdn,ecrs,eczet,zlfc )
        call DF_PW91C( rs,zet,g,ec,ecrs,eczet,t,uu,vv,ww,
     $   h,dvcup,dvcdn, dht )
c
        vxc = vx + ( vcup + vcdn + dvcup + dvcdn ) / two
        exc = ex + ec + h
c
c       Gradient terms for stress calculation:
        dexs = dexs * ( rho / twokf**2 )
        dht  = dht  * ( rho / twoksg**2 )
        dexc = dexs + dht
        dexc2 = zero
c
c       Convert to Rydberg units:
        vxc = two*vxc
        exc = two*exc
        dexc = two*dexc
c
      elseif( i_fcnal .eq. i_pbe )then
c
c       >>>>>>>>> PBE/GGA <<<<<<<<<
c
c       lgga: 0=local exchange; 1=gga exchange
c       lpot: 0=only energy;    1=potential, too.
        lgga = 1
        lpot = 1
        call DF_PBEX( rho, s,u,v, lgga,lpot, ex,vx, dexs )
        call DF_PBEC( rs,zet,t,uu,vv,ww,lgga,lpot, ec,vcup,vcdn,
     $   h,dvcup,dvcdn, dht )
c
        vxc = vx + ( vcup + vcdn + dvcup + dvcdn ) / two
        exc = ex + ec + h
c
c       Gradient terms for stress calculation:
        dexs = dexs * ( rho / twokf**2 )
        dht  = dht  * ( rho / twoksg**2 )
        dexc = dexs + dht
        dexc2 = zero
c
c       Convert to Rydberg units:
        vxc = two*vxc
        exc = two*exc
        dexc = two*dexc
c
      elseif( i_fcnal .eq. i_blyp )then
c
c       >>>>>>>>> BLYP/GGA <<<<<<<<<
c
        call DF_B1X( rho, s, u, v, ex, vx, dexs )
c
c
c       Note: LYP use the bare gradient and laplacians, grho and g2rho.
c       It also gives the full correlation out in ec, in contrast
c       to other correlation functionals.
c
        rhup = rho / two
        rhdn = rho / two
        g2rhup = g2rho / two
        g2rhdn = g2rho / two
        call DF_LYP( rhup, rhdn, grho, grhup, grhdn, g2rhup, g2rhdn,
     $                 ec, vcup, vcdn, dcgrho, dcgrhup, dcgrhdn,
     $                 dcg2rhup, dcg2rhdn )
c
c
        vxc = vx + ( vcup + vcdn ) / two
        exc = ex + ec
c
c       Gradient terms for stress calculation:
        dexs = dexs * ( rho / twokf**2 )
        dexc = dexs + dcgrho + ( dcgrhup + dcgrhdn ) / ( two*two )
        dexc2 = ( dcg2rhup + dcg2rhdn ) / two
c
c       Convert to Rydberg units:
        vxc = two*vxc
        exc = two*exc
        dexc = two*dexc
        dexc2 = two*dexc2
c
      elseif( i_fcnal .eq. i_am05 )then 
c
c       >>>>>>>>> AM05 <<<<<<<<<
c
        call DF_AM05( rho, s, u, v, ex, ec, vx, vc, dexs, decs )
c
c
c       Note: Implemented without separate exchange
c       and correlation. Only spinless am05
c       available yet, thus I use vc, not vcup and vcdn.
c
        vxc = vx + vc
        exc = ex + ec
        dexc = ( dexs + decs ) * ( rho / twokf**2 )
        dexc2 = zero
c        
c       Convert to Rydberg units:
        vxc = two*vxc
        exc = two*exc
        dexc = two*dexc
c
      elseif( i_fcnal .eq. i_hse )then
c
c       >>>>>>>>> HSE/GGA <<<<<<<<<
c
c       NOTE: HSE06 parameters, omega=0.11, fraction=0.25
c       lgga: 0=local exchange; 1=gga exchange
c       lpot: 0=only energy;    1=potential, too.
        lgga = 1
        lpot = 1
        call DF_HSEX(rho,s,u,v,0.11d0,ex,vx,dexs)
        call DF_PBEC( rs,zet,t,uu,vv,ww,lgga,lpot, ec,vcup,vcdn,
     $   h,dvcup,dvcdn, dht )
        ex = ex * 0.25d0
        vx = vx * 0.25d0
        dexs = dexs * 0.25d0
c
        vxc = vx + ( vcup + vcdn + dvcup + dvcdn ) / two
        exc = ex + ec + h
c
c       Gradient terms for stress calculation:
        dexs = dexs * ( rho / twokf**2 )
        dht  = dht  * ( rho / twoksg**2 )
        dexc = dexs + dht
        dexc2 = zero
c
c       Convert to Rydberg units:
        vxc = two*vxc
        exc = two*exc
        dexc = two*dexc
c
      elseif( i_fcnal .eq. i_hse2 )then
c
c       >>>>>>>>> HSE/GGA <<<<<<<<<
c
c       NOTE: HSE12 parameters, omega=0.098, fraction=0.313
c       lgga: 0=local exchange; 1=gga exchange
c       lpot: 0=only energy;    1=potential, too.
        lgga = 1
        lpot = 1
        call DF_HSEX(rho,s,u,v,0.098d0,ex,vx,dexs)
        call DF_PBEC( rs,zet,t,uu,vv,ww,lgga,lpot, ec,vcup,vcdn,
     $   h,dvcup,dvcdn, dht )
        ex = ex * 0.313d0
        vx = vx * 0.313d0
        dexs = dexs * 0.313d0
c
        vxc = vx + ( vcup + vcdn + dvcup + dvcdn ) / two
        exc = ex + ec + h
c
c       Gradient terms for stress calculation:
        dexs = dexs * ( rho / twokf**2 )
        dht  = dht  * ( rho / twoksg**2 )
        dexc = dexs + dht
        dexc2 = zero
c
c       Convert to Rydberg units:
        vxc = two*vxc
        exc = two*exc
        dexc = two*dexc
c
      else
c
c       >>>>>>>>> ERROR <<<<<<<<<
c
        call STOPXERR( 'gga-type/unknown gga type' )
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
