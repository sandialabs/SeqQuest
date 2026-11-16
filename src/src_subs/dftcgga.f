c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTCGGA
c
c
      subroutine DFTCGGA( rhup, grho,g2rho,gpggp, vxcup,excup,
     $                    rhdn, grhup,grhdn,g2rhup, vxcdn,excdn,
     $                    dcgrho, dcgrhup, dcgrhdn,
     $                    dcg2rhup, dcg2rhdn )
c---------------------------------------------------------------
c Purpose: compute spin-polarized gga correlation potential
c
c Written: Peter A. Schultz, 10-August-2001, for v2.49 (spin)
c
c Revision history:
c  13Nov12-JEM/2.65: HSE added, piggybacked on PBE
C  17Spt08-RJM/2.62: added call to am05 spin  
c   8Oct03-AEM/2.57: 5 stress arguments added
c   8Oct03-AEM/2.57: argument deleted in DF_PW91C and DF_PBEC
c   7Oct03-AEM/2.57: added 5 arguments to DF_LYP
c  26Jun03-AEM/2.56: added g2rhup as argument, needed for blypsp.
c   8Apr02-PAS/2.53: gga stress
c---------------------------------------------------------------
c
c Input quantities
c   rhup = up-electron density
c   grho  = | GRAD[rho] |
c   g2rho = GRAD**2[rho]
c   gpggp = GRAD[rho] DOT GRAD[ | GRAD[rho] | ]
c   rhdn = dn-electron density
c   grhup = | GRAD[rhup] |
c   grhdn = | GRAD[rhdn] |
c   g2rhup = GRAD**2[rhup]
c Output quantities [up/dn]:
c   vxc = xc-potential (in rydbergs)
c   exc = xc energy density
c   dcgrho = stress term due to gradients (grho)
c   dcgrhup = stress term due to gradients up (grhup)
c   dcgrhdn = stress term due to gradients down (grhdn)
c   dcg2rhup = stress term due to laplacian up (g2rhup)
c   dcg2rhdn = stress term due to laplacian down (g2rhdn)
c
c Note:
c  Have to convert to Rydbergs for the GGA functionals
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Set parameters needed to construct various terms:
      PARAMETER  ( zero = 0.d0, one = 1.d0, two = 2.d0 )
      PARAMETER  ( third = 1.d0/3.d0, fourthrd = 4.d0/3.d0 )
c
      PARAMETER  ( xcrhocut = 1.d-21 )
c

      REAL*8  decsup,decsdn,decsupr,decsdnr
      REAL*8  twokfup,twokfdn

      COMMON  /DFTTYPS/ i_fcnal, i_gga,i_spin,
     $ i_capz,i_pw91,i_pbe,i_blyp,i_am05,i_hse,i_hse2
c
c >>>> EXECUTABLE CODE:
c
      rho = rhup + rhdn
      if( rho .lt. xcrhocut )then
        excup = zero
        vxcup = zero
        excdn = zero
        vxcdn = zero
        RETURN
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c           Need to assemble inputs to gga routines:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Relative spin polarization
      if( rhdn .gt. zero )then
        zet = ( rhup - rhdn ) / rho
      else
        zet = one
      endif
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
      if( i_fcnal .eq. i_pw91 )then
c
c       >>>>>>>>> PW91/GGA <<<<<<<<<
c
        call DF_LSDC( rs,zet,ec,vcup,vcdn,ecrs,eczet,zlfc )
        call DF_PW91C( rs,zet,g,ec,ecrs,eczet,t,uu,vv,ww,
     $   h,dvcup,dvcdn, dht )
c
        vxcup = vcup + dvcup
        vxcdn = vcdn + dvcdn
        exc = ec + h
c
c       Gradient terms for stress calculation:
        dht  = dht  * ( rho / twoksg**2 )
        dcgrho = dht
        dcgrhup = zero
        dcgrhdn = zero
        dcg2rhup = zero
        dcg2rhdn = zero
c
c       Convert to Rydberg units:
        vxcup = two*vxcup
        vxcdn = two*vxcdn
        excup = two*exc
        excdn = two*exc
        dcgrho = two*dcgrho
c
      elseif( i_fcnal .eq. i_pbe .or. i_fcnal .eq. i_hse
     $         .or. i_fcnal .eq. i_hse2 )then
c
c       >>>>>>>>> PBE/GGA <<<<<<<<<
c
c       lgga: 0=local exchange; 1=gga exchange
c       lpot: 0=only energy;    1=potential, too.
        lgga = 1
        lpot = 1
        call DF_PBEC( rs,zet,t,uu,vv,ww,lgga,lpot, ec,vcup,vcdn,
     $   h,dvcup,dvcdn, dht )
c
        vxcup = vcup + dvcup
        vxcdn = vcdn + dvcdn
        exc = ec + h
c
c       Gradient terms for stress calculation:
        dht  = dht  * ( rho / twoksg**2 )
        dcgrho = dht
        dcgrhup = zero
        dcgrhdn = zero
        dcg2rhup = zero
        dcg2rhdn = zero
c
c       Convert to Rydberg units:
        vxcup = two*vxcup
        vxcdn = two*vxcdn
        excup = two*exc
        excdn = two*exc
        dcgrho = two*dcgrho
c
      elseif( i_fcnal .eq. i_blyp )then
c
c       >>>>>>>>> BLYP/GGA <<<<<<<<<
c
        g2rhdn = g2rho - g2rhup
        call DF_LYP( rhup, rhdn, grho, grhup, grhdn, g2rhup, g2rhdn,
     $                 ec, vcup, vcdn, dcgrho, dcgrhup, dcgrhdn,
     $                 dcg2rhup, dcg2rhdn )
c
c       Convert to Rydberg units:
        vxcup = two*vcup
        vxcdn = two*vcdn
        excup = two*ec
        excdn = two*ec
        dcgrho = two*dcgrho
        dcgrhup = two*dcgrhup
        dcgrhdn = two*dcgrhdn
        dcg2rhup = two*dcg2rhup
        dcg2rhdn = two*dcg2rhdn
c

      elseif( i_fcnal .eq. i_am05 )then
c
c       >>>>>>>>>AM05/GGA <<<<<<<<<
c
        call DF_AM05s2( rhup, rhdn,grhup,grhdn,grho,ec,
     &                  vcup,vcdn,decsup,decsdn)
c       Gradient terms for stress calculation:
        twokfup=twokf*(two*rhup/(rhup+rhdn))**fourthrd
        twokfdn=twokf*(two*rhdn/(rhup+rhdn))**fourthrd
c
        decsup = decsup * two*( two*rhup / twokfup**two )
        decsdn = decsdn * two*( two*rhdn / twokfdn**two )
c
        dcgrho = zero
        dcgrhup = decsup
        dcgrhdn = decsdn
        dcg2rhup = zero
        dcg2rhdn = zero
c       Convert to Rydberg units:
        vxcup = two*vcup
        vxcdn = two*vcdn
        excup = two*ec
        excdn = two*ec
        dcgrhup = two*dcgrhup
        dcgrhdn = two*dcgrhdn
c
      else
c
c       >>>>>>>>> ERROR <<<<<<<<<
c
        call STOPXERR( 'gga-type/unknown gga type 1' )
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
