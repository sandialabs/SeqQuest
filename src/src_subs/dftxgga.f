c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTXGGA
c
c
      subroutine DFTXGGA( rho, rop,grho,g2rho,gpggp, vxr,exr, dexsr )
c---------------------------------------------------------------
c Purpose: compute spin-polarized gga exchange potentials
c
c Written: Peter A. Schultz, 10-August-2001, for v2.49 (spin)
c
c Revision history:
c  13Nov12-JEM/2.65: HSE added
c  17Spt08-RJM/2.62: added opp spin passing for AM05 spin stress 
c  14Oct03-AEM/2.57: stress terms doubled to give correct value
c   7Oct03-AEM/2.57: stress term added
c  25Jun03-AEM/2.55: blyp added
c   8Apr02-PAS/2.53: gga stress
c---------------------------------------------------------------
c
c Input quantities
c   rho   = total electron number density
c   grho  = | GRAD[rho] |
c   g2rho = GRAD**2[rho]
c   gpggp = GRAD[rho] DOT GRAD[ | GRAD[rho] | ]
c Output quantities
c   vxr = xc-potential (in rydbergs)
c   exr = xc energy density
c   dexs = stress term
c
c Note:
c  Have to convert to Rydbergs for the GGA functionals
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
        exr = zero
        vxr = zero
        dexsr = zero
        RETURN
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c           Need to assemble inputs to gga routines:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      call DF_GGA0X( rho, grho,g2rho,gpggp, twokf, s,u,v )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                      Flavors of DFT:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( i_fcnal .eq. i_pw91 )then
c
c       >>>>>>>>> PW91/GGA <<<<<<<<<
c
        call DF_PW91X( rho, s,u,v, ex,vx, dexs )
c
c       Gradient terms for stress calculation:
        dexs = dexs * ( rho / twokf**2 )
c       This only gives half the required term, so double dexs
        dexs = two*dexs
c
c       Convert to Rydberg units:
        vxr = two*vx
        exr = two*ex
        dexsr = two*dexs
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
c
c       Gradient terms for stress calculation:
        dexs = dexs * ( rho / twokf**2 )
c       This only gives half the required term, so double dexs
        dexs = two*dexs
c
c       Convert to Rydberg units:
        vxr = two*vx
        exr = two*ex
        dexsr = two*dexs
c
      elseif( i_fcnal .eq. i_blyp )then
c
c       >>>>>>>>> BLYP/GGA <<<<<<<<<
c
        call DF_B1X( rho, s, u, v, ex,vx, dexs )
c
c       Gradient terms for stress calculation:
        dexs = dexs * ( rho / twokf**2 )
c       This only gives half the required term, so double dexs
        dexs = two*dexs
c
c       Convert to Rydberg units:
        vxr = two*vx
        exr = two*ex
        dexsr = two*dexs
c
      elseif( i_fcnal .eq. i_am05 )then
c
c       >>>>>>>>>AM05/GGA <<<<<<<<<
c
c       Note: rho=2*spin a, rop=spin b
c
        call DF_AM05s1( rho, rop,s, u, v, ex,vx, dexs )
c
c       Gradient terms for stress calculation:
        dexs = dexs * ( rho / twokf**2 )
c       This only gives half the required term, so double dexs
        dexs = two*dexs
c       Convert to Rydberg units:
        vxr = two*vx
        exr = two*ex
        dexsr = two*dexs
c
      elseif( i_fcnal .eq. i_hse )then
c
c       >>>>>>>>> HSE/GGA <<<<<<<<<
c
c       NOTE: HSE06 parameters, omega=0.11, fraction=0.25
c
        call DF_HSEX(rho,s,u,v,0.11d0,ex,vx,dexs)
        ex = ex * 0.25d0
        vx = vx * 0.25d0
        dexs = dexs * 0.25d0
c
c       Gradient terms for stress calculation:
        dexs = dexs * ( rho / twokf**2 )
c       This only gives half the required term, so double dexs
        dexs = two*dexs
c
c       Convert to Rydberg units:
        vxr = two*vx
        exr = two*ex
        dexsr = two*dexs
c
      elseif( i_fcnal .eq. i_hse2 )then
c
c       >>>>>>>>> HSE/GGA <<<<<<<<<
c
c       NOTE: HSE12 parameters, omega=0.098, fraction=0.313
c
        call DF_HSEX(rho,s,u,v,0.098d0,ex,vx,dexs)
        ex = ex * 0.313d0
        vx = vx * 0.313d0
        dexs = dexs * 0.313d0
c
c       Gradient terms for stress calculation:
        dexs = dexs * ( rho / twokf**2 )
c       This only gives half the required term, so double dexs
        dexs = two*dexs
c
c       Convert to Rydberg units:
        vxr = two*vx
        exr = two*ex
        dexsr = two*dexs
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
