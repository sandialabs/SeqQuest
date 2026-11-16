c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKSHL2D
c
c
      subroutine MKSHL2D( lshl, nrs, rngr,rsq, a0e,a1e,a2e,
     $ rho,drho,ddrh )
c---------------------------------------------------------------
c Purpose: produce radial derivatives of reference shells
c          from sums of Occ*Ca*Cb*(a+b)**N*exp(-(a+b)*r**2)
c
c Written: Peter A. Schultz, 11-January-2001, for v2.46 (gga)
c
c Revision history:
c  21Jun01-PAS/2.48: replace STOPs
c---------------------------------------------------------------
c
c  We have accumulated these sums over the gaussians of a shell:
c   aNe(r) = SUM[a;b] [ Occ*Ca*Cb * (a+b)**N * exp(-(a+b)*r**2) ]
c  for N=0,1,2.  Want to assemble density and derivatives from
c  these sums, for a shell of angular momentum l=lshl.
c
c  rho(r)      = r**(2l) * a0e(r)
c  d(rho)/dr   = [l=0] -2*r*a1e(r)
c              = [l>0] -2*r**(2l+1)*a1e(r) + 2l*r**(2l-1)*a0e(r)
c              = r**(2l-1) * ( -2*r2*a1e(r) + 2l*a0e(r) )
c  d2(rho)/dr2 = [l=0] 4*r2*a2e(r) - 2*a1e(r)
c              = [l>0] 4*r**(2l+2)*a2e(r) - (8l+2)*r**2l*a1e(r)
c                    + (2l)*(2l-1)*a0e(r)
c              = r**(2l-2) * ( 4*r4*a2e(r) - (8l+2)*r2*a1e(r) +
c                              (2l)*(2l-1)*a0e )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c input arrays: (r,r2, sums)
      DIMENSION  rngr(nrs), rsq(nrs)
      DIMENSION  a0e(nrs), a1e(nrs), a2e(nrs)
c output arrays:
      DIMENSION  rho(nrs), drho(nrs), ddrh(nrs)
c
c local declarations:
      PARAMETER ( zero=0.d0, one=1.d0, two=2.d0, four=4.d0 )
c
c >>>> EXECUTABLE CODE:
c
      twol = DBLE( 2*lshl )
c
      twolm1 = twol - one
      cab0 = twol*twolm1
      cab1 = two*( two*twol + one )
c
      if( lshl.eq.0 )then
        do  nr=1,nrs
          rho(nr) = rho(nr) + a0e(nr)
          drho(nr) = drho(nr)
     $     - a1e(nr)*two*rngr(nr)
          ddrh(nr) = ddrh(nr) +
     $     a2e(nr)*four*rsq(nr) - a1e(nr)*two
        enddo
      elseif( lshl.eq.1 )then
        do  nr=1,nrs
          r = rngr(nr)
          r2 = rsq(nr)
          rho(nr) = rho(nr) + a0e(nr)*r2
          drho(nr) = drho(nr) + r*
     $     ( a0e(nr)*twol - a1e(nr)*two*r2 )
          ddrh(nr) = ddrh(nr) +
     $     ( a2e(nr)*four*r2*r2 - a1e(nr)*cab1*r2 + a0e(nr)*cab0 )
         enddo
      elseif( lshl.eq.2 )then
        do  nr=1,nrs
          r = rngr(nr)
          r2 = rsq(nr)
          rho(nr) = rho(nr) + a0e(nr)*r2*r2
          drho(nr) = drho(nr) + r*r2*
     $     ( a0e(nr)*twol - a1e(nr)*two*r2 )
          ddrh(nr) = ddrh(nr) + r2*
     $     ( a2e(nr)*four*r2*r2 - a1e(nr)*cab1*r2 + a0e(nr)*cab0 )
        enddo
      elseif( lshl.gt.2 )then
        do  nr=1,nrs
          r = rngr(nr)
          r2 = rsq(nr)
          r2lm2 = r2**(lshl-1)
          rho(nr) = rho(nr) + a0e(nr)*r2*r2lm2
          drho(nr) = drho(nr) + r*r2lm2*
     $     ( a0e(nr)*twol - a1e(nr)*two*r2 )
          ddrh(nr) = ddrh(nr) + r2lm2*
     $     ( a2e(nr)*four*r2*r2 - a1e(nr)*cab1*r2 + a0e(nr)*cab0 )
        enddo
      else
c       This should never happen but test is free.
        call STOPXERR( 'SH2-lshl/ bad lshl in MKSHL2D' )
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKSHL1D
c
c
      entry MKSHL1D( lshl, nrs, rngr,rsq, a0e,a1e,
     $ rho,drho )
c---------------------------------------------------------------
c Purpose: shell density and first derivative only.
c---------------------------------------------------------------
c
c >>>> EXECUTABLE CODE:
c
      twol = DBLE( 2*lshl )
c
      if( lshl.eq.0 )then
        do  nr=1,nrs
          rho(nr) = rho(nr) + a0e(nr)
          drho(nr) = drho(nr)
     $     - a1e(nr)*two*rngr(nr)
        enddo
      elseif( lshl.eq.1 )then
        do  nr=1,nrs
          r = rngr(nr)
          r2 = rsq(nr)
          rho(nr) = rho(nr) + a0e(nr)*r2
          drho(nr) = drho(nr) + r*
     $     ( a0e(nr)*twol - a1e(nr)*two*r2 )
         enddo
      elseif( lshl.eq.2 )then
        do  nr=1,nrs
          r = rngr(nr)
          r2 = rsq(nr)
          rho(nr) = rho(nr) + a0e(nr)*r2*r2
          drho(nr) = drho(nr) + r*r2*
     $     ( a0e(nr)*twol - a1e(nr)*two*r2 )
         enddo
      elseif( lshl.gt.2 )then
        do  nr=1,nrs
          r = rngr(nr)
          r2 = rsq(nr)
          r2lm2 = r2**(lshl-1)
          rho(nr) = rho(nr) + a0e(nr)*r2*r2lm2
          drho(nr) = drho(nr) + r*r2lm2*
     $     ( a0e(nr)*twol - a1e(nr)*two*r2 )
        enddo
      else
c       This should never happen but test is free.
        call STOPXERR( 'SH1-lshl/ bad lshl in MKSHL1D' )
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKSHL
c
c
      entry MKSHL( lshl, nrs, rsq, a0e, rho )
c---------------------------------------------------------------
c Purpose: shell density only, no derivatives.
c---------------------------------------------------------------
c
c >>>> EXECUTABLE CODE:
c
      if( lshl.eq.0 )then
        do  nr=1,nrs
          rho(nr) = rho(nr) + a0e(nr)
        enddo
      elseif( lshl.eq.1 )then
        do  nr=1,nrs
          rho(nr) = rho(nr) + a0e(nr) * rsq(nr)
        enddo
      elseif( lshl.eq.2 )then
        do  nr=1,nrs
          rho(nr) = rho(nr) + a0e(nr) * rsq(nr)*rsq(nr)
       enddo
      elseif( lshl.gt.2 )then
        do  nr=1,nrs
          rho(nr) = rho(nr) + a0e(nr) * rsq(nr)**lshl
        enddo
      else
c       This should never happen but test is free.
        call STOPXERR( 'SHL-lshl/ bad lshl in MKSHL' )
      endif
c
      RETURN
c
c    That's all Folks!
c
      END
