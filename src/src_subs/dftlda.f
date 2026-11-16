c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTLDA
c
c
      subroutine DFTLDA( rho, vxc,exc )
c---------------------------------------------------------------
c Purpose: compute lda xc-potentials
c
c Written: Peter A. Schultz, Dec-2000, for v2.46 (gga)
c
c Revision history:
c  23Jun15-PAS/2.67: fix dfttyps declare
c   7Mar02-PAS/2.52: all explicit constants to d0
c   7Aug01-PAS/2.49: spin-polarization, and mimic PZ formulas
c  21Jun01-PAS/2.48: replace STOPs
c  19Dec00-PAS/2.46: to accomodate gga, renamed to "dftlda"
c---------------------------------------------------------------
c
c Input quantity:
c   rho = total electron number density
c Output quantities:
c   vxc = xc-potential (in rydbergs)
c   exc = xc energy density
c
c Notes:
c   For CAPZ = Perdew/Zunger param of Ceperley/Alder:
c   Potential not fully consistent with pseudo-potential generator.
c   Most notably, relativistic factors are NOT in the pseudopotential.
c   Otherwise minor numerics (it works out, trust me).  5Aug94-PAS
c
c   Fixed constants for vxc to at least 8 digits, noting that
c     vxc = exc - (rs/3)*(d[exc]/d[rs]) - 9Dec94-PAS
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0, one = 1.d0, two = 2.d0 )
      PARAMETER  ( third = 0.3333333333333333 d0 )
c  volfac=3/(4*pi)
      PARAMETER  ( volfac = 0.23873241463784 d0 )
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
        RETURN
      endif
c
c Flavors of LDA (only one right now):
c
      if( i_fcnal .eq. i_capz )then
c
c       CAPZ: Ceperley/Alder, parameterized by Perdew/Zunger:
c       Phys. Rev. B 23, 5048 (1981)
c
        rs = (volfac/rho)**third
c
crelativistics turned on:
        beta = 0.014d0 / rs
        rtb = SQRT( one + beta**2 )
        blog = LOG( beta + rtb )
        vrelfac = -one/two + 1.5d0*blog/(beta*rtb)
        erelfac = one - 1.5d0*((beta*rtb-blog)/beta**2)**2
c
crelativistics turned off:
crel        vrelfac = -one/two
crel        erelfac = one
c
c       Compute exchange part:
        Ex  =  -0.4582d0*erelfac/rs
        Vx  =  -0.61093333d0*vrelfac/rs
c
c       Compute correlation part:
        if( rs.lt.one )then
          rslog = LOG( rs )
          Ec  = 0.0311d0*rslog - 0.0480d0
     $        + 0.0020d0*rs*rslog       - 0.0116d0*rs
          Vc  = 0.0311d0*rslog - 0.058366665d0
     $        + 0.0013333333d0*rs*rslog - 0.0084d0*rs
        else
          rtrs = SQRT( rs )
          denom = 1.d0 + 1.0529d0*rtrs + 0.3334d0*rs
          Ec  = -0.1423d0 / denom
          Vc  = Ec *(1.d0+1.22838333d0*rtrs+0.44453333d0*rs) / denom
        endif
c
        exc = Ex + Ec
c
        vxc = Vx + Vc
c
c       End CAPZ potential
c
      else
c
        call STOPXERR( 'lda-type unknown' )
c
      endif
c
c Convert to Rydberg:
      exc = two*exc
      vxc = two*vxc
c
c    That's all Folks!
c
      RETURN
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
