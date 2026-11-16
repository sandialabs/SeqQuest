c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DFTSLDA
c
c
      subroutine DFTSLDA( rhup, vxcup,excup, rhdn, vxcdn,excdn )
c---------------------------------------------------------------
c Purpose: compute spin-polarized lda xc-potentials
c
c Written: Peter A. Schultz, 7-August-2001, for v2.49 (spin)
c
c Revision history:
c  23Jun15-PAS/2.67: fix dfttyps declare
c   7Mar02-PAS/2.52: all explicit constants to d0
c---------------------------------------------------------------
c
c Input quantities:
c   rh[up/dn] = [up/dn]-spin electron number density
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
      PARAMETER  ( r4d3  = 1.333333333333333 d0 )
c  volfac=3/(4*pi)
      PARAMETER  ( volfac = 0.23873241463784 d0 )
      PARAMETER  ( xcrhocut = 1.d-21 )
c
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
c zet = relative spin polarization
c fzet = spin interpolation factor (PZ's, from von Barth/Hedin)
c dfzet = d(fzet)/d(zet)
c
      if( rhdn .ne. zero )then
        zet = ( rhup - rhdn ) / rho
        fzet = ( one + zet )**r4d3 + ( one - zet )**r4d3 - two
        fdenom = two**r4d3 - two
        fzet = fzet / fdenom
        dfzet = r4d3*( (one + zet)**third - (one - zet)**third )
        dfzet = dfzet / fdenom
      else
c       Get the hydrogen atom right
        zet = one
        fzet = one
        dfzet = two * third * ( one + one / ( two**third - one ) )
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
crel        vrelfac = one
crel        erelfac = one
c
c U/P=Unpolarized/Polarized limits
c
c       Compute exchange part:
        ExU =  -0.4582d0*erelfac/rs
        VxU =  -0.61093333d0*vrelfac/rs
        ExP =  ExU * two**third
        VxP =  VxU * two**third
c
c       Compute correlation part:
        if( rs.lt.one )then
          rslog = LOG( rs )
          EcU = 0.0311d0*rslog - 0.0480d0
     $        + 0.0020d0*rs*rslog       - 0.0116d0*rs
          VcU = 0.0311d0*rslog - 0.058366666d0
     $        + 0.0013333333d0*rs*rslog - 0.0084d0*rs
          EcP = 0.01555d0*rslog - 0.0269d0
     $        + 0.0007d0*rs*rslog       - 0.0048d0*rs
          VcP = 0.01555d0*rslog - 0.032083333d0
     $        + 0.0004666666d0*rs*rslog - 0.0034333333d0*rs
        else
          rtrs = SQRT( rs )
          denom = 1.d0 + 1.0529d0*rtrs + 0.3334d0*rs
          EcU = -0.1423d0 / denom
          VcU = EcU*(1.d0+1.22838333d0*rtrs+0.44453333d0*rs) / denom
          denom = 1.d0 + 1.3981d0*rtrs + 0.2611d0*rs
          EcP = -0.0843d0 / denom
          VcP = EcP*(1.d0+1.63111666d0*rtrs+0.34813333d0*rs) / denom
        endif
c
        ExcU = ExU + EcU
        ExcP = ExP + EcP
        Exc = ExcU + fzet * ( ExcP - ExcU )
        excup = Exc
        excdn = Exc
c
        VxcU = VxU + VcU
        VxcP = VxP + VcP
        Vxc0 = VxcU + fzet * ( VxcP - VxcU )
        Vxc1 = ( ExcP - ExcU )*dfzet
        facup =  one - zet
        facdn = -one - zet
        vxcup = Vxc0 + facup*Vxc1
        vxcdn = Vxc0 + facdn*Vxc1
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
      excup = two*excup
      vxcup = two*vxcup
      excdn = two*excdn
      vxcdn = two*vxcdn
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
