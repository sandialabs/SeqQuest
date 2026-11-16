c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FERMIPOP
c
c
      subroutine FERMIPOP( efermi,etemp, eltotl,nlvl,nocc,
     $ occfull, occlvl, elvl, wlvl )
c---------------------------------------------------------------
c Purpose: find occupations and elec# for given efermi and temp.
c          Use a Fermi-Dirac filling function.
c
c Written: Peter A. Schultz, 15-February-2005, for v2.59
c          (extracted from solver/pop routines
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c  Routine assumes eigenvalues are ordered in increasing order.
c  Filling factor of state with eigenvalue e(i), fermi level Ef,
c   and electronic temperature et is given by the Fermi-Dirac
c   distribution:
c         1 / ( 1 + exp[ {e(i)-Ef}/et ] )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER ( zero=0.d0, one=1.d0 )
c
c Input:
      DIMENSION  elvl(nlvl), wlvl(nlvl)
c Output:
      DIMENSION  occlvl(nlvl)
c
c Local:
c   cutocc=argument cutoff for exp, to take exp(-tol)=0
c   Probably tol=19. is good enough (r*4 filling precision)
c   but we give it a little extra here just to be safe
c   We could also use a limit expansion computing the
c   filling factor for the larger exponents, but why?
      DATA  cutocc  / 30.d0 /
c
c >>>> EXECUTABLE CODE:
c
c  Standard Fermi-Dirac filling function
c
      elsum = zero
      do 210 ilvl=1,nlvl
        tol = cutocc
        xp = ( elvl(ilvl) - efermi ) / etemp
        if( xp .lt. -tol )then
          ff = one
        elseif( xp .lt. tol )then
          ff = one / ( one + EXP( xp ) )
        else
          ff = zero
          goto 211
        endif
        occlvl(ilvl) = ff*occfull
        elsum = elsum + occlvl(ilvl)*wlvl(ilvl)
        nocc = ilvl
  210 continue
  211 continue
c
      eltotl = elsum
c
      RETURN
c
c    That's all Folks!
c
      END
