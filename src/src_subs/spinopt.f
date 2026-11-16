c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SPINOPT
c
c
      subroutine SPINOPT( nspin, fermilvl, spinold,spinnew,
     $ etemp, elecno, nstate,nk, wtk,eigval,
     $ wklvl )
c---------------------------------------------------------------
c Purpose: compute fermi level across spin, associated spinpol
c
c Written: Peter A. Schultz, 11-January-2008, for 2.61
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER ( zero=0.d0, one=1.d0, two=2.d0 )
c
c Input:
      DIMENSION  wtk(nk), eigval(nstate,nk,*)
      DIMENSION  elecno(2)
c Scratch:
      DIMENSION  wklvl(*)
c
c Local:
      DIMENSION  elspin(2)
      LOGICAL    xch,pos,neg
      DATA  elcerr,elconv / 1.d-6,1.d-11 /
      DATA  occnil / 1.d-9 /
c Local:
      DIMENSION  occlvl(2)
      DATA       occlvl / 2.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      if( nspin .ne. 2 ) call STOPXERR( 'spinopt/no-spin' )
c
      nstk = nstate*nk
      nstall = nstk*2
c
      occfull = occlvl(nspin)
      elnumbr = elecno(1) + elecno(2)
c
      call POP( nspin, spindif,
     $ fermilvl, egap,
     $ etemp, elnumbr, esumsp, occfull,
     $ nstate,nk, wtk, eigval,wklvl(1),npop,
     $ wklvl(1+1*nstall),wklvl(1+2*nstall),wklvl(1+3*nstall),
     $ wklvl(1*4*nstall),wklvl(1+5*nstall) )
c -->  elvl wlvl klvl nlvl flvl
c
      spinold = elecno(1) - elecno(2)
      spinnew = spindif
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SPUPDATE
c
c
      subroutine SPUPDATE( istat_sp, istep_sp, spindata,
     $ elecno, spinold,spinnew,spinges,
     $ nstep_sp,meth_sp,conv_sp,blend_sp )
c---------------------------------------------------------------
c Purpose: compute update for spin polarization
c
c Written: Peter A. Schultz, 11-January-2008, for 2.61
c
c Revision history:
c   2Dec10-PAS/2.62: prevent jump to negative spin on first step
c---------------------------------------------------------------
c
c I have two schemes here.
c One "simple" is just to use the up/dn electron count ("spinnew")
c  from a fermi level equilibration with the eigenspectrum with
c  the current up/dn electron count ("spinold").
c The other, called "linear", is to take the intersection of the
c  lines given by the input up/dn electron count ("old") and the
c  output up/dn electron count ("new"), with constraints to make
c  sure that things behave reasonably.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER ( zero=0.d0, one=1.0d0, two=2.d0 )
c
c Input/Output:
      DIMENSION  elecno(2)
      DIMENSION  spindata(4,2)
c Local:
      DIMENSION  occlvl(2)
      DATA       occlvl / 2.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      elnumbr = elecno(1) + elecno(2)
      spindel = spinnew - spinold
      spinchg = ABS( spindel )
c
c For the moment, only implement the SIMPLE scheme.
c
      istat_sp = 0
      if( spinchg .lt. conv_sp )then
c       Spin is converged
        istat_sp = 2
        spinges = spinold
      else
c       Spin is not converged
        if( istep_sp .eq. nstep_sp ) istat_sp = 1
c
        if( istep_sp .eq. 1 .or. meth_sp .eq. 0 )then
c         Use simple blend
          spinsmpl = spinold + spindel*blend_sp
          spinges = spinsmpl
          spinslope = blend_sp
        elseif( istep_sp .gt. 1 )then
c         Use linear fit for guess, unless numerical issues
          spinold1 = spindata(1,1)
          spinnew1 = spindata(2,1)
          spindel1 = spinnew1 - spinold1
          spinold2 = spinold
          spinnew2 = spinnew
          spindel2 = spindel
          spinfac = (spinold1-spinold2) / (spindel2-spindel1)
          if( spinfac .lt. zero )then
c           Then this is not convergent, revert to simple blend
            spinslope = blend_sp
          elseif( spinfac .gt. two )then
c           This looks dangerously large, limit slope
            spinslope = two
          else
c           We have a reasonable slope, intersect the two lines
            spinslope = spinfac
          endif
          spinges = spinold2 + spindel2*spinslope
c
c         Shift up data from previous last step
          spindata(1,2) = spindata(1,1)
          spindata(2,2) = spindata(2,1)
          spindata(3,2) = spindata(3,1)
          spindata(4,2) = spindata(4,1)
        else
c         Cannot be, and check is free
          call STOPXERR( 'istep_sp=0 in SPUPDATE' )
        endif
c
        if( spinges .lt. zero )then
c         We overshot from above into a (spurious) negative spin, stifle
          spinges = spinnew
          spinslope = one
c         If by some numerical noise, spin<0, set to zero:
          if( spinges .lt. zero ) spinges = zero
        endif
c
c       Pack up this step:
        spindata(1,1) = spinold
        spindata(2,1) = spinnew
        spindata(3,1) = spinges
        spindata(4,1) = spinslope
c
c       Set up electron counts for next spin step:
        elecno(2) = ( elnumbr - spinges ) / two
        elecno(1) = elecno(2) + spinges
      endif
c
c    That's all, Folks!
c
      RETURN
      END
