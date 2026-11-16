c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EBLENDASD
c
c
      subroutine EBLENDASD( IWR, sufx, itstep, memblnde,
     $ nrlx, dt,
     $ fi,ri, vi,wt, savblnd )
c---------------------------------------------------------------
c Purpose: Do ASD-type SCF blend
c
c Written: Peter A. Schultz, Decmeber-2014, for 2.66
c
c Revision History:
c   14Dec14-PAS/2.66: adapted from DMD routine
c---------------------------------------------------------------
c
C ***** IMPORTANT: THIS DOES NOT WORK YET *****
c Variables:
c  nrlx = total number of independent components to be relaxed
c   (this will usually be 3-vectors, but *could* be potentially
c    less if one introduces constraints on the geometry, e.g.
c    fixing the z-coordinate of an atom to be a fixed number)
c  ri(),fi() = (input) positions, gradients
c  vi()      = (scratch) velocities
c  wt() = inertia factor for velocity update
c  savblnd() = saved velocity and eight factors (should overlap vi and wt <sigh>)
c
c  NB: output ri() is copied into vi() tor compatibility with broyden blend
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input/output variables
      DIMENSION  ri(*),vi(*),fi(*), wt(*)
c file stuff:
      CHARACTER  sufx*(*), dmd_tag*40
c scratch declaration:
      DIMENSION  savblnd(*)
c
c local declarations
      DATA  lstout / 2 /
      DATA  zero,one / 0.d0,1.d0 /
      DATA  two / 2.d0 /
c Defines blending type
      DATA    methbl_broy,methbl_asd / 1,2 /
c
c >>>> EXECUTABLE CODE:
c
c To fix to ASD blend
      methbl = methbl_asd
c
      dmd_tag = 'dmd_' // sufx
      if( itstep.eq.1 )then
c       Initialize velocities and inertia factors on first step
        call MKZERO( nrlx, vi )
        call DCOPY( nrlx , one,0, wt,1 )
      else
c       Recover velocities, and inertia factors (if applicable):
        if( memblnde .lt. 2*nrlx )then
          call FLOPENB( idmdfl, dmd_tag )
          REWIND( unit=idmdfl )
          call READBIG( idmdfl, nrlx, vi )
          call READBIG( idmdfl, nrlx , wt )
        else
          call DCOPY( nrlx, savblnd(1),1, vi,1 )
          call DCOPY( nrlx, savblnd(nrlx+1),1, wt,1 )
        endif
      endif
c
c Update velocity with specified quenching scheme:
      call EASDVELOC( nrlx, vi,fi,wt, dt )
c
c Update the atomic coordinates along the velocity vector:
      call DMDSTEP( nrlx, ri,vi, dt )
c
c And output some diagnostics:
      if( lstout.gt.0 )
     $ write(IWR,*) 'e-asd blend step #', itstep
c
      if( lstout.gt.1 )then
        iwtmn = 1
        iwtmx = 1
        wtmin = ABS( wt(1) )
        wtmax = ABS( wt(1) )
        do  irlx=1,nrlx
          wti = ABS( wt(irlx) )
            if( wti .lt. wtmin )then
            wtmin = wti
            iwtmn = irlx
          endif
          if( wti .gt. wtmax )then
            wtmax = wti
            iwtmx = irlx
          endif
        enddo
        write(IWR,'(a,f12.6,1x,i12)')'EASD: min mass fac=',wtmin,iwtmn
        write(IWR,'(a,f12.6,1x,i12)')'EASD: max mass fac=',wtmax,iwtmx
c
        ivvmn = 1
        ivvmx = 1
        vvmin = ABS( vi(1) )
        vvmax = ABS( vi(1) )
        do  irlx=1,nrlx
          vvi = ABS( vi(irlx) )
          if( vvi .lt. vvmin )then
            vvmin = vvi
            ivvmn = irlx
          endif
          if( vvi .gt. vvmax )then
            vvmax = vvi
            ivvmx = irlx
          endif
        enddo
        write(IWR,'(a,f12.6,1x,i12)')'EASD: min velocity=',vvmin,ivvmn
        write(IWR,'(a,f12.6,1x,i12)')'EASD: max velocity=',vvmax,ivvmx
      endif
c
c Put out the current velocity vector:
      if( memblnde .lt. 2*nrlx )then
        REWIND( unit=idmdfl )
        call WRITBIG( idmdfl, nrlx, vi )
        call WRITBIG( idmdfl, nrlx, wt )
      else
        call DCOPY( nrlx, vi,1, savblnd(1),1 )
        call DCOPY( nrlx, wt,1, savblnd(nrlx+1),1 )
      endif
c
      if( memblnde .lt. 2*nrlx )
     $call FLCLOSE( idmdfl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EASDVELOC
c
c
      subroutine EASDVELOC( nrlx, vi,fi,wt, dt )
c---------------------------------------------------------------
c Purpose: (damped) dynamics update of particle velocities
c
c Written: Peter A. Schultz, 12-October-2001, for v2.50
c
c Revision History:
c  12Oct01-PAS/2.50: developed new damped dynamics scheme
c---------------------------------------------------------------
c
c Parameters for modified damped dynamics:
c  dvup = boosting factor for rate of acceleration
c         must be .ge.1, should be .le.dvdn
c  dvdn = reducing factor (if rate deemed too fast)
c         must be .ge.1, should be between 2 and 4(?)
c  cutdn = v-par/dv-update ratio that triggers deceleration
c         must be .ge.1 (if dv bigger than v-par, must slow down)
c         should be ~1.5-3.0 (dv and v-par comparable, should slow)
c  Setting dvdn=dvup=1.0 reduces to simple Jonsson damped dynamics.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0 )
c input/output arrays:
      DIMENSION  vi(*),fi(*), wt(*)
c
c >>>> EXECUTABLE CODE:
c
c Set parameters for the damped acceleration dynamics ...
c  ... the acceleration factor (cube root of 2.0):
      dvup = 1.25992105d0
c      dvup = 1.414213562d0
c  ... the deceleration factor (2.0 seems good number):
      dvdn = 2.d0
c  ... and cutoff ratio for deceleration (2.0-3.0 appears to work):
      cutdn = 2.d0
c  ... and limit weight boosting, and velocity update
      wtlimit = 2.d0/dt
      filimit = 2.0d0
cxxx:
      call FLGETIWR( IWR )
      write(IWR,*) 'DBG/EASDVELOC: dt,filimit=',dt,filimit
c
c       >>>> accel. steepest desc.
c
        do  irlx =1,nrlx
c
c         Get size of force&velocity vectors, and dot product
          vdotf = vi(irlx)*fi(irlx)
          fdotf = fi(irlx)*fi(irlx)
          vdotv = vi(irlx)*vi(irlx)
c
          if( vdotf .lt. zero )then
c           Velocity and force reversed: quench, and update
c
            dvfac = ABS( wt(irlx) )
            dvfacx = dvfac
c
            dvmag = (fi(irlx)* dvfac*dt)**2
            fac = vdotf/fdotf
            vpar = fac*fac*fdotf
            dvmag = SQRT( dvmag )
            vpar = SQRT( vpar )
c
c           If we have a big reversal ...
            if( cutdn*dvmag .gt. vpar .or.
c              ... or reversed twice in a row ...
     $                     wt(irlx) .gt. zero )then
c              ... decrease acceleration rate:
              dvfac = dvfac / dvdn
              dvfacx = dvfac
c
c             If we still overshoot after damping ...
              if( (dvmag/dvdn) .gt. vpar )then
c                ... limit update size to previous step size:
                dvfacx = dvfac*( vpar / (dvmag/dvdn) )
              endif
c
            endif
c
            vi(irlx) = fi(irlx)* dvfacx*dt
c
c           We just reversed, set as boost-ineligible:
            wt(irlx) = dvfac
c
          else
c           Velocity and force aligned, update velocity along force
c
            dvfac = ABS( wt(irlx) )
c
            if( fdotf .ne. zero )then
              fac = vdotf/fdotf
              dvmag = (fi(irlx)* dvfac*dt)**2
              vpar = fac*fac*fdotf
              dvmag = SQRT( dvmag )
              vpar = SQRT( vpar )
            else
              fac = zero
              dvmag = zero
              vpar = zero
            endif
c
            if( wt(irlx) .lt. zero )then
c             We are repeatedly heading in the right direction ...
c
cpas: probably should do unconditional boost.
              if( two*dvmag .lt. vpar )then
c                ... and not accelerating fast enough - boost
                dvfac = dvfac * dvup
              endif
c             Note boost-eligibility in sign of mass factor:
c               ... after limiting weighting:
              if( dvfac .gt. wtlimit ) dvfac = wtlimit
              wt(irlx) = -dvfac
c
            else
c             We recently (but not this time) reversed
c
              if( two*dvmag .lt. vpar )then
c               Looks like we are converging, quench velocity ...
                fac = zero
c                ... and leave this as boost-ineligible
                wt(irlx) = dvfac
              else
c               We have healthly acceleration ...
                if( vpar .gt. zero )then
c                  ... and constructive (non-zero) alignment
                  if( two*dvmag .gt. vpar )then
c                    ... and the acceleration is hefty - boost
c                    (even though we were nominally ineligible)
                    dvfac = dvfac * dvup
                  endif
c                 Now boost-eligible, note in sign of mass factor:
                   if( dvfac .gt. wtlimit ) dvfac = wtlimit
c               ... after limiting weighting:
                  wt(irlx) = -dvfac
                endif
              endif
c
            endif
c
C            vi(irlx) = fac*fi(irlx) + fi(irlx)*dvfac*dt
            facfi = fac + dvfac*dt
            if( facfi .gt. filimit )then
c             Jumpino too far, limit step, and unboost eligible
              wt(irlx) = ABS( wt(irlx) )
              facfi = filimit
            endif
            vi(irlx) = facfi*fi(irlx)
c
          endif
c
        enddo
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
c End
