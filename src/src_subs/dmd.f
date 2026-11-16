c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DMD
c
c
      subroutine DMD( IWR, sufx, igstep, igdmdmod,
     $ nrlx, natmg,nxatm, dt,
     $ am, ri,vi,fi, wt )
c---------------------------------------------------------------
c Purpose: damped dynamics update of particle positions
c
c Written: Aidan Thompson (SNL), September 2001
c
c Revision History:
c   28Aug28-PAS/2.61f: set ASD atomic masses to one (vary just inertial)
c    4Feb02-PAS/2.51c: added restart (v,m) file
c   11Oct01-PAS/2.50: adapted to SeqQuest conventions; more options
c---------------------------------------------------------------
c
c Variables:
c  nrlx = total number of independent components to be relaxed
c  natmg = number of "atoms" to be relaxed
c  nxatm = number of components to be relaxed for each atom
c   (this will usually be 3-vectors, but *could* be potentially
c    less if one introduces constraints on the geometry, e.g.
c    fixing the z-coordinate of an atom to be a fixed number)
c  am() = masses for each of those components (NB: not by atom)
c  ri(),vi(),fi() = positions, velocities, and forces
c  wt() = inertia factor for velocity update
c
c igdmdmod  function
c    0      simple unquenched MD step
c    1      Jonsson's damped dynamics (component-wise)
c    2      Jonsson's damped dynamics (atom-by-atom)
c    3      PAS-modified damped dynamics: accelerated steepest descent
c   Split into unquenched and quenched if igdmdmod < 0
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input/output variables
      DIMENSION  am(*), ri(*),vi(*),fi(*), wt(*)
      DIMENSION  nxatm(*)
c file stuff:
      CHARACTER  sufx*(*), dmd_tag*40
c
c local declarations
      DATA  lstout / 2 /
      DATA  zero,one / 0.d0,1.d0 /
      DATA  two / 2.d0 /
c
c >>>> EXECUTABLE CODE:
c
      dmd_tag = 'dmd_' // sufx
      call FLOPENB( idmdfl, dmd_tag )
      REWIND( unit=idmdfl )
      if( igstep.eq.1 )then
c       Initialize velocities and inertia factors on first step
        call MKZERO( nrlx, vi )
        call DCOPY( natmg, one,0, wt,1 )
      else
c       Recover velocities, and inertia factors (if applicable):
        call READBIG( idmdfl, nrlx, vi )
        if( IABS( igdmdmod ) .ge. 3 )then
          call READBIG( idmdfl, natmg, wt )
        else
          call DCOPY( natmg, one,0, wt,1 )
        endif
      endif
c
      iquench = igdmdmod
      dtvel = dt
c
      if( igdmdmod .lt. 0 )then
c       Split the time interval into unquenched and quenched ...
        iquench = 0
        dtvel = dt / two
c        ... the unquenched velocity update in the first half-step:
        call DMDVELOC( nrlx, natmg,nxatm, vi,am,fi,wt, dtvel, iquench )
c
c        ... set up specified quench for second half of time-step:
        iquench = -igdmdmod
      endif
c
c Update velocity with specified quenching scheme:
      call DMDVELOC( nrlx, natmg,nxatm, vi,am,fi,wt, dtvel, iquench )
c
c Update the atomic coordinates along the velocity vector:
      call DMDSTEP( nrlx, ri,vi, dt )
c
c And output some diagnostics:
      if( lstout.gt.0 )
     $ write(IWR,*) 'dynamics update step #', igstep
c
c Put out the current velocity vector:
      REWIND( unit=idmdfl )
      call WRITBIG( idmdfl, nrlx, vi )
      if( IABS( igdmdmod ) .ge. 3 )then
        call WRITBIG( idmdfl, nrlx, wt )
        wtmin = ABS( wt(1) )
        wtmax = ABS( wt(1) )
        do  iatmg=1,natmg
          wti = ABS( wt(iatmg) )
          if( wti .lt. wtmin ) wtmin = wti
          if( wti .gt. wtmax ) wtmax = wti
        enddo
        if( lstout.gt.1 )
     $  write(IWR,'(a,2f9.4)') 'ASD: min/max mass factors=',wtmin,wtmax
c
        if( lstout.gt.2 )then
 9013     format(1x,'ASDWT:',2i4,2x,8f9.4)
          do  iatmg=1,natmg,8
            jatmg = MIN( iatmg+7, natmg )
            write(IWR,9013) iatmg,jatmg, (wt(i),i=iatmg,jatmg)
          enddo
        endif
      endif
c
      call FLCLOSE( idmdfl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DMDVELOC
c
c
      subroutine DMDVELOC( nrlx, natmg,nxatm, vi,am,fi,wt, dt, iquench )
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
      DIMENSION  am(*), vi(*),fi(*), wt(*)
      DIMENSION  nxatm(*)
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
cpas: set masses to **one** for ASD and related (this needs to cleaned up)
      if( iquench.gt.0 )then
c       set all mass elements to one
        do  irlx=1,nrlx
          am(irlx) = one
        enddo
      endif
c
      if( iquench .eq. 0 )then
c
c       >>>> No quench, strict dynamics update
c
        do  irlx=1,nrlx
          vi(irlx) = vi(irlx) + fi(irlx)* dt/am(irlx)
        enddo
c
      elseif( iquench .eq. 1 )then
c
c       >>>> Jonsson's damped dynamics (component-wise)
c
        do  irlx=1,nrlx
          if(  vi(irlx)*fi(irlx) .lt. zero )then
c           Quench velocity, and update
            vi(irlx) =            fi(irlx)* dt/am(irlx)
          else
c           Add to velocity ...
            vi(irlx) = vi(irlx) + fi(irlx)* dt/am(irlx)
          endif
        enddo
c
      elseif( iquench .eq. 2 )then
c
c       >>>> Jonsson's damped dynamics (atom-by-atom)
c
        irlx0 = 0
        do  iatmg=1,natmg
c         Get number of components on this atom to be relaxed
          lenvec = nxatm(iatmg)
c
          vdotf = zero
          fdotf = zero
          do  irlx=irlx0+1,irlx0+lenvec
            vdotf = vdotf + vi(irlx)*fi(irlx)
            fdotf = fdotf + fi(irlx)*fi(irlx)
          enddo
c
          if( vdotf .lt. zero )then
c           Velocity and force reversed: quench v, and update
            do  irlx=irlx0+1,irlx0+lenvec
              vi(irlx) =                fi(irlx)* dt/am(irlx)
            enddo
          else
c           Velocity and force aligned, update velocity along force
            if( fdotf .ne. zero )then
              fac = vdotf/fdotf
            else
              fac = zero
            endif
            do  irlx=irlx0+1,irlx0+lenvec
              vi(irlx) = fac*fi(irlx) + fi(irlx)* dt/am(irlx)
            enddo
          endif
c
        enddo
c
      elseif( iquench .eq. 3 )then
c
c       >>>> PAS modified damped dynamics, accel. steepest desc.
c
        irlx0 = 0
        do  iatmg=1,natmg
c         Get number of components on this atom to be relaxed
          lenvec = nxatm(iatmg)
c
c         Get size of force&velocity vectors, and dot product
          vdotf = zero
          fdotf = zero
          vdotv = zero
          do  irlx=irlx0+1,irlx0+lenvec
            vdotf = vdotf + vi(irlx)*fi(irlx)
            fdotf = fdotf + fi(irlx)*fi(irlx)
            vdotv = vdotv + vi(irlx)*vi(irlx)
          enddo
c
          if( vdotf .lt. zero )then
c           Velocity and force reversed: quench, and update
c
            dvfac = ABS( wt(iatmg) )
            dvfacx = dvfac
c
            dvmag = zero
            do  irlx=irlx0+1,irlx0+lenvec
              dvmag = dvmag + (fi(irlx)* dvfac*dt/am(irlx))**2
            enddo
            fac = vdotf/fdotf
            vpar = fac*fac*fdotf
            dvmag = SQRT( dvmag )
            vpar = SQRT( vpar )
c
c           If we have a big reversal ...
            if( cutdn*dvmag .gt. vpar .or.
c              ... or reversed twice in a row ...
     $                     wt(iatmg) .gt. zero )then
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
            do  irlx=irlx0+1,irlx0+lenvec
              vi(irlx) = fi(irlx)* dvfacx*dt/am(irlx)
            enddo
c
c           We just reversed, set as boost-ineligible:
            wt(iatmg) = dvfac
c
          else
c           Velocity and force aligned, update velocity along force
c
            dvfac = ABS( wt(iatmg) )
c
            if( fdotf .ne. zero )then
              fac = vdotf/fdotf
              dvmag = zero
              do  irlx=irlx0+1,irlx0+lenvec
                dvmag = dvmag + (fi(irlx)* dvfac*dt/am(irlx))**2
              enddo
              vpar = fac*fac*fdotf
              dvmag = SQRT( dvmag )
              vpar = SQRT( vpar )
            else
              fac = zero
              dvmag = zero
              vpar = zero
            endif
c
            if( wt(iatmg) .lt. zero )then
c             We are repeatedly heading in the right direction ...
c
cpas: probably should do unconditional boost.
              if( two*dvmag .lt. vpar )then
c                ... and not accelerating fast enough - boost
                dvfac = dvfac * dvup
              endif
c             Note boost-eligibility in sign of mass factor:
              wt(iatmg) = -dvfac
c
            else
c             We recently (but not this time) reversed
c
              if( two*dvmag .lt. vpar )then
c               Looks like we are converging, quench velocity ...
                fac = zero
c                ... and leave this as boost-ineligible
                wt(iatmg) = dvfac
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
                  wt(iatmg) = -dvfac
                endif
              endif
c
            endif
c
            do  irlx=irlx0+1,irlx0+lenvec
              vi(irlx) = fac*fi(irlx) + fi(irlx)*dvfac*dt/am(irlx)
            enddo
c
          endif
c
          irlx0 = irlx0 + lenvec
        enddo
c
      else
c
c       >>>> Invalid quenching option
c
        call STOPXERR( 'DMDVEL  - coding error/invalid iquench' )
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DMDSTEP
c
c
      subroutine DMDSTEP( nrlx, ri,vi, dt )
c
      implicit none
c
c     --- variables passed to and from subroutine
      DIMENSION  ri(*),vi(*)
      integer  nrlx
      real*8   ri,vi,   dt
c
c     --- local variables
      integer  irlx
c
c >>>> EXECUTABLE CODE:
c
      do  irlx=1,nrlx
        ri(irlx) = ri(irlx) + dt*vi(irlx)
      enddo
c
      RETURN
      END
