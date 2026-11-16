c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMRMAT
c
c
      subroutine SYMRMAT( nsymi, isymtyp,syminfo,
     $ nsymd, nsym, rmatsym,  rprim )
c---------------------------------------------------------------
c Purpose: generate cartesian reps of 3D symmetry operations
c
c Written: Peter A. Schultz, 6-January-1998, for v2.30
c          Adapted from symr3mat, and also replacing symr2mat
c
c Revision history:
c  21Jun01-PAS/2.48: replace STOPs, writes
c  21Jul99-PAS/2.37: project translation vector+right handedness
c  26Apr99-PAS/2.34: fix rotation+inversion generation
c---------------------------------------------------------------
c
c  Variables:
c   nsymi - number of input non-identity symmetry definitions
c   nsym  - number of output non-identity symmetry operations
c   isymtyp() - type of symmetry operation
c     order of rotation(+inversion if less than 0)
c     Valid inputs: -6,-4,-3,-2,-1,1,2,3,4,6
c     NB: 1=reflection equivalent to "-2" 2-fold rotation+inv
c   syminfo() - defining axis of those operations (normalized)
c     syminfo(1-3) = normalized axis of rotation (Cartesian)
c     syminfo(4-6) = translation vector (lattice units)
c   rmatsym() - output matrix of cartesian reps of symmetries
c   rprim()   - primitive lattice vectors
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
c input arrays:
      DIMENSION  rprim(3,3)
      DIMENSION  isymtyp(*), syminfo(6,*)
c output array:
      DIMENSION  rmatsym(3,4,*)
c local declarations:
      LOGICAL  doinvrs
      DIMENSION  u(3),off(3)
      DATA  zero,one,two / 0.d0,1.d0,2.d0 /
c
c >>>> EXECUTABLE CODE:
c
      nsym = 0
      if( nsymi.lt.1 ) RETURN
c
      isym = 0
      do 100 isymi=1,nsymi
c
        doinvrs = .false.
        istyp = isymtyp(isymi)
c
c  Load rotation axis:
        u(1) = syminfo(1,isymi)
        u(2) = syminfo(2,isymi)
        u(3) = syminfo(3,isymi)
c
c  Project translation vector into primary unit cell:
c  (thanks to Dave Raczkowski for pointing out that
c   translation vectors of the input generators need to be
c   projected into the primary unit cell, or the group
c   generation code might fail - 21Jul99-PAS)
        do 20 j=4,6
c  I want the lattice projection to be 0.le.x.lt.1
          tr = syminfo(j,isymi)
          if( tr.ge.one )then
            itr = tr
            ishift = -itr
          elseif( tr.lt.zero )then
            trm = -tr + 0.999
            ishift = trm
          else
            goto 20
          endif
          syminfo(j,isymi) = syminfo(j,isymi) + ishift
   20   continue
c
c  Generate translation vector:
        do  j=1,3
          off(j) = syminfo(4,isymi)*rprim(j,1)
     $           + syminfo(5,isymi)*rprim(j,2)
     $           + syminfo(6,isymi)*rprim(j,3)
        enddo
        offsq = off(1)**2 + off(2)**2 + off(3)**2
c
        rotang = zero
c I will generate all rotations of a generator ...
        nrepeat = IABS(istyp) - 1
c  ... unless it has a translation, then let group code do it elsewhere
        if( offsq .gt. 1.d-4 ) nrepeat = 1
c
   40   isym = isym + 1
        if( isym.gt.nsymd ) goto 90
c
c Load translation vector:
c
        rmatsym(1,4,isym) = off(1)
        rmatsym(2,4,isym) = off(2)
        rmatsym(3,4,isym) = off(3)
c
c Compute transformation matrices:
c
        if( istyp .eq. 1 )then
c
c         REFLECTION
          do  j=1,3
            do  i=1,3
              rmatsym(i,j,isym) = -two*u(i)*u(j)
            enddo
          enddo
          rmatsym(1,1,isym) = rmatsym(1,1,isym) + one
          rmatsym(2,2,isym) = rmatsym(2,2,isym) + one
          rmatsym(3,3,isym) = rmatsym(3,3,isym) + one
c
        elseif( IABS(istyp) .gt. 1 )then
c
c         ROTATION
          nrot = IABS(istyp)
          rotang = rotang + two*pi/DBLE(nrot)
          cosa = COS(rotang)
          sina = SIN(rotang)
          fac = one - cosa
          do  j=1,3
            do  i=1,3
              rmatsym(i,j,isym) = rmatsym(i,j,isym) + fac*u(i)*u(j)
            enddo
          enddo
          rmatsym(1,1,isym) = rmatsym(1,1,isym) + cosa
          rmatsym(2,2,isym) = rmatsym(2,2,isym) + cosa
          rmatsym(3,3,isym) = rmatsym(3,3,isym) + cosa
          rmatsym(1,2,isym) = rmatsym(1,2,isym) + sina*u(3)
          rmatsym(2,1,isym) = rmatsym(2,1,isym) - sina*u(3)
          rmatsym(2,3,isym) = rmatsym(2,3,isym) + sina*u(1)
          rmatsym(3,2,isym) = rmatsym(3,2,isym) - sina*u(1)
          rmatsym(3,1,isym) = rmatsym(3,1,isym) + sina*u(2)
          rmatsym(1,3,isym) = rmatsym(1,3,isym) - sina*u(2)
c
c decide whether to toggle inversion:
          if( istyp.lt.0 ) doinvrs = .not. doinvrs
c
        elseif( istyp .eq. -1 .or. istyp.eq.0 )then
c
          doinvrs = .true.
c         INVERSION - wait, install identity
          do  j=1,3
            do  i=1,3
              rmatsym(i,j,isym) = zero
            enddo
          enddo
          rmatsym(1,1,isym) = one
          rmatsym(2,2,isym) = one
          rmatsym(3,3,isym) = one
c
        endif
c
        if( doinvrs )then
c         INVERSION - actually do it
          do  j=1,3
            do  i=1,3
              rmatsym(i,j,isym) = -rmatsym(i,j,isym)
            enddo
          enddo
        endif
c
   90   continue
        nrepeat = nrepeat - 1
c Go back and generate next rotation in this set:
        if( nrepeat .gt. 0 ) goto 40
c
  100 continue
c
c  Record number of symmetry operations generated here:
c
      nsym = isym
      if( nsym.gt.nsymd )then
        call FLGETIERR( IERRFL )
        write(IERRFL,*) '>>>>> total so far=',nsym,', max=',nsymd
        call STOPXERR( 'symrmatd/ too many sym ops in group' )
      endif
c
c    That's all Folks!
c
      RETURN
      END
