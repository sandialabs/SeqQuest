c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FIXFRAME
c
c
      subroutine FIXFRAME( natm, iatframe, ratm, fatm )
c---------------------------------------------------------------
c Purpose: Adjust the forces on a geometry in order to fix
c          the frame of reference.
c
c Written: Peter A. Schultz, 28-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input: frame-defining atoms, and atomic coordinates
      DIMENSION  iatframe(3)
      DIMENSION  ratm(3,natm)
c Input/Output: forces on atoms
      DIMENSION  fatm(3,natm)
c
c Local:
      DIMENSION    vecij(3), vecjk(3), vecijk(3)
      DATA         zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Identify the atoms which define the frame
c   iatmfix = is the atom which is fixed in place
c   iatmvec = atom fixed along vector from iatmfix
c   iatmpln = with iatmfix and iatmvec, defines fixed plane of atoms
c             (if iatmpln is negative, enforce collinearity)
c
      iatmfix = iatframe(1)
      if( iatmfix.eq.0 ) RETURN
      iatmvec = iatframe(2)
      iatmpln = iatframe(3)
c
c Enforce fixed atom by zeroing forces on iatmfix:
c
      do  i=1,3
        fatm(i,iatmfix) = zero
      enddo
c
c Enforce vector of iatmvec from iatmfix ...
c
c  ... compute normalized vector connecting atoms iatmfix amd iatmvec
      vecnrm = zero
      do  i=1,3
        vecij(i) = ratm(i,iatmvec) - ratm(i,iatmfix)
        vecnrm = vecnrm + vecij(i)**2
      enddo
      vecnrm = one / SQRT( vecnrm )
      do  i=1,3
        vecij(i) = vecnrm*vecij(i)
      enddo
c  ... project iatmvec force along vector
      fdotv = zero
      do  i=1,3
        fdotv = fdotv + fatm(i,iatmvec)*vecij(i)
      enddo
c  ... and keep only force along this vector
      do  i=1,3
        fatm(i,iatmvec) = fdotv*vecij(i)
      enddo
c
c Enforce iatmfix-iatmvec-iatmpln constraint
c
      if( iatmpln.gt.0 )then
c
c       Enforce plane ...
c
c        ... compute iatmvec-iatmpln normalized vector
        vecnrm = zero
        do  i=1,3
          vecjk(i) = ratm(i,iatmpln) - ratm(i,iatmvec)
          vecnrm = vecnrm + vecjk(i)**2
        enddo
        vecnrm = one / SQRT( vecnrm )
        do  i=1,3
          vecjk(i) = vecnrm*vecjk(i)
        enddo
c        ... compute normalized normal to plane
        call CROSS(  vecij,vecjk, vecijk )
        vecnrm = zero
        do  i=1,3
          vecnrm = vecnrm + vecijk(i)**2
        enddo
        vecnrm = one / SQRT( vecnrm )
        do  i=1,3
          vecijk(i) = vecnrm*vecijk(i)
        enddo
c        ... project iatmpln force along normal to plane ...
        fdotv = zero
        do  i=1,3
          fdotv = fdotv + fatm(i,iatmpln)*vecijk(i)
        enddo
c        ... and remove from force vector on this atom:
        do  i=1,3
          fatm(i,iatmpln) = fatm(i,iatmpln) - fdotv*vecijk(i)
        enddo
c
      elseif( iatmpln.lt.0 )then
c
c       Enforce collinearity of third atom with first two ...
c
        jatmpln = -iatmpln
c        ... project iatmpln force along iatmfix-iatmvec vector
        fdotv = zero
        do  i=1,3
          fdotv = fdotv + fatm(i,jatmpln)*vecij(i)
        enddo
c        ... and keep only that component of the force:
        do  i=1,3
          fatm(i,jatmpln) = fdotv*vecij(i)
        enddo
c
c     else
c       Nothing to enforce.
      endif
c
c    That's all Folks!
c
      RETURN
      END
