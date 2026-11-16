c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ROTMAT
c
c
      subroutine ROTMAT( nsymi,nsyms, isymtyp,syminfo, dsym )
c---------------------------------------------------------------
c Purpose: compute inverse rotation and reflection matrices,
c          for two-dimensional symmetries.
c
c Revision history:
c  14Jun05-PAS/2.59: File units
c  21Jun01-PAS/2.48: replace STOPs
c   9Jan99-PAS/2.30: redo symmetry scheme
c  10Dec98-PAS/2.29: promote source to explicit double precision
c   3Sep96-PAS/2.15: fix mkzero dimensions
c---------------------------------------------------------------
c
c Note: the following code adapted from PJF's original, which
c was restricted to 2D symmetries and used symmetry info saved
c in different format.  To keep the working code here, have to
c back out the original info from the current data scheme.
c At some point, it would be a good idea to retool this routine
c to use the cartesian symmetry matrices, and maybe even
c generalize to 3D calculations, but not today. - 9Jan99-PAS
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c input symmetry angles:
      DIMENSION  isymtyp(*), syminfo(6,*)
c output symmetry transformation matrices:
      DIMENSION  dsym(5,5,3,*)
c local declarations:
      DATA  zero,one,two / 0.d0,1.d0,2.d0 /
c
c >>>> EXECUTABLE CODE:
c
      radians = ASIN(one)/90.d0
      call MKZERO( 5*5*3*nsyms, dsym )
c
c     IDENTITY (isym=1)
c
      isym = 1
      do  l=1,3
        do  li=1,5
          dsym(li,li,l,isym) = one
        enddo
      enddo
c
      do 100 isymi=1,nsymi
c
c Need to back out what the operation is:
        istyp = isymtyp(isymi)
        if( istyp.gt.1 )then
c         ROTATION
          dang = two*pi / DBLE( istyp )
          ang = zero
          nrepeat = istyp - 1
        elseif( istyp.eq.-2 )then
c         REFLECTION
c         Note: the angle is based on indata assumption
          sy = -syminfo(1,isymi)
          sx =  syminfo(2,isymi)
          smag = SQRT( sx*sx + sy*sy )
          if( sy.lt.zero )then
            sy = -sy
            sx = -sx
          endif
          arg = sx/smag
c          arg = syminfo(2,isymi)
          ang = ACOS( arg )
          nrepeat = 1
        else
c         ERROR
          call STOPXERR( 'rotmat2d' )
        endif
c
   50   if( nrepeat.le.0 ) goto 100
        nrepeat = nrepeat - 1
        isym = isym + 1
c
        if( istyp.gt.1 )then
c         ROTATION
          ang = ang + dang
          cosa = COS(ang)
          sina = SIN(ang)
          cos2a = (cosa-sina)*(cosa+sina)
          sin2a = two*cosa*sina
c l=0 matrix
          dsym(1,1,1,isym) =  one
c l=1 matrix
          dsym(1,1,2,isym) =  cosa
          dsym(2,2,2,isym) =  cosa
          dsym(1,2,2,isym) = -sina
          dsym(2,1,2,isym) =  sina
          dsym(3,3,2,isym) =  one
c l=2 matrix
          dsym(1,1,3,isym) =  cos2a
          dsym(2,2,3,isym) =  cosa
          dsym(3,3,3,isym) =  cosa
          dsym(4,4,3,isym) =  cos2a
          dsym(5,5,3,isym) =  one
          dsym(1,4,3,isym) =  sin2a
          dsym(2,3,3,isym) =  sina
          dsym(3,2,3,isym) = -sina
          dsym(4,1,3,isym) = -sin2a
        elseif( istyp.eq.-2 )then
c         REFLECTION
          cos2a = COS(ang+ang)
          sin2a = SIN(ang+ang)
          cos4a = (cos2a-sin2a)*(cos2a+sin2a)
          sin4a = two*cos2a*sin2a
c l=0 matrix
          dsym(1,1,1,isym) =  one
c l=1 matrix
          dsym(1,1,2,isym) =  cos2a
          dsym(2,2,2,isym) = -cos2a
          dsym(1,2,2,isym) =  sin2a
          dsym(2,1,2,isym) =  sin2a
          dsym(3,3,2,isym) =  one
c l=2 matrix
          dsym(1,1,3,isym) = -cos4a
          dsym(2,2,3,isym) = -cos2a
          dsym(3,3,3,isym) =  cos2a
          dsym(4,4,3,isym) =  cos4a
          dsym(5,5,3,isym) =  one
          dsym(1,4,3,isym) =  sin4a
          dsym(2,3,3,isym) =  sin2a
          dsym(3,2,3,isym) =  sin2a
          dsym(4,1,3,isym) =  sin4a
        endif
c
        goto 50
c
  100 continue
c
      if( isym.ne.nsyms )then
        call FLGETIERR( IERRFL )
        write(IERRFL,*) '>>>>> must input all reflections for 2D'
        write(IERRFL,*) '>>>>> nsym(rotmat),nsym(code)=',isym,nsyms
        call STOPXERR( '2d-syms - 2D symmetry input incomplete' )
      endif
c
      RETURN
      END
