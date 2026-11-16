c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMVEC
c
c
      subroutine SYMVEC( nvec, vec, rsym, nsym, rmatsym )
c---------------------------------------------------------------
c Purpose: symmetrize vectors
c
c Written: Peter A. Schultz, 28-September-1999, v2.40
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input/output arrays:
      DIMENSION  vec(3,nvec)
      DIMENSION  rsym(3), rmatsym(3,4,*)
c local declarations
      DIMENSION  vi(3),svi(3),vis(3)
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      if( nsym.le.0 ) RETURN
c
c  Loop over all vectors, accumulate symmetry vector
c
      symscal = one/ DBLE( nsym+1 )
      do 300 ivec=1,nvec
        vi(1) = vec(1,ivec) - rsym(1)
        vi(2) = vec(2,ivec) - rsym(2)
        vi(3) = vec(3,ivec) - rsym(3)
        vis(1) = vi(1)
        vis(2) = vi(2)
        vis(3) = vi(3)
        do  isym=1,nsym
c
          svi(1) = rmatsym(1,1,isym)*vi(1) + rmatsym(2,1,isym)*vi(2)
     $           + rmatsym(3,1,isym)*vi(3)
          svi(2) = rmatsym(1,2,isym)*vi(1) + rmatsym(2,2,isym)*vi(2)
     $           + rmatsym(3,2,isym)*vi(3)
          svi(3) = rmatsym(1,3,isym)*vi(1) + rmatsym(2,3,isym)*vi(2)
     $           + rmatsym(3,3,isym)*vi(3)
c
          vis(1) = vis(1) + svi(1)
          vis(2) = vis(2) + svi(2)
          vis(3) = vis(3) + svi(3)
        enddo
        vec(1,ivec) = rsym(1) + symscal*vis(1)
        vec(2,ivec) = rsym(2) + symscal*vis(2)
        vec(3,ivec) = rsym(3) + symscal*vis(3)
c
  300 continue
c
c    That's all Folks!
c
      RETURN
      END
