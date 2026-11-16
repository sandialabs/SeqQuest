c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RINTRNL
c
c
      subroutine RINTRNL( ndim,nvec, dolvatm, rprim,orig,rscale,
     $ rextr,rintr )
c---------------------------------------------------------------
c Purpose: take atom coordinates, shift, scale to/from user input
c
c Written: Peter A. Schultz, 24-September-1997, v2.21
c
c Revision history:
c   7Jul00-PAS/2.45: unify scaling functions; lattice coords
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c Note: For lattice coordinates, I *enforce* the assumption,
c       *not* checked, that for reduced dimension systems,
c       non-periodic axes are in special directions.
c       E.g., slab normal in a 2D calculation is the z-axis
c Note: I might consider pulling ruc() up
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input declarations:
      LOGICAL    dolvatm
      DIMENSION  rprim(3,4),orig(3),rscale(3)
c input/output arrays:
      DIMENSION  rextr(3,nvec),rintr(3,nvec)
c local declarations:
      DIMENSION  ruc(3,3), vect(3)
      DATA  zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      if( dolvatm )then
c
c       Use reduced dimension unit cell vectors:
        do  j=1,3
          do  i=1,3
            ruc(i,j) = zero
          enddo
          ruc(j,j) = one
        enddo
        do  j=1,ndim
          do  i=1,ndim
            ruc(i,j) = rprim(i,j) / rprim(j,4)
          enddo
        enddo
c       I choose to install scaling factor here.
        do  i=ndim+1,3
          ruc(i,i) = ruc(i,i) * rscale(i)
        enddo
c
c       Expand the vectors into coordinate space:
        do  ivec=1,nvec
          call LVXPND( ruc, rintr(1,ivec), rextr(1,ivec) )
        enddo
c
      else
c
        do  ivec=1,nvec
          rintr(1,ivec) = rextr(1,ivec)*rscale(1)
          rintr(2,ivec) = rextr(2,ivec)*rscale(2)
          rintr(3,ivec) = rextr(3,ivec)*rscale(3)
        enddo
c
      endif
c
      do  ivec=1,nvec
        rintr(1,ivec) = rintr(1,ivec) - orig(1)
        rintr(2,ivec) = rintr(2,ivec) - orig(2)
        rintr(3,ivec) = rintr(3,ivec) - orig(3)
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> REXTRNL
c
c
      entry REXTRNL( ndim,nvec, dolvatm, rprim,orig,rscale,
     $ rextr,rintr )
c
c >>>> EXECUTABLE CODE:
c
      do  ivec=1,nvec
        rextr(1,ivec) = rintr(1,ivec) + orig(1)
        rextr(2,ivec) = rintr(2,ivec) + orig(2)
        rextr(3,ivec) = rintr(3,ivec) + orig(3)
      enddo
c
      if( dolvatm )then
c
c       Use reduced dimension unit cell vectors:
        do  j=1,3
          do  i=1,3
            ruc(i,j) = zero
          enddo
          ruc(j,j) = one
        enddo
        do  j=1,ndim
          do  i=1,ndim
            ruc(i,j) = rprim(i,j) / rprim(j,4)
          enddo
        enddo
c       I choose to install scaling factor here.
        do  i=ndim+1,3
          ruc(i,i) = ruc(i,i) * rscale(i)
        enddo
c
c       Project coordinates onto reduced unit cell vectors:
        do  ivec=1,nvec
          call LVPROJ( ruc, rextr(1,ivec), vect )
          rextr(1,ivec) = vect(1)
          rextr(2,ivec) = vect(2)
          rextr(3,ivec) = vect(3)
        enddo
c
      else
c
        do  ivec=1,nvec
          rextr(1,ivec) = rextr(1,ivec)/rscale(1)
          rextr(2,ivec) = rextr(2,ivec)/rscale(2)
          rextr(3,ivec) = rextr(3,ivec)/rscale(3)
        enddo
c
      endif
c
      RETURN
c
      END
