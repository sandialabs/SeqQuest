c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMMAP2
c
c
      subroutine SYMMAP2( badsym, nsym, n1r,n2r, isymctr,matsym,
     $ isymv,nsymv )
c---------------------------------------------------------------
c Purpose: find mesh points related by 2D symmetry operations.
c
c Written: Peter A. Schultz, 20-February-1990
c
c Revision history:
c   6Jan98-PAS/2.30: redo symmetry scheme
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Sep98-PAS/2.27: treat non-origin center of symmetry
c  22Jul94-PAS/catch incomplete-sym bug
c---------------------------------------------------------------
c
c  The "nsym" symmetry operators, represented as mesh index
c  transformation matrices, are given in matsym.  Note, however,
c  that the (m',n') = S(m,n) obtained by performing the symmetry
c  operation S on point (m,n) will, in general, lie outside the
c  unit cell.  To get it back in the unit cell, one must add
c  a lattice vector.  I do this by using the MOD operation on
c  (m',n') to get the (m",n") that lies within the unit cell.
c  Note that for this transformation, the mesh points are indexed
c  starting with 0 rather than 1 [so that the origin lies at (0,0)
c  rather than (1,1)].
c
c  Variables:
c    nsym      - number of symmetry operations
c    isymctr() - grid index of center of symmetry
c    matsym()  - array containing symmetry operations
c    isymv()   - index of first of symmetry-related mesh points
c    nsymv()   - number of symmetry related mesh points in that set
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input array:
      DIMENSION  matsym(3,4,nsym), isymctr(3)
c output arrays:
      DIMENSION  isymv(*),nsymv(*)
      LOGICAL    badsym
c
c >>>> EXECUTABLE CODE:
c
      badsym = .false.
c
      do  ir=1,n1r*n2r
        isymv(ir) = 0
        nsymv(ir) = 0
      enddo
c
c  Determine the symmetry-related sets of points:
c
      ir = 0
      do 102 i2s=0,n2r-1
        i2 = i2s - isymctr(2)
        do 101 i1s=0,n1r-1
          i1 = i1s - isymctr(1)
          ir = ir + 1
c
          if( isymv(ir).eq.0 )then
c           we have ourselves a virgin point
            isymv(ir) = ir
c           generate all of the symmetry-related points:
            do 50 isym=1,nsym
c             generate symmetry point ...
              j1 = matsym(1,1,isym)*i1 + matsym(2,1,isym)*i2
     $           + isymctr(1)
              j2 = matsym(1,2,isym)*i1 + matsym(2,2,isym)*i2
     $           + isymctr(2)
c              ... bring into primary unit cell:
              j1 = MOD( j1, n1r )
              j2 = MOD( j2, n2r )
c              ... 'cause MOD gives negative numbers:
              if( j1.lt.0 ) j1 = j1 + n1r
              if( j2.lt.0 ) j2 = j2 + n2r
              jr = n1r*j2 + j1 + 1
              if( isymv(jr).ne.0 .and. isymv(jr).ne.ir )then
                badsym = .true.
                RETURN
              endif
              isymv(jr) = ir
   50       continue
          endif
c
  101   continue
  102 continue
c
c  Count up the number of points in each symmetry set:
c
      do  ir=1,n1r*n2r
        ifst = isymv(ir)
        nsymv(ifst) = nsymv(ifst) + 1
      enddo
c
c    That's all Folks!
c
      RETURN
      END
