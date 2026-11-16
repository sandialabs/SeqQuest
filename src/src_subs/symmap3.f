c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMMAP3
c
c
      subroutine SYMMAP3( badsym, nsym, n1r,n2r,n3r, isymctr,matsym,
     $ isymv,nsymv )
c---------------------------------------------------------------
c Purpose: generate image map of 3D symmetries
c
c Written: Peter A. Schultz, 3-November-1994
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Sep98-PAS/2.27: treat non-origin center of symmetry
c---------------------------------------------------------------
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
      do  ir=1,n1r*n2r*n3r
        isymv(ir) = 0
        nsymv(ir) = 0
      enddo
c
c  Determine the symmetry-related sets of points:
c
      ir = 0
      do 103 i3s=0,n3r-1
        i3 = i3s - isymctr(3)
        do 102 i2s=0,n2r-1
          i2 = i2s - isymctr(2)
          do 101 i1s=0,n1r-1
            i1 = i1s - isymctr(1)
            ir = ir + 1
c
            if( isymv(ir).eq.0 )then
c             we have ourselves a virgin point
              isymv(ir) = ir
c             generate all of the symmetry-related points:
              do 50 isym=1,nsym
c               generate symmetry point ...
                j1 = matsym(1,1,isym)*i1 + matsym(2,1,isym)*i2
     $             + matsym(3,1,isym)*i3 + matsym(1,4,isym)
     $             + isymctr(1)
                j2 = matsym(1,2,isym)*i1 + matsym(2,2,isym)*i2
     $             + matsym(3,2,isym)*i3 + matsym(2,4,isym)
     $             + isymctr(2)
                j3 = matsym(1,3,isym)*i1 + matsym(2,3,isym)*i2
     $             + matsym(3,3,isym)*i3 + matsym(3,4,isym)
     $             + isymctr(3)
c                ... bring into primary unit cell:
                j1 = MOD( j1, n1r )
                j2 = MOD( j2, n2r )
                j3 = MOD( j3, n3r )
c                ... 'cause MOD gives negative numbers:
                if( j1.lt.0 ) j1 = j1 + n1r
                if( j2.lt.0 ) j2 = j2 + n2r
                if( j3.lt.0 ) j3 = j3 + n3r
                jr = n1r*( n2r*j3 + j2 ) + j1 + 1
                if( isymv(jr).ne.0 .and. isymv(jr).ne.ir )then
                  badsym = .true.
                  RETURN
                endif
                isymv(jr) = ir
   50         continue
            endif
c
  101     continue
  102   continue
  103 continue
c
c  Count up the number of points in each symmetry set:
c
      do  ir=1,n1r*n2r*n3r
        ifst = isymv(ir)
        nsymv(ifst) = nsymv(ifst) + 1
      enddo
c
c    That's all Folks!
c
      RETURN
      END
