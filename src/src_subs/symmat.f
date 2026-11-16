c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMMAT
c
c
      subroutine SYMMAT( subsym, ibadsym, nsym,nsymg,
     $ rprim, n1r,n2r,n3r, hh, rsym, rmatsym,isymgrd,
     $ isymctr,matsym )
c---------------------------------------------------------------
c Purpose: generate matrix reps of 3D symmetry operations
c
c Written: Peter A. Schultz, 6-January-1999
c          Adapted from symmat3, and also replacing symmat2
c
c Revision history:
c   1May02-PAS/2.53: reconcile sym and cell-opt
c   5Jan99-PAS/2.30: redo symmetry scheme
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Sep98-PAS/2.27: treat non-origin center of symmetry
c---------------------------------------------------------------
c
c  Variables:
c   subsym    - true  = ok to reduce symmetry
c               false = trigger error upon symmetry failure
c   nsym      - # input non-identity symmetry operations
c   nsymg     - # non-identity symmetry operations used by grid
c   hh()      - the mesh lattice point basis vectors
c   rsym()    - cartesian coordinates of symmetry center
c   rmatsym() - the symmetry operations in cartesian form
c   isymctr() - output grid index of symmetry center
c   matsym()  - the output mesh index transformation matrix
c   isymgrd() - list of which of all symmetries ok on grid
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c input symmetry data:
      LOGICAL    subsym
      DIMENSION  rsym(3), rmatsym(3,4,*)
c grid-generating primitive vectors:
      DIMENSION  hh(3,3)
      DIMENSION  rprim(3,*)
c matrix reps of symmetry ops:
      DIMENSION  isymctr(3), matsym(3,4,*)
      DIMENSION  isymgrd(*)
c local declarations:
      DIMENSION  rsymctr(3)
      DIMENSION  a(3),b(3),c(3)
      DIMENSION  u(3), off(3),of2(3)
      DIMENSION  sa(3),sb(3),sc(3), sa2(3),sb2(3),sc2(3)
c
      DATA  zero,one,two, extra,dmax / 0.d0,1.d0,2.d0, 1.0001,1.d-7 /
c
c >>>> EXECUTABLE CODE:
c
      ibadsym = 0
c
      do  i=1,3
        a(i) = hh(i,1)
        b(i) = hh(i,2)
        c(i) = hh(i,3)
      enddo
c
c Copy symmetry center:
      do  i=1,3
        rsymctr(i) = rsym(i)
      enddo
c
      isymg = 1
      nsymg = nsym
      do 100 isym=1,nsym
c
c       Transform grid generating vectors by this symmetry operation:
        do  i=1,3
          sa(i) = rmatsym(1,i,isym)*a(1) + rmatsym(2,i,isym)*a(2)
     $          + rmatsym(3,i,isym)*a(3)
          sb(i) = rmatsym(1,i,isym)*b(1) + rmatsym(2,i,isym)*b(2)
     $          + rmatsym(3,i,isym)*b(3)
          sc(i) = rmatsym(1,i,isym)*c(1) + rmatsym(2,i,isym)*c(2)
     $          + rmatsym(3,i,isym)*c(3)
        enddo
        call LVXPND( rprim, off, rmatsym(1,4,isym) )
c
c       Project these vectors onto original primitive vectors:
c        (multiplied by "extra" to prevent truncation error)
c
        call PROJLV( sa, a,b,c, sa2 )
        matsym(1,1,isymg) = extra*sa2(1)
        matsym(1,2,isymg) = extra*sa2(2)
        matsym(1,3,isymg) = extra*sa2(3)
        call PROJLV( sb, a,b,c, sb2 )
        matsym(2,1,isymg) = extra*sb2(1)
        matsym(2,2,isymg) = extra*sb2(2)
        matsym(2,3,isymg) = extra*sb2(3)
        call PROJLV( sc, a,b,c, sc2 )
        matsym(3,1,isymg) = extra*sc2(1)
        matsym(3,2,isymg) = extra*sc2(2)
        matsym(3,3,isymg) = extra*sc2(3)
c        ... and load offset vector:
        call PROJLV( off, a,b,c, of2 )
        matsym(1,4,isymg) = extra*of2(1)
        matsym(2,4,isymg) = extra*of2(2)
        matsym(3,4,isymg) = extra*of2(3)
c
c       Check to be sure they have the right symmetries:
c
        do  j=1,3
          sa2(j) = a(j)*matsym(1,1,isymg) + b(j)*matsym(1,2,isymg)
     $           + c(j)*matsym(1,3,isymg)
          sb2(j) = a(j)*matsym(2,1,isymg) + b(j)*matsym(2,2,isymg)
     $           + c(j)*matsym(2,3,isymg)
          sc2(j) = a(j)*matsym(3,1,isymg) + b(j)*matsym(3,2,isymg)
     $           + c(j)*matsym(3,3,isymg)
          of2(j) = a(j)*matsym(1,4,isymg) + b(j)*matsym(2,4,isymg)
     $           + c(j)*matsym(3,4,isymg)
        enddo
c
        difa = SQRT( (sa(1)-sa2(1))**2 + (sa(2)-sa2(2))**2
     $             + (sa(3)-sa2(3))**2 )
        difb = SQRT( (sb(1)-sb2(1))**2 + (sb(2)-sb2(2))**2
     $             + (sb(3)-sb2(3))**2 )
        difc = SQRT( (sc(1)-sc2(1))**2 + (sc(2)-sc2(2))**2
     $             + (sc(3)-sc2(3))**2 )
        difo = SQRT( (off(1)-of2(1))**2 +
     $               (off(2)-of2(2))**2 +
     $               (off(3)-of2(3))**2 )
c
        if( difa .gt. dmax .or. difb .gt. dmax .or. difc .gt. dmax
     $     .or. difo .gt. dmax )then
          nsymg = nsymg - 1
          if( .not. subsym )then
c           We have a hard check, and this symmetry failed
            ibadsym = isym
            RETURN
          endif
        else
c         Have a good new symmetry element
          isymgrd(isymg) = isym
          isymg = isymg + 1
        endif
c
  100 continue
c
c And take care of symmetry center vector:
c
      call PROJLV( rsymctr, a,b,c, u )
      isymctr(1) = extra*u(1)
      isymctr(2) = extra*u(2)
      isymctr(3) = extra*u(3)
      do  j=1,3
        of2(j) = a(j)*isymctr(1) + b(j)*isymctr(2)
     $         + c(j)*isymctr(3)
      enddo
      difs = SQRT( (rsymctr(1)-of2(1))**2 +
     $             (rsymctr(2)-of2(2))**2 +
     $             (rsymctr(3)-of2(3))**2 )
      if( difs.gt.dmax ) ibadsym = -1
c
c    That's all Folks!
c
      RETURN
      END
