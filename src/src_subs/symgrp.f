c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMGRP
c
c
      subroutine SYMGRP( nsymd,nsym0,nsym, rmatsym, rprim )
c---------------------------------------------------------------
c Purpose: produce a set of symmetry elements that form a group
c          from a given collection of symmetry elements, all
c          represented in cartesian form.
c          (i.e., compute all products, eliminate duplicates)
c
c Written: Peter A. Schultz, 4-January-1999, for v2.30
c
c Revision history:
c  21Jun01-PAS/2.48: replace STOPs
c---------------------------------------------------------------
c
c  Variables:
c   nsymd - maximum number of allowed symmetry elements
c   nsym0 - # of non-identity (e) symmetry elements on input
c   nsym  - # of output (non-e) symmetry elements in group
c   rmatsym() - cartesian representations of symmetry elements
c   rprim(3,3) - primitive lattice vectors
c
c  Note:
c   The algorithm employed to generate and validate the group
c   is gruesomely wasteful, but it is not worth the effort to
c   cook up something more clever.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input/output array:
      DIMENSION  rmatsym(3,4,*)
c input array:
      DIMENSION  rprim(3,3)
c local declarations:
      DIMENSION  ematsym(3,4), tmatsym(3,4)
c  set symmetry identity:
      DATA  ematsym / 1.d0, 0.d0, 0.d0,
     $                0.d0, 1.d0, 0.d0,
     $                0.d0, 0.d0, 1.d0,
     $                0.d0, 0.d0, 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      nsym = 0
      if( nsym0.lt.1 ) goto 999
c
c  Clean out duplicates from input set of symmetry elements:
c
      do 100 isym=1,nsym0
c       Check that candidate is not identity ...
        call DUPEMAT( 12, ksym0, rmatsym(1,1,isym), 1,ematsym )
c          ... and not the same as earlier element:
        ksym1 = 0
        if( ksym0.eq.0 .and. nsym.gt.0 )
     $   call DUPEMAT( 12, ksym1, rmatsym(1,1,isym), nsym,rmatsym )
c
        if( ksym0.le.0 .and. ksym1.le.0 )then
          nsym = nsym + 1
          if( nsym.ne.isym )then
            do  j=1,4
              rmatsym(1,j,nsym) = rmatsym(1,j,isym)
              rmatsym(2,j,nsym) = rmatsym(2,j,isym)
              rmatsym(3,j,nsym) = rmatsym(3,j,isym)
            enddo
          endif
        endif
  100 continue
c
c  Now have nsym = # unique non-identity symmetry elements from input
c  set.  Want to expand this into a full group of symmetry elements.
c  I take all possible products (both sides), and check against the
c  identity and all existing elements, and keep the unique elements.
c  Can do much smarter than this pairwise products with full check,
c  but that would require tedious preprocessing and bookkeeping, and
c  I have forgotten too much of my group theory (sigh).  The current
c  scheme is simpler to code and, even though order nsym**3, is good
c  enough because the number of elements nsym will always be small.
c
      isym = 0
  200 isym = isym + 1
      if( isym.gt.nsym ) goto 999
c
      do 300 jsym=1,isym
c
        js = jsym
        is = isym
c      Compute product of symmetry elements:
  210   call SYMPROD( rmatsym(1,1,js),rmatsym(1,1,is), tmatsym, rprim )
c
c      Check that product is not identity ...
        call DUPEMAT( 12, ksym, tmatsym, 1,ematsym )
        if( ksym.gt.0 ) goto 300
c
c       ... and not the same as existing group element:
        call DUPEMAT( 12, ksym, tmatsym, nsym,rmatsym )
        if( ksym.gt.0 ) goto 290
c
c  We have a unique new group element, add it to set:
        if( nsym.gt.nsymd )then
          call STOPXERR( 'symgrp-d/too many sym elements in group' )
        endif
c
        nsym = nsym + 1
        do  j=1,4
          rmatsym(1,j,nsym) = tmatsym(1,j)
          rmatsym(2,j,nsym) = tmatsym(2,j)
          rmatsym(3,j,nsym) = tmatsym(3,j)
        enddo
c
  290   continue
        if( js.lt.isym )then
c         Check reverse product (could be non-commutative group)
          js = isym
          is = jsym
          goto 210
        endif
c
  300 continue
      goto 200
c
  999 continue
c
c    That's all Folks!
c
      RETURN
      END
