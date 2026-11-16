c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMBUILD
c
c
      subroutine SYMBUILD( IWR, nsym, nsymi,nsymd,  ndimsy,
     $ isymtyp,syminfo,  rprim, rmatsym )
c---------------------------------------------------------------
c Purpose: build symmetry stuff
c
c Written: Peter A. Schultz, 30-May-2001, for 2.47
c
c Revision history:
c   1May02-PAS/2.53: reconcile sym and cell-opt
c  21Jun01-PAS/2.48: replace STOPs
c  30May01-PAS/2.47: extracted from main routine
c---------------------------------------------------------------
c
c  nsymi = number of input symmetry elements
c  nsym  = number of (non-identity) symmetry elements in group
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  isymtyp(*), syminfo(*)
      DIMENSION  rprim(3,3)
c output array:
      DIMENSION  rmatsym(*)
c local declaration:
      DATA       zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      if( nsymi.gt.0 )then
c        ... cast symmetry operations in cartesian matrix form:
        call SYMRMAT( nsymi,isymtyp,syminfo,
     $   nsymd, nsymj, rmatsym, rprim )
c        ... generate group from these elements:
        call SYMGRP( nsymd, nsymj,nsym, rmatsym, rprim )
c        ... put translations back into lattice units:
        call SYMLAT( nsym, rmatsym, rprim )
c
        write(IWR,*) 'Number of symmetry generators       =',nsymi
        write(IWR,*) 'Expands to number of elements(w/id) =',nsymj+1
        write(IWR,*) 'Size of output symmetry group (w/id)=',nsym+1
        if( ndimsy.eq.2 )then
c          ... check for requirement that slab be x-y + z-axis
          if( rprim(1,3).ne.zero .or. rprim(2,3).ne.zero .or.
     $        rprim(3,1).ne.zero .or. rprim(3,2).ne.zero )then
            write(IWR,*) '>>>>> slab normal must be z-axis for 2D'
            call STOPXERR( '2D sym workable iff z-axis is normal' )
          endif
        endif
      else
        nsym = 0
      endif
c
c    That's all Folks!
c
      RETURN
      END
