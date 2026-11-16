c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMRHO
c
c
      subroutine SYMRHO( IWR, subsym,ndimsy, nsym,nsymg,
     $ rsym, rmatsym, isymgrd,  rprim, n1r,n2r,n3r, hh,
     $ rho, iwrk )
c---------------------------------------------------------------
c Purpose: driving routine to symmetrize density on the grid
c
c Written: Peter A. Schultz, 3-November-1994
c
c Revision history:
c   1May02-PAS/2.53: reconcile sym and cell-opt
c  21Jun01-PAS/2.48: replace STOPs
c   6Jan98-PAS/2.30: redo symmetry scheme
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Sep98-PAS/2.27: treat non-origin center of symmetry
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c symmetry information:
      LOGICAL    subsym
      DIMENSION  rsym(3), rmatsym(12,*), isymgrd(*)
      DIMENSION  rprim(3,*)
c grid-lattice primitive vectors:
      DIMENSION  hh(3,3)
c density to be symmetrized and work space:
      DIMENSION  rho(*), iwrk(*)
c local:
      DIMENSION  isymctr(3)
      LOGICAL  badsym
c
c >>>> EXECUTABLE CODE:
c
c skip if nothing to symmetrize:
      if( nsym.le.0 ) RETURN
c
c Allocate storage for symmetry operations:
      isymops = 1
      is1 = isymops + 12*nsym
c
c  Generate index symmetry operations from cartesian operations:
c
      call SYMMAT( subsym, ibadsym, nsym,nsymg,
     $ rprim, n1r,n2r,n3r, hh, rsym, rmatsym, isymgrd,
     $ isymctr,iwrk(isymops) )
c
      if( nsym.ne.nsymg )then
        write(IWR,*) '***** NOTE: symmetry of grid is reduced'
        write(IWR,*) '  # of system symmetry operations(w/id)=',nsym+1
        write(IWR,*) '  # of viable symmetries on grid(w/id)=',nsymg+1
      endif
c
      if( ibadsym.gt.0 )then
        write(IWR,*) '>>>>> grid does not transform by sym#',ibadsym
        write(IWR,*) 'Cartesian representation of symmetry operation:'
        write(IWR,'(10x,3f15.10)') (rmatsym(ij,ibadsym),ij=1,12)
        call STOPXERR( 'symrho-r/grid does not xform by symmetry' )
      elseif( ibadsym.lt.0 )then
        call STOPXERR( 'symctrgr/sym center not on grid point' )
      endif
c
      if( ndimsy.ne.2 )then
c
c >>>> 3D SYMMETRIZATION (general case)
c
        is2 = is1 + n1r*n2r*n3r
c
        call SYMMAP3( badsym, nsymg,n1r,n2r,n3r,isymctr,iwrk(isymops),
     $   iwrk(is1),iwrk(is2) )
c
        if( badsym )then
          call STOPXERR( 'symset3d/ bad 3D symmetry group' )
        endif
c
        call SYMRHO3( n1r,n2r,n3r, iwrk(is1),iwrk(is2), rho )
c
      else
c
c >>>> 2D SYMMETRIZATION (special case)
c
        is2 = is1 + n1r*n2r
c
        call SYMMAP2( badsym, nsymg, n1r,n2r, isymctr,iwrk(isymops),
     $   iwrk(is1),iwrk(is2) )
c
        if( badsym )then
          call STOPXERR( 'symset2d/ bad 2D symmetry group' )
        endif
c
        call SYMRHO2( n1r,n2r,n3r, iwrk(is1),iwrk(is2), rho )
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
