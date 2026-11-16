c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KCART
c
c
      subroutine KCART( ndim, rprim, nk,veckg,veck )
c---------------------------------------------------------------
c Purpose: convert k-vectors to Cartesian(KCART)/lattice(KLATT)
c
c Written: Peter A. Schultz, 6-August-2003, for v2.55a
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
c Input: real space primitive vectors
      DIMENSION  rprim(3,3)
c
c Input<->Output: k-vectors: Cartesian and lattice units
      DIMENSION  veck(3,*),veckg(3,*)
c
c Local: k-point cell primitive vectors
      DIMENSION  gvec(3,3)
c
c >>>> EXECUTABLE CODE:
c
c If molecular, or just the gamma point, load zeroes:
      if( ndim.eq.0 .or.
     $    nk.eq.1 .and. veckg(1,1).eq.zero .and. veckg(2,1).eq.zero
     $            .and. veckg(3,1).eq.zero )then
        veck(1,1) = zero
        veck(2,1) = zero
        veck(3,1) = zero
        RETURN
      endif
c
c Generate reciprocal space primitive vectors:
      call KCELL( ndim, rprim, gvec )
c
c Expand Cartesian k-vectors from lattice-unit k-vectors:
      do  ik=1,nk
        call LVXPND( gvec, veck(1,ik), veckg(1,ik) )
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KLATT
c
c
      entry KLATT( ndim, rprim, nk,veck,veckg )
c
c If molecular, or just the gamma point, load zeroes:
      if( ndim.eq.0 .or.
     $    nk.eq.1 .and. veck(1,1).eq.zero .and. veck(2,1).eq.zero
     $            .and. veck(3,1).eq.zero )then
        veckg(1,1) = zero
        veckg(2,1) = zero
        veckg(3,1) = zero
        RETURN
      endif
c
c Generate reciprocal space primitive vectors:
      call KCELL( ndim, rprim, gvec )
c
c Project Cartesian k-vectors onto primitive k-vectors lattice:
      do  ik=1,nk
        call LVPROJ( gvec, veck(1,ik), veckg(1,ik) )
      enddo
c
      RETURN
c
c    That's all Folks!
c
      END
