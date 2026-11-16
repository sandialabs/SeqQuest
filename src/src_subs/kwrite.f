c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KWRITE
c
c
      subroutine KWRITE( IWR, ndim, rprim, nk,nk0,wtk,wtkscal,veckg )
c---------------------------------------------------------------
c Purpose: write out k-vectors to output file
c
c Written: Peter A. Schultz, 6-August-2003, for v2.55a
c
c Revision history:
c   8Aug03-PAS/2.56: use LVPROJ instead of PROJLV
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input: real space primitive vectors; k-point weights, vectors
      DIMENSION  rprim(3,3)
      DIMENSION  wtk(*), veckg(3,*)
c Local:
      DIMENSION  gvec(3,3), vkr(3),vkg(3)
      DATA  zero / 0.d0 /
      DATA  lstout / 2 /
c
c >>>> EXECUTABLE CODE:
c
c If molecular, or just the gamma point, skip this:
      if( ndim.eq.0 ) RETURN
      if( nk.eq.1 .and. veckg(1,1).eq.zero .and. veckg(2,1).eq.zero
     $            .and. veckg(3,1).eq.zero ) RETURN
c
      if( IWR.le.0 ) RETURN
c
      if( lstout.gt.0 )then
c       Get reciprocal space primitive vectors:
        call KCELL( ndim, rprim, gvec )
c
c       Summarize output k-point set:
        write(IWR,*) 'IBZ Bloch vectors, total=',nk,nk0
        write(IWR,9101)
 9101   format(1x,'  K# ',2x,4x,'wt',4x,2x,9x,'k-vector:(r);(g)')
 9102   format(1x,  i5,   2x,  f10.6,   2x, 3f10.6, 2x,3f8.4 )
        if( lstout.gt.1 )then
          do  k=1,nk
            vkg(1) = veckg(1,k)
            vkg(2) = veckg(2,k)
            vkg(3) = veckg(3,k)
            call LVXPND( gvec, vkr, vkg )
            write(IWR,9102)  k, wtkscal*wtk(k), vkr, vkg
          enddo
        endif
      endif
c
c    That's all Folks!
c
      RETURN
      END
