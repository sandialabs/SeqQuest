c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BSIALL
c
c
      subroutine BSIALL( nrs, llmxi, barg, bess, radmsh,
     $  x,w1,w2,w3,w4, xx )
c---------------------------------------------------------------
c Purpose: to compute half odd order bessel functions
c          bess(nr,k)=i(k-0.5,x(nr)), k=1,8 for x>1.d-40
c
c Written: Peter A. Schultz, 26-April-1989,  based on PJF's
c          adaptation of Don Amos' original bessel routines
c
c Revision history:
c  11Mar02-PAS/2.52: last literal reals to d0
c   8Aug00-PAS/2.45: renamed "chbsv" to "bsicbsv"
c  22Oct99-PAS/2.40: simple optimizations with chbsv
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Aug98-PAS/2.25: promote constants
c  22May98-PAS/2.22: cleanup
c---------------------------------------------------------------
c
c Arrays:
c    bess()    - output bessel functions on mesh
c    radmsh()  - input points on radial mesh
c    w1()-w4() - scratch arrays for chbsv routine
c    x()       - bessel argument on mesh=barg*radmsh()
c    xx()      - input parameters for chbsv routine
c
c Caution: w1(),w2(); w3,w4 need to be consecutive in mem, w1() and
c          w3() are used double sized in chbsv and bsixgt4 routines
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
      DIMENSION  bess(nrs,*),radmsh(nrs)
c scratch arrays:
      DIMENSION  x(nrs),xx(nrs)
      DIMENSION  w1(*),w2(nrs),w3(*),w4(nrs)
c
c local declarations:
      DATA  zero,two / 0.d0,2.d0 /
      DATA  four,five,fivp5 / 4.d0,5.d0,5.5d0 /
c
c >>>> EXECUTABLE CODE:
c
      rthafpi = SQRT( pi/two )
      nra = 0
      do  nr=1,nrs
        x(nr) = barg*radmsh(nr)
      enddo
c
c    Determine which range x's belong in:
c
      if( x(nrs).le.four )then
        nra = nrs
        goto 100
      elseif( x(1).gt.four )then
        goto 200
      endif
      nra = 1
      do  nr=2,nrs-1
        if( x(nr).gt.four ) goto 100
        nra = nra + 1
      enddo
c
c  * * * * *  x le 4  * * * * *
c
  100 continue
      do  nr=1,nra
        xx(nr) = x(nr) - two
      enddo
c
      call BSICBSV( nrs,nra, bess(1,6), xx, w1,w3 )
c
c pulling the next two loops separately out of the third
c loop was just a touch faster than combined loop(s)
c The cost here appears to be taken up in the EXP and LOG.
      do  nr=1,nra
        w2(nr) = five* LOG( x(nr) ) - x(nr)
      enddo
      do  nr=1,nra
        w3(nr) = rthafpi*EXP( w2(nr) )
      enddo
      do  nr=1,nra
        bess(nr,6) = bess(nr,6)*w3(nr)
        bess(nr,7) = bess(nr,7)*w3(nr)*x(nr)
        w1(nr) = two/x(nr)
        w2(nr) = fivp5*w1(nr)
      enddo
c get the rest by recursion:
c (note: unrolling the k-loop yields no performance gain).
      do  k=5,1,-1
        do  nr=1,nra
          bess(nr,k) = w2(nr)*bess(nr,k+1) + bess(nr,k+2)
          w2(nr) = w2(nr) - w1(nr)
        enddo
      enddo
c fwd recursion to get extra bessel fcn (if needed)
      if( llmxi.eq.8 )then
        do  nr=1,nra
          bess(nr,8) = bess(nr,6) - 13.d0*bess(nr,7)/x(nr)
        enddo
      endif
c
      if( nra.eq.nrs ) RETURN
c
c  * * * * *  x > 4  * * * * *
c
  200 continue
c
c Caution: w1() used double sized inside bsixgt4
c     This routine dependent on w1,w2 being consecutive
c     in memory (or at least that w1 can be safely double-sized).
c Note: array x() is destroyed in bsixgt4
c
      call BSIXGT4( nra+1,nrs, llmxi, bess, x,xx,w1,w3,w4 )
c
c    That's all Folks!
c
      RETURN
      END
