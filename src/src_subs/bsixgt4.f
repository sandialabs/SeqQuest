c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BSIXGT4
c
c
      subroutine BSIXGT4( nrsmin,nrs, ltop, bessi,x,
     $ xinv,sum1,xtok,expm2r )
c---------------------------------------------------------------
c Revision history:
c   7Mar02-PAS/2.52: explicitly promote constants to r** (d0)
c  16Nov99-PAS/2.42: cleanup and commentary
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c Notes: there is no easy optimization available here.  Tried
c        several ideas, but none were helpful. PAS/16Nov99
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
      DIMENSION  bessi(nrs,ltop),x(nrs)
      DIMENSION  xinv(nrs),sum1(nrs,2),xtok(nrs),expm2r(nrs)
c
c the following data array has form ( (j-1)+(i-1) )! / ( (j-1)!*(i-1)! )
      DIMENSION  a(8,8)
      DATA  (a(i,1),i=1,8) / 8*1.d0 /
      DATA  (a(i,2),i=2,8)
     $ /  2.d0,  6.d0,  12.d0,  20.d0,  30.d0,  42.d0,  56.d0/
      DATA  (a(i,3),i=3,8)
     $ /  12.d0,  60.d0,  180.d0,  420.d0,  840.d0, 1512.d0/
      DATA  (a(i,4),i=4,8)
     $ /  120.d0,  840.d0,  3360.d0,  10080.d0, 25200.d0/
      DATA  (a(i,5),i=5,8)
     $ /  1680.d0,  15120.d0, 75600.d0,  277200.d0/
      DATA  (a(i,6),i=6,8)
     $ /  30240.d0,  332640.d0, 1995840.d0/
      DATA  (a(i,7),i=7,8)
     $ /  665280.d0,  8648640.d0/
      DATA   a(8,8)
     $ / 17297280.d0/
c
      DATA  zero,half,one / 0.d0,0.5d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
c  cutexp=largest neg. argument before exp underflow
      argmx  = half*cutexp
      argtst = half*argmx
      nrsmax = nrs
      do  nr=nrs,nrsmin,-1
        if( x(nr).lt.argtst ) goto 6
        nrsmax = nr - 1
        xinv(nr) = half/x(nr)
        expm2r(nr) = zero
      enddo
    6 continue
c
      do  nr=nrsmin,nrsmax
        xinv(nr) = half/x(nr)
        arg = x(nr) + x(nr)
        expm2r(nr) = EXP( -arg )
      enddo
c
      do  nr=nrsmin,nrs
        bessi(nr,1) = xinv(nr)*a(1,1)*(one-expm2r(nr))
      enddo
      if( ltop.ge.2 )then
        do  nr=nrsmin,nrs
          bessi(nr,2) = xinv(nr)*( a(2,1)*(one+expm2r(nr))
     $                  - xinv(nr)*a(2,2)*(one-expm2r(nr)) )
        enddo
      endif
c
      do 100 l=3,ltop
c
        do  nr=nrsmin,nrs
          expm2r(nr) = -expm2r(nr)
          xtok(nr) = xinv(nr)*xinv(nr)
          sum1(nr,1) = a(l,1)*xinv(nr)
          sum1(nr,2) = a(l,2)*xtok(nr)
        enddo
c
c  You'd think pulling this loop apart would help, but
c  it does not seem to.
        do 40 k=3,l
          is = 2 - MOD(k,2)
          alk = a(l,k)
          do  nr=nrsmin,nrs
            xtok(nr) = xtok(nr)*xinv(nr)
            sum1(nr,is) = sum1(nr,is) + alk*xtok(nr)
          enddo
   40   continue
c
        do  nr=nrsmin,nrs
          bessi(nr,l) = (sum1(nr,1)-sum1(nr,2))
     $                + (sum1(nr,1)+sum1(nr,2))*expm2r(nr)
        enddo
c
  100 continue
c
c    That's all Folks!
c
      RETURN
      END
