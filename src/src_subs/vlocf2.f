c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLOCF2
c
c
      subroutine VLOCF2(
     $ fesi0,fesi, fesj0,fesj, fxci0,fxci, fxcj0,fxcj,
     $ li,lj, ci,cj, gradci,gradcj,
     $ nx, vloces,vlocxc, twoali,twoalj, p,
     $ pr0, dpri0,dpri,pridp, dprj0,dprj,prjdp,
     $ psq1,psq2,psq3,  b1,b12,b123,bfac )
c---------------------------------------------------------------
c Purpose: compute quintuple sum of symmetry coefficients times
c          the radial integrals "vloces" and "vlocxc" to
c          computed derivatives at the orbital sites.
c
c Written: Peter A. Schultz, 30-September-1998 (for v2.27)
c          based on PJF's original "engausb" routine
c
c Revision history:
c   9Nov99-PAS/2.42: bugfix/cleanup
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c Notes:
c  This routine (along with its xc-less little sister vlocf1)
c  is the most expensive part of the local three-center force
c  calculation. In a couple of tests, it amounts to up to 2/3
c  of the cost of the local force integrals. Small changes in
c  this routine can have big consequences in efficiency.  The
c  trick is that a lot of loops are *very* short, and figuring
c  out how to shift work around to reduce op count and branching
c  is tricky.  I'm going to settle for what I've already done,
c  but extra work tuning the vlocfX routines would be worthwhile.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output arrays:
      DIMENSION  fesi(3),fesj(3), fxci(3),fxcj(3)
c input arrays:
      DIMENSION  ci(10),cj(10), gradci(3,10),gradcj(3,10)
c
      DIMENSION  vloces(6,nx),vlocxc(6,nx)
      DIMENSION  twoali(nx),twoalj(nx), p(3,nx)
c scratch arrays:
      DIMENSION  pr0(nx), dpri0(nx),dpri(3,nx),pridp(3,nx)
      DIMENSION           dprj0(nx),dprj(3,nx),prjdp(3,nx)
      DIMENSION  psq1(nx),psq2(nx),psq3(nx)
      DIMENSION  b1(nx),b12(nx),b123(nx),bfac(nx)
c local declarations:
      EQUIVALENCE  (nq(1),nq1),(nqmn(1),nqmn1)
      COMMON  /nqs/   nq1,nq2,nq3
      COMMON  /nqmns/ nqmn1,nqmn2,nqmn3
      DIMENSION  nq(3),nqmn(3)
      LOGICAL    EVEN
c
      DIMENSION  pow(3),powzer(3)
      DIMENSION  ns(9),np(10,3),bpq(5,5),nqmin(5)
      DATA  ns / 1, 2,3,4, 5,6,7,9,10 /
      DATA  np / 0, 1,0,0, 1,0,1,2,0,0,
     $           0, 0,1,0, 1,1,0,0,2,0,
     $           0, 0,0,1, 0,1,1,0,0,2 /
      DATA  bpq / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     $            0.d0, 1.d0, 1.d0, 0.d0, 0.d0,
     $            0.d0, 0.d0, 1.d0, 3.d0, 3.d0,
     $            0.d0, 0.d0, 0.d0, 1.d0, 6.d0,
     $            0.d0, 0.d0, 0.d0, 0.d0, 1.d0 /
      DATA  nqmin / 1,2,2,3,3 /
c
      DATA  zero,one,two / 0.d0,1.d0,2.d0 /
c statement function:
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
      nsimx = ns(li)
      nsjmx = ns(lj)
c
      do  ix=1,nx
        psq1(ix) = p(1,ix)**2
        psq2(ix) = p(2,ix)**2
        psq3(ix) = p(3,ix)**2
      enddo
c
      do 102 nsi=1,nsimx
       do 101 nsj=1,nsjmx
c
        cicj = ci(nsi)*cj(nsj)
        do  ix=1,nx
          pr0(ix) = cicj
          dpri(1,ix) = gradci(1,nsi)*cj(nsj)
          dpri(2,ix) = gradci(2,nsi)*cj(nsj)
          dpri(3,ix) = gradci(3,nsi)*cj(nsj)
          dprj(1,ix) = ci(nsi)*gradcj(1,nsj)
          dprj(2,ix) = ci(nsi)*gradcj(2,nsj)
          dprj(3,ix) = ci(nsi)*gradcj(3,nsj)
          dpri0(ix)   = -twoali(ix)*cicj
          dprj0(ix)   = -twoalj(ix)*cicj
          pridp(1,ix) = -twoali(ix)*cicj
          pridp(2,ix) = -twoali(ix)*cicj
          pridp(3,ix) = -twoali(ix)*cicj
          prjdp(1,ix) = -twoalj(ix)*cicj
          prjdp(2,ix) = -twoalj(ix)*cicj
          prjdp(3,ix) = -twoalj(ix)*cicj
        enddo
c
c  "pridp" is used for the case that pow(n)=1. and p(n)=0.  In this
c  case we have to be careful in taking derivative wrt p(n), since
c  have function equal zero whose derivative isn't.  In original code,
c  the zero/zero that arises in this case was evaluated incorrectly.
c
        do 20 n=1,3
          powzer(n) = zero
          nq(n) = np(nsi,n) + np(nsj,n) + 1
          nqmn(n) = nqmin(nq(n))
          if( EVEN(nq(n)) )then
            do  ix=1,nx
              pr0(ix) = pr0(ix)*p(n,ix)
              dpri0(ix) = dpri0(ix)*p(n,ix)
              dprj0(ix) = dprj0(ix)*p(n,ix)
            enddo
            do 18 nn=1,3
              do  ix=1,nx
                dpri(nn,ix) = dpri(nn,ix)*p(n,ix)
                dprj(nn,ix) = dprj(nn,ix)*p(n,ix)
              enddo
              if( nn.ne.n )then
                do  ix=1,nx
                  pridp(nn,ix) = pridp(nn,ix)*p(n,ix)
                  prjdp(nn,ix) = prjdp(nn,ix)*p(n,ix)
                enddo
              endif
   18       continue
            powzer(n) = powzer(n) + one
          else
            if( nq(n).gt.2 )then
              do  ix=1,nx
                if( p(n,ix).ne.zero )then
                  pridp(n,ix) = pridp(n,ix)/p(n,ix)
                  prjdp(n,ix) = prjdp(n,ix)/p(n,ix)
                endif
              enddo
            endif
          endif
   20   continue
c
        pow(1) = powzer(1)
        pow(2) = powzer(2)
        pow(3) = powzer(3)
c
        do  ix=1,nx
          b1(ix) = one
        enddo
c
        do 50 iq1=nqmn1,nq1
          do  ix=1,nx
            b12(ix) = b1(ix)*bpq(nq1,iq1)
          enddo
c
          do 40 iq2=nqmn2,nq2
            do  ix=1,nx
              b123(ix) = b12(ix)*bpq(nq2,iq2)
            enddo
c
            do 30 iq3=nqmn3,nq3
              do  ix=1,nx
                bfac(ix) = b123(ix)*bpq(nq3,iq3)
              enddo
              iqsum = iq1 + iq2 + iq3 - 2
c
              do  ix=1,nx
                ves = bfac(ix)*vloces(iqsum,ix)
                vxc = bfac(ix)*vlocxc(iqsum,ix)
                fesi(1) = fesi(1) - ves*dpri(1,ix)
                fesi(2) = fesi(2) - ves*dpri(2,ix)
                fesi(3) = fesi(3) - ves*dpri(3,ix)
                fxci(1) = fxci(1) - vxc*dpri(1,ix)
                fxci(2) = fxci(2) - vxc*dpri(2,ix)
                fxci(3) = fxci(3) - vxc*dpri(3,ix)
                fesj(1) = fesj(1) - ves*dprj(1,ix)
                fesj(2) = fesj(2) - ves*dprj(2,ix)
                fesj(3) = fesj(3) - ves*dprj(3,ix)
                fxcj(1) = fxcj(1) - vxc*dprj(1,ix)
                fxcj(2) = fxcj(2) - vxc*dprj(2,ix)
                fxcj(3) = fxcj(3) - vxc*dprj(3,ix)
              enddo
c
              if( cicj.ne.zero )then
                do  ix=1,nx
                  if( bfac(ix).ne.zero )then
                    ves = bfac(ix)*vloces(iqsum,ix)
                    vxc = bfac(ix)*vlocxc(iqsum,ix)
                    tfaci = bfac(ix)*dpri0(ix)
                    dvesi = vloces(iqsum+1,ix)*tfaci
                    dvxci = vlocxc(iqsum+1,ix)*tfaci
                    fesi0   = fesi0   - ves*dpri0(ix)
                    fesi(1) = fesi(1) + dvesi*p(1,ix)
                    fesi(2) = fesi(2) + dvesi*p(2,ix)
                    fesi(3) = fesi(3) + dvesi*p(3,ix)
                    fxci0   = fxci0   - vxc*dpri0(ix)
                    fxci(1) = fxci(1) + dvxci*p(1,ix)
                    fxci(2) = fxci(2) + dvxci*p(2,ix)
                    fxci(3) = fxci(3) + dvxci*p(3,ix)
                    tfacj = bfac(ix)*dprj0(ix)
                    dvesj = vloces(iqsum+1,ix)*tfacj
                    dvxcj = vlocxc(iqsum+1,ix)*tfacj
                    fesj0   = fesj0   - ves*dprj0(ix)
                    fesj(1) = fesj(1) + dvesj*p(1,ix)
                    fesj(2) = fesj(2) + dvesj*p(2,ix)
                    fesj(3) = fesj(3) + dvesj*p(3,ix)
                    fxcj0   = fxcj0   - vxc*dprj0(ix)
                    fxcj(1) = fxcj(1) + dvxcj*p(1,ix)
                    fxcj(2) = fxcj(2) + dvxcj*p(2,ix)
                    fxcj(3) = fxcj(3) + dvxcj*p(3,ix)
                    if( pow(1).ne.zero )then
                      pfaci = pow(1)*pridp(1,ix)
                      fesi(1) = fesi(1) + ves*pfaci
                      fxci(1) = fxci(1) + vxc*pfaci
                      pfacj = pow(1)*prjdp(1,ix)
                      fesj(1) = fesj(1) + ves*pfacj
                      fxcj(1) = fxcj(1) + vxc*pfacj
                    endif
                    if( pow(2).ne.zero )then
                      pfaci = pow(2)*pridp(2,ix)
                      fesi(2) = fesi(2) + ves*pfaci
                      fxci(2) = fxci(2) + vxc*pfaci
                      pfacj = pow(2)*prjdp(2,ix)
                      fesj(2) = fesj(2) + ves*pfacj
                      fxcj(2) = fxcj(2) + vxc*pfacj
                    endif
                    if( pow(3).ne.zero )then
                      pfaci = pow(3)*pridp(3,ix)
                      fesi(3) = fesi(3) + ves*pfaci
                      fxci(3) = fxci(3) + vxc*pfaci
                      pfacj = pow(3)*prjdp(3,ix)
                      fesj(3) = fesj(3) + ves*pfacj
                      fxcj(3) = fxcj(3) + vxc*pfacj
                    endif
                  endif
                enddo
              endif
c
              pow(3) = pow(3) + two
              do  ix=1,nx
                b123(ix) = b123(ix)*psq3(ix)
              enddo
   30       continue
c
            pow(3) = powzer(3)
            pow(2) = pow(2) + two
            do  ix=1,nx
              b12(ix) = b12(ix)*psq2(ix)
            enddo
   40     continue
c
          pow(2) = powzer(2)
          pow(1) = pow(1) + two
          do  ix=1,nx
            b1(ix) = b1(ix)*psq1(ix)
          enddo
   50   continue
c
  101  continue
  102 continue
c
c    That's all Folks!
c
      RETURN
      END
