c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLOCF1
c
c
      subroutine VLOCF1(
     $ fesi0,fesi, fesj0,fesj,
     $ li,lj, ci,cj, gradci,gradcj,
     $ nx, vloces, twoali,twoalj, p,
     $ tves,dpr,   xx1,xx2,   tdves,pridp,prjdp,
     $ psq1,psq2,psq3,  b1,b12,b123,bfac )
c---------------------------------------------------------------
c Purpose: compute quintuple sum of symmetry coefficients times
c          the radial integrals "vloces" to
c          computed derivatives at the orbital sites.
c
c Written: Peter A. Schultz, 30-September-1998 (for v2.27)
c          based on PJF's original "engausb" routine
c
c Revision history:
c  11Feb00-PAS/2.43: some simple optimizations
c   9Nov99-PAS/2.42: bugfix/cleanup
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c Notes:
c  This routine is just the vlocf2 routine with the xc stuff
c  ripped out of it.  Fully identical otherwise.
c  NB: vlocf1 tuned for performance, but same tuning NOT carried
c  over to vlocf0 and vlocf2. PAS/14Feb00
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output arrays:
      DIMENSION  fesi(3),fesj(3)
c input arrays:
      DIMENSION  ci(10),cj(10), gradci(3,10),gradcj(3,10)
c
      DIMENSION  vloces(6,nx)
      DIMENSION  twoali(nx),twoalj(nx), p(3,nx)
c scratch arrays:
      DIMENSION  tves(nx), dpr(nx)
      DIMENSION  xx1(3,nx),xx2(3,nx)
      DIMENSION  tdves(nx), pridp(3,nx),prjdp(3,nx)
      DIMENSION  psq1(nx),psq2(nx),psq3(nx)
      DIMENSION  b1(nx),b12(nx),b123(nx),bfac(nx)
c local declarations:
      EQUIVALENCE  (nq(1),nq1),(nqmn(1),nqmn1)
      COMMON  /nqs/   nq1,nq2,nq3
      COMMON  /nqmns/ nqmn1,nqmn2,nqmn3
      DIMENSION  nq(3),nqmn(3)
      LOGICAL    ci_eq_0
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
       cii = ci(nsi)
       ci_eq_0 = .false.
       if( cii.eq.zero ) ci_eq_0 = .true.
c
       do 101 nsj=1,nsjmx
        cjj = cj(nsj)
        if( ci_eq_0 .and. cjj.eq.zero ) goto 101
c
        cicj = cii*cjj
        dci1cj = gradci(1,nsi)*cjj
        dci2cj = gradci(2,nsi)*cjj
        dci3cj = gradci(3,nsi)*cjj
        cidcj1 = cii*gradcj(1,nsj)
        cidcj2 = cii*gradcj(2,nsj)
        cidcj3 = cii*gradcj(3,nsj)
c
        if( cicj.eq.zero )then
c
          do  ix=1,nx
            dpr(ix) = one
          enddo
c
          do 20 n=1,3
            nq(n) = np(nsi,n) + np(nsj,n) + 1
            nqmn(n) = nqmin(nq(n))
            if( EVEN( nq(n) ) )then
              do  ix=1,nx
                dpr(ix) = dpr(ix)*p(n,ix)
              enddo
            endif
   20     continue
c
          do  ix=1,nx
            b1(ix) = dpr(ix)
            tves(ix) = zero
          enddo
c
          do 80 iq1=nqmn1,nq1
            do  ix=1,nx
              b12(ix) = b1(ix)*bpq(nq1,iq1)
            enddo
c
            do 60 iq2=nqmn2,nq2
              do  ix=1,nx
                b123(ix) = b12(ix)*bpq(nq2,iq2)
              enddo
c
              do 40 iq3=nqmn3,nq3
                cb = bpq(nq3,iq3)
                iqsum = iq1 + iq2 + iq3 - 2
c
                do  ix=1,nx
                  tves(ix) = tves(ix) + cb*b123(ix)*vloces(iqsum,ix)
                enddo
c
                do  ix=1,nx
                  b123(ix) = b123(ix)*psq3(ix)
                enddo
   40         continue
c
              do  ix=1,nx
                b12(ix) = b12(ix)*psq2(ix)
              enddo
   60       continue
c
            do  ix=1,nx
              b1(ix) = b1(ix)*psq1(ix)
            enddo
   80     continue
c
          fint = zero
          do  ix=1,nx
            fint = fint + tves(ix)
          enddo
          fesi(1) = fesi(1) - fint*dci1cj
          fesi(2) = fesi(2) - fint*dci2cj
          fesi(3) = fesi(3) - fint*dci3cj
          fesj(1) = fesj(1) - fint*cidcj1
          fesj(2) = fesj(2) - fint*cidcj2
          fesj(3) = fesj(3) - fint*cidcj3
c
        else
c
          do  ix=1,nx
            dpr(ix) = one
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
          do 120 n=1,3
            powzer(n) = zero
            nq(n) = np(nsi,n) + np(nsj,n) + 1
            nqmn(n) = nqmin(nq(n))
            if( EVEN( nq(n) ) )then
              do  ix=1,nx
                dpr(ix) = dpr(ix)*p(n,ix)
              enddo
              do 118 nn=1,3
                if( nn.ne.n )then
                  do  ix=1,nx
                    pridp(nn,ix) = pridp(nn,ix)*p(n,ix)
                    prjdp(nn,ix) = prjdp(nn,ix)*p(n,ix)
                  enddo
                endif
  118         continue
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
  120     continue
c
          pow(1) = powzer(1)
          pow(2) = powzer(2)
          pow(3) = powzer(3)
c
          do  ix=1,nx
            b1(ix) = one
            tves(ix) = zero
            tdves(ix) = zero
          enddo
c
          do 180 iq1=nqmn1,nq1
            do  ix=1,nx
              b12(ix) = b1(ix)*bpq(nq1,iq1)
            enddo
c
            do 160 iq2=nqmn2,nq2
              do  ix=1,nx
                b123(ix) = b12(ix)*bpq(nq2,iq2)
              enddo
c
              do 140 iq3=nqmn3,nq3
                iqsum = iq1 + iq2 + iq3 - 2
c
                cb = bpq(nq3,iq3)
                do  ix=1,nx
                  bfac(ix) = b123(ix)*cb
                enddo
c
                do  ix=1,nx
                  tves(ix) = tves(ix) + bfac(ix)*vloces(iqsum,ix)
                  tdves(ix) = tdves(ix) + bfac(ix)*vloces(iqsum+1,ix)
                enddo
c
                do  ix=1,nx
                  ves = bfac(ix)*vloces(iqsum,ix)
                  fesi(1) = fesi(1) + ves*pow(1)*pridp(1,ix)
                  fesi(2) = fesi(2) + ves*pow(2)*pridp(2,ix)
                  fesi(3) = fesi(3) + ves*pow(3)*pridp(3,ix)
                  fesj(1) = fesj(1) + ves*pow(1)*prjdp(1,ix)
                  fesj(2) = fesj(2) + ves*pow(2)*prjdp(2,ix)
                  fesj(3) = fesj(3) + ves*pow(3)*prjdp(3,ix)
                enddo
c
                pow(3) = pow(3) + two
                do  ix=1,nx
                  b123(ix) = b123(ix)*psq3(ix)
                enddo
  140         continue
c
              pow(3) = powzer(3)
              pow(2) = pow(2) + two
              do  ix=1,nx
                b12(ix) = b12(ix)*psq2(ix)
              enddo
  160       continue
c
            pow(2) = powzer(2)
            pow(1) = pow(1) + two
            do  ix=1,nx
              b1(ix) = b1(ix)*psq1(ix)
            enddo
  180     continue
c
          do  ix=1,nx
            tves(ix) = tves(ix)*dpr(ix)
            tdves(ix) = tdves(ix)*dpr(ix)
          enddo
          fint   = zero
          finti0 = zero
          fintj0 = zero
          finti1 = zero
          finti2 = zero
          finti3 = zero
          fintj1 = zero
          fintj2 = zero
          fintj3 = zero
          do  ix=1,nx
            fint   = fint   - tves(ix)
            finti0 = finti0 + tves(ix)*twoali(ix)
            fintj0 = fintj0 + tves(ix)*twoalj(ix)
            finti1 = finti1 - tdves(ix)*twoali(ix)*p(1,ix)
            finti2 = finti2 - tdves(ix)*twoali(ix)*p(2,ix)
            finti3 = finti3 - tdves(ix)*twoali(ix)*p(3,ix)
            fintj1 = fintj1 - tdves(ix)*twoalj(ix)*p(1,ix)
            fintj2 = fintj2 - tdves(ix)*twoalj(ix)*p(2,ix)
            fintj3 = fintj3 - tdves(ix)*twoalj(ix)*p(3,ix)
          enddo
          fesi(1) = fesi(1) + fint*dci1cj + finti1*cicj
          fesi(2) = fesi(2) + fint*dci2cj + finti2*cicj
          fesi(3) = fesi(3) + fint*dci3cj + finti3*cicj
          fesj(1) = fesj(1) + fint*cidcj1 + fintj1*cicj
          fesj(2) = fesj(2) + fint*cidcj2 + fintj2*cicj
          fesj(3) = fesj(3) + fint*cidcj3 + fintj3*cicj
          fesi0   = fesi0   + finti0*cicj
          fesj0   = fesj0   + fintj0*cicj
c
        endif
c
  101  continue
  102 continue
c
c    That's all Folks!
c
      RETURN
      END
