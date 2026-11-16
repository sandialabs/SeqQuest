c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLOCF0
c
c
      subroutine VLOCF0( fxez,fyez,fzez,
     $ li,lj, ci,cj, gradci,gradcj,
     $ nx, vlocez, twoals, p,
     $ pr0, dpr0,pr0dp,
     $ psq1,psq2,psq3, b1,b12,b123,bfac )
c---------------------------------------------------------------
c Purpose: compute quintuple sum of symmetry coefficients times
c          the radial integrals "vlocez" to get es force integrals
c          for reference atom occupations
c
c Written: Peter A. Schultz, 30-September-1998 (for v2.27)
c          based on PJF's original "engausf" routine
c
c Revision history:
c   9Nov99-PAS/2.42: bugfix/cleanup
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  ci(10),cj(10), gradci(3,10),gradcj(3,10)
c
      DIMENSION  vlocez(6,nx)
      DIMENSION  twoals(nx), p(3,nx)
c scratch arrays:
      DIMENSION  pr0(nx),dpr0(3,nx),pr0dp(3,nx)
      DIMENSION  psq1(nx),psq2(nx),psq3(nx)
      DIMENSION  b1(49),b12(49),b123(49), bfac(49)
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
      fxez = zero
      fyez = zero
      fzez = zero
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
          pr0dp(1,ix) = -twoals(ix)*pr0(ix)
          pr0dp(2,ix) = -twoals(ix)*pr0(ix)
          pr0dp(3,ix) = -twoals(ix)*pr0(ix)
          dpr0(1,ix) = pr0(ix)*p(1,ix)
     $     - gradci(1,nsi)*cj(nsj) - ci(nsi)*gradcj(1,nsj)
          dpr0(2,ix) = pr0(ix)*p(2,ix)
     $     - gradci(2,nsi)*cj(nsj) - ci(nsi)*gradcj(2,nsj)
          dpr0(3,ix) = pr0(ix)*p(3,ix)
     $     - gradci(3,nsi)*cj(nsj) - ci(nsi)*gradcj(3,nsj)
        enddo
c
c  "pr0dp" is used for the case that pow(n)=1. and p(n)=0.  In this
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
            enddo
            do 18 nn=1,3
              do  ix=1,nx
                dpr0(nn,ix) = dpr0(nn,ix)*p(n,ix)
              enddo
              if( nn.ne.n )then
                do  ix=1,nx
                  pr0dp(nn,ix) = pr0dp(nn,ix)*p(n,ix)
                enddo
              endif
   18       continue
            powzer(n) = powzer(n) + one
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
              iqsum = iq1 + iq2 + iq3 - 2
              do  ix=1,nx
                bfac(ix) = b123(ix)*bpq(nq3,iq3)
              enddo
c
              do  ix=1,nx
                vezfac = bfac(ix)*vlocez(iqsum,ix)
                tfac = bfac(ix)*twoals(ix)*pr0(ix)
                dvezfac = vlocez(iqsum+1,ix)*tfac
                fxez = fxez + vezfac*dpr0(1,ix) - dvezfac*p(1,ix)
                fyez = fyez + vezfac*dpr0(2,ix) - dvezfac*p(2,ix)
                fzez = fzez + vezfac*dpr0(3,ix) - dvezfac*p(3,ix)
              enddo
              do  ix=1,nx
                vezfac = bfac(ix)*vlocez(iqsum,ix)
                if( pow(1).eq.one .and. p(1,ix).eq.zero )
     $            fxez = fxez + vezfac*pr0dp(1,ix)
                if( pow(2).eq.one .and. p(2,ix).eq.zero )
     $            fyez = fyez + vezfac*pr0dp(2,ix)
                if( pow(3).eq.one .and. p(3,ix).eq.zero )
     $            fzez = fzez + vezfac*pr0dp(3,ix)
                if( pr0(ix).ne.zero )then
                  vezfac = vezfac*twoals(ix)*pr0(ix)
                  if( p(1,ix).ne.zero .and. pow(1).ne.zero )
     $              fxez = fxez - vezfac*pow(1)/p(1,ix)
                  if( p(2,ix).ne.zero .and. pow(2).ne.zero )
     $              fyez = fyez - vezfac*pow(2)/p(2,ix)
                  if( p(3,ix).ne.zero .and. pow(3).ne.zero )
     $              fzez = fzez - vezfac*pow(3)/p(3,ix)
                endif
              enddo
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
