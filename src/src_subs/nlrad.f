c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NLRAD
c
c
      subroutine NLRAD( ven, alsum,bargi,bargj,bargs,
     $ exparg, llmxi,llmxj,lmxp1,lstmx,lpowmx,
     $ nrs,nrd,nr0, radmsh,vpsrad,
     $ bessi,bessj,radfac,w1,w2,w3,w4,w5,w6 )
c---------------------------------------------------------------
c Purpose: radial integrals needed to evaluate matrix elements
c   of non-local short-range part of effective nuclear potential
c
c Revision history:
c  14Mar00-_AS/2.43: rename "radnloc" to "nlrad"
c   7Dec99-PAS/2.42: clean-up loops
c  10Dec98-PAS/2.29: promote source to explicit double precision
c    Aug94-PAS/2.11: cleanup, external underflow control
c    May89-PAS: install optimized bessel function routines
c---------------------------------------------------------------
c
c  Assume vpsrad() contains l-dependent potential times "weight"
c  necessary for numerical quadrature on the radial mesh "radmsh"
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output array:
      DIMENSION  ven(llmxi,llmxj,lstmx)
c input arrays:
      DIMENSION  radmsh(*),vpsrad(nrd,4)
c scratch arrays:
      DIMENSION  bessi(nrs,*),bessj(nrs,*),radfac(nrs,*),
     $   w1(nrs),w2(nrs),w3(nrs),w4(nrs),w5(nrs),w6(nrs)
c local stuff:
      DATA  zero,one,two / 0.d0,1.d0,2.d0 /
c
c >>>> EXECUTABLE CODE:
c
      call MKZERO( llmxi*llmxj*lstmx, ven )
c
      do  nr=1,nrs
        r = radmsh(nr0+nr)
        w2(nr) = r*EXP(-alsum*r**2 + bargs*r - exparg)
      enddo
c
c    Multiply potential by l-independent factors
c
      do  llpot=1,lmxp1
        do  nr=1,nrs
          radfac(nr,llpot) = vpsrad(nr0+nr,llpot)*w2(nr)
        enddo
      enddo
c
c    Compute spherical bessel functions of imaginary argument
c
      if( bargi.gt.zero )then
        call BSIALL( nrs, llmxi, bargi, bessi, radmsh(nr0+1),
     $   w1,w2,w3,w4,w5,w6 )
      else
        do  nr=1,nrs
          bessi(nr,1) = one
        enddo
      endif
c
      if( bargj.gt.zero )then
        call BSIALL( nrs, llmxj, bargj, bessj, radmsh(nr0+1),
     $   w1,w2,w3,w4,w5,w6)
      else
        do  nr=1,nrs
          bessj(nr,1) = one
        enddo
      endif
c
c    Loop on sum of l-primes and lpot
c
      lstore = 0
      do 102 llpot=1,lmxp1
        do 101 lpower=1,lpowmx
          lstore = lstore + 1
          do  nr=1,nrs
            radfac(nr,llpot) = radfac(nr,llpot)*radmsh(nr0+nr)
          enddo
c
          do 92 llj=1,llmxj
            do 91 lli=1,llmxi
              do  nr=1,nrs
                ven(lli,llj,lstore) = ven(lli,llj,lstore) +
     $           radfac(nr,llpot)*bessi(nr,lli)*bessj(nr,llj)
              enddo
   91       continue
   92     continue
c
  101   continue
  102 continue
c
c    That's all Folks!
c
      RETURN
      END
