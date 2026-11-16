c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NLINT
c
c
      subroutine NLINT( vnlr,vnli, vnlint,
     $ limxp1,ljmxp1,lmxp1,lstmx,lpowmx, rimag,rjmag,
     $ li,lj, ai,aj, ylmj,angsav,
     $ nk,nlat, coskr,sinkr, jlat )
c---------------------------------------------------------------
c Purpose: evaluate lm-sums for non-local part of
c          effective potential matrix elements
c
c Revision history:
c   2Oct02-PAS/2.54: single jlat version, extracted from NLINTEQ
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output arrays:
      DIMENSION  vnlr(nk),vnli(nk)
c input arrays:
      DIMENSION  vnlint(limxp1,ljmxp1,lstmx)
      DIMENSION  ai(9),aj(9), ylmj(2,25)
      DIMENSION  angsav(limxp1,ljmxp1*ljmxp1,lmxp1,9,9)
      DIMENSION  coskr(nk,nlat),sinkr(nk,nlat)
c local declarations:
      DIMENSION  lang(9)
      DATA  lang / 0, 1,1,1, 2,2,2,2,2 /
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      call MKZERO( nk, vnlr )
      call MKZERO( nk, vnli )
c
      vnltmp = zero
      do 200 lmpi=1,li
        if( ai(lmpi).eq.zero ) goto 200
c
        lpi = lang(lmpi)
c
        do 100 lmpj=1,lj
          if( aj(lmpj).eq.zero ) goto 100
c
          lpj = lang(lmpj)
          aij = ai(lmpi)*aj(lmpj)
c
          do 90 llpot=1,lmxp1
            lstore = lpi + lpj + 1 + lpowmx*(llpot-1)
c
            do 70 lli=1,limxp1
              jlm = 0
              jlx = 0
c
              mmj = 1
                do 41 llj=mmj,ljmxp1
                  jlm = jlm + 1
c                   jpprl = 1
                    jlx = jlx + 1
                    vnltmp = vnltmp + aij*vnlint(lli,llj,lstore)*
     $               ylmj(1,jlm)*angsav(lli,jlx,llpot,lmpi,lmpj)
   41           continue
c
              do 52 mmj=2,ljmxp1
                do 51 llj=mmj,ljmxp1
                  jlm = jlm + 1
c                   jpprl = 1
                    jlx = jlx + 1
                    vnltmp = vnltmp + aij*vnlint(lli,llj,lstore)*
     $               ylmj(1,jlm)*angsav(lli,jlx,llpot,lmpi,lmpj)
c                   jpprl = 2
                    jlx = jlx + 1
                    vnltmp = vnltmp + aij*vnlint(lli,llj,lstore)*
     $               ylmj(2,jlm)*angsav(lli,jlx,llpot,lmpi,lmpj)
c
   51           continue
   52         continue
c
   70       continue
c
   90     continue
c
  100   continue
  200 continue
c
        vint = vnltmp
        do  k=1,nk
          vnlr(k) = vnlr(k) + vint*coskr(k,jlat)
          vnli(k) = vnli(k) + vint*sinkr(k,jlat)
        enddo
c
c    That's all Folks!
c
      RETURN
      END
