c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NLINTEQ
c
c
      subroutine NLINTEQ( vnlr,vnli, vnlint,
     $ limxp1,ljmxp1,lmxp1,lstmx,lpowmx, rimag,rjmag,
     $ li,lj,izero,ai,aj,ylmj,angsav,
     $ nk,nlat, coskr,sinkr,
     $ neqlat,jlateq,  vnltmp )
c---------------------------------------------------------------
c Purpose: evaluate lm-sums for non-local part of
c          effective potential matrix elements
c
c Revision history:
c   2Oct02-PAS/2.54: optimization via unrolls, variable retasks
c   7Mar00-PAS/2.43: rename "lloops" to "nlint"
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Apr97-PAS/2.17: naming conventions
c  29Sep95-PAS/2.13; extract k-loop from bottom of loops
c  12Aug94-PAS/2.11: cleanup/dimensioning
c  19Jun89-PAS: move ANGSET/extract ylmi prod/pack angsav
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output arrays:
      DIMENSION  vnlr(nk),vnli(nk)
c input arrays:
      DIMENSION  vnlint(limxp1,ljmxp1,lstmx)
      LOGICAL    izero(9,9)
      DIMENSION  ai(9),aj(9,9,*), ylmj(2,25,*)
      DIMENSION  angsav(limxp1,ljmxp1*ljmxp1,lmxp1,9,9)
      DIMENSION  coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  jlateq(neqlat)
c scratch arrays:
      DIMENSION  vnltmp(neqlat)
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
      call MKZERO( neqlat, vnltmp )
      do 200 lmpi=1,li
        lpi = lang(lmpi)
c
        do 100 lmpj=1,lj
          if( izero(lmpi,lmpj) ) goto 100
c
          aii = ai(lmpi)
          lpj = lang(lmpj)
c
          do 90 llpot=1,lmxp1
            lstore = lpi + lpj + 1 + lpowmx*(llpot-1)
c
            do 70 lli=1,limxp1
              jlx = 0
              jlm = 0
c
              mmj = 1
                do 41 llj=mmj,ljmxp1
                  jlm = jlm + 1
c                 jpprl = 1
                  jlx = jlx + 1
                  angfac = angsav(lli,jlx,llpot,lmpi,lmpj)
                  if( angfac.ne.zero )then
                    vt = aii*angfac*vnlint(lli,llj,lstore)
                    do  jeql=1,neqlat
                      vnltmp(jeql) = vnltmp(jeql) +
     $                 vt*aj(lmpj,lj,jeql)*ylmj(    1,jlm,jeql)
                    enddo
                  endif
   41           continue
c
              do 52 mmj=2,ljmxp1
                do 51 llj=mmj,ljmxp1
                  jlm = jlm + 1
c                   jpprl = 1
                    jlx = jlx + 1
                    angfac = angsav(lli,jlx,llpot,lmpi,lmpj)
                    if( angfac.ne.zero )then
                      vt = aii*angfac*vnlint(lli,llj,lstore)
                      do  jeql=1,neqlat
                        vnltmp(jeql) = vnltmp(jeql) +
     $                   vt*aj(lmpj,lj,jeql)*ylmj(    1,jlm,jeql)
                      enddo
                    endif
c                   jpprl = 2
                    jlx = jlx + 1
                    angfac = angsav(lli,jlx,llpot,lmpi,lmpj)
                    if( angfac.ne.zero )then
                      vt = aii*angfac*vnlint(lli,llj,lstore)
                      do  jeql=1,neqlat
                        vnltmp(jeql) = vnltmp(jeql) +
     $                   vt*aj(lmpj,lj,jeql)*ylmj(    2,jlm,jeql)
                      enddo
                    endif
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
      do  jeql=1,neqlat
        vint = vnltmp(jeql)
        jlat = jlateq(jeql)
        do  k=1,nk
          vnlr(k) = vnlr(k) + vint*coskr(k,jlat)
          vnli(k) = vnli(k) + vint*sinkr(k,jlat)
        enddo
      enddo
c
c    That's all Folks!
c
      RETURN
      END
