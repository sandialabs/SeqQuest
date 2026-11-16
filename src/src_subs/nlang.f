c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NLANG
c
c
      subroutine NLANG( lmi,lmj,limxp1,ljmxp1,lmxp1,rimag,rjmag,
     $ ylmi,gntcomb,angsav,gntci,cgntcj )
c---------------------------------------------------------------
c Purpose: compute sums of three-j , for non-local integrals
c
c Revision history:
c   6Mar00-PAS/2.43: rename "angset" to "nlang"
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c   5Jun89-PAS: extracted ANGSET from orbital loops
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c input arrays:
      DIMENSION  ylmi(2,25)
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
c output array:
      DIMENSION  angsav(limxp1,ljmxp1*ljmxp1,lmxp1,9,9)
c scratch arrays:
      DOUBLE COMPLEX    gntci(limxp1,lmi,lmxp1*lmxp1),
     $          cgntcj(ljmxp1*ljmxp1,lmj,lmxp1*lmxp1)
c local declarations:
      DIMENSION  ckm(9)
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      fourpi = pi + pi + pi + pi
c zero angsav:
      call MKZERO( limxp1*ljmxp1*ljmxp1*lmxp1*9*9, angsav )
c
c    Collect scale factors:
c
      do  llpot=1,lmxp1
        k = llpot - 1
        km0 = k*k
        numkm = 2*k + 1
        scale = fourpi* DBLE( numkm )
        do  km=km0+1,km0+numkm
          ckm(km) = scale
        enddo
      enddo
      maxkm = lmxp1*lmxp1
c
c    Pack gaunt coef * ylm's for i:
c
      call MKZERO( 2*limxp1*lmi*maxkm, gntci )
c
      ilm = 0
      iplst = 1
      do 100 mmi=1,limxp1
        do  lli=mmi,limxp1
          ilm = ilm + 1
          do  ipprl=1,iplst
c
            do  km=1,maxkm
              scale = ckm(km)*ylmi(ipprl,ilm)
              do  lmpi=1,lmi
                gntci(lli,lmpi,km) = gntci(lli,lmpi,km) +
     $           scale*gntcomb(lmpi,lli,mmi,ipprl,km)
              enddo
            enddo
c
          enddo
        enddo
  100 iplst = 2
c
c    Pack gaunt coef's for j:
c
      jlx = 0
      jplst = 1
      do 200 mmj=1,ljmxp1
        do  llj=mmj,ljmxp1
          do  jpprl=1,jplst
            jlx = jlx + 1
c
            do  km=1,maxkm
              do  lmpj=1,lmj
                cgntcj(jlx,lmpj,km) =
     $           CONJG( gntcomb(lmpj,llj,mmj,jpprl,km) )
              enddo
            enddo
c
          enddo
        enddo
  200 jplst=2
c
c    Do prods of gaunt coef's:
c
      do 302 lmpj=1,lmj
        do 301 lmpi=1,lmi
c
          do  llpot=1,lmxp1
            k = llpot - 1
            km0 = k*k
            numkm = 2*k + 1
c
            do  km=km0+1,km0+numkm
              do  lli=1,limxp1
                do  jlx=1,ljmxp1*ljmxp1
                  angsav(lli,jlx,llpot,lmpi,lmpj) =
     $             angsav(lli,jlx,llpot,lmpi,lmpj) +
     $             gntci(lli,lmpi,km)*cgntcj(jlx,lmpj,km)
                enddo
              enddo
            enddo
c
          enddo
c
  301   continue
  302 continue
c
c    That's all Folks!
c
      RETURN
      END
