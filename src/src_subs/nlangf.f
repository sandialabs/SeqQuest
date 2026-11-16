c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NLANGF
c
c
      subroutine NLANGF( lmi,lmj,limxp1,ljmxp1,lmxp1,lisq,ljsq,
     $ ylmi,ylmj,grylmi,grylmj,gntcomb,
     $ angsav,angtmp,gntci,cgntcj )
c---------------------------------------------------------------
c Purpose: do angular products and sums necessary for evaluation
c          of contribution to the force from the non-local part
c          of the potential.
c
c Written: Peter A. Schultz, 14-June-89,  based on PJF's
c          original ANGSETF,LLOOPSF,VLOOPSF routines.
c
c Revision history:
c  21Jun01-PAS/2.48: replace STOPs
c  15Oct99-PAS/2.41: cleaned up and minor tuning
c   5Oct99-DBR/2.41: optimization (unroll loops, reord dims)
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c  For later usage, need a variety of angular sums.  They are
c  stored in angsav() as follows:
c
c   angsav(1...) = sum of 3-j symbols * ylmj * ylmi
c   angsav(2...) = "" * ylmj * d/dx(ylmi) [=grylmi(i,1,ii)]
c   angsav(3...) = "" * ylmj * d/dy(ylmi) [=grylmi(i,2,ii)]
c   angsav(4...) = "" * ylmj * d/dz(ylmi) [=grylmi(i,3,ii)]
c   angsav(5...) = "" * d/dx(ylmj) * ylmi
c   angsav(6...) = "" * d/dy(ylmj) * ylmi
c   angsav(7...) = "" * d/dz(ylmj) * ylmi
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c input arrays:
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
c output array:
      DIMENSION  angsav(7,limxp1,ljmxp1,lmxp1,9,9)
c scratch arrays:
      DIMENSION  angtmp(lisq,ljsq)
      DOUBLE COMPLEX    gntci(lisq,lmi,lmxp1*lmxp1),
     $                 cgntcj(ljsq,lmj,lmxp1*lmxp1)
c local declarations:
      DIMENSION  ckm(9)
      DIMENSION  angij(25,25)
      DIMENSION  ang_yi(25,25), ang_yj(25,25)
      DIMENSION  ylmj(2,5,5), ylmi(2,5,5)
      DIMENSION  grylmj(3,2,5,5), grylmi(3,2,5,5)
      DIMENSION  ylmj_t(25), ylmi_t(25)
      DIMENSION  grylmj_t(25,3), grylmi_t(25,3)
      DIMENSION  angi(9),angj(9)
      DIMENSION  ime(5)
      DATA  ime/ 1,2,2,2,2/
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      if( limxp1 .gt. 5 .or. ljmxp1 .gt. 5 )then
        call STOPXERR( 'nlangf-0/NLANGFlmax exceeded' )
      endif
c
      fourpi = pi + pi + pi + pi
c zero angsav():
      call MKZERO( 7*limxp1*ljmxp1*lmxp1*9*9, angsav )
c
c    Collect scale factors:
c
      m = 0
      do  lli=1,limxp1
        do  mmi=1,lli
          do  ipprl=1,ime(mmi)
            m = m + 1
            ylmi_t(m) = ylmi(ipprl,mmi,lli)
            grylmi_t(m,1) = grylmi(1,ipprl,mmi,lli)
            grylmi_t(m,2) = grylmi(2,ipprl,mmi,lli)
            grylmi_t(m,3) = grylmi(3,ipprl,mmi,lli)
          enddo
        enddo
      enddo
c
      m = 0
      do  llj=1,ljmxp1
        do  mmj=1,llj
          do  jpprl=1,ime(mmj)
            m = m + 1
            ylmj_t(m) = ylmj(jpprl,mmj,llj)
            grylmj_t(m,1) = grylmj(1,jpprl,mmj,llj)
            grylmj_t(m,2) = grylmj(2,jpprl,mmj,llj)
            grylmj_t(m,3) = grylmj(3,jpprl,mmj,llj)
          enddo
        enddo
      enddo
c
      do  llpot=1,lmxp1
        k = llpot - 1
        km0 = k*k
        numkm = 2*k + 1
        scale = fourpi* DBLE(numkm)
        do  km=km0+1,km0+numkm
          ckm(km) = scale
        enddo
      enddo
      maxkm = lmxp1*lmxp1
c
      ilx = 0
      do  lli=1,limxp1
        do  mmi=1,lli
          do  ipprl=1,ime(mmi)
            ilx = ilx + 1
c
            do  km=1,maxkm
              scale = ckm(km)
              do  lmpi=1,lmi
                gntci(ilx,lmpi,km) = scale*
     $           gntcomb(lmpi,lli,mmi,ipprl,km)
              enddo
            enddo
c
          enddo
        enddo
      enddo
c
c zero out cgntcj():
      call MKZERO( 2*ljsq*lmj*maxkm, cgntcj )
c
      jlx = 0
      do  llj=1,ljmxp1
        do  mmj=1,llj
          do  jpprl=1,ime(mmj)
            jlx = jlx + 1
c
            do  km=1,maxkm
              do  lmpj=1,lmj
                cgntcj(jlx,lmpj,km) =
     $           CONJG(gntcomb(lmpj,llj,mmj,jpprl,km))
              enddo
            enddo
c
          enddo
        enddo
      enddo
c
c    Generate products of gnt coef's and factor in ylm's,grylm's:
c
      do 302 lmpj=1,lmj
        do 301 lmpi=1,lmi
c
          do 290 llpot=1,lmxp1
            k = llpot - 1
            km0 = k*k
            numkm = 2*k + 1
c
            call MKZERO( 25*25, angij )
            do  km=km0+1,km0+numkm
              do  jlx=1,ljsq
                do  ilx=1,lisq
                  angij(ilx,jlx) = angij(ilx,jlx) +
     $             gntci(ilx,lmpi,km)*cgntcj(jlx,lmpj,km)
                enddo
              enddo
            enddo
c
            if( lisq .ge. ljsq )then
              do  jlx=1,ljsq
                ylj = ylmj_t(jlx)
                do  ilx=1,lisq
                  ang_yj(jlx,ilx) = angij(ilx,jlx)*ylj
                  ang_yi(ilx,jlx) = angij(ilx,jlx)*ylmi_t(ilx)
                enddo
              enddo
            else
              do  ilx=1,lisq
                yli = ylmi_t(ilx)
                do  jlx=1,ljsq
                  ang_yj(jlx,ilx) = angij(ilx,jlx)*ylmj_t(jlx)
                  ang_yi(ilx,jlx) = angij(ilx,jlx)*yli
                enddo
              enddo
            endif
c
            do  llj=1,ljmxp1
c
              if( llj .eq. 1 )then
c
                do  lli=1,limxp1
c
                  if( lli .eq. 1 )then
c
                    angi(1) = ang_yi(1,1)
c
                    angj(1) = ang_yj(1,1)
c
                    angsav(1,1,1,llpot,lmpi,lmpj) =
     $               ylmi_t(1)*angj(1)
c
                    angsav(2,1,1,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,1)
                    angsav(3,1,1,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,2)
                    angsav(4,1,1,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,3)
c
                    angsav(5,1,1,llpot,lmpi,lmpj) =
     $               angi(1)*grylmj_t(1,1)
                    angsav(6,1,1,llpot,lmpi,lmpj) =
     $               angi(1)*grylmj_t(1,2)
                    angsav(7,1,1,llpot,lmpi,lmpj) =
     $               angi(1)*grylmj_t(1,3)
c
                  elseif( lli .eq. 2 )then
c
                    angi(1) = ang_yi(2,1)
     $                      + ang_yi(3,1)
     $                      + ang_yi(4,1)
c
                    angsav(1,2,1,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(1)
c
                    angsav(2,2,1,llpot,lmpi,lmpj) =
     $               ang_yj(1,2)*grylmi_t(2,1) +
     $               ang_yj(1,3)*grylmi_t(3,1) +
     $               ang_yj(1,4)*grylmi_t(4,1)
                    angsav(3,2,1,llpot,lmpi,lmpj) =
     $               ang_yj(1,2)*grylmi_t(2,2) +
     $               ang_yj(1,3)*grylmi_t(3,2) +
     $               ang_yj(1,4)*grylmi_t(4,2)
                    angsav(4,2,1,llpot,lmpi,lmpj) =
     $               ang_yj(1,2)*grylmi_t(2,3) +
     $               ang_yj(1,3)*grylmi_t(3,3) +
     $               ang_yj(1,4)*grylmi_t(4,3)
c
                    angsav(5,2,1,llpot,lmpi,lmpj) =
     $               grylmj_t(1,1)*angi(1)
                    angsav(6,2,1,llpot,lmpi,lmpj) =
     $               grylmj_t(1,2)*angi(1)
                    angsav(7,2,1,llpot,lmpi,lmpj) =
     $               grylmj_t(1,3)*angi(1)
c
                  elseif( lli .eq. 3 )then
c
                    angi(1) = ang_yi(5,1)
     $                      + ang_yi(6,1)
     $                      + ang_yi(7,1)
     $                      + ang_yi(8,1)
     $                      + ang_yi(9,1)
c
                    angsav(1,3,1,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(1)
c
                    angsav(2,3,1,llpot,lmpi,lmpj) =
     $               ang_yj(1,5)*grylmi_t(5,1) +
     $               ang_yj(1,6)*grylmi_t(6,1) +
     $               ang_yj(1,7)*grylmi_t(7,1) +
     $               ang_yj(1,8)*grylmi_t(8,1) +
     $               ang_yj(1,9)*grylmi_t(9,1)
                    angsav(3,3,1,llpot,lmpi,lmpj) =
     $               ang_yj(1,5)*grylmi_t(5,2) +
     $               ang_yj(1,6)*grylmi_t(6,2) +
     $               ang_yj(1,7)*grylmi_t(7,2) +
     $               ang_yj(1,8)*grylmi_t(8,2) +
     $               ang_yj(1,9)*grylmi_t(9,2)
                    angsav(4,3,1,llpot,lmpi,lmpj) =
     $               ang_yj(1,5)*grylmi_t(5,3) +
     $               ang_yj(1,6)*grylmi_t(6,3) +
     $               ang_yj(1,7)*grylmi_t(7,3) +
     $               ang_yj(1,8)*grylmi_t(8,3) +
     $               ang_yj(1,9)*grylmi_t(9,3)
c
                    angsav(5,3,1,llpot,lmpi,lmpj) =
     $               grylmj_t(1,1)*angi(1)
                    angsav(6,3,1,llpot,lmpi,lmpj) =
     $               grylmj_t(1,2)*angi(1)
                    angsav(7,3,1,llpot,lmpi,lmpj) =
     $               grylmj_t(1,3)*angi(1)
c
                  elseif( lli .eq. 4 )then
c
                    angi(1) = ang_yi(10,1)
     $                      + ang_yi(11,1)
     $                      + ang_yi(12,1)
     $                      + ang_yi(13,1)
     $                      + ang_yi(14,1)
     $                      + ang_yi(15,1)
     $                      + ang_yi(16,1)
c
                    angsav(1,4,1,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(1)
c
                    angsav(2,4,1,llpot,lmpi,lmpj) =
     $               ang_yj(1,10)*grylmi_t(10,1) +
     $               ang_yj(1,11)*grylmi_t(11,1) +
     $               ang_yj(1,12)*grylmi_t(12,1) +
     $               ang_yj(1,13)*grylmi_t(13,1) +
     $               ang_yj(1,14)*grylmi_t(14,1) +
     $               ang_yj(1,15)*grylmi_t(15,1) +
     $               ang_yj(1,16)*grylmi_t(16,1)
                    angsav(3,4,1,llpot,lmpi,lmpj) =
     $               ang_yj(1,10)*grylmi_t(10,2) +
     $               ang_yj(1,11)*grylmi_t(11,2) +
     $               ang_yj(1,12)*grylmi_t(12,2) +
     $               ang_yj(1,13)*grylmi_t(13,2) +
     $               ang_yj(1,14)*grylmi_t(14,2) +
     $               ang_yj(1,15)*grylmi_t(15,2) +
     $               ang_yj(1,16)*grylmi_t(16,2)
                    angsav(4,4,1,llpot,lmpi,lmpj) =
     $               ang_yj(1,10)*grylmi_t(10,3) +
     $               ang_yj(1,11)*grylmi_t(11,3) +
     $               ang_yj(1,12)*grylmi_t(12,3) +
     $               ang_yj(1,13)*grylmi_t(13,3) +
     $               ang_yj(1,14)*grylmi_t(14,3) +
     $               ang_yj(1,15)*grylmi_t(15,3) +
     $               ang_yj(1,16)*grylmi_t(16,3)
c
                    angsav(5,4,1,llpot,lmpi,lmpj) =
     $               grylmj_t(1,1)*angi(1)
                    angsav(6,4,1,llpot,lmpi,lmpj) =
     $               grylmj_t(1,2)*angi(1)
                    angsav(7,4,1,llpot,lmpi,lmpj) =
     $               grylmj_t(1,3)*angi(1)
c
                  elseif( lli .eq. 5 )then
c
                    angi(1) = ang_yi(17,1)
     $                      + ang_yi(18,1)
     $                      + ang_yi(19,1)
     $                      + ang_yi(20,1)
     $                      + ang_yi(21,1)
     $                      + ang_yi(22,1)
     $                      + ang_yi(23,1)
     $                      + ang_yi(24,1)
     $                      + ang_yi(25,1)
c
                    angsav(1,5,1,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(1)
c
                    angsav(2,5,1,llpot,lmpi,lmpj) =
     $               ang_yj(1,17)*grylmi_t(17,1) +
     $               ang_yj(1,18)*grylmi_t(18,1) +
     $               ang_yj(1,19)*grylmi_t(19,1) +
     $               ang_yj(1,20)*grylmi_t(20,1) +
     $               ang_yj(1,21)*grylmi_t(21,1) +
     $               ang_yj(1,22)*grylmi_t(22,1) +
     $               ang_yj(1,23)*grylmi_t(23,1) +
     $               ang_yj(1,24)*grylmi_t(24,1) +
     $               ang_yj(1,25)*grylmi_t(25,1)
                    angsav(3,5,1,llpot,lmpi,lmpj) =
     $               ang_yj(1,17)*grylmi_t(17,2) +
     $               ang_yj(1,18)*grylmi_t(18,2) +
     $               ang_yj(1,19)*grylmi_t(19,2) +
     $               ang_yj(1,20)*grylmi_t(20,2) +
     $               ang_yj(1,21)*grylmi_t(21,2) +
     $               ang_yj(1,22)*grylmi_t(22,2) +
     $               ang_yj(1,23)*grylmi_t(23,2) +
     $               ang_yj(1,24)*grylmi_t(24,2) +
     $               ang_yj(1,25)*grylmi_t(25,2)
                    angsav(4,5,1,llpot,lmpi,lmpj) =
     $               ang_yj(1,17)*grylmi_t(17,3) +
     $               ang_yj(1,18)*grylmi_t(18,3) +
     $               ang_yj(1,19)*grylmi_t(19,3) +
     $               ang_yj(1,20)*grylmi_t(20,3) +
     $               ang_yj(1,21)*grylmi_t(21,3) +
     $               ang_yj(1,22)*grylmi_t(22,3) +
     $               ang_yj(1,23)*grylmi_t(23,3) +
     $               ang_yj(1,24)*grylmi_t(24,3) +
     $               ang_yj(1,25)*grylmi_t(25,3)
c
                    angsav(5,5,1,llpot,lmpi,lmpj) =
     $               grylmj_t(1,1)*angi(1)
                    angsav(6,5,1,llpot,lmpi,lmpj) =
     $               grylmj_t(1,2)*angi(1)
                    angsav(7,5,1,llpot,lmpi,lmpj) =
     $               grylmj_t(1,3)*angi(1)
c
                  else
c
                    call STOPXERR( 'nlangf 1' )
c
                  endif
c
                enddo
c
              elseif( llj .eq. 2 )then
c
                do  lli=1,limxp1
c
                  if( lli .eq. 1 )then
c
                    angj(1) = ang_yj(2,1)
     $                      + ang_yj(3,1)
     $                      + ang_yj(4,1)
c
                    angsav(1,1,2,llpot,lmpi,lmpj) =
     $               ylmi_t(1)*angj(1)
c
                    angsav(2,1,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,1)
                    angsav(3,1,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,2)
                    angsav(4,1,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,3)
c
                    angsav(5,1,2,llpot,lmpi,lmpj) =
     $               ang_yi(1,2)*grylmj_t(2,1) +
     $               ang_yi(1,3)*grylmj_t(3,1) +
     $               ang_yi(1,4)*grylmj_t(4,1)
                    angsav(6,1,2,llpot,lmpi,lmpj) =
     $               ang_yi(1,2)*grylmj_t(2,2) +
     $               ang_yi(1,3)*grylmj_t(3,2) +
     $               ang_yi(1,4)*grylmj_t(4,2)
                    angsav(7,1,2,llpot,lmpi,lmpj) =
     $               ang_yi(1,2)*grylmj_t(2,3) +
     $               ang_yi(1,3)*grylmj_t(3,3) +
     $               ang_yi(1,4)*grylmj_t(4,3)
c
                  elseif( lli .eq. 2 )then
c
                    angi(1) = ang_yi(2,2)
     $                      + ang_yi(3,2)
     $                      + ang_yi(4,2)
                    angi(2) = ang_yi(2,3)
     $                      + ang_yi(3,3)
     $                      + ang_yi(4,3)
                    angi(3) = ang_yi(2,4)
     $                      + ang_yi(3,4)
     $                      + ang_yi(4,4)
c
                    angj(1) = ang_yj(2,2)
     $                      + ang_yj(3,2)
     $                      + ang_yj(4,2)
                    angj(2) = ang_yj(2,3)
     $                      + ang_yj(3,3)
     $                      + ang_yj(4,3)
                    angj(3) = ang_yj(2,4)
     $                      + ang_yj(3,4)
     $                      + ang_yj(4,4)
c
                    angsav(1,2,2,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(2) +
     $               angi(2)*ylmj_t(3) +
     $               angi(3)*ylmj_t(4)
c
                    angsav(2,2,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(2,1) +
     $               angj(2)*grylmi_t(3,1) +
     $               angj(3)*grylmi_t(4,1)
                    angsav(3,2,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(2,2) +
     $               angj(2)*grylmi_t(3,2) +
     $               angj(3)*grylmi_t(4,2)
                    angsav(4,2,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(2,3) +
     $               angj(2)*grylmi_t(3,3) +
     $               angj(3)*grylmi_t(4,3)
c
                    angsav(5,2,2,llpot,lmpi,lmpj) =
     $               grylmj_t(2,1)*angi(1) +
     $               grylmj_t(3,1)*angi(2) +
     $               grylmj_t(4,1)*angi(3)
                    angsav(6,2,2,llpot,lmpi,lmpj) =
     $               grylmj_t(2,2)*angi(1) +
     $               grylmj_t(3,2)*angi(2) +
     $               grylmj_t(4,2)*angi(3)
                    angsav(7,2,2,llpot,lmpi,lmpj) =
     $               grylmj_t(2,3)*angi(1) +
     $               grylmj_t(3,3)*angi(2) +
     $               grylmj_t(4,3)*angi(3)
c
                  elseif( lli .eq. 3 )then
c
                    angi(1) = ang_yi(5,2)
     $                      + ang_yi(6,2)
     $                      + ang_yi(7,2)
     $                      + ang_yi(8,2)
     $                      + ang_yi(9,2)
                    angi(2) = ang_yi(5,3)
     $                      + ang_yi(6,3)
     $                      + ang_yi(7,3)
     $                      + ang_yi(8,3)
     $                      + ang_yi(9,3)
                    angi(3) = ang_yi(5,4)
     $                      + ang_yi(6,4)
     $                      + ang_yi(7,4)
     $                      + ang_yi(8,4)
     $                      + ang_yi(9,4)
c
                    angj(1) = ang_yj(2,5)
     $                      + ang_yj(3,5)
     $                      + ang_yj(4,5)
                    angj(2) = ang_yj(2,6)
     $                      + ang_yj(3,6)
     $                      + ang_yj(4,6)
                    angj(3) = ang_yj(2,7)
     $                      + ang_yj(3,7)
     $                      + ang_yj(4,7)
                    angj(4) = ang_yj(2,8)
     $                      + ang_yj(3,8)
     $                      + ang_yj(4,8)
                    angj(5) = ang_yj(2,9)
     $                      + ang_yj(3,9)
     $                      + ang_yj(4,9)
c
                    angsav(1,3,2,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(2) +
     $               angi(2)*ylmj_t(3) +
     $               angi(3)*ylmj_t(4)
c
                    angsav(2,3,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(5,1) +
     $               angj(2)*grylmi_t(6,1) +
     $               angj(3)*grylmi_t(7,1) +
     $               angj(4)*grylmi_t(8,1) +
     $               angj(5)*grylmi_t(9,1)
                    angsav(3,3,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(5,2) +
     $               angj(2)*grylmi_t(6,2) +
     $               angj(3)*grylmi_t(7,2) +
     $               angj(4)*grylmi_t(8,2) +
     $               angj(5)*grylmi_t(9,2)
                    angsav(4,3,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(5,3) +
     $               angj(2)*grylmi_t(6,3) +
     $               angj(3)*grylmi_t(7,3) +
     $               angj(4)*grylmi_t(8,3) +
     $               angj(5)*grylmi_t(9,3)
c
                    angsav(5,3,2,llpot,lmpi,lmpj) =
     $               grylmj_t(2,1)*angi(1) +
     $               grylmj_t(3,1)*angi(2) +
     $               grylmj_t(4,1)*angi(3)
                    angsav(6,3,2,llpot,lmpi,lmpj) =
     $               grylmj_t(2,2)*angi(1) +
     $               grylmj_t(3,2)*angi(2) +
     $               grylmj_t(4,2)*angi(3)
                    angsav(7,3,2,llpot,lmpi,lmpj) =
     $               grylmj_t(2,3)*angi(1) +
     $               grylmj_t(3,3)*angi(2) +
     $               grylmj_t(4,3)*angi(3)
c
                  elseif( lli .eq. 4 )then
c
                    angi(1) = ang_yi(10,2)
     $                      + ang_yi(11,2)
     $                      + ang_yi(12,2)
     $                      + ang_yi(13,2)
     $                      + ang_yi(14,2)
     $                      + ang_yi(15,2)
     $                      + ang_yi(16,2)
                    angi(2) = ang_yi(10,3)
     $                      + ang_yi(11,3)
     $                      + ang_yi(12,3)
     $                      + ang_yi(13,3)
     $                      + ang_yi(14,3)
     $                      + ang_yi(15,3)
     $                      + ang_yi(16,3)
                    angi(3) = ang_yi(10,4)
     $                      + ang_yi(11,4)
     $                      + ang_yi(12,4)
     $                      + ang_yi(13,4)
     $                      + ang_yi(14,4)
     $                      + ang_yi(15,4)
     $                      + ang_yi(16,4)
c
                    angj(1) = ang_yj(2,10)
     $                      + ang_yj(3,10)
     $                      + ang_yj(4,10)
                    angj(2) = ang_yj(2,11)
     $                      + ang_yj(3,11)
     $                      + ang_yj(4,11)
                    angj(3) = ang_yj(2,12)
     $                      + ang_yj(3,12)
     $                      + ang_yj(4,12)
                    angj(4) = ang_yj(2,13)
     $                      + ang_yj(3,13)
     $                      + ang_yj(4,13)
                    angj(5) = ang_yj(2,14)
     $                      + ang_yj(3,14)
     $                      + ang_yj(4,14)
                    angj(6) = ang_yj(2,15)
     $                      + ang_yj(3,15)
     $                      + ang_yj(4,15)
                    angj(7) = ang_yj(2,16)
     $                      + ang_yj(3,16)
     $                      + ang_yj(4,16)
c
                    angsav(1,4,2,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(2) +
     $               angi(2)*ylmj_t(3) +
     $               angi(3)*ylmj_t(4)
c
                    angsav(2,4,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(10,1) +
     $               angj(2)*grylmi_t(11,1) +
     $               angj(3)*grylmi_t(12,1) +
     $               angj(4)*grylmi_t(13,1) +
     $               angj(5)*grylmi_t(14,1) +
     $               angj(6)*grylmi_t(15,1) +
     $               angj(7)*grylmi_t(16,1)
                    angsav(3,4,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(10,2) +
     $               angj(2)*grylmi_t(11,2) +
     $               angj(3)*grylmi_t(12,2) +
     $               angj(4)*grylmi_t(13,2) +
     $               angj(5)*grylmi_t(14,2) +
     $               angj(6)*grylmi_t(15,2) +
     $               angj(7)*grylmi_t(16,2)
                    angsav(4,4,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(10,3) +
     $               angj(2)*grylmi_t(11,3) +
     $               angj(3)*grylmi_t(12,3) +
     $               angj(4)*grylmi_t(13,3) +
     $               angj(5)*grylmi_t(14,3) +
     $               angj(6)*grylmi_t(15,3) +
     $               angj(7)*grylmi_t(16,3)
c
                    angsav(5,4,2,llpot,lmpi,lmpj) =
     $               grylmj_t(2,1)*angi(1) +
     $               grylmj_t(3,1)*angi(2) +
     $               grylmj_t(4,1)*angi(3)
                    angsav(6,4,2,llpot,lmpi,lmpj) =
     $               grylmj_t(2,2)*angi(1) +
     $               grylmj_t(3,2)*angi(2) +
     $               grylmj_t(4,2)*angi(3)
                    angsav(7,4,2,llpot,lmpi,lmpj) =
     $               grylmj_t(2,3)*angi(1) +
     $               grylmj_t(3,3)*angi(2) +
     $               grylmj_t(4,3)*angi(3)
c
                  elseif( lli .eq. 5 )then
c
                    angi(1) = ang_yi(17,2)
     $                      + ang_yi(18,2)
     $                      + ang_yi(19,2)
     $                      + ang_yi(20,2)
     $                      + ang_yi(21,2)
     $                      + ang_yi(22,2)
     $                      + ang_yi(23,2)
     $                      + ang_yi(24,2)
     $                      + ang_yi(25,2)
                    angi(2) = ang_yi(17,3)
     $                      + ang_yi(18,3)
     $                      + ang_yi(19,3)
     $                      + ang_yi(20,3)
     $                      + ang_yi(21,3)
     $                      + ang_yi(22,3)
     $                      + ang_yi(23,3)
     $                      + ang_yi(24,3)
     $                      + ang_yi(25,3)
                    angi(3) = ang_yi(17,4)
     $                      + ang_yi(18,4)
     $                      + ang_yi(19,4)
     $                      + ang_yi(20,4)
     $                      + ang_yi(21,4)
     $                      + ang_yi(22,4)
     $                      + ang_yi(23,4)
     $                      + ang_yi(24,4)
     $                      + ang_yi(25,4)
c
                    angj(1) = ang_yj(2,17)
     $                      + ang_yj(3,17)
     $                      + ang_yj(4,17)
                    angj(2) = ang_yj(2,18)
     $                      + ang_yj(3,18)
     $                      + ang_yj(4,18)
                    angj(3) = ang_yj(2,19)
     $                      + ang_yj(3,19)
     $                      + ang_yj(4,19)
                    angj(4) = ang_yj(2,20)
     $                      + ang_yj(3,20)
     $                      + ang_yj(4,20)
                    angj(5) = ang_yj(2,21)
     $                      + ang_yj(3,21)
     $                      + ang_yj(4,21)
                    angj(6) = ang_yj(2,22)
     $                      + ang_yj(3,22)
     $                      + ang_yj(4,22)
                    angj(7) = ang_yj(2,23)
     $                      + ang_yj(3,23)
     $                      + ang_yj(4,23)
                    angj(8) = ang_yj(2,24)
     $                      + ang_yj(3,24)
     $                      + ang_yj(4,24)
                    angj(9) = ang_yj(2,25)
     $                      + ang_yj(3,25)
     $                      + ang_yj(4,25)
c
                    angsav(1,5,2,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(2) +
     $               angi(2)*ylmj_t(3) +
     $               angi(3)*ylmj_t(4)
c
                    angsav(2,5,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(17,1) +
     $               angj(2)*grylmi_t(18,1) +
     $               angj(3)*grylmi_t(19,1) +
     $               angj(4)*grylmi_t(20,1) +
     $               angj(5)*grylmi_t(21,1) +
     $               angj(6)*grylmi_t(22,1) +
     $               angj(7)*grylmi_t(23,1) +
     $               angj(8)*grylmi_t(24,1) +
     $               angj(9)*grylmi_t(25,1)
                    angsav(3,5,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(17,2) +
     $               angj(2)*grylmi_t(18,2) +
     $               angj(3)*grylmi_t(19,2) +
     $               angj(4)*grylmi_t(20,2) +
     $               angj(5)*grylmi_t(21,2) +
     $               angj(6)*grylmi_t(22,2) +
     $               angj(7)*grylmi_t(23,2) +
     $               angj(8)*grylmi_t(24,2) +
     $               angj(9)*grylmi_t(25,2)
                    angsav(4,5,2,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(17,3) +
     $               angj(2)*grylmi_t(18,3) +
     $               angj(3)*grylmi_t(19,3) +
     $               angj(4)*grylmi_t(20,3) +
     $               angj(5)*grylmi_t(21,3) +
     $               angj(6)*grylmi_t(22,3) +
     $               angj(7)*grylmi_t(23,3) +
     $               angj(8)*grylmi_t(24,3) +
     $               angj(9)*grylmi_t(25,3)
c
                    angsav(5,5,2,llpot,lmpi,lmpj) =
     $               grylmj_t(2,1)*angi(1) +
     $               grylmj_t(3,1)*angi(2) +
     $               grylmj_t(4,1)*angi(3)
                    angsav(6,5,2,llpot,lmpi,lmpj) =
     $               grylmj_t(2,2)*angi(1) +
     $               grylmj_t(3,2)*angi(2) +
     $               grylmj_t(4,2)*angi(3)
                    angsav(7,5,2,llpot,lmpi,lmpj) =
     $               grylmj_t(2,3)*angi(1) +
     $               grylmj_t(3,3)*angi(2) +
     $               grylmj_t(4,3)*angi(3)
c
                  else
c
                    call STOPXERR( 'nlangf-2' )
c
                  endif
c
                enddo
c
              elseif( llj .eq. 3 )then
c
                do  lli=1,limxp1
c
                  if( lli .eq. 1 )then
c
                    angj(1) = ang_yj(5,1)
     $                      + ang_yj(6,1)
     $                      + ang_yj(7,1)
     $                      + ang_yj(8,1)
     $                      + ang_yj(9,1)
c
                    angsav(1,1,3,llpot,lmpi,lmpj) =
     $               ylmi_t(1)*angj(1)
c
                    angsav(2,1,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,1)
                    angsav(3,1,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,2)
                    angsav(4,1,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,3)
c
                    angsav(5,1,3,llpot,lmpi,lmpj) =
     $               ang_yi(1,5)*grylmj_t(5,1) +
     $               ang_yi(1,6)*grylmj_t(6,1) +
     $               ang_yi(1,7)*grylmj_t(7,1) +
     $               ang_yi(1,8)*grylmj_t(8,1) +
     $               ang_yi(1,9)*grylmj_t(9,1)
                    angsav(6,1,3,llpot,lmpi,lmpj) =
     $               ang_yi(1,5)*grylmj_t(5,2) +
     $               ang_yi(1,6)*grylmj_t(6,2) +
     $               ang_yi(1,7)*grylmj_t(7,2) +
     $               ang_yi(1,8)*grylmj_t(8,2) +
     $               ang_yi(1,9)*grylmj_t(9,2)
                    angsav(7,1,3,llpot,lmpi,lmpj) =
     $               ang_yi(1,5)*grylmj_t(5,3) +
     $               ang_yi(1,6)*grylmj_t(6,3) +
     $               ang_yi(1,7)*grylmj_t(7,3) +
     $               ang_yi(1,8)*grylmj_t(8,3) +
     $               ang_yi(1,9)*grylmj_t(9,3)
c
                  elseif( lli .eq. 2 )then
c
                    angi(1) = ang_yi(2,5)
     $                      + ang_yi(3,5)
     $                      + ang_yi(4,5)
                    angi(2) = ang_yi(2,6)
     $                      + ang_yi(3,6)
     $                      + ang_yi(4,6)
                    angi(3) = ang_yi(2,7)
     $                      + ang_yi(3,7)
     $                      + ang_yi(4,7)
                    angi(4) = ang_yi(2,8)
     $                      + ang_yi(3,8)
     $                      + ang_yi(4,8)
                    angi(5) = ang_yi(2,9)
     $                      + ang_yi(3,9)
     $                      + ang_yi(4,9)
c
                    angj(1) = ang_yj(5,2)
     $                      + ang_yj(6,2)
     $                      + ang_yj(7,2)
     $                      + ang_yj(8,2)
     $                      + ang_yj(9,2)
                    angj(2) = ang_yj(5,3)
     $                      + ang_yj(6,3)
     $                      + ang_yj(7,3)
     $                      + ang_yj(8,3)
     $                      + ang_yj(9,3)
                    angj(3) = ang_yj(5,4)
     $                      + ang_yj(6,4)
     $                      + ang_yj(7,4)
     $                      + ang_yj(8,4)
     $                      + ang_yj(9,4)
c
                    angsav(1,2,3,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(5) +
     $               angi(2)*ylmj_t(6) +
     $               angi(3)*ylmj_t(7) +
     $               angi(4)*ylmj_t(8) +
     $               angi(5)*ylmj_t(9)
c
                    angsav(2,2,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(2,1) +
     $               angj(2)*grylmi_t(3,1) +
     $               angj(3)*grylmi_t(4,1)
                    angsav(3,2,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(2,2) +
     $               angj(2)*grylmi_t(3,2) +
     $               angj(3)*grylmi_t(4,2)
                    angsav(4,2,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(2,3) +
     $               angj(2)*grylmi_t(3,3) +
     $               angj(3)*grylmi_t(4,3)
c
                    angsav(5,2,3,llpot,lmpi,lmpj) =
     $               grylmj_t(5,1)*angi(1) +
     $               grylmj_t(6,1)*angi(2) +
     $               grylmj_t(7,1)*angi(3) +
     $               grylmj_t(8,1)*angi(4) +
     $               grylmj_t(9,1)*angi(5)
                    angsav(6,2,3,llpot,lmpi,lmpj) =
     $               grylmj_t(5,2)*angi(1) +
     $               grylmj_t(6,2)*angi(2) +
     $               grylmj_t(7,2)*angi(3) +
     $               grylmj_t(8,2)*angi(4) +
     $               grylmj_t(9,2)*angi(5)
                    angsav(7,2,3,llpot,lmpi,lmpj) =
     $               grylmj_t(5,3)*angi(1) +
     $               grylmj_t(6,3)*angi(2) +
     $               grylmj_t(7,3)*angi(3) +
     $               grylmj_t(8,3)*angi(4) +
     $               grylmj_t(9,3)*angi(5)
c
                  elseif( lli .eq. 3 )then
c
                    angi(1) = ang_yi(5,5)
     $                      + ang_yi(6,5)
     $                      + ang_yi(7,5)
     $                      + ang_yi(8,5)
     $                      + ang_yi(9,5)
                    angi(2) = ang_yi(5,6)
     $                      + ang_yi(6,6)
     $                      + ang_yi(7,6)
     $                      + ang_yi(8,6)
     $                      + ang_yi(9,6)
                    angi(3) = ang_yi(5,7)
     $                      + ang_yi(6,7)
     $                      + ang_yi(7,7)
     $                      + ang_yi(8,7)
     $                      + ang_yi(9,7)
                    angi(4) = ang_yi(5,8)
     $                      + ang_yi(6,8)
     $                      + ang_yi(7,8)
     $                      + ang_yi(8,8)
     $                      + ang_yi(9,8)
                    angi(5) = ang_yi(5,9)
     $                      + ang_yi(6,9)
     $                      + ang_yi(7,9)
     $                      + ang_yi(8,9)
     $                      + ang_yi(9,9)
c
                    angj(1) = ang_yj(5,5)
     $                      + ang_yj(6,5)
     $                      + ang_yj(7,5)
     $                      + ang_yj(8,5)
     $                      + ang_yj(9,5)
                    angj(2) = ang_yj(5,6)
     $                      + ang_yj(6,6)
     $                      + ang_yj(7,6)
     $                      + ang_yj(8,6)
     $                      + ang_yj(9,6)
                    angj(3) = ang_yj(5,7)
     $                      + ang_yj(6,7)
     $                      + ang_yj(7,7)
     $                      + ang_yj(8,7)
     $                      + ang_yj(9,7)
                    angj(4) = ang_yj(5,8)
     $                      + ang_yj(6,8)
     $                      + ang_yj(7,8)
     $                      + ang_yj(8,8)
     $                      + ang_yj(9,8)
                    angj(5) = ang_yj(5,9)
     $                      + ang_yj(6,9)
     $                      + ang_yj(7,9)
     $                      + ang_yj(8,9)
     $                      + ang_yj(9,9)
c
                    angsav(1,3,3,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(5) +
     $               angi(2)*ylmj_t(6) +
     $               angi(3)*ylmj_t(7) +
     $               angi(4)*ylmj_t(8) +
     $               angi(5)*ylmj_t(9)
c
                    angsav(2,3,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(5,1) +
     $               angj(2)*grylmi_t(6,1) +
     $               angj(3)*grylmi_t(7,1) +
     $               angj(4)*grylmi_t(8,1) +
     $               angj(5)*grylmi_t(9,1)
                    angsav(3,3,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(5,2) +
     $               angj(2)*grylmi_t(6,2) +
     $               angj(3)*grylmi_t(7,2) +
     $               angj(4)*grylmi_t(8,2) +
     $               angj(5)*grylmi_t(9,2)
                    angsav(4,3,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(5,3) +
     $               angj(2)*grylmi_t(6,3) +
     $               angj(3)*grylmi_t(7,3) +
     $               angj(4)*grylmi_t(8,3) +
     $               angj(5)*grylmi_t(9,3)
c
                    angsav(5,3,3,llpot,lmpi,lmpj) =
     $               grylmj_t(5,1)*angi(1) +
     $               grylmj_t(6,1)*angi(2) +
     $               grylmj_t(7,1)*angi(3) +
     $               grylmj_t(8,1)*angi(4) +
     $               grylmj_t(9,1)*angi(5)
                    angsav(6,3,3,llpot,lmpi,lmpj) =
     $               grylmj_t(5,2)*angi(1) +
     $               grylmj_t(6,2)*angi(2) +
     $               grylmj_t(7,2)*angi(3) +
     $               grylmj_t(8,2)*angi(4) +
     $               grylmj_t(9,2)*angi(5)
                    angsav(7,3,3,llpot,lmpi,lmpj) =
     $               grylmj_t(5,3)*angi(1) +
     $               grylmj_t(6,3)*angi(2) +
     $               grylmj_t(7,3)*angi(3) +
     $               grylmj_t(8,3)*angi(4) +
     $               grylmj_t(9,3)*angi(5)
c
                  elseif( lli .eq. 4 )then
c
                    angi(1) = ang_yi(10,5)
     $                      + ang_yi(11,5)
     $                      + ang_yi(12,5)
     $                      + ang_yi(13,5)
     $                      + ang_yi(14,5)
     $                      + ang_yi(15,5)
     $                      + ang_yi(16,5)
                    angi(2) = ang_yi(10,6)
     $                      + ang_yi(11,6)
     $                      + ang_yi(12,6)
     $                      + ang_yi(13,6)
     $                      + ang_yi(14,6)
     $                      + ang_yi(15,6)
     $                      + ang_yi(16,6)
                    angi(3) = ang_yi(10,7)
     $                      + ang_yi(11,7)
     $                      + ang_yi(12,7)
     $                      + ang_yi(13,7)
     $                      + ang_yi(14,7)
     $                      + ang_yi(15,7)
     $                      + ang_yi(16,7)
                    angi(4) = ang_yi(10,8)
     $                      + ang_yi(11,8)
     $                      + ang_yi(12,8)
     $                      + ang_yi(13,8)
     $                      + ang_yi(14,8)
     $                      + ang_yi(15,8)
     $                      + ang_yi(16,8)
                    angi(5) = ang_yi(10,9)
     $                      + ang_yi(11,9)
     $                      + ang_yi(12,9)
     $                      + ang_yi(13,9)
     $                      + ang_yi(14,9)
     $                      + ang_yi(15,9)
     $                      + ang_yi(16,9)
c
                    angj(1) = ang_yj(5,10)
     $                      + ang_yj(6,10)
     $                      + ang_yj(7,10)
     $                      + ang_yj(8,10)
     $                      + ang_yj(9,10)
                    angj(2) = ang_yj(5,11)
     $                      + ang_yj(6,11)
     $                      + ang_yj(7,11)
     $                      + ang_yj(8,11)
     $                      + ang_yj(9,11)
                    angj(3) = ang_yj(5,12)
     $                      + ang_yj(6,12)
     $                      + ang_yj(7,12)
     $                      + ang_yj(8,12)
     $                      + ang_yj(9,12)
                    angj(4) = ang_yj(5,13)
     $                      + ang_yj(6,13)
     $                      + ang_yj(7,13)
     $                      + ang_yj(8,13)
     $                      + ang_yj(9,13)
                    angj(5) = ang_yj(5,14)
     $                      + ang_yj(6,14)
     $                      + ang_yj(7,14)
     $                      + ang_yj(8,14)
     $                      + ang_yj(9,14)
                    angj(6) = ang_yj(5,15)
     $                      + ang_yj(6,15)
     $                      + ang_yj(7,15)
     $                      + ang_yj(8,15)
     $                      + ang_yj(9,15)
                    angj(7) = ang_yj(5,16)
     $                      + ang_yj(6,16)
     $                      + ang_yj(7,16)
     $                      + ang_yj(8,16)
     $                      + ang_yj(9,16)
c
                    angsav(1,4,3,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(5) +
     $               angi(2)*ylmj_t(6) +
     $               angi(3)*ylmj_t(7) +
     $               angi(4)*ylmj_t(8) +
     $               angi(5)*ylmj_t(9)
c
                    angsav(2,4,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(10,1) +
     $               angj(2)*grylmi_t(11,1) +
     $               angj(3)*grylmi_t(12,1) +
     $               angj(4)*grylmi_t(13,1) +
     $               angj(5)*grylmi_t(14,1) +
     $               angj(6)*grylmi_t(15,1) +
     $               angj(7)*grylmi_t(16,1)
                    angsav(3,4,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(10,2) +
     $               angj(2)*grylmi_t(11,2) +
     $               angj(3)*grylmi_t(12,2) +
     $               angj(4)*grylmi_t(13,2) +
     $               angj(5)*grylmi_t(14,2) +
     $               angj(6)*grylmi_t(15,2) +
     $               angj(7)*grylmi_t(16,2)
                    angsav(4,4,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(10,3) +
     $               angj(2)*grylmi_t(11,3) +
     $               angj(3)*grylmi_t(12,3) +
     $               angj(4)*grylmi_t(13,3) +
     $               angj(5)*grylmi_t(14,3) +
     $               angj(6)*grylmi_t(15,3) +
     $               angj(7)*grylmi_t(16,3)
c
                    angsav(5,4,3,llpot,lmpi,lmpj) =
     $               grylmj_t(5,1)*angi(1) +
     $               grylmj_t(6,1)*angi(2) +
     $               grylmj_t(7,1)*angi(3) +
     $               grylmj_t(8,1)*angi(4) +
     $               grylmj_t(9,1)*angi(5)
                    angsav(6,4,3,llpot,lmpi,lmpj) =
     $               grylmj_t(5,2)*angi(1) +
     $               grylmj_t(6,2)*angi(2) +
     $               grylmj_t(7,2)*angi(3) +
     $               grylmj_t(8,2)*angi(4) +
     $               grylmj_t(9,2)*angi(5)
                    angsav(7,4,3,llpot,lmpi,lmpj) =
     $               grylmj_t(5,3)*angi(1) +
     $               grylmj_t(6,3)*angi(2) +
     $               grylmj_t(7,3)*angi(3) +
     $               grylmj_t(8,3)*angi(4) +
     $               grylmj_t(9,3)*angi(5)
c
                  elseif( lli .eq. 5 )then
c
                    angi(1) = ang_yi(17,5)
     $                      + ang_yi(18,5)
     $                      + ang_yi(19,5)
     $                      + ang_yi(20,5)
     $                      + ang_yi(21,5)
     $                      + ang_yi(22,5)
     $                      + ang_yi(23,5)
     $                      + ang_yi(24,5)
     $                      + ang_yi(25,5)
                    angi(2) = ang_yi(17,6)
     $                      + ang_yi(18,6)
     $                      + ang_yi(19,6)
     $                      + ang_yi(20,6)
     $                      + ang_yi(21,6)
     $                      + ang_yi(22,6)
     $                      + ang_yi(23,6)
     $                      + ang_yi(24,6)
     $                      + ang_yi(25,6)
                    angi(3) = ang_yi(17,7)
     $                      + ang_yi(18,7)
     $                      + ang_yi(19,7)
     $                      + ang_yi(20,7)
     $                      + ang_yi(21,7)
     $                      + ang_yi(22,7)
     $                      + ang_yi(23,7)
     $                      + ang_yi(24,7)
     $                      + ang_yi(25,7)
                    angi(4) = ang_yi(17,8)
     $                      + ang_yi(18,8)
     $                      + ang_yi(19,8)
     $                      + ang_yi(20,8)
     $                      + ang_yi(21,8)
     $                      + ang_yi(22,8)
     $                      + ang_yi(23,8)
     $                      + ang_yi(24,8)
     $                      + ang_yi(25,8)
                    angi(5) = ang_yi(17,9)
     $                      + ang_yi(18,9)
     $                      + ang_yi(19,9)
     $                      + ang_yi(20,9)
     $                      + ang_yi(21,9)
     $                      + ang_yi(22,9)
     $                      + ang_yi(23,9)
     $                      + ang_yi(24,9)
     $                      + ang_yi(25,9)
c
                    angj(1) = ang_yj(5,17)
     $                      + ang_yj(6,17)
     $                      + ang_yj(7,17)
     $                      + ang_yj(8,17)
     $                      + ang_yj(9,17)
                    angj(2) = ang_yj(5,18)
     $                      + ang_yj(6,18)
     $                      + ang_yj(7,18)
     $                      + ang_yj(8,18)
     $                      + ang_yj(9,18)
                    angj(3) = ang_yj(5,19)
     $                      + ang_yj(6,19)
     $                      + ang_yj(7,19)
     $                      + ang_yj(8,19)
     $                      + ang_yj(9,19)
                    angj(4) = ang_yj(5,20)
     $                      + ang_yj(6,20)
     $                      + ang_yj(7,20)
     $                      + ang_yj(8,20)
     $                      + ang_yj(9,20)
                    angj(5) = ang_yj(5,21)
     $                      + ang_yj(6,21)
     $                      + ang_yj(7,21)
     $                      + ang_yj(8,21)
     $                      + ang_yj(9,21)
                    angj(6) = ang_yj(5,22)
     $                      + ang_yj(6,22)
     $                      + ang_yj(7,22)
     $                      + ang_yj(8,22)
     $                      + ang_yj(9,22)
                    angj(7) = ang_yj(5,23)
     $                      + ang_yj(6,23)
     $                      + ang_yj(7,23)
     $                      + ang_yj(8,23)
     $                      + ang_yj(9,23)
                    angj(8) = ang_yj(5,24)
     $                      + ang_yj(6,24)
     $                      + ang_yj(7,24)
     $                      + ang_yj(8,24)
     $                      + ang_yj(9,24)
                    angj(9) = ang_yj(5,25)
     $                      + ang_yj(6,25)
     $                      + ang_yj(7,25)
     $                      + ang_yj(8,25)
     $                      + ang_yj(9,25)
c
                    angsav(1,5,3,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(5) +
     $               angi(2)*ylmj_t(6) +
     $               angi(3)*ylmj_t(7) +
     $               angi(4)*ylmj_t(8) +
     $               angi(5)*ylmj_t(9)
c
                    angsav(2,5,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(17,1) +
     $               angj(2)*grylmi_t(18,1) +
     $               angj(3)*grylmi_t(19,1) +
     $               angj(4)*grylmi_t(20,1) +
     $               angj(5)*grylmi_t(21,1) +
     $               angj(6)*grylmi_t(22,1) +
     $               angj(7)*grylmi_t(23,1) +
     $               angj(8)*grylmi_t(24,1) +
     $               angj(9)*grylmi_t(25,1)
                    angsav(3,5,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(17,2) +
     $               angj(2)*grylmi_t(18,2) +
     $               angj(3)*grylmi_t(19,2) +
     $               angj(4)*grylmi_t(20,2) +
     $               angj(5)*grylmi_t(21,2) +
     $               angj(6)*grylmi_t(22,2) +
     $               angj(7)*grylmi_t(23,2) +
     $               angj(8)*grylmi_t(24,2) +
     $               angj(9)*grylmi_t(25,2)
                    angsav(4,5,3,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(17,3) +
     $               angj(2)*grylmi_t(18,3) +
     $               angj(3)*grylmi_t(19,3) +
     $               angj(4)*grylmi_t(20,3) +
     $               angj(5)*grylmi_t(21,3) +
     $               angj(6)*grylmi_t(22,3) +
     $               angj(7)*grylmi_t(23,3) +
     $               angj(8)*grylmi_t(24,3) +
     $               angj(9)*grylmi_t(25,3)
c
                    angsav(5,5,3,llpot,lmpi,lmpj) =
     $               grylmj_t(5,1)*angi(1) +
     $               grylmj_t(6,1)*angi(2) +
     $               grylmj_t(7,1)*angi(3) +
     $               grylmj_t(8,1)*angi(4) +
     $               grylmj_t(9,1)*angi(5)
                    angsav(6,5,3,llpot,lmpi,lmpj) =
     $               grylmj_t(5,2)*angi(1) +
     $               grylmj_t(6,2)*angi(2) +
     $               grylmj_t(7,2)*angi(3) +
     $               grylmj_t(8,2)*angi(4) +
     $               grylmj_t(9,2)*angi(5)
                    angsav(7,5,3,llpot,lmpi,lmpj) =
     $               grylmj_t(5,3)*angi(1) +
     $               grylmj_t(6,3)*angi(2) +
     $               grylmj_t(7,3)*angi(3) +
     $               grylmj_t(8,3)*angi(4) +
     $               grylmj_t(9,3)*angi(5)
c
                  else
c
                    call STOPXERR( 'nlangf-3' )
c
                  endif
c
                enddo
c
              elseif( llj .eq. 4 )then
c
                do  lli=1,limxp1
c
                  if( lli .eq. 1 )then
c
                    angj(1) = ang_yj(10,1)
     $                      + ang_yj(11,1)
     $                      + ang_yj(12,1)
     $                      + ang_yj(13,1)
     $                      + ang_yj(14,1)
     $                      + ang_yj(15,1)
     $                      + ang_yj(16,1)
c
                    angsav(1,1,4,llpot,lmpi,lmpj) =
     $               ylmi_t(1)*angj(1)
c
                    angsav(2,1,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,1)
                    angsav(3,1,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,2)
                    angsav(4,1,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,3)
c
                    angsav(5,1,4,llpot,lmpi,lmpj) =
     $               ang_yi(1,10)*grylmj_t(10,1) +
     $               ang_yi(1,11)*grylmj_t(11,1) +
     $               ang_yi(1,12)*grylmj_t(12,1) +
     $               ang_yi(1,13)*grylmj_t(13,1) +
     $               ang_yi(1,14)*grylmj_t(14,1) +
     $               ang_yi(1,15)*grylmj_t(15,1) +
     $               ang_yi(1,16)*grylmj_t(16,1)
                    angsav(6,1,4,llpot,lmpi,lmpj) =
     $               ang_yi(1,10)*grylmj_t(10,2) +
     $               ang_yi(1,11)*grylmj_t(11,2) +
     $               ang_yi(1,12)*grylmj_t(12,2) +
     $               ang_yi(1,13)*grylmj_t(13,2) +
     $               ang_yi(1,14)*grylmj_t(14,2) +
     $               ang_yi(1,15)*grylmj_t(15,2) +
     $               ang_yi(1,16)*grylmj_t(16,2)
                    angsav(7,1,4,llpot,lmpi,lmpj) =
     $               ang_yi(1,10)*grylmj_t(10,3) +
     $               ang_yi(1,11)*grylmj_t(11,3) +
     $               ang_yi(1,12)*grylmj_t(12,3) +
     $               ang_yi(1,13)*grylmj_t(13,3) +
     $               ang_yi(1,14)*grylmj_t(14,3) +
     $               ang_yi(1,15)*grylmj_t(15,3) +
     $               ang_yi(1,16)*grylmj_t(16,3)
c
                  elseif( lli .eq. 2 )then
c
                    angi(1) = ang_yi(2,10)
     $                      + ang_yi(3,10)
     $                      + ang_yi(4,10)
                    angi(2) = ang_yi(2,11)
     $                      + ang_yi(3,11)
     $                      + ang_yi(4,11)
                    angi(3) = ang_yi(2,12)
     $                      + ang_yi(3,12)
     $                      + ang_yi(4,12)
                    angi(4) = ang_yi(2,13)
     $                      + ang_yi(3,13)
     $                      + ang_yi(4,13)
                    angi(5) = ang_yi(2,14)
     $                      + ang_yi(3,14)
     $                      + ang_yi(4,14)
                    angi(6) = ang_yi(2,15)
     $                      + ang_yi(3,15)
     $                      + ang_yi(4,15)
                    angi(7) = ang_yi(2,16)
     $                      + ang_yi(3,16)
     $                      + ang_yi(4,16)
c
                    angj(1) = ang_yj(10,2)
     $                      + ang_yj(11,2)
     $                      + ang_yj(12,2)
     $                      + ang_yj(13,2)
     $                      + ang_yj(14,2)
     $                      + ang_yj(15,2)
     $                      + ang_yj(16,2)
                    angj(2) = ang_yj(10,3)
     $                      + ang_yj(11,3)
     $                      + ang_yj(12,3)
     $                      + ang_yj(13,3)
     $                      + ang_yj(14,3)
     $                      + ang_yj(15,3)
     $                      + ang_yj(16,3)
                    angj(3) = ang_yj(10,4)
     $                      + ang_yj(11,4)
     $                      + ang_yj(12,4)
     $                      + ang_yj(13,4)
     $                      + ang_yj(14,4)
     $                      + ang_yj(15,4)
     $                      + ang_yj(16,4)
c
                    angsav(1,2,4,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(10) +
     $               angi(2)*ylmj_t(11) +
     $               angi(3)*ylmj_t(12) +
     $               angi(4)*ylmj_t(13) +
     $               angi(5)*ylmj_t(14) +
     $               angi(6)*ylmj_t(15) +
     $               angi(7)*ylmj_t(16)
c
                    angsav(2,2,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(2,1) +
     $               angj(2)*grylmi_t(3,1) +
     $               angj(3)*grylmi_t(4,1)
                    angsav(3,2,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(2,2) +
     $               angj(2)*grylmi_t(3,2) +
     $               angj(3)*grylmi_t(4,2)
                    angsav(4,2,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(2,3) +
     $               angj(2)*grylmi_t(3,3) +
     $               angj(3)*grylmi_t(4,3)
c
                    angsav(5,2,4,llpot,lmpi,lmpj) =
     $               grylmj_t(10,1)*angi(1) +
     $               grylmj_t(11,1)*angi(2) +
     $               grylmj_t(12,1)*angi(3) +
     $               grylmj_t(13,1)*angi(4) +
     $               grylmj_t(14,1)*angi(5) +
     $               grylmj_t(15,1)*angi(6) +
     $               grylmj_t(16,1)*angi(7)
                    angsav(6,2,4,llpot,lmpi,lmpj) =
     $               grylmj_t(10,2)*angi(1) +
     $               grylmj_t(11,2)*angi(2) +
     $               grylmj_t(12,2)*angi(3) +
     $               grylmj_t(13,2)*angi(4) +
     $               grylmj_t(14,2)*angi(5) +
     $               grylmj_t(15,2)*angi(6) +
     $               grylmj_t(16,2)*angi(7)
                    angsav(7,2,4,llpot,lmpi,lmpj) =
     $               grylmj_t(10,3)*angi(1) +
     $               grylmj_t(11,3)*angi(2) +
     $               grylmj_t(12,3)*angi(3) +
     $               grylmj_t(13,3)*angi(4) +
     $               grylmj_t(14,3)*angi(5) +
     $               grylmj_t(15,3)*angi(6) +
     $               grylmj_t(16,3)*angi(7)
c
                  elseif( lli .eq. 3 )then
c
                    angi(1) = ang_yi(5,10)
     $                      + ang_yi(6,10)
     $                      + ang_yi(7,10)
     $                      + ang_yi(8,10)
     $                      + ang_yi(9,10)
                    angi(2) = ang_yi(5,11)
     $                      + ang_yi(6,11)
     $                      + ang_yi(7,11)
     $                      + ang_yi(8,11)
     $                      + ang_yi(9,11)
                    angi(3) = ang_yi(5,12)
     $                      + ang_yi(6,12)
     $                      + ang_yi(7,12)
     $                      + ang_yi(8,12)
     $                      + ang_yi(9,12)
                    angi(4) = ang_yi(5,13)
     $                      + ang_yi(6,13)
     $                      + ang_yi(7,13)
     $                      + ang_yi(8,13)
     $                      + ang_yi(9,13)
                    angi(5) = ang_yi(5,14)
     $                      + ang_yi(6,14)
     $                      + ang_yi(7,14)
     $                      + ang_yi(8,14)
     $                      + ang_yi(9,14)
                    angi(6) = ang_yi(5,15)
     $                      + ang_yi(6,15)
     $                      + ang_yi(7,15)
     $                      + ang_yi(8,15)
     $                      + ang_yi(9,15)
                    angi(7) = ang_yi(5,16)
     $                      + ang_yi(6,16)
     $                      + ang_yi(7,16)
     $                      + ang_yi(8,16)
     $                      + ang_yi(9,16)
c
                    angj(1) = ang_yj(10,5)
     $                      + ang_yj(11,5)
     $                      + ang_yj(12,5)
     $                      + ang_yj(13,5)
     $                      + ang_yj(14,5)
     $                      + ang_yj(15,5)
     $                      + ang_yj(16,5)
                    angj(2) = ang_yj(10,6)
     $                      + ang_yj(11,6)
     $                      + ang_yj(12,6)
     $                      + ang_yj(13,6)
     $                      + ang_yj(14,6)
     $                      + ang_yj(15,6)
     $                      + ang_yj(16,6)
                    angj(3) = ang_yj(10,7)
     $                      + ang_yj(11,7)
     $                      + ang_yj(12,7)
     $                      + ang_yj(13,7)
     $                      + ang_yj(14,7)
     $                      + ang_yj(15,7)
     $                      + ang_yj(16,7)
                    angj(4) = ang_yj(10,8)
     $                      + ang_yj(11,8)
     $                      + ang_yj(12,8)
     $                      + ang_yj(13,8)
     $                      + ang_yj(14,8)
     $                      + ang_yj(15,8)
     $                      + ang_yj(16,8)
                    angj(5) = ang_yj(10,9)
     $                      + ang_yj(11,9)
     $                      + ang_yj(12,9)
     $                      + ang_yj(13,9)
     $                      + ang_yj(14,9)
     $                      + ang_yj(15,9)
     $                      + ang_yj(16,9)
c
                    angsav(1,3,4,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(10) +
     $               angi(2)*ylmj_t(11) +
     $               angi(3)*ylmj_t(12) +
     $               angi(4)*ylmj_t(13) +
     $               angi(5)*ylmj_t(14) +
     $               angi(6)*ylmj_t(15) +
     $               angi(7)*ylmj_t(16)
c
                    angsav(2,3,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(5,1) +
     $               angj(2)*grylmi_t(6,1) +
     $               angj(3)*grylmi_t(7,1) +
     $               angj(4)*grylmi_t(8,1) +
     $               angj(5)*grylmi_t(9,1)
                    angsav(3,3,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(5,2) +
     $               angj(2)*grylmi_t(6,2) +
     $               angj(3)*grylmi_t(7,2) +
     $               angj(4)*grylmi_t(8,2) +
     $               angj(5)*grylmi_t(9,2)
                    angsav(4,3,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(5,3) +
     $               angj(2)*grylmi_t(6,3) +
     $               angj(3)*grylmi_t(7,3) +
     $               angj(4)*grylmi_t(8,3) +
     $               angj(5)*grylmi_t(9,3)
c
                    angsav(5,3,4,llpot,lmpi,lmpj) =
     $               grylmj_t(10,1)*angi(1) +
     $               grylmj_t(11,1)*angi(2) +
     $               grylmj_t(12,1)*angi(3) +
     $               grylmj_t(13,1)*angi(4) +
     $               grylmj_t(14,1)*angi(5) +
     $               grylmj_t(15,1)*angi(6) +
     $               grylmj_t(16,1)*angi(7)
                    angsav(6,3,4,llpot,lmpi,lmpj) =
     $               grylmj_t(10,2)*angi(1) +
     $               grylmj_t(11,2)*angi(2) +
     $               grylmj_t(12,2)*angi(3) +
     $               grylmj_t(13,2)*angi(4) +
     $               grylmj_t(14,2)*angi(5) +
     $               grylmj_t(15,2)*angi(6) +
     $               grylmj_t(16,2)*angi(7)
                    angsav(7,3,4,llpot,lmpi,lmpj) =
     $               grylmj_t(10,3)*angi(1) +
     $               grylmj_t(11,3)*angi(2) +
     $               grylmj_t(12,3)*angi(3) +
     $               grylmj_t(13,3)*angi(4) +
     $               grylmj_t(14,3)*angi(5) +
     $               grylmj_t(15,3)*angi(6) +
     $               grylmj_t(16,3)*angi(7)
c
                  elseif( lli .eq. 4 )then
c
                    angi(1) = ang_yi(10,10)
     $                      + ang_yi(11,10)
     $                      + ang_yi(12,10)
     $                      + ang_yi(13,10)
     $                      + ang_yi(14,10)
     $                      + ang_yi(15,10)
     $                      + ang_yi(16,10)
                    angi(2) = ang_yi(10,11)
     $                      + ang_yi(11,11)
     $                      + ang_yi(12,11)
     $                      + ang_yi(13,11)
     $                      + ang_yi(14,11)
     $                      + ang_yi(15,11)
     $                      + ang_yi(16,11)
                    angi(3) = ang_yi(10,12)
     $                      + ang_yi(11,12)
     $                      + ang_yi(12,12)
     $                      + ang_yi(13,12)
     $                      + ang_yi(14,12)
     $                      + ang_yi(15,12)
     $                      + ang_yi(16,12)
                    angi(4) = ang_yi(10,13)
     $                      + ang_yi(11,13)
     $                      + ang_yi(12,13)
     $                      + ang_yi(13,13)
     $                      + ang_yi(14,13)
     $                      + ang_yi(15,13)
     $                      + ang_yi(16,13)
                    angi(5) = ang_yi(10,14)
     $                      + ang_yi(11,14)
     $                      + ang_yi(12,14)
     $                      + ang_yi(13,14)
     $                      + ang_yi(14,14)
     $                      + ang_yi(15,14)
     $                      + ang_yi(16,14)
                    angi(6) = ang_yi(10,15)
     $                      + ang_yi(11,15)
     $                      + ang_yi(12,15)
     $                      + ang_yi(13,15)
     $                      + ang_yi(14,15)
     $                      + ang_yi(15,15)
     $                      + ang_yi(16,15)
                    angi(7) = ang_yi(10,16)
     $                      + ang_yi(11,16)
     $                      + ang_yi(12,16)
     $                      + ang_yi(13,16)
     $                      + ang_yi(14,16)
     $                      + ang_yi(15,16)
     $                      + ang_yi(16,16)
c
                    angj(1) = ang_yj(10,10)
     $                      + ang_yj(11,10)
     $                      + ang_yj(12,10)
     $                      + ang_yj(13,10)
     $                      + ang_yj(14,10)
     $                      + ang_yj(15,10)
     $                      + ang_yj(16,10)
                    angj(2) = ang_yj(10,11)
     $                      + ang_yj(11,11)
     $                      + ang_yj(12,11)
     $                      + ang_yj(13,11)
     $                      + ang_yj(14,11)
     $                      + ang_yj(15,11)
     $                      + ang_yj(16,11)
                    angj(3) = ang_yj(10,12)
     $                      + ang_yj(11,12)
     $                      + ang_yj(12,12)
     $                      + ang_yj(13,12)
     $                      + ang_yj(14,12)
     $                      + ang_yj(15,12)
     $                      + ang_yj(16,12)
                    angj(4) = ang_yj(10,13)
     $                      + ang_yj(11,13)
     $                      + ang_yj(12,13)
     $                      + ang_yj(13,13)
     $                      + ang_yj(14,13)
     $                      + ang_yj(15,13)
     $                      + ang_yj(16,13)
                    angj(5) = ang_yj(10,14)
     $                      + ang_yj(11,14)
     $                      + ang_yj(12,14)
     $                      + ang_yj(13,14)
     $                      + ang_yj(14,14)
     $                      + ang_yj(15,14)
     $                      + ang_yj(16,14)
                    angj(6) = ang_yj(10,15)
     $                      + ang_yj(11,15)
     $                      + ang_yj(12,15)
     $                      + ang_yj(13,15)
     $                      + ang_yj(14,15)
     $                      + ang_yj(15,15)
     $                      + ang_yj(16,15)
                    angj(7) = ang_yj(10,16)
     $                      + ang_yj(11,16)
     $                      + ang_yj(12,16)
     $                      + ang_yj(13,16)
     $                      + ang_yj(14,16)
     $                      + ang_yj(15,16)
     $                      + ang_yj(16,16)
c
                    angsav(1,4,4,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(10) +
     $               angi(2)*ylmj_t(11) +
     $               angi(3)*ylmj_t(12) +
     $               angi(4)*ylmj_t(13) +
     $               angi(5)*ylmj_t(14) +
     $               angi(6)*ylmj_t(15) +
     $               angi(7)*ylmj_t(16)
c
                    angsav(2,4,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(10,1) +
     $               angj(2)*grylmi_t(11,1) +
     $               angj(3)*grylmi_t(12,1) +
     $               angj(4)*grylmi_t(13,1) +
     $               angj(5)*grylmi_t(14,1) +
     $               angj(6)*grylmi_t(15,1) +
     $               angj(7)*grylmi_t(16,1)
                    angsav(3,4,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(10,2) +
     $               angj(2)*grylmi_t(11,2) +
     $               angj(3)*grylmi_t(12,2) +
     $               angj(4)*grylmi_t(13,2) +
     $               angj(5)*grylmi_t(14,2) +
     $               angj(6)*grylmi_t(15,2) +
     $               angj(7)*grylmi_t(16,2)
                    angsav(4,4,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(10,3) +
     $               angj(2)*grylmi_t(11,3) +
     $               angj(3)*grylmi_t(12,3) +
     $               angj(4)*grylmi_t(13,3) +
     $               angj(5)*grylmi_t(14,3) +
     $               angj(6)*grylmi_t(15,3) +
     $               angj(7)*grylmi_t(16,3)
c
                    angsav(5,4,4,llpot,lmpi,lmpj) =
     $               grylmj_t(10,1)*angi(1) +
     $               grylmj_t(11,1)*angi(2) +
     $               grylmj_t(12,1)*angi(3) +
     $               grylmj_t(13,1)*angi(4) +
     $               grylmj_t(14,1)*angi(5) +
     $               grylmj_t(15,1)*angi(6) +
     $               grylmj_t(16,1)*angi(7)
                    angsav(6,4,4,llpot,lmpi,lmpj) =
     $               grylmj_t(10,2)*angi(1) +
     $               grylmj_t(11,2)*angi(2) +
     $               grylmj_t(12,2)*angi(3) +
     $               grylmj_t(13,2)*angi(4) +
     $               grylmj_t(14,2)*angi(5) +
     $               grylmj_t(15,2)*angi(6) +
     $               grylmj_t(16,2)*angi(7)
                    angsav(7,4,4,llpot,lmpi,lmpj) =
     $               grylmj_t(10,3)*angi(1) +
     $               grylmj_t(11,3)*angi(2) +
     $               grylmj_t(12,3)*angi(3) +
     $               grylmj_t(13,3)*angi(4) +
     $               grylmj_t(14,3)*angi(5) +
     $               grylmj_t(15,3)*angi(6) +
     $               grylmj_t(16,3)*angi(7)
c
                  elseif( lli .eq. 5 )then
c
                    angi(1) = ang_yi(17,10)
     $                      + ang_yi(18,10)
     $                      + ang_yi(19,10)
     $                      + ang_yi(20,10)
     $                      + ang_yi(21,10)
     $                      + ang_yi(22,10)
     $                      + ang_yi(23,10)
     $                      + ang_yi(24,10)
     $                      + ang_yi(25,10)
                    angi(2) = ang_yi(17,11)
     $                      + ang_yi(18,11)
     $                      + ang_yi(19,11)
     $                      + ang_yi(20,11)
     $                      + ang_yi(21,11)
     $                      + ang_yi(22,11)
     $                      + ang_yi(23,11)
     $                      + ang_yi(24,11)
     $                      + ang_yi(25,11)
                    angi(3) = ang_yi(17,12)
     $                      + ang_yi(18,12)
     $                      + ang_yi(19,12)
     $                      + ang_yi(20,12)
     $                      + ang_yi(21,12)
     $                      + ang_yi(22,12)
     $                      + ang_yi(23,12)
     $                      + ang_yi(24,12)
     $                      + ang_yi(25,12)
                    angi(4) = ang_yi(17,13)
     $                      + ang_yi(18,13)
     $                      + ang_yi(19,13)
     $                      + ang_yi(20,13)
     $                      + ang_yi(21,13)
     $                      + ang_yi(22,13)
     $                      + ang_yi(23,13)
     $                      + ang_yi(24,13)
     $                      + ang_yi(25,13)
                    angi(5) = ang_yi(17,14)
     $                      + ang_yi(18,14)
     $                      + ang_yi(19,14)
     $                      + ang_yi(20,14)
     $                      + ang_yi(21,14)
     $                      + ang_yi(22,14)
     $                      + ang_yi(23,14)
     $                      + ang_yi(24,14)
     $                      + ang_yi(25,14)
                    angi(6) = ang_yi(17,15)
     $                      + ang_yi(18,15)
     $                      + ang_yi(19,15)
     $                      + ang_yi(20,15)
     $                      + ang_yi(21,15)
     $                      + ang_yi(22,15)
     $                      + ang_yi(23,15)
     $                      + ang_yi(24,15)
     $                      + ang_yi(25,15)
                    angi(7) = ang_yi(17,16)
     $                      + ang_yi(18,16)
     $                      + ang_yi(19,16)
     $                      + ang_yi(20,16)
     $                      + ang_yi(21,16)
     $                      + ang_yi(22,16)
     $                      + ang_yi(23,16)
     $                      + ang_yi(24,16)
     $                      + ang_yi(25,16)
c
                    angj(1) = ang_yj(10,17)
     $                      + ang_yj(11,17)
     $                      + ang_yj(12,17)
     $                      + ang_yj(13,17)
     $                      + ang_yj(14,17)
     $                      + ang_yj(15,17)
     $                      + ang_yj(16,17)
                    angj(2) = ang_yj(10,18)
     $                      + ang_yj(11,18)
     $                      + ang_yj(12,18)
     $                      + ang_yj(13,18)
     $                      + ang_yj(14,18)
     $                      + ang_yj(15,18)
     $                      + ang_yj(16,18)
                    angj(3) = ang_yj(10,19)
     $                      + ang_yj(11,19)
     $                      + ang_yj(12,19)
     $                      + ang_yj(13,19)
     $                      + ang_yj(14,19)
     $                      + ang_yj(15,19)
     $                      + ang_yj(16,19)
                    angj(4) = ang_yj(10,20)
     $                      + ang_yj(11,20)
     $                      + ang_yj(12,20)
     $                      + ang_yj(13,20)
     $                      + ang_yj(14,20)
     $                      + ang_yj(15,20)
     $                      + ang_yj(16,20)
                    angj(5) = ang_yj(10,21)
     $                      + ang_yj(11,21)
     $                      + ang_yj(12,21)
     $                      + ang_yj(13,21)
     $                      + ang_yj(14,21)
     $                      + ang_yj(15,21)
     $                      + ang_yj(16,21)
                    angj(6) = ang_yj(10,22)
     $                      + ang_yj(11,22)
     $                      + ang_yj(12,22)
     $                      + ang_yj(13,22)
     $                      + ang_yj(14,22)
     $                      + ang_yj(15,22)
     $                      + ang_yj(16,22)
                    angj(7) = ang_yj(10,23)
     $                      + ang_yj(11,23)
     $                      + ang_yj(12,23)
     $                      + ang_yj(13,23)
     $                      + ang_yj(14,23)
     $                      + ang_yj(15,23)
     $                      + ang_yj(16,23)
                    angj(8) = ang_yj(10,24)
     $                      + ang_yj(11,24)
     $                      + ang_yj(12,24)
     $                      + ang_yj(13,24)
     $                      + ang_yj(14,24)
     $                      + ang_yj(15,24)
     $                      + ang_yj(16,24)
                    angj(9) = ang_yj(10,25)
     $                      + ang_yj(11,25)
     $                      + ang_yj(12,25)
     $                      + ang_yj(13,25)
     $                      + ang_yj(14,25)
     $                      + ang_yj(15,25)
     $                      + ang_yj(16,25)
c
                    angsav(1,5,4,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(10) +
     $               angi(2)*ylmj_t(11) +
     $               angi(3)*ylmj_t(12) +
     $               angi(4)*ylmj_t(13) +
     $               angi(5)*ylmj_t(14) +
     $               angi(6)*ylmj_t(15) +
     $               angi(7)*ylmj_t(16)
c
                    angsav(2,5,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(17,1) +
     $               angj(2)*grylmi_t(18,1) +
     $               angj(3)*grylmi_t(19,1) +
     $               angj(4)*grylmi_t(20,1) +
     $               angj(5)*grylmi_t(21,1) +
     $               angj(6)*grylmi_t(22,1) +
     $               angj(7)*grylmi_t(23,1) +
     $               angj(8)*grylmi_t(24,1) +
     $               angj(9)*grylmi_t(25,1)
                    angsav(3,5,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(17,2) +
     $               angj(2)*grylmi_t(18,2) +
     $               angj(3)*grylmi_t(19,2) +
     $               angj(4)*grylmi_t(20,2) +
     $               angj(5)*grylmi_t(21,2) +
     $               angj(6)*grylmi_t(22,2) +
     $               angj(7)*grylmi_t(23,2) +
     $               angj(8)*grylmi_t(24,2) +
     $               angj(9)*grylmi_t(25,2)
                    angsav(4,5,4,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(17,3) +
     $               angj(2)*grylmi_t(18,3) +
     $               angj(3)*grylmi_t(19,3) +
     $               angj(4)*grylmi_t(20,3) +
     $               angj(5)*grylmi_t(21,3) +
     $               angj(6)*grylmi_t(22,3) +
     $               angj(7)*grylmi_t(23,3) +
     $               angj(8)*grylmi_t(24,3) +
     $               angj(9)*grylmi_t(25,3)
c
                    angsav(5,5,4,llpot,lmpi,lmpj) =
     $               grylmj_t(10,1)*angi(1) +
     $               grylmj_t(11,1)*angi(2) +
     $               grylmj_t(12,1)*angi(3) +
     $               grylmj_t(13,1)*angi(4) +
     $               grylmj_t(14,1)*angi(5) +
     $               grylmj_t(15,1)*angi(6) +
     $               grylmj_t(16,1)*angi(7)
                    angsav(6,5,4,llpot,lmpi,lmpj) =
     $               grylmj_t(10,2)*angi(1) +
     $               grylmj_t(11,2)*angi(2) +
     $               grylmj_t(12,2)*angi(3) +
     $               grylmj_t(13,2)*angi(4) +
     $               grylmj_t(14,2)*angi(5) +
     $               grylmj_t(15,2)*angi(6) +
     $               grylmj_t(16,2)*angi(7)
                    angsav(7,5,4,llpot,lmpi,lmpj) =
     $               grylmj_t(10,3)*angi(1) +
     $               grylmj_t(11,3)*angi(2) +
     $               grylmj_t(12,3)*angi(3) +
     $               grylmj_t(13,3)*angi(4) +
     $               grylmj_t(14,3)*angi(5) +
     $               grylmj_t(15,3)*angi(6) +
     $               grylmj_t(16,3)*angi(7)
c
                  else
c
                    call STOPXERR( 'nlangf-4' )
c
                  endif
c
                enddo
c
              elseif( llj .eq. 5 )then
c
                do  lli=1,limxp1
c
                  if( lli .eq. 1 )then
c
                    angj(1) = ang_yj(17,1)
     $                      + ang_yj(18,1)
     $                      + ang_yj(19,1)
     $                      + ang_yj(20,1)
     $                      + ang_yj(21,1)
     $                      + ang_yj(22,1)
     $                      + ang_yj(23,1)
     $                      + ang_yj(24,1)
     $                      + ang_yj(25,1)
c
                    angsav(1,1,5,llpot,lmpi,lmpj) =
     $               ylmi_t(1)*angj(1)
c
                    angsav(2,1,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,1)
                    angsav(3,1,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,2)
                    angsav(4,1,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(1,3)
c
                    angsav(5,1,5,llpot,lmpi,lmpj) =
     $               ang_yi(1,17)*grylmj_t(17,1) +
     $               ang_yi(1,18)*grylmj_t(18,1) +
     $               ang_yi(1,19)*grylmj_t(19,1) +
     $               ang_yi(1,20)*grylmj_t(20,1) +
     $               ang_yi(1,21)*grylmj_t(21,1) +
     $               ang_yi(1,22)*grylmj_t(22,1) +
     $               ang_yi(1,23)*grylmj_t(23,1) +
     $               ang_yi(1,24)*grylmj_t(24,1) +
     $               ang_yi(1,25)*grylmj_t(25,1)
                    angsav(6,1,5,llpot,lmpi,lmpj) =
     $               ang_yi(1,17)*grylmj_t(17,2) +
     $               ang_yi(1,18)*grylmj_t(18,2) +
     $               ang_yi(1,19)*grylmj_t(19,2) +
     $               ang_yi(1,20)*grylmj_t(20,2) +
     $               ang_yi(1,21)*grylmj_t(21,2) +
     $               ang_yi(1,22)*grylmj_t(22,2) +
     $               ang_yi(1,23)*grylmj_t(23,2) +
     $               ang_yi(1,24)*grylmj_t(24,2) +
     $               ang_yi(1,25)*grylmj_t(25,2)
                    angsav(7,1,5,llpot,lmpi,lmpj) =
     $               ang_yi(1,17)*grylmj_t(17,3) +
     $               ang_yi(1,18)*grylmj_t(18,3) +
     $               ang_yi(1,19)*grylmj_t(19,3) +
     $               ang_yi(1,20)*grylmj_t(20,3) +
     $               ang_yi(1,21)*grylmj_t(21,3) +
     $               ang_yi(1,22)*grylmj_t(22,3) +
     $               ang_yi(1,23)*grylmj_t(23,3) +
     $               ang_yi(1,24)*grylmj_t(24,3) +
     $               ang_yi(1,25)*grylmj_t(25,3)
c
                  elseif( lli .eq. 2 )then
c
                    angi(1) = ang_yi(2,17)
     $                      + ang_yi(3,17)
     $                      + ang_yi(4,17)
                    angi(2) = ang_yi(2,18)
     $                      + ang_yi(3,18)
     $                      + ang_yi(4,18)
                    angi(3) = ang_yi(2,19)
     $                      + ang_yi(3,19)
     $                      + ang_yi(4,19)
                    angi(4) = ang_yi(2,20)
     $                      + ang_yi(3,20)
     $                      + ang_yi(4,20)
                    angi(5) = ang_yi(2,21)
     $                      + ang_yi(3,21)
     $                      + ang_yi(4,21)
                    angi(6) = ang_yi(2,22)
     $                      + ang_yi(3,22)
     $                      + ang_yi(4,22)
                    angi(7) = ang_yi(2,23)
     $                      + ang_yi(3,23)
     $                      + ang_yi(4,23)
                    angi(8) = ang_yi(2,24)
     $                      + ang_yi(3,24)
     $                      + ang_yi(4,24)
                    angi(9) = ang_yi(2,25)
     $                      + ang_yi(3,25)
     $                      + ang_yi(4,25)
c
                    angj(1) = ang_yj(17,2)
     $                      + ang_yj(18,2)
     $                      + ang_yj(19,2)
     $                      + ang_yj(20,2)
     $                      + ang_yj(21,2)
     $                      + ang_yj(22,2)
     $                      + ang_yj(23,2)
     $                      + ang_yj(24,2)
     $                      + ang_yj(25,2)
                    angj(2) = ang_yj(17,3)
     $                      + ang_yj(18,3)
     $                      + ang_yj(19,3)
     $                      + ang_yj(20,3)
     $                      + ang_yj(21,3)
     $                      + ang_yj(22,3)
     $                      + ang_yj(23,3)
     $                      + ang_yj(24,3)
     $                      + ang_yj(25,3)
                    angj(3) = ang_yj(17,4)
     $                      + ang_yj(18,4)
     $                      + ang_yj(19,4)
     $                      + ang_yj(20,4)
     $                      + ang_yj(21,4)
     $                      + ang_yj(22,4)
     $                      + ang_yj(23,4)
     $                      + ang_yj(24,4)
     $                      + ang_yj(25,4)
c
                    angsav(1,2,5,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(17) +
     $               angi(2)*ylmj_t(18) +
     $               angi(3)*ylmj_t(19) +
     $               angi(4)*ylmj_t(20) +
     $               angi(5)*ylmj_t(21) +
     $               angi(6)*ylmj_t(22) +
     $               angi(7)*ylmj_t(23) +
     $               angi(8)*ylmj_t(24) +
     $               angi(9)*ylmj_t(25)
c
                    angsav(2,2,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(2,1) +
     $               angj(2)*grylmi_t(3,1) +
     $               angj(3)*grylmi_t(4,1)
                    angsav(3,2,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(2,2) +
     $               angj(2)*grylmi_t(3,2) +
     $               angj(3)*grylmi_t(4,2)
                    angsav(4,2,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(2,3) +
     $               angj(2)*grylmi_t(3,3) +
     $               angj(3)*grylmi_t(4,3)
c
                    angsav(5,2,5,llpot,lmpi,lmpj) =
     $               grylmj_t(17,1)*angi(1) +
     $               grylmj_t(18,1)*angi(2) +
     $               grylmj_t(19,1)*angi(3) +
     $               grylmj_t(20,1)*angi(4) +
     $               grylmj_t(21,1)*angi(5) +
     $               grylmj_t(22,1)*angi(6) +
     $               grylmj_t(23,1)*angi(7) +
     $               grylmj_t(24,1)*angi(8) +
     $               grylmj_t(25,1)*angi(9)
                    angsav(6,2,5,llpot,lmpi,lmpj) =
     $               grylmj_t(17,2)*angi(1) +
     $               grylmj_t(18,2)*angi(2) +
     $               grylmj_t(19,2)*angi(3) +
     $               grylmj_t(20,2)*angi(4) +
     $               grylmj_t(21,2)*angi(5) +
     $               grylmj_t(22,2)*angi(6) +
     $               grylmj_t(23,2)*angi(7) +
     $               grylmj_t(24,2)*angi(8) +
     $               grylmj_t(25,2)*angi(9)
                    angsav(7,2,5,llpot,lmpi,lmpj) =
     $               grylmj_t(17,3)*angi(1) +
     $               grylmj_t(18,3)*angi(2) +
     $               grylmj_t(19,3)*angi(3) +
     $               grylmj_t(20,3)*angi(4) +
     $               grylmj_t(21,3)*angi(5) +
     $               grylmj_t(22,3)*angi(6) +
     $               grylmj_t(23,3)*angi(7) +
     $               grylmj_t(24,3)*angi(8) +
     $               grylmj_t(25,3)*angi(9)
c
                  elseif( lli .eq. 3 )then
c
                    angi(1) = ang_yi(5,17)
     $                      + ang_yi(6,17)
     $                      + ang_yi(7,17)
     $                      + ang_yi(8,17)
     $                      + ang_yi(9,17)
                    angi(2) = ang_yi(5,18)
     $                      + ang_yi(6,18)
     $                      + ang_yi(7,18)
     $                      + ang_yi(8,18)
     $                      + ang_yi(9,18)
                    angi(3) = ang_yi(5,19)
     $                      + ang_yi(6,19)
     $                      + ang_yi(7,19)
     $                      + ang_yi(8,19)
     $                      + ang_yi(9,19)
                    angi(4) = ang_yi(5,20)
     $                      + ang_yi(6,20)
     $                      + ang_yi(7,20)
     $                      + ang_yi(8,20)
     $                      + ang_yi(9,20)
                    angi(5) = ang_yi(5,21)
     $                      + ang_yi(6,21)
     $                      + ang_yi(7,21)
     $                      + ang_yi(8,21)
     $                      + ang_yi(9,21)
                    angi(6) = ang_yi(5,22)
     $                      + ang_yi(6,22)
     $                      + ang_yi(7,22)
     $                      + ang_yi(8,22)
     $                      + ang_yi(9,22)
                    angi(7) = ang_yi(5,23)
     $                      + ang_yi(6,23)
     $                      + ang_yi(7,23)
     $                      + ang_yi(8,23)
     $                      + ang_yi(9,23)
                    angi(8) = ang_yi(5,24)
     $                      + ang_yi(6,24)
     $                      + ang_yi(7,24)
     $                      + ang_yi(8,24)
     $                      + ang_yi(9,24)
                    angi(9) = ang_yi(5,25)
     $                      + ang_yi(6,25)
     $                      + ang_yi(7,25)
     $                      + ang_yi(8,25)
     $                      + ang_yi(9,25)
c
                    angj(1) = ang_yj(17,5)
     $                      + ang_yj(18,5)
     $                      + ang_yj(19,5)
     $                      + ang_yj(20,5)
     $                      + ang_yj(21,5)
     $                      + ang_yj(22,5)
     $                      + ang_yj(23,5)
     $                      + ang_yj(24,5)
     $                      + ang_yj(25,5)
                    angj(2) = ang_yj(17,6)
     $                      + ang_yj(18,6)
     $                      + ang_yj(19,6)
     $                      + ang_yj(20,6)
     $                      + ang_yj(21,6)
     $                      + ang_yj(22,6)
     $                      + ang_yj(23,6)
     $                      + ang_yj(24,6)
     $                      + ang_yj(25,6)
                    angj(3) = ang_yj(17,7)
     $                      + ang_yj(18,7)
     $                      + ang_yj(19,7)
     $                      + ang_yj(20,7)
     $                      + ang_yj(21,7)
     $                      + ang_yj(22,7)
     $                      + ang_yj(23,7)
     $                      + ang_yj(24,7)
     $                      + ang_yj(25,7)
                    angj(4) = ang_yj(17,8)
     $                      + ang_yj(18,8)
     $                      + ang_yj(19,8)
     $                      + ang_yj(20,8)
     $                      + ang_yj(21,8)
     $                      + ang_yj(22,8)
     $                      + ang_yj(23,8)
     $                      + ang_yj(24,8)
     $                      + ang_yj(25,8)
                    angj(5) = ang_yj(17,9)
     $                      + ang_yj(18,9)
     $                      + ang_yj(19,9)
     $                      + ang_yj(20,9)
     $                      + ang_yj(21,9)
     $                      + ang_yj(22,9)
     $                      + ang_yj(23,9)
     $                      + ang_yj(24,9)
     $                      + ang_yj(25,9)
c
                    angsav(1,3,5,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(17) +
     $               angi(2)*ylmj_t(18) +
     $               angi(3)*ylmj_t(19) +
     $               angi(4)*ylmj_t(20) +
     $               angi(5)*ylmj_t(21) +
     $               angi(6)*ylmj_t(22) +
     $               angi(7)*ylmj_t(23) +
     $               angi(8)*ylmj_t(24) +
     $               angi(9)*ylmj_t(25)
c
                    angsav(2,3,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(5,1) +
     $               angj(2)*grylmi_t(6,1) +
     $               angj(3)*grylmi_t(7,1) +
     $               angj(4)*grylmi_t(8,1) +
     $               angj(5)*grylmi_t(9,1)
                    angsav(3,3,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(5,2) +
     $               angj(2)*grylmi_t(6,2) +
     $               angj(3)*grylmi_t(7,2) +
     $               angj(4)*grylmi_t(8,2) +
     $               angj(5)*grylmi_t(9,2)
                    angsav(4,3,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(5,3) +
     $               angj(2)*grylmi_t(6,3) +
     $               angj(3)*grylmi_t(7,3) +
     $               angj(4)*grylmi_t(8,3) +
     $               angj(5)*grylmi_t(9,3)
c
                    angsav(5,3,5,llpot,lmpi,lmpj) =
     $               grylmj_t(17,1)*angi(1) +
     $               grylmj_t(18,1)*angi(2) +
     $               grylmj_t(19,1)*angi(3) +
     $               grylmj_t(20,1)*angi(4) +
     $               grylmj_t(21,1)*angi(5) +
     $               grylmj_t(22,1)*angi(6) +
     $               grylmj_t(23,1)*angi(7) +
     $               grylmj_t(24,1)*angi(8) +
     $               grylmj_t(25,1)*angi(9)
                    angsav(6,3,5,llpot,lmpi,lmpj) =
     $               grylmj_t(17,2)*angi(1) +
     $               grylmj_t(18,2)*angi(2) +
     $               grylmj_t(19,2)*angi(3) +
     $               grylmj_t(20,2)*angi(4) +
     $               grylmj_t(21,2)*angi(5) +
     $               grylmj_t(22,2)*angi(6) +
     $               grylmj_t(23,2)*angi(7) +
     $               grylmj_t(24,2)*angi(8) +
     $               grylmj_t(25,2)*angi(9)
                    angsav(7,3,5,llpot,lmpi,lmpj) =
     $               grylmj_t(17,3)*angi(1) +
     $               grylmj_t(18,3)*angi(2) +
     $               grylmj_t(19,3)*angi(3) +
     $               grylmj_t(20,3)*angi(4) +
     $               grylmj_t(21,3)*angi(5) +
     $               grylmj_t(22,3)*angi(6) +
     $               grylmj_t(23,3)*angi(7) +
     $               grylmj_t(24,3)*angi(8) +
     $               grylmj_t(25,3)*angi(9)
c
                  elseif( lli .eq. 4 )then
c
                    angi(1) = ang_yi(10,17)
     $                      + ang_yi(11,17)
     $                      + ang_yi(12,17)
     $                      + ang_yi(13,17)
     $                      + ang_yi(14,17)
     $                      + ang_yi(15,17)
     $                      + ang_yi(16,17)
                    angi(2) = ang_yi(10,18)
     $                      + ang_yi(11,18)
     $                      + ang_yi(12,18)
     $                      + ang_yi(13,18)
     $                      + ang_yi(14,18)
     $                      + ang_yi(15,18)
     $                      + ang_yi(16,18)
                    angi(3) = ang_yi(10,19)
     $                      + ang_yi(11,19)
     $                      + ang_yi(12,19)
     $                      + ang_yi(13,19)
     $                      + ang_yi(14,19)
     $                      + ang_yi(15,19)
     $                      + ang_yi(16,19)
                    angi(4) = ang_yi(10,20)
     $                      + ang_yi(11,20)
     $                      + ang_yi(12,20)
     $                      + ang_yi(13,20)
     $                      + ang_yi(14,20)
     $                      + ang_yi(15,20)
     $                      + ang_yi(16,20)
                    angi(5) = ang_yi(10,21)
     $                      + ang_yi(11,21)
     $                      + ang_yi(12,21)
     $                      + ang_yi(13,21)
     $                      + ang_yi(14,21)
     $                      + ang_yi(15,21)
     $                      + ang_yi(16,21)
                    angi(6) = ang_yi(10,22)
     $                      + ang_yi(11,22)
     $                      + ang_yi(12,22)
     $                      + ang_yi(13,22)
     $                      + ang_yi(14,22)
     $                      + ang_yi(15,22)
     $                      + ang_yi(16,22)
                    angi(7) = ang_yi(10,23)
     $                      + ang_yi(11,23)
     $                      + ang_yi(12,23)
     $                      + ang_yi(13,23)
     $                      + ang_yi(14,23)
     $                      + ang_yi(15,23)
     $                      + ang_yi(16,23)
                    angi(8) = ang_yi(10,24)
     $                      + ang_yi(11,24)
     $                      + ang_yi(12,24)
     $                      + ang_yi(13,24)
     $                      + ang_yi(14,24)
     $                      + ang_yi(15,24)
     $                      + ang_yi(16,24)
                    angi(9) = ang_yi(10,25)
     $                      + ang_yi(11,25)
     $                      + ang_yi(12,25)
     $                      + ang_yi(13,25)
     $                      + ang_yi(14,25)
     $                      + ang_yi(15,25)
     $                      + ang_yi(16,25)
c
                    angj(1) = ang_yj(17,10)
     $                      + ang_yj(18,10)
     $                      + ang_yj(19,10)
     $                      + ang_yj(20,10)
     $                      + ang_yj(21,10)
     $                      + ang_yj(22,10)
     $                      + ang_yj(23,10)
     $                      + ang_yj(24,10)
     $                      + ang_yj(25,10)
                    angj(2) = ang_yj(17,11)
     $                      + ang_yj(18,11)
     $                      + ang_yj(19,11)
     $                      + ang_yj(20,11)
     $                      + ang_yj(21,11)
     $                      + ang_yj(22,11)
     $                      + ang_yj(23,11)
     $                      + ang_yj(24,11)
     $                      + ang_yj(25,11)
                    angj(3) = ang_yj(17,12)
     $                      + ang_yj(18,12)
     $                      + ang_yj(19,12)
     $                      + ang_yj(20,12)
     $                      + ang_yj(21,12)
     $                      + ang_yj(22,12)
     $                      + ang_yj(23,12)
     $                      + ang_yj(24,12)
     $                      + ang_yj(25,12)
                    angj(4) = ang_yj(17,13)
     $                      + ang_yj(18,13)
     $                      + ang_yj(19,13)
     $                      + ang_yj(20,13)
     $                      + ang_yj(21,13)
     $                      + ang_yj(22,13)
     $                      + ang_yj(23,13)
     $                      + ang_yj(24,13)
     $                      + ang_yj(25,13)
                    angj(5) = ang_yj(17,14)
     $                      + ang_yj(18,14)
     $                      + ang_yj(19,14)
     $                      + ang_yj(20,14)
     $                      + ang_yj(21,14)
     $                      + ang_yj(22,14)
     $                      + ang_yj(23,14)
     $                      + ang_yj(24,14)
     $                      + ang_yj(25,14)
                    angj(6) = ang_yj(17,15)
     $                      + ang_yj(18,15)
     $                      + ang_yj(19,15)
     $                      + ang_yj(20,15)
     $                      + ang_yj(21,15)
     $                      + ang_yj(22,15)
     $                      + ang_yj(23,15)
     $                      + ang_yj(24,15)
     $                      + ang_yj(25,15)
                    angj(7) = ang_yj(17,16)
     $                      + ang_yj(18,16)
     $                      + ang_yj(19,16)
     $                      + ang_yj(20,16)
     $                      + ang_yj(21,16)
     $                      + ang_yj(22,16)
     $                      + ang_yj(23,16)
     $                      + ang_yj(24,16)
     $                      + ang_yj(25,16)
c
                    angsav(1,4,5,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(17) +
     $               angi(2)*ylmj_t(18) +
     $               angi(3)*ylmj_t(19) +
     $               angi(4)*ylmj_t(20) +
     $               angi(5)*ylmj_t(21) +
     $               angi(6)*ylmj_t(22) +
     $               angi(7)*ylmj_t(23) +
     $               angi(8)*ylmj_t(24) +
     $               angi(9)*ylmj_t(25)
c
                    angsav(2,4,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(10,1) +
     $               angj(2)*grylmi_t(11,1) +
     $               angj(3)*grylmi_t(12,1) +
     $               angj(4)*grylmi_t(13,1) +
     $               angj(5)*grylmi_t(14,1) +
     $               angj(6)*grylmi_t(15,1) +
     $               angj(7)*grylmi_t(16,1)
                    angsav(3,4,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(10,2) +
     $               angj(2)*grylmi_t(11,2) +
     $               angj(3)*grylmi_t(12,2) +
     $               angj(4)*grylmi_t(13,2) +
     $               angj(5)*grylmi_t(14,2) +
     $               angj(6)*grylmi_t(15,2) +
     $               angj(7)*grylmi_t(16,2)
                    angsav(4,4,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(10,3) +
     $               angj(2)*grylmi_t(11,3) +
     $               angj(3)*grylmi_t(12,3) +
     $               angj(4)*grylmi_t(13,3) +
     $               angj(5)*grylmi_t(14,3) +
     $               angj(6)*grylmi_t(15,3) +
     $               angj(7)*grylmi_t(16,3)
c
                    angsav(5,4,5,llpot,lmpi,lmpj) =
     $               grylmj_t(17,1)*angi(1) +
     $               grylmj_t(18,1)*angi(2) +
     $               grylmj_t(19,1)*angi(3) +
     $               grylmj_t(20,1)*angi(4) +
     $               grylmj_t(21,1)*angi(5) +
     $               grylmj_t(22,1)*angi(6) +
     $               grylmj_t(23,1)*angi(7) +
     $               grylmj_t(24,1)*angi(8) +
     $               grylmj_t(25,1)*angi(9)
                    angsav(6,4,5,llpot,lmpi,lmpj) =
     $               grylmj_t(17,2)*angi(1) +
     $               grylmj_t(18,2)*angi(2) +
     $               grylmj_t(19,2)*angi(3) +
     $               grylmj_t(20,2)*angi(4) +
     $               grylmj_t(21,2)*angi(5) +
     $               grylmj_t(22,2)*angi(6) +
     $               grylmj_t(23,2)*angi(7) +
     $               grylmj_t(24,2)*angi(8) +
     $               grylmj_t(25,2)*angi(9)
                    angsav(7,4,5,llpot,lmpi,lmpj) =
     $               grylmj_t(17,3)*angi(1) +
     $               grylmj_t(18,3)*angi(2) +
     $               grylmj_t(19,3)*angi(3) +
     $               grylmj_t(20,3)*angi(4) +
     $               grylmj_t(21,3)*angi(5) +
     $               grylmj_t(22,3)*angi(6) +
     $               grylmj_t(23,3)*angi(7) +
     $               grylmj_t(24,3)*angi(8) +
     $               grylmj_t(25,3)*angi(9)
c
                  elseif( lli .eq. 5 )then
c
                    angi(1) = ang_yi(17,17)
     $                      + ang_yi(18,17)
     $                      + ang_yi(19,17)
     $                      + ang_yi(20,17)
     $                      + ang_yi(21,17)
     $                      + ang_yi(22,17)
     $                      + ang_yi(23,17)
     $                      + ang_yi(24,17)
     $                      + ang_yi(25,17)
                    angi(2) = ang_yi(17,18)
     $                      + ang_yi(18,18)
     $                      + ang_yi(19,18)
     $                      + ang_yi(20,18)
     $                      + ang_yi(21,18)
     $                      + ang_yi(22,18)
     $                      + ang_yi(23,18)
     $                      + ang_yi(24,18)
     $                      + ang_yi(25,18)
                    angi(3) = ang_yi(17,19)
     $                      + ang_yi(18,19)
     $                      + ang_yi(19,19)
     $                      + ang_yi(20,19)
     $                      + ang_yi(21,19)
     $                      + ang_yi(22,19)
     $                      + ang_yi(23,19)
     $                      + ang_yi(24,19)
     $                      + ang_yi(25,19)
                    angi(4) = ang_yi(17,20)
     $                      + ang_yi(18,20)
     $                      + ang_yi(19,20)
     $                      + ang_yi(20,20)
     $                      + ang_yi(21,20)
     $                      + ang_yi(22,20)
     $                      + ang_yi(23,20)
     $                      + ang_yi(24,20)
     $                      + ang_yi(25,20)
                    angi(5) = ang_yi(17,21)
     $                      + ang_yi(18,21)
     $                      + ang_yi(19,21)
     $                      + ang_yi(20,21)
     $                      + ang_yi(21,21)
     $                      + ang_yi(22,21)
     $                      + ang_yi(23,21)
     $                      + ang_yi(24,21)
     $                      + ang_yi(25,21)
                    angi(6) = ang_yi(17,22)
     $                      + ang_yi(18,22)
     $                      + ang_yi(19,22)
     $                      + ang_yi(20,22)
     $                      + ang_yi(21,22)
     $                      + ang_yi(22,22)
     $                      + ang_yi(23,22)
     $                      + ang_yi(24,22)
     $                      + ang_yi(25,22)
                    angi(7) = ang_yi(17,23)
     $                      + ang_yi(18,23)
     $                      + ang_yi(19,23)
     $                      + ang_yi(20,23)
     $                      + ang_yi(21,23)
     $                      + ang_yi(22,23)
     $                      + ang_yi(23,23)
     $                      + ang_yi(24,23)
     $                      + ang_yi(25,23)
                    angi(8) = ang_yi(17,24)
     $                      + ang_yi(18,24)
     $                      + ang_yi(19,24)
     $                      + ang_yi(20,24)
     $                      + ang_yi(21,24)
     $                      + ang_yi(22,24)
     $                      + ang_yi(23,24)
     $                      + ang_yi(24,24)
     $                      + ang_yi(25,24)
                    angi(9) = ang_yi(17,25)
     $                      + ang_yi(18,25)
     $                      + ang_yi(19,25)
     $                      + ang_yi(20,25)
     $                      + ang_yi(21,25)
     $                      + ang_yi(22,25)
     $                      + ang_yi(23,25)
     $                      + ang_yi(24,25)
     $                      + ang_yi(25,25)
c
                    angj(1) = ang_yj(17,17)
     $                      + ang_yj(18,17)
     $                      + ang_yj(19,17)
     $                      + ang_yj(20,17)
     $                      + ang_yj(21,17)
     $                      + ang_yj(22,17)
     $                      + ang_yj(23,17)
     $                      + ang_yj(24,17)
     $                      + ang_yj(25,17)
                    angj(2) = ang_yj(17,18)
     $                      + ang_yj(18,18)
     $                      + ang_yj(19,18)
     $                      + ang_yj(20,18)
     $                      + ang_yj(21,18)
     $                      + ang_yj(22,18)
     $                      + ang_yj(23,18)
     $                      + ang_yj(24,18)
     $                      + ang_yj(25,18)
                    angj(3) = ang_yj(17,19)
     $                      + ang_yj(18,19)
     $                      + ang_yj(19,19)
     $                      + ang_yj(20,19)
     $                      + ang_yj(21,19)
     $                      + ang_yj(22,19)
     $                      + ang_yj(23,19)
     $                      + ang_yj(24,19)
     $                      + ang_yj(25,19)
                    angj(4) = ang_yj(17,20)
     $                      + ang_yj(18,20)
     $                      + ang_yj(19,20)
     $                      + ang_yj(20,20)
     $                      + ang_yj(21,20)
     $                      + ang_yj(22,20)
     $                      + ang_yj(23,20)
     $                      + ang_yj(24,20)
     $                      + ang_yj(25,20)
                    angj(5) = ang_yj(17,21)
     $                      + ang_yj(18,21)
     $                      + ang_yj(19,21)
     $                      + ang_yj(20,21)
     $                      + ang_yj(21,21)
     $                      + ang_yj(22,21)
     $                      + ang_yj(23,21)
     $                      + ang_yj(24,21)
     $                      + ang_yj(25,21)
                    angj(6) = ang_yj(17,22)
     $                      + ang_yj(18,22)
     $                      + ang_yj(19,22)
     $                      + ang_yj(20,22)
     $                      + ang_yj(21,22)
     $                      + ang_yj(22,22)
     $                      + ang_yj(23,22)
     $                      + ang_yj(24,22)
     $                      + ang_yj(25,22)
                    angj(7) = ang_yj(17,23)
     $                      + ang_yj(18,23)
     $                      + ang_yj(19,23)
     $                      + ang_yj(20,23)
     $                      + ang_yj(21,23)
     $                      + ang_yj(22,23)
     $                      + ang_yj(23,23)
     $                      + ang_yj(24,23)
     $                      + ang_yj(25,23)
                    angj(8) = ang_yj(17,24)
     $                      + ang_yj(18,24)
     $                      + ang_yj(19,24)
     $                      + ang_yj(20,24)
     $                      + ang_yj(21,24)
     $                      + ang_yj(22,24)
     $                      + ang_yj(23,24)
     $                      + ang_yj(24,24)
     $                      + ang_yj(25,24)
                    angj(9) = ang_yj(17,25)
     $                      + ang_yj(18,25)
     $                      + ang_yj(19,25)
     $                      + ang_yj(20,25)
     $                      + ang_yj(21,25)
     $                      + ang_yj(22,25)
     $                      + ang_yj(23,25)
     $                      + ang_yj(24,25)
     $                      + ang_yj(25,25)
c
                    angsav(1,5,5,llpot,lmpi,lmpj) =
     $               angi(1)*ylmj_t(17) +
     $               angi(2)*ylmj_t(18) +
     $               angi(3)*ylmj_t(19) +
     $               angi(4)*ylmj_t(20) +
     $               angi(5)*ylmj_t(21) +
     $               angi(6)*ylmj_t(22) +
     $               angi(7)*ylmj_t(23) +
     $               angi(8)*ylmj_t(24) +
     $               angi(9)*ylmj_t(25)
c
                    angsav(2,5,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(17,1) +
     $               angj(2)*grylmi_t(18,1) +
     $               angj(3)*grylmi_t(19,1) +
     $               angj(4)*grylmi_t(20,1) +
     $               angj(5)*grylmi_t(21,1) +
     $               angj(6)*grylmi_t(22,1) +
     $               angj(7)*grylmi_t(23,1) +
     $               angj(8)*grylmi_t(24,1) +
     $               angj(9)*grylmi_t(25,1)
                    angsav(3,5,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(17,2) +
     $               angj(2)*grylmi_t(18,2) +
     $               angj(3)*grylmi_t(19,2) +
     $               angj(4)*grylmi_t(20,2) +
     $               angj(5)*grylmi_t(21,2) +
     $               angj(6)*grylmi_t(22,2) +
     $               angj(7)*grylmi_t(23,2) +
     $               angj(8)*grylmi_t(24,2) +
     $               angj(9)*grylmi_t(25,2)
                    angsav(4,5,5,llpot,lmpi,lmpj) =
     $               angj(1)*grylmi_t(17,3) +
     $               angj(2)*grylmi_t(18,3) +
     $               angj(3)*grylmi_t(19,3) +
     $               angj(4)*grylmi_t(20,3) +
     $               angj(5)*grylmi_t(21,3) +
     $               angj(6)*grylmi_t(22,3) +
     $               angj(7)*grylmi_t(23,3) +
     $               angj(8)*grylmi_t(24,3) +
     $               angj(9)*grylmi_t(25,3)
c
                    angsav(5,5,5,llpot,lmpi,lmpj) =
     $               grylmj_t(17,1)*angi(1) +
     $               grylmj_t(18,1)*angi(2) +
     $               grylmj_t(19,1)*angi(3) +
     $               grylmj_t(20,1)*angi(4) +
     $               grylmj_t(21,1)*angi(5) +
     $               grylmj_t(22,1)*angi(6) +
     $               grylmj_t(23,1)*angi(7) +
     $               grylmj_t(24,1)*angi(8) +
     $               grylmj_t(25,1)*angi(9)
                    angsav(6,5,5,llpot,lmpi,lmpj) =
     $               grylmj_t(17,2)*angi(1) +
     $               grylmj_t(18,2)*angi(2) +
     $               grylmj_t(19,2)*angi(3) +
     $               grylmj_t(20,2)*angi(4) +
     $               grylmj_t(21,2)*angi(5) +
     $               grylmj_t(22,2)*angi(6) +
     $               grylmj_t(23,2)*angi(7) +
     $               grylmj_t(24,2)*angi(8) +
     $               grylmj_t(25,2)*angi(9)
                    angsav(7,5,5,llpot,lmpi,lmpj) =
     $               grylmj_t(17,3)*angi(1) +
     $               grylmj_t(18,3)*angi(2) +
     $               grylmj_t(19,3)*angi(3) +
     $               grylmj_t(20,3)*angi(4) +
     $               grylmj_t(21,3)*angi(5) +
     $               grylmj_t(22,3)*angi(6) +
     $               grylmj_t(23,3)*angi(7) +
     $               grylmj_t(24,3)*angi(8) +
     $               grylmj_t(25,3)*angi(9)
c
                  else
c
                    call STOPXERR( 'nlangf-5' )
c
                  endif
c
                enddo
c
              endif
c
            enddo
c
  290     continue
c
  301   continue
  302 continue
c
c    That's all Folks!
c
      RETURN
      END
