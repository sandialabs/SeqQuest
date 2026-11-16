c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NLINTF
c
c
      subroutine NLINTF( enpri,enprj, engri,engrj,
     $ vnln,vnli,vnlj, rimag,rjmag, lmi,lmj,
     $ limxp1,ljmxp1,lmxp1,lpowmx,lpowi,lpowj,lstmx,
     $ ai,gradai, aj,gradaj, angsav )
c---------------------------------------------------------------
c Purpose: do angular sums NLANGF has not already done
c
c Written: Peter A. Schultz, 28-August-1998, for v2.25
c
c Revision History:
c  15Oct99-PAS/2.41: cleaned up, and reorganized for further opt
c   5Oct99-DBR/2.41: optimization: reordered angv
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c output arrays:
      DIMENSION  engri(3),engrj(3)
c input arrays:
      DIMENSION  vnln(limxp1+lpowi,ljmxp1+lpowj,lstmx)
      DIMENSION  vnli(limxp1+lpowi,ljmxp1+lpowj,lstmx)
      DIMENSION  vnlj(limxp1+lpowi,ljmxp1+lpowj,lstmx)
      DIMENSION  ai(9),gradai(3,9), aj(9),gradaj(3,9)
      DIMENSION  angsav(7,limxp1,ljmxp1,lmxp1,9,9)
c local declarations:
      LOGICAL    ri_gt_0, rj_gt_0
      DIMENSION  aigraj(3),ajgrai(3)
      DIMENSION  lang(9)
      DATA  lang / 0, 1,1,1, 2,2,2,2,2 /
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      enpri = zero
      enprj = zero
c
      engri(1) = zero
      engri(2) = zero
      engri(3) = zero
      engrj(1) = zero
      engrj(2) = zero
      engrj(3) = zero
c
      eneli = zero
      enelj = zero
c
      ri_gt_0 = .false.
      rj_gt_0 = .false.
      if( rimag .ge. 1.d-12 ) ri_gt_0 = .true.
      if( rjmag .ge. 1.d-12 ) rj_gt_0 = .true.
c
      do 200 lmpj=1,lmj
        sumj = DABS( aj(lmpj) ) + DABS( gradaj(1,lmpj) )
     $   + DABS( gradaj(2,lmpj) ) + DABS( gradaj(3,lmpj) )
        if( sumj .eq. zero ) goto 200
c
        lpj = lang(lmpj)
c
        do 100 lmpi=1,lmi
          sumi = DABS( ai(lmpi) ) + DABS( gradai(1,lmpi) )
     $     + DABS( gradai(2,lmpi) ) + DABS( gradai(3,lmpi) )
          if( sumi .eq. zero ) goto 100
c
          lpi = lang(lmpi)
c
          aiaj = ai(lmpi)*aj(lmpj)
c
          ajgrai(1) = gradai(1,lmpi)*aj(lmpj)
          ajgrai(2) = gradai(2,lmpi)*aj(lmpj)
          ajgrai(3) = gradai(3,lmpi)*aj(lmpj)
          aigraj(1) = ai(lmpi)*gradaj(1,lmpj)
          aigraj(2) = ai(lmpi)*gradaj(2,lmpj)
          aigraj(3) = ai(lmpi)*gradaj(3,lmpj)
c
          sum1i  = zero
          sum1i1 = zero
          sum1j  = zero
          sum1j1 = zero
          sum1n  = zero
          sum2n  = zero
          sum3n  = zero
          sum4n  = zero
          sum5n  = zero
          sum6n  = zero
          sum7n  = zero
c
          do  llpot=1,lmxp1
            lstore = lpi + lpj + 1 + lpowmx*(llpot-1)
c
            do  llj=1,ljmxp1
              do  lli=1,limxp1
c
                ang1 = angsav(1,lli,llj,llpot,lmpi,lmpj)
                venij = vnln(lli,llj,lstore)
c
                sum1i  = sum1i  + ang1*vnli(lli,llj,lstore)
                sum1i1 = sum1i1 + ang1*vnli(lli+1,llj,lstore+1)
c
                sum1j  = sum1j  + ang1*vnlj(lli,llj,lstore)
                sum1j1 = sum1j1 + ang1*vnlj(lli,llj+1,lstore+1)
c
                sum1n  = sum1n  + ang1*venij
c
                sum2n  = sum2n  + venij*
     $           angsav(2,lli,llj,llpot,lmpi,lmpj)
                sum3n  = sum3n  + venij*
     $           angsav(3,lli,llj,llpot,lmpi,lmpj)
                sum4n  = sum4n  + venij*
     $           angsav(4,lli,llj,llpot,lmpi,lmpj)
c
                sum5n  = sum5n  + venij*
     $           angsav(5,lli,llj,llpot,lmpi,lmpj)
                sum6n  = sum6n  + venij*
     $           angsav(6,lli,llj,llpot,lmpi,lmpj)
                sum7n  = sum7n  + venij*
     $           angsav(7,lli,llj,llpot,lmpi,lmpj)
c
              enddo
            enddo
c
          enddo
c
          eneli = eneli + aiaj*sum1i
          if( ri_gt_0 ) enpri = enpri - aiaj*sum1i1
c
          enelj = enelj + aiaj*sum1j
          if( rj_gt_0 ) enprj = enprj - aiaj*sum1j1
c
          engri(1) = engri(1) - ( aiaj*sum2n + ajgrai(1)*sum1n )
          engri(2) = engri(2) - ( aiaj*sum3n + ajgrai(2)*sum1n )
          engri(3) = engri(3) - ( aiaj*sum4n + ajgrai(3)*sum1n )
c
          engrj(1) = engrj(1) - ( aiaj*sum5n + aigraj(1)*sum1n )
          engrj(2) = engrj(2) - ( aiaj*sum6n + aigraj(2)*sum1n )
          engrj(3) = engrj(3) - ( aiaj*sum7n + aigraj(3)*sum1n )
c
  100   continue
  200 continue
c
      if( ri_gt_0 )then
        enpri = enpri/rimag + eneli
      else
        enpri = eneli
      endif
c
      if( rj_gt_0 )then
        enprj = enprj/rjmag + enelj
      else
        enprj = enelj
      endif
c
c    That's all Folks!
c
      RETURN
      END
