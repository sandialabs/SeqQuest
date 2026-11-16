c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EQTRIP
c
c
      subroutine EQTRIP( ntest, iateq,jateq,iatm,jatm,nuc,
     $ natm, itypa,ratm )
c---------------------------------------------------------------
c Purpose: look through previous atom triplets to see if earlier
c          triplet of same atom types separated by same vectors.
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION itypa(natm),ratm(3,natm)
      LOGICAL  ntest
c NB: this term must be identical in eqtrip and gttrip and it is a
c     good idea to use some "odd" number to prevent accidental match
      DATA  rtol / 9.96340732d-7 /
c
      ntest = .false.
      if( iatm.ge.jatm )then
        ntrip = jatm + (iatm-1)*natm + (nuc-1)*natm**2
      else
        ntrip = iatm + (jatm-1)*natm + (nuc-1)*natm**2
      endif
      do 100 nuceq=1,natm
        if( itypa(nuceq).ne.itypa(nuc) ) goto 100
        do 50 iateq=1,natm
          if( itypa(iateq).ne.itypa(iatm) ) goto 50
          do 10 jr=1,3
            dr = ratm(jr,iatm)  - ratm(jr,nuc)
     $         - ratm(jr,iateq) + ratm(jr,nuceq)
            if( ABS(dr).gt.rtol ) goto 50
   10     continue
          do 40 jateq=1,iateq
            neqtrip = jateq + (iateq-1)*natm + (nuceq-1)*natm**2
            if( neqtrip.ge.ntrip ) RETURN
            if( itypa(jateq).ne.itypa(jatm) ) goto 40
            do 20 jr=1,3
              dr = ratm(jr,jatm)  - ratm(jr,nuc)
     $           - ratm(jr,jateq) + ratm(jr,nuceq)
              if( ABS(dr).gt.rtol ) goto 40
   20       continue
            ntest = .true.
            RETURN
   40     continue
   50   continue
  100 continue
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GTTRIP
c
c
      entry GTTRIP( ntest, iateq,jateq,nateq,iatm,jatm,nuc,
     $ natm, itypa,ratm )
c---------------------------------------------------------------
c Purpose: look through later atom pairs to see if later
c          triplet of same atom types separated by same vectors.
c---------------------------------------------------------------
c
      ntest = .false.
      if( iatm.ge.jatm )then
        ntrip = jatm + (iatm-1)*natm + (nuc-1)*natm**2
      else
        ntrip = iatm + (jatm-1)*natm + (nuc-1)*natm**2
      endif
      if( itypa(nateq).ne.itypa(nuc) )RETURN
      do 150 iateq=1,natm
        if( itypa(iateq).ne.itypa(iatm) ) goto 150
        do 110 jr=1,3
          dr = ratm(jr,iatm)  - ratm(jr,nuc)
     $       - ratm(jr,iateq) + ratm(jr,nateq)
          if( ABS(dr).gt.rtol ) goto 150
  110   continue
        do 140 jateq=1,iateq
          neqtrip = jateq + (iateq-1)*natm + (nateq-1)*natm**2
          if( neqtrip.le.ntrip ) goto 140
          if( itypa(jateq).ne.itypa(jatm) ) goto 140
          do 120 jr=1,3
            dr = ratm(jr,jatm)  - ratm(jr,nuc)
     $         - ratm(jr,jateq) + ratm(jr,nateq)
            if( ABS(dr).gt.rtol ) goto 140
  120     continue
          ntest = .true.
          RETURN
  140   continue
  150 continue
c
      RETURN
c
      END
