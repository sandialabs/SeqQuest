c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRHODM2
c
c
      subroutine GRHODM2( ibox,
     $ norb,nk, natm,ntyp, nshld,  itypa,numshl,lshel,
     $ nbox,nrec,inrec, nrbox, nobox,inboxs,
     $ wvfcns, rhobox, wkbi,
     $ dmat )
c---------------------------------------------------------------
c Purpose: computes electron density in coordinate space,
c          though a density matrix scheme (complex)
c
c Written: Peter A. Schultz, 11-May-1999, for v2.35
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c  Note that dmat() and wvfcns() are complex quantities, but
c  are represented as double-length (1:2,*) real quantities.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  wvfcns, wkbi
      REAL  dmat, dmatr,dmati,dmatv
      COMPLEX  dmatc
      REAL  rhobox
c
c input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp)
      DIMENSION  inboxs(nrec,nbox)
      DIMENSION  wvfcns(2,nrbox,nk*nobox)
      DIMENSION  dmat(2,nk,norb,norb)
c output array:
      DIMENSION  rhobox(nrbox)
c scratch arrays:
      DIMENSION  wkbi(2,nrbox,nk*5)
c local declarations:
      DIMENSION  dmatv(2)
      EQUIVALENCE  (dmatv(1),dmatr),(dmatv(2),dmati)
      EQUIVALENCE  (dmatv,dmatc)
      LOGICAL    EVEN
c statement function:
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
      inbox = 0
      jnbox = 0
c
c Loop over atoms/shells
c
        irec = 0
        iatshl = 0
        iwvshl = 0
        inext = 0
        do 200 iatm=1,natm
          ityp = itypa(iatm)
c
          do 180 ishl=1,numshl(ityp)
            norbshi = 2*lshel(ishl,ityp) + 1
            inbox = inbox/2
            if( inrec*(iatshl/inrec).eq.iatshl )then
              irec = irec + 1
              inbox = inboxs(irec,ibox)
            endif
            iatshl = iatshl + 1
            i0 = inext
            inext = inext + norbshi
            if( EVEN(inbox) ) goto 180
c
c Clear workspace:
            call MKSZERO( 2*nk*norbshi*nrbox, wkbi )
c
            jrec = 0
            jatshl = 0
            jwvshl = 0
            jnext = 0
            do 100 jatm=1,iatm
              jtyp = itypa(jatm)
              nshlj = numshl(jtyp)
              if( jatm.eq.iatm ) nshlj = ishl
c
              do 80 jshl=1,nshlj
                norbshj = 2*lshel(jshl,jtyp) + 1
c
                jnbox = jnbox/2
                if( inrec*(jatshl/inrec).eq.jatshl )then
                  jrec = jrec + 1
                  jnbox = inboxs(jrec,ibox)
                endif
                jatshl = jatshl + 1
                j0 = jnext
                jnext = jnext + norbshj
                if( EVEN(jnbox) ) goto 80
c
                joff = 0
                if( jatm.eq.iatm .and. jshl.eq.ishl )then
                  norbshj = 0
                  joff = 1
                endif
c
                iwvxorb = 0
c                iwvorb = iwvshl
                i = i0
                do 170 iorbsh=1,norbshi
                  i = i + 1
c
                  norbshj = norbshj + joff
c
                  jwvorb = jwvshl
                  j = j0
                  do 70 jorbsh=1,norbshj
                    j = j + 1
c
                    iwvx = iwvxorb
c                    iwv = iwvorb
                    jwv = jwvorb
c
                    do 60 k=1,nk
                      iwvx = iwvx + 1
c                      iwv = iwv + 1
                      jwv = jwv + 1
                      dmatr =  dmat(1,k,i,j)
                      dmati = -dmat(2,k,i,j)
c
c use blas to substitute for explicit loop ...
                      call CAXPY( nrbox, dmatc, wvfcns(1,1,jwv),1,
     $                 wkbi(1,1,iwvx),1 )
c
c  ... and this is the reduced op code it replaces:
c                      do  nr=1,nrbox
c                        wkbi(1,nr,iwvx) = wkbi(1,nr,iwvx)
c     $                   + dmatr*wvfcns(1,nr,jwv)
c     $                   - dmati*wvfcns(2,nr,jwv)
c                        wkbi(2,nr,iwvx) = wkbi(2,nr,iwvx)
c     $                   + dmatr*wvfcns(2,nr,jwv)
c     $                   + dmati*wvfcns(1,nr,jwv)
c                      enddo
c  ... and this was the original full code:
c                      do  nr=1,nrbox
c                        rhobox(nr) = rhobox(nr) + dmatr*
c     $                   ( wvfcns(1,nr,iwv)*wvfcns(1,nr,jwv)
c     $                   + wvfcns(2,nr,iwv)*wvfcns(2,nr,jwv) )
c                        rhobox(nr) = rhobox(nr) + dmati*
c     $                   ( wvfcns(2,nr,iwv)*wvfcns(1,nr,jwv)
c     $                   - wvfcns(1,nr,iwv)*wvfcns(2,nr,jwv) )
c                      enddo
c
   60               continue
c
                    jwvorb = jwvorb + nk
   70             continue
                  iwvxorb = iwvxorb + nk
c                  iwvorb = iwvorb + nk
  170           continue
c
                jwvshl = jwvshl + nk*norbshj
   80         continue
  100       continue
c
c Now collect all the terms:
            iwvx = 0
            iwv = iwvshl
            do  iorbsh=1,norbshi
              do  k=1,nk
                iwvx = iwvx + 1
                iwv = iwv + 1
                do  nr=1,nrbox
                  rhobox(nr) = rhobox(nr)
     $             + wvfcns(1,nr,iwv)*wkbi(1,nr,iwvx)
     $             + wvfcns(2,nr,iwv)*wkbi(2,nr,iwvx)
                enddo
              enddo
            enddo
c
            iwvshl = iwvshl + nk*norbshi
  180     continue
  200   continue
c
c    That's all Folks!
c
      RETURN
      END
