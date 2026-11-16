c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VSLOM1
c
c
      subroutine VSLOM1( ibox,
     $ norb, natm,ntyp, nshld,  itypa,numshl,lshel,
     $ nbox,nrec,inrec, nrbox,nobox,inboxs,
     $ wvfcns, vslobox, wkbi,
     $ vkij )
c---------------------------------------------------------------
c Purpose: compute matrix elements of grid potentials
c          over real grid basis orbitals
c
c Written: Peter A. Schultz, 12-May-1999, for v2.35
c
c Revision history:
c  13Jan12-PAS/2.62: fix bug (declare sdot as real!)
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  wvfcns, wkbi
      REAL  vkij, vintr
      REAL  vslobox
      REAL  r4zero
      REAL  SDOT
c
c output arrays:
      DIMENSION  vkij(norb,norb)
c input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp)
      DIMENSION  inboxs(nrec,nbox)
      DIMENSION  wvfcns(nrbox,nobox)
      DIMENSION  vslobox(nrbox)
c scratch arrays:
      DIMENSION  wkbi(nrbox,5)
c local declarations:
      DATA  r4zero / 0.e0 /
      LOGICAL    EVEN
c statement function:
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
      inbox = 0
      jnbox = 0
c
c    Loop over atoms/shells
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
            iwvx = 0
            iwv = iwvshl
            do  iorbsh=1,norbshi
              iwvx = iwvx + 1
              iwv = iwv + 1
              do  nr=1,nrbox
                wkbi(nr,iwvx) = vslobox(nr)*wvfcns(nr,iwv)
              enddo
            enddo
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
                    iwvx = iwvxorb + 1
c                    iwv = iwvorb + 1
                    jwv = jwvorb + 1
c
c Use blas to substitute for explicit loop ...
                    vintr = SDOT( nrbox, wkbi(1,iwvx),1,
     $               wvfcns(1,jwv),1 )
c
c  ... and this is the reduced op code it replaces:
c                    vintr = r4zero
c                    do  nr=1,nrbox
c                      vintr = vintr + wkbi(nr,iwvx)*wvfcns(nr,jwv)
c                    enddo
c  ... and this was the original full code:
c                    do  nr=1,nrbox
c                    vint =+ wvfcns(nr,iwv)*CONJG(wvfcns(nr,jwv))*v
c                      vintr = vintr + vslobox(nr)*
c     $                 wvfcns(nr,iwv)*wvfcns(nr,jwv)
c                    enddo
c
                    vkij(i,j) = vkij(i,j) + vintr
c
                    jwvorb = jwvorb + 1
   70             continue
c
                  iwvxorb = iwvxorb + 1
c                  iwvorb = iwvorb + 1
  170           continue
c
                jwvshl = jwvshl + norbshj
   80         continue
  100       continue
c
            iwvshl = iwvshl + norbshi
  180     continue
  200   continue
c
c    That's all Folks!
c
      RETURN
      END
