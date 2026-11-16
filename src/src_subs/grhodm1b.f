c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRHODM1B
c
c
      subroutine GRHODM1B( ibox,
     $ norb, natm,ntyp, nshld,  itypa,numshl,lshel,
     $ nbox,nrec,inrec, nrbox, nobox,inboxs,
     $ wvfcns, rhobox, wkbi,
     $ dmat, dmatbx )
c---------------------------------------------------------------
c Purpose: computes electron density in coordinate space,
c          through a density matrix scheme (real)
c          Version using level-3 blas
c
c Written: Peter A. Schultz,  4-Jun-1999, for v2.36
c          Based on a suggestion from Dave Raczkowski
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Notes:
c  This routine first takes the full dense density matrix in
c  dmat, and packs the part of it used in this box into dmatbx.
c  This packed array is then sent into the level-3 blas routine
c  sgemm for the cpu-heavy arithmetic.
c
c  Tried substituting ssymm, the symmetric form of the level-3
c  blas, for sgemm (note code commented by csym), and the code
c  was slower by 1/3.  Hence, use sgemm instead.  4jun99-pas
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  wvfcns, wkbi
      REAL  dmat, dmatbx
      REAL  rhobox
      REAL  r4zero,r4one
csym      REAL  r4half
c
c input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp)
      DIMENSION  inboxs(nrec,nbox)
      DIMENSION  wvfcns(nrbox,nobox)
      DIMENSION  dmat(norb,norb)
c output array:
      DIMENSION  rhobox(nrbox)
c scratch arrays:
      DIMENSION  wkbi(nrbox,nobox)
      DIMENSION  dmatbx(nobox,nobox)
c local declarations:
      DATA  r4zero,r4one / 0.e0,1.e0 /
csym      DATA  r4half / 0.d0 /
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
        inext = 0
        ibnext = 0
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
            ib0 = ibnext
            ibnext = ibnext + norbshi
c
            jrec = 0
            jatshl = 0
            jnext = 0
            jbnext = 0
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
                jb0 = jbnext
                jbnext = jbnext + norbshj
c
                joff = 0
                if( jatm.eq.iatm .and. jshl.eq.ishl )then
                  norbshj = 0
                  joff = 1
                endif
c
                i = i0
                ib = ib0
                do 170 iorbsh=1,norbshi
                  i = i + 1
                  ib = ib + 1
c
                  norbshj = norbshj + joff
c
                  j = j0
                  jb = jb0
                  do 70 jorbsh=1,norbshj
                    j = j + 1
                    jb = jb + 1
c
c Supposed to do triangle, actually do full block
                    dmatbx(jb,ib) = r4zero
                    dmatbx(ib,jb) = dmat(i,j)
csymc
csymc Use symmetric form of blas routine, lower triangle
csym                    if( ib.eq.jb )then
csym                      dmatbx(ib,jb) = dmat(i,j)
csym                    else
csym                      dmatbx(ib,jb) = r4half*dmat(i,j)
csym                    endif
c
   70             continue
  170           continue
c
   80         continue
  100       continue
c
  180     continue
  200   continue
c
c  Now have dmat in this box packed in dmatbx.
c
c Do heavy lifting inside level-3 blas routine:
        call SGEMM( 'n','n', nrbox,nobox,nobox, r4one,
     $   wvfcns,nrbox, dmatbx,nobox, r4zero, wkbi,nrbox )
csymc
csym        call SSYMM( 'r','l', nrbox,nobox, r4one,
csym     $   dmatbx,nobox, wvfcns,nrbox, r4zero, wkbi,nrbox )
c
c Complete density calculation:
        do  ib=1,nobox
          do  nr=1,nrbox
            rhobox(nr) = rhobox(nr)
     $       + wvfcns(nr,ib)*wkbi(nr,ib)
          enddo
        enddo
c
c    That's all Folks!
c
      RETURN
      END
