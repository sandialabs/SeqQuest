c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRHOEF1B
c
c
      subroutine GRHOEF1B( ibox, numpop,
     $ norb, natm,ntyp, nshld,  itypa,numshl,lshel,
     $ nbox,nrec,inrec, nrbox, nobox,inboxs,
     $ wvfcns, rhobox, wkbi,
     $ eigvec, evbx )
c---------------------------------------------------------------
c Purpose: computes electron density in coordinate space,
c          via grid eigenfunctions (real)
c          Version using level-3 blas
c
c Written: Peter A. Schultz,  4-Jun-1999, for v2.36
c          Based on a suggestion by Dave Raczkowski
c
c Revision history:
c  23Jul01-PAS/2.49: eigvec indices
c---------------------------------------------------------------
c
c Notes:
c  This routine first takes the full dense eigenvector matrix in
c  eigvec, and packs the part of it used in this box into evbx.
c  This packed array is then sent into the level-3 blas routine
c  sgemm for the cpu-heavy arithmetic.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  wvfcns, wkbi
      REAL  eigvec, evbx
      REAL  rhobox
      REAL  r4zero,r4one
c
c input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp)
      DIMENSION  inboxs(nrec,nbox)
      DIMENSION  wvfcns(nrbox,nobox)
      DIMENSION  eigvec(numpop,norb)
c output array:
      DIMENSION  rhobox(nrbox)
c scratch arrays:
      DIMENSION  wkbi(nrbox,numpop)
      DIMENSION  evbx(numpop,nobox)
c local declarations:
      DATA  r4zero,r4one / 0.e0,1.e0 /
      LOGICAL    EVEN
c statement function:
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
      inbox = 0
c
c  Loop over atoms/shells
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
            i = i0
            ib = ib0
            do 170 iorbsh=1,norbshi
              i = i + 1
              ib = ib + 1
c
              do 70 n=1,numpop
c
c Pack up box-specific eigenvector array:
                evbx(n,ib) = eigvec(n,i)
c
   70         continue
  170       continue
c
  180     continue
  200   continue
c
c Do heavy lifting inside level-3 blas routine:
        call SGEMM( 'n','t', nrbox,numpop,nobox, r4one,
     $   wvfcns,nrbox, evbx,numpop, r4zero, wkbi,nrbox )
c
c Have occup-weighted eigenfunctions, build density:
        do  n=1,numpop
          do  nr=1,nrbox
            rhobox(nr) = rhobox(nr) + wkbi(nr,n)**2
          enddo
        enddo
c
c    That's all Folks!
c
      RETURN
      END
