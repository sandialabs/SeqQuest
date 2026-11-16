c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRHOEF2
c
c
      subroutine GRHOEF2( ibox, numpop,
     $ norb,nk, natm,ntyp, nshld,  itypa,numshl,lshel,
     $ nbox,nrec,inrec, nrbox, nobox,inboxs,
     $ wvfcns, rhobox, wkbi,
     $ eigvec )
c---------------------------------------------------------------
c Purpose: computes electron density in coordinate space,
c          via grid eigenfunctions (complex)
c
c Written: Peter A. Schultz, 11-May-1999, for v2.35,
c          Based on a suggestion by Dave Raczkowski
c
c Revision history:
c  23Jul01-PAS/2.49: eigvec indices
c---------------------------------------------------------------
c
c  Note that eigvec() and wvfcns() are complex quantities, but
c  are represented as double-length (1:2,*) real quantities.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  wvfcns, wkbi
      REAL  eigvec, cni,cnir,cnii
      COMPLEX  cnic
      REAL  rhobox
c
c input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp)
      DIMENSION  inboxs(nrec,nbox)
      DIMENSION  wvfcns(2,nrbox,nk*nobox)
      DIMENSION  eigvec(2,numpop,nk,norb)
c output array:
      DIMENSION  rhobox(nrbox)
c scratch arrays:
      DIMENSION  wkbi(2,nrbox,nk*numpop)
c local declarations:
      DIMENSION  cni(2)
      EQUIVALENCE  (cni(1),cnir),(cni(2),cnii)
      EQUIVALENCE  (cni,cnic)
      LOGICAL    EVEN
c statement function:
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
c Clear out grid eigenvector wavefunctions
      call MKSZERO( 2*nk*numpop*nrbox, wkbi )
      inbox = 0
c
c  Loop over atoms/shells
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
            iwvorb = iwvshl
            i = i0
            do 170 iorbsh=1,norbshi
              i = i + 1
c
              iwv = iwvorb + 1
              nwv0 = 0
c
              do 60 k=1,nk
c
                nwv = nwv0
                do 70 n=1,numpop
                  nwv = nwv + 1
                  cnir = eigvec(1,n,k,i)
                  cnii = eigvec(2,n,k,i)
c
c Use blas to substitute for explicit loop ...
                  call CAXPY( nrbox, cnic, wvfcns(1,1,iwv),1,
     $             wkbi(1,1,nwv),1 )
c
c  ...and this is the reduced op code it replaces:
c                  do  nr=1,nrbox
c                    wkbi(1,nr,nwv) = wkbi(1,nr,nwv)
c     $               + cnir*wvfcns(1,nr,iwv)
c     $               - cnii*wvfcns(2,nr,iwv)
c                    wkbi(2,nr,nwv) = wkbi(2,nr,nwv)
c     $               + cnir*wvfcns(2,nr,iwv)
c     $               + cnii*wvfcns(1,nr,iwv)
c                  enddo
c
   70           continue
c
                nwv0 = nwv0 + numpop
                iwv = iwv + 1
   60         continue
c
              iwvorb = iwvorb + nk
  170       continue
c
            iwvshl = iwvshl + nk*norbshi
  180     continue
  200   continue
c
c  Have occup-weighted eigenfunctions, build density:
c
        do  nwv=1,nk*numpop
          do  nr=1,nrbox
            rhobox(nr) = rhobox(nr) + wkbi(1,nr,nwv)**2
     $                              + wkbi(2,nr,nwv)**2
          enddo
        enddo
c
c    That's all Folks!
c
      RETURN
      END
