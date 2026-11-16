c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GRHOEF1
c
c
      subroutine GRHOEF1( ibox, numpop,
     $ norb, natm,ntyp, nshld,  itypa,numshl,lshel,
     $ nbox,nrec,inrec, nrbox, nobox,inboxs,
     $ wvfcns, rhobox, wkbi,
     $ eigvec )
c---------------------------------------------------------------
c Purpose: computes electron density in coordinate space,
c          via grid eigenfunctions (real)
c
c Written: Peter A. Schultz, 11-May-1999, for v2.35,
c          Based on a suggestion by Dave Raczkowski
c
c Revision history:
c  23Jul01-PAS/2.49: eigvec indices
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  wvfcns, wkbi
      REAL  eigvec, cni
      REAL  rhobox
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
c local declarations:
      LOGICAL    EVEN
c statement function:
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
c Clear out grid eigenvector wavefunctions
      call MKSZERO( numpop*nrbox, wkbi )
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
              iwv = iwvorb + 1
              i = i + 1
c
              do 70 n=1,numpop
                cni = eigvec(n,i)
c
c Use blas to build occupation-weighted grid eigenfunctions ...
                call SAXPY( nrbox, cni, wvfcns(1,iwv),1,
     $           wkbi(1,n),1 )
c
   70         continue
c
              iwvorb = iwvorb + 1
  170       continue
c
            iwvshl = iwvshl + norbshi
  180     continue
  200   continue
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
