c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VSLOF1
c
c
      subroutine VSLOF1( ibox, ndim, frcslo,strslo,
     $ norb, natm,ntyp, nshld,nald, nlat,
     $ itypa,numshl,lshel,nala,ala,cala, alamin,
     $ ratm, rlat, coskr,sinkr,
     $ hh, nbox,nrec,inrec, nrbox,nobox,inboxs, boxrad,boxctr,
     $ wf, rmesh, vbox, wkbox,wkbi,
     $ dphix,dphiy,dphiz,  dphiuv,
     $ dmat )
c---------------------------------------------------------------
c Purpose: computes force/stress from grid potentials in a box
c          for real basis orbitals
c
c Written: Peter A. Schultz, 10-June-1999, for v2.36
c
c Revision history:
c  23Jan09-PAS?2.63: k-parallel
c  27Sep99-PAS/2.40: full three-dimensional stress tensor
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  wf, wkbi
      REAL  dphix,dphiy,dphiz,  dphiuv
      REAL  vbox
      REAL  dmat, dmatr
      REAL  fxint,fyint,fzint,  strbox
      REAL  r4zero
c
c output arrays:
      DIMENSION  frcslo(3,natm), strslo(3,3)
c input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp),
     $ nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp)
      DIMENSION  ratm(3,natm), rlat(3,nlat)
      DIMENSION  coskr(nlat),sinkr(nlat)
      DIMENSION  hh(3),inboxs(nrec,nbox)
      DIMENSION  wf(nrbox,nobox)
      DIMENSION  rmesh(3,nrbox)
      DIMENSION  vbox(nrbox)
      DIMENSION  dmat(norb,norb)
c scratch arrays:
      DIMENSION  wkbi(nrbox), wkbox(nrbox,*)
      DIMENSION  dphix(nrbox),dphiy(nrbox),dphiz(nrbox)
      DIMENSION  dphiuv(ndim,ndim,nrbox)
c local declarations:
      DIMENSION  strbox(3,3)
      DIMENSION  boxctr(3)
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
c Loop over atoms/shells
c
       irec = 0
       iatshl = 0
       inext = 0
       ibnext = 0
       do 200 iatm=1,natm
         ityp = itypa(iatm)
         almina = alamin(ityp)
c
         do 180 ishl=1,numshl(ityp)
           lshli = lshel(ishl,ityp)
           norbshi = 2*lshli + 1
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
           almins = ala(1,ishl,ityp)
           do  ial=2,nala(ishl,ityp)
             ali = ala(ial,ishl,ityp)
             if( ali.lt.almins ) almins = ali
           enddo
           almin = almina
c
           i = i0
           ib = ib0
           do 170 iorbsh=1,norbshi
            i = i + 1
            ib = ib + 1
c
c  For force calculation:
c   compute derivatives of i-th wave function
c  For stress calculation:
c   compute (x-X) times x-derivative of i-th wave function
c   compute (y-Y) times y-derivative of i-th wave function
c
            call GRADWF( ibox,     1, ndim, iatm,ishl,iorbsh, almin,
     $        1,0,1, natm,ntyp, nshld,nald, nlat,
     $       itypa, numshl,lshel,nala,ala,cala, ratm, rlat, coskr,sinkr,
     $       nrbox, boxrad, boxctr,
     $       rmesh,  wkbox(1,1), wkbox(1,2),wkbox(1,3),
c -->        rv-3io  rsq-s       sumg-s     sum2ag-s
     $       wkbox(1,4),wkbox(1,5),wkbox(1,6),
c -->        dpdx-s     dpdy-s     dpdz-s
     $       dphix,dphiy,dphiz,  dphiuv )
c
c    Clear SUM(j) of dmij*phi-j
c
            call MKSZERO( nrbox, wkbi )
c
c    Atom j loop:
c
            jrec = 0
            jatshl = 0
            jnext = 0
            jbnext = 0
            do 100 jatm=1,natm
              jtyp = itypa(jatm)
c
              do 80 jshl=1,numshl(jtyp)
               lshlj = lshel(jshl,jtyp)
               norbshj = 2*lshlj + 1
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
               j = j0
               jb = jb0
               do 70 jorbsh=1,norbshj
                j = j + 1
                jb = jb + 1
c
                dmatr = dmat(i,j)
                call SAXPY( nrbox, dmatr, wf(1,jb),1, wkbi,1)
Cc
C                do  nr=1,nrbox
C                  wkbi(nr) = wkbi(nr) + dmatr*wf(nr,jb)
C                enddo
c
c    Close atom j orbital, atom loops:
c
   70          continue
c
   80         continue
  100       continue
c
c  Complete force integrals ...
c
            fxint = r4zero
            fyint = r4zero
            fzint = r4zero
            call MKSZERO( 9, strbox )
c
            do  nr=1,nrbox
              wkbi(nr) = vbox(nr)*wkbi(nr)
            enddo
c
            do  nr=1,nrbox
              fxint = fxint + dphix(nr)*wkbi(nr)
              fyint = fyint + dphiy(nr)*wkbi(nr)
              fzint = fzint + dphiz(nr)*wkbi(nr)
            enddo
c
            if( ndim.eq.1 )then
              do  nr=1,nrbox
                strbox(1,1) = strbox(1,1) + dphiuv(1,1,nr)*wkbi(nr)
              enddo
            elseif( ndim.eq.2 )then
              do  nr=1,nrbox
                strbox(1,1) = strbox(1,1) + dphiuv(1,1,nr)*wkbi(nr)
                strbox(2,1) = strbox(2,1) + dphiuv(2,1,nr)*wkbi(nr)
                strbox(1,2) = strbox(1,2) + dphiuv(1,2,nr)*wkbi(nr)
                strbox(2,2) = strbox(2,2) + dphiuv(2,2,nr)*wkbi(nr)
              enddo
            elseif( ndim.eq.3 )then
              do  nr=1,nrbox
                strbox(1,1) = strbox(1,1) + dphiuv(1,1,nr)*wkbi(nr)
                strbox(2,1) = strbox(2,1) + dphiuv(2,1,nr)*wkbi(nr)
                strbox(3,1) = strbox(3,1) + dphiuv(3,1,nr)*wkbi(nr)
                strbox(1,2) = strbox(1,2) + dphiuv(1,2,nr)*wkbi(nr)
                strbox(2,2) = strbox(2,2) + dphiuv(2,2,nr)*wkbi(nr)
                strbox(3,2) = strbox(3,2) + dphiuv(3,2,nr)*wkbi(nr)
                strbox(1,3) = strbox(1,3) + dphiuv(1,3,nr)*wkbi(nr)
                strbox(2,3) = strbox(2,3) + dphiuv(2,3,nr)*wkbi(nr)
                strbox(3,3) = strbox(3,3) + dphiuv(3,3,nr)*wkbi(nr)
              enddo
            endif
c
            frcslo(1,iatm) = frcslo(1,iatm) + fxint
            frcslo(2,iatm) = frcslo(2,iatm) + fyint
            frcslo(3,iatm) = frcslo(3,iatm) + fzint
c
            do  jd=1,ndim
              do  id=1,ndim
                strslo(id,jd) = strslo(id,jd) + strbox(id,jd)
              enddo
            enddo
c
c    Close atom i orbital, atom loops:
c
  170      continue
  180     continue
  200   continue
c
c    That's all Folks!
c
      RETURN
      END
c         
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:        
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
