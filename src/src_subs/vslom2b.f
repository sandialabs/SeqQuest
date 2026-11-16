c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VSLOM2B
c
c
      subroutine VSLOM2B( ibox,
     $ norb,nk, natm,ntyp, nshld,  itypa,numshl,lshel,
     $ nbox,nrec,inrec, nrbox, nobox,inboxs,
     $ wvfcns, vslobox, wkbi,
     $ vkij,vkbx )
c---------------------------------------------------------------
c Purpose: compute matrix elements of grid potentials
c          over complex grid basis orbitals
c          Version using level-3 blas
c
c Written: Peter A. Schultz,  9-Jun-1999, for v2.36
c          Based on suggestion by Dave Raczkowski
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Notes:
c  The matrix elements are computed over the grid for wave fcns
c  in this box using the level-3 blas routine for the cpu-heavy
c  arithmetic, and are loaded into the full dense matrix after.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL  wvfcns, wkbi
      REAL  vkij, vkbx
      REAL  vslobox
      REAL  c4zero,c4one
c
c output arrays:
      DIMENSION  vkij(nk,norb,norb)
c input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp)
      DIMENSION  inboxs(nrec,nbox)
      DIMENSION  wvfcns(2,nrbox,nk,nobox)
      DIMENSION  vslobox(nrbox)
c scratch arrays:
      DIMENSION  wkbi(2,nrbox,nobox), vkbx(2,nobox,nobox,nk)
c local declarations:
      DIMENSION  c4zero(2),c4one(2)
      DATA  c4zero / 0.e0,0.e0 /, c4one / 1.e0,0.e0 /
      LOGICAL    EVEN
c statement function:
      EVEN(int) = (AND(int,1).eq.0)
c
c >>>> EXECUTABLE CODE:
c
      do  k=1,nk
c
c Put potentials into left grid orbitals
        do  ib=1,nobox
          do  nr=1,nrbox
            wkbi(1,nr,ib) =  vslobox(nr)*wvfcns(1,nr,k,ib)
            wkbi(2,nr,ib) = -vslobox(nr)*wvfcns(2,nr,k,ib)
          enddo
        enddo
c
c Do heavy lifting inside level-3 blas routine, k-by-k:
        call CGEMM( 't','n', nobox,nobox,nrbox, c4one,
     $   wvfcns(1,1,k,1),nk*nrbox,  wkbi, nrbox,
     $   c4zero, vkbx(1,1,1,k), nobox )
c
      enddo
c
c  Put matrix elements into full dense matrix:
c
      inbox = 0
      jnbox = 0
c
c    Loop over atoms/shells
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
                    do 60 k=1,nk
c
                      vkij(k,i,j) = vkij(k,i,j) + vkbx(1,jb,ib,k)
                      if( i.ne.j )
     $                 vkij(k,j,i) = vkij(k,j,i) - vkbx(2,jb,ib,k)
c
   60               continue
c
   70             continue
c
  170           continue
c
   80         continue
  100       continue
c
  180     continue
  200   continue
c
c    That's all Folks!
c
      RETURN
      END
