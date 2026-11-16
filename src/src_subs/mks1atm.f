c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKS1ATM
c
c
      subroutine MKS1ATM( s1atom,
     $ norb,nk,ncplx, natm,ntyp, noad, itypa, norba,
     $ smat )
c---------------------------------------------------------------
c Purpose: distill single atom overlaps from total overlap matrix
c
c Written: P.A. Schultz, 15-May-2001 for v2.47
c
c Revision history:
c  10Sep08-PAS/2.62: clean and adapt for striped-matrices
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input arrays:
      DIMENSION  itypa(natm), norba(ntyp)
      DIMENSION  smat(norb,norb,nk)
c Output array:
      DIMENSION  s1atom(noad,noad,nk,natm)
c
c >>>> EXECUTABLE CODE:
c
C      write(6,*) 'SDEV/mks1atm: ijkmat'
      jnext = 0
      do 1000 jatm=1,natm
        jtyp = itypa(jatm)
        norbatj = norba(jtyp)
        jstart = jnext
        jnext = jstart + norbatj
c
        j0 = jstart
c
        do  k=1,nk
c
          do  jorbat=1,norbatj
            do  iorbat=1,norbatj
              s1atom(iorbat,jorbat,k,jatm) = smat(j0+iorbat,j0+jorbat,k)
            enddo
          enddo
c
        enddo
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKS1STR
c
c
      subroutine MKS1STR( s1atom,
     $ norb,nk,ncplx, natm,ntyp, noad, itypa, norba,
     $ myj0,myatm1,myatm2,
     $ nprocstr, master,iproc,icomm,
     $ smat )
c---------------------------------------------------------------
c Purpose: distill single atom overlaps from total overlap matrix
c
c Written: P.A. Schultz, 15-May-2001 for v2.47
c          based on P.J. Feibelman's original
c
c Revision history:
c  25Sep08-PAS/2.62: pass ncplx into this routine
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  itypa(natm), norba(ntyp)
      DIMENSION  smat(nk,*)
c output array:
      DIMENSION  s1atom(noad,noad,nk,natm)
c
c >>>> EXECUTABLE CODE:
c
C      write(6,*) 'SDEV/mks1str: stripe-parallel'
      myatm = 0
      jjnext = 0
      jnext = myj0
      do 1000 jatm=myatm1,myatm2
        myatm = myatm + 1
        jtyp = itypa(jatm)
        norbatj = norba(jtyp)
        jstart = jnext
        jnext = jstart + norbatj
c
        matij = norbatj*norbatj
        njstripe = norb - jstart
        matjstr = norbatj*njstripe
        jjstart = jjnext
        jjnext = jjstart + matjstr*ncplx
c
        j0 = jstart
c
        ij0 = jjstart
        ij0 = ij0 - norbatj
        ijc0 = ij0 + matij
c
        do 420 jorbat=1,norbatj
c
          if( ncplx.eq.2 .and. jorbat.gt.1 )then
c           Need to flip i-j counting in stripe to pick up complex
            ijc = ijc0 + jorbat
            do 412 iorbat=1,jorbat-1
              ijc = ijc + norbatj
              do  k=1,nk
                s1atom(iorbat,jorbat,k,myatm) = smat(k,ijc)
              enddo
 412        continue
          endif
c
          ij0 = ij0 + norbatj
          do 410 iorbat=jorbat,norbatj
            ij = ij0 + iorbat
            do  k=1,nk
              s1atom(iorbat,jorbat,k,myatm) = smat(k,ij)
            enddo
  410     continue
c
  420   continue
c
c       Close Atom i loop:
 1000 continue
c
      if( nprocstr.gt.1 )then
        nvatm = noad*noad*nk
        nvtotl = nvatm*natm
        nvlocl = nvatm*myatm
        call STRMERGE( nprocstr, master,iproc,icomm, nvtotl,nvlocl,
     $   s1atom )
      endif
c
c    That's all Folks!
c
      RETURN
      END
