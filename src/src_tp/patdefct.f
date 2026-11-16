c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATDEFCT
c
c
      subroutine PATDEFCT( IWR, idmatfl, ncplx, znuc,rhoatom,defatom,
     $ norb,nk, natm,ntyp, noad, itypa, norba,
     $ near2c, s1atom, orbint,orbofk,
     $ dmat)
c---------------------------------------------------------------
c Purpose: compute nearby-construction local atom charge defects
c
c Written: Peter A. Schultz, June-1993, for v2.10
c
c Revision history:
c   6Jul06-APT/    : merged task parallel and image parallel
c  30Jul05-PAS/2.59: master-only matrix files
c  23Jul01-PAS/2.49: compact dmat; spin-polarized dft
c   3May99-PAS/2.35: use flag to set use of two-center terms
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Apr97-PAS/2.17: naming conventions
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 ,two = 2.d0 )
c
      DIMENSION  znuc(ntyp), rhoatom(natm), defatom(natm,*)
      DIMENSION  itypa(natm), norba(ntyp)
      DIMENSION  s1atom(noad,noad,nk,natm)
      DIMENSION  orbint(2,norb,nk),orbofk(2,natm,norb,nk)
c scratch array:
      DIMENSION  dmat(norb,norb)
c local declarations:
      DATA  lstout / 1 /
c
c >>>> EXECUTABLE CODE:
c
c  Get Local MP info
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
      do 1000 k=1,nk
c
c       Read in density matrix for this k-vector
        if( iproc.eq.master ) call READBIG( idmatfl, norb*norb, dmat )
        if( nprocs.gt.0 ) call MPBCAST8( master, norb*norb, dmat, 
     $       icomm )
c
        i = 0
        do 200 iatm=1,natm
c
          do 190 iorbat=1,norba(itypa(iatm))
            i = i + 1
c
            j = 0
            do 100 jatm=1,iatm
              norbatj = norba(itypa(jatm))
              if( jatm.eq.iatm ) norbatj = iorbat
c
              do 90 jorbat=1,norbatj
                j = j + 1
c
                if( iatm.eq.jatm )then
c
                  si1r = s1atom(iorbat,jorbat,k,iatm)
                  if( iorbat.ne.jorbat )then
                    si1i = s1atom(jorbat,iorbat,k,iatm)
                    defatom(iatm,1) = defatom(iatm,1)
     $               + two*( si1r*dmat(i,j) + si1i*dmat(j,i) )
                  else
                    defatom(iatm,1) = defatom(iatm,1) + dmat(i,i)*si1r
                  endif
c
                elseif( near2c.eq.2 )then
c
C                  si2r = REAL(orbint(i,k)*CONJG(orbofk(iatm,j,k)))
C                  si2i = IMAG(orbint(i,k)*CONJG(orbofk(iatm,j,k)))
C                  sj2r = REAL(orbofk(jatm,i,k)*CONJG(orbint(j,k)))
C                  sj2i = IMAG(orbofk(jatm,i,k)*CONJG(orbint(j,k)))
                  si2r = orbint(1,i,k)*orbofk(1,iatm,j,k)
     $                 + orbint(2,i,k)*orbofk(2,iatm,j,k)
                  si2i = orbint(2,i,k)*orbofk(1,iatm,j,k)
     $                 - orbint(1,i,k)*orbofk(2,iatm,j,k)
                  sj2r = orbofk(1,jatm,i,k)*orbint(1,j,k)
     $                 + orbofk(2,jatm,i,k)*orbint(2,j,k)
                  sj2i = orbofk(2,jatm,i,k)*orbint(1,j,k)
     $                 - orbofk(1,jatm,i,k)*orbint(2,j,k)
                  defatom(iatm,2) = defatom(iatm,2)
     $             + two*( si2r*dmat(i,j) + si2i*dmat(j,i) )
                  defatom(jatm,2) = defatom(jatm,2)
     $             + two*( sj2r*dmat(i,j) + sj2i*dmat(j,i) )
c
                endif
c
   90         continue
c
  100       continue
c
  190     continue
c
  200   continue
c
 1000 continue
c
      if( lstout.gt.0 )then
        write(IWR,*)
        if( lstout.gt.1 )then
          write(IWR,*)'ATDEFCT: nearby-construction atom charge defects'
          write(IWR,900)
        endif
  900   format(2x,'atom',9x,'z-rho',9x,'total',9x,'1-ctr',9x,'2-ctr')
  901   format(2x,i3,3x   ,4f14.8)
        defatm = zero
        def1atm = zero
        do 300 iatm=1,natm
          defgr = znuc(itypa(iatm)) - rhoatom(iatm)
          defnr = defatom(iatm,1) + defatom(iatm,2)
          def1atm = def1atm + defatom(iatm,1)
          defatm = defatm + defnr
          if( lstout.gt.1 )
     $    write(IWR,901)iatm,defgr,defnr,defatom(iatm,1),defatom(iatm,2)
  300   continue
        if( near2c.lt.2 )then
          write(IWR,'(1x,a,2f18.8)')
     $     'one-center local charge defect = ',def1atm
        else
          write(IWR,'(1x,a,2f18.8)')
     $     '1ctr/total local charge defect = ',def1atm,defatm
        endif
      endif
c
c    That's all Folks!
c
      RETURN
      END
