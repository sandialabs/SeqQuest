c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BFSETUP
c
c
      subroutine BFSETUP( IWR, natm,ntyp,nshld,nald,
     $ norb,norbd,
     $ numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ znuc,occsh,occij,alamin )
c---------------------------------------------------------------
c Purpose: process some atom information
c
c Written: Peter A. Schultz,  9-May-2000, for v2.45
c
c Revision history:
c  10Jan08-PAS/2.61: eliminate unused passed arguments
c  21Jun01-PAS/2.48: replace STOPs
c   7Jun01-PAS/2.47: eliminate nstd
c---------------------------------------------------------------
c Note: this routine copies the basis info onto the shell info,
c       getting ready for a possible split basis tranformation
c       It also sets up the occupation arrays and finds min
c       alphas for each type to be used later in the code.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input/output arrays:
c  ... basis description
      DIMENSION  numshl(ntyp),lshel(nshld,ntyp),
     $ nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
c  ... reference shell description
      DIMENSION  nocc(ntyp),locc(nshld,ntyp),
     $ nalsh(nshld,ntyp),alsh(nald,nshld,ntyp),calsh(nald,nshld,ntyp)
      DIMENSION  znuc(ntyp), occsh(nshld,ntyp), occij(nshld,nshld,ntyp)
      DIMENSION  alamin(ntyp)
c
c local declarations:
      DATA  zero / 0.d0 /
      DATA  cutocc / 0.9d-8 /
      LOGICAL  booboo
c
c >>>> EXECUTABLE CODE:
c
      booboo = .false.
c
      do 1000 ityp=1,ntyp
        nshl = numshl(ityp)
c
c  Copy basis into shell info:
        nocc(ityp) = numshl(ityp)
        do  ishl=1,nshl
          locc(ishl,ityp) = lshel(ishl,ityp)
          nalsh(ishl,ityp) = nala(ishl,ityp)
          do  ial=1,nalsh(ishl,ityp)
            alsh(ial,ishl,ityp) = ala(ial,ishl,ityp)
            calsh(ial,ishl,ityp) = cala(ial,ishl,ityp)
          enddo
        enddo
c
c  Take care of occupation information:
        occatm = zero
        call MKZERO( nshld*nshld, occij(1,1,ityp) )
        do  ishl=1,nshl
          occshl = occsh(ishl,ityp)
          if( ABS( occshl ).lt.cutocc ) occshl = zero
          occatm = occatm + occshl
          occsh(ishl,ityp) = occshl
          norbsh = 2*lshel(ishl,ityp) + 1
          occij(ishl,ishl,ityp) = occshl / DBLE( norbsh )
        enddo
c
        if( ABS( znuc(ityp)-occatm ) .gt. cutocc )then
          write(IWR,*) '***** ERROR: occups .ne. z on atom type',ityp
          write(IWR,*) '>>>>> z,occatm=',znuc(ityp), occatm
          booboo = .true.
        endif
c
c  Find and store minimum alpha:
        almina = ala(1,1,ityp)
        do  ishl=1,nshl
          lshl = lshel(ishl,ityp)
          nal = nala(ishl,ityp)
          call GNORM( lshl, nal, ala(1,ishl,ityp),cala(1,ishl,ityp) )
          do  ial=1,nal
            almina = MIN( almina, ala(ial,ishl,ityp) )
          enddo
        enddo
        alamin(ityp) = almina
c
 1000 continue
c
      if( booboo ) call STOPXERR( 'bfsetup-errors' )
c
c    That's all Folks!
c
      RETURN
      END
