c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EF2DM
c
c
      subroutine EF2DM( IWR,lvlout,doblas3,npop0k,npopdm,nbox,noboxs )
c---------------------------------------------------------------
c Purpose: determine number of occupied states at which the cost
c          of the density matrix (DM) scheme for generating the
c          density on the grid becomes cheaper the eigenfunction
c          (EF) scheme.
c
c Written: Peter A. Schultz, 22-May-2001, for v2.47
c
c Revision history:
c   8Jan08-PAS/2.61: min/max orb-in-box diagnostic added
c  28Jun02-PAS/2.54: add comment to document, fix output diagnostics
c   3Dec01-PAS/2.51: purge uninitialized variable
c  18May01-PAS/2.47: extract from main routine
c---------------------------------------------------------------
c
c  Summed over boxes, the DM scheme costs SUM[ nobox*(nobox+1)/2 ]
c  and the EF scheme costs SUM[ nobox*numpop ], where nobox is the
c  number of orbitals having non-negligible density in a box, and
c  numpop is the number of occupied levels.
c  The rub: the explicit code takes advantage of the symmetric/
c  Hermitian density matrix to do only triangle of matrix, while
c  the code employing blas3 does a full square because the sym-blas3
c  is slower than full-blas3 (why? who knows). Hence, the crossover
c  depends on whether we have blas3 on or off.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  noboxs(nbox)
      LOGICAL    doblas3
c
c >>>> EXECUTABLE CODE:
c
c Silence output if invalid write unit
      lvlprnt = lvlout
      if( IWR .lt. 0 ) lvlprnt = -1
c
      minobox = noboxs(1)
      maxobox = noboxs(1)
      nsum = 0
      n2sum = 0
      n3sum = 0
      do  ibox=1,nbox
        nobox = noboxs(ibox)
        if( nobox .lt. minobox ) minobox = nobox
        if( nobox .gt. maxobox ) maxobox = nobox
        if( nobox.gt.0 )then
          nsum = nsum + nobox
c         In the following arithmetic, I remove nobox*npop0k
c         from the sum in order to keep n?sum small so that
c         this arithmetic does not overflow the integer variable.
c         I add SUM[nobox*npop0k] = SUM[nobox]*npop0k later
c         when I add to npopdm: npop0k = SUM[nobox*npop0k]/nsum
          n2sum = n2sum + (nobox*(nobox+1-2*npop0k)) / 2
          n3sum = n3sum + nobox*(nobox-npop0k)
        endif
      enddo
c
      if( doblas3 )then
c       Using full rectangle in dm density routines with blas3
        npopdm = npop0k + (n3sum+nsum-1)/nsum
      else
c       Using triangle only in dm compiled density routines
        npopdm = npop0k + (n2sum+nsum-1)/nsum
      endif
c
      if( lvlprnt.gt.1 )then
        write(IWR,*) 'Turnover EF->DM npopdm, npop0k=',npopdm,npop0k
        if( lvlout.gt.2 )then
          write(IWR,'(1x,a,2i7)')  'EF2DM: grid boxes(nbox)=', nbox
          aveobox = DBLE( nsum ) / DBLE( nbox )
          write(IWR,'(1x,a,f8.2,2(1x,i6))') 'Ave/min/max orbitals/box=',
     $     aveobox, minobox, maxobox
        endif
        if( lvlprnt.gt.4 )then
          if( doblas3 )then
            write(IWR,'(1x,a,2i12)') 'EF: n3sum,nsum=',n3sum,nsum
          else
            write(IWR,'(1x,a,2i12)') 'EF: n2sum,nsum=',n2sum,nsum
          endif
        endif
      endif
c
c    That's all Folks!
c
      RETURN
      END
