c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EVPRNT
c
c
      subroutine EVPRNT( IWR, lprnt, nspin, efermi,egap, npop,nocc0k,
     $ nstate,nk, wtk,eigval,eigpop )
c---------------------------------------------------------------
c Purpose: output computed levels and occupations
c
c Written: Peter A. Schultz, 9-February-1999, for v2.31
c
c Revision history:
c   8Apr15-PAS/2.66: snipped printout of occupied states
c  12Nov08-PAS/2.62: check nstate vs. npop for sufficiency
c   8Mar02-PAS/2.52: clean unused variables
c  19Dec01-PAS/2.51: side-by-side spin up/dn output
c  20Jul01-PAS/2.49: spin-polarized dft
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input arrays:
      DIMENSION  efermi(2),egap(2)
      DIMENSION  wtk(nk)
      DIMENSION  eigval(nstate,nk,nspin),eigpop(nstate,nk,nspin)
c
c Local declarations:
      CHARACTER*7  spinlbl, spinless,spinup,spindn
      DATA  spinless /'       '/
      DATA  spinup /'spin-up'/, spindn /'spin-dn'/
      DIMENSION  tlvl(4),tpop(4)
      CHARACTER*(14) hdglbl(3)
      CHARACTER      headg*80,elips*21
c
      DIMENSION  occlvl(2)
      DATA       occlvl / 2.d0,1.d0 /
c
      DATA  toev / 13.605805 d0 /
      DATA  nloprnt, nhiprnt, nmtprnt / 4,20,8 /
c
c >>>> EXECUTABLE CODE:
c
c Set spin-type label for eigenspectrum
      occfull = occlvl(nspin)
      if( nspin.eq.2 )then
        spinlbl = spinup
      else
        spinlbl = spinless
      endif
      hdglbl(1) =        'level/DN/occup'
      hdglbl(2) =        'level    occup'
      hdglbl(3) =        'level/UP/occup'
c
      mxset = 4
  911 format(/1x,'k-vector:', i8,3i21 )
  912 format( 1x,' weight:',3x,f9.6,3(12x,f9.6) )
  914 format( 9x,      a)
  915 format( i5,      4(f13.6,f8.4) )
  916 format( 4x, a, a, a, a )
c
  921 format(/1x,'k-vector:',i19,i42)
  922 format( 1x,' weight:',14x,f9.6,  33x,f9.6 )
      elips = '         ...      ...'
c
c >>>>> Write out fermi level/gap info:
c
      write(IWR,*)
      do 100 ispin=1,nspin
        evgap = egap(ispin)*toev
c
        write(IWR,910)  spinlbl, efermi(ispin),egap(ispin),evgap
  910   format(1x,a,' fermi level =',f12.6,' Ry,',
     $   3x,'gap(0K)=',f9.6,' Ry =',f6.3,' eV' )
c
c       Switch to dn-spin label
        spinlbl = spindn
  100 continue
c
c Check for sufficiency of nstate vs. npop
c
      nufstate = 0
      npopmx = npop
      if( npopmx .ge. nstate )then
c       We might have a problem here.
        nufstate = 1
c
        popzero = 0.00004
        do  ispin=1,nspin
          do  k=1,nk
            if( eigpop(nstate,k,ispin) .gt. popzero ) nufstate = 2
          enddo
        enddo
      endif
c
      if( lprnt.lt.2 .and. nufstate.lt.2 ) goto 900
c
c  Determine is skipping over occupied in print, and how many
        nskip = nocc0k - nloprnt - nhiprnt
        if( nskip .lt. 8 )  nskip = 0
        kset = mxset/nspin
c
c >>>>> Write out levels and occupations:
c
        mtprnt0 = nmtprnt
        if( lprnt.gt.2 ) mtprnt0 = nstate
c
        icn =  -4
        headg = ' '
        jspin = 0
        if( nspin.eq.2 ) jspin = 1
        do  k=1,4
          ic1 = icn + 8
          icn = icn + 21
          mspin = 2 + jspin
          headg(ic1:icn) = hdglbl(mspin)
          jspin = -jspin
        enddo
c
c Set up partial prints of occupied states
        do 300 k1=1,nk,kset
          kn = MIN( k1-1+kset, nk )
          k0 = k1 - 1
          mk = kn - k1 + 1
c
          nksp = mk*nspin
          if( nspin.eq.1 )then
            write(IWR,911) ( k,      k=k1,kn )
            write(IWR,912) ( wtk(k), k=k1,kn )
          else
            write(IWR,921) ( k,      k=k1,kn )
            write(IWR,922) ( wtk(k), k=k1,kn )
          endif
          nc = 21*nksp - 4
          write(IWR,914)  headg(1:nc)
          write(IWR,*) ' '
c
c         Limit number of empty orbitals that are printed out:
          iskip0 = 1
          mtprnt = mtprnt0
          partocc = occfull - 0.00002
          do 200 n=1,nstate
            noccs = 0
            iskip = iskip0
            iksp = 0
            do  k=k1,kn
              do  ispin=1,nspin
                iksp = iksp + 1
                tlvl(iksp) = eigval(n,k,ispin)
                tpop(iksp) = eigpop(n,k,ispin)
                if( tpop(iksp) .gt. 0.000001 ) noccs = noccs + 1
                if( tpop(iksp) .lt. partocc ) iskip = 0
              enddo
            enddo
            iskip0 = iskip
c
c           Determine if we skip this occupied state
            if( nskip .gt. 0 .and. n.gt.nloprnt .and. n.lt.nocc0k )then
              if( iskip.eq.1 .and. n.eq.(nloprnt+1) )then
c               Our first skipped occupied level, write out our skip
                write(IWR,916) (elips,iksp=1,nksp)
                goto 200
              elseif( iskip.eq.1 .and. n .le. (nocc0k-nhiprnt) )then
                goto 200
              endif
            endif
c
            if( noccs.eq.0 ) mtprnt = mtprnt - 1
            if( mtprnt.lt.0 )then
              if( n.lt.nstate )then
c               Need to put out highest computed (empty) level ...
                write(IWR,916) (elips,iksp=1,nksp)
                write(IWR,915) nstate,
     $           (( eigval(nstate,k,ispin), eigpop(nstate,k,ispin),
     $             ispin=1,nspin ), k=k1,kn )
              endif
c              ... and jump to new set of k-points:
              goto 300
            endif
            write(IWR,915) n,( tlvl(iksp),tpop(iksp), iksp=1,nksp )
  200     continue
c
  300   continue
c
  900 continue
c
      if( nufstate .eq. 1 )then
        write(IWR,*) '***** CAUTION: number of solved states=',nstate
        write(IWR,*) '>>>>> This is only marginally sufficient'
      elseif( nufstate .eq. 2 )then
        write(IWR,*) '***** ERROR: number of solved states=',nstate
        write(IWR,*) '>>>>> Not enough to catch all populated states'
        write(IWR,*) '>>>>> Specify larger number of states',
     $              ' using >states< keyword in run phase input'
        call STOPXERR( 'Insufficient states for populations' )
      endif
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
c End
