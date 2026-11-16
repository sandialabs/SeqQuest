c************************************************************************
c Small module for i/o of 1e-eigvals and occups for Koopmans/occdft
c 
c  WREIGOCC - controller to write out eigvals and occups into ascii files
c  WREO     - writes the data to a file
c  RDEO     - reads  the data from a file
c  RDEIGS   - reads the eigvalues
c  RDOCCS   - reads the occupations
c  These reads must reconcile the base spin/states/kpts with the new run
c 
c************************************************************************
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WREIGOCC
c
c
      subroutine WREIGOCC( nspin, nk, npop,nocc0k,nstate, wtk,
     $  eigval,eigpop )
c---------------------------------------------------------------
c Purpose: output eigenvalues/occupations into files (for koopmans)
c
c Written: Peter A. Schultz, 5-December-2022, for v2.69
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input arrays:
C      DIMENSION  efermi(2),egap(2)
      DIMENSION  wtk(nk)
      DIMENSION  eigval(nstate,nk,nspin),eigpop(nstate,nk,nspin)
c
c >>>> EXECUTABLE CODE:
c
c >>>>> First, write out the eigenvalues
c
      call FLOPENA( ioccfl, 'eigs' )
c
      REWIND( ioccfl )
      write(ioccfl,'(a6)') 'eigval'
      call WREO( ioccfl, nspin,nk,nstate, 
     $                eigval )
c
      call FLCLOSE( ioccfl )
c
c >>>>> Second, write out the occupations
c
      call FLOPENA( ioccfl, 'occs' )
c
      REWIND( ioccfl )
      write(ioccfl,'(a6)') 'occups'
      call WREO( ioccfl, nspin,nk,nstate, 
     $                eigpop )
c
      call FLCLOSE( ioccfl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WREO
c
c
      subroutine WREO( ieofile, nspin, nk, nstate,
     $  eigocc )
c---------------------------------------------------------------
c Purpose: write eigenvalues/occupations to file (for koopmans)
c
c Written: Peter A. Schultz, 5-December-2022, for v2.69
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input arrays:
C      DIMENSION  efermi(2),egap(2)
      DIMENSION  wtk(nk)
      DIMENSION  eigocc(nstate,nk,nspin)
c
c >>>> EXECUTABLE CODE:
c
c Set up formatting for output data:
  901 format(1x,i5,a )
  910 format(a)
  911 format('k-vector',i7,7(7x,i7))
  915 format(1x,i5,2x,8f14.8)
      kset = 8
c
      write(ieofile,901) nspin ,  ' = nspin'
      write(ieofile,901) nk    ,  ' = kpoints'
      write(ieofile,901) nstate,  ' = nstate'
c
      if( nspin.eq.1 )then
        write(ieofile,910) 'spinless'
      else
        write(ieofile,910) 'spinup'
      endif
c
      do 500 ispin=1,nspin
        if( ispin.eq.2 ) write(ieofile,910) 'spindn'
c
        do 400 k1=1,nk,kset
c         Cycle over sets of k-points (sets of 8)
          kn = MIN( k1-1+kset, nk )
c
c         Write out the k-vector heading
          write(ieofile,911) (k, k=k1,kn)
c
c         Write out the full list of eigvals/occups
          do 300 n=1,nstate
            write(ieofile,915) n, ( eigocc(n,k,ispin), k=k1,kn )
  300     continue
c
  400   continue
  500 continue
c
      write(ieofile,910) 'EOF'
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RDEO
c
c
      subroutine RDEO( ieofile, eolbl, ioccdft, 
     $ mspin,nspin, nk, nstate,
     $ eigocc )
c---------------------------------------------------------------
c Purpose: read eigvals/occups for a koopman/occdft file
c
c Written: Peter A. Schultz, 5-December-2022, for v2.69
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero=0.d0 )
      PARAMETER  (idebug=0)
c
c Input data:
      CHARACTER*6  eolbl
c Output data:
      DIMENSION  eigocc(nstate,nk,nspin)
c Local declarations:
      CHARACTER*6  wantlbl, readlbl, errlbl
c     The eigval increment to mark extra states above baseline state
      DATA  eincr / 10.0 /
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
      if( idebug.gt.0 ) write(IWR,*) 'RDEO: Read eig/occ unit=', ieofile
c
c Set up formatting for input data:
  910 format(a6)
  911 format('k-vector',i7,7(7x,i7))
  915 format(1x,i5,2x,8f14.8)
  916 format(1x,5x,2x,8f14.8)
      kset = 8
c
      REWIND( unit=ieofile )
      errlbl = 'fltype'
C      wantlbl = 'eigval'
      wantlbl = eolbl
      read(ieofile,910,err=1313)  readlbl
      write(IWR   ,910)           readlbl
      if( readlbl .ne. wantlbl )then
        write(IWR,*) 'ERROR: eigval/occdft file has wrong heading'
        write(IWR,*) 'Seeking ',wantlbl,' found ',readlbl
        goto 1313
      endif
c
      read(ieofile,*) mspin
      read(ieofile,*) mk
      read(ieofile,*) mstate
      if( idebug.gt.1 )
     $write(IWR,'(a,3i8)') 'Eigval/occ file input: nspin,nk,nstate=',
     $ mspin,mk,mstate
c
      if( mspin.ne.nspin .and. eolbl.eq. 'occups' )then
c       Spin mis-match not allowed for occups (twinned for eigvals)
        write(IWR,*) 'ERROR: occups file spin .ne. problem spin'
        errlbl = 'sp-occs'
        goto 1313
      endif
c
      if( mk.ne.nk )then
        if( eolbl.eq. 'eigval' )then
c         The eigenvalue file must match k-points
          write(IWR,*) 'ERROR: eigval file k-pts .ne. nk=',mk,nk
          errlbl = 'eigkpt'
          goto 1313
        else
c         Check the occups in more detail
          if( mk.ne.1 .or.
     $      ( mk.eq.1 .and. ioccdft.ne.1 ) )then
            write(IWR,*) 'ERROR: occups file k-pts .ne. nk=',mk,nk
            errlbl = 'occkpt'
            goto 1313
c         else
c           input mk=1 and ioccdft=1 ==> copy k=1 into all nk
          endif
        endif
      endif
c
      nstrd = MIN( nstate, mstate )
      do 500 ispin=1,mspin
        errlbl = 'spinxx'
        read(ieofile,910)  readlbl
        write(IWR   ,910)  readlbl
        if( readlbl(1:4) .ne. 'spin' ) goto 1313
c
        do 400 k1=1,nk,kset
c         Cycle over sets of k-points (sets of 8)
          kn = MIN( k1-1+kset, nk )
c
c         Write out the k-vector heading
         
          errlbl = 'klabel'
          wantlbl = 'k-vect'
          read(ieofile,910)  readlbl
          write(IWR   ,910)  readlbl
          if( readlbl .ne. wantlbl )then
            write(IWR,*) 'ERROR: format error in eig/occ file'
            write(IWR,*) 'Seeking ',wantlbl,' found ',readlbl
            goto 1313
          endif
c
c         Read in the full list of eigvals/occups
          do 300 n=1,nstrd
            read(ieofile,916)       ( eigocc(n,k,ispin), k=k1,kn )
  300     continue
          if( mstate .gt. nstate )then
c           Blank read through the extras in the base state
            do  n=nstate+1,mstate
              read(ieofile,*)
            enddo
          endif
c
  400   continue
c
c       If states in file is fewer than states in job ...
        if( nstate .gt. mstate )then
c         ... Fill in entries for missing states at the tail end:
          if( eolbl.eq. 'eigval' )then
c           Eigval: Arbitrary, fixed offset (higher) from top base state
            do  k=1,mk
              do  n=mstate+1,nstate
                eigocc(n,k,ispin) = eigocc(mstate,k,ispin) + eincr
              enddo
            enddo
          elseif( eolbl.eq. 'occups' )then
c           Occups: higher states must be empty, set to zero
            do  k=1,mk
              do  n=mstate+1,nstate
                eigocc(n,k,ispin) = zero
              enddo
            enddo
          endif
        endif
c
  500 continue
c
      if( nspin.eq.2 .and. mspin.eq.1 )then
c       From spinless base to spin-polarized ... copy up-to-dn
        if( idebug.gt.1 )
     $  write(IWR,*) 'DEV/rdeo: twin eigoccs sp0 into spin (eigval)'
        do  k=1,nk
          do  n=1,nstate
            eigocc(n,k,2) = eigocc(n,k,1)
          enddo
        enddo
      endif
      if( mk.eq.1 .and. nk.gt.1 )then
c       Expand (copy) from k=1 into all the other k-points
        if( idebug.gt.1 )
     $  write(IWR,*) 'DEV/rdeo: expand eigocc k=1 into nk (occups)'
        do  isp=1,nspin
          do  k=2,nk
            call DCOPY( nstate, eigocc(1,1,isp),1, eigocc(1,k,isp),1 )
          enddo
        enddo
      endif
c
      goto 999
c
 1313 continue
c     Error handling
      write(IWR,*) 'ERROR in eigval/occups input file=',errlbl
      call STOPXERR( 'ERROR koopmans/occdft file' )
c
  999 continue
c
c    That's all, Folks!
c
      RETURN
      END
c Customization of Fortran mode in Emacs. Must be in last 3000
characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
