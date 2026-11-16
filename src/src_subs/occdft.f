c************************************************************************
c OCCDFT routines
c   SETOCCDFT - sets up the occupations
c   OCCDFT0   - set up default aufbau
c************************************************************************
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SETOCCDFT
c
c
      subroutine SETOCCDFT( IDAT,IWR, ioccdft,iflcdft,
     $ nspin, nk, nstate,  elecno,
     $  eigpop )
c---------------------------------------------------------------
c Purpose: read eigvals for a koopmans runo file (for koopmans)
c
c Written: Peter A. Schultz, 6-December-2022, for v2.69 (koopmans/occdft)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( idbgkoop=0 )
c
c Input array:
      DIMENSION  elecno(nspin)
c Output array:
      DIMENSION  eigpop(nstate,nk,nspin)
c Local declarations:
      DIMENSION    noccs(2)
      CHARACTER*6  wantlbl, readlbl, errlbl
      CHARACTER*128 linerd, filenm
      DATA  nline/128/
c
c >>>> EXECUTABLE CODE:
c
  913 format( a, 5i8 )
      if( idbgkoop.ne.0 )then
        write(IWR,913) 'DBG/OCCDFT: nspin,nk,nstate=',nspin,nk,nstate
        write(IWR,*  ) 'DBG/OCCDFT: elecno=',(elecno(isp),isp=1,nspin)
      endif
      call MPNODE( iproc )
      call MPNODE0( master )
c
c Set up formatting for input data:
 8000 format(a128)
  901 format(1x,i5,a )
  910 format(a6)
  911 format('k-vector',i7,7(7x,i7))
  915 format(1x,i5,2x,8f14.8)
  916 format(1x,5x,2x,8f14.8)
      kset = 8
c
c Initialize with standard Aufbau occupations
c
      ioccdft = 0
      iflcdft = 0
      call OCCDFT0( nspin,nk,nstate, noccs, elecno, eigpop )
c
      if( idbgkoop.gt.0 )
     $write(IWR,*) 'DBG:AFT OCCDFT0, noccs=',noccs(1),noccs(2)
c
      read(IDAT,910)  readlbl
      write(IWR,910)  readlbl
      if( readlbl .eq. 'occfil' )then
c       Attempt to open an existing occdft file
        read(IDAT,8000)   linerd
        filenm = ' '
C        write(IDAT,8000)  linerd
        call LNWRIT( IWR, linerd, nline )
        call STRPARS(  linerd, 1,nline, n1,n2 )
        nfile = n2 - n1 + 1
        nfile2 = nfile
        filenm(1:nfile) = linerd(n1:n2)
c
  100 continue
c       Try to open an existing file
C        if( iproc .eq. master )then
          call FLOPEN( iflcdft, filenm, 'OLD', 'FOR', iflerr )
          if( iflerr .eq. 0 )then
c           We found an existing file, use it, ignore rest of occdft input
c           This way, can leave the "excite" setup input in the file.
            write(IWR,*) 'Found occdft file=',filenm(1:nfile2)
c           From occdft file over (presumably) full set of nk k-points
            ioccdft = 2
            goto 700
          else
c           Try again with an assumed '.occdft' suffix
            write(IWR,*) 'Failed to find occfile=',filenm
            if( nfile2 .eq. nfile )then
              nfile1 = nfile + 1
              nfile2 = nfile + 7
              filenm(nfile1:nfile2) = '.occdft'
              goto 100
            endif
            call STOPXERR( 'occdft file missing' )
          endif
C        endif
        read(IDAT,910)  readlbl
        write(IWR,910)  readlbl
      endif
c
      ierr = 0
c
  200 continue
c
 9200 format(1x,2i6,6x,a)
c
      if( readlbl .eq. 'excite' )then
c       spinless --> spinless only
        if( nspin.eq.2 ) call STOPXERR( 'excite only with spinless' )
        ispin = 1
        read(IDAT,*)  iocc,imt 
        write(IWR,9200) iocc,imt ,'excites nocc-i --> nocc+j spinless'
        call OCCEXCITE( ispin,iocc,imt, ierr,
     $   nspin, nk, nstate,  noccs,
     $   eigpop )
c
      elseif( readlbl .eq. 'upexci' )then
c       --> spinfull up-spin excitation
        if( nspin.eq.1 ) call STOPXERR( 'upexcite implies spin' )
        ispin = 1
        read(IDAT,*)  iocc,imt 
        write(IWR,9200) iocc,imt ,'excites nocc+i --> nocc+j spin-up'
        call OCCEXCITE( ispin,iocc,imt, ierr,
     $   nspin, nk, nstate,  noccs,
     $   eigpop )
c
      elseif( readlbl .eq. 'dnexci' )then
c       --> spinfull dn-spin excitation
        if( nspin.eq.1 ) call STOPXERR( 'dnexcite implies spin' )
        ispin = 2
        read(IDAT,*)  iocc,imt 
        write(IWR,9200) iocc,imt ,'excites nocc+i --> nocc+j spin-dn'
        call OCCEXCITE( ispin,iocc,imt, ierr,
     $   nspin, nk, nstate,  noccs,
     $   eigpop )
c
      elseif( readlbl(1:4) .eq. 'save' .or. readlbl .eq. 'end_oc' )then
c       Save this to a file, only want the master doing this
c       Set flag to indicate occdft file has a single k-point
        ioccdft = 1
        if( iproc .eq. master )then
          call FLOPENA( iflcdft, 'occdft' )
          REWIND( unit=iflcdft )
          write( iflcdft, '(a6)' ) 'occups'
          call WREO( iflcdft, nspin,nk,nstate, eigpop )
        endif
c
        if( readlbl .eq. 'end_oc' ) goto 999
c       Ignore the rest of the occdft input after the "save"
        goto 700
c
c
c ***** Set flag that there is an active occdft file
c
      else
        write(IWR,*) '**** ERROR: unknown occdft keyword >',readlbl,'<'
        call STOPXERR( 'unknown occdft input keyword' )
c
      endif
c
      if( ierr.gt.0 )then
c       An error in specifying the excitations
        write(IWR,'(a,i2)') 'ERROR unallowed excitation, ierr=',ierr
        errlbl = 'EXCITE'
        goto 1313
      endif
      read(IDAT,910)  readlbl
      write(IWR,910)  readlbl
c
c     Go parse more instructions
      goto 200
c
  700 continue
      wantlbl = 'end_oc'
      read(IDAT,910)  readlbl
      write(IWR,910)  readlbl
      if( readlbl .ne. wantlbl ) goto 700
      goto 999
c
 1313 continue
c     Error handling
      write(IWR,*) 'ERROR in occdft input, flag=',errlbl
      call STOPXERR( 'occdft input' )
c
  999 continue
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> OCCDFT0
c
c
      subroutine OCCDFT0( nspin, nk, nstate,  noccs, elecno,
     $  eigpop )
c---------------------------------------------------------------
c Purpose: initialize occups to closed aufbau
c
c Written: Peter A. Schultz, 6-December-2022, for v2.69 (koopmans/occdft)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero=0d0, one=1.d0, two=2.d0 )
      PARAMETER  ( idbgkoop=0 )
c
c Input array:
      DIMENSION  elecno(*)
c Output array:
      DIMENSION  eigpop(nstate,nk,*), noccs(*)
c Local declarations:
      DIMENSION    occfull(2)
      DATA         occfull /2.0d0,1.0d0/
      DATA         delocc / 1.d-8 /
      CHARACTER*6  wantlbl, readlbl, errlbl
      CHARACTER*128 linerd, filenm
c     The eigval increment for extra states above base states
c
c >>>> EXECUTABLE CODE:
c
  913 format( a, 5i8 )
      if( idbgkoop.gt.0 )then
        call FLGETIWR( IWR )
        write(IWR,913) 'DBG/OCCDFT0: nspin,nk,nstate=',nspin,nk,nstate
        write(IWR,913) 'DEB/OCCDFT0: elecno=',(elecno(isp),isp=1,nspin)
      endif
c
c  Start off with everything empty
      call MKZERO( nstate*nk*nspin, eigpop )
c
      full = occfull(nspin)
      do 1000 ispin=1,nspin
        nocc = ( (elecno(ispin)+delocc) / full )
        elecsp = nocc*full
        elecn1 = full
        if( (elecsp+delocc) .lt. elecno(ispin) )then
          nocc = nocc + 1
          elecn1 = elecno(ispin) - elecsp
        endif
        noccs(ispin) = nocc
        do  k=1,nk
          do  n=1,nocc
            eigpop(n,k,ispin) = full
          enddo
          eigpop(nocc,k,ispin) = elecn1
        enddo
 1000 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> OCCEXCITE
c
c
      subroutine OCCEXCITE( ispin,iocc,imt, ierr,
     $ nspin, nk, nstate,  noccs,
     $ eigpop )
c---------------------------------------------------------------
c Purpose: installs desired 1e-excitation into occupations
c
c Written: Peter A. Schultz, 6-December-2022, for v2.69 (koopmans/occdft)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c eigpop() records the desired occupations that will be used later
c
c   eigpop had originally been initialized to standard Aufbau,
c    but this routine will operate on a previously excited eigpop,
c    allowing for repeated excitations.
c   Here, the current eigpop() will be updated with an excitation iocc->imt
c    in the occupations of the specified spin=ispin,
c    where iocc and imt are referenced to the nocc(ispin) (top orbital),
c    such that the excitation is from state (nocc+iocc) to state (nocc+imt)
c    E.g. (iocc,,imt) - excitation
c         (0,1) - top occupied state to first empty state (nocc-->nocc+1)
c        (-1,2) - first below top occupied to second empty (nocc-1)->(nocc+2)
c   This routine error-checks ...
c     ... ensure excitation in-bounds of eigpop: 0 < #state <= nstate
c     ... is possible: iocc has an electron, imt has room for an electron
c     ... and returns an error flag if these fail
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero=0d0, one=1.d0, two=2.d0 )
      PARAMETER  ( idbgkoop=0 )
c
c Input array:
      DIMENSION  noccs(*)
c Input/Output array:
      DIMENSION  eigpop(nstate,nk,*)
c Local declarations:
      DIMENSION    occfull(2)
      DATA         occfull /1.0d0,2.0d0/
      DATA         delocc / 1.d-8 /
c
c >>>> EXECUTABLE CODE:
c
      ierr = 0
      call FLGETIWR( IWR )
      if( idbgkoop.gt.0 )then
  913   format( a, 5i8 )
c       Stamp arrival in occexcite
        write(IWR,913) 'DBG/OCCEXCITE: nspin,nk,nstate=',nspin,nk,nstate
        write(IWR,913) 'DBG/OCCEXCITE: isp,nocc,nst, iocc,imt=',
     $                  ispin,noccs(ispin),nstate, iocc,imt
      endif
c
      nocc = noccs(ispin)
      kocc = nocc + iocc
      kmt  = nocc + imt
c Check if these are allowed values, return with error if the are not:
      if( kocc.lt.1  .or. kocc.gt.nstate .or.
     $    kmt .lt.1  .or. kmt .gt.nstate )then
        write(IWR,*) 'ERROR: kocc,kmt,nocc,nstate=',kocc,kmt,nocc,nstate
c       ierr=2: specified occ or empty excited state is out-of-bounds
        ierr = 2
        goto 1313
      endif
c     
      full = occfull(nspin) + delocc
c
      do 200 k=1,nk
        if( (eigpop(kocc,k,ispin) - one + delocc ) .lt. zero .or.
     $      (eigpop(kmt ,k,ispin) + one - delocc ) .gt. full )then
c         ierr=11: Source state lacks an electron, or target state is not open
          ierr = 11
          goto 1313
        endif
        eigpop(kocc,k,ispin) = eigpop(kocc,k,ispin) - one
        eigpop(kmt ,k,ispin) = eigpop(kmt ,k,ispin) + one
  200 continue
      goto 999
c
c ERROR handling
 1313 continue
      write(IWR,*) 'ERROR: in occdft excite specification'
      RETURN
c
c Successful exit
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
