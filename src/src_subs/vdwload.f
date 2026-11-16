c*******************************************************************************
c VDWLOAD - routines to populate VDW data
c
C Routines:
c  VDWNULL  - initialize vdw parameters to a null state
c  VDWCHECK - check if any atom types are missing vdW parameters
c  VDWLOADZ - from deduced atomic Z values, populate vdW ff parameters
c  VDWDATA  - from input file (within setdata), read user input for vdW ff parameters
c
c Revision history:
c   9Apr15-PAS/2.66: LN line/character manipulation routine xfers out to strpars
c*******************************************************************************
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWNULL
c
c
      subroutine VDWNULL( ntyp, c6vdw,r0vdw, c6ijvdw,r0ijvdw )
c---------------------------------------------------------------
c Purpose: Initialize vdw parameters to a null state
c Written: Peter A. Schultz, April-2013, for 2.65
c---------------------------------------------------------------
c
c The vdw parameters are all positive definite.
c A negative value denotes a parameter that is uninitialized
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( vdwna = -1.d0 )
c Output arrays:
      DIMENSION  c6vdw(ntyp), r0vdw(ntyp)
      DIMENSION  c6ijvdw(ntyp,ntyp), r0ijvdw(ntyp,ntyp)
c
c >>>> EXECUTABLE CODE:
c
      do  ityp=1,ntyp
        c6vdw(ityp) = vdwna
        r0vdw(ityp) = vdwna
        do jtyp=1,ntyp
          c6ijvdw(ityp,jtyp) = vdwna
          r0ijvdw(ityp,jtyp) = vdwna
        enddo
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWCHECK
c
c
      subroutine VDWCHECK( IWR, ntyp, typnm, c6vdw,r0vdw, iunkn )
c---------------------------------------------------------------
c Purpose: check for initialized atom types for vdW
c Written: Peter A. Schultz, April-2013, for 2.65
c---------------------------------------------------------------
c
c  ierr returns number of atom types without valid vdW parameters
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero = 0.d0 )
c Input arrays:
      DIMENSION  c6vdw(ntyp), r0vdw(ntyp)
      CHARACTER*(*)  typnm(ntyp)
c
c >>>> EXECTUABLE CODE:
c
      iunkn = 0
      do  ityp=1,ntyp
        if( c6vdw(ityp).lt.zero .or. r0vdw(ityp).lt.zero )then
c         We have an unaccounted-for atom type
          if( IWR.gt.0 )
     $    write(IWR,*) 'VDW param unknown, type=',ityp,'  ',typnm(ityp)
          iunkn = iunkn + 1
        endif
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWLOADZ
c
c
      subroutine VDWLOADZ( IWRF, ntyp, znuc,typnm,natmnm, c6vdw,r0vdw )
c---------------------------------------------------------------
c Purpose: load vdw parameters from internal data, using deduced atomic-Z
c
c Written: Peter A. Schultz, xx-February-2013, for 2.65
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c  Two forms of pair-potential vdW corrections are currently implemented:
c  (1) DFT-D2, per Stefan Grimme, J. Comput. Chem. 27, 1787 (2006).
c  (2) ULG, per H. Kim, J.-M. Choi, and WAG III, JPhChL 3, 360 (2012).
c
c  Note that meanings of d_vdw and rt0vdw are different in these two methods.
c   d_vdw is used in different functional forms for damp short-range interations.
c   r0vdw in DFT-D@ is an atomic vdW radius, in ULG it is a dimer bond minimum.
c   rc_vdw os a global long-range cutoff for evaluation of these terms
c
c  Input vdW parameters are assumed to be Ry and bohr.
c  Output units are also in Ry atomic units
c  
c  Note that I will use overlap cutoffs to limit range of interacions,
c  as the lattice vectors are set up for this (and <10^-6 terms will not matter)
c
c  This version is non-parallel, all processors do whole thing
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON /VDWCTRL / i_vdw, i_dftd2,i_ulg
      COMMON /VDWSTUFF/ vdw_s,vdw_d,vdw_rc
c
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0, half=0.5d0 )
c
c Output arrays:
      DIMENSION  c6vdw(ntyp), r0vdw(ntyp)
c
c Input arrays:
      DIMENSION  znuc(ntyp)
      CHARACTER*(*)  typnm(ntyp)
c
c Local declarations:
      CHARACTER*80  strtmp
      CHARACTER*2   symbty
      LOGICAL       printon
c
c >>>> EXECUTABLE CODE:
c
      printon = ( IWRF .ge. 0 )
      if( i_vdw.eq.0 ) RETURN
c
      if( printon )
     $write(IWRF,*) '#>>>> VDWLOADZ: load internal c6, ivdw=',i_vdw
c
      do 1000 ityp=1,ntyp
c
        if( znuc(ityp).le.zero )then
c         Floating orbital or other potential-less site
          c6vdw(ityp) = zero
          r0vdw(ityp) = zero
          goto 1000
        endif
c
c       Attempt to get element type from type name
        symbty = '  '
        call STRPARS( typnm(ityp), 1,natmnm, n1,n2 )
        if( n1 .eq. 0 ) call STOPXERR( 'nameless type in vdw' )
        strtmp(1:natmnm) = typnm(ityp)
        lennm = MIN( (n2+1-n1), 2 )
        n2 = n1 + lennm - 1
        symbty(1:lennm) =  strtmp(n1:n2)
c
        call ELEMSYMB2Z( iatZ, symbty )
        if( i_vdw .eq. i_dftd2 )then
          call ELEM_DFTD2( iatZ, c6i,r0i )
        elseif( i_vdw .eq. i_ulg )then
          call ELEM_ULG12( iatZ, c6i,r0i )
        else
          if( printon )
     $     write(IWRF,*)'#>>>> VDWLOADZ: c6 not found type=',ityp,symbty
          goto 1000
        endif
c
        c6vdw(ityp) = c6i
        r0vdw(ityp) = r0i
c
        if( printon )     
     $   write(IWRF,'(a,i2,3x,a2,2f16.4)')
     $    '# VDWLOADZ: ityp,c6(Ry),r0(Bohr)=', ityp,symbty,c6i,r0i
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VDWDATA
c
c
      subroutine VDWDATA( IWRF,IDAT, ntyp, typnm,natmnm, c6vdw,r0vdw,
     $                    c6ijvdw,r0ijvdw )
c---------------------------------------------------------------
c Purpose: read user input that modifies the vdW C6 potentials
c
c Written: Peter A. Schultz, xx-February-2013, for 2.65
c
c Revision history:
c  17May14-PAS/2.65c: allow direct input of c6(pair) [w/o r0(pr)]
c   9May14-PAS/2.65c: fix bugs-ijvdw args, and misplaced goto
c---------------------------------------------------------------
c
c  Units are Rydberg atomic units: energy in Rydberg, distance in Bohr.
c
c  This routine designed to allow the entire vdw specification to be specified,
c  or only any part of it that the user wants to modify.
c  With this design, this specification could be located in a separate vdw file.
c   which would require only a little work (in setdata) to implement.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0, half=0.5d0 )
c
c Output arrays:
      DIMENSION  c6vdw(ntyp), r0vdw(ntyp)
      DIMENSION  c6ijvdw(ntyp,ntyp), r0ijvdw(ntyp,ntyp)
c
c Input arrays:
      CHARACTER*(*)  typnm(ntyp)
c
c Local declarations:
      PARAMETER  (nline=80)
      CHARACTER*(nline)  line, tmpline
      CHARACTER*6  wantlbl,readlbl
      CHARACTER*1  ktype
      EQUIVALENCE  (line,readlbl)
      CHARACTER*2   symbty
      LOGICAL      printon
c
c >>>> EXECUTABLE CODE:
c
      printon = ( IWRF.ge.0 )
 8000 format(a80)
 9000 format(1x,a)
 9002 format(1x, 2(5x,f16.8) )
 9012 format(1x, a, 2(5x,f16.8) )
 9022 format(1x, a, 2x, a, 2(5x,f16.8) )
c
  100 continue
      read(IDAT,8000)  line
      if(printon) call LNWRIT( IWRF, line, nline )
c
      if( readlbl .eq. 'vdw_da' )then
c       File heading
        goto 100
      endif
c
      if( readlbl(1:5) .eq. 'vdwff' .or.
     $    readlbl .eq. 'fftype' )then
c       Read the vdw correcion type ...
c       ... and check that it is consistent with earlier input
        read(IDAT,8000)  line
        if(printon) call LNWRIT( IWRF, line, nline )
c       
        call STRPARS( line , 1,nline, n1,n2 )
        if( n1 .eq. 0 )then
           write(IWRF,*) 'no vdw type  label input is blank'
           goto 1300
        endif
        wantlbl = '      '
        nlbl = MIN( n2-n1+1, 6 )
c       Note that cannot use "readlbl" here, as readlbl==line:
        wantlbl(1:nlbl) = line(n1:n2)
c
        call VDWGETNAME( readlbl )
        if( wantlbl .ne. readlbl )then
c         This file and calculation vdw ff disagree - kill
          call FLGETIERR( IERRFL )
          write(IERRFL,*) '**** QM input and vdw file ff mismatch'
          write(IERRFL,*) '>>>> internal  ff=', readlbl
          write(IERRFL,*) '>>>> vdw-input ff=', wantlbl
          call STOPXERR( 'QM and VDW mismatch ff type' )
        endif
      endif
c
c >>>>> Get global factors for vdw pair potentials
c
  300 continue
      read(IDAT,8000)  line
      if(printon) call LNWRIT( IWRF, line, nline )
c
      if( readlbl .eq. 'sfacto' )then
c       Set the s-factor (fcnal-dependent scaling)
        read(IDAT,*) s_vdw
        if(printon)  write(IWRF,9002)  s_vdw
        call VDWSET_S( s_vdw )
        goto 300
      endif
c
      if( readlbl .eq. 'dfacto' .or. 
     $    readlbl(1:3) .eq. 'blg' )then
c       Set the d-factor (damping function factor)
        read(IDAT,*) d_vdw
        if(printon)  write(IWRF,9002)  d_vdw
        if( d_vdw .gt. zero )then
          call VDWSET_D( d_vdw )
        else
c         d-factor must be greater than zero
          call STOPXERR( 'invalid value for d_vdw' )
        endif
        goto 300
      endif
c
      if( readlbl(1:6) .eq. 'rcutof' )then
c       Set the range cutoff for vdW (in bohr)
        read(IDAT,*) rc_vdw
        if(printon)  write(IWRF,9002)  rc_vdw
        if( rc_vdw .gt. zero )then
          call VDWSET_RC( rc_vdw )
        else
c         vdw cutoff distance must be positive
          call STOPXERR( 'invalid value for rc_vdw' )
        endif
        goto 300
      endif
c
      if( readlbl .eq. 'vdw_co' )then
c       Go to the section to get parameters
        goto 400
      endif
c
c     If we lacked it, insert the keyword for coeff section:
      if(printon)  write(IWRF,9000) 'vdw_coefficients'
      goto 401
c
c >>>>> Get single-atom vdw parameters
c
  400 continue
      read(IDAT,8000)  line
      if(printon) call LNWRIT( IWRF, line, nline )
  401 continue
c
      if( readlbl .eq. 'atomc6' .or. readlbl .eq. 'c6atom' )then
c       Read(IDAT,*)     typ  c6
        read(IDAT,8000)  line
        n1 = 1
        call LNREADTYP( line,n1,nline, ityp, ntyp,typnm,natmnm )
        if( ityp.gt.0 .and. ityp.le.ntyp )then
c         Extract c6 coefficient
          call LNREADR8( line,n1,nline, c6i, ierr )
          if( ierr.eq.0 )then
            if(printon) write(IWRF,9012)  typnm(ityp), c6i
            c6vdw(ityp) = c6i
          else
            if(printon) write(IWRF,*) '***** ERROR reading c6 from:'
            if(printon) call LNWRIT( IWRF, line, nline )
            call STOPXERR( 'error reading c6 coeff' )
          endif
        else
          if(printon) call LNWRIT( IWRF, line, nline )
        endif
        goto 400
c
      elseif( readlbl .eq. 'atomr0' .or. readlbl .eq. 'r0atom' )then
C       Read(IDAT,*)     typ  c6
        read(IDAT,8000)  line
        n1 = 1
        call LNREADTYP( line,n1,nline, ityp, ntyp,typnm,natmnm )
        if( ityp.gt.0 .and. ityp.le.ntyp )then
c         Extract r0 coefficient
          call LNREADR8( line,n1,nline, r0i, ierr )
          if( ierr.eq.0 )then
            if(printon) write(IWRF,9012)  typnm(ityp), r0i
            r0vdw(ityp) = r0i
          else
            if(printon) write(IWRF,*) 'ERROR reading r0 from:'
            if(printon) call LNWRIT( IWRF, line, nline )
            call STOPXERR( 'error reading r0 coeff' )
          endif
        else
          if(printon) call LNWRIT( IWRF, line, nline )
        endif
        goto 400
c
      elseif( readlbl .eq. 'atomff' .or. readlbl.eq. 'ffatom' )then
C       Read(IDAT,*)     typ  c6  r0
        read(IDAT,8000)  line
        n1 = 1
        call LNREADTYP( line,n1,nline, ityp, ntyp,typnm,natmnm )
        if( ityp.gt.0 .and. ityp.le.ntyp )then
c         Extract vdw (c6 and r0) coefficients
          call LNREADR8( line,n1,nline, c6i, ierr1 )
          call LNREADR8( line,n1,nline, r0i, ierr2 )
          if( ierr1.eq.0 .and. ierr2.eq.0 )then
            if(printon) write(IWRF,9012)  typnm(ityp), c6i, r0i
            c6vdw(ityp) = c6i
            r0vdw(ityp) = r0i
          else
            if(printon) write(IWRF,*) 'ERROR reading c6 and r0 from:'
            if(printon) call LNWRIT( IWRF, line, nline )
            call STOPXERR( 'error reading vdw coeffs' )
          endif
        else
          if(printon) call LNWRIT( IWRF, line, nline )
        endif
        goto 400
      endif
c
c     End input of atomic c6/r0 parameters, check that dance card is full
      call VDWCHECK( IWR, ntyp, typnm, c6vdw,r0vdw, ierr )
      if( ierr .ne. 0 )then
        write(IWR,*) 'Total atom types lacking vdW parameters=',ierr
        call STOPXERR( 'Some atoms missing vdW parameters' )
      endif
c    
c     Build (default) cross-terms, and check if have 2-center input
      call VDWCROSS( ntyp, c6vdw,r0vdw, c6ijvdw,r0ijvdw )
c
      if( readlbl .ne. 'vdw_2c' ) goto 800
c
c >>>>> Modify atom-pair vdw potentials
c
  500 continue
      read(IDAT,8000)  line
      if(printon) call LNWRIT( IWRF, line, nline )
  501 continue
c
      if( readlbl .eq. 'pairff' .or. readlbl .eq. 'ffpair' )then
c       Read c6/r0 for atom pair
C       read(IDAT)       typ1  typ2  c6ij  r0ij
        read(IDAT,8000)  line
        n1 = 1
        call LNREADTYP( line,n1,nline, ityp, ntyp,typnm,natmnm )
        call LNREADTYP( line,n1,nline, jtyp, ntyp,typnm,natmnm )
        if( ityp.gt.0 .and. ityp.le.ntyp  .and.
     $      jtyp.gt.0 .and. jtyp.le.ntyp )then
c         Valid entries, extract vdw (c6 and r0) coefficients
          call LNREADR8( line,n1,nline, c6ij, ierr1 )
          call LNREADR8( line,n1,nline, r0ij, ierr2 )
          if( ierr1.eq.0 .and. ierr2.eq.0 )then
            if(printon)
     $      write(IWRF,9022)  typnm(ityp), typnm(jtyp), c6ij, r0ij
            c6ijvdw(ityp,jtyp) = c6ij
            r0ijvdw(ityp,jtyp) = r0ij
            c6ijvdw(jtyp,ityp) = c6ij
            r0ijvdw(jtyp,ityp) = r0ij
          else
            if(printon) write(IWRF,*) 'ERROR reading c6 and r0 from:'
            if(printon) call LNWRIT( IWRF, line, nline )
            call STOPXERR( 'error reading vdw coeffs' )
          endif
        else
          if(printon) call LNWRIT( IWRF, line, nline )
        endif
        goto 500
c
      elseif( readlbl .eq. 'pairc6' .or. readlbl .eq. 'c6pair' )then
c       Read c6/r0 for atom pair
C       read(IDAT)       typ1  typ2  c6ij
        read(IDAT,8000)  line
        n1 = 1
        call LNREADTYP( line,n1,nline, ityp, ntyp,typnm,natmnm )
        call LNREADTYP( line,n1,nline, jtyp, ntyp,typnm,natmnm )
        if( ityp.gt.0 .and. ityp.le.ntyp  .and.
     $      jtyp.gt.0 .and. jtyp.le.ntyp )then
c         Valid entries, extract c6 coefficients
          call LNREADR8( line,n1,nline, c6ij, ierr1 )
          if( ierr1.eq.0 )then
            if(printon)
     $      write(IWRF,9022)  typnm(ityp), typnm(jtyp), c6ij
            c6ijvdw(ityp,jtyp) = c6ij
            c6ijvdw(jtyp,ityp) = c6ij
          else
            if(printon) write(IWRF,*) 'ERROR reading c6(pair) from:'
            if(printon) call LNWRIT( IWRF, line, nline )
            call STOPXERR( 'error reading vdw coeffs' )
          endif
        else
          if(printon) call LNWRIT( IWRF, line, nline )
        endif
        goto 500
c
      elseif( readlbl .eq. 'zeroff' .or. readlbl .eq. 'ffzero' )then
c       Read ityp pair for which vdw interaction to be off(zero)
C       Read(IDAT)       typ1  typ2
        read(IDAT,8000)  line
        n1 = 1
        call LNREADTYP( line,n1,nline, ityp, ntyp,typnm,natmnm )
        call LNREADTYP( line,n1,nline, jtyp, ntyp,typnm,natmnm )
        if( ityp.gt.0 .and. ityp.le.ntyp  .and.
     $      jtyp.gt.0 .and. jtyp.le.ntyp )then
          if(printon)
     $    write(IWRF,9022)  typnm(ityp), typnm(jtyp)
          c6ijvdw(ityp,jtyp) = zero
          r0ijvdw(ityp,jtyp) = zero
          c6ijvdw(jtyp,ityp) = zero
          r0ijvdw(jtyp,ityp) = zero
        else
          if(printon) call LNWRIT( IWRF, line, nline )
        endif
        goto 500
c
      endif
c
c     End of 2-center cross-potential input
  800 continue
c
      if( readlbl .eq. 'end_vd' ) goto 900
c
c     We had a problem in the vdw input somewhere
      call FLGETIERR( IERRFL )
      write(IERRFL,*) '**** ERROR in vdw data'
      write(IERRFL,*) '>>>> keyword not recognized:', readlbl
      write(IERRFL,*) '>>>> NB: vdw section must end in: >end_vdwdata<'
 1300 continue
      call STOPXERR( 'error in vdw data input section' )
c
  900 continue
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
