c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATFMT
c
c
      subroutine ATFMTIN( iatmfmt, line )
c---------------------------------------------------------------
c Purpose: connect type of atom-line format to atom-line label
c
c Written: Peter A. Schultz, 11-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c iatmfmt = format of atom description
c         = 0  "coordinates"
c         = 1  "atom, type, positions"
c
c Warning: routine assumes (does not check) that input "line"
c is at least six characters long, and that the output "line"
c has enough space to load the identifying label.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER*(*)  line
      CHARACTER*(*)  readlbl*6
c
c >>>> EXECUTABLE CODE:
c
      iatmfmt = -1
c
      readlbl = line(1:6)
      if( readlbl.eq. 'coordi' )then
c       This will be position-only input (x,y,z)
        iatmfmt = 0
      elseif( readlbl.eq. 'atom, ' )then
c       This denotes atom,type,position (at,ty,x,y,z) input
        iatmfmt = 2
      endif
      RETURN
c
      entry ATFMTOUT( iatmfmt, line )
c
      if( iatmfmt.eq.0 )then
        line = 'coordinates (x y z) for atoms'
      elseif( iatmfmt.eq.1 .or. iatmfmt.eq.2 )then
        line = 'atom, type, position (at,ty,x,y,z) for atoms'
      else
        call STOPXERR( 'ATFMTOUT: unknown atom output format' )
      endif
      RETURN
c
c >>>> That's all, Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CONFIGRD
c
c
      subroutine CONFIGRD( IDT,IWR, iatmfmt, natmnm, ntyp,typnm,
     $ natm, atmnm, itypa, ratm )
c---------------------------------------------------------------
c Purpose: input atomic configuration
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  10Jan02-PAS/2.54: bare coord read enabled
c---------------------------------------------------------------
c
c iatmfmt = format of atom description
c         = 0  just x,y,z
c         = 1  atom, type(as number), coords
c         = 2  atom name, type(as name), coords
c
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input: (type names)
      CHARACTER*(*)  typnm
      DIMENSION  typnm(*)
c output: atom names, atom types, atom coordinates
      CHARACTER*(*)  atmnm
      DIMENSION  atmnm(*)
      DIMENSION  itypa(*)
      DIMENSION  ratm(3,natm)
c
c >>>> EXECUTABLE CODE:
c
      do  iatm=1,natm
        call ATOMRD( IDT,IWR, iatmfmt, natmnm, ntyp, typnm,
     $   atmnm(iatm), iatm,itypa, ratm(1,iatm) )
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CONFIGWR
c
c
      entry CONFIGWR( IDT, iatmfmt, natmnm, ntyp,typnm,
     $ natm, atmnm, itypa, ratm )
c---------------------------------------------------------------
c Purpose: output atomic configuration
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      do  iatm=1,natm
        call ATOMWR( IDT, iatmfmt, natmnm, ntyp, typnm,
     $   atmnm(iatm), iatm,itypa, ratm(1,iatm) )
      enddo
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COORDRD
c
c
      subroutine COORDRD( IDT, iatmfmt, natm,ratm )
c---------------------------------------------------------------
c Purpose: extract atomic coordinates from geometry file
c
c Written: Peter A. Schultz, 10-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c atomic coordinates:
      DIMENSION  ratm(3,natm)
c Local
      DIMENSION  ratmi(3)
c
c >>>> EXECUTABLE CODE:
c
      jatmfmt = iatmfmt
      do  iatm=1,natm
        call RATMRD( IDT, jatmfmt, ratmi )
        ratm(1,iatm) = ratmi(1)
        ratm(2,iatm) = ratmi(2)
        ratm(3,iatm) = ratmi(3)
      enddo
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COORDWR
c
c
      entry COORDWR( IDT,  natm,ratm )
c---------------------------------------------------------------
c Purpose: dump raw atomic coordinates to geometry file
c
c Written: Peter A. Schultz, 10-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      do  iatm=1,natm
        ratmi(1) = ratm(1,iatm)
        ratmi(2) = ratm(2,iatm)
        ratmi(3) = ratm(3,iatm)
        call RATMWR( IDT, ratmi )
      enddo
c
      RETURN
c
c    That's all Folks!
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATOMRD
c
c
      subroutine ATOMRD( IDT,IWR, iatmfmt, natmnm, ntyp,typnm,
     $ atmnm, iatm,itypa, ratm )
c---------------------------------------------------------------
c Purpose: input atom coordinates
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Oh, the ineffable joys of character manipulation in fortran ...
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input: (type names)
      CHARACTER*(*)  typnm
      DIMENSION  typnm(*)
c output: atom names, atom types, atom coordinates
      CHARACTER*(*)  atmnm
      DIMENSION  itypa(*)
      DIMENSION  ratm(3)
c
c local declarations:
      PARAMETER  ( nline=80 )
      CHARACTER  line*(nline), strtmp*(nline)
      LOGICAL    printon
c
c >>>> EXECUTABLE CODE:
c
      printon = ( IWR .ge. 0 )
c Global format for output of atomic coordinates:
 9100 format(1x,3f16.10)
c
      if( iatmfmt .eq. 0 )then
c       Read the raw coordinates
        read(IDT,*)  ratm
        if( printon ) write(IWR,9100)  ratm
c
      elseif( iatmfmt .eq. 1 .or. iatmfmt .eq. 2 )then
c
c       Read the line with the atomic information:
c
        read(IDT,'(a80)',end=1301,err=1301 )  line
        if( printon ) write(IWR,'(a)' )   line
        nl1 = 1
        nl2 = nline
c
c       Parse the atom number/name (first string on line):
c
        call STRPARS( line, nl1,nl2, n1,n2 )
        if( n1 .eq. 0 ) call STOPXERR( 'no atom name found' )
        nmlen = n2 - n1 + 1
        if( nmlen .gt. natmnm ) call STOPXERR( 'atom name too long' )
        strtmp = ' '
        strtmp(1:nmlen) = line(n1:n2)
        atmnm = strtmp(1:natmnm)
        nl1 = n2 + 1
c
c       First, check if atoms input as number sequence.
c       If they are, check that atoms are indeed in sequence.
        read( atmnm, * , err=120 ) jatm
        if( jatm .ne. iatm ) call STOPXERR( 'atom out of sequence' )
c       Shift atom name to right aligned 6-char deep
        atmnm = ' '
        write( atmnm, '(i6)' ) iatm
c
c       Parse the atom type:
c
  120   continue
        call STRPARS( line, nl1,nl2, n1,n2 )
        if( n1 .eq. 0 ) call STOPXERR( 'no atom type found on line' )
        nmlen = n2 - n1 + 1
        if( nmlen .gt. natmnm ) call STOPXERR( 'type name too long' )
        nl1 = n2 + 1
c
c       And now try and get the atom type:
c
        if( iatmfmt .eq. 2 )then
          itypa(iatm) = 0
          do  ityp=1,ntyp
            strtmp = typnm(ityp)
            call STRPARS( strtmp, 1,nline, nt1,nt2 )
            ntlen = nt2 - nt1 + 1
            if( ntlen .eq. nmlen .and.
     $          line(n1:n2) .eq. strtmp(nt1:nt2) )then
              itypa(iatm) = ityp
              goto 200
            endif
          enddo
        endif
c       Did not find type in list of type names, perhaps as type #?
        read( line(n1:n2), * , err=130 ) ityp
        if( ityp.lt.1 .or. ityp.gt.ntyp ) call STOPXERR( 'bad type' )
        itypa(iatm) = ityp
c
c       Atom type are input as numbers, not names:
        iatmfmt = 1
c
        goto 200
c
  130   call STOPXERR( 'atom type unknown' )
c
  200   continue
c
c       Get the atom coordinates:
        read( line(nl1:nl2), * )  ratm
c
      else
c       Unknown geometry format
        call STOPXERR( 'unknown geometry (input) format' )
      endif
c
      RETURN
Cc
Cc >>>> ERROR HANDLING
Cc
C1301 call STOPXERR( 'error (or eof) reading atom coordinates' )
C     RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATOMWR
c
c
       entry ATOMWR( IDT, iatmfmt, natmnm, ntyp,typnm,
     $ atmnm, iatm,itypa, ratm )
c---------------------------------------------------------------
c Purpose: output atom coordinates
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( iatmfmt .eq. 0 )then
c       Write out raw coordinates without atomic id's
        write(IDT,9100)  ratm
c
      elseif( iatmfmt .eq. 1 .or. iatmfmt .eq. 2 )then
c       Build line with the atomic information:
        line = ' '
        nl1 = 2
        nl2 = nline
c
c       Parse and incorporate the atom name:
        call STRPARS( atmnm, 1,natmnm, n1,n2 )
        if( n1 .eq. 0 ) call STOPXERR( 'nameless atom in atom output' )
        nmlen = n2
        nl2 = nl1 + nmlen - 1
        line(nl1:nl2) = atmnm(1:n2)
        nl1 = nl1 + natmnm + 1
c
c       Parse and incorporate the atom type name:
        ityp = itypa(iatm)
        strtmp = ' '
        if( iatmfmt .eq. 1 )then
c         Output type as type #, not type name
          n2 = 6
          write( strtmp(1:6), '(i6)' ) ityp
        else
c         Output type as type name
          call STRPARS( typnm(ityp), 1,natmnm, n1,n2 )
          if( n1 .eq. 0 ) call STOPXERR( 'nameless type in output' )
          strtmp(1:natmnm) = typnm(ityp)
        endif
        nmlen = n2
        nl2 = nl1 + nmlen - 1
        line(nl1:nl2) = strtmp(1:n2)
        nl1 = nl1 + natmnm + 1
c
c       Incorporate the atom coordinates:
        nl2 = nl1 + 1 + 3*16
        if( nl2 .gt. nline ) call STOPXERR( 'atom type/coord too long' )
        write( line(nl1:nl2), 9100 )  ratm
c
        write(IDT,'(a)')  line(1:nl2)
c
      else
c       Unknown geometry format
        call STOPXERR( 'unknown geometry (output) format' )
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RATMRD
c
c
      entry RATMRD( IDT, iatmfmt, ratm )
c---------------------------------------------------------------
c Purpose: output atom coordinates
c
c Written: Peter A. Schultz, 10-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      if( iatmfmt .eq. 0 )then
c       Read the raw coordinates
        read(IDT,*)  ratm
c
      elseif( iatmfmt .eq. 1 .or. iatmfmt .eq. 2 )then
c
c       Read the line with the atomic information:
c
        read(IDT,'(a80)',end=1301,err=1301 )  line
        nl1 = 1
        nl2 = nline
c
c       Skip the atom number/name (first string on line):
c
        call STRPARS( line, nl1,nl2, n1,n2 )
        if( n1 .eq. 0 ) call STOPXERR( 'no atom name found' )
        nl1 = n2 + 1
c
c       Skip the atom type (second string on line):
c
        call STRPARS( line, nl1,nl2, n1,n2 )
        if( n1 .eq. 0 ) call STOPXERR( 'no atom type found on line' )
        nl1 = n2 + 1
c
c       Get the atom coordinates:
        read( line(nl1:nl2), * )  ratm
c
      else
c       Unknown geometry format
        call STOPXERR( 'unknown geometry (input) format' )
c
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RATMWR
c
c
      entry RATMWR( IDT, ratm )
c---------------------------------------------------------------
c Purpose: output atom coordinates
c
c Written: Peter A. Schultz, 10-January-2002, for v2.54
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      write(IDT,9100)  ratm
c
      RETURN
c
c >>>> ERROR HANDLING
c
 1301 call STOPXERR( 'error (or eof) reading atom coordinates' )
      RETURN
c
c    That's all Folks!
c
      END
