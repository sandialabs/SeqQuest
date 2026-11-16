c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c STRTOOLS
c
c Purpose: Tools for manipulating strings
c Written: Peter A. Schultz, Decemeber-2014, for 2.65
c
c Routines:
c  STRPARS    - extracts the first (bland-delimited) character string
c  LNENDS     - finds index in string of last non-blank character
c  LNWRIT     - writes a character line, out to last non-blank character
c  LNREADTYP  - reads an atom type out of a character strong line
c  LNREADR8   - reads a REAL*8 value out of character string line
c
c Revision history:
c   5Jun15-PAS/2.67: special cases for atom type reads
c   9Apr15-PAS/2.66: LN line routines added into this module (from vdwload)
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRPARS
c
c
      subroutine STRPARS( line, nl1,nl2, n1,n2 )
c---------------------------------------------------------------
c Purpose: parse bounds of first character string within line
c
c Written: Peter A. Schultz, 26-November-2001, for v2.51
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER*(*)  line
c
c local declarations:
      CHARACTER*(1)  blank
      DATA           blank / ' ' /
c
c >>>> EXECUTABLE CODE:
c
      n1 = 0
      n2 = nl2
      do  il=nl1,nl2
        if( line(il:il) .ne. blank )then
c         Found the first non-blank character ...
          n1 = il
          n2 = n1
c          ... and find the last character in string:
          do  jl=n1+1,nl2
            if( line(jl:jl) .eq. blank ) goto 999
            n2 = jl
          enddo
          goto 999
        endif
      enddo
c
  999 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LNENDS
c
c
      subroutine LNENDS( line, nline, n1line, ndline )
c
      IMPLICIT NONE
c
      INTEGER    nline, n1line, ndline
      CHARACTER*(*)  line
c Local declarations:
      INTEGER    ic
c
c Find last non-blank character in line:
      ndline = nline
      do 10 ic=nline,1,-1
        if( line(ic:ic) .ne. ' ' )then
          ndline = ic - 1
          goto 11
        endif
   10 continue
c     Line is blank
      ndline = 0
      n1line = 0
      goto 99
   11 continue
c
c Find first non-blank character in line:
      n1line = 1
      do 20 ic=1,ndline
        if( line(ic:ic) .ne. ' ' )then
          n1line = ic
          goto 21
        endif
   20 continue
   21 continue
c
   99 continue
c     
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LNWRIT
c
c
      subroutine LNWRIT( IWRLN, line, nline )
c---------------------------------------------------------------
c Purpose: write out line, only out to last non-blank char
c Written: Peter A. Schultz, February-2013, for v2.65 (vdw)
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c input
      INTEGER        IWRLN
      CHARACTER*(*)  line
      INTEGER        nline
c Local
      INTEGER        ic, lenlin
      LOGICAL        printon
c
      printon = ( IWRLN .ge. 0 )
c
      if( .not. printon ) RETURN
c
      lenlin = 0
c Determine length of line to write
      do  ic=nline,1,-1
        if( line(ic:ic) .ne. ' ' )then
          lenlin = ic
          goto 120
        endif
      enddo
c
  120 continue
      if( lenlin .gt. 0 )then
c       Write the line out
        write(IWRLN,'(a)')  line(1:lenlin)
      else
c       Write a blank line
        write(IWRLN,*)
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LNREADTYP
c
c
      subroutine LNREADTYP( line,n1,nline, ityp, ntyp,typnm,natmnm )
c---------------------------------------------------------------
c Purpose: determine typ # from input
c Written: Peter A. Schultz, February-2013, for v2.65 (vdw)
c Revised:
c  29Apr16-PAS/2.67b - removed length-limimit on integer read
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c input/Output
      CHARACTER*(*)  line
      INTEGER        n1,nline
      CHARACTER*(*)  typnm
      DIMENSION  typnm(*)
      INTEGER        ityp, ntyp,natmnm
c Local:
c     atmnm needs to be at least as long as typnm (I.e. natmnm)
      CHARACTER*(12)  atmnm
      INTEGER         is1,isn,islast, nmlen, jtyp
c
c >>>> EXECUTABLE CODE:
c
      ityp = -2
      if( n1 .gt. nline ) RETURN
      call STRPARS( line, n1,nline, is1,isn )
c
c     CHeck if this line is blank
      ityp = -1
      if( is1 .le. 0 ) RETURN
      n1 = isn + 1
c
      atmnm = ' '
      islast = isn
      nmlen = islast - is1 + 1
      if( nmlen .gt. natmnm ) nmlen = natmnm
      islast = is1 + nmlen - 1
      atmnm(1:natmnm) = line(is1:islast)
c
C      if( nmlen.lt.3 )then
c       Attempt to read an atom type number out of this:
        read( atmnm, * , err=100 ) jtyp
        if( jtyp.gt.0 .and. jtyp.le.ntyp )then
          ityp = jtyp
          RETURN
        endif
C      endif
c
  100 continue
      do  jtyp=1,ntyp
        if( atmnm(1:natmnm) .eq. typnm(jtyp) )then
          ityp = jtyp
          RETURN
        endif
      enddo
c
c  Did not find the type, check some special cases:
      if( atmnm(1:5) .eq. 'FIRST' )then
        ityp = 1
      elseif( atmnm(1:4) .eq. 'LAST' )then
        ityp = ntyp
      else 
c       We could not identify the atom
        ityp = 0
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LNREADR8
c
c
      subroutine LNREADR8( line,n1,nline, value, ierr )
c---------------------------------------------------------------
c Purpose: determine real*8 from input
c Written: Peter A. Schultz, February-2013, for v2.65 (vdw)
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c input/Output
      CHARACTER*(*)  line
      INTEGER        n1,nline,ierr
      REAL*8         value
c Local
      INTEGER        is1,isn
c
c >>>> EXECUTABLE CODE:
c
      ierr = -2
      if( n1 .gt. nline ) RETURN
      call STRPARS( line, n1,nline, is1,isn )
c
c     CHeck if this line is blank
      ierr = -1
      if( is1 .le. 0 ) RETURN
c
c Read the value from the line, and update the char counter
      read( line(is1:isn), * , err=130 )  value
      n1 = isn + 1
      ierr = 0
c
      RETURN
c
c Error encountered
  130 continue
      value = 0.d0
      ierr = -1
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
