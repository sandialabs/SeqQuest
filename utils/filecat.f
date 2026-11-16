      PROGRAM  filecat
c---------------------------------------------------------------
c Program: filecat
c Purpose: combine multiple files into single sh-separable file
c Written: 23-March-1998, Peter A. Schultz
c Last revised: 15Jun03-PAS, for cleanup
c---------------------------------------------------------------
c
      PARAMETER  (lenflnm=80)
      CHARACTER  infile*(lenflnm), outfil*(lenflnm)
      CHARACTER  blank*1
      DATA  blank / ' ' /
c
      DATA  IRD,IWR, IFIL, ILIB / 5,6, 8, 10 /
c
c >>>> EXECUTABLE CODE:
c
      write(IWR,*) 'Enter output archive name:'
      read(IRD,'(a)')  outfil
      if( outfil.eq.blank .or. outfil.eq.'quit' ) goto 999
c
      OPEN( unit=ILIB, file=outfil, form='formatted', status='new',
     $ err=1301 )
c
      write(IWR,*)
   10 continue
c
      write(IWR,*) 'Enter file to add to archive (<CR>,"quit"=end):'
      read(IRD,'(a)',end=900)  infile
      if( infile.eq.blank .or. infile.eq.'quit' ) goto 900
c
      OPEN( unit=IFIL, file=infile, form='formatted', status='old',
     $ err=1302 )
c
  901 format("PATH=/bin:/usr/bin")
  902 format("cat > ",a," <<'QED_FIL'")
  910 format("QED_FIL")
c
      nchar = lenflnm
      call NUMCHAR( nchar, infile )
      write(ILIB,901)
      write(ILIB,902)  infile(1:nchar)
c
      call COPYFL( IFIL, ILIB, llines )
      write(IWR,*) 'Number of lines copied=',llines
c
      write(ILIB,910)
c
      CLOSE(unit=IFIL)
c
      goto 10
c
 1301 write(IWR,*) 'ERROR: output file already exists'
      goto 999
c
 1302 write(IWR,*)  'ERROR: Could not open file: ',infile
      write(ILIB,910)
      goto 10
c
  900 continue
      CLOSE( unit=ILIB )
      goto 999
c
  999 continue
      STOP 'bye guy'
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COPYFL
c
      subroutine COPYFL( IRD, IWR, llines )
      CHARACTER*128 line
c
      REWIND(unit=IRD)
c
      llines = 0
   10 continue
      read(IRD,'(a128)',end=990,err=990) line
      lenl =      128
      call WRITLN( IWR, lenl, line )
      llines = llines + 1
c
      goto 10
c
  990 continue
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NUMCHAR
c
      subroutine NUMCHAR( nlen, str )
      CHARACTER*(*)  str
      CHARACTER*1  blank
      DATA  blank / ' ' /
c
      do 10 ilen=nlen,1,-1
        if( str(ilen:ilen) .ne. blank )then
          jlen = ilen
          goto 11
        endif
   10 continue
      jlen = 0
   11 continue
c put string length back into length variable
      nlen = jlen
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WRITLN
c
      subroutine WRITLN( IDT, lenl, line )
      CHARACTER*(*)  line
c
      n = lenl
      do 100 i=lenl,1,-1
        if( line(i:i) .eq. ' ' ) goto 100
        n = i
        goto 102
  100 continue
      lenl = 0
      write(IDT,*)
      RETURN
  102 write(IDT,'(a)' ) line(1:n)
c
      RETURN
      END
