      PROGRAM  mktapes
c---------------------------------------------------------------
c Program: mktapes
c Purpose: generate data deck for hcp, varying a
c Written: 13-February-1998, Peter A. Schultz
c---------------------------------------------------------------
c
      CHARACTER  filenm*25,headfl*25,tailfl*25,filtyp*20,seq*3
      CHARACTER  ascale*6,cscale*6
      CHARACTER*1  blank
      DATA  blank / ' ' /
      DATA  IRD,IWR, IDT, IHD,IBD,ITL / 5,6, 8, 10,11,12 /
c
C      write(*,*) 'Enter sequence label (3 char max):'
C      read(*,'(a3)')  seq
C      nseq = 3
C      call NUMCHAR( nseq, seq )
c
      write(*,*) 'Enter file type to make:'
      read(*,'(a20)')  filtyp
c
      nfil = 20
      call NUMCHAR( nfil, filtyp )
c
      filenm = 'axxxxcxxxx.tape5'
C      if1 = nseq + 1
C      ifn = nseq + 4
      headfl = filtyp(1:nfil)//'.head'
      tailfl = filtyp(1:nfil)//'.tail'
c
   23 write(*,*) 'Enter range for c (x.xxx-y.yyy) xxxx,yyyy,del:'
      read(*,*)  l1,l2,ldel
      if( l1.lt.1000 .or. l1.gt.l2 .or. l2.gt.2000 ) goto 23
      if( ldel.lt.1 ) goto 23
c
   21 write(*,*) 'Enter range for a (x.xxx-y.yyy) xxxx,yyyy,del:'
      read(*,*)  m1,m2,mdel
      if( m1.lt.1000 .or. m1.gt.m2 .or. m2.ge.10000 ) goto 21
      if( mdel.lt.1 ) goto 21
c
      OPEN(unit=IHD,file=headfl,status='old',form='formatted')
      OPEN(unit=ITL,file=tailfl,status='old',form='formatted')
c
      cscale(1:6) = ' x.xxx'
      do 300 ic=l1,l2,ldel
        write( filenm(7:10), '(i4)' ) ic
        c = ic
        cval = c / 1000.d0
        write( cscale(2:6), '(f5.3)' ) cval
c
        ascale(1:6) = ' x.xxx'
        do 200 ia=m1,m2,mdel
          write( filenm(2:5), '(i4)' ) ia
          a = ia
          aval = a / 1000.d0
          write( ascale(2:6), '(f5.3)' ) aval
c
            OPEN(unit=IDT,file=filenm,form='formatted')
c
            call COPYFL( IHD, IDT )
            write(IDT, '(a,a6)' ) 'scale factor, a=',ascale
            write(IDT, '(a6)' ) ascale
            write(IDT, '(a,a6)' ) 'scalez factor, c/a=',cscale
            write(IDT, '(a6)' ) cscale
            call COPYFL( ITL, IDT )
c
            CLOSE(unit=IDT)
c
  100     continue
c
  200   continue
c
  300 continue
c
      STOP 'bye guy'
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COPYFL
c
      subroutine COPYFL( IRD, IWR )
      CHARACTER*80 line
c
      REWIND(unit=IRD)
c
   10 continue
      read(IRD,'(a80)',end=990,err=990) line
      lenl = 80
      call WRITLN( IWR, lenl, line )
      if( lenl.lt.1 ) goto 990
c
      goto 10
c
  990 RETURN
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
      RETURN
  102 write(IDT,'(a)' ) line(1:n)
c
      RETURN
      END
        
