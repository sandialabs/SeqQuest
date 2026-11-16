      PROGRAM  mktapes
c---------------------------------------------------------------
c Program: mktapes
c Purpose: generate data deck for td, varying a
c Written: 13-February-1998, Peter A. Schultz
c---------------------------------------------------------------
c
      CHARACTER  filenm*25,headfl*25,tailfl*25,filtyp*20,seq*4
      CHARACTER  ascale*7
      CHARACTER*1  blank
      DATA  blank / ' ' /
      DATA  IRD,IWR, IDT, IHD,IBD,ITL / 5,6, 8, 10,11,12 /
c
      write(*,*) 'Enter sequence label (4 char max):'
      read(*,'(a4)')  seq
c
      write(*,*) 'Enter file type to make:'
      read(*,'(a20)')  filtyp
c
      nseq = 4
      call NUMCHAR( nseq, seq )
      nfil = 20
      call NUMCHAR( nfil, filtyp )
c
      filenm = seq(1:nseq)//'xxxx.tape5'
      if1 = nseq + 1
      ifn = nseq + 4
      headfl = filtyp(1:nfil)//'.head'
      tailfl = filtyp(1:nfil)//'.tail'
c
   20 write(*,*) 'Enter range for a (x.xx-y.yy) xxx,yyy,del:'
      read(*,*)  m1,m2,mdel
      if( m1.lt.100 .or. m1.gt.m2 .or. m2.gt.2000 ) goto 20
      if( mdel.lt.1 ) goto 20
c
      OPEN(unit=IHD,file=headfl,status='old',form='formatted')
      OPEN(unit=ITL,file=tailfl,status='old',form='formatted')
c
      ascale(1:7) = ' xx.xx0'
      do 300 ia=m1,m2,mdel
        write( filenm(if1:ifn), '(i4)' ) ia
        if( ia .lt. 1000 ) filenm(if1:if1) = '0'
        a = ia
        aval = a/100.d0
        write( ascale(2:6), '(f5.2)' ) aval
c
            OPEN(unit=IDT,file=filenm,form='formatted')
c
            call COPYFL( IHD, IDT )
            write(IDT, '(a,a7)' ) 'scale factor, a=',ascale
            write(IDT, '(a7)' ) ascale
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
        
