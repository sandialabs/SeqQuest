      PROGRAM MKEMIN
c-----------------------------------------------------------
c Program: EMIN
c Written: Peter A. Schultz
c Revised: 12Jun96-PAS/graceful multi-pt/multi-del analysis
c-----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  e(50)
      CHARACTER  line*100
      DATA  zero / 0.d0 /
c
   10 continue
      write(*,*) 'Enter delta-r:'
      read(*,*) dr
      if( dr.le. 0.d0 ) goto 999
c
   20 write(*,*) 'Enter r1,e1,e2,e3,e4,... :'
   21 read(*,'(a100)')  line
      nval = 0
      l1 = 1
      l2 = 100
c
      call NEXTNUM( r1, line, l1,l2 )
      if( l1.eq.0 ) goto 21
      if( r1.le.zero )then
        if( nval.gt.0 )then
          goto 40
        elseif( r1.eq.zero )then
          goto 10
        elseif( r1.lt.zero )then
          goto 999
        endif
      endif
c
      nval = 0
   30 call NEXTNUM( enval, line, l1,l2 )
      if( l1.gt.0 )then
c       Got a valid pt.
        nval = nval + 1
        e(nval) = enval
        goto 30
      elseif( l1.eq.0 )then
c       Continuation
        read(*,'(a100)')  line
        l1 = 1
        l2 = 100
        goto 30
      endif
c
   40 mxdel = (nval-1)/2
      if( mxdel.lt.1 )then
c       Not enough pts to find min, get new interval
        goto 20
      elseif( mxdel.eq.1 )then
c       Only 3 or 4 pts, 1 int fit only
        ndel = 1
        goto 61
      else
        goto 60
      endif
c
   59 write(*,*) 'Max jump intervals =',(nval-1)/2
   60 write(*,*) 'Enter number of intervals to jump:'
      read(*,*)  ndel
      if( ndel.lt.0 ) goto 999
      if( ndel.eq.0 ) goto 20
      if( ndel.gt.mxdel ) goto 59
c
   61 npt = nval - 2*ndel
      rmid = r1 + (ndel-1)*dr
      rdel = ndel*dr
      do 100 i1=1,npt
        rmid = rmid + dr
        e1 = e(i1)
        e2 = e(i1+ndel)
        e3 = e(i1+2*ndel)
        call MKMIN( rmin,emin, rmid,rdel,e1,e2,e3 )
        write(*,910)  rmid,rmin, rmin, emin
  910   format(1x,'About r =',f7.4,',  rfit = ',f8.5,' = ',f7.4,
     $   '     efit=',f14.6)
  100 continue
      if( mxdel.eq.1 ) goto 20
      goto 60
c
  999 STOP 'bye guy'
      END
c
c
c
      SUBROUTINE MKMIN( rmin,emin, r2,rdel, e1,e2,e3 )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      DATA  two / 2.d0 /
c
      dd = two*e2 - e3 - e1
      d31 = e3 - e1
c
      roff = rdel*(d31/dd)/two
      rmin = r2 + roff
c
      eoff = ((d31/two)**2)/(two*dd)
      emin = e2 + eoff
c
      RETURN
      END
c
      SUBROUTINE NEXTNUM( val, str, l1,l2 )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      CHARACTER  str*(*)
      CHARACTER  str2*12, blank*1, xtnd*1
      DATA  blank / ' ' /, xtnd / '\\' /
      str2 = ' '
c go through leading blanks
      do 10 ll=l1,l2
        lead = ll
        if( str(ll:ll).ne.blank )then
          if( str(ll:ll).eq.xtnd )then
            l1 = 0
            goto 99
          else
            goto 20
          endif
        endif
   10 continue
      l1 = -1
      goto 99
c
c find trailing blank
   20 lend = lead + 1
      do 40 ll=lead+1,l2
        lend = ll
        if( str(ll:ll).eq.blank ) goto 50
   40 continue
c I require a trailing blank to guarantee end of input
      l1 = -1
      goto 99
c
c load and read string
   50 len = lend - lead + 1
      str2(1:len) = str(lead:lend)
      read(str2,*) val
      l1 = lend + 1
c
   99 RETURN
      END
