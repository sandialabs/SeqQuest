c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATOMNM
c
c
      subroutine ATOMNM( atmfil,nflnm, filenm, atmnm,natnm )
c---------------------------------------------------------------
c Purpose: parse atom file input to get atom name AND file name
c
c Written: Peter A. Schultz, 31-May-2001 for v2.47
c          Based on an idea by Art Edwards
c
c Revision history:
c  R21Jun01-PAS/2.48: replace STOPs
c---------------------------------------------------------------
c
c Looking to parse an atom name and a file name out of an input
c file type name.  Two types of input:
c  (1) <atmnm> = <filenm>
c     - where aliasing is explicit, and the presence of '='
c       uniquely signifies use of an alias.  We allow no spaces
c       within <atmnm> or <filenm>, but will allow spaces on
c       either side of the '='
c  (2) <filenm>
c     - where we will attempt to parse out a default atom name
c       from the atom file name, the atom name being XXXX, where
c       the atom file is XXXX.atm.  Must ignore directory name.
c I take care to be sure that the length of the character string
c for the atom name does not exceed that alloted in 'atmnm'
c
c input:
      CHARACTER*(*)  atmfil
c output:
      CHARACTER*(*)  filenm
      CHARACTER*(*)  atmnm
c
c local declarations:
      LOGICAL      got_name
      CHARACTER*4  atmtag
      DATA         atmtag /'.atm'/
      CHARACTER*1  blank,       equals,       slash
      DATA         blank /' '/, equals /'='/, slash /'/'/
c
c >>>> EXECUTABLE CODE:
c
      got_name = .false.
c Set default atom name to no-name:
      atmnm = blank
      filenm = blank
c
      ic1 = 1
c
   10 continue
c
c first step: find first non-blank character
c
      do 20 ic=ic1,nflnm
        i1 = ic
        if( atmfil(ic:ic) .ne. blank ) goto 21
   20 continue
      call STOPXERR( 'blankatm/atom name or file parsed blank' )
   21 continue
c
c next step: find end of string (blank), and look for '=' which
c            would signify an alias.
c
      i2 = nflnm
      do 50 ic=i1+1,nflnm
c
        if( atmfil(ic:ic) .eq. blank )then
          i2 = ic - 1
c         If already have an atom name, this string is file name:
          if( got_name ) goto 100
c
c         Look for first non-blank character and check for '='
c         that would signify an alias:
          do 40 jc=ic+1,nflnm
c           End string, we either have atom name or file name:
            if( atmfil(jc:jc) .ne. blank ) then
              if( atmfil(jc:jc) .eq. equals )then
c               Found an alias
                got_name = .true.
                inm1 = i1
                inm2 = i2
                ic1 = jc + 1
c               Got atom name, now go get file name:
                goto 10
              else
c               Not an alias, end of file name instead:
                goto 100
              endif
            endif
   40     continue
c         Have nothing but blanks left, must be file name:
          goto 100
        elseif( atmfil(ic:ic) .eq. equals )then
c         Make sure not *two* equals around:
          if( got_name ) call STOPXERR( 'at=file=/two equals in atom' )
          i2 = ic - 1
c         Found an alias
          got_name = .true.
          inm1 = i1
          inm2 = i2
          ic1 = ic + 1
          goto 10
        endif
c
   50 continue
c
  100 continue
c  We have parsed an atom file name, load it:
      lenfl = i2 - i1 + 1
      filenm(1:lenfl) = atmfil(i1:i2)
c
      if( .not. got_name .and. lenfl.gt.4 )then
c       Try to get atom name out of atom file name:
        it1 = i2 - 3
        it2 = i2
        inm2 = it1 - 1
        if( atmfil(it1:it2) .eq. atmtag .and.
     $      atmfil(inm2:inm2) .ne. slash )then
          got_name = .true.
c         Go looking backwards for directory slash:
          inm1 = i1
          do 140 ic=inm2-1,i1,-1
            if( atmfil(ic:ic) .eq. slash )then
              inm1 = ic + 1
              goto 150
            endif
  140     continue
        endif
      endif
c
  150 continue
      if( got_name )then
c       Guarantee that atom name fits in space alloted:
        if( inm2 .gt. (inm1+natnm-1 ) ) inm2 = inm1+natnm-1
        lennm = inm2 - inm1 + 1
        atmnm(1:lennm) = atmfil(inm1:inm2)
      endif
c      print*, "Debug printing in atomnm (delete):"
c      print*, atmtag,blank,equals,slash
c      print*, nflnm,natnm
c      print*, lennm,inm1,inm2
c      print*, atmnm
c      print*, filenm
c
c    That's all Folks!
c
      RETURN
      END
