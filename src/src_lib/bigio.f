c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BIGIO
c
c
c===============================================================
c Purpose: these routines manage big i/o to files
c
c Written: Peter A. Schultz, 22-July-2003, for 2.56
c          (expanded from original big i/o routines)
c
c Revision history:
c  26Mar20-PAS/2.68e add READB1 to read partial vector (ivecfl)
c   9Apr15-PAS/2.66: install MPBC routine (master-only i/o + bcasts)
c===============================================================
c
c Notes:
c  These routines prompted by problem encountered on some dec
c  alphas where the OS limited the length of the a single record
c  in a read (but, apparently, not in a write - go figure).
c
c External dependencies:
c  none
c
c Public routines:
c  (administrative)
c    BIGIOINIT( maxreclen ) - Initializes module
c    IOSETRECL( maxreclen ) - Sets max record length
c    IOGETRECL( maxreclen ) - Gets max record length
c  (functional)
c    WRITBIG( ifile, lenvec, vec ) - write vector to unit
c    READBIG( ifile, lenvec, vec ) - read vector from unit
c    READB1 ( ifile, lendim, lenvec, vec ) - read part of vector from unit
c    BACKBIG( ifile, lenvec )      - backspace over vector on unit
c    READSKP( ifile, lenvec )      - skip over vector on unit
c  (master-parallel - i.e. master reads and broadcasts)
c    MPBCREWIND( ifile )               - rewind unit
c    MPBCREADBIG( ifile, lenvec, vec ) - read vector from unit
c    MPBCREADSKP( ifile, lenvec, vec ) - skip over vector on unit
c
c Variables:
c  ifile  - file unit number being written to (read from)
c  vec()  - the vector to be written (read)
c  lenvec - length of the vector
c  maxreclen = >0 maximum length of a single write/read record (r8)
c            =  0 uses default for max length of record
c            = <0 unlimited record length
c
c===============================================================
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BIGIOINIT
c
c
      subroutine BIGIOINIT( maxrecl )
c---------------------------------------------------------------
c Purpose: initialize big i/o stuff
c
c Written: Peter A. Schultz, 22-July-2003
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      maxlen = maxrecl
      call IOSETRECL( maxlen )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> IOSETRECL
c
c
      subroutine IOSETRECL( maxrecl )
c---------------------------------------------------------------
c Purpose: Set/get max record length
c
c Written: Peter A. Schultz, 22-July-2003
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /BIGIO/ maxreclen
c
      maxreclen = maxrecl
c
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> IOGETRECL
c
      entry IOGETRECL( maxrecl )
c
      maxrecl = maxreclen
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WRITBIG
c
c
      subroutine WRITBIG( ifile, lenvec, vec )
c---------------------------------------------------------------
c Purpose: write/read big line within max size records
c          ***** MACHINE DEPENDENT *****
c
c Written: Peter A. Schultz, 15-February-1999, for v2.31
c
c Revision history:
c  22Jul03-PAS/2.56: modularized, and removed from 'utl'
c  21Jun01-PAS/2.47a: moved to utl, enddo
c---------------------------------------------------------------
c
c Variables:
c  ifile  - file unit number being written to (read from)
c  vec()  - the vector to be written (read)
c  lenvec - length of the vector
c  maxlen - maximum length of a single write/read record
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  (maxrecdef= 1 000 000)
c
      COMMON  /BIGIO/ maxreclen
c
      DIMENSION  vec(lenvec)
c
c >>>> EXECUTABLE CODE:
c
      maxlen = maxreclen
      if( maxlen.eq.0 )then
c       If undefined, use internal default
        maxlen = maxrecdef
      elseif( maxlen.lt.0 )then
c       If turned OFF (maxlen<0), write full vectors
        maxlen = lenvec
      endif
c
      nvrec = lenvec/maxlen
      if( nvrec.eq.0 )then
c       The vector will fit in a single record:
        write(ifile)  vec
      else
c       The vector needs to be broken up:
        iv0 = 0
        do  ivrec=1,nvrec
          write(ifile)  (vec(iv),iv=iv0+1,iv0+maxlen)
          iv0 = iv0 + maxlen
        enddo
        if( iv0.lt.lenvec )  write(ifile) (vec(iv),iv=iv0+1,lenvec)
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> READBIG
c
c
      entry READBIG( ifile, lenvec, vec )
c
      maxlen = maxreclen
      if( maxlen.eq.0 )then
        maxlen = maxrecdef
      elseif( maxlen.lt.0 )then
        maxlen = lenvec
      endif
c
      nvrec = lenvec/maxlen
      if( nvrec.eq.0 )then
c       The vector will fit in a single record:
        read(ifile)  vec
      else
c       The vector needs to be broken up:
        iv0 = 0
        do  ivrec=1,nvrec
          read(ifile)  (vec(iv),iv=iv0+1,iv0+maxlen)
          iv0 = iv0 + maxlen
        enddo
        if( iv0.lt.lenvec ) read(ifile)  (vec(iv),iv=iv0+1,lenvec)
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> READB1
c
c
      entry READB1 ( ifile, lendim,lenvec, vec )
c
      maxlen = maxreclen
      if( maxlen.eq.0 )then
        maxlen = maxrecdef
      elseif( maxlen.lt.0 )then
        maxlen = lenvec
      endif
c
      nvrec = lenvec / maxlen
      ndrec = lendim / maxlen
      if( ndrec*maxlen .lt. lendim ) ndrec = ndrec + 1
c
      if( nvrec.eq.0 )then
c       The vector will fit in a single record:
        ndrec = ndrec - 1
        read(ifile)  vec
      else
c       The vector needs to be broken up:
        iv0 = 0
        do  ivrec=1,nvrec
          ndrec = ndrec - 1
          read(ifile)  (vec(iv),iv=iv0+1,iv0+maxlen)
          iv0 = iv0 + maxlen
        enddo
        if( iv0.lt.lenvec )then
          ndrec = ndrec - 1
          read(ifile)  (vec(iv),iv=iv0+1,lenvec)
        endif
      endif
c
      if( ndrec.gt.0 )then
c       Skip through (unread) rest of the whole vector
        do  idrec=1,ndrec
          read(ifile)
        enddo
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BACKBIG
c
c
      entry BACKBIG( ifile, lenvec )
c
      maxlen = maxreclen
      if( maxlen.eq.0 )then
        maxlen = maxrecdef
      elseif( maxlen.lt.0 )then
        maxlen = lenvec
      endif
c
      nvrec = lenvec/maxlen
      if( nvrec.eq.0 )then
c       The vector will fit in a single record:
        BACKSPACE( ifile )
      else
c       The vector needs to be broken up:
        iv0 = 0
        do  ivrec=1,nvrec
          BACKSPACE( ifile )
          iv0 = iv0 + maxlen
        enddo
        if( iv0.lt.lenvec ) BACKSPACE( ifile )
      endif
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> READSKP
c
c
      entry READSKP( ifile, lenvec )
c
      maxlen = maxreclen
      if( maxlen.eq.0 )then
        maxlen = maxrecdef
      elseif( maxlen.lt.0 )then
        maxlen = lenvec
      endif
c
      nvrec = lenvec/maxlen
      if( nvrec.eq.0 )then
c       The vector will fit in a single record:
        read(ifile)
      else
c       The vector needs to be broken up:
        iv0 = 0
        do  ivrec=1,nvrec
          read(ifile)
          iv0 = iv0 + maxlen
        enddo
        if( iv0.lt.lenvec ) read(ifile)
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBCREWIND
c
c
      subroutine MPBCREWIND( ifile )
c---------------------------------------------------------------
c Purpose: master rewind file
c
c Written: Peter A. Schultz, 28-January-2008, for 2.61
c
c Revision history:
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c >>>> EXECUTABLE CODE:
c
      call MPNODES( nprocs )
      call MPNODE0( master )
      call MPNODE( iproc )
      call MPCOMM( icomm )
c
      if( iproc .eq. master )then
        REWIND( ifile )
      endif
c
      RETURN
      END
c
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBCREADBIG
c
c
      subroutine MPBCREADBIG( ifile, nvec, vec )
c---------------------------------------------------------------
c Purpose: master read vector and broadcast to all
c
c Written: Peter A. Schultz, 27-July-2007, for 2.60
c
c Revision history:
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  vec(*)
c
c >>>> EXECUTABLE CODE:
c
      call MPNODES( nprocs )
      call MPNODE0( master )
      call MPNODE( iproc )
      call MPCOMM( icomm )
c
      if( iproc .eq. master )then
        call READBIG( ifile, nvec, vec )
      endif
c
      if( nprocs .gt. 1 )then
        lenmsg = nvec
        call MPBCAST8( master, lenmsg, vec, icomm )
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBCREADSKP
c
c
      subroutine MPBCREADSKP( ifile, nvec, vec )
c---------------------------------------------------------------
c Purpose: master read vector and broadcast to all
c
c Written: Peter A. Schultz, 27-July-2007, for 2.60
c
c Revision history:
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  vec(*)
c
c >>>> EXECUTABLE CODE:
c
      call MPNODES( nprocs )
      call MPNODE0( master )
      call MPNODE( iproc )
      call MPCOMM( icomm )
c
      if( iproc .eq. master )then
        call READSKP( ifile, nvec )
      endif
c
c    That's all Folks!
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
