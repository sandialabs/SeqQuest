c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TIMSTUFF
c
c
c===============================================================
c Purpose: this module handles timing instrumentation for the code
c
c Written: Peter A. Schultz, 17-July-2003, for version 2.56
c          (expanded from original timing routine)
c
c Revision history:
c  28Oct04-PAS/2.59: to install wall time routine in timer
c  30Apr08-PAS/2.62: to deal with daily reset in wall timer
c===============================================================
c
c Routines contained:
c  TIMINIT( itimer) - initializes timing module
c  TIMFINI( IO_UNIT ) - final timing summary
c  TIMER( msg ) - timing statement
c
c External dependencies:
c  GETTLFT( t ) - returns time-left from system
c  GETTUSE( ttotl, tuser ) - gets time usage from system
c  GETWALLT( twall ) - gets wall time from system
c  FLGETIWR( iwr ) - gets output unit
c
c===============================================================
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TIMINIT
c
c
      subroutine TIMINIT( jtimer )
c---------------------------------------------------------------
c Purpose: initialize timing variables
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /TIMSTUFF/  timetotl,timeuser,timesys,timelim,
     $ timewall,tim0wall,timewoff, itimer
c
c >>>> EXECUTABLE CODE:
c
      timetotl = 0.
      timeuser = 0.
      timesys  = 0.
      timelim  = -1.0
c
      tim0wall = 0.
      timewall = 0.
      timewoff = 0.
c
      if( jtimer .eq. 0 )then
c       Turn timer off entirely
        itimer = 0
      elseif( jtimer .eq. 1 )then
c       System time
        itimer = 1
      elseif( jtimer .eq. 2 )then
c       Wall timer
        itimer = 2
      else
c       Default is system timer
        itimer = 1
      endif
c
c     Initialize the wall timer
      call GETWALLT( twall )
      tim0wall = twall
      timewall = twall
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TIMFINI
c
c
      subroutine TIMFINI( IWRTIM )
c---------------------------------------------------------------
c Purpose: handle closing down of final timing instrumentation
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /TIMSTUFF/  timetotl,timeuser,timesys,timelim,
     $ timewall,tim0wall,timewoff, itimer
c
      if( IWRTIM.gt.0 )then
c
        if( itimer.eq.1 )then
c         Without wall timer:
          write(IWRTIM,901) timetotl,timeuser,timesys
  901     format(1x,'Total TIME used=',f10.1,' seconds (user=',f10.1,
     $     '  sys=',f10.1,')')
c
        elseif( itimer.eq.2 )then
c         With wall timer:
          totlwall = timewall - tim0wall
          timelost = totlwall - timetotl
          write(IWRTIM,902) totlwall,timeuser,timesys,timelost
  902     format(1x,'Total TIME=',f10.1,' seconds (usr=',f10.1,
     $     '  sys=',f10.1,'  lost=',f10.1,')')
c
        endif
c
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TIMER
c
c
      subroutine TIMER( sublabl )
c---------------------------------------------------------------
c Purpose: track running time for code, with labels
c
c Revision history:
c  30Apr08-PAS/2.62: account for day rollover in wall clock
c  17Jul03-PAS/2.56: get IWR from FL module, expand timing info
c  28Aug02-PAS/2.54: increase format for time used
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /TIMSTUFF/  timetotl,timeuser,timesys,timelim,
     $ timewall,tim0wall,timewoff, itimer
c
      CHARACTER*(*)  sublabl
      CHARACTER*(4)  timlbl
c
c >>>> EXECUTABLE CODE:
c
c Retrieve total running time, and user time from system:
      call GETTUSE( ttotl, tuser )
c
c Compute elapsed (tot,usr,sys) times since last record
      oldttotl = timetotl
      timetotl = ttotl
      delttotl = timetotl - oldttotl
c
      oldtuser = timeuser
      timeuser = tuser
      deltuser = timeuser - oldtuser
c
      timesys = timetotl - timeuser
      deltsys = delttotl - deltuser
c
      ioffday = 0
      if( itimer .eq. 2 )then
c       Elapsed wall time since last entry:
        call GETWALLT( twall )
        twall = twall + timewoff
c       Check that we have not had a roll-over of the clock
        if( twall .lt. timewall )then
c         Our wall clock is a day clock, and reset by 86400sec.
c         Adjust our time offset by a day.
cpas: temporary stamp of timer rollover
          write(6,*) 'TIMER rollover, twall,timewoff=',twall,timewoff
          ioffday = ioffday + 1
          delclock = 86400.d0
          timewoff = timewoff + delclock
          twall = twall + delclock
c         In principle, more than one day could have passed between
c         timer calls, and offset error could be more than one day.
c         However, with runs that long, the timer is irrelevant and
c         the time that comes out of that timer unimportant.
        endif
        oldwall = timewall
        timewall = twall
        delwall = timewall - oldwall
        totlwall = timewall - tim0wall
c
        deltime = delwall
        timetot = totlwall
      endif
c
c Retrieve output unit, and write out data:
      call FLGETIWR( IWR )
      if( IWR .gt. 0 )then
c
C       write(IWR,910)  sublabl, delttotl,timetotl
C  910  format(/1x,a,'+',f9.3,' =',f10.1,' seconds used')
c
        if( itimer .eq. 1 )then
          write(IWR,911)  sublabl, delttotl,timetotl
     $                           , deltuser,deltsys
  911     format(/1x,a,'+',f8.2,' =',f10.1,' seconds used'
     $               ,' (',f8.2,' usr',f8.2,' sys)')
c
        elseif( itimer .eq. 2 )then
          write(IWR,912)  sublabl, deltime, timetot
     $                           , deltuser,deltsys
  912     format(/1x,a,'+',f8.2,' =',f10.1,' seconds wall'
     $               ,' (',f8.2,' usr',f8.2,' sys)')
          if( ioffday.ne.0 )then
            write(IWR,9121)
 9121       format(1x,' (after offset of 86400 seconds rollover)' )
          endif
        endif
c
        call MTBUFF
      endif
c
c    That's all Folks!
c
      RETURN
      END
