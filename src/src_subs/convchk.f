c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CONVCHK
c
c
      subroutine CONVCHK( IWR, chglbl, istat, iter,itstop, convrg,
     $ chgmax,chgrms,  nd1,nd2, dvec )
c-->                           (nd1*nd2)
c---------------------------------------------------------------
c Purpose: test for convergence/completion
c
c Written: Peter A. Schultz, 25-September-1997, v2.21
c
c Revision history:
c   5Dec22-PAS/2.69: tweak itstop convergence
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
      CHARACTER*(*)  chglbl
      DIMENSION  dvec(nd1,*)
C     DIMENSION  dvec(nd1,nd2)
c
c >>>> EXECUTABLE CODE:
c
      nmx = 1
      chgmax = zero
      chgrms = zero
      do  id2=1,nd2
        chg = zero
        do  id1=1,nd1
          chg = chg + dvec(id1,id2)**2
        enddo
        chgrms = chgrms + chg
        if( chg .gt. chgmax )then
          nmx = id2
          chgmax = chg
        endif
      enddo
      chgmax =  SQRT( chgmax )
      chgrms = chgrms / DBLE( nd2 )
      chgrms = SQRT( chgrms )
c
      write(IWR,9005)  chglbl, chgrms,chgmax, nmx
 9005 format(1x,'rms,max ',a,'=',1p2d12.4,' for component#',i9)
c
      istat = 0
      if( iter .ge. itstop ) istat = 1
      if( chgmax .lt. convrg ) istat = 2
c
      RETURN
      END
c
c
c  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SCFCHECK
c
c
      subroutine SCFCHECK( IWR, chglbl, iscfstat,
     $ iter,itstop, convrg,  chgmax,chgrms,  ndvloc,
     $ dvscf )
c --> (ndvloc)
c---------------------------------------------------------------
c Purpose: test for convergence/completion, k-parallel distributed vector
c
c Written: Andrew C. Pineda and PAS, Nov-2009 - Apr 2012, for 2.63
c          Adapted from convchk
c
c Revision history:
c---------------------------------------------------------------
c Notes:
c  On entry and exit:  dvscf=dv
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
      CHARACTER*(*)  chglbl
      DIMENSION  dvscf(*)  
C     DIMENSION  dvscf(ndvloc)
c
c >>>> EXECUTABLE CODE:
c
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
c  Set local variables: 
      idvmx = 1
      chgmax = zero
      chgrms = zero
c
c  The next block is only done by processors having part of this vector
      if( ndvloc .gt. 0 )then
c       Accumulate chgmax & idvmx, and chgrms on active processors
        do 100 idv=1,ndvloc
          chg = dvscf(idv)**2
          chgrms = chgrms + chg
          if( chg .gt. chgmax )then
            idvmx = idv
            chgmax = chg
          endif
  100   continue
        chgmax = SQRT( chgmax )
      endif
c
c  Reduce ndv, chgmax & idvmx, and chgrms back to master.
c
      ndvtot = 0
      chgrmstot = zero
      lenmsg = 1
c
      do 500 jproc=0,nprocs-1
c
        if( jproc .eq. master )then
c         Local accumulation on the master
c
          if( iproc .eq. master .and. ndvloc .gt. 0 )then
            ndvtot = ndvtot + ndvloc
            chgrmstot = chgrmstot + chgrms
          endif
c
        else
c         Accumulate results from procs at master
c
          if( iproc .eq. master )then
c
c           Image master receives size of vector on jproc ...
            call MPRECVI( jproc, lenmsg, ndvjpr, icomm )
c            ... and receives results if jproc had part of vector
            if( ndvjpr .gt. 0 ) then
              call MPRECVI(  jproc, lenmsg, idvmxjpr,  icomm )
              call MPRECVR8( jproc, lenmsg, chgmaxjpr, icomm )
              call MPRECVR8( jproc, lenmsg, chgrmsjpr, icomm )
              if( chgmaxjpr .gt. chgmax )then
                idvmx = ndvtot + idvmxjpr
                chgmax = chgmaxjpr
              endif
              chgrmstot = chgrmstot + chgrmsjpr
              ndvtot = ndvtot + ndvjpr
            endif
c
          elseif( iproc .eq. jproc )then
c
c           This proc sends size of its vector to image master
            call MPSENDI( master, lenmsg, ndvloc, icomm )
c           Send convergnce analysis to master, if jproc had part of vector
            if( ndvloc .gt. 0 )then
              call MPSENDI(  master, lenmsg, idvmx,  icomm )
              call MPSENDR8( master, lenmsg, chgmax, icomm )
              call MPSENDR8( master, lenmsg, chgrms, icomm )
            endif
c
          endif
        endif
  500 continue
c
      if( iproc .eq. master )then
c       Only master dumps a record to the listing file
        chgrms = chgrmstot / DBLE( ndvtot )
        chgrms = SQRT( chgrms )
c
        write(IWR,9005)  chglbl, chgrms,chgmax, idvmx
 9005   format(1x,'rms,max ',a,'=',1p2d12.4,' for component#',i10)
c
        iscfstat = 0
        if( iter .ge. itstop ) iscfstat = 1
        if( chgmax .lt. convrg ) iscfstat = 2
      endif
c
c Broadcast status to all processors:
      if( nprocs .gt. 1 ) call MPBCASTI( master, 1, iscfstat, icomm )
c
c    That's all, Folks!
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
c End:
