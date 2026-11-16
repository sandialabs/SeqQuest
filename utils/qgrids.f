      PROGRAM QGRIDS
c###############################################################
c Program: QGRIDS
c Purpose: Retrieve and do some limited processing of grid data
c          produced by a DFT Quest calculation.
c Written: Peter A. Schultz, 13-August-2003
c          Design based on SeqQuest 2.56(alpha)
c Revision history:
c 17Aug10/PAS/1.0e: bugfix in PLOTGRID call
c 23Apr10/PAS/1.0d: a couple more cosmetic changes
c  4Dec04-PAS/1.0c: fix decl of filenm
c 27Apr04-PAS/1.0b: fix bug in bugfix
c  7Oct03-PAS/1.0a: fix bug constraining output plot grid
c###############################################################
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER (maxwkd= 8 000 000)
c
      COMMON  /BIGIO/ maxreclen
c
      CHARACTER  lastup*11
      DATA       lastup /'17Aug10-PAS'/
      DIMENSION  wk(maxwkd)
c
      DIMENSION  rprim(3,3)
      DIMENSION  gridvec(3,3)
      DIMENSION  orig(3)
      DIMENSION  vec3(3)
c
      LOGICAL    qdebug
      CHARACTER  opt*(1), yesno*(1)
      CHARACTER  jobnm*(32), filenm*(80)
      DATA       IRD,IWR / 5,6 /
      DATA       IGRDFL,IOUTFL / 7,8 /
      DATA       zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      qdebug = .false.
      write(IWR,9000)  lastup
 9000 format(
     $  1x,'Program QGRIDS - last update: ',a11,
     $/ 1x,'Purpose: retrieve/manipulated Quest grid output'/
     $/ 1x,'Needs access to Quest files with suffixes .grd and .grd0'/
     $/ 1x,'All Quest jobs in this analysis must have same unit cell',
     $/ 1x,'This code works by constructing a work grid field W',
     $/ 1x,'It updates/adds-to this workspace W as instructed, from',
     $/ 1x,'either Quest output grid fields (which must have been',
     $/ 1x,'previously Loaded), or saved grid fields that had been',
     $/ 1x,'Read into a temp field earlier.  The update of W does NOT'
     $/ 1x,'occur when the grids are Loaded or Read.  The workspace W'
     $/ 1x,'is only updated when it is instructed to do so after the'
     $/ 1x,'fields have been Read in or Loaded from disk.',
     $/ 1x,'Once you have constructed the work grid field you desire',
     $/ 1x,'through load/read and update, you o(utput) the result'
     $/ )
c
c     Initialize maximum record-length for I/O:
c     (this must match the size used in the Quest calculation)
c     
      maxreclen = 1 000 000
   10 continue
      call BIGIOINIT( maxreclen )
      write(IWR,*) 'maxreclen default=',maxreclen
c
   20 continue
      write(IWR,*) 'Enter Quest primitive lattice vectors A,B,C (bohr):'
      read(IRD,*)   rprim
      call CROSS( rprim(1,1),rprim(1,2), vec3 )
      volume = ABS( rprim(1,3)*vec3(1) + rprim(2,3)*vec3(2)
     $            + rprim(3,3)*vec3(3) )
      write(IWR,'(a,f16.6)') 'Cell volume (bohr^3)=',volume
c
   30 continue
      write(IWR,*) 'Enter grid dimensions of Quest calculation:'
      read(IRD,*)   n1r,n2r,n3r
      if( n1r.lt.1 .or. n2r.lt.1 .or. n3r.lt.1 )then
c       Impossible grid dimensions
        write(IWR,*) 'Grid dimensions out-of-bounds'
        goto 30
      endif
      npts = n1r*n2r*n3r
c Grid integration weight is simply the volume of single element:
      weight = volume / DBLE( npts )
c
c Generate the grid primitive vectors from the lattice primitive
c  vectors and the grid dimensions:
   40 continue
      do  i=1,3
        gridvec(i,1) = rprim(i,1) / n1r
        gridvec(i,2) = rprim(i,2) / n2r
        gridvec(i,3) = rprim(i,3) / n3r
      enddo
c
c Define workspaces for grid fields:
c   i_fld  = the field we are constructing
c   i_rho0 = the reference density on the grid
c   i_drho = the delta-density (scf-ref) on the grid
c   i_dves = the delta-electrostatic potential
c   i_tmp  = scratch space
      i_start = 1
      i_fld  = i_start
      i_rho0 = i_start + 1*npts
      i_drho = i_start + 2*npts
      i_dves = i_start + 3*npts
      i_tmp  = i_start + 4*npts
      i_end  = i_start + 5*npts - 1
      if( i_end .gt. maxwkd )then
        write(IWR,*) '***** ERROR: not enough memory dimensioned'
        write(IWR,*) '>>>>> maxwkd; needed=',maxwkd,i_end
      endif
c
c Clear out grid field:
   80 continue
      call MKZERO( npts, wk(i_fld) )
      write(IWR,9080) n1r,n2r,n3r, npts
 9080 format(1x,'Work grid field empty and ready, dims=',
     $  i4, ' x', i4, ' x', i4, ' =' ,i10 )
c
c Record status of grid field as empty, and no Quest job loaded yet:
      ifldstat = 0
      ijobstat = 0
c
c **** Help message ****
c
   90 continue
      write(IWR,9090)
 9090 format(/1x,'Options:'
     $/5x,'l(oad)       - load Quest output grid fields into program',
     $/5x,'m(axreclen)  - reset max record length for quest grid i/o',
     $/5x,' ',
     $/5x,'u(pdate)     - update work grid field W from Quest fields',
     $/5x,' ',
     $/5x,'o(utput)     - output (ASCII) grid field W for plotting',
     $/5x,' ',
     $/5x,'i(nitialize) - clear work grid field W',
     $/5x,'s(cale)      - scale work grid field W by constant',
     $/5x,'w(rite)      - write work grid field W to disk for later',
     $/5x,'r(ead)       - read saved work grid field W from disk',
     $/5x,' ',
     $/5x,'q(uit)       - quit',
     $/5x,' ',
     $/5x,'h(elp)       - write this listing'/)
c
c **** Main options menu ****
c
  100 continue
      write(IWR,*) 'Enter option:'
      read(IRD,'(a)')  opt
c
      if( opt.eq.'i' .or. opt.eq.'I' )then
c       Initialize work grid field to zero
c
        if( ifldstat.eq.1 )then
          write(IWR,*) 'Your constructed grid field has not been saved'
          write(IWR,*) 'Do you really want to destroy it (def=n)?'
          read(IRD,'(a)')  yesno
          if( yesno.eq.'y' .or. yesno.eq.'Y' ) goto 80
          goto 100
        endif
        goto 80
c
      elseif( opt.eq.'o' .or. opt.eq.'O' )then
c       Call routine to decide how to output the data:
c
        if( ifldstat.eq.0 ) write(IWR,*) 'Note: your field is blank'
        call PLOTGRID( IRD,IWR, IGRDFL,
     $   n1r,n2r,n3r,npts, gridvec,
     $   wk(i_fld), wk(i_tmp) )
c       Change field status to "saved" (2)
        ifldstat = 2
        goto 100
c
      elseif( opt.eq.'w' .or. opt.eq.'W' )then
c       Write full work grid field to disk
c
        write(IWR,*) 'Enter (new) file name to write grid field to:'
        read(IRD,'(a)')  filenm
        if( filenm.eq. ' ' ) goto 100
        OPEN( unit=IGRDFL, file=filenm, status='NEW',
     $   form='unformatted',err=100)
c       Write it using standard (record-size limited) binary:
        call WRITBIG( IGRDFL, npts, wk(i_fld) )
        write(IWR,*) 'Full grid written to binary file'
c       Change field status to "saved" (2)
        ifldstat = 2
        goto 100
c
      elseif( opt.eq.'r' .or. opt.eq.'R' )then
c       Read saved grid field from disk, and do something with it.
c
        write(IWR,*) 'Enter (old) file name to read grid field from:'
        read(IRD,'(a)')  filenm
        if( filenm.eq. ' ' ) goto 100
        OPEN( unit=IGRDFL, file=filenm, status='OLD',
     $   form='unformatted',err=100)
c       Read it using standard (record-size limited) binary:
        call READBIG( IGRDFL, npts, wk(i_tmp) )
c
        write(IWR,*) 'Add into working grid: W(r) = W(r) + c*G(r)'
        write(IWR,*) 'Enter constant c (0=skip):'
        read(IRD,*)  scale
        if( scale.eq.zero ) goto 100
        call FPLUSCG( scale, npts, wk(i_fld), wk(i_tmp) )
        ifldstat = 1
        goto 100
c
      elseif( opt.eq.'s' .or. opt.eq.'S' )then
c       Scale the work grid field by a constant:
c
        write(IWR,*) 'Enter constant to scale work field (0=skip):'
        read(IRD,*)  scale
        if( scale.eq.zero ) goto 100
        call FSCALE( scale, npts, wk(i_fld) )
        ifldstat = 1
        goto 100
c
      elseif( opt.eq.'m' .or. opt.eq.'M' )then
c       Reset maxreclen to be compatible with quest output
        write(IWR,*) 'Enter max record length for binary i/o:'
        read(IRD,*)  maxrecl
        call IOSETRECL( maxrecl )
        goto 100
c
      elseif( opt.eq.'l' .or. opt.eq.'L' )then
c       Load Quest grid outputs:
c
        if( ijobstat.eq.1 )then
          write(IWR,*) 'Not yet used (UPDATE) previous Quest data'
          write(IWR,*) 'Do you really want to overload it (def=n)?'
          read(IRD,'(a)') yesno
          if( yesno.ne.'y' .and. yesno.ne.'Y' ) goto 100
        endif
        write(IWR,*) 'To load Quest grid output -'
        write(IWR,*) 'Enter base name XXXX (e.g., XXXX.grd):'
        read(IRD,'(a)')  jobnm
        call STRPARS( jobnm, 1,32, j1,jn )
        if( j1.eq.0 ) goto 100
c
c       Retrieve data from Quest files:
        call GETFIELDS( IRD,IWR, IGRDFL, jobnm(j1:jn), npts, ierr,
     $   wk(i_rho0),wk(i_drho),wk(i_dves),wk(i_tmp) )
c
        if( ierr.eq.0 )then
          write(IWR,*) 'Successfully loaded rho0/drho/dves from files'
          call FSUM( rho0tot, npts, wk(i_rho0) )
          rho0tot = weight*rho0tot
          call FSUM( drhotot, npts, wk(i_drho) )
          drhotot = weight*drhotot
          rhotot = rho0tot + drhotot
          call FABSSUM( absdrho, npts, wk(i_drho) )
          absdrho = weight*absdrho
          if( qdebug )then
            call FSUM( vestot , npts, wk(i_dves) )
            call FABSSUM( vestot0 , npts, wk(i_dves) )
            write(IWR,9213) 'EVES SUM DEBUG=',vestot,vestot0
          endif
          call FDOT( eves, npts, wk(i_drho), wk(i_dves) )
          eves = weight*eves / 2.d0
 9213     format(1x,a,2f20.12)
          write(IWR,9213) ' Grid-integrated reference density=',rho0tot
          write(IWR,9213) ' Grid-integrated total density=',rhotot
          write(IWR,9213) ' Grid-integrated delta movement=',absdrho
          write(IWR,9213) ' non-spherical es ENERGY=',eves
          write(IWR,*) 'Use **update** to add into work field W'
          ijobstat = 1
        else
          write(IWR,*) 'Failed - did not load any data'
          ijobstat = 0
        endif
        goto 100
c
      elseif( opt.eq.'u' .or. opt.eq.'U' )then
c       Update work grid field with data from Quest job:
c
        if( ijobstat.eq.0 )then
          write(IWR,*) 'No Quest grid data currently loaded'
          write(IWR,*) 'Must *load* Quest data before updating'
          goto 100
        endif
        write(IWR,*) 'Update work grid field W from Quest data:'
        write(IWR,*) ' W(r) = W(r) + c*Q(r)'
        write(IWR,*) 'You will be prompted for factor c, but first'
        write(IWR,*) 'Choose field Q to add (rho0=1,drho=2,dves=3):'
        read(IRD,*) iqfld
        if(     iqfld.eq.1 )then
          i_update = i_rho0
          write(IWR,*) '    Chosen: Reference density rho0'
        elseif( iqfld.eq.2 )then
          i_update = i_drho
          write(IWR,*) '    Chosen: Delta-density drho'
        elseif( iqfld.eq.3 )then
          i_update = i_dves
          write(IWR,*) '    Chosen: Delta electrostatic potential'
        else
          write(IWR,*) '***** ERROR: field out-of-bounds'
          write(IWR,*) '>>>>> No action taken'
          goto 100
        endif
        write(IWR,*) 'Enter addition coefficient c:'
        read(IRD,*)  scale
        if( scale.eq.zero )then
          write(IWR,*) '***** WARNING: W is unchanged, coeff was zero'
          goto 100
        endif
        call FPLUSCG( scale, npts, wk(i_fld), wk(i_update) )
c       Note Quest data as used, and work grid as modified
        ijobstat = 2
        ifldstat = 1
        goto 100
c
      elseif( opt.eq.'q' .or. opt.eq.'Q' )then
c       Quit
c
        if( ifldstat.eq.1 )then
          write(IWR,*) 'Your constructed grid field not saved/output'
          write(IWR,*) 'Do you really want to quit (def=n)?'
          read(IRD,'(a)')  yesno
          if( yesno.eq.'y' .or. yesno.eq.'Y' ) goto 999
          goto 100
        elseif( ijobstat.eq.1 )then
          write(IWR,*) 'You have not used the Quest data you loaded'
          write(IWR,*) 'Do you really want to quit (def=n)?'
          read(IRD,'(a)')  yesno
          if( yesno.eq.'y' .or. yesno.eq.'Y' ) goto 999
          goto 100
        endif
        goto 999
c
      endif
c
      write(IWR,*) 'Option unknown: ',opt
      goto 90
c
  999 continue
      write(IWR,*) 'All done'
c
      STOP
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PLOTGRID
c
c
      subroutine PLOTGRID( IRD,IWR, IGRDFL,
     $    n1r,n2r,n3r,npts, gridvec,
     $    wfld, wplot )
c---------------------------------------------------------------
c Purpose: manage output of grid field
c Written: Peter A. Schultz, 13-August-2003 from 2.56(alpha)
c Revision history:
c  none
c---------------------------------------------------------------
c
c The grid field has been constructed in "wfld()"
c The grid intervals in the three dimensions are specified
c   in "gridvec(3,3)".
c
c This routine has lots of options to make it very flexible
c   to change what is output and the format it is output in.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input grid primitive vectors:
      DIMENSION  gridvec(3,3)
c Input work grid field to be processed:
      DIMENSION  wfld(n1r,n2r,n3r)
c Scratch space
      DIMENSION  wplot(npts)
c Local declarations:
      DIMENSION  igrdim(3)
      DIMENSION  ipldim(3), iplorig(3), iplend(3,3), iplint(3,3)
      DIMENSION  rpt(3), wpts(10)
      LOGICAL    r_write
      CHARACTER  filenm*(80)
      CHARACTER  opt*(1)
      CHARACTER  rfmt*(80),vfmt*(80)
      DATA       zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
      write(IWR,*) 'Entered output routine PLOTGRID'
c Load total working grid dimensions into vector
      igrdim(1) = n1r
      igrdim(2) = n2r
      igrdim(3) = n3r
c
   20 continue
      write(IWR,*)
      write(IWR,*) 'Choose dimensionality of grid desired for output'
      write(IWR,*) 'Write 1D line (1), 2D grid (2), or 3D grid (3)?:'
      read(IRD,*,err=20)   ndplot
      if( ndplot.lt.1 .or. ndplot.gt.3 ) goto 20
c
c Select output grid:
   30 continue
      write(IWR,9020)
 9020 format(1x,'Define plotting grid of points for output',
     $/1x,'You will be asked for (as integer 3-vectors):',
     $/1x,'  (1) origin of plot grid',
     $/1x,'  (2) end point of each axis',
     $/1x,'  (2) number of points from origin to end of axis')
      write(IWR,9021)
 9021 format(1x,'2D Example: for 3D Quest grid that is 8x8x20',
     $/1x,'To define an 8x20 2D grid along the 110:'
     $/1x,'  Origin = (1 1 1)',
     $/1x,'                  endpoint    points',
     $/1x,'  First axis=     (8 8 1)       8 ',
     $/1x,'  Second axis=    (1 1 20)      20',
     $/1x,'This will result in a 8x20 grid of points being written')
c
   40 continue
      write(IWR,'(/a,4i4)') 'Dimension of total grid:',n1r,n2r,n3r
      write(IWR,*) 'Enter origin (3 integers):'
      read(IRD,*) iplorig
      call CHECKDIM( igrdim, iplorig, ierr )
      if( ierr.ne.0 )then
        write(IWR,*) 'ERROR: Origin is outside of grid'
        goto 40
      endif
c
   45 continue
c     Initialize output plot dimensions to zero.
      do  idim=1,3
        ipldim(idim) = 1
        do  i=1,3
          iplint(i,idim) = 0
        enddo
      enddo
c
c     Get endpoints of exes, and dimension of output grid
      do 100 idplot=1,ndplot
c
        write(IWR,'(a,i2,a)')'Enter endpoint (3-int) of axis',idplot,':'
        read(IRD,*)  (iplend(i,idplot),i=1,3)
        call CHECKDIM( igrdim, iplend(1,idplot), ierr )
        if( ierr.ne.0 )then
          write(IWR,*) 'ERROR: Endpoint is outside of grid'
          goto 45
        endif
c
        write(IWR,'(a)') 'Enter number of points along axis:'
        read(IRD,*)  ipld
        if( ipld .lt. 1 )then
          write(IWR,*) 'ERROR: must be at least one'
          goto 45
        endif
        do  idim=1,3
c         Check that origin-endpoint vector is divisible by #pts
          idimlen = IABS( iplorig(idim) - iplend(idim,idplot) ) + 1
          if( idimlen.gt.1 )then
            idimint = idimlen / ipld
            if( idimlen .ne. 1 .and. idimlen.ne.(idimint*ipld) )then
              write(IWR,*) 'ERROR: orig-end not divisible into #pts'
              goto 40
            endif
            if( iplorig(idim).gt.iplend(idim,idplot) ) idimint =-idimint
          else
            idimint = 0
          endif
          iplint(idim,idplot) = idimint
        enddo
        ipldim(idplot) = ipld
c
  100 continue
c
c We have now successfully defined a plotting grid
      write(IWR,*) 'Plotting grid defined'
      write(IWR,'(a,3i4)') 'Dims of plot grid:',(ipldim(i),i=1,ndplot)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c        Write the result to an output file
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Decide whether to write positions with data or not:
  200 continue
      r_write = .false.
      write(IWR,9200)
 9200 format(/1x,'Two choices for output:',
     $/1x,'  Output position with point, 1pt/line: x y z f',
     $/1x,'  Output grid values only',
     $//1x,'Write positions with data (y/n)?')
      read(IRD,'(a)')  opt
      if( opt.eq.'y' .or. opt.eq.'Y' )then
        r_write = .true.
        write(IWR,*) 'Will write field, 1 pt/line, as: {x y z f}'
      elseif( opt.eq.'n' .or. opt.eq.'N' )then
        r_write = .false.
        write(IWR,*) 'Will write grid field values only'
      else
        goto 200
      endif
c
  210 continue
      write(IWR,*) 'Enter name of (new) output file:'
      read(IRD,'(a)')  filenm
      if( filenm.eq. ' ' ) goto 200
      OPEN( unit=IGRDFL, file=filenm, status='NEW', form='FORMATTED',
     $ err=210)
c
c Format for {x y z f}:
      rfmt = '(1x,3f13.6,f15.8)'
c Format for grid field alone:
      vfmt = '(1x,5f15.8)'
      if( r_write )then
c       We are writing the positions AND values, 1pt/line
        i0 = iplorig(1)
        j0 = iplorig(2)
        k0 = iplorig(3)
        ijkpl = 0
        do  kp=0,ipldim(3)-1
          do  jp=0,ipldim(2)-1
            do  ip=0,ipldim(1)-1
c             Locate pt i,j,k on full grid:
              i = i0 + ip*iplint(1,1) + jp*iplint(1,2) + kp*iplint(1,3)
              j = j0 + ip*iplint(2,1) + jp*iplint(2,2) + kp*iplint(2,3)
              k = k0 + ip*iplint(3,1) + jp*iplint(3,2) + kp*iplint(3,3)
c             Construct the position of this point:
              do  ir=1,3
                rpt(ir) = DBLE(i-1)*gridvec(ir,1)
     $                  + DBLE(j-1)*gridvec(ir,2)
     $                  + DBLE(k-1)*gridvec(ir,3)
              enddo
c             Write out position, and value of field at this point
              ijkpl = ijkpl + 1
              write(IGRDFL,rfmt)  rpt, wfld(i,j,k)
            enddo
          enddo
        enddo
        nptpl = ijkpl
      else
c       We are writing just the grid field values:
        i0 = iplorig(1)
        j0 = iplorig(2)
        k0 = iplorig(3)
        ijkpl = 0
        do  kp=0,ipldim(3)-1
          do  jp=0,ipldim(2)-1
            do  ip=0,ipldim(1)-1
c             Locate pt i,j,k on full grid:
              i = i0 + ip*iplint(1,1) + jp*iplint(1,2) + kp*iplint(1,3)
              j = j0 + ip*iplint(2,1) + jp*iplint(2,2) + kp*iplint(2,3)
              k = k0 + ip*iplint(3,1) + jp*iplint(3,2) + kp*iplint(3,3)
              ijkpl = ijkpl + 1
c             Load value at this position
              wplot(ijkpl) = wfld(i,j,k)
            enddo
          enddo
        enddo
        nptpl = ijkpl
c       Write the whole thing:
        write(IGRDFL,vfmt)  (wplot(ijk),ijk=1,nptpl)
      endif
c
      CLOSE( unit=IGRDFL )
c
      write(IWR,*) 'Number of pts written to output file=',nptpl
c
  999 continue
      write(IWR,*) 'Exiting output routine'
c
c >>>> That's all Folks: 
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CHECKDIM
c
c
      subroutine CHECKDIM( idimvec, ivec, ierr )
c---------------------------------------------------------------
c Purpose: check that ivec(3) is inside dims given by idimvec(3)
c Written: Peter A. Schultz, 13-August-2003 from 2.56(alpha)
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  idimvec(3),ivec(3)
c
      ierr = 0
      do  idim=1,3
        if( ivec(idim).lt.1 .or. ivec(idim).gt.idimvec(idim) )then
          ierr = idim
          RETURN
        endif
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FPLUSCG
c
c
       subroutine FPLUSCG( c, npts, f, g )
c---------------------------------------------------------------
c Purpose: F = F + cG, add scaled vector into first vector
c Written: Peter A. Schultz, 13-August-2003 from 2.56(alpha)
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  f(npts), g(npts)
c
      do  i=1,npts
        f(i) = f(i) + c*g(i)
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FSCALE
c
c
       subroutine FSCALE( scale, npts, f )
c---------------------------------------------------------------
c Purpose: F = cF, multiply input array by a constant
c Written: Peter A. Schultz, 13-August-2003 from 2.56(alpha)
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  f(npts)
c
      do  i=1,npts
        f(i) = scale*f(i)
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FSUM
c
c
       subroutine FSUM( sum, npts, f )
c---------------------------------------------------------------
c Purpose: sum of elements of array
c Written: Peter A. Schultz, 14-August-2003 from 2.56(alpha)
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  f(npts)
c
      sum = 0.d0
      do  i=1,npts
        sum = sum + f(i)
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FABSSUM
c
c
       subroutine FABSSUM( sum, npts, f )
c---------------------------------------------------------------
c Purpose: sum of absolute values of elements of array
c Written: Peter A. Schultz, 14-August-2003 from 2.56(alpha)
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  f(npts)
c
      sum = 0.d0
      do  i=1,npts
        sum = sum + ABS( f(i) )
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FDOT
c
c
       subroutine FDOT( sum, npts, f, g )
c---------------------------------------------------------------
c Purpose: take do product of two vectors
c Written: Peter A. Schultz, 14-August-2003 from 2.56(alpha)
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  f(npts),g(npts)
c
      sum = 0.d0
      do  i=1,npts
        sum = sum + f(i)*g(i)
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CROSS
c
c
      subroutine CROSS( a,b, axb )
c---------------------------------------------------------------
c Purpose: Perform cross product axb() = a() X b()
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  a(3),b(3), axb(3)
c
      axb(1) = a(2)*b(3) - a(3)*b(2)
      axb(2) = a(3)*b(1) - a(1)*b(3)
      axb(3) = a(1)*b(2) - a(2)*b(1)
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GETFIELDS
c
c
      subroutine GETFIELDS( IRD,IWR, IGRDFL, jobnm, npts, ierr,
     $ rho0,drho,dves,tmp )
c---------------------------------------------------------------
c Purpose: extract grid data from Quest output grid files
c Written: Peter A. Schultz, 12-August-2003 from 2.56(alpha)
c Revision history:
c  none
c---------------------------------------------------------------
c
c The contents of the reference grid files ".grd0":
c  (1) rho0 - the reference density on the grid
c  (2) rho_core - core density for non-linear correction
c                 NB: this is absent if no atoms have core corrections
c  (3) vxc0 - sum of atomic cutoff xc-potentials
c  (4) exc0 - sum of atomic cutoff xc energy densities
c
c The content of the scf grid files:
c   NO spin:
c    (1) drho (2) dves (3) vslo
c   WITH spin:
c    (1) drho[up] (2) drho[dn] (3) dves (4) vslo[up] (5) vslo[dn]
c   where
c    "drho" is the delta-density: rho(scf)-rho0
c    "dves" is delta-Ves: the electostatic potential from "drho"
c    "vslo" is the sum of the delta-Ves and delta-Vxc, the total
c           delta potential on the grid
c     and drho and vslo are potentially spin-resolved.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input: job/file name:
      CHARACTER  jobnm*(*)
      LOGICAL    do_spin
c Output:
      DIMENSION  rho0(npts), drho(npts), dves(npts)
c Scratch:
      DIMENSION  tmp(npts)
c Local:
      CHARACTER  opt*(1)
      CHARACTER  filenm*(80)
      DATA       zero,one / 0.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Preset flag to "failure"
      ierr = 1
c
c Need to know if spin to decode grid output files:
   10 continue
      write(IWR,*) 'Spin polarized (y-es) or not (n-o)?'
      read(IRD,'(a)') opt
      if( opt.eq.'y' .or. opt.eq.'Y' )then
        do_spin = .true.
      elseif( opt.eq.'n' .or. opt.eq.'N' )then
        do_spin = .false.
      else
        goto 10
      endif
c
c First, open the reference grid file ".grd0", ...
      filenm = jobnm // '.grd0'
      OPEN( unit=IGRDFL,file=filenm,status='OLD',form='UNFORMATTED',
     $ err = 1301 )
c
c      ... and retrieve the reference density:
      call READBIG( IGRDFL, npts, rho0 )
      CLOSE( unit=IGRDFL )
c
c Second, open the scf grid file ".grd" ...
      filenm = jobnm // '.grd'
      OPEN( unit=IGRDFL,file=filenm,status='OLD',form='UNFORMATTED',
     $ err = 1301 )
c
c      .. and retrieve delta-density:
      call READBIG( IGRDFL, npts, drho )
      if( do_spin )then
c       Add the down-spin delta density into the up-spin
        call READBIG( IGRDFL, npts, tmp )
        do  i=1,npts
          drho(i) = drho(i) + tmp(i)
        enddo
      endif
c
c      ... and retrieve the electrostatic potential:
      call READBIG( IGRDFL, npts, dves )
c
c      ... and that is all we want:
      CLOSE( unit=IGRDFL )
c
c Set flag to "success"
      ierr = 0
      goto 999
c
 1301 continue
      call STRPARS( filenm, 1,80, j1,jn )
      write(IWR,*) '***** ERROR: missing file:',filenm(j1:jn)
c
  999 continue
      RETURN
      END     
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKZERO
c
c
      subroutine MKZERO( n, a )
c---------------------------------------------------------------
c Purpose: zero out a double precision real array
c
c Revision history:
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      DIMENSION  a(n)
      DATA  zero / 0.d0 /
c
      do  i=1,n
        a(i) = zero
      enddo
c
      RETURN
      END
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
c  none
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
c    BACKBIG( ifile, lenvec )      - skip over vector on unit
c    READSKP( ifile, lenvec )      - backspace over vector on unit
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
