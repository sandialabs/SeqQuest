c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DEFBOX
c
c
      subroutine DEFBOX( IWR, ndim, n1r,n2r,n3r,
     $ nbox,mxrbox, ndbox )
c---------------------------------------------------------------
c Purpose: carve up grid unit cell into smaller boxes
c
c Written: Peter A. Schultz,  3-June-1997, for v2.20 (boxes)
c
c Revision history:
c  21Jun01-PAS/2.48: replace STOPs
c  18Dec00-PAS/2.46: tune to limit box sizes to 216 or less.
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c Notes:
c   What we would like to do here is divide the grid points
c   in a unit cell into boxes, in such a way as to:
c     1.  keep grid calcns within cache, requiring smaller boxes.
c     2.  keep vectorization, requiring larger boxes
c     3.  keep overhead low, also requiring larger boxes
c     4.  get maximum benefit of cutoffs, requiring smaller
c         boxes in those directions that cutoffs can take effect.
c     5.  keep boxes relatively uniformly sized/shaped.
c   What I adopt as a default is boxes that have edges of length
c   four or five grid points, defining boxes with 64-125 total
c   grid points in them.  However, in lower dimension calcs,
c   e.g. slabs, this could be wasteful - you would like smaller
c   divisions in the z-directions.  So the following code does
c   some ad hoc division for breaking up the unit cell that
c   hopefully is more efficient.  A 6/5 rather than a 5/4 scheme
c   might be more efficient, but will depend on machine cache
c   sizes, and number of orbitals and k-points.
c
c   The difficulty with using reduced dimension slicing is that
c   main code may not be smart enough to take advantage of the
c   thinner slicing in one direction, and you end up with a more
c   expensive calculation rather than cheaper one.
c   "reduce_d" controls whether effective dimension changed.
c   I may just purge this special case treatment in the future.
c
c   I impose a limit of nrboxd=216 for the biggest possible box.
c   DO NOT CHANGE this value from 216!!!!!!
c   This enables me to make viable (but checked) assumptions in
c   main routine about amount of work space needed in some rtns
c   that are called from main (e.g. "grid0xc").
c   The scheme implemented here will produce boxes up to 216,
c   for either default box edges of 5 or 6.  This rtn has been
c   tested thoroughly for both default=5 and default=6.
c
c Variables:
c  ndim         - input dimension (periodicity) of problem
c  n1r,n2r,n3r  - input grid dimensions of unit cell
c  nbox         - output number of boxes
c  mxrbox       - output size of largest box
c  ndbox(2,2,3) - output breakup of cell into boxes
c  Along direction "iu":
c    (1,1,iu),(2,1,iu)  is number of boxes with ...
c    (1,2,iu),(2,2,iu)   ... points along the side of the box
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c local parameters:
c   Restrict biggest box to 216 pts (DO NOT CHANGE!):
      PARAMETER  (nrboxd=216)
c   Set default box edge length (must be five or six!!!)
c    idefr1=5 -> (4:5)*(4:5)*(4:5)=[ 64:125] pts/box default
c    idefr1=6 -> (5:6)*(5:6)*(5:6)=[125:216] pts/box default
      PARAMETER  (idefr1=5)
c
c output array:
      DIMENSION  ndbox(2,2,3)
c
c local declarations:
      DIMENSION  nrv(3),idrv(3)
      LOGICAL    reduce_d
      LOGICAL    booboo
c
c >>>> EXECUTABLE CODE:
c
c Our box breakdown only works for: default=5 or default=6.
      if( idefr1.lt.5 .or. idefr1.gt.6 ) call STOPXERR( 'defbox-d' )
c
c Set up default box dimensions:
c
      nrv(1) = n1r
      nrv(2) = n2r
      nrv(3) = n3r
      do  idim=1,3
c       Set (large) box edge to default value ...
        idrv1 = idefr1
c        ... and fix small cell dimension cases:
        nvr = nrv(idim)
        if( nvr.eq.19 .or. nvr.eq.14 .or. nvr.eq.13 .or. nvr.eq.9 )then
c         This fixes bad cases specific for default=6
          idrv1 = 5
        elseif( nvr.eq.11 )then
c         This fixes bad case specific for default=5
          idrv1 = 6
        elseif( nvr.eq.8 .or. nvr.eq.7 )then
          idrv1 = 4
        elseif( nvr.lt.7 )then
c         Should *never* have cell this small, but for completeness:
          idrv1 = nvr
        endif
        idrv(idim) = idrv1
      enddo
      idra1 = idrv(1)
      idrb1 = idrv(2)
      idrc1 = idrv(3)
c
c The reduce-dimensions special cases can be skipped (reduce_d=F)
c
      reduce_d = .true.
      if( reduce_d )then
c
c       Determine effective system dimension:
c
        mdim = ndim
c       Try to reduce dimension:
        if( mdim.eq.3 .and. n3r.gt.64 ) mdim = 2
        if( mdim.eq.2 .and. n2r.gt.64 ) mdim = 1
        if( mdim.eq.1 .and. n1r.gt.64 ) mdim = 0
c       Try to increase dimension (this to catch pathologies):
        if( mdim.eq.1 .and. n2r.lt.32 ) mdim = 2
        if( mdim.eq.1 .and. n3r.lt.32 ) mdim = 3
        if( mdim.eq.2 .and. n3r.lt.32 ) mdim = 3
c
        if( mdim.eq.1 )then
c
c         Have a 1-D unit cell (n1r.le.64)
c         Presume: (n2r.gt.64).and.(n3r.gt.64), but unknown
c
          if( n1r.lt.12 )then
c           n1r=( 1:11) -> ( 1:11)x4x4=[ 16:176]
            idra1 = n1r
            idrb1 = 4
            idrc1 = 4
          elseif( n1r.lt.20 )then
c           n1r=(12:19) -> (12:19)x3x3=[108:171]
            idra1 = n1r
            idrb1 = 3
            idrc1 = 3
          elseif( n1r.lt.55 )then
c           n1r=(20:54) -> (20:54)x2x2=[ 80:216]
            idra1 = n1r
            idrb1 = 2
            idrc1 = 2
          elseif( n1r.le.64 )then
c           n1r=(55:64) -> (28:32)x2x2=[112:128]
            idra1 = (n1r-1)/2 + 1
            idrb1 = 2
            idrc1 = 2
          else
c           cannot be
            call STOPXERR( 'defbox1d' )
          endif
c
        elseif( mdim.eq.2 )then
c
c         Have 2-D unit cell (n2r.le.64)
c         Presume: (n3r.gt.64).and.(n1r.le.64), but unknown
c         Want: thin slabs.  Get: no more than 14 on a side
c         with 14**2=196 in each layer of slice.
c         n1r,n2r    box edge
c         (57+  ) ->  (6:7)      -> ( 6: 7)
c         (43:56) -> (n-1)/4 + 1 -> (11:14)
c         (29:42) -> (n-1)/3 + 1 -> (10:14)
c         (15:28) -> (n-1)/2 + 1 -> ( 8:14)
c         ( 1:14) ->   n         -> ( 1:14)
c
c         Divide up x-vector
          if( n1r.gt.56 )then
            idra1 = 7
          elseif( n1r.gt.42 )then
            idra1 = (n1r-1)/4 + 1
          elseif( n1r.gt.28 )then
            idra1 = (n1r-1)/3 + 1
          elseif( n1r.gt.14 )then
            idra1 = (n1r-1)/2 + 1
          else
            idra1 = n1r
          endif
c         Divide up y-vector
          if( n2r.gt.56 )then
            idrb1 = 7
          elseif( n2r.gt.42 )then
            idrb1 = (n2r-1)/4 + 1
          elseif( n2r.gt.28 )then
            idrb1 = (n2r-1)/3 + 1
          elseif( n2r.gt.14 )then
            idrb1 = (n2r-1)/2 + 1
          else
            idrb1 = n2r
          endif
c
          n12d = idra1*idrb1
c         Divide up z-vector, with n12d=(1:14)*(1:14)=(1:196)
c
          if( n12d.lt.26 )then
c           ( 1:25)*5 = [  5:125]
            idrc1 = 5
          elseif( n12d.lt.41 )then
c           (26:40)*4 = [104:160]
            idrc1 = 4
          elseif( n12d.lt.56 )then
c           (41:55)*3 = [123:165]
            idrc1 = 3
          elseif( n12d.lt.91 )then
c           (56:90)*2 = [112:180]
            idrc1 = 2
          else
c           (91:196)*1 = [91:196]
            idrc1 = 1
          endif
c
        endif
c
      endif
c
c  Box edge lengths determined, now get number of boxes ...
c
c   ... first define smaller box edge length
      idra2 = idra1 - 1
      idrb2 = idrb1 - 1
      idrc2 = idrc1 - 1
c   ... get number of short lengths
      nboxa2 = ( (n1r+idra2) / idra1 )*idra1 - n1r
      nboxb2 = ( (n2r+idrb2) / idrb1 )*idrb1 - n2r
      nboxc2 = ( (n3r+idrc2) / idrc1 )*idrc1 - n3r
c   ... and get number of longer lengths
      nboxa1 = (n1r - nboxa2*idra2)/idra1
      nboxb1 = (n2r - nboxb2*idrb2)/idrb1
      nboxc1 = (n3r - nboxc2*idrc2)/idrc1
c   ... shift if have zero larger boxes
      if( nboxa1.eq.0 )then
        nboxa1 = nboxa2
        nboxa2 = 0
        idra1 = idra2
        idra2 = idra1 - 1
      endif
      if( nboxb1.eq.0 )then
        nboxb1 = nboxb2
        nboxb2 = 0
        idrb1 = idrb2
        idrb2 = idrb1 - 1
      endif
      if( nboxc1.eq.0 )then
        nboxc1 = nboxc2
        nboxc2 = 0
        idrc1 = idrc2
        idrc2 = idrc1 - 1
      endif
c
c Do all sorts of idiot double checks ...
c   ... we only get failure if code design is messed up, default ok:
      booboo = .false.
c   ... see if we get nonsense breakdown:
      if( nboxa1.lt.1 .or. nboxa2.lt.0 .or. idra1.lt.1 .or.
     $    nboxb1.lt.1 .or. nboxb2.lt.0 .or. idrb1.lt.1 .or.
     $    nboxc1.lt.1 .or. nboxc2.lt.0 .or. idrc1.lt.1 ) booboo =.true.
c   ... and check to be sure it hangs together:
      m1r = nboxa1*idra1 + nboxa2*idra2
      m2r = nboxb1*idrb1 + nboxb2*idrb2
      m3r = nboxc1*idrc1 + nboxc2*idrc2
      if( m1r.ne.n1r .or. m2r.ne.n2r .or. m3r.ne.n3r ) booboo = .true.
c   ... and double-check that box size conforms to design:
      mxrbox = idra1*idrb1*idrc1
      if( mxrbox.gt.nrboxd ) booboo = .true.
c
      if( booboo )then
c       Code is designed to work right, yet somehow we missed
c       something.  DO NOT simply change dimensioning (nrboxd)
c       or hack code to make this error go away.  This branch
c       is only hit if above code has a very bad error in logic
c       somewhere, or if defaults changed to unacceptable range.
        write(IWR,*) '>>>>> coding failure! (this is VERY bad)'
        write(IWR,*) 'Report grid dimensions to code architect:'
        write(IWR,*) 'n1r,n2r,n3r=',n1r,n2r,n3r
        write(IWR,*) ' nboxa1,idra1;[b];[c]=',
     $   nboxa1,idra1,nboxb1,idrb1,nboxc1,idrc1
        write(IWR,*) ' nboxa2,idra2;[b];[c]=',
     $   nboxa2,idra2,nboxb2,idrb2,nboxc2,idrc2
        write(IWR,*) 'yielding n1r,n2r,n3r=',m1r,m2r,m3r
        write(IWR,*) 'box size,max=',mxrbox,nrboxd
        call STOPXERR( 'defbox - error defining box sizes' )
      endif
c
c  Load up the unit cell box break-up
c
      nbox = (nboxa1+nboxa2)*(nboxb1+nboxb2)*(nboxc1+nboxc2)
      ndbox(1,1,1) = nboxa1
      ndbox(2,1,1) = nboxa2
      ndbox(1,2,1) = idra1
      ndbox(2,2,1) = idra2
      ndbox(1,1,2) = nboxb1
      ndbox(2,1,2) = nboxb2
      ndbox(1,2,2) = idrb1
      ndbox(2,2,2) = idrb2
      ndbox(1,1,3) = nboxc1
      ndbox(2,1,3) = nboxc2
      ndbox(1,2,3) = idrc1
      ndbox(2,2,3) = idrc2
c
c    That's all Folks!
c
      RETURN
      END
