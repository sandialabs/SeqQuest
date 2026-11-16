c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KGRID
c
c
      subroutine KGRID( ndim, rprim,
     $ iktype, nka,nkb,nkc, nkfull, veck )
c---------------------------------------------------------------
c Purpose: produce a regular grid of k-points in full BZ
c
c Written: Peter A. Schultz, 19-October-1999, for v2.40
c
c Revision history:
c  14Jun05-PAS/2.59: File units
c  18Feb05-PAS/2.59: changes needed for closeshl or defect sampling
c  23Jul03-PAS/2.56: add special hexagonal offset grids
c   2Jul01-PAS/2.48: replace STOPs
c---------------------------------------------------------------
c
c This gets a little messy and makes bunches of assumptions that
c would be easy to violate, because it tries to take care of all
c sorts of special cases simultaneously.  It should work for all
c common cases and might fail if the user gets really cute.
c
c Variables:
c   ndim    - periodic dimensionality
c   rprim() - primitve lattice vectors
c   iktype  - type of kgrid
c           = 0,1,2,3 number of 0-offsets (rather than 1/2 offset)
c           = 20 1x1 hex grid offset to triangle center
c   nkU     - dimension of k-grid in dimensions
c   nkfull  - total number of k-points in grid
c   veck(3,*) - resulting k-point grid
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
c Input: primitive lattice vectors
      DIMENSION  rprim(3,3)
c Output: k-vectors
      DIMENSION  veck(3,*)
c Local declarations:
      DIMENSION  nkx(3), nkt(3)
      DIMENSION  rvec(3,3), gvec(3,3)
      DIMENSION  vk0(3), dk(3,3)
      DATA  zero,one,two,three / 0.d0,1.d0,2.d0,3.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Determine type of k-grid to generate:
      if( iktype.ge.0 .and. iktype.le.3 )then
c       Number of offsets to set to 0 ( rather than 1/2 )
        n0k = iktype
      elseif( iktype.eq.20 )then
c       Hex grid, offset of 1/3 (center of triangle)
        n0k = 2
      else
c       Default: offset to body-center ( all 1/2 )
        n0k = 0
      endif
c
c Load k-dimensions, setting non-periodic directions to gamma:
      nkx(1) = nka
      nkx(2) = nkb
      nkx(3) = nkc
      do  jd=ndim+1,3
        nkx(jd) = 0
      enddo
c
c Check for gamma point calculation:
      if( nka.eq.0 .and. nkb.eq.0 .and. nkc.eq.0 )then
c       Trivial gamma, point, load and exit
        nkfull = 1
        veck(1,1) = zero
        veck(2,1) = zero
        veck(3,1) = zero
        RETURN
      endif
c
c Determine k-space dimensionality (number of non-zero dims)
c   kdim = number of non-zero dimensions
c   mdim = index of last non-zero dimension
      kdim = 0
      do  jd=1,3
        nkt(jd) = nkx(jd)
        if( nkx(jd).ne.0 )then
          kdim = kdim + 1
          mdim = jd
        endif
      enddo
c
c At this point have:
c  ndim = dimension of problem (molecule=0, ... , 3=bulk)
c  kdim = number of non-zero dimensions of k-grid
c  mdim = dimension index of last non-zero dimension of grid
c
c We have choices how to generate the k-grid, corresponding to:
c  First, how to set the effective dimension of the problem.
c   (1) Choose the actual problem dimension: ndim
c       -straightforward
c   (2) Choose the default k-space dimension: mdim
c       - also straightforward, leaving the problem in the
c         "standard" state of periodic dimension first, simply
c         reducing the k-grid dimension from the r-space dimension
c   (3) Choose the number of non-zero k-dims: kdim
c       - perhaps the "true" dimension of the grid and is likely
c         the closest thing to an "ideal" sample.  However, since
c         it is not standard order (mixed periodic/non-periodic)
c         it requires reordering and careful bookkeeping.
c  Second, how to use the reduced dimension:
c   (a) Use full 3-dimensional R-vectors to get G-vec
c       - straightforward
c   (b) Use reduced dimension space, i.e., assuming
c       non-periodic dimension are normal to periodic dimensions.
c       - some modest effort to created "reduced" R-vectors.
c
c We will choose (1)(b). (Specified by setting kdim=ndim)
        kdim = ndim
c   We choose (b), because non-periodic directions should NOT
c   be used in generating the reciprocal space.  Note this only
c   has an effect if non-periodic R-vectors are not orthogonal
c   to periodic vectors, which should not usually happen.
c   We choose (1) because otherwise we run into problems in
c   other places in the code, i.e. cell-optimization, which
c   must do transformations of the G-vectors, but will not
c   have knowledge of the reduced dimension and thereby
c   cause errors unless this information (kdim,mdim) is
c   propagated.  Sending this information around is prone
c   to coding errors and therefore we will skip it. <sigh>
c
c Standard usage:
      call DCOPY( 9, rprim,1, rvec,1 )
c
      if( mdim.gt.kdim )then
c       Have non-periodics between/before periodics, reorder vecs
        if( n0k .ne. 0 )then
          call STOPXERR( 'kgrid dims and offsets incompatible' )
        endif
        jdk = 0
        jdx = kdim
        do  jd=1,3
          if( nkt(jd).ne.0 )then
            jdk = jdk + 1
            kd = jdk
          else
            jdx = jdx + 1
            kd = jdx
          endif
          nkx(kd) = nkt(jd)
          rvec(1,kd) = rprim(1,jd)
          rvec(2,kd) = rprim(2,jd)
          rvec(3,kd) = rprim(3,jd)
        enddo
      endif
c
c Produce reduced primitive reciprocal lattice vectors:
      call KCELL( kdim, rvec, gvec )
c
c Define k-point intervals:
      do  jd=1,3
        if( nkx(jd).ne.0 )then
          nkt(jd) = nkx(jd)
        else
          nkt(jd) = 1
        endif
        do  i=1,3
          dk(i,jd) = gvec(i,jd)/ DBLE( nkt(jd) )
        enddo
      enddo
c
c Define k-grid offset:
      vk0(1) = zero
      vk0(2) = zero
      vk0(3) = zero
      if( iktype.ge.0 .and. iktype.lt.3 )then
c       Simple 1/2 (body-center) offset
        do  jd=1,kdim
          if( jd.gt.n0k .and. nkx(jd).ne.0 )then
            do  i=1,3
              vk0(i) = vk0(i) + dk(i,jd)/two
            enddo
          endif
        enddo
      elseif( iktype.eq.20 )then
c       2D Hex 1/3 offset
c       First, another idiot check:
        aklen = SQRT( dk(1,1)**2 + dk(2,1)**2 + dk(3,1)**2 )
        bklen = SQRT( dk(1,2)**2 + dk(2,2)**2 + dk(3,2)**2 )
        adotb = dk(1,1)*dk(1,2) + dk(2,1)*dk(2,2) + dk(3,1)*dk(3,2)
        diflen = ABS( one - bklen/aklen )
        difdot = ABS( one - two*ABS( adotb ) / (aklen*bklen) )
        if( diflen.gt.(1.d-6) .or. difdot.gt.(1.d-6) )then
          call FLGETIERR( IERRFL )
          write(IERRFL,'(a,3f16.8)') 'dkA=',(dk(i,1),i=1,3)
          write(IERRFL,'(a,3f16.8)') 'dkB=',(dk(i,2),i=1,3)
          write(IERRFL,'(a,3f16.8)') 'A,B,A.B=',aklen,bklen,adotb
          write(IERRFL,'(a,3f16.8)') 'scaled=',one,bklen/aklen,
     $     adotb/(aklen*bklen)
          write(IERRFL,'(a,3f16.8)') 'diflen,difdot=',diflen,difdot
          call STOPXERR( 'KGRID: hex k-grid vectors are not hex' )
        endif
        bsign = one
        if( adotb.lt.zero ) bsign = -one
        do  i=1,3
          vk0(i) = ( dk(i,1) + bsign*dk(i,2) ) / three
        enddo
        jd = 3
        if(  jd.gt.n0k .and. nkx(jd).ne.0 )then
          do  i=1,3
            vk0(i) = vk0(i) + dk(i,3)/two
          enddo
        endif
      endif
c
c Generate resulting k-grid:
      nkfull = nkt(1)*nkt(2)*nkt(3)
      call MKRBOX( dk, vk0, nkt(1),nkt(2),nkt(3), veck )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KGAMMA
c
c
      subroutine KGAMMA( IWR, ikgamma, nk,nkd,wtk, veck,veckg )
c---------------------------------------------------------------
c Purpose: if missing, add the gamma-pt to the k-sample
c Written: P.A. Schultz, 14Feb05, for 2.59/defect k-sampling
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0, one=1.d0 )
c
c Input/Output: k-vectors
      DIMENSION  wtk(*)
      DIMENSION  veck(3,*), veckg(3,*)
c
c >>>> EXECUTABLE CODE:
c
c First, try to find the gamma pt in the list:
      do 100 ik=1,nk
        veckmag = veck(1,ik)**2 + veck(2,ik)**2 + veck(3,ik)**3
        if( veckmag .eq. zero )then
          ikgamma = ik
          if( IWR.gt.0 )
     $    write(IWR,*) 'KGAMMA: found gamma-pt at k=',ikgam
          goto 900
        endif
  100 continue
c
      if( IWR.gt.0 )
     $write(IWR,*) 'KGAMMA: will add gamma-pt to k-sample'
      nk = nk + 1
      if( nk .gt. nkd ) call STOPXERR( 'too many kpts' )
c
c Load zero-weighted gamma pt at end of list
      ikgam = nk
c     Denote special point outsider status by sending
c      defect k-point as negative
      ikgamma = -ikgam
      wtk(ikgam) = one
      do  i=1,3
        veck(i,ikgam) = zero
        veckg(i,ikgam) = zero
      enddo
c
  900 continue
c
      RETURN
c
c    That's all Folks!
c
      END
