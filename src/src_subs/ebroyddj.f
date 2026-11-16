c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EBROYDDJ
c
c
      subroutine EBROYDDJ( IWR, sufx, iter,ithist, nhist, blendfac,
     $ chg, nvec, nvecloc, membroy,
     $ vv1,vv2,vv3,vv4, wbroy )
c---------------------------------------------------------------
c Purpose: blend guess vectors using Johnson modified Broyden
c          [D.D. Johnson, PRB 38,12807 (1988)]
c
c Written: Peter A. Schultz
c          Based on original code by D.D. Johnson
c          Via code modified by P.J. Feibelman
c
c Revision history:
c  26Apr12-ACP&PAS/2.63: for kparallel, lift distribution out of here
c  20Jun07-PAS/2.60: merge serial and tp
c  13Feb07-PAS/2.60: make tp/serial more consistent;
c                    extract file-handling to separate file;
c                    in-core broyden history; purge Bcast ges
c   6Jul06-APT/    : merged task parallel and image parallel
c   3Nov04-PAS/2.59-tp: parallelized (split form)
c   5Feb03-PAS/2.57a: fix bug in diagnostic output (prune),
c     and fix bug that sets number of Broyden files.
c  11Sep03-PAS/2.57: push itmaxd into broyden
c   8Jan03-PAS/2.55: file size limits, limited history added
c   8Mar02-PAS/2.52: literal constants to d0
c   6Dec01-PAS/2.51: file mgmt
c  11Oct01-PAS/2.50: purge 'inquire' usage
c  21Dec99-PAS/2.42: cosmetic clean-up
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  22Sep97-PAS/2.21: pass dv in vv1, take out convrg tests,
c                    use for geo opt, and rename.
c   9May95-PAS/2.12: shrink memory usage 2X
c---------------------------------------------------------------
c
c Notes:
c  This routine protects itself like a module.
c  It is initialized with input of iter.le.1, but afterwards
c  maintains its state internally.  It should be able to
c  protect itself from all eventualities short of malicious
c  tampering with its files.
c
c Variables:
c  sufx  = character string/tag identifying optimized quantity
c  iter  = .le.1 - restart broyden; else - *ignored*
c          (broyden step numbers are maintained internally)
c  ithist = number of iterations in current history
c  nhist = if>0, max number of steps used in Broyden history
c  blendfac = mixing parameter for first step
c          (recommend: electronic: ~.1->.3; geometric: ~0.5->6.0)
c  chg   = for weighting set, measure of goodnees of this set
c          (e.g. max force/change, or rms, etc)
c
c  nvec  = length of total vector to be optimized
c  nvecloc  = length of vector on local process
c  vv2() = v(in) [input] -> new guess v() [output]
c  vv1() = dv()  [input]
c  vv3,vv4 = scratch space of dimensioned length
c  wbroy() = in-core saved broyden history
c  membroy = amount of space available for in-core
c  
c
c  itmaxd = maximum number of active Broyden steps
c
c Notes:
c  u() and dft() = johnson's u() and df(transpose)
c    These are continually updated and all iterations stored on
c    unit ibroyfl to avoid memory cost of holding entire jacobian
c  dv()  = v(out) - v(in) OR d(F)/d(v())
c  ddv() = dv(m+1) - dv(m)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.0d0, one=1.0d0 )
c
c Input/Output:
      CHARACTER  sufx*(*)
      DIMENSION  vv1(*),vv2(*)
      DIMENSION  wbroy(*)
c Scratch:
      DIMENSION  vv3(*),vv4(*)
c
c Local declarations:
c    itmaxd = max steps in history for blend (30 is likely overkill)
      PARAMETER  ( itmaxd = 30 )
      DIMENSION  a(itmaxd,itmaxd),b(itmaxd,itmaxd),cm(itmaxd)
      DIMENSION  d(itmaxd,itmaxd),w(itmaxd)
      DIMENSION  td(itmaxd),ad(itmaxd),bd(itmaxd)
c    maxfl_broy = maximum files to divide Broyden data into
      PARAMETER  ( maxfl_broy = 20 )
      DIMENSION  ibrfls(maxfl_broy)
      CHARACTER  broyfl*80
c
      DATA  lstout / 0 /
c
c >>>> EXECUTABLE CODE:
c
c Routine now uses input size for distribution
      nveclocl = nvecloc
c
c Obtain parallel environment:
      call MPNODES( nprocl )
      call MPNODE0( masterl )
      call MPNODE( iprocl )
      call MPCOMM( icomml )
c
c Open binary file for broyden:
      broyfl = 'broy0_' // sufx
      call FLOPENB( ibrfl0, broyfl )
c
      iterb = iter
c     Decide on storage (disk=0, or in-core>0) of broyden history.
      incorebr = 1
      call BROYINIT( incorebr, iterb,ithist,nhist,nveclocl,membroy )
      broyfl = 'broy0v_' // sufx
C      call FLOPENB( ibrfl0v, broyfl )
      call BROY0OPEN( ibrfl0v, broyfl )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c        Broyden blend: first iteration is simple blend
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( iter .le. 1 .or. iterb .eq. 0 )then
c       If either signaled to (iter.le.1), or lacking history (iterb=0),
c       then restart blend
c
        blendbr = blendfac
        write(IWR,'(a,a,a,f9.6)') '>>>>> Broyden: first step (', sufx,
     $      ') simple blend ratio=', blendbr
        iters = 1
        iterbm1 = 0
        REWIND( unit=ibrfl0 )
        write(ibrfl0)  blendbr, iters,iterbm1, nprocl
c
c       Store first vin/vv2 and dv/vv1:
        call BROY0RWND( ibrfl0v )
        ithist = 1
        call BROY0PUT( ibrfl0v, nveclocl, vv2, wbroy )
        call BROY0PUT( ibrfl0v, nveclocl, vv1, wbroy )
c
c       On first iteration, simple blend to get new vector.
        do  n=1,nveclocl
          vv2(n) = vv2(n) + blendbr*vv1(n)
        enddo
c
        goto 900
c
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c        Broyden blending method for later iterations
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write(IWR,'(a,a)') '>>>>> Broyden blend to get next step - ',sufx
c
c  w0 = weighting factor for the zeroth iteration
      w0 = 0.01d0
      REWIND( unit=ibrfl0 )
      read(ibrfl0,end=1301,err=1301)  blendbr, itersm1,iterbm2, nprocbr
      if( nprocbr .ne. nprocl ) call STOPXERR( 'BROY = #procs changed' )
c
c  Retrieve vv3<=previous vin ; vv4<=previous dv
      call BROY0RWND( ibrfl0v )
      call BROY0GET( ibrfl0v, nveclocl, vv3    , wbroy )
      call BROY0GET( ibrfl0v, nveclocl, vv4    , wbroy )
c
      if( iterbm2.gt.0 .and. iterbm2.lt.itmaxd )then
        read(ibrfl0)  ( ( a(jt,it),jt=1,iterbm2 ), it=1,iterbm2 )
        read(ibrfl0)  ( w(it), it=1,iterbm2 )
      endif
c
c Check limits on blend counts:
      itmaxbr = itmaxd
      if( nhist.gt.0 .and. nhist.lt.itmaxd ) itmaxbr = nhist
      if( iterbm2.eq.itmaxbr )then
c       Drop old iterations in limited-history blend
        write(IWR,'(a)') ' Broyden: limit-history, dropping oldest'
        do  jt=2,iterbm2
          w(jt-1) = w(jt)
          do it=2,iterbm2
            a(it-1,jt-1) = a(it,jt)
          enddo
        enddo
        iterbm2 = iterbm2 - 1
      elseif( iterbm2.gt.itmaxbr )then
        write(IWR,*) 'BROYDEN: more iterations than dimensioned for'
        write(IWR,*) 'This event indicates a coding error, sys failure'
        call STOPXERR( 'BROY dimensions' )
      endif
c Set up blend counts:
c  iter (as input) = total number of steps in overall scf attempt
c  iters = total number of steps in this broyden scf
c  iterb = number of steps in active broyden history
      iters   = itersm1 + 1
      iterbm1 = iterbm2 + 1
      iterb   = iterbm1 + 1
      itprune = iters - iterb
c
c Have vin and dv of current iteration, write to ibrfl0 for use later
      REWIND( unit=ibrfl0 )
      write(ibrfl0)  blendbr, iters,iterbm1, nprocl
c
c  Store current vv2=vin ; vv1=dv
      call BROY0RWND( ibrfl0v )
      call BROY0PUT( ibrfl0v, nveclocl, vv2    , wbroy )
      call BROY0PUT( ibrfl0v, nveclocl, vv1    , wbroy )
c
c  vv1=current dv ; vv2=current vin ; vv3=previous vin ; vv4=previous dv
c  vv3=>delta-vin ; vv4=>ddv
      do  n=1,nveclocl
        vv3(n) = vv2(n) - vv3(n)
        vv4(n) = vv1(n) - vv4(n)
      enddo
c
c For i-th iter, dfnorm = | dv[i] - dv[i-1] |, for normalization
      dfnorm = zero
      do  n=1,nveclocl
        dfnorm = dfnorm + vv4(n)*vv4(n)
      enddo
c
      if( nprocl .gt. 1 )then
c       Reduce and brodcast dfnorm
        lenred = 1
        call MPREDUC8( masterl, lenred, dfnorm,dfnormtp, icomml )
C        call MPI_Reduce( dfnorm,dfnormtp,1,MPI_DOUBLE_PRECISION,
C     $   MPI_SUM, masterl, icomml, ierr )
        dfnorm = dfnormtp
        call MPBCAST8( masterl, lenred, dfnorm, icomml )
C        call MPI_Bcast ( dfnorm, 1, MPI_DOUBLE_PRECISION,
C     $   masterl, icomml, ierr )
      endif
      dfnorm = SQRT( dfnorm )
      if( lstout.gt.0 ) write(IWR,'(a,f16.8)')  ' dfnorm=',dfnorm
c
      fac2 = one / dfnorm
      fac1 = blendbr*fac2
c  vv3=>u ; vv4=>dft
      do  n=1,nveclocl
        vv3(n) = fac1*vv4(n) + fac2*vv3(n)
        vv4(n) = fac2*vv4(n)
      enddo
c
c*********** calculation of coefficient matrices *************
c***********    and the sum for corrections      *************
c
c recall: a(i,j) is a symmetric matrix
c       : b(i,j) is the inverse of [ w0**2 I + a ]
c
c Initialize history storage:
      call BROY2OPEN( itersm1,iterb,itmaxbr, sufx,
     $ ibroyfl, nveclocl )
c
c Rewind history storage:
      call BROY2RWND( ibroyfl )
c
      if( iterb .gt. 2 )then
c
        if( itprune .gt. 0 )then
c         Skip itprune steps, to begin of active Broyden history:
          call BROY2SKIP( itprune, ibroyfl, nveclocl )
        endif
c
        do  it=1,iterbm2
c         Get dft[it] in vv2(), skipping u[it]
          call BROY2GETD( ibroyfl, nveclocl, vv2, wbroy )
c
          aij = zero
          cmj = zero
          do  n=1,nveclocl
            aij = aij + vv2(n)*vv4(n)
            cmj = cmj + vv2(n)*vv1(n)
          enddo
          a(iterbm1,it) = aij
          a(it,iterbm1) = aij
          cm(it)       = cmj
        enddo
      endif
c
      cmj = zero
      aij = zero
      do  n=1,nveclocl
        cmj = cmj + vv4(n)*vv1(n)
        aij = aij + vv4(n)*vv4(n)
      enddo
      cm(iterbm1)        = cmj
      a(iterbm1,iterbm1) = aij
c
      if( nprocl .gt. 1 )then
c       Reduce/broadcast a/cm matrices (bd() is free for scratch)
        call MPREDUC8( masterl, iterbm1, a(1,iterbm1),bd, icomml )
C        call MPI_Reduce( a(1,iterbm1),bd, iterbm1, MPI_DOUBLE_PRECISION,
C     $   MPI_SUM, masterl, icomml, ierr )
        call DCOPY( iterbm1, bd,1, a(1,iterbm1),1 )
        call MPBCAST8( masterl, iterbm1, a(1,iterbm1), icomml )
C        call MPI_Bcast ( a(1,iterbm1), iterbm1, MPI_DOUBLE_PRECISION,
C     $   masterl, icomml, ierr )
        do  it=1,iterbm1
          a(iterbm1,it) = a(it,iterbm1)
        enddo
c
        call MPREDUC8( masterl, iterbm1, cm,bd, icomml )
C        call MPI_Reduce( cm,bd, iterbm1, MPI_DOUBLE_PRECISION,
C     $   MPI_SUM, masterl, icomml, ierr)
        call DCOPY( iterbm1, bd,1, cm,1 )
        call MPBCAST8( masterl, iterbm1, cm, icomml )
C        call MPI_Bcast ( cm, iterbm1, MPI_DOUBLE_PRECISION,
C     $   masterl, icomml, ierr)
      endif
c
c  Store latest vv3=u ; vv4=dft
      call BROY2PUT( ibroyfl, nveclocl, vv3  ,vv4  , wbroy )
c
c Weighting factors for each iteration set to one over sqrt of max error
c This need not be the case, e.g., could be one over rms instead.
      wtmp = zero
      if( chg.gt. 1.0d-09 ) wtmp = 2.0d0*SQRT( 0.010d0/chg )
      if( wtmp.lt.one ) wtmp = one
      w(iterbm1) = wtmp
      if( lstout.gt.0 ) write(IWR,'(a,f20.6)') 'weighting set=',wtmp
c
      write(ibrfl0)  ( ( a(jt,it), jt=1,iterbm1 ), it=1,iterbm1 )
      write(ibrfl0)  ( w(it), it=1,iterbm1 )
c
c Set up and calculate beta matrix
      do  it=1,iterbm1
        do  jt=1,iterbm1
          d(jt,it) = a(jt,it)*w(jt)*w(it)
          b(jt,it) = zero
        enddo
        b(it,it) = one
        d(it,it) = w0**2 + d(it,it)
      enddo
c
      call INVERT( d,b,iterbm1, td,ad,bd, itmaxd )
c
c Compute input vector for new iteration
c
      REWIND( unit=ibrfl0 )
      read(ibrfl0)  blendbr, itm1_t,itm2_t, nprocbr
      if( nprocbr .ne. nprocl ) call STOPXERR( 'BROY - #procs changed' )
c
      call BROY0RWND( ibrfl0v )
      call BROY0GET( ibrfl0v, nveclocl, vv2    , wbroy )
      do  n=1,nveclocl
        vv2(n) = vv2(n) + blendbr*vv1(n)
      enddo
c
      call BROY2RWND( ibroyfl )
c
      if( itprune .gt. 0 )then
c       Skip itprune steps, to beginning of active Broyden history:
        call BROY2SKIP( itprune, ibroyfl, nveclocl )
      endif
c
      do  it=1,iterbm1
c       Get u[it] in vv3(), skipping dft[it]
        call BROY2GETU( ibroyfl, nveclocl, vv3, wbroy )
c
        gmi = zero
        do  jt=1,iterbm1
          gmi = gmi + cm(jt)*b(jt,it)*w(jt)
        enddo
        do  n=1,nveclocl
          vv2(n) = vv2(n) - gmi*vv3(n)*w(it)
        enddo
      enddo
c
c Close the Broyden history file(s):
C      call FLCLOSE( ibroyfl )
      call BROY2CLOSE( ibroyfl )
c
      goto 900
c
 1301 continue
      call STOPXERR( 'BROYDDJ - broyden files missing/mangled' )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c        Clean up and get out:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  900 continue
c
      call FLCLOSE( ibrfl0 )
      call BROY0CLOSE( ibrfl0v )
c
      write(IWR,'(a,i6)') '<<<<< Broyden step complete, scf step=',iter
      if( iter.ne.iters )then
        write(IWR,'(a,i6)') '<<<<< Broyden blend step=',iters
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
c End:
