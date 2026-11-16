c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BFSPLIT
c
c
      subroutine BFSPLIT( ntyp,nshld,nald, cbfsh,occsh,occij,
     $ numshl,lshel,nala,ala,cala )
c---------------------------------------------------------------
c Purpose: construct split valence basis from full shell basis
c
c Written: Peter A. Schultz, 27-July-1998, for v2.23 (split basis)
c
c Revision history:
c  21Jun01-PAS/2.48: replace STOPs
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c  On input, the bf and occ shells are copies of each other
c  On output:
c    1) The bf shells (numshl-cala) are in split basis form.
c    2) cbfsh gives orig shells from split basis shells:
c                sh(j) = SUM[i]( cbfsh(i,j)*bf(i) )
c    3) occij(i,j) = matrix coefs for i.gt.j (shell-wise)
c                  = effective z for i.le.j
c
c  I assume that bf coefficients are all non-zero on input.
c  Otherwise the splitting algorithm *might* have a problem.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  (lmax=4)
c output/input arrays:
      DIMENSION  cbfsh(nshld,nshld,ntyp)
      DIMENSION  occsh(nshld,ntyp), occij(nshld,nshld,ntyp)
      DIMENSION  numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
c local declarations:
      DIMENSION  numl(0:lmax)
      DATA  zero,one,two / 0.d0,1.d0,2.d0 /, vsmall / 1.d-6 /
c
c >>>> EXECUTABLE CODE
c
      call MKZERO( nshld*nshld*ntyp, cbfsh )
      call MKZERO( nshld*nshld*ntyp, occij )
c
c  Cycle over all atom types:
c
      do 1000 ityp=1,ntyp
        nshl = numshl(ityp)
c
c  Start with identity, i.e. bf=occ shells:
        do  ishl=1,nshl
          cbfsh(ishl,ishl,ityp) = one
        enddo
c
c  Count up number of shells in each angular momentum:
        do  l=0,lmax
          numl(l) = 0
        enddo
        do  ishl=1,nshl
          lshl = lshel(ishl,ityp)
          if( lshl.gt.lmax .or. lshl.lt.0 ) call STOPXERR( 'split-l' )
          numl(lshl) = numl(lshl) + 1
        enddo
c
c Loop over shells of a given angular momentum
        do 500 l=0,lmax
c
          alnext = zero
          numspl = numl(l)
c
c  Return here to split out next shell:
c
  100     numspl = numspl - 1
c
c  Setting this cutoff at 1 would restrict kind of splitting allowed
          if( numspl.lt.0 ) goto 500
c
c  Find next shell/alpha to be split off
          ishspl = 0
          aldone = alnext
          do 200 ishl=1,nshl
            if( lshel(ishl,ityp).ne.l ) goto 200
c
c  This shell has right l, find its range of alphas:
            ispl = 1
            almins = ala(1,ishl,ityp)
            almaxs = ala(1,ishl,ityp)
            do 110 ial=2,nala(ishl,ityp)
              ali = ala(ial,ishl,ityp)
              if( ali.lt.almins )then
                ispl = ial
                almins = ali
              endif
              if( ali.gt.almaxs ) almaxs = ali
  110       continue
c
            if( ishspl.eq.0 )then
              if( almins.gt.aldone )then
c  First shell with l, load up
                ishspl = ishl
                ialspl = ispl
                alnext = almins
                almax = almaxs
              endif
            elseif( almins.gt.aldone .and. almins.lt.alnext )then
c  Has smaller unsplit alpha than previous, load it
              ishspl = ishl
              ialspl = ispl
              alnext = almins
              almax = almaxs
            elseif( almins.eq.alnext .and. almaxs.lt.almax )then
c  Has same splitting alpha, but smaller range, load it instead
              ishspl = ishl
              ialspl = ispl
              almax = almaxs
            endif
c
  200     continue
          if( ishspl.eq.0 ) call STOPXERR( 'split-n' )
c
c  Project out this splitting shell from other shells:
c
          do 300 ishl=1,nshl
            if( lshel(ishl,ityp).ne.l ) goto 300
            if( ishl.eq.ishspl ) goto 300
c
c  Right, at this point I have some choices:
c   A0) project out of *every* shell having this alpha
c   A1) project only those where projected number of non-zero
c       alphas does not increase (i.e., at most one new alpha)
c   A2) project only those where all splitting shell gaussians
c       are also present in target shell
c      (Hence *strictly* reducing in number of gaussians)
c   B0) restrict to those where target has same min alpha
c       (with or without conditions 1 and 2)
c
c  I am going to choose A2.  The bookkeeping is simpler, and
c  gives best set of min alphas for all but pathological cases
c  (internal alphas could be better optimized), and easier
c  formulas for getting bf:occ shell projections.
c
c  This shell has right l, check the projecting alphas:
            jalspl = 0
            nalj = nala(ishl,ityp)
            do 220 ial=1,nala(ishspl,ityp)
              ali = ala(ial,ishspl,ityp)
              do 210 jal=1,nalj
                alj = ala(jal,ishl,ityp)
                if( alj.eq.ali )then
                  if( ial.eq.ialspl ) jalspl = jal
                  goto 220
                endif
  210         continue
c  We did not find a match for a splitting shell gaussian:
              goto 300
  220       continue
c
c  I have assumed a non-zero coefficient.  Idiot-check this here:
            if( ABS(cala(ialspl,ishspl,ityp)).lt.vsmall )
     $       call STOPXERR( 'split-c' )
c
c  Do the arithmetic:
            coef = cala(jalspl,ishl,ityp)/cala(ialspl,ishspl,ityp)
c
c  Keep track of how to back-transform:
            do 240 jshl=1,nshl
              cbfsh(ishspl,jshl,ityp) = cbfsh(ishspl,jshl,ityp) +
     $                             coef*cbfsh(ishl,jshl,ityp)
  240       continue
c
c  And do the work on the coefficients:
            do 260 ial=1,nala(ishspl,ityp)
              ali = ala(ial,ishspl,ityp)
              cali = coef*cala(ial,ishspl,ityp)
              do 250 jal=1,nalj
                alj = ala(jal,ishl,ityp)
                if( alj.eq.ali )then
                  cala(jal,ishl,ityp) = cala(jal,ishl,ityp) - cali
                  goto 260
                endif
  250         continue
  260       continue
            cala(jalspl,ishl,ityp) = zero
c
c  And compact the result (i.e. clean out alphas with zero weight):
            do 280 jal=nalj,1,-1
              if( ABS(cala(jal,ishl,ityp)) .lt. vsmall )then
                nala(ishl,ityp) = nala(ishl,ityp) - 1
                do 270 jj=jal+1,nalj
                  ala(jj-1,ishl,ityp) = ala(jj,ishl,ityp)
                  cala(jj-1,ishl,ityp) = cala(jj,ishl,ityp)
  270           continue
              endif
  280       continue
            if( nala(ishl,ityp).eq.0 ) call STOPXERR( 'split-x' )
c
c  Close shell loop:
c
  300      continue
c
c  This splitting shell/alpha now taken out of all relevant orbs.
c  Go find the next splitting shell/alpha for this angular momentum:
c
           goto 100
c
c  Close angular momentum loop:
c
  500   continue
c
c  Normalize the split basis functions:
c
        do 550 ishl=1,nshl
          lshl = lshel(ishl,ityp)
          nal  = nala(ishl,ityp)
c
          call GOVLP( s,lshl, nal,ala(1,ishl,ityp),cala(1,ishl,ityp),
     $                        nal,ala(1,ishl,ityp),cala(1,ishl,ityp) )
c
c  Both normalize the basis functions *and* scale transformation matrix:
          scale = SQRT(s)
          do  jshl=1,nshl
            cbfsh(ishl,jshl,ityp) = scale*cbfsh(ishl,jshl,ityp)
          enddo
          snorm = one/scale
          do  ial=1,nal
            cala(ial,ishl,ityp) = snorm*cala(ial,ishl,ityp)
          enddo
  550   continue
c
c  Now construct reference occupation numbers for split basis:
        do 650 kshl=1,nshl
          lshl = lshel(kshl,ityp)
          norbsh = 2*lshl + 1
          occfac = occsh(kshl,ityp)/DBLE(norbsh)
          if( occfac.ne.zero )then
            do  jshl=1,nshl
              do  ishl=1,nshl
                occij(ishl,jshl,ityp) = occij(ishl,jshl,ityp) +
     $           occfac*cbfsh(ishl,kshl,ityp)*cbfsh(jshl,kshl,ityp)
              enddo
            enddo
          endif
  650   continue
c
c  Now get z numbers of cross terms, double, and put into
c  ishl<jshl triangle (for use in "vlocmii")
        do 700 jshl=2,nshl
          l = lshel(jshl,ityp)
          nalj = nala(jshl,ityp)
          do 690 ishl=1,jshl-1
            if( occij(ishl,jshl,ityp).ne.zero )then
              nali = nala(ishl,ityp)
              call GOVLP(s,l,nali,ala(1,ishl,ityp),cala(1,ishl,ityp),
     $                       nalj,ala(1,jshl,ityp),cala(1,jshl,ityp) )
              occij(ishl,jshl,ityp) = two*s*occij(ishl,jshl,ityp)
            endif
  690     continue
  700   continue
c
c  Close loop over types:
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
