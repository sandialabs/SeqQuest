c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RH1ORB
c
c
      subroutine RH1ORB( irh1orb,
     $ natm,ntyp,nshld,nald, n1r,n2r,n3r,nptr,nlat,maxwf,
     $ itypa, numshl,lshel,nala,ala,cala,  ratm,hh,rlat,
     $ lmx1ctr,almnv1c,  mngr3,mxgr3,ngr,
     $ r, dfn,rsq,gsum, wf )
c---------------------------------------------------------------
c Purpose: computes on the mesh the products of the orbitals
c          of each atom and summing these products over the
c          lattice vectors.  Needed for the subtraction of the
c          "1-atom terms" from the density in RHO1CTR
c
c Written: Peter A. Schultz, based on PJF's original routine
c
c Revision history:
c  24Sep08-PAS/2.62: fix unit variable in this unused routine
c   7Mar02-PAS/2.52: ang constants extracted
c  21Jun01-PAS/2.48: replace STOPs
c  19Dec99-PAS/2.42: clean-up; read/write big records
c  10Dec98-PAS/2.29: promote source to explicit double precision
c   1May97-PAS/2.17: naming conventions, on/off criterion,
c                    new grid coord rtn
c  15Dec95-PAS/2.14: shell-by-shell grid ranges, cleanup
c  17Sep93-PAS: more stringent cutoffs installed, some opt.
c   5Apr90-PAS: made more efficient and smaller storage
c---------------------------------------------------------------
c
c   On input, the coordinate mesh is assumed to lie in r()
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  ratm(3,natm),hh(3,3),rlat(3,nlat)
      DIMENSION  lmx1ctr(ntyp),almnv1c(ntyp)
      DIMENSION  mngr3(*),mxgr3(*),ngr(*)
      DIMENSION  r(3,nptr)
c scratch arrays:
      DIMENSION  dfn(5,nptr),rsq(nptr),gsum(nptr)
      DIMENSION  wf(*)
c local declarations:
      DIMENSION  r0(3),ratom(3),rlast(3)
c
      DIMENSION  vnorm(5)
      DATA  zero,half, three / 0.d0,0.5d0, 3.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Get angular normalizations:
      call VNORMS( vnorm )
c
      argmx = half*cutgaus
      do  jr=1,3
        rlast(jr) = zero
      enddo
c
c  Compute grid coordinates
c
      r0(1) = zero
      r0(2) = zero
      r0(3) = zero
      call MKRBOX( hh,r0, n1r,n2r,n3r, r )
c
c    Center loop:
c
      i1next = 1
      do 1000 iatm=1,natm
        ityp = itypa(iatm)
        nshl = numshl(ityp)
c I need this for grid indexing arrays:
        iorb1 = i1next
        do  ishl=1,nshl
          i1next = i1next + 2*lshel(ishl,ityp) + 1
        enddo
        if( lmx1ctr(ityp).lt.0 ) goto 1000
c
        al1cmin = almnv1c(ityp)
c
        almin = ala(1,1,ityp)
        almax = ala(1,1,ityp)
        do  ishl=1,nshl
          do  ial=1,nala(ishl,ityp)
            if( ala(ial,ishl,ityp).lt.almin ) almin = ala(ial,ishl,ityp)
            if( ala(ial,ishl,ityp).gt.almax ) almax = ala(ial,ishl,ityp)
          enddo
        enddo
c       rsqtop = cutgrid/(almin+almin)
c find range-limiting gaussian (smallest of (largest alf of 1c pair)):
        al1c = almax
        do  ishl=1,nshl
          do  jshl=1,ishl
            do  ial=1,nala(ishl,ityp)
              ali = ala(ial,ishl,ityp)
              do  jal=1,nala(jshl,ityp)
                alj = ala(jal,jshl,ityp)
                als = ali + alj
                if( als.ge.al1cmin )then
                  algt = MAX(ali,alj)
                  al1c = MIN(al1c,algt)
                endif
              enddo
            enddo
          enddo
        enddo
        twoalc = al1c + almin
        rsqtop = cutgrid/twoalc
c
        nlatm = nlat
        iorbat0 = 0
        do 850 ishl=1,nshl
          lshli = lshel(ishl,ityp)
          nali = nala(ishl,ityp)
          norbshi = 2*lshli + 1
c
          nrshi = ngr(iorb1+iorbat0)
          imin3sh = mngr3(iorb1+iorbat0)
c
          alimax = ala(nali,ishl,ityp)
          do  ial=1,nali
            ali = ala(ial,ishl,ityp)
            if( ali.gt.alimax ) alimax = ali
          enddo
c
          jorbat0 = 0
          do 800 jshl=1,ishl
            lshlj = lshel(jshl,ityp)
            nalj = nala(jshl,ityp)
            norbshj = 2*lshlj + 1
            lsum = lshli + lshlj
c
            nrshj = ngr(iorb1+jorbat0)
            jmin3sh = mngr3(iorb1+jorbat0)
c
            nrs = MIN( nrshi,nrshj )
            ijmin3 = MAX( imin3sh,jmin3sh )
            nr0 = n1r*n2r*(ijmin3-1)
c
            aljmax = ala(nalj,jshl,ityp)
            do  jal=1,nalj
              alj = ala(jal,jshl,ityp)
              if( alj.gt.aljmax )  aljmax = alj
            enddo
            alsmx = alimax + aljmax
c
            if( lsum.gt.lmx1ctr(ityp) .or. alsmx.lt.al1cmin)goto 800
c
            al1c = almax
            do  ial=1,nala(ishl,ityp)
              ali = ala(ial,ishl,ityp)
              do  jal=1,nala(jshl,ityp)
                alj = ala(jal,jshl,ityp)
                als = ali + alj
                if( als.ge.al1cmin )then
                  algt = MAX(ali,alj)
                  al1c = MIN(al1c,algt)
                endif
              enddo
            enddo
            twoalc = al1c + almin
            r2cutsh = cutgrid/twoalc
c
c this stuff below used to ensure only i<j calculated in lower loops:
            norbij = norbshi*norbshj
            norbj0 = norbshj
            joff = 0
            if( ishl.eq.jshl )then
              norbij = (norbij + norbshi)/2
              norbj0 = 0
              joff = 1
            endif
            nrset = norbij*nrs
            if( nrset.gt.2*maxwf ) call STOPXERR( 'rh1orb1' )
c
            do  nr=1,nrset
              wf(nr) = zero
            enddo
c
            nlatmx = 1
            do 900 ilat=1,nlatm
c
              do  jr=1,3
                ratom(jr) = ratm(jr,iatm) + rlat(jr,ilat)
              enddo
c
              do  nr=1,nrs
                r(1,nr0+nr) = r(1,nr0+nr) - ratom(1) + rlast(1)
                r(2,nr0+nr) = r(2,nr0+nr) - ratom(2) + rlast(2)
                r(3,nr0+nr) = r(3,nr0+nr) - ratom(3) + rlast(3)
                rsq(nr) = r(1,nr0+nr)**2 + r(2,nr0+nr)**2
     $                  + r(3,nr0+nr)**2
                gsum(nr) = zero
              enddo
c
              do  jr=1,3
                rlast(jr) = ratom(jr)
              enddo
c
              rsqmin = rsq(1)
              rsqmax = rsq(1)
              do  nr=1,nrs
                if( rsq(nr).lt.rsqmin ) rsqmin = rsq(nr)
                if( rsq(nr).gt.rsqmax ) rsqmax = rsq(nr)
              enddo
              if( rsqmin.lt.rsqtop .and. ilat.gt.nlatmx ) nlatmx = ilat
              if( rsqmin.gt.r2cutsh ) goto 900
c
              if( lshli.ge.2 .or. lshlj.ge.2 )then
                ang4 = vnorm(3)/vnorm(4)
                ang5 = vnorm(3)/vnorm(5)
                do  nr=1,nrs
                  dfn(1,nr) = r(1,nr0+nr)*r(2,nr0+nr)
                  dfn(2,nr) = r(2,nr0+nr)*r(3,nr0+nr)
                  dfn(3,nr) = r(3,nr0+nr)*r(1,nr0+nr)
                  dfn(4,nr) = (r(1,nr0+nr)**2 - r(2,nr0+nr)**2)*ang4
                  dfn(5,nr) = (three*r(3,nr0+nr)**2 - rsq(nr))*ang5
                enddo
              endif
c
              do 290 ial=1,nali
                ali = ala(ial,ishl,ityp)
                cfi = cala(ial,ishl,ityp)/vnorm(lshli+1)
                do 280 jal=1,nalj
                  alj = ala(jal,jshl,ityp)
                  alsum = ali + alj
                  if( alsum.lt.al1cmin ) goto 280
c
                  algt = ali
                  if( alj.gt.ali ) algt = alj
                  twoalgt = algt + algt
                  twoalc = algt + almin
                  rsqcut = cutgrid/twoalc
c  Skip if alphas too large
                  if( rsqmin.gt.rsqcut ) goto 280
c
                  cij = cfi*cala(jal,jshl,ityp)/vnorm(lshlj+1)
                  if( twoalgt*rsqmax.lt.argmx )then
                    do  nr=1,nrs
                      gsum(nr) = gsum(nr) + cij*EXP( -alsum*rsq(nr) )
                    enddo
                  else
                    do  nr=1,nrs
                      if( rsq(nr).lt.rsqcut )
     $                gsum(nr) = gsum(nr) + cij*EXP( -alsum*rsq(nr) )
                    enddo
                  endif
  280           continue
  290         continue
c
              if( lshli+lshlj.eq.0 )then
c
                do  nr=1,nrs
                  wf(nr) = wf(nr) + gsum(nr)
                enddo
c
              elseif( lshli+lshlj.eq.1 )then
c
                do  nr=1,nrs
                  wf(      nr) = wf(      nr) + gsum(nr)*r(1,nr0+nr)
                  wf(  nrs+nr) = wf(  nrs+nr) + gsum(nr)*r(2,nr0+nr)
                  wf(2*nrs+nr) = wf(2*nrs+nr) + gsum(nr)*r(3,nr0+nr)
                enddo
c
              elseif( lshli+lshlj.eq.2.and.lshli*lshlj.eq.0 )then
c
                do  nr=1,nrs
                  wf(      nr) = wf(      nr) + gsum(nr)*dfn(1,nr)
                  wf(  nrs+nr) = wf(  nrs+nr) + gsum(nr)*dfn(2,nr)
                  wf(2*nrs+nr) = wf(2*nrs+nr) + gsum(nr)*dfn(3,nr)
                  wf(3*nrs+nr) = wf(3*nrs+nr) + gsum(nr)*dfn(4,nr)
                  wf(4*nrs+nr) = wf(4*nrs+nr) + gsum(nr)*dfn(5,nr)
                enddo
c
              elseif( lshli.eq.1.and.lshlj.eq.1 )then
c
                norbshj = norbj0
                nr00 = 0
                do  iorbsh=1,norbshi
                  norbshj = norbshj + joff
                  do  jorbsh=1,norbshj
                    do  nr=1,nrs
                      wf(nr00+nr) = wf(nr00+nr) + gsum(nr)*
     $                 r(iorbsh,nr0+nr)*r(jorbsh,nr0+nr)
                    enddo
                    nr00 = nr00 + nrs
                  enddo
                enddo
c
              elseif( lshli+lshlj.eq.3.and.lshli*lshlj.eq.2 )then
c
                if(lshli.eq.1)then
                  n1 = 5
                  n2 = 1
                else
                  n1 = 1
                  n2 = 3
                endif
                do  k1=1,3
                  do  k2=1,5
                    nr00 = (n1*(k1-1)+n2*(k2-1))*nrs
                    do  nr=1,nrs
                      wf(nr00+nr) = wf(nr00+nr) + gsum(nr)*
     $                 r(k1,nr0+nr) * dfn(k2,nr)
                    enddo
                  enddo
                enddo
c
              elseif( lshli.eq.2.and.lshlj.eq.2 )then
c
                norbshj = norbj0
                nr00 = 0
                do  iorbsh=1,norbshi
                  norbshj = norbshj + joff
                  do  jorbsh=1,norbshj
                    do  nr=1,nrs
                      wf(nr00+nr) = wf(nr00+nr) + gsum(nr)*
     $                 dfn(iorbsh,nr)*dfn(jorbsh,nr)
                    enddo
                    nr00 = nr00 + nrs
                  enddo
                enddo
c
              else
                call STOPXERR( 'rh1orb2' )
              endif
c
c    Close lattice vector loop:
c
  900       continue
            nlatm = nlatmx
c
c    Write out the results:
c
            nrset0 = 0
            do  ijprod=1,norbij
C             write(irh1orb)  (wf(nr),nr=nrset0+1,nrset0+nrs)
              call WRITBIG( irh1orb, nrs, wf(nrset0+1) )
              nrset0 = nrset0 + nrs
            enddo
c
c    Reset the coordinate mesh:
c
            do  nr=1,nrs
              r(1,nr0+nr) = r(1,nr0+nr) + rlast(1)
              r(2,nr0+nr) = r(2,nr0+nr) + rlast(2)
              r(3,nr0+nr) = r(3,nr0+nr) + rlast(3)
            enddo
            do  jr=1,3
              rlast(jr) = zero
            enddo
c
c    Close radial function loops:
c
  800     jorbat0 = jorbat0 + norbshj
  850   iorbat0 = iorbat0 + norbshi
c
c    Close center loop:
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
