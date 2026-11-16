c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BOXRNG
c
c
      subroutine BOXRNG( IWR,lvlout, maxrec,    mngr3,mxgr3,ngr,
     $ ndim, n1r,n2r,n3r, natm,ntyp, nshld,nald,
     $ itypa,numshl,lshel,nala,norba,ala, alamin, ratm,rprim,
     $ hh, ndbox,nbox,nrec,inrec,noboxs,inboxs, boxrad0,
     $ wkbox )
c---------------------------------------------------------------
c Purpose: to determine extent of orbitals on grid, and record
c   which shells are in which boxes, radius of boxes.
c
c Written: Peter A. Schultz,  3-June-1997, for 2.20 (boxes)
c
c Revision history:
c  27Feb07-PAS/2.60: adjust output flags
c   2Dec03-PAS/2.58: add atm/shl record when atom not in cell
c   8Mar02-PAS/2.52: clean unused variables & literal const d0
c  21Jun01-PAS/2.48: replace STOPs
c  24May01-PAS/2.47: transfer out wf storage requirements code
c   2May01-PAS/2.47: remove disk storage of grid orbitals
c  11Apr00-PAS/2.43: pass out length of grid wf storage
c  21Dec99-PAS/2.42: cosmetic clean-up
c  22Jul99-PAS/2.37: filter parallelpiped of boxes, chg diagnostic
c   1May99-PAS/2.35: pass max single box mem.
c   9Jan99-PAS/2.30: redo symmetry scheme
c   9Nov98-PAS/2.28: single precision grid orbitals/wf4
c  10Sep98-PAS/2.26: cosmetic - documentation
c  14May98-PAS/2.22: redid
c---------------------------------------------------------------
c
c Notes:
c  The filter of box occupancy within parallelpiped is done via
c  a crude box-center to atom distance scheme.  This should be
c  good enough.  You'll miss a couple of boxes you should keep
c  where a few points at a corner will be within range of an
c  atom, and you will include some boxes where all the points
c  will be out of range, but not all that many of either.
c  You could do a finer point-wise check for the suspect boxes,
c  but it is not worth the effort and cost.  26Jul99-PAS
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c input arrays:
      DIMENSION  itypa(natm),numshl(ntyp),lshel(nshld,ntyp),
     $ nala(nshld,ntyp),ala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp)
      DIMENSION  ratm(3,natm), rprim(3,3), hh(3,3)
      DIMENSION  ndbox(2,2,3),noboxs(nbox),inboxs(*)
cpas:1 these are for non-standard usage
      DIMENSION  mngr3(*),mxgr3(*),ngr(*)
c scratch arrays:
      DIMENSION  wkbox(4,nbox)
c local declarations:
      DIMENSION  rlnn(3,27), r0box(3)
      DIMENSION  g(3,3),rg(3),nrv(3),rnv(3), v(3),irngs(2,3)
      DIMENSION  ibox1(3),iboxn(3)
      LOGICAL    booboo
c
      DATA  zero,half,one / 0.d0,0.5d0,1.d0 /
      DATA  lstout / 1 /
c
c >>>> EXECUTABLE CODE:
c
      booboo = .false.
      lstout = lvlout
      if( lstout.gt.0 ) write(IWR,*)
c
c  Find minimum radius of a sphere that contains largest box:
      do  idim=1,3
        g(1,idim) = (ndbox(1,2,idim)-1)*hh(1,idim)
        g(2,idim) = (ndbox(1,2,idim)-1)*hh(2,idim)
        g(3,idim) = (ndbox(1,2,idim)-1)*hh(3,idim)
      enddo
      rg(1) = half*( g(1,1) + g(1,2) + g(1,3) )
      rg(2) = half*( g(2,1) + g(2,2) + g(2,3) )
      rg(3) = half*( g(3,1) + g(3,2) + g(3,3) )
      sqpp = zero
      sqpm = zero
      sqmp = zero
      sqmm = zero
      do  j=1,3
        sqpp = sqpp +   rg(j)**2
        sqpm = sqpm + ( rg(j) - g(j,3) )**2
        sqmp = sqmp + ( rg(j) - g(j,2) )**2
        sqmm = sqmm + ( rg(j) - g(j,2) - g(j,3) )**2
      enddo
      boxrad0 = MAX( sqpp, sqpm, sqmp, sqmm )
      boxrad0 = SQRT( boxrad0 )
      if( lstout.gt.4 ) write(IWR,*) 'BOXRNG: boxrad=',boxrad0
c
c  Actually, using this true extreme radius without point-wise
c  filter will catch boxes that have very few, even zero, points
c  within range of the atom, unless you do a point-wise check.
c  I sacrifice boxes with a few in-range corner points to kick
c  out boxes with zero in-range points.  I do this by using a
c  box radius that is smaller that the extreme.  Assuming a cube
c  and box radius, the radius is the 111-diag, and the smallest
c  radius is a 100-face, with ratio rt3:1.  Try an intermediate
c  ratio such as .75*rt3:1 (~1.3:1).
      boxrad = 0.75*boxrad0
c
c Collect coordinates of box centers
c
      irpot = 0
      do 30 ibox=1,nbox
c
c  Get basic box info
        call INDBOX( ibox, ndbox, iboxa,iboxb,iboxc,
     $   ira0,irb0,irc0, idra,idrb,idrc )
        nrbox = idra*idrb*idrc
c
        iv0 = irpot
        irpot = irpot + nrbox
c
c  Generate coordinates
        call BOXR( hh,ira0,irb0,irc0, idra,idrb,idrc, r0box,v )
c  Load coordinates of center of this box:
        wkbox(1,ibox) = v(1)
        wkbox(2,ibox) = v(2)
        wkbox(3,ibox) = v(3)
c        call MKRBOX( hh,r0box, idra,idrb,idrc, r(1,iv0+1) )
c
   30 continue
c
c Set up nearby lattice vectors
c
      nla = 0
      nlb = 0
      nlc = 0
      if( ndim.gt.0 ) nla = 1
      if( ndim.gt.1 ) nlb = 1
      if( ndim.gt.2 ) nlc = 1
c  Put origin first
      il = 1
      rlnn(1,1) = zero
      rlnn(2,1) = zero
      rlnn(3,1) = zero
      do  kc=-nlc,nlc
        do  kb=-nlb,nlb
          do  ka=-nla,nla
            if( ka.ne.0 .or. kb.ne.0 .or. kc.ne.0 )then
              il = il + 1
              rlnn(1,il) = ka*rprim(1,1) + kb*rprim(1,2) + kc*rprim(1,3)
              rlnn(2,il) = ka*rprim(2,1) + kb*rprim(2,2) + kc*rprim(2,3)
              rlnn(3,il) = ka*rprim(3,1) + kb*rprim(3,2) + kc*rprim(3,3)
            endif
          enddo
        enddo
      enddo
      nlnn = il
c
      nrv(1) = n1r
      nrv(2) = n2r
      nrv(3) = n3r
c generate orthogonal projection of grid vectors:
      call CROSS( hh(1,2),hh(1,3), g(1,1) )
      call CROSS( hh(1,3),hh(1,1), g(1,2) )
      call CROSS( hh(1,1),hh(1,2), g(1,3) )
      do  id=1,3
        gmag = SQRT( g(1,id)**2+g(2,id)**2+g(3,id)**2 )
        g(1,id) = g(1,id)/gmag
        g(2,id) = g(2,id)/gmag
        g(3,id) = g(3,id)/gmag
        rnv(id) = DBLE( nrv(id) )
        rg(id) = hh(1,id)*g(1,id) + hh(2,id)*g(2,id) + hh(3,id)*g(3,id)
        if( rg(id).lt.zero )then
          rg(id) = -rg(id)
          g(1,id) = -g(1,id)
          g(2,id) = -g(2,id)
          g(3,id) = -g(3,id)
        endif
      enddo
c
c  Make sure have enough box records
      natshl = 0
      do  iatm=1,natm
        ityp = itypa(iatm)
        natshl = natshl + numshl(ityp)
      enddo
      nrec = (natshl+inrec-1)/inrec
      if( nrec.gt.(maxrec-1) )then
        write(IWR,*) '>>>>> nrec,nbox=',nrec,nbox
        write(IWR,*) '>>>>> need nrecd =',(nrec+1)*nbox
        call STOPXERR( 'box-recs/not enough box records (nrecd)' )
      else
c       Clear out box records
        do  ibox=1,nbox
          noboxs(ibox) = 0
        enddo
        do  ibr=1,nrec*nbox
          inboxs(ibr) = 0
        enddo
      endif
c
 9000 format(/1x,'ATOMIC ORBITAL RANGES:', /1x,
     $   ' Atom Shell  Min-alf     range',14x,'grid range')
 9001 format( 1x,i4,3x,2x,1x,f10.6,f11.6 )
 9002 format( 1x,4x,3x,i2,1x,f10.6,f11.6,3x,
     $ 2('(',i3,'-',i3,')x'),'(',i3,'-',i3,')' )
 9012 format( 2x,4x,3x,i2,1x,f10.6,f11.6)
 9003 format( 1x,4x,3x,2x,1x, 10x ,'nbox = ',i4,' = ',
     $ 2('(',i3,'-',i3,')x'),'(',i3,'-',i3,')' )
c
      if( lstout.gt.0 )then
        write(IWR,*) ' Atom'
        write(IWR,*) ' Type  Shell  Min-alf     range'
        do  ityp=1,ntyp
          almina = alamin(ityp)
          range = SQRT( cutgrid / (almina+almina) )
          write(IWR,9001)  ityp,almina,range
          if( lstout.gt.2 )then
            do  ishl=1,numshl(ityp)
              almins = ala(1,ishl,ityp)
              do  ial=1,nala(ishl,ityp)
                ali = ala(ial,ishl,ityp)
                if( ali.lt.almins ) almins = ali
              enddo
              almin = almina
              range = SQRT( cutgrid / (almin +almins) )
              write(IWR,9012)  ishl,almins,range
            enddo
          endif
        enddo
      endif
c
      if( lstout.gt.3 ) write(IWR,9000)
cpas:1 for non-standard code
      iorb = 0
      nobchk = 0
      nobocc = 0
      nsbchk = 0
      nsbocc = 0
      irec = 0
      iatshl = 0
      inbox = 0
      do 300 iatm=1,natm
        ityp = itypa(iatm)
c
        do  idim=1,3
          v(idim) = ratm(1,iatm)*g(1,idim)
     $            + ratm(2,iatm)*g(2,idim)
     $            + ratm(3,iatm)*g(3,idim)
        enddo
c
        almina = alamin(ityp)
        range = SQRT( cutgrid / (almina+almina) )
c
        if( lstout.gt.3 ) write(IWR,9001)  iatm,almina,range
c
c  Do grid orbitals ranges, shell-by-shell:
c  On the first pass through, create cutoff alpha using longest
c  alpha on atom.  This is a bad hack, but other choices are bad, too.
        do 280 ishl=0,numshl(ityp)
c
          if( ishl.eq.0 )then
c
c  Pass 0: Compute cutoff alphas by box, use longest atom alpha
c
            almins = almina
c
          else
c
c  Otherwise, go ahead and compute ranges
c
            norbsh = 2*lshel(ishl,ityp) + 1
            inbox = inbox*2
            if( inrec*(iatshl/inrec).eq.iatshl )then
              irec = irec + 1
              inbox = 1
            endif
            iatshl = iatshl + 1
c
            almins = ala(1,ishl,ityp)
            do  ial=2,nala(ishl,ityp)
              ali = ala(ial,ishl,ityp)
              if( ali.lt.almins ) almins = ali
            enddo
c
          endif
          almin = almina
          range = SQRT( cutgrid / (almin +almins) )
c
          do 150 idim=1,3
            irdim = nrv(idim)
            rdim = rg(idim)
c  Get orbital range ...
            rmn = (v(idim) - range)/rdim
            rmx = (v(idim) + range)/rdim
            if( rmn.gt.-half )then
              irmn = ( rmn + half + one )
            else
              irmn = -ABS( rmn - half - one )
            endif
            if( rmx.gt.-half )then
              irmx = ( rmx + half + one )
            else
              irmx = -ABS( rmx - half - one )
            endif
c
            if( idim.gt.ndim .and. (irmn.lt.1 .or. irmx.gt.irdim) )then
c             See how bad it is
              ioutdn = MAX( 0, 1 - irmn )
              ioutup = MAX( 0, irmx - irdim )
              routdn = ioutdn*rg(idim)
              routup = ioutup*rg(idim)
  913         format(1x,a,4i6)
  914         format(1x,a,2f12.3)
              write(IWR,913)
     $         '***** WARNING: orbitals extend outside unit cell'
              write(IWR,913)
     $         '>>>>> atom,shell,direction=',iatm,ishl,idim
              write(IWR,913)
     $         '>>>>> Grid: irmn,irmx, limits=',irmn,irmx, 1,irdim
              write(IWR,914)
     $         '>>>>> outside cell (dn/up) by:', routdn,routup
              rngerr = MAX( routup, routdn )
              if( rngerr/range .gt. 0.25d0 )then
                booboo = .true.
              endif
              irmn = MAX( irmn, 1 )
              irmx = MIN( irmx, irdim )
            endif
c
c  ... and convert to range over boxes
            nbox1 = ndbox(1,1,idim)
            nbox2 = ndbox(2,1,idim)
            nboxs = nbox1 + nbox2
            idr1 = ndbox(1,2,idim)
            idr2 = ndbox(2,2,idim)
            if( (irmx-irmn+idr1) .gt. irdim ) then
              ibox1(idim) = 1
              iboxn(idim) = nboxs
              jrmn = 1
              jrmx = irdim
            else
              jrmn = irmn
              if( irmn.lt.1 ) jrmn = irmn + irdim
              if( jrmn.gt.(nbox1*idr1) )then
                ibox1(idim) = nbox1 + (jrmn-nbox1*idr1+idr2-1) / idr2
              else
                ibox1(idim) = (jrmn+idr1-1) / idr1
              endif
              jrmx = irmx
              if( irmx.gt.irdim ) jrmx =  irmx - irdim
              if( jrmx.gt.(nbox1*idr1) )then
                iboxn(idim) = nbox1 + (jrmx-nbox1*idr1+idr2-1) / idr2
              else
                iboxn(idim) = (jrmx+idr1-1) / idr1
              endif
              if( irmn.lt.  1   ) ibox1(idim) = ibox1(idim) - nboxs
              if( irmx.gt.irdim ) iboxn(idim) = iboxn(idim) + nboxs
c    ... check if have all boxes
              if( (iboxn(idim)-ibox1(idim)+1) .ge. nboxs ) then
                ibox1(idim) = 1
                iboxn(idim) = nboxs
              elseif( ibox1(idim).lt.1 )then
                ibox1(idim) = ibox1(idim) + nboxs
                iboxn(idim) = iboxn(idim) + nboxs
              endif
            endif
            irngs(1,idim) = jrmn
            irngs(2,idim) = jrmx
  150     continue
c
c At this point, use range to box to determine cutoff alpha,
c and use this cutoff alpha to determine if shell overlaps box.
c Only want to compute cutoffs once, hence need to branch here
c according to: first pass, compute cutoff alphas; second
c passes for shells all cutoff alphas already computed, load
c box occupancies.
c
          nboxa = ndbox(1,1,1) + ndbox(2,1,1)
          nboxb = ndbox(1,1,2) + ndbox(2,1,2)
          nboxc = ndbox(1,1,3) + ndbox(2,1,3)
          jboxa1 = ibox1(1)
          jboxb1 = ibox1(2)
          jboxc1 = ibox1(3)
          jboxan = iboxn(1)
          jboxbn = iboxn(2)
          jboxcn = iboxn(3)
          nboxat = (jboxan-jboxa1+1)*(jboxbn-jboxb1+1)*(jboxcn-jboxc1+1)
c
          if( ishl.eq.0 )then
c
c  Cycle over boxes potentially within range of any shell of this atom:
            do 173 jboxc=jboxc1,jboxcn
              iboxc = jboxc
              if( jboxc.gt.nboxc ) iboxc = jboxc - nboxc
              do 172 jboxb=jboxb1,jboxbn
                iboxb = jboxb
                if( jboxb.gt.nboxb ) iboxb = jboxb - nboxb
                do 171 jboxa=jboxa1,jboxan
                  iboxa = jboxa
                  if( jboxa.gt.nboxa ) iboxa = jboxa - nboxa
                  ibox = iboxa + nboxa*( iboxb - 1 + nboxb*(iboxc-1) )
c  Compute cutoff alphas for all boxes w/in range of longest atom shell
                  ilmn = 1
                  rsqmn = (ratm(1,iatm)-wkbox(1,ibox))**2
     $                  + (ratm(2,iatm)-wkbox(2,ibox))**2
     $                  + (ratm(3,iatm)-wkbox(3,ibox))**2
                  do 170 il=2,nlnn
                    rsq = (ratm(1,iatm)-wkbox(1,ibox)-rlnn(1,il))**2
     $                  + (ratm(2,iatm)-wkbox(2,ibox)-rlnn(2,il))**2
     $                  + (ratm(3,iatm)-wkbox(3,ibox)-rlnn(3,il))**2
                    if( rsq.lt.rsqmn )then
                      ilmn = il
                      rsqmn = rsq
                    endif
  170             continue
c  We have a range**2 of box, convert to max alpha that overlaps it:
                  rngbox = SQRT( rsqmn ) - boxrad
                  if( rngbox.lt. 1.d-8 ) rngbox = 1.d-8
                  cutalf = cutgrid/(rngbox**2) - almin
                  wkbox(4,ibox) = cutalf
  171           continue
  172         continue
  173       continue
c
          else
c
c  Cycle over boxes potentially within range of this shell ...
            do 183 jboxc=jboxc1,jboxcn
              iboxc = jboxc
              if( jboxc.gt.nboxc ) iboxc = jboxc - nboxc
              do 182 jboxb=jboxb1,jboxbn
                iboxb = jboxb
                if( jboxb.gt.nboxb ) iboxb = jboxb - nboxb
                do 181 jboxa=jboxa1,jboxan
                  iboxa = jboxa
                  if( jboxa.gt.nboxa ) iboxa = jboxa - nboxa
                  ibox = iboxa + nboxa*( iboxb - 1 + nboxb*(iboxc-1) )
c   ... and check shell alpha against cutoff alpha to overlap box:
                  if( almins.lt.wkbox(4,ibox) )then
                    nirec = nrec*(ibox-1) + irec
                    noboxs(ibox) = noboxs(ibox) + norbsh
                    inboxs(nirec) = inboxs(nirec) + inbox
c
                    nsbocc = nsbocc + 1
                    nobocc = nobocc + norbsh
                  endif
c
                    nsbchk = nsbchk + 1
                    nobchk = nobchk + norbsh
  181           continue
  182         continue
  183       continue
c
cpas:start/oldstuff
            imin3 = irngs(1,3)
            imax3 = irngs(2,3)
            nr = n1r*n2r*(imax3-imin3+1)
            do  iorbsh=1,norbsh
              iorb = iorb + 1
              ngr(iorb) = nr
              mngr3(iorb) = imin3
              mxgr3(iorb) = imax3
            enddo
cpas:end/oldstuff
            if( lstout.gt.4 )then
              write(IWR,9002) ishl,almins,range,irngs
              if( lstout.gt.5 )
     $          write(IWR,9003)   nboxat,
     $           jboxa1,jboxan,jboxb1,jboxbn,jboxc1,jboxcn
            endif
c
          endif
c
  280   continue
c
  300 continue
c
      if( lstout.gt.3 )then
        write(IWR,'(a,2i10)')'BOXRNG: shell/orb box chks=',nsbchk,nobchk
        write(IWR,'(a,2i10)')'BOXRNG: shell/orb box occs=',nsbocc,nobocc
      endif
c
      if( booboo )then
        call STOPXERR( 'orb-rng /orbitals range too far outside cell' )
      endif
c
c    That's all Folks!
c
      RETURN
      END
