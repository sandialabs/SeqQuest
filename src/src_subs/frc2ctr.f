c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FRC2CTR
c
c
      subroutine FRC2CTR( ndim,
     $ nafrc,iforce, frcs,frct, strs,strt,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, nlat,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba,
     $ ratm, rlat,coskr,sinkr,
     $ dmat,emat, frctp )
c---------------------------------------------------------------
c Purpose: compute contribution to force and stress from the
c          kinetic and overlap matrices
c
c Written: Peter A. Schultz, based on original PJF routine stress2c
c
c Revision history:
c  15Oct08-PAS/2.62: stripe-parallel
c  21Jun07-PAS/2.60: task-parallelism
c   7Mar02-PAS/2.52: ang constants extracted
c  23Jul01-PAS/2.49: compact dmat, spin-polarized dft
c  21Dec99-PAS/2.42: cosmetic clean-up
c  27Sep99-PAS/2.40: full three-dimensional stress tensor;
c                    and upper triangle; renamed to frc2ctr
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  17Dec98-PAS/2.30: clean force/stress arrays
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  13Aug98-PAS/2.24: cutoffs, names, restructure
c  16Jul98-PAS/2.23: fix label problem (thx D.Raczkowski)
c  28Apr97-PAS/2.17: naming conventions
c  10Feb96-PAS/2.13: standardize structure
c  14Apr94-PAS/2.10: cleanup and optimize
c  14Mar94-PAS/2.04: cleanup, match setup routines
c---------------------------------------------------------------
c
c This routine combines the kinetic energy and overlap matrix calcs
c  of the bulk setup phase with the corresponding impurity code
c  force calcs to give the kinetic energy and overlap contributions
c  to the stress of the slab, i.e., the derivative of the k.e.
c  and overlap trace with respect to expansion of the lattice in
c  the surface plane, as well as the corresponding normal forces on
c  various inequivalent layers. Since the calculations are based on
c  the Hellmann-Feynman theorem, the derivative of the density matrix
c  matrix is not included. Terms involving this derivative add to
c  zero, according to the theorem.
c
c  Calculates kinetic/overlap derivatives analytically
c  This routine is an adaptation of "enij"
c
c  Radial wave functions are assumed to have been normalized (this is
c  done in routine "orbbulk", where the cala coefficients have been
c  rescaled). Normalization of the angular wave functions is done in
c  the present routine via the division by the quantities "vnorm"
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
      PARAMETER  ( zero=0.d0, two=2.d0, half=0.5d0 )
c
c Output(/input) arrays:
      DIMENSION  iforce(*)
      DIMENSION  frcs(3,natm),frct(3,natm), strs(3,3),strt(3,3)
c Input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat),coskr(nk,nlat),sinkr(nk,nlat)
      DIMENSION  dmat(nk,*),emat(nk,*)
c Scratch/parallel:
      DIMENSION  frctp(3,natm), strtp(3,3)
c
c Local declarations:
      LOGICAL     diagshl
      DIMENSION   rij(3), p(3), ci(10,9),cj(10,9)
      DIMENSION   engrs(3),engrt(3), sgr(3),tgr(3)
c
      DIMENSION   vnorm(9)
c
c >>>> EXECUTABLE CODE:
c
C      write(6,*) 'FDEV/FRC2CTR: stripe-parallel'
c Get angular normalizations:
      call VNORMS9( vnorm )
      argmx = half*cutgaus
      emaxf = cutgwfn
      emaxsl = cutslow
c
      fourpi = pi + pi + pi + pi
c
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
c    Initialize stress and force
c
      call MKZERO( 3*natm, frcs )
      call MKZERO( 3*natm, frct )
      call MKZERO( 9, strs )
      call MKZERO( 9, strt )
c
c  Determine size of full striped matrix
      call MATTOTSTR( matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
c
      nmatstr = ncplx*nk*matstr
c
c  Figure out how much belongs on local processor:
      call MATLOCSTR( iproc,nprocs,nprocstr,
     $ myatm1, myatm2, myj0, maxmatpr, mymatpr,
     $ matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
c
      nmatpr = mymatpr*ncplx*nk
c
c     Atom j loop:
c
      jjnext = 0
      jnext = myj0
      do 2000 jatm=myatm1,myatm2
        jtyp = itypa(jatm)
        nshlj = numshl(jtyp)
        norbatj = norba(jtyp)
        jstart = jnext
        jnext = jstart + norbatj
c
        njstripe = norb - jstart
        matjstr = norbatj*njstripe
        jjstart = jjnext
        jjnext = jjstart + matjstr*ncplx
Cc     
Cc       Each processor will do every nprocs-th atom in this loop.
C        if( MOD( jatm, nprocs ) .ne. iproc ) goto 2000
c
        alaminj = alamin(jtyp)
        call CINIT( cj )
c     
c       Atom j shell loop:
c
        jorbatnxt = 0
        j = jstart
        do 1800 jshl=1,nshlj
          nalj = nala(jshl,jtyp)
          lshlj = lshel(jshl,jtyp)
          norbshj = 2*lshlj + 1
          lj0 = lshlj*lshlj
c
          j0 = j
          jorbat0 = jorbatnxt
          jorbatnxt = jorbat0 + norbshj
c
          alminj = ala(1,jshl,jtyp)
          almaxj = ala(nalj,jshl,jtyp)
          do  jal=1,nalj
            alj = ala(jal,jshl,jtyp)
            if( alj.lt.alminj ) alminj = alj
            if( alj.gt.almaxj ) almaxj = alj
          enddo
          alaminj = alminj
c
          emaxij = emaxf
          if( alminj.lt.alfast ) emaxij = emaxsl
c
c         Atom i loop:
c
          ijnext = jjstart
          inext = jstart
          do 1000 iatm=jatm,natm
            ityp = itypa(iatm)
            nshli = numshl(ityp)
            norbati = norba(ityp)
            istart = inext
            inext = istart + norbati
c
            matij = norbati*norbatj
            ijstart = ijnext
            ijnext = ijstart + matij*ncplx
c
c           Following helps ensure only i.ge.j done:
            iorbatn0 = 0
            ishl1 = 1
            if( iatm.eq.jatm )then
              ishl1 = jshl
              do  ishl=1,ishl1-1
                lshli = lshel(ishl,ityp)
                norbshi = 2*lshli + 1
                istart = istart + norbshi
                iorbat0 = iorbatn0
                iorbatn0 = iorbat0 + norbshi
              enddo
            endif
c
            alamini = alamin(ityp)
            emax = emaxij
            rijmax = emax*(alamini+alaminj) / (alamini*alaminj)
c
c           Lattice vector loop for Atom i:
c
            do 900 ilat=1,nlat
c
              rijsq = zero
              do  ir=1,3
                rij(ir) =  ratm(ir,jatm) - ratm(ir,iatm) - rlat(ir,ilat)
                rijsq = rijsq + rij(ir)*rij(ir)
              enddo
              if( rijsq.gt.rijmax ) goto 900
c
              call CINIT( ci )
              call COFR( cj,rij )
c         
c             Atom i shell loop:
c 
              iorbatnxt = iorbatn0
              i = istart
              do 800 ishl=ishl1,nshli
                nali = nala(ishl,ityp)
                lshli = lshel(ishl,ityp)
                norbshi = 2*lshli + 1
                li0 = lshli*lshli
c
                i0 = i
                iorbat0 = iorbatnxt
                iorbatnxt = iorbat0 + norbshi
c
                diagshl = .false.
                if( i0.eq.j0 ) diagshl = .true.
c
                almini = ala(1,ishl,ityp)
                almaxi = ala(nali,ishl,ityp)
                do  ial=1,nali
                  ali = ala(ial,ishl,ityp)
                  if( ali.lt.almini ) almini = ali
                  if( ali.gt.almaxi ) almaxi = ali
                enddo
c
                emax = emaxf
                if( almaxi.lt.alfast .and. almaxj.lt.alfast )
     $           emax = emaxsl
                rijmx = emax*(alaminj+almini) / (alaminj*almini)
                if( rijsq.gt.rijmx ) goto 799
c
c               Angular loops:
c
                ij0 = ijstart + jorbat0*norbati + iorbat0
                ij0 = ij0 - norbati
                mli1 = 1
                do 420 mlj=1,norbshj
c                 lj=1,....,9   is orbital symmetry index of orbital j
                  lj = lj0 + mlj
                  jorbat = jorbat0 + mlj
                  j = j0 + mlj
c
                  ij0 = ij0 + norbati
                  if( diagshl ) mli1 = mlj
                  do 410 mli=mli1,norbshi
c                   li=1,....,9   is orbital symmetry index of orbital i
                    li = li0 + mli
                    iorbat = iorbat0 + mli
                    i = i0 + mli
                    ij = ij0 + mli
c
c                   Alpha loops:
c
                    nx = 0
                    do  n=1,3
                      sgr(n) = zero
                      tgr(n) = zero
                    enddo
c
                    do 320 jal=1,nalj
                      alj = ala(jal,jshl,jtyp)
c
                      do 310 ial=1,nali
                        ali = ala(ial,ishl,ityp)
c
                        rijmx = emax*(ali+alj) / (ali*alj)
                        if( rijsq.gt.rijmx ) goto 310
c
                        nx = nx + 1
                        alsum = ali + alj
                        do  ir=1,3
                          p(ir) = two*alj*rij(ir)
                        enddo
                        pmag = SQRT( p(1)*p(1)+p(2)*p(2)+p(3)*p(3) )
c
c                       Evaluate contribution to stress and forces
c
                        call FINTSUM( engrt,engrs,p,lj,li,
     $                   cj(1,lj),ci(1,li), alj,ali,alsum, pmag )
c
                        ffac = EXP( -ali*alj*rijsq/alsum )*
     $                   cala(ial,ishl,ityp)*cala(jal,jshl,jtyp)
                        do  n=1,3
                          sgr(n) = sgr(n) + engrs(n)*ffac
                          tgr(n) = tgr(n) + engrt(n)*ffac
                        enddo
c
c                       Close alpha loops:
  310                 continue
  320               continue
c
c                   Multiply by Bloch phase factors and add to stress and forces
c
                    emfac = zero
                    dmfac = zero
                    if( i.eq.j )then
                      do  k=1,nk
C                        emfac = emfac + coskr(k,ilat)*emat(i,i,k)
C                        dmfac = dmfac + coskr(k,ilat)*dmat(i,i,k)
                        emfac = emfac + coskr(k,ilat)*emat(k,ij)
                        dmfac = dmfac + coskr(k,ilat)*dmat(k,ij)
                      enddo
                    else
c                     i.ne.j
                      do  k=1,nk
C                        emfac = emfac + coskr(k,ilat)*emat(i,j,k)
C                        dmfac = dmfac + coskr(k,ilat)*dmat(i,j,k)
                        emfac = emfac + coskr(k,ilat)*emat(k,ij)
                        dmfac = dmfac + coskr(k,ilat)*dmat(k,ij)
                      enddo
                      if( ncplx.eq.2 )then
                        ijc = ij + matij
                        do  k=1,nk
C                          emfac = emfac + sinkr(k,ilat)*emat(j,i,k)
C                          dmfac = dmfac + sinkr(k,ilat)*dmat(j,i,k)
                          emfac = emfac + sinkr(k,ilat)*emat(k,ijc)
                          dmfac = dmfac + sinkr(k,ilat)*dmat(k,ijc)
                        enddo
                      endif
c
c                     Double because only doing one triangle
                      emfac = emfac + emfac
                      dmfac = dmfac + dmfac
                    endif
c
                    facs = emfac*fourpi / ( vnorm(li)*vnorm(lj) )
                    fact = dmfac*fourpi / ( vnorm(li)*vnorm(lj) )
c
                    fxs = facs*sgr(1)
                    fys = facs*sgr(2)
                    fzs = facs*sgr(3)
                    fxt = fact*tgr(1)
                    fyt = fact*tgr(2)
                    fzt = fact*tgr(3)
c
                    frcs(1,iatm) = frcs(1,iatm) - fxs
                    frcs(2,iatm) = frcs(2,iatm) - fys
                    frcs(3,iatm) = frcs(3,iatm) - fzs
                    frcs(1,jatm) = frcs(1,jatm) + fxs
                    frcs(2,jatm) = frcs(2,jatm) + fys
                    frcs(3,jatm) = frcs(3,jatm) + fzs
                    strs(1,1) = strs(1,1) - fxs*rij(1)
                    strs(2,1) = strs(2,1) - fys*rij(1)
                    strs(3,1) = strs(3,1) - fzs*rij(1)
                    strs(1,2) = strs(1,2) - fxs*rij(2)
                    strs(2,2) = strs(2,2) - fys*rij(2)
                    strs(3,2) = strs(3,2) - fzs*rij(2)
                    strs(1,3) = strs(1,3) - fxs*rij(3)
                    strs(2,3) = strs(2,3) - fys*rij(3)
                    strs(3,3) = strs(3,3) - fzs*rij(3)
                    frct(1,iatm) = frct(1,iatm) + fxt
                    frct(2,iatm) = frct(2,iatm) + fyt
                    frct(3,iatm) = frct(3,iatm) + fzt
                    frct(1,jatm) = frct(1,jatm) - fxt
                    frct(2,jatm) = frct(2,jatm) - fyt
                    frct(3,jatm) = frct(3,jatm) - fzt
                    strt(1,1) = strt(1,1) + fxt*rij(1)
                    strt(2,1) = strt(2,1) + fyt*rij(1)
                    strt(3,1) = strt(3,1) + fzt*rij(1)
                    strt(1,2) = strt(1,2) + fxt*rij(2)
                    strt(2,2) = strt(2,2) + fyt*rij(2)
                    strt(3,2) = strt(3,2) + fzt*rij(2)
                    strt(1,3) = strt(1,3) + fxt*rij(3)
                    strt(2,3) = strt(2,3) + fyt*rij(3)
                    strt(3,3) = strt(3,3) + fzt*rij(3)
c
c                   Close orbital loops:
  410             continue
  420           continue
c
  799           continue
c               Close Atom i shell loop:
                i = i0 + norbshi
  800         continue
c
c             Close Atom i lattice loop:
  900       continue
c
c           Close Atom i loop:
 1000     continue
c
 1799     continue
c         Close Atom j shell loop:
          j = j0 + norbshj
 1800   continue
c
c       Close Atom j loop:
 2000 continue
c
c >>>> Collect final force/stress
c
      if( nprocs .gt. 1 )then
c       Each local processor has forces over its share of work.
c       Reduce(sum) these to master (node 0).
c       Separate reduces or combined reduce - no diff in performance
c       For now, also broadcast total to all processors.
c
        nreduc = 3*natm
c
        call MPREDUC8( master, nreduc, frcs, frctp , icomm )
        call DCOPY( 3*natm, frctp,1, frcs,1 )
        call MPBCAST8( master, nreduc, frcs, icomm )
c
        call MPREDUC8( master, nreduc, frct, frctp , icomm )
        call DCOPY( 3*natm, frctp,1, frct,1 )
        call MPBCAST8( master, nreduc, frct, icomm )
c       
        nreduc = 9
c       
        call MPREDUC8( master, nreduc, strs, strtp, icomm )
        call DCOPY( 9, strtp,1, strs,1 )
        call MPBCAST8( master, nreduc, strs, icomm )
c       
        call MPREDUC8( master, nreduc, strt, strtp, icomm )
        call DCOPY( 9, strtp,1, strt,1 )
        call MPBCAST8( master, nreduc, strt, icomm )
c       
c       End tp
      endif
c
c    That's all Folks!
c
      RETURN
      END
