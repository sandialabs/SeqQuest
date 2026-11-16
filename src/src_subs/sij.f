c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SIJ
c
c
      subroutine SIJ( norb,nk,ncplx, natm,ntyp,nshld,nald, nlat,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba,
     $ ratm, rlat,coskr,sinkr, wksm,
     $ smat )
c---------------------------------------------------------------
c Purpose: compute overlap matrix for a band structure
c          calculation. The real part is stored in the lower
c          triangle of matrix, the imaginary part in the upper
c
c Written: P.A. Schultz, based on P.J. Feibelman's original
c
c Revision history:
c  17Oct08-PAS/2.62: stripe-distributed-parallel
c   7Mar02-PAS/2.52: ang constants extracted
c  14May01-PAS/2.47: extract s1atom
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  13Aug98-PAS/2.24: cutoffs, names, restructure
c  16Jul98-PAS/2.23: fix label problem (thx D.Raczkowski)
c  28Apr97-PAS/2.17: naming conventions
c  12Apr95-PAS/2.12: stretched alpha cutoffs installed
c  18Jan94-PAS/2.04: full cleanup, match sij,tij,enij,ennloc
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
      PARAMETER  ( zero=0.d0, two=2.d0, half=0.5d0 )
c
c Output array:
      DIMENSION  smat(nk,*)
c Input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat),coskr(nk,nlat),sinkr(nk,nlat)
c Scratch space:
      DIMENSION  wksm(*)
c
c Local declarations:
      LOGICAL    diagshl
      DIMENSION  rij(3), p(3), ci(10,9),cj(10,9)
c
      DIMENSION  vnorm(9)
c
c >>>> EXECUTABLE CODE:
c
C      write(6,*) 'SDEV/SIJ: stripe-parallel'
c Get angular normalizations:
      call VNORMS9( vnorm )
c
      argmx = half*cutgaus
      emaxf = cutgwfn
      emaxsl = cutslow
      fourpi = pi + pi + pi + pi
c
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
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
c  Clear overlap matrix
c
      call MKZERO( nmatpr, smat )
c
c     Atom j loop:
c
      nlatmx = 0
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
                    ssum = zero
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
c                       Evaluate contribution to matrix element
                        call INTSUM( selem, p,lj,li,cj(1,lj),ci(1,li),
     $                               alsum )
                        ssum = ssum + EXP( -ali*alj*rijsq/alsum )*
     $                   selem*cala(ial,ishl,ityp)*cala(jal,jshl,jtyp)
c
c                       Close alpha loops:
  310                 continue
  320               continue
                    ssum = fourpi*ssum / ( vnorm(li)*vnorm(lj) )
c
c                   Multiply by Bloch phase factors and load output array
c
                    do  k=1,nk
C                      smat(i,j,k) = smat(i,j,k) + ssum*coskr(k,ilat)
                      smat(k,ij) = smat(k,ij) + ssum*coskr(k,ilat)
                    enddo
                    if( ncplx.eq.2 .and. i.ne.j )then
                      ijc = ij + matij
                      do  k=1,nk
C                        smat(j,i,k) = smat(j,i,k) + ssum*sinkr(k,ilat)
                        smat(k,ijc) = smat(k,ijc) + ssum*sinkr(k,ilat)
                      enddo
                    endif
c
                    if( ilat.gt.nlatmx ) nlatmx = ilat
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
c         Close Atom j shell loop:
          j = j0 + norbshj
 1800   continue
c
c       Close Atom j loop:
 2000 continue
c
      if( nprocs .gt. 1 )then
        call TIMER('  sij compute    ')
        if( iproc.eq.master )then
          do  jproc=1,nprocstr-1
            lenmsg = 1
            call MPRECVI( jproc, lenmsg, nlatpr, icomm )
            if( nlatpr.gt.nlatmx ) nlatmx = nlatpr
          enddo
        elseif( iproc.lt.nprocstr )then
          lenmsg = 1
          call MPSENDI( master, lenmsg, nlatmx, icomm )
        endif
c
c       Broadcast the result to everyone:
c       (do not use openmpi Bcast cuz of interconnect init cost)
        if( iproc.eq.master )then
          do  jproc=1,nprocs-1
            lenmsg = 1
            call MPSENDI( jproc, lenmsg, nlatmx, icomm )
          enddo
        else
          lenmsg = 1
          call MPRECVI( master, lenmsg, nlatmx, icomm )
        endif
        call TIMER('  sij/nlat-bcast ')
        nlat = nlatmx
      endif
c
c  Shrink nlat if list end not used (ie, unit cells not interact):
      nlat = nlatmx
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SIJ1
c
c
      subroutine SIJ1( norb,nk,ncplx, natm,ntyp,nshld,nald, nlat,
     $ itypa, numshl,lshel,nala,ala,cala, alamin, norba,
     $ ratm, rlat,coskr,sinkr,
     $ smat )
c---------------------------------------------------------------
c Purpose: compute overlap matrix for a band structure
c          calculation. The real part is stored in the lower
c          triangle of matrix, the imaginary part in the upper
c
c Written: P.A. Schultz, based on P.J. Feibelman's original
c
c Revision history:
c   6Sep08-PAS/2.62: switched ij-ji
c   7Mar02-PAS/2.52: ang constants extracted
c  14May01-PAS/2.47: extract s1atom
c  21Dec99-PAS/2.42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  13Aug98-PAS/2.24: cutoffs, names, restructure
c  16Jul98-PAS/2.23: fix label problem (thx D.Raczkowski)
c  28Apr97-PAS/2.17: naming conventions
c  12Apr95-PAS/2.12: stretched alpha cutoffs installed
c  18Jan94-PAS/2.04: full cleanup, match sij,tij,enij,ennloc
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
      PARAMETER  ( zero=0.d0, two=2.d0, half=0.5d0 )
c
c Output array:
      DIMENSION  smat(norb,norb,nk)
c Input arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  alamin(ntyp), norba(ntyp)
      DIMENSION  ratm(3,natm)
      DIMENSION  rlat(3,nlat),coskr(nk,nlat),sinkr(nk,nlat)
c Local declarations:
      LOGICAL    diagshl
      DIMENSION  rij(3), p(3), ci(10,9),cj(10,9)
c
      DIMENSION  vnorm(9)
c
c >>>> EXECUTABLE CODE:
c
C      write(6,*) 'SDEV/SIJ1: ijkmat'
c Get angular normalizations:
      call VNORMS9( vnorm )
c
      argmx = half*cutgaus
      emaxf = cutgwfn
      emaxsl = cutslow
      fourpi = pi + pi + pi + pi
      nlatmx = 0
c
c  Clear overlap matrix
c
      call MKZERO( norb*norb*nk, smat )
c
c     Atom j loop:
c
      jnext = 0
      do 2000 jatm=1,natm
        jtyp = itypa(jatm)
        nshlj = numshl(jtyp)
        norbatj = norba(jtyp)
        jstart = jnext
        jnext = jstart + norbatj
c
        alaminj = alamin(jtyp)
        call CINIT( cj )
c
c       Atom j shell loop:
c
        j = jstart
        do 1800 jshl=1,nshlj
          nalj = nala(jshl,jtyp)
          lshlj = lshel(jshl,jtyp)
          norbshj = 2*lshlj + 1
          lj0 = lshlj*lshlj
c
          j0 = j
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
          inext = jstart
          do 1000 iatm=jatm,natm
            ityp = itypa(iatm)
            nshli = numshl(ityp)
            norbati = norba(ityp)
            istart = inext
            inext = istart + norbati
c
c           Following helps ensure only i.ge.j done:
            ishl1 = 1
            if( iatm.eq.jatm )then
              ishl1 = jshl
              do  ishl=1,ishl1-1
                lshli = lshel(ishl,ityp)
                norbshi = 2*lshli + 1
                istart = istart + norbshi
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
              i = istart
              do 800 ishl=ishl1,nshli
                nali = nala(ishl,ityp)
                lshli = lshel(ishl,ityp)
                norbshi = 2*lshli + 1
                li0 = lshli*lshli
c
                i0 = i
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
                mli1 = 1
                do 420 mlj=1,norbshj
c                 lj=1,....,9   is orbital symmetry index of orbital j
                  lj = lj0 + mlj
                  j = j0 + mlj
c
                  if( diagshl ) mli1 = mlj
                  do 410 mli=mli1,norbshi
c                   li=1,....,9   is orbital symmetry index of orbital i
                    li = li0 + mli
                    i = i0 + mli
c
c                   Alpha loops:
c
                    nx = 0
                    ssum = zero
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
c                       Evaluate contribution to matrix element
                        call INTSUM( selem, p,lj,li,cj(1,lj),ci(1,li),
     $                               alsum )
                        ssum = ssum + EXP( -ali*alj*rijsq/alsum )*
     $                   selem*cala(ial,ishl,ityp)*cala(jal,jshl,jtyp)
c
c                       Close alpha loops:
  310                 continue
  320               continue
                    ssum = fourpi*ssum / ( vnorm(li)*vnorm(lj) )
c
c                   Multiply by Bloch phase factors and load output array
c
                    do  k=1,nk
                      smat(i,j,k) = smat(i,j,k) + ssum*coskr(k,ilat)
                    enddo
                    if( ncplx.eq.2 .and. i.ne.j )then
                      do  k=1,nk
                        smat(j,i,k) = smat(j,i,k) + ssum*sinkr(k,ilat)
                      enddo
                    endif
c
                    if( ilat.gt.nlatmx ) nlatmx = ilat
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
c         Close Atom j shell loop:
          j = j0 + norbshj
 1800   continue
c
c       Close Atom j loop:
 2000 continue
c
c  Shrink nlat if list end not used (ie, unit cells not interact):
      nlat = nlatmx
c
c    That's all Folks!
c
      RETURN
      END
