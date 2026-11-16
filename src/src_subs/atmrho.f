c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ATMRHO
c
c
      subroutine ATMRHO( ilocrho,
     $ natm,ntyp, nshld,nald, noad, nrd, nlat, nang,
     $ itypa, numshl,lshel,nala,ala,cala, occ,nfityp, nrad,radmsh,
     $ ratm,rlat, angpts,
     $ rhorad,rsph,rpt,rsq,rsqtol )
c---------------------------------------------------------------
c Purpose: to compute local radial charge about atoms due to all
c          OTHER atoms for the purpose of getting 0th order
c          correction for the local density used to generate the
c          local xc-potential/energy density for 1-center fit.
c
c Written: Peter A. Schultz, 19-February-1993
c
c Revision history:
c  11Mar02-PAS/2.52: literal reals to d0
c   1Mar00-PAS/2.43: bugfix for non-xcfit
c  21Dec99-PAS/2,42: cosmetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Apr97-PAS/2.17: naming conventions
c---------------------------------------------------------------
c
c  Note: could use orbital occupations to give directionality
c        to the density.  I am guessing using spherical atoms
c        with shell occupations is good enough. 19Feb93-PAS
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c data arrays:
      DIMENSION  itypa(natm), numshl(ntyp),lshel(nshld,ntyp),
     $  nala(nshld,ntyp),ala(nald,nshld,ntyp),cala(nald,nshld,ntyp)
      DIMENSION  occ(noad,*), nfityp(*),nrad(*),radmsh(nrd,*)
      DIMENSION  ratm(3,natm),rlat(3,nlat)
      DIMENSION  angpts(3,*)
c scratch arrays:
      DIMENSION  rhorad(*)
      DIMENSION  rsph(3,*),rpt(3,*),rsq(*),rsqtol(*)
c local declarations:
      DIMENSION  rnuc(3),ratom(3)
      DATA  zero,half / 0.d0,0.5d0 /
c
c >>>> EXECUTABLE CODE:
c
      fourpi = pi + pi + pi + pi
      argmx = half*cutexp
c
c      Loop over sites, include only those with potentials:
c
      do 2000 nuc=1,natm
        ktyp = itypa(nuc)
        if( nfityp(ktyp).eq.0 ) goto 2000
c
        do  j=1,3
          rnuc(j) = ratm(j,nuc)
        enddo
        nrxc = nrad(ktyp)
c  Blank density array, compute positions of site's cluster of points
        ir = 0
        do  iang=1,nang
          do  nr=1,nrxc
            r = radmsh(nr,ktyp)
            ir = ir + 1
            rhorad(ir) = zero
            rsph(1,ir) = angpts(1,iang)*r
            rsph(2,ir) = angpts(2,iang)*r
            rsph(3,ir) = angpts(3,iang)*r
          enddo
        enddo
        nrsph = ir
c
c      Loop over atoms, find those with occupied shells
c
        do 1000 iatm=1,natm
          ityp = itypa(iatm)
c         Make sure the site's own density is excluded
          ilat1 = 1
          if( iatm.eq.nuc ) ilat1 = 2
c
          iorbat = 0
          do 900 ishl=1,numshl(ityp)
            lshl = lshel(ishl,ityp)
c
            occshl = zero
            do  ll=1,2*lshl+1
              iorbat = iorbat + 1
              occshl = occshl + occ(iorbat,ityp)
            enddo
c           Check if shell is occupied
            if( ABS( occshl ).lt.1.d-10 ) goto 900
c
            nal = nala(ishl,ityp)
            almini = ala(1,ishl,ityp)
            do  ial=1,nal
              if(ala(ial,ishl,ityp).lt.almini)almini=ala(ial,ishl,ityp)
            enddo
            alsmin = almini + almini
c
c      Loop over this atom in all unit cells:
c
            do 500 ilat=ilat1,nlat
c
c             Get nuclear position wrt this atom
              do  j=1,3
                ratom(j) =  rnuc(j) - ratm(j,iatm) - rlat(j,ilat)
              enddo
              rsqat = ratom(1)**2 + ratom(2)**2 + ratom(3)**2
c
c             Check if atom shell density overlaps site
              if( alsmin*rsqat.gt.argmx ) goto 500
c
c             Compute positions of site cluster of points wrt atom
              do  ir=1,nrsph
                rpt(1,ir) = ratom(1) + rsph(1,ir)
                rpt(2,ir) = ratom(2) + rsph(2,ir)
                rpt(3,ir) = ratom(3) + rsph(3,ir)
                rsq(ir) = rpt(1,ir)**2 + rpt(2,ir)**2 + rpt(3,ir)**2
              enddo
c
              if( lshl.gt.0 )then
                do  ir=1,nrsph
                  rsqtol(ir) = rsq(ir)**lshl
                enddo
              endif
c
              do 420 ial=1,nal
                ali = ala(ial,ishl,ityp)
                cali = occshl*cala(ial,ishl,ityp)
                cali0 = cali
c
                do 410 jal=ial,1,-1
                  alsum = ali + ala(jal,ishl,ityp)
                  if( alsum*rsqat.gt.argmx ) goto 410
c
                  prod = cali*cala(jal,ishl,ityp)
                  cali = cali0 + cali0
                  rhoprod = prod/fourpi
c
                  ir0 = 0
                  do 300 iang=1,nang
                    rendsq = rsq(ir0+nrxc)
                    nrxcn = nrxc
c
                    if( alsum*rendsq.gt.argmx )then
c                     Outer end of points too large, cut em off
                      rsqcut = argmx/alsum
                      do 200 ir=ir0+nrxc,ir0+1,-1
                        if( rsq(ir).lt.rsqcut ) goto 201
                        nrxcn = nrxcn - 1
  200                 continue
  201                 continue
                    endif
c
c      Compute contribution of density into this angular slice of points
c
                    if( lshl.eq.0 )then
c
                      do  ir=ir0+1,ir0+nrxcn
                        rhorad(ir) = rhorad(ir) + rhoprod*
     $                                EXP(-alsum*rsq(ir))
                      enddo
c
                    else
c
                      do  ir=ir0+1,ir0+nrxcn
                        rhorad(ir) = rhorad(ir) + rhoprod*rsqtol(ir)*
     $                                EXP(-alsum*rsq(ir))
                      enddo
c
                    endif
c
                    ir0 = ir0 + nrxc
  300             continue
c
c    Close alpha loops:
c
  410           continue
  420         continue
c
c    Close lattice loop:
c
  500       continue
c
c    Close shell loop:
c
  900     continue
c
c    Close atom loop:
c
 1000   continue
c
c    Write local density about this site to disk:
c
      write(ilocrho) nrxc,nang
      write(ilocrho) (rhorad(ir),ir=1,nrsph)
c
c    Close nuclear loop:
c
 2000 continue
c
c    That's all Folks!
c
      RETURN
      END
