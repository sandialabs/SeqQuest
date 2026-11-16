c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SYMAGES
c
c
      subroutine SYMAGES( IWR, subsym, rsym,rprim, natm,ratm,itypa,
     $ natmnm,atmnm,
     $ nsym,nsyma, naofsym, rmatsym )
c---------------------------------------------------------------
c Purpose: find symmetry images of atoms
c
c Written: Peter A. Schultz, 26-January-1996, for v2.15
c
c Revision history:
c  23Jun15-PAS/2.67: atom labels, and some extra cosmetics
c   1May02-PAS/2.53: reconcile sym and cell-opt; cosmetics
c  21Jun01-PAS/2.48: replace STOPs
c   6Jan99-PAS/2.30: redo symmetry scheme
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Sep98-PAS/2.27: control level for output
c  24Sep97-PAS/2.21: move commentary on atom mapping here
c  24Apr97-PAS/2.17: naming conventions
c---------------------------------------------------------------
c
c  Use tranformation matrices to produce images of each atom,
c  and map onto original atoms
c
c  Variables:
c   subsym    - true  = ok to reduce symmetry
c               false = trigger error upon symmetry failure
c   nsym      - # input non-identity symmetry operations
c   nsyma     - # non-identity symmetries obeyed by atoms
c   rsym()    - offset from origin of center of symmetry
c   naofsym() - mapping of atoms to other atoms by symmetry
c   rmatsym() - cartesian representations of symmetry operators
c
c  Note: that naofsym begins with the identity, while the
c  cartesian symmetry operations skip the identity, so that
c  the indexing for naofsym and rmatsym differ by 1.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c
c input arrays:
      LOGICAL    subsym
      DIMENSION  rsym(3),rprim(3,3), ratm(3,*), itypa(*)
      CHARACTER*(*)  atmnm
      DIMENSION  atmnm(*)
      DIMENSION  rmatsym(3,4,*)
c output array:
      DIMENSION  naofsym(natm,*)
c local declarations:
      DIMENSION  rsymctr(3),rsyveci(3)
      DIMENSION  ri(3),sri(3),tri(3),vij(3),clat(3)
      LOGICAL  outuc,shift, badgrp,goodsym
      DATA  rngtol / 0.0001d0 /
      DATA  zero,one / 0.d0,1.d0 /
      DATA  lstout / 0 /
c
c >>>> EXECUTABLE CODE
c
      badgrp = .false.
c  All atoms map to self under identity:
      do  iatm=1,natm
        naofsym(iatm,1) = iatm
      enddo
      if( nsym.le.0 ) RETURN
c
c Symmetry center offset:
      rsymctr(1) = rsym(1)
      rsymctr(2) = rsym(2)
      rsymctr(3) = rsym(3)
c
      nsyma = nsym
      isyma = 0
      do 300 isym=1,nsym
c
        goodsym = .true.
        do  iatm=1,natm
          naofsym(iatm,isym+1) = -1
        enddo
c
c       Convert symmetry element translation to Cartesian coords:
        call LVXPND( rprim, rsyveci, rmatsym(1,4,isym) )
        rsymag = rsyveci(1)**2 + rsyveci(2)**2 + rsyveci(3)**2
c
        do 200 iatm=1,natm
          ityp = itypa(iatm)
c
c         Do transformations about center of symmetry:
          ri(1) = ratm(1,iatm) - rsymctr(1)
          ri(2) = ratm(2,iatm) - rsymctr(2)
          ri(3) = ratm(3,iatm) - rsymctr(3)
c
c         Transform coords of this atom by this symmetry operation
          sri(1) = rmatsym(1,1,isym)*ri(1) + rmatsym(2,1,isym)*ri(2)
     $           + rmatsym(3,1,isym)*ri(3) + rsyveci(1)
     $           + rsymctr(1)
          sri(2) = rmatsym(1,2,isym)*ri(1) + rmatsym(2,2,isym)*ri(2)
     $           + rmatsym(3,2,isym)*ri(3) + rsyveci(2)
     $           + rsymctr(2)
          sri(3) = rmatsym(1,3,isym)*ri(1) + rmatsym(2,3,isym)*ri(2)
     $           + rmatsym(3,3,isym)*ri(3) + rsyveci(3)
     $           + rsymctr(3)
c
c         Project image atom coordinates into primary unit cell ...
          call PROJUC( outuc,shift, 3, sri,tri, rprim )
c
c          ... and match against atom list, checking first against
c         itself (most common), and then against all atoms.
c         If this fails, try local search by lattice vector as we
c         are guaranteed to be close.
c
          itry = 1
   90     jatm = iatm
          jnext = 1
  100     jtyp = itypa(jatm)
c
c         Check against atom type first
          if( jtyp.ne.ityp) goto 140
c
          vij(1) = ratm(1,jatm) - tri(1)
          vij(2) = ratm(2,jatm) - tri(2)
          vij(3) = ratm(3,jatm) - tri(3)
          if( itry.eq.1 )then
            rij = SQRT( vij(1)**2 + vij(2)**2 + vij(3)**2 )
          else
c           Try a local more expensive search on second go-round
            call LVPROJ( rprim, vij, clat )
            rij = zero
            do  idim=1,3
              cl0 = ABS( clat(idim) )
              cl1 = ABS( cl0-one )
              rij = rij + MIN( cl0, cl1 )
            enddo
          endif
c
          if( rij.lt.rngtol )then
c           Have jatm = S(iatm), or iatm=(S**-1)(jatm)
            naofsym(jatm,isyma+2) = iatm
            goto 190
          endif
c
  140     jatm = jnext
          jnext = jnext + 1
          if( jatm.le.natm ) goto 100
c
c         Failed to find symmetry image, try alternate search ...
          itry = itry + 1
          if( itry.eq.2 ) goto 90
c
c          ... and alternate search failed, too:
          if( .not. subsym .or. lstout.gt.2 )then
c           Have hard check, trigger error.
 9131       format(1x,a,3f15.8)
            if( goodsym )then
c             First failure found in this symmetry, describe operator:
              write(IWR,*) '***** NOTE: Atomic symmetry failure'
              write(IWR,*)'Cartesian representation of failed symmetry:'
              write(IWR,'(5x,3f12.8)') ((rmatsym(i,j,isym),i=1,3),j=1,3)
              if( rsymag .gt. 0.00001d0 )then
                write(IWR,*) ' ... and translation vector:'
                write(IWR,'(5x,3f12.8)') rsyveci
              endif
              write(IWR,9131) 'Symmetry center offset = ', rsymctr
              if(  .not. subsym .or. lstout.gt.4 )
     $        write(IWR,*) 'Coords (internal - scaled and shifted):'
            endif
c
            write(IWR,9130)  iatm,  atmnm(iatm)
 9130       format(1x,'>>>>> Atom#',i6,' = ',a,' fails to have image')
            if( .not. subsym .or. lstout.gt.4 )then
              write(IWR,9131) '    Internal coords   =',
     $                        ratm(1,iatm),ratm(2,iatm),ratm(3,iatm)
              write(IWR,9131) '    Symmetry center   =', rsymctr
              write(IWR,9131) '    Coordinates       =', ri
              write(IWR,9131) '    Sym transformation=',sri
              write(IWR,9131) '    Map into unit cell=',tri
              write(IWR,*)   '     ... for which there is no match'
            endif
          endif
c         Mark this symmetry as failed, so do not repeat operator
          goodsym = .false.
c
  190     continue
c
c         Close loop over atoms:
c
  200   continue
c
        if( goodsym )then
          isyma = isyma + 1
          if( isyma.ne.isym )then
c           Reorder matrix ops to put active ops to front:
            do  j=1,4
              do  i=1,3
                tmp = rmatsym(i,j,isyma)
                rmatsym(i,j,isyma) = rmatsym(i,j,isym)
                rmatsym(i,j,isym)  = tmp
              enddo
            enddo
          endif
        else
          if( .not. subsym )then
            isyma = isyma + 1
            badgrp = .true.
          endif
          nsyma = nsyma - 1
        endif
c
c       Close symmetry operation loop:
c
  300 continue
c
      if( nsyma.ne.nsym )then
        write(IWR,*) '***** NOTE: atomic symmetry is reduced'
        write(IWR,*) '  # of input symmetry operations (w/id)=',nsym+1
        write(IWR,*) '  # of valid atomic symmetries (w/id)  =',nsyma+1
      endif
c
c  Use isyma rather than nsyma as it puts out full mapping, even
c  with an error condition
      if( isyma.gt.0 .and. ( lstout.gt.4 .or. badgrp ) )then
        write(IWR,'(1h0,a)')'atom, atom mapping under inverse sym ops'
        do  iatm=1,natm
          write(IWR,'(1x,24(i4))')  (naofsym(iatm,isym),isym=1,isyma)
        enddo
      endif
c
      if( badgrp ) call STOPXERR( 'bad-sym - SYMAGES atom map fails' )
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
c End
