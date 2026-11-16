c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GUPDATE
c
c
      subroutine GUPDATE( ierrg, ratm0,
     $ natm, iforce, ndim,nsym,nsyma, rprim, rmatsym, naofsym,
     $ nslics,islics,islicats,vslics,    iatframe,
     $ ratmP,
     $ ratmg, ratmx, dratm0, dratmsy )
c---------------------------------------------------------------
c Purpose: take new geometry, check symmetry, and update the
c          master coordinate vector.
c
c Written: Peter A. Schultz, 18-May-2001, for v2.47
c
c Revision history:
c   1Aug17-PAS/2.68: bonds/slab constraints; unified SLICS
c   5Aug15-PAS/2.67: return error flag rather than stop upon failure
c  20Feb09-PAS/2.62: add options for vgfix atoms along vector
c  14Jun05-PAS/2.59: File units
c  11Feb02-PAS/2.54: install frame-constraint, better documentation
c  30Oct01-PAS/2.50: documentation typo
c  21Jun01-PAS/2.48: replace STOPs
c  18May01-PAS/2.47: extract from main routine
c---------------------------------------------------------------
c
cRRR add in stuff for coordinate de-projection
cRRR needed:  check max allowed angular motions
cRRR                project internal coordinates back to cartesian
cRRR                update cartesian positions in master coord array
c
c Notes:
c  The packed (by relax v. fixed) new coordinates are in ratmg().
c  They are unpacked here, and the difference from the original
c  coordinates computed.  That difference is checked for symmetry,
c  and then used to update the master coordinate array in ratm0().
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero= 0.d0,half=0.5d0 )
c Input/Output array:
      DIMENSION  ratm0(3,natm)
c Input arrays:
      DIMENSION  iforce(natm)
      DIMENSION  rprim(3,3), rmatsym(*), naofsym(*)
      DIMENSION  ratmg(*)
c
      DIMENSION  islics(4,*),islicats(*),vslics(3,*)
      DIMENSION  iatframe(3)
cRMV for slic
      DIMENSION  ratmP(3,natm)
c
c Scratch space:
      DIMENSION  ratmx(3,natm),dratm0(3,natm),dratmsy(3,natm)
c
c Local declarations:
      DIMENSION  dum9(3,3)
      DATA       dum9 / 9*0.d0 /
c Diagnostic level
      DATA  iprint / 2 /
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
      ierrg = 0
c
      if( iprint.gt.0 )
     $ write(IWR,*) '>>>>> GUPDATE begins, nslics=', nslics
c
c Expand relaxed coordinates into full vector:
c
      irlx = 0
      do  iatm=1,natm
        if( iforce(iatm) .eq. 0 )then
c         Frozen atoms: Copy in old coordinates for unmodified atoms
          ratmx(1,iatm) = ratm0(1,iatm)
          ratmx(2,iatm) = ratm0(2,iatm)
          ratmx(3,iatm) = ratm0(3,iatm)
        elseif( iforce(iatm) .ne. 0 )then
c         Install new coordinates for relaxed atoms
          do  idf=1,3
            irlx = irlx + 1
            ratmx(idf,iatm) = ratmg(irlx)
          enddo
        endif
      enddo
c
c Compute change vector ...
c
      do 1000 islic=1,nslics
        islictyp = islics(1,islic)
        if( iprint.gt.1 ) 
     $   write(IWR,*) '  islic,islictyp=',islic,islictyp
c
        if(     islictyp.eq.1 .or. islictyp.eq.2 )then
c         SLIC projections - 1-fold sphere, 2-fold cylinder
          iatmproj = islics(2,islic)
          if( iprint.gt.1 ) write(IWR,*) 'SLIC/projection' 
          if( iprint.gt.2 )
     $    write(IWR,*) 'Active update slic-projected atom: ',iatmproj
          if( iprint.gt.3 )then
            write(IWR,*) 'old proj coords ',(ratmP(i,iatmproj),i=1,3)
            write(IWR,*) 'new proj coords ',(ratmx(i,iatmproj),i=1,3)
          endif
c****** DEV Should reverse cyl/sph order to 1-fold(sph) and 2-fold (cyl)
          if( islictyp.eq.2 )then
            write(IWR,*) 'SLIC/slsphb: de-project 1-fold spherical'
            call SLSPHB( IWR,natm, islics,islic,
     $                   ratm0,ratmx,ratmP )
          elseif( islictyp.eq.1 )then
            write(IWR,*) 'SLIC/slcylb: de-project 2-fold cylinder'
            call SLCYLB( IWR,natm, islics,islic,
     $                   ratm0,ratmx,ratmP )
          endif
c
        elseif( islictyp.eq.4 )then
c         BONDS - do "trivial" thing first, simply scale along bond
          if( iprint.gt.1 )  write(IWR,*) 'SLICS/BOND-length'
          call SLBONDFR( IWR, islic,islics,
     $     ratm0,ratmx,ratmP )
c
        elseif( islictyp.ge.5 .and. islictyp.le.8 )then
c       . VGFIX: purify the vector-projected atom geometry update:
          if( iprint.gt.1 ) 
     $     write(IWR,*) 'SLICS/VGFIX: slictyp=',islictyp
          nvgfix = islics(2,islic)
          if( nvgfix .gt. 0 )then
            call FLGETIWR( IWR )
cpas:       (vgfix should ultimately share SLIC data structures)
            isl1 = islics(3,islic)
            iatmvg1 = islicats(isl1)
            islicv = islics(4,islic)
c           Get center of mass motion along vector from lead vgfix atom
            dvec = zero
            do  id=1,3
C              dvec = dvec + vslics(id,islicv)*dratm0(id,iatmvg1)
              dvec = dvec + vslics(id,islicv)*
     $                     (ratmx(id,iatmvg1)-ratm0(id,iatmvg1))
            enddo
            if( iprint.gt.3 )
     $      write(IWR,*) 'VGFIX: purify coordinate update, dvec=',dvec
c
            do  ivgatm=isl1,isl1+nvgfix-1
              iatm = islicats(ivgatm)
              dveci = zero
              do  id=1,3
C                dveci = dveci + vslics(id,islicv)*dratm0(id,iatm)
                dveci = dveci + vslics(id,islicv)*
     $               (ratmx(id,iatm)-ratm0(id,iatm))
              enddo
c
              if( islictyp .eq. 7 )then
c               Constrain geometry change for atom to be along axis
                if( iprint.gt.3 )
     $          write(IWR,*) 'vector project R-update on atom',iatm
c               Geometry change is strictly motion along vector:
                do  id=1,3
C                  dratm0(id,iatm) = dveci*vslics(id,islicv)
                  ratmx(id,iatm) = ratm0(id,iatm) 
     $                           + dveci*vslics(id,islicv)
                enddo
              else
c               Allow motion normal to axis, adjust motion along axis
                if( islictyp .ne. 5 )then
c                 Install center-of-mass motion along vector
c                 Note that lead atom has cm-vec motion already installed,
c                  and that non-lead atom motion along vec had been suppressed.
                  dveci = dveci - dvec
                endif
                do  id=1,3
C                 dratm0(id,iatm) = dratm0(id,iatm)
                 ratmx(id,iatm) = ratmx(id,iatm)
     $                           - dveci*vslics(id,islicv)
                enddo
              endif
            enddo
c
          endif
c
        elseif( islictyp.eq.9 )then
c         SLINKY?LAYERS slab option
          if( iprint.gt.1 ) write(IWR,*) 'SLICS/SLINKY-LAYERS'
            call SLINKUPD( IWR, islic, islics,islicats,vslics,
     $         ratm0,ratmx,ratmP )
c         End SLICS options
        endif
 1000 continue
c
      do  iatm=1,natm
        dratm0(1,iatm) = ratmx(1,iatm) - ratm0(1,iatm)
        dratm0(2,iatm) = ratmx(2,iatm) - ratm0(2,iatm)
        dratm0(3,iatm) = ratmx(3,iatm) - ratm0(3,iatm)
      enddo
c
      do  iatm=1,natm
        dratmsy(1,iatm) = dratm0(1,iatm)
        dratmsy(2,iatm) = dratm0(2,iatm)
        dratmsy(3,iatm) = dratm0(3,iatm)
      enddo
c
c  ... and symmetrize the change vector
c
c (finished with ratmx and ratmg, and can use as scratch here)
      call SYMFRC( 0, 'move at', dum9, dratmsy, ratmx, ratmg,
     $ ndim, natm, nsyma, rprim, rmatsym, naofsym )
c
c Impose fixed-frame constraints on the change vector ...
c
      if( iprint.gt.1 .and. iatframe(1).ne.0 )
     $  write(IWR,*) 'GUPDATE/FRAME constraint applied'
      call FIXFRAME( natm, iatframe,
     $       ratm0,       dratmsy )
c -->        ratm(3a)-io  fatm(3a)-io
c
c  ... and update coordinates with sym'd change for new geometry:
c
      mxrd = 0
      rdmx = zero
      dratmax = zero
      idratmx = 0
      do  iatm=1,natm
c
c       Update coordinates with symmetrized update vectors
        ratm0(1,iatm) = ratm0(1,iatm) + dratmsy(1,iatm)
        ratm0(2,iatm) = ratm0(2,iatm) + dratmsy(2,iatm)
        ratm0(3,iatm) = ratm0(3,iatm) + dratmsy(3,iatm)
        dratmi = dratmsy(1,iatm)**2 + dratmsy(2,iatm)**2
     $         + dratmsy(3,iatm)**2
        if( dratmi .gt. dratmax )then
          dratmax = dratmi
          idratmx = iatm
        endif
c
c       Check deviation of symmetrized from raw update vectors
        rdiff = ( dratmsy(1,iatm) - dratm0(1,iatm) )**2
     $        + ( dratmsy(2,iatm) - dratm0(2,iatm) )**2
     $        + ( dratmsy(3,iatm) - dratm0(3,iatm) )**2
        if( rdiff.gt.rdmx )then
          rdmx = rdiff
          mxrd = iatm
        endif
c
      enddo
      dratmax = SQRT( dratmax )
      write(IWR,'(1x,a,i6,f12.6)')
     $  'GUPDATE: Largest relax: atom#, dr(au)=',idratmx,dratmax
c
      rdmx =  SQRT( rdmx )
      if( rdmx .gt. 0.0001 )then
        call FLGETIERR( IERRFL )
        write(IERRFL,*) '***** ERROR: Geometry update violates symmetry'
        write(IERRFL,*) '>>>>> Atom#, deviation distance:',mxrd,rdmx
        write(IERRFL,'(a/(1x,i4,3f15.8/5x,3f15.8))')
     $   'Change in coordinates, and symmetrized change vector:',
     $   ( iatm, dratm0(1,iatm), dratm0(2,iatm), dratm0(3,iatm),
     $           dratmsy(1,iatm),dratmsy(2,iatm),dratmsy(3,iatm),
     $     iatm=1,natm )
C        call STOPXERR( 'geom-sym/geometry update violates symmetry' )
        ierrg = 2
      endif
c
  999 continue
c
      if( iprint.gt.0 ) write(IWR,*) '<<<<< GUPDATE ends'
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
