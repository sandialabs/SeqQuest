c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GSELECT
c
c
      subroutine GSELECT( natm, nafrc, iforce, nrlxdim,
     $ nslics,islics,islicats,vslics,
     $ ratmP,frcP,
     $ itypa,atmass, natmg,nxatm,
     $ ratm,frc, ratmg,frcg,amassg )
c---------------------------------------------------------------
c Purpose: select and pack force components of atoms to be
c          relaxed in geometry update
c
c Written: Peter A. Schultz, 18-May-2001, for v2.47
c
c Revision history:
c   1Aug17-PAS/2.68: add bond, unified slics
c  20Feb09-PAS/2.62: refined constraints along vector
c   9Sep08-PAS/2.62: option to freeze atoms along vector
c    Oct06-RMV/2.60: initial SLIC projection
c   3Oct01-PAS/2.50: damped dynamics
c  21Jun01-PAS/2.48: replace STOPs
c  18May01-PAS/2.47: extract from main routine
c  14Jun05-PAS/2.59: File unit
c---------------------------------------------------------------
c
cRMV add in stuff for coordinate de-projection
cRMV needed:  check max 
cRMV                project selected coords cartesian from internal
cRMV                to correct type of internal
cRMV              
c
c Notes:
c  For the moment, we are simply picking frozen and relaxed atoms.
c  This could easily be generalized to enforce internal
c  restrictions on atoms (such as keeping an atom a given plane)
c  by manipulating the forces and coordinates.
c
c On output:
c  natmg   = number of atoms with relaxed coordinates
c  nxatm() = number of degress of freedom per atom
c  nrlxdim = the length of the vector being relaxed
c  ratmg() = packed coordinates of components to be relaxed
c  frcg()  = packed forces of componentes to be relaxed
c  amassg()= packed particle masses for the components to be relaxed
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero=0.d0, one=1.d0 )
c
c Input arrays:
      DIMENSION  iforce(natm)
      DIMENSION  itypa(natm), atmass(*)
      DIMENSION  ratm(3,natm), frc(3,natm)
cRMV constraints/slics
      DIMENSION  islics(4,*), islicats(*), vslics(3,*)
      DIMENSION  ratmP(3,natm), frcP(3,natm)
c Output arrays:
      DIMENSION  nxatm(*)
      DIMENSION  ratmg(*), frcg(*), amassg(*)
c
c Local declarations:
      DIMENSION  ftmp(3)
c Diagnostics/debugging 
      DATA  iprint / 2 /
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
      nrlxdim = 3*nafrc
      natmg = 0
c
      if( iprint.gt.1 )
     $write(IWR,*) '>>>>> GSELECT begins, nslics=',nslics
c
c  Set up SLIC options
c
      do  iatm=1,natm
c       Copy raw forces into "projected" forces, as default
        frcP(1,iatm) = frc(1,iatm)
        frcP(2,iatm) = frc(2,iatm)
        frcP(3,iatm) = frc(3,iatm)
      enddo
c
c Take care of preparing Unified SLICS
c
      do 200 islic=1,nslics
        islictyp = islics(1,islic)
        if( iprint.gt.1 )
     $   write(IWR,*) 'GSELECT: islic,islictyp=',islic,islictyp
c
        if( islictyp.eq.1 .or. islictyp.eq.2 )then
c         SLIC projections: 1-fold or 2-fold
          if( iprint.gt.1 )
     $    write(IWR,*) 'SLIC/PROJ: projection type=',islictyp
          iatmproj = islics(2,islic)
          iatmank1 = islics(3,islic)
          iatmank2 = islics(4,islic)
 9201     format(1x,'SLIC/',a,': projected atom, anchors=',3i5)
c****** DEV Should reverse cyl/sph order to 1-fold(sph) and 2-fold (cyl)
          if(     islictyp.eq.2 )then
            if( iprint.gt.1 )
     $      write(IWR,9201) 'SLIC/slsph:',iatmproj,iatmank1,iatmank2
            call SLSPH( IWR,natm, islics,islic,
     $                  ratm,frc, ratmP,frcP )
          elseif( islictyp.eq.1 )then
            if( iprint.gt.1 )
     $      write(IWR,9201) 'SLIC/slcyl:',iatmproj,iatmank1,iatmank2
            call SLCYL( IWR,natm, islics,islic,
     $                  ratm,frc, ratmP,frcP )
          endif
c
        elseif( islictyp.eq.3 )then
c         FRAME
          if( iprint.gt.0 )then
            write(IWR,*) 'SLICS/FRAME; slictyp=',islictyp
          endif
          write(IWR,*) 'DEV/FRAME still done outside of SLICS'
c
        elseif( islictyp.eq.4 )then
c         BONDS
          if( iprint.gt.1 )
     $    write(IWR,*) 'SLICS/BOND - length constrained'
          call SLBONDTO( IWR,islic,islics, ratm,frc, ratmP,frcP )
c
        elseif( islictyp.ge.5 .and. islictyp.le.8 )then
c         VGFIX slictyp
          if( iprint.gt.1 )
     $    write(IWR,*) 'SLICS/VGFIX: slictyp=',islictyp
          nvgfix = islics(2,islic)
          if( nvgfix .gt. 0 )then
            isl1 = islics(3,islic)
            islicv = islics(4,islic)
c
            frcvcm = zero
            do  ivgatm=isl1,isl1+nvgfix-1
              iatm = islicats(ivgatm)
              fdotv = zero
              do  id=1,3
                fdotv = fdotv + vslics(id,islicv)*frc(id,iatm)
              enddo
              frcvcm = frcvcm + fdotv
            enddo
            if( iprint.gt.3 )
     $      write(IWR,*) 'VGFIX: CM force along vector=',frcvcm
c
            do  ivgatm=isl1,isl1+nvgfix-1
              iatm = islicats(ivgatm)
              fdotv = zero
              do  id=1,3
                fdotv = fdotv + vslics(id,islicv)*frc(id,iatm)
              enddo
              if( ivgatm.eq.isl1 .and.
     $          ( islictyp.eq.6 .or. islictyp.eq.8 ) )then
c               Allow CM motion along vector
c               Install CM-vec forces into lead atom, suppress other cm
                fdotv = fdotv - frcvcm
              endif
              do  id=1,3
                ratmP(id,iatm) = fdotv*vslics(id,islicv) 
                if( islictyp.eq. 7 )then
c                 Allow motion along vector, suppress normal
                  frcP(id,iatm) = fdotv*vslics(id,islicv) 
                else
c                 islic=3,4,6 
c                 Allow motion norma to vector, remove along vector
                  frcP(id,iatm)= frc(id,iatm) - fdotv*vslics(id,islicv) 
                endif
              enddo
            enddo
          endif
c
        elseif( islictyp.eq.9 )then
c         SLINKY/LAYERS slictyp
          if( iprint.gt.1 )
     $    write(IWR,*) 'SLICS/SLINKY: slictyp=',islictyp
          call SLINKSEL( IWR, islic, islics,islicats,vslics,
     $      ratm,frc, ratmP,frcP )
c
c         End of all SLICS options
        endif

  200 continue
c
c  Pack the relaxation vectors
c
      irlx = 0
      do 300 iatm=1,natm
c
        if( iforce(iatm) .ne. 0 )then
c         This atom will have some sort of relaxation
c
          ityp = itypa(iatm)
          amass = atmass(ityp)
          natmg = natmg + 1
          nxatm(natmg) = 3
c
          if(     iforce(iatm) .eq. 1 )then
c           Standard cartesian relaxation
            do  id=1,3
              irlx = irlx + 1
              ratmg(irlx) = ratm(id,iatm)
c             Yes, frcP is deliberate here, to enable refining anchors
              frcg(irlx)  = frcP(id,iatm)
              amassg(irlx)= amass
            enddo
c
          elseif( iforce(iatm) .eq. 2 )then
cRMV        SLIC-projected - 1-fold sphere or 2-fold cylinder
            write(IWR,*) 'slic projected atom: ',iatm
            do  id=1,3
c             NB: for new slic1/sphere, will also need frcP on ank1
              irlx = irlx + 1
              ratmg(irlx) = ratmP(id,iatm)
              frcg(irlx)  = frcP(id,iatm)
              if( iprint.gt.4 )then
                write(IWR,'(a,2i6,f14.8)') 'coord ',iatm,id, ratmg(irlx)
                write(IWR,'(a,12x,f14.8)') 'force ',frcg(irlx)
              endif
              amassg(irlx)= amass
            enddo
c
          elseif( iforce(iatm) .eq. 4 )then
c           BOND - length constraint
c
            if( iprint.gt.3 )
     $      write(IWR,*) 'BOND: part of fixed bond, iatm= ',iatm
            do  id=1,3
              irlx = irlx + 1
              ratmg(irlx) = ratmP(id,iatm)
              frcg(irlx) = frcP(id,iatm)
              amassg(irlx)= amass
            enddo
c
          elseif( iforce(iatm) .ge. 5 .and. iforce(iatm) .le. 8 )then
c           VGFIX - Constrain atom along vector
c
            if( iprint.gt.3 )
     $      write(IWR,*) 'VGFIX: vector project force on atom ',iatm
            do  id=1,3
              irlx = irlx + 1
              ratmg(irlx) = ratm(id,iatm)
              frcg(irlx) = frcP(id,iatm)
              amassg(irlx)= amass
            enddo
c
          elseif( iforce(iatm) .ge. 9 )then
c           SLINKY/LAYERS - constraint
c
            if( iprint.gt.3 )
     $      write(IWR,*) 'SLINKY/LAYERS optimization atom ',iatm
            do  id=1,3
              irlx = irlx + 1
              ratmg(irlx) = ratmP(id,iatm)
              frcg(irlx) = frcP(id,iatm)
              amassg(irlx)= amass
            enddo
          endif
c
        endif
c
  300 continue
c
      if( iprint.gt.0 )
     $write(IWR,*) '<<<<< GSELECT ends'
c Idiot-check:
      if( nafrc.ne.natmg )then
        call FLGETIERR( IERRFL )
        write(IERRFL,*) '>>>>> nafrc, sum(iforce)=',nafrc,natmg
        call STOPXERR( 'gselect1/number of relaxed atoms mismatch' )
      elseif( irlx .ne. nrlxdim )then
        call FLGETIERR( IERRFL )
        write(IERRFL,*) '>>>>> nrlxdim,irlx=',nrlxdim,irlx
        call STOPXERR( 'gselect2/dim of relaxation vector mismatch' )
      endif
c
c    That's all Folks!
c
      RETURN
      END
