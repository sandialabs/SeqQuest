c***********************************************************************
c  SLICS - selection of slic-related routines
c
c List of routines:
c
c    SLSPH    - spherical projection routine for 1-fold atoms
c    SLSPHB   - spherical de-projection routine for 1-fold atoms
c    SLCYL    - cylindrical projection routine for 2-fold atoms
c    SLCYL    - cylindrical de-projection routine for 2-fold atoms
c
c    SLBONDTO - bond-length fixing projection selection routine
c    SLBONTFR - bond-length fixing de-projection update routine
c
c    SLINKSEL - slinky/layers selection routine for slab relaxation
c    SLINKUPD - slinky/layers update routine for slab relaxation
c
c    UVWCOORD - generate orthogonal u.v.w vectors, from single input vector
c               Presumably z-axis, mapped to w.
c               Then orthgonalize least overlapping x & y to get u,v,w
c
c  ***** WARNING ***** projection routines are NOT periodic-safe
c
c***********************************************************************
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLSPH
c
c
      subroutine SLSPH( IWR,natm, islics,islic,
     $ ratm,frc, ratmP,frcP )
c---------------------------------------------------------------
c Purpose: projects 1-fold coordinates and forces to spherical coordinates
c
c Written: Renee M. Van Ginhoven
c
c Revision history:
c  15Aug17-PAS/2.68: integration into full slics
c  28Sep08-PAS/2.62: cosmetic standardization
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
c
      PARAMETER  ( zero=0.d0 )
cRMV  inputs from gselect
      INTEGER    islics(4,*)
      DIMENSION  ratm(3,natm), frc(3,natm), ratmP(3,natm), frcP(3,natm)
c
c Local declarations:
      DIMENSION  Uvec(3),Vvec(3),Wvec(3), rvec(3),bondvec(3)
      DIMENSION  Unorm(3),Vnorm(3),Wnorm(3),rnorm(3)
      DIMENSION  Fuvec(3),Fvvec(3),Fwvec(3),cartnew(3),Fuwperp(3)
      DIMENSION  Fuwperpnorm(3)
c    
      DATA  Pi / 3.141592654 /
      DATA  degperrad / 57.29577951 /
c
c spherical coordinate projection
c
c First, define new cartesian coordinate system, x,y,z -> u,v,w
c based on positions of projected and anchor atoms
c origin is set to anchor position 1
c Wvec direction is vector from a2 to a1
c Wvec = R(anchor1) -R(anchor2)  -> norm
c Rvec = bondvec = R(atomproj) - R(anchor1)
c Uvec = Uvec' = bondvec - ((bondvec dot Wvec')/|Wvec'|**2)(Wvec')  -> norm
c Vvec' = Wvec' cross Uvec' -> norm
c  initial theta is zero
c  there will be a problem if the a2,a1,pa angle is near 180!! (phi = 0)
c   
c 
c then transform projected atom to spherical coords
c r, phi and theta
c
c R = |rvec|
c W = (rvec dot Wvec')   (assume Wvec' is normalized)
c phi = cos^-1[W/R] 
c theta = sin^-1[V/rsin(phi)] 
c initially, Vvec' = 0, so theta = 0
c
c Project force vector onto new cartesian coordinates
c component along Wvec'  Fw' = (Fatom dot Wvec')
c component along Uvec'  Fu' = (Fatom dot Uvec')
c Fv = F - compU -compW
c transform to spherical
c Fr = sqrt(Fu**2+Fv**2 + Fw**2)
c       or Fr = (Fatom dot rvec)/|rvec|**2
c
c Fphi = sin^-1(FuwcompPerp*sign/Rcoord)
c phi force = magnitude of force in UW plane
c perpendicular to the bond (R) times the direction
c
c Ftheta = sin^-1[Fv/(Fr sin(Fphi))]
c
c update velocity, position
c output files (esp velocity)
c backtransform new coordinates
c
c Unew' = rsin(phi)cos(theta)
c Vnew' = rsin(phi)sin(theta)
c Wnew' = rcos(theta)
c
c return cartesians to original orientation in UNSPHERE
c
c Add option to add 1-fold atom force into anchor atom - 12Aug17
c   Adding the force component projected along the bond, worked well
c   In one test (methanol), adding the full force converged poorly
c   Not clear why the 'full" force propagation would be badly worse
c
      DATA  scalephi, scaletheta / 1.5, 1.5 /
      DATA  ifrcank / 2 /
c ifrcank = 0 raw unmodified anchor force
c         = 1 add full 1-fold force into anchor force
c         = 2 add bond-vector projected 1-fold force into anchor
      DATA  iprint / 2 /
c
c >>>> EXECUTABLE CODE:
c
 9103 format(1x,a,3f14.8)
      iatmproj = islics(2,islic)
      ia1 = islics(3,islic)
      ia2 = islics(4,islic)
      if( iprint.gt.0 )then
        write(IWR,*)'>>>>> slic/slsph: 1-fold sphere projection'
        write(IWR,*)'slic projected atom (1-fold spherical)= #',iatmproj
        write(IWR,*)'  first anchor atom (sphere center) is #',ia1
        write(IWR,*)'  second anchor atom is # ',ia2
      endif
c
c  Retrieve setup parameters
      call SLICPGET( 'FANK  ', pslic, islfank )
      if( islfank.ge.0 .and. islfank.le.2 .and. ifrcank.ne.islfank )then
        if( iprint.gt.0 ) write(IWR,*) 'f-anchor option to ',islfank,
     $    '(0=raw,1=full force,2=bond-projected)'
        ifrcank = islfank
      endif
      call SLICPGET( 'PAR1  ', pslic, isl )
      if( pslic.gt.zero )then
        if( iprint.gt.0 ) write(IWR,9103) 'Change scalephi to',pslic
        scalephi = pslic
      endif
      call SLICPGET( 'PAR2  ', pslic, isl )
      if( pslic.gt.zero )then
        if( iprint.gt.0 ) write(IWR,9103) 'Change scaletheta to',pslic
        scaletheta = pslic
      endif
c set new z->W direction, based on anchor atoms
c z direction is from a2 to a1
      wsq = zero
      do id = 1,3
         Wvec(id) = ratm(id,ia1)-ratm(id,ia2)
         wsq = wsq + (Wvec(id))**2
      enddo
      Wmag = SQRT( wsq )
      if( iprint.gt.3 ) write(IWR,*)'Wmag (reference a1-a2 bond)= ',Wmag
      do id = 1,3
         Wnorm(id) = Wvec(id)/Wmag
      enddo
c
c-------
c set rvec, based on ratom(iatmproj) and anchor1
c first, get bondvec = rvec
c get Rcoord = |bondvec|=bondmag
      bondsq = zero
      do id = 1,3
         bondvec(id)=ratm(id,iatmproj)-ratm(id,ia1)
         bondsq = bondsq +(bondvec(id))**2
      enddo
      bondmag = SQRT( bondsq )
      Rcoord = bondmag
c
c also useful, the r unit vector
      do id = 1,3
         Rnorm(id) = bondvec(id)/Rcoord
      enddo
c     
c     get Wcoord = bond dot Wnorm
      Wcoord = zero
      do id = 1,3
         Wcoord = Wcoord + Wnorm(id)*bondvec(id)
      enddo
      if( iprint.gt.3 ) write(IWR,*) 'Wcoord (b4 step)= ',Wcoord
c
c     now get uvec =normalized component of bondvec perp to Wvec
      Uvecsq = zero
      do id = 1,3
         Uvec(id) = bondvec(id)-Wcoord*Wnorm(id)
         Uvecsq = Uvecsq + (Uvec(id))**2
      enddo
      Ucoord = SQRT( Uvecsq )
      if( iprint.gt.3 ) write(IWR,*) 'Ucoord = ',Ucoord
      do id = 1,3
         Unorm(id) =  Uvec(id)/Ucoord
      enddo
c---- 
c     as a check ,also get Vcoord of iatmprojroj (should be zero!)
c use cross product  Wcrossr
      Vvec(1)=Wnorm(2)*Unorm(3)-Wnorm(3)*Unorm(2)
      Vvec(2)=Wnorm(3)*Unorm(1)-Wnorm(1)*Unorm(3)
      Vvec(3)=Wnorm(1)*Unorm(2)-Wnorm(2)*Unorm(1)
      if( iprint.gt.3 ) write(IWR,9103) 'Vvec =',(Vvec(i),i=1,3)
      Vprmagsq = zero
      do id = 1,3
         Vprmagsq = Vprmagsq +(Vvec(id))**2
      enddo
      Vprmag = SQRT( Vprmagsq )
      if( iprint.gt.3 ) write(IWR,*) 'Vprmag = ',Vprmag
      do id = 1,3
         Vnorm(id)=Vvec(id)/Vprmag
      enddo
c get Vcoord = bond dot Vnorm (should be zero!)
      Vcoord = zero
      do id = 1,3
         Vcoord = Vcoord + bondvec(id)*Vnorm(id)
      enddo
      if( iprint.gt.3 ) write(IWR,*) 'Vcoord (should=0)= ',Vcoord
c and for later we want the "Q" direction, which is just the 
c propjection of R into the UV plane
c  as a check, this initial Q should be the same as U
c
      Qcoordsq = Vcoord**2 + Ucoord**2
      Qcoord = SQRT( Qcoordsq )
      if( iprint.gt.4 )
     $write(IWR,9103) 'PROJ: Qcoord, Ucoord = ',Qcoord, Ucoord
c
c  finish spherical projection of atom coordinates
c
      phi = ACOS( Wcoord/Rcoord )
      if( iprint.gt.2 ) write(IWR,*) 'phi(radians) from acos(z/r)= ',phi
c
c test that the initial theta is zero
      sinphi = SIN( phi )
      if( iprint.gt.2 ) write(IWR,*) 'sinphi = ',sinphi
      ratio = Vcoord/Rcoord
c     if( iprint.gt.4 ) write(IWR,*)'ratio = ',ratio
      theta = ASIN( ratio/sinphi )
c     if( iprint.gt.4 ) write(IWR,*)'theta (should=0) = ',theta
c
c     ***** PROJECT FORCES *****
c
      if( ifrcank.gt.0 )then
c       Propagate outer atom forces into anchor atom
        if(     ifrcank.eq.1 )then
c         Add full 1-fold atom force to anchor (central) atom force
          if( iprint.gt.0 )
     $    write(IWR,*) 'SPHERE: full 1-fold force into central atom'
          do  id=1,3
            frcP(id,ia1) = frc(id,ia1) + frc(id,iatmproj)
          enddo
        elseif( ifrcank.eq.2 )then
c         Add bond-projected 1-fold atom force to anchor atom force
          if( iprint.gt.0 )
     $    write(IWR,*) 'SPHERE: bond 1-fold force into central atom'
          frc12 = Rnorm(1)*frc(1,iatmproj) + Rnorm(2)*frc(2,iatmproj)
     $          + Rnorm(3)*frc(3,iatmproj)
          do  id=1,3
            frcP(id,ia1) = frc(id,ia1) + frc12*Rnorm(id)
          enddo
        else
          call STOPXERR( 'Unknown anchor force update' )
        endif
        if( iprint.gt.1 )then
          write(IWR,9103) 'Raw frc 1-fold =',(frc(i,iatmproj),i=1,3)
          write(IWR,9103) 'Raw frc anchor =',(frc(i,ia1     ),i=1,3)
          write(IWR,9103) 'New frc anchor =',(frcP(i,ia1    ),i=1,3)
        endif
      endif
c
c     first onto rotated cartesians
      Fu = zero
      Fv = zero
      Fw = zero
      do id = 1,3
         Fu = Fu + frc(id,iatmproj)*Unorm(id)
         Fv = Fv + frc(id,iatmproj)*Vnorm(id)
         Fw = Fw + frc(id,iatmproj)*Wnorm(id)
      enddo
c
cFr is the easiest
      Fr = zero
      do id = 1,3
         Fr = Fr + frc(id,iatmproj)*rnorm(id)
      enddo
c
c next get projection of uw plane tangental force
c Fu+Fw component perp to Fr  (No Fv component needed)
c need vector magnitude as input to get Fphi
cget vector
c
      do i=1,3
         Fuwperp(i) = Fu*Unorm(i)+Fw*Wnorm(i)-Fr*Rnorm(i)
      enddo
      Fuwperpmagsq = zero
      do id = 1,3
         Fuwperpmagsq = Fuwperpmagsq +Fuwperp(id)**2
      enddo
      Fuwperpmag = SQRT( Fuwperpmagsq )
      if( iprint.gt.4 ) write(IWR,*) 'Perp to R force = ',Fuwperpmag
c   
c normalize Fuwperp
      do id=1,3
         Fuwperpnorm(id)= Fuwperp(id)/Fuwperpmag
      enddo
c check math...
c
c find sign of phi force Fuwperp 
c Fuwperp dot Wvec
      FperpdotW = zero
      do id=1,3
         FperpdotW = FperpdotW + Fuwperp(id)*Wvec(id)
      enddo
c set direction of force
      fsign = 1.0
      if (FperpdotW.gt.zero)then
         fsign = -1.0
      endif
      if( iprint.gt.3 ) write(IWR,*) 'sign of force = ',fsign
c set Phi force
      Fphi = ASIN( Fuwperpmag*fsign/Rcoord )
      if( iprint.gt.3 )write(IWR,*)'Fphi from asin(Fuwcomp/Rcoord)',Fphi
c if phi too close to 180, revert to cartesian
c      if (abs(phi).gt.3.0)then
c         Fphi = Fuwperpmag
c      endif
      sinphi = SIN( phi )
      ratio = Fv/(Qcoord)
      if( iprint.gt.3 ) write(IWR,9103)'ratio/Fv/Qcoor=',ratio,Fv,Qcoord
c      ratiocheck = Fv/(Rcoord*sinphi)
      Ftheta = ASIN( ratio )
c      Fthetacheck = ASIN( ratiocheck )
      if( iprint.gt.1 )
     $write(IWR,9103)'PROJ: Fr, Fphi, Ftheta =',Fr,Fphi,Ftheta
c,Fthetacheck
c    
c assign spherical forces to force vector
c return projected coordinates and forces in arrays ratmP frcP
      ratmP(1,iatmproj) = Rcoord
      ratmP(2,iatmproj) = phi
      ratmP(3,iatmproj) = theta
c forces, scaled for directions
      frcP(1,iatmproj) = Fr
      frcP(2,iatmproj) = scalephi*Fphi
      frcP(3,iatmproj) = scaletheta*Ftheta
c
      if( iprint.gt.4 )
     $write(IWR,9103)'check forces for r,phi,theta =',Fr,Fphi,Ftheta
c 
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLSPHB
c
c
      subroutine SLSPHB( IWR,natm, islics,islic, 
     $ ratm0,ratmx,ratmP )
c---------------------------------------------------------------
c Purpose: de-projects 1-fold coordinates from spherical back to cartesian
c
c Written: Renee M. Van Ginhoven
c
c Revision history:
c  24Sep08-PAS/2.62: standardization, some minor cleanup
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
c
cRMV  inputs from gselect
      DIMENSION  ratm0(3,natm),ratmx(3,natm),ratmP(3,natm)
      INTEGER    islics(4,*)
c
c Local declarations:
      DIMENSION  Uvec(3),Vvec(3),Wvec(3), rvec(3),bondvec(3)
      DIMENSION  Unorm(3),Vnorm(3),Wnorm(3),rnorm(3)
      DIMENSION  Fuvec(3),Fvvec(3),Fwvec(3),cartnew(3),Fuwperp(3)
      DIMENSION  unrot(3),unshift(3)
c
      DATA  Pi / 3.141592654 /
      DATA  degperrad / 57.29577951 /
      DATA  zero / 0.d0 /
c
c spherical coordinate projection
c
c First, define new cartesian coordinate system,
c based on positions of projected and anchor atoms
c origin is set to anchor position 1
c zvec direction is vector from a2 to a1
c zvec' = R(anchor1) -R(anchor2)  -> norm
c rvec = bondvec = R(atomproj) - R(anchor1)
c xvec' = xvec' = bondvec - ((bondvec dot zvec')/|zvec'|**2)(zvec')  -> norm
c yvec' = zvec' cross xvec' -> norm
c  initial theta is zero
c  there will be a problem if the a2,a1,pa angle is near 180!! (phi = 0)
c
c
c then transform projected atom to spherical coords
c r, phi and theta
c
c r = |rvec|
c z = (rvec dot zvec')   (assume zvec' is normalized)
c phi = cos^-1[z/r] 
c theta = sin^-1[y/rsin(phi)] 
c initially, yvec' = 0, so theta = 0
c
c Project force vector onto new cartesian coordinates
c component along zvec'  Fz' = (Fatom dot zvec')
c component along xvec'  Fx' = (Fatom dot xvec')
c Fy = F - compX -compZ
c transform to spherical
c Fr = sqrt(Fx**2+Fy**2 + Fx**2)
c       or Fr = (Fatom dot rvec)/|rvec|**2
c Fphi = cos^-1(Fz/Fr)
c Ftheta = sin^-1[Fy/(Fr sin(Fphi))]
c
c update velocity, position
c output files (esp velocity)
c backtransform new coordinates
c
c xnew' = rsin(phi)cos(theta)
c ynew' = rsin(phi)sin(theta)
c znew' = rcos(theta)
cc
c
c return cartesians to original orientation
c
      DATA  iprint / 3 /
c
c >>>> EXECUTABLE CODE:
c
 9103 format(1x,a,3f14.8)
      iatmproj = islics(2,islic)
      ia1 = islics(3,islic)
      ia2 = islics(4,islic)
      if( iprint.gt.0 )then
        write(IWR,*)'>>>>> slic/slsphb: 1-fold sphere de-projection'
        write(IWR,*)'slic projected atom (1-fold spherical)= #',iatmproj
        write(IWR,*)'  first anchor atom (sphere center) is #',ia1
        write(IWR,*)'  second anchor atom is # ',ia2
      endif
c set reference W direction, based  on anchor atoms
c W direction is from a2 to a1
      wsq = zero
      do id = 1,3
         Wvec(id) = ratm0(id,ia1)-ratm0(id,ia2)
         wsq = wsq + (Wvec(id))**2
      enddo
      Wmag = SQRT( wsq )
      do id = 1,3
         Wnorm(id) = Wvec(id)/Wmag
      enddo
c write out Wprime vector as check
      if( iprint.gt.3 ) write(IWR,*)'Wmag (reference a1-a2 bond)= ',Wmag
      if( iprint.gt.3 ) write(IWR,9103)'Wvec = Wmag * ',(Wnorm(i),i=1,3)
c
c-------
c set rvec, based on ratom(iatmproj) and anchor1
c first, get bondvec = rvec
c get Rcoord = |bondvec|=bondmag
c NOTE for nested scheme, the bond vector will have to
c be kept from the original coordinate transformation
c rather than re-calculated here
      bondsq = zero
      do id = 1,3
         bondvec(id)=ratm0(id,iatmproj)-ratm0(id,ia1)
         bondsq = bondsq +(bondvec(id))**2
      enddo
      bondmag = SQRT( bondsq )
      Rcoord0 = bondmag
      if( iprint.gt.3 ) write(IWR,9103)'Rcoord0 = ',Rcoord0
      if( iprint.gt.3 ) write(IWR,9103)'bondvec =',(bondvec(i),i=1,3)
c also useful, the r unit vector
      do id = 1,3
         rnorm(id) = bondvec(id)/Rcoord0
      enddo
c
c      
c     get Wcoord = bond dot Wnorm
      Wcoord0 = zero
      do id = 1,3
         Wcoord0 = Wcoord0 + Wnorm(id)*bondvec(id)
      enddo
c
c     now get xvec =normaliWed component of bondvec perp to Wvec
      Uvecsq = zero
      do id = 1,3
         Uvec(id) = bondvec(id)-Wcoord0*Wnorm(id)
         Uvecsq = Uvecsq + (Uvec(id))**2
      enddo
      Ucoord0 = SQRT( Uvecsq )
c
      do id = 1,3
         Unorm(id) =  Uvec(id)/Ucoord0
      enddo
c
c---- 
c     as a check ,also get yprimecoord of iatmprojroj (should be zero!)
c use cross product  Wcrossr
      Vvec(1)=Wnorm(2)*Unorm(3)-Wnorm(3)*Unorm(2)
      Vvec(2)=Wnorm(3)*Unorm(1)-Wnorm(1)*Unorm(3)
      Vvec(3)=Wnorm(1)*Unorm(2)-Wnorm(2)*Unorm(1)
      if( iprint.gt.3 ) write(IWR,9103) 'Vvec =',(Vvec(i),i=1,3)
      Vmagsq = zero
      do id = 1,3
         Vmagsq = Vmagsq +(Vvec(id))**2
      enddo
      Vmag = SQRT( Vmagsq )
      if( iprint.gt.3 ) write(IWR,*) 'Vmag = ',Vmag
      do id = 1,3
         Vnorm(id)=Vvec(id)/Vmag
      enddo
c get Vcoord = bond dot Vnorm (should be zero!)
      Vcoord0 = zero
      do id = 1,3
         Vcoord0 = Vcoord0 + bondvec(id)*Vnorm(id)
      enddo
      if( iprint.gt.3 ) write(IWR,*) 'Vcoord0 = ', Vcoord0
c
c check dot products
      VdotW = zero
      VdotU = zero
      WdotU = zero
      do id = 1,3
         VdotW = VdotW + Vvec(id)*Wvec(id)
         VdotU = VdotU + Vvec(id)*Uvec(id)
         WdotU = WdotU + Wvec(id)*Uvec(id)
      enddo
      if( iprint.gt.3 ) write(IWR,*) 'VdotW = ',VdotW
      if( iprint.gt.3 ) write(IWR,*) 'VdotU = ',VdotU
      if( iprint.gt.3 ) write(IWR,*) 'WdotU = ',WdotU
c-----
c     check maximum step for theta and phi
c    need this to prevent excessive twirling
c     
c      angvelmax = 10.
c      angvelmaxn = -10.
c      if (theta .gt. angvelmax) then
c         theta = angvelmax
c      elseif (theta .lt. angvelmaxn)then
c         theta = angvelmaxn
c      endif
c      if (phi .gt. angvelmax) then
c         phi = angvelmax
c      elseif (phi .lt. angvelmaxn)then
c         phi = angvelmaxn
c      endif
c---- 
c     Backtransformation of coordinates
c     to cartesian
c
      Rcoord = ratmx(1,iatmproj)
      phi    = ratmx(2,iatmproj)
      theta  = ratmx(3,iatmproj)
      if( iprint.gt.1 )
     $  write(IWR,9103) 'Spherical coordinates: R,phi,theta =',
     $                   Rcoord,phi,theta
      cartnew(1) = Rcoord*SIN( phi )*COS( theta )
      cartnew(2) = Rcoord*SIN( phi )*SIN( theta )
      cartnew(3) = Rcoord*COS( phi )
      if( iprint.gt.2 )then
        write(IWR,9103) 'new coordinates                   ',
     $                  (cartnew(i),i=1,3)
        write(IWR,9103) 'original shifted coordinates      ',
     $                   Ucoord0,Vcoord0,Wcoord0
        write(IWR,9103) 'original unshifted coordinates',
     $                  (Ratm0(i,iatmproj),i=1,3)
      endif
c rotate to original reference frame
      unrot(1) = cartnew(1)*Unorm(1) + cartnew(2)*Vnorm(1)  
     &         + cartnew(3)*Wnorm(1)
      unrot(2) = cartnew(1)*Unorm(2) + cartnew(2)*Vnorm(2)  
     &         + cartnew(3)*Wnorm(2)
      unrot(3) = cartnew(1)*Unorm(3) + cartnew(2)*Vnorm(3)  
     &         + cartnew(3)*Wnorm(3)
c check rotation
      if( iprint.gt.2 )
     $write(IWR,9103)'unrotated coords',(unrot(i),i=1,3)
c
c reset origin
      do i = 1,3
         unshift(i)=unrot(i)+Ratm0(i,ia1)
         if( iprint.gt.4 ) write(IWR,*)ia1,Ratm0(i,ia1)
      enddo
      if( iprint.gt.4 )
     $write(IWR,9103)'new coords',(unshift(i),i=1,3)
c
c return cartesian coordinates to ratmx array
      do id = 1,3
         ratmx(id,iatmproj) = unshift(id)
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLCYL
c
c
      subroutine SLCYL( IWR,natm, islics,islic, 
     $ ratm,frc, ratmP,frcP )
c---------------------------------------------------------------
c Purpose: projects 2-fold coordinates and forces to cylindrical coordinates
c
c Written: Renee M. Van Ginhoven
c
c Revision history:
c  24Sep08-PAS/2.62: cosmetics and some minor cleanip
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
c
cRMV  inputs from gselect
      DIMENSION  ratm(3,natm), frc(3,natm), ratmP(3,natm),frcP(3,natm)
      INTEGER    islics(4,*)
c
c islics contains the projection type, and projected atom and anchor atoms
c 1 = cylindrical projection
c 2 = spherical projection
c
c Local declarations:
      DIMENSION  Uvec(3),Vvec(3),Wvec(3), rvec(3),bondvec(3)
      DIMENSION  Unorm(3),Vnorm(3),Wnorm(3),rnorm(3)
      DIMENSION  Fuvec(3),Fvvec(3),Fwvec(3),cartnew(3),Fuwperp(3)
      DIMENSION  Fuwperpnorm(3)
      DIMENSION  ypr(3)
      DIMENSION  Frvec(3)
      DIMENSION  unrot(3),unshift(3)
c
      DATA  Pi / 3.141592654 /
      DATA  degperrad / 57.29577951 /
      DATA  zero /0.d0/
c
c cylindrical coordinate projection
c
c first, define new cartesian coordinate system,
c based on positions of projected and anchor atoms
c  Wvec' = R(anchor2)-R(anchor1)  -> norm
c  bondvec = R(atomproj)-R(anchor1)
c  Uvec' = bondvec - ((bondvec dot Wvec')/|Wvec'|**2)(zvec')  -> norm
c  Vvec' = Wvec' cross Uvec'  -> norm
cdon't actually need Vvec here
c
c then transform projected atom to polar coordinates
c W = (bondvec dot zvec')/|Wvec'|**2
c r = |Uvec'|
c r = U'/cos(theta),  theta = 0; r = U'
c theta = tan^-1(v/u)   or  theta = cos^-1(u/r)
c initially, Vvec' = 0, so theta = 0
c 
c Project force vector onto new cartesian coordinates
c component along W, component along U (=r0)
c component along y; Fy = F- compU - compW
c  transform to polar
c Ftheta = sin^-1(V/r) 
c update velocity, position
c output files (esp velocity)
c backtransform new coordinates
c
      DATA  iprint / 2 /
c
c >>>> EXECUTABLE CODE:
c
 9103 format(1x,a,3f14.8)
      iatmproj = islics(2,islic)
      ia1 = islics(3,islic)
      ia2 = islics(4,islic)
      if( iprint.gt.0 )then
        write(IWR,*)'>>>>> slic/slcyl: 2-fold cylinder projection'
        write(IWR,*)'slic projected atom (2-fold cylinder)= #',iatmproj
        write(IWR,*)'  first anchor atom is #',ia1
        write(IWR,*)'  second anchor atom is # ',ia2
      endif
c set new z direction, based  on anchor atoms
      wsq = zero
      do id = 1,3
         Wvec(id) = ratm(id,ia2)-ratm(id,ia1)
         wsq = wsq + (Wvec(id))**2
      enddo
      Wmag = SQRT( wsq )
      if( iprint.gt.3 ) write(IWR,*) 'Wmag = ' ,Wmag
      do id = 1,3
         Wnorm(id) = Wvec(id)/Wmag
      enddo
c write out W vector as check
      if( iprint.gt.3 )write(IWR,9103) 'Wvec = Wmag *',(Wnorm(i),i=1,3)
c-------
c set rvec, based on ratom(iatmproj) and new w
c first, get bondvec
      bondsq = zero
      do id = 1,3
         bondvec(id)=ratm(id,iatmproj)-ratm(id,ia1)
         bondsq = bondsq +(bondvec(id))**2
      enddo
      bondmag = SQRT( bondsq )
      if( iprint.gt.3 ) write(IWR,9103) 'bondmag =',bondmag
      if( iprint.gt.3 ) write(IWR,9103) 'bondvec =',(bondvec(i),i=1,3)
c     get Wcoord = bond dot wnorm
      Wcoord = zero
      do id = 1,3
         Wcoord = Wcoord + Wnorm(id)*bondvec(id)
      enddo
      if( iprint.gt.3 ) write(IWR,*) 'Wcoord (b4 step) = ',Wcoord
      rvecsq = zero
c     now get rvec and Rcoord
      do id = 1,3
         rvec(id) = bondvec(id)-Wcoord*Wnorm(id)
         rvecsq = rvecsq + (rvec(id))**2
      enddo
c     rvecmag = Rcoord = xprime
      Rcoord = SQRT( rvecsq )
      do id = 1,3
         rnorm(id) = rvec(id)/Rcoord
      enddo
      if( iprint.gt.3 ) write(IWR,9103) 'Rcoord =',Rcoord
      if( iprint.gt.3 ) write(IWR,9103) 'rnorm =',(rnorm(i),i=1,3)
      theta = zero
c---- 
c     as a check ,also get yprimecoord of iatmprojroj (should be zero!)
c use cross product  zcrossr
      Vvec(1)=Wcoord*Wnorm(2)*rvec(3)-Wcoord*wnorm(3)*rvec(2)
      Vvec(2)=Wcoord*Wnorm(3)*rvec(1)-Wcoord*wnorm(1)*rvec(3)
      Vvec(3)=Wcoord*wnorm(1)*rvec(2)-Wcoord*wnorm(2)*rvec(1)
      if( iprint.gt.3 ) write(IWR,9103) 'Vvec =',(Vvec(i),i=1,3)
      Vmagsq = zero
      do id = 1,3
         Vmagsq = Vmagsq +(Vvec(id))**2
      enddo
      vmag = SQRT( vmagsq )
      if( iprint.gt.3 ) write(IWR,*) 'Vmag = ',Vmag
      do id = 1,3
         Vnorm(id)=Vvec(id)/Vmag
      enddo
c get ycoord = bond dot yprnorm (should be zero!)
      Vcoord = zero
      do id = 1,3
         Vcoord = Vcoord + bondvec(id)*Vnorm(id)
      enddo
      if( iprint.gt.3 ) write(IWR,*) 'Vcoord = ', Vcoord
      theta = ASIN( Vcoord/rcoord )
      if( iprint.gt.3 ) write(IWR,*) 'THETA = ',theta
c
c check dot products
      VdotW = zero
      VdotU = zero
      WdotU = zero
      do id = 1,3
         VdotW = VdotW + Vvec(id)*Wvec(id)
         VdotU = VdotU + Vvec(id)*rvec(id)
         WdotU = WdotU + Wvec(id)*rvec(id)
      enddo
      if( iprint.gt.3 ) write(IWR,*)'VdotW = ',VdotW
      if( iprint.gt.3 ) write(IWR,*)'VdotU = ',VdotU
      if( iprint.gt.3 ) write(IWR,*)'WdotU = ',WdotU
c-------
c     project forces
c     
      Fw = zero
      do id = 1,3
         Fw =Fw + Frc(id,iatmproj)*Wnorm(id)
      enddo
      do id=1,3
         Fwvec(id) =  Fw*wnorm(id) 
      enddo
c      
      Fr = zero
      do id = 1,3
         Fr =Fr + Frc(id,iatmproj)*rnorm(id)
      enddo      
      do id=1,3
         Frvec(id) =  Fr*rnorm(id) 
      enddo
c     
c     
      Fv = zero
      Ftheta = zero
      do id = 1,3
        Fv = Fv + Frc(id,iatmproj)*Vnorm(id)
      enddo
c 
c
      FthetaCH = ASIN( Fv/Rcoord )
      if( iprint.gt.3 ) write(IWR,*) 'check force magnitudes'
      cartfsq = zero
      do i=1,3
         cartfsq=cartfsq+frc(i,iatmproj)*frc(i,iatmproj)
      enddo
      cartF = SQRT( cartfsq )
      polfsq = zero
      rotcartfsq = Fr*Fr +Fv*Fv+ Fw*Fw
      rotcartf = SQRT( rotcartfsq )
      polfsq = Fr*Fr + FthetaCH*FthetaCH+Fw*Fw
      polf = SQRT( polfsq )
      if( iprint.gt.4 )
     $write(IWR,9103)'cart and pol forces',cartf,rotcartf,polf
c
c      
c assign polar forces to force vector
c----------------------------------------------
c option 1: Fr,Ftheta,Fz with no modifications
c----------------------------------------------
c option 2: use Fy instead of Ftheta -  good if steps small enough
c      Ftheta = Fy
c----------------------------------------------
c option 3: limit maximum Fy/Rcoord ratio
c     Fy/Rcoord =ratio
c     ratiomax = 0.173
c  
c      if (ratio.gt. ratiomax)then
c         ratio = ratiomax
c      endif
c      Ftheta = ASIN( ratio )   
c-----------------------------------------------
c option 4: use a constant for Rcoord (i.e. unit vector = 1)
c and use a max cutoff for Fy to restrict dtheta to <30deg
c Therefore Fy/rfixed < 0.5       Fy < 0.008
      rfixed = 0.05
      Fvmax = 0.03
      Fvmaxn =-0.03
c
c
c-----------------------------------------------------
c option5: rescalable rfixed:  initial value based on initial Fy
c  think some more about this
c some of this function will be taked care of by asde
c
c
c
      Ftheta = 2.0*ASIN( Fv/Rcoord )
c      if (Fy.gt. Fymax)then
c         Fy = Fymax
c      endif
c      if (Fy .lt. Fymaxn) then
c         Fy = Fymaxn
c      endif
c
c      tinyforce = 0.001
c      if (abs(Fy) .lt. tinyforce)then
c         rfixed = 0.005
c      endif
c       Ftheta = ASIN( Fy/rfixed )  
c-----------------------------------------------------
c set forces
c
      if( iprint.gt.2 ) write(IWR,9103)'check cyl forces=',Fr,Ftheta,Fw
      if( iprint.gt.2 )
     $write(IWR,9103) 'check Fxyz',(frc(i,iatmproj),i=1,3)
c set atom forces to projected forces
      frcP(1,iatmproj)=Fr
      frcP(2,iatmproj)=Ftheta
      frcP(3,iatmproj)=Fw
C set atom coordinates to projected coordinates
      ratmP(1,iatmproj)=Rcoord
      ratmP(2,iatmproj)=theta
      ratmP(3,iatmproj)=Wcoord
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLCYLB
c
c
      subroutine SLCYLB( IWR,natm, islics,islic,
     $ ratm0,ratmx,ratmP )
c---------------------------------------------------------------
c Purpose: deproject 2-fold curvelinear coodrdinates back to cartesian
c
c Written: Renee M. Van Ginhoven
c
c Revision history:
c  24Sep08-PAS/2.62: some cosmetic standardization
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
c
cRMV  inputs from gupdate
      DIMENSION  ratm0(3,natm), ratmx(3,natm), ratmP(3,natm)
      INTEGER    islics(4,*)
c 
c Local declarations:
      DIMENSION  Uvec(3),Vvec(3),Wvec(3), rvec(3),bondvec(3)
      DIMENSION  Unorm(3),Vnorm(3),Wnorm(3),rnorm(3)
      DIMENSION  Fuvec(3),Fvvec(3),Fwvec(3),cartnew(3),Fuwperp(3)
      DIMENSION  Fuwperpnorm(3)
      DIMENSION  Frvec(3)
      DIMENSION  unrot(3),unshift(3)
c    
      DATA  Pi / 3.141592654 /
      DATA  degperrad / 57.29577951 /
      DATA  zero / 0.d0 /
c
c cylindrical coordinate projection
c
c first, define new cartesian coordinate system,
c based on positions of projected and anchor atoms
c  zvec' = R(anchor2)-R(anchor1)  -> norm
c  bondvec = R(atomproj)-R(anchor1)
c  xvec' = bondvec - ((bondvec dot zvec')/|zvec'|**2)(zvec')  -> norm
c  yvec' = zvec' cross xvec'  -> norm
cdon't actually need yvec here
c
c then transform projected atom to polar coordinates
c z = (bondvec dot zvec')/|zvec'|**2
c r = |xvec'|
c r = x'/cos(theta),  theta = 0; r = x'
c theta = tan^-1(y/x)   or  theta = cos^-1(x/r)
c initially, yvec' = 0, so theta = 0
c
c Project force vector onto new cartesian coordinates
c component along z, component along x (=r0)
c component along y; Fy = F- compX - compZ
c  transform to polar
c Ftheta = sin^-1(y/r) 
c update velocity, position
c output files (esp velocity)
c backtransform new coordinates
c
c--------
c set new z direction, based  on anchor atoms
c nested update is taked care of iin gupdate.  The projection vectors
c in this routine are based on the ratm0 (master array) coordinates
c If this is a nested update, the master anchor atoms
c will have been updated previously.
c and yes, this messes with the symmetry options
c
      DATA  iprint / 2 /
c
c >>>> EXECUTABLE CODE:
c
 9103 format(1x,a,3f14.8)
      iatmproj = islics(2,islic)
      ia1 = islics(3,islic)
      ia2 = islics(4,islic)
      if( iprint.gt.0 )then
        write(IWR,*)'>>>>> slic/slcylb: 2-fold cylinder de-projection'
        write(IWR,*)'slic projected atom (2-fold cylinder)= #',iatmproj
        write(IWR,*)'  first anchor atom is #',ia1
        write(IWR,*)'  second anchor atom is # ',ia2
      endif
c
c set new z direction, based  on anchor atoms
      wsq = zero
      do id = 1,3
         Wvec(id) = ratm0(id,ia2)-ratm0(id,ia1)
         wsq = wsq + (Wvec(id))**2
      enddo
      Wmag = SQRT( wsq )
      if( iprint.gt.3 ) write(IWR,9103)'Wmag = ' ,Wmag
      do id = 1,3
         Wnorm(id) = Wvec(id)/Wmag
      enddo
c write out zprime vector as check
      if( iprint.gt.3 ) write(IWR,9103) 'Wvec = Wnag *',(Wnorm(i),i=1,3)
c-------
c set rvec, based on ratom(iatmproj) and new W
c first, get bondvec
      bondsq = zero
      do id = 1,3
         bondvec(id)=ratm0(id,iatmproj)-ratm0(id,ia1)
         bondsq = bondsq +(bondvec(id))**2
      enddo
      bondmag = SQRT( bondsq )
      if( iprint.gt.3 ) write(IWR,9103) 'bondmag =',bondmag
      if( iprint.gt.3 ) write(IWR,9103) 'bondvec =',(bondvec(i),i=1,3)
c     get Wcoord = bond dot znorm
      Wcoord = zero
      do id = 1,3
         Wcoord = Wcoord + Wnorm(id)*bondvec(id)
      enddo
      if( iprint.gt.3 ) write(IWR,*) 'Wcoord (b4 step)= ',Wcoord
      if( iprint.gt.3 ) write(IWR,*) 'check coordinate transfer'
      if( iprint.gt.3 )
     $write(IWR,*)'Wcoord from previous projection',ratmP(3,iatmproj)
c
      rvecsq = zero
c     now get rvec and Rcoord
      do id = 1,3
         rvec(id) = bondvec(id)-Wcoord*Wnorm(id)
         rvecsq = rvecsq + (rvec(id))**2
      enddo
c     rvecmag = Rcoord = xprime
      Rcoord = SQRT( rvecsq )
      do id = 1,3
         rnorm(id) = rvec(id)/Rcoord
      enddo
      if( iprint.gt.3 ) write(IWR,9103) 'Rcoord =',Rcoord
      if( iprint.gt.3 ) write(IWR,9103) 'rnorm =',(rnorm(i),i=1,3)
      theta = zero
c---- 
c     as a check ,also get yprimecoord of iatmprojroj (should be zero!)
c use cross product  zcrossr
      Vvec(1)=Wcoord*Wnorm(2)*rvec(3)-Wcoord*Wnorm(3)*rvec(2)
      Vvec(2)=Wcoord*Wnorm(3)*rvec(1)-Wcoord*Wnorm(1)*rvec(3)
      Vvec(3)=Wcoord*Wnorm(1)*rvec(2)-Wcoord*Wnorm(2)*rvec(1)
      if( iprint.gt.3 ) write(IWR,9103) 'Vvec =',(Vvec(i),i=1,3)
      Vmagsq = zero
      do id = 1,3
         Vmagsq = Vmagsq +(Vvec(id))**2
      enddo
      Vmag = SQRT( Vmagsq )
      if( iprint.gt.3 ) write(IWR,9103) 'Vmag = ',Vmag
      do id = 1,3
         Vnorm(id)=Vvec(id)/Vmag
      enddo
c get ycoord = bond dot yprnorm (should be zero!)
      Vcoord = zero
      do id = 1,3
         Vcoord = Vcoord + bondvec(id)*Vnorm(id)
      enddo
      if( iprint.gt.3 ) write(IWR,*) 'Vcoord = ', Vcoord
      theta = ASIN( Vcoord/rcoord )
      if( iprint.gt.3 ) write(IWR,*) 'THETA (should=0)= ',theta
c
c check dot products
      VdotW = zero
      VdotU = zero
      WdotU = zero
      do id = 1,3
         VdotW = VdotW + Vvec(id)*Wvec(id)
         VdotU = VdotU + Vvec(id)*rvec(id)
         WdotU = WdotU + Wvec(id)*rvec(id)
      enddo
      if( iprint.gt.3 ) write(IWR,9103) 'VdotW = ',VdotW
      if( iprint.gt.3 ) write(IWR,9103) 'VdotU = ',VdotU
      if( iprint.gt.3 ) write(IWR,9103) 'WdotU = ',WdotU
c-----
c     update polar coordinates of iatmprojroj
      if( iprint.gt.3 ) write(IWR,*) 'wcoord (b4 update)= ',Wcoord
c     
c       theta1 = theta
      Rcoord = ratmx(1,iatmproj)
      theta  = ratmx(2,iatmproj)
      Wcoord = ratmx(3,iatmproj)
      if( iprint.gt.1 )
     $write(IWR,9103) 'Rcoord,theta,wcoord =',rcoord,theta,wcoord
c
c     check maximum step for theta
      angvelmax = 0.25
      angvelmaxn = -0.25
      if (theta .gt. angvelmax) then
         theta = angvelmax
      elseif (theta .lt. angvelmaxn)then
         theta = angvelmaxn
      endif  
c
c---- 
c     Backtransformation of coordinates
c     to cartesian
      if( iprint.gt.1 )
     $write(IWR,9103) 'polar coordinates, R,theta,W=',Rcoord,theta,Wcoord
      cartnew(1) = Rcoord*COS( theta )
      cartnew(2) = Rcoord*SIN( theta )
c      cartnew(1) = Rcoord
c      cartnew(2) = theta
      cartnew(3) = Wcoord
      if( iprint.gt.2 )then
        write(IWR,9103) 'new coordinates',(cartnew(i),i=1,3)
        write(IWR,9103) 'original shifted coordinates',
     $                   Rcoord,Vcoord,Wcoord
        write(IWR,9103) 'original unshifted coordinates',
     $                  (Ratm0(i,iatmproj),i=1,3)
      endif
c rotate to original reference frame
      unrot(1) = cartnew(1)*rnorm(1) + cartnew(2)*Vnorm(1)  
     &         + cartnew(3)*Wnorm(1)
      unrot(2) = cartnew(1)*rnorm(2) + cartnew(2)*Vnorm(2)  
     &         + cartnew(3)*Wnorm(2)
      unrot(3) = cartnew(1)*rnorm(3) + cartnew(2)*Vnorm(3)  
     &         + cartnew(3)*Wnorm(3)
c check rotation
      if( iprint.gt.3 )
     $ write(IWR,9103) 'unrotated coordinates',(unrot(i),i=1,3)
c
c reset origin
      do i = 1,3
         unshift(i)=unrot(i)+Ratm0(i,ia1)
c         write(IWR,*)ia1,Ratom(ia1,i)
      enddo
      if( iprint.gt.2 )
     $ write(IWR,9103) 'new coordinates',(unshift(i),i=1,3)
c
c return cartesian coordinates to ratmx array
      do id = 1,3
         ratmx(id,iatmproj) = unshift(id)
      enddo
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLBONDTO
c
c
      subroutine SLBONDTO( IWR, islic, islics,
     $ ratm,frc, ratmP,frcP )
c---------------------------------------------------------------
c Purpose: transform into bond coordinates
c Written: 1-Aug-2017, for v2.68
c Revisions: none
c---------------------------------------------------------------
c
c ratm = original cartesian coordinates
c ratmP = original projected coordinates (CM/12)
c frc = original cartesian forces
c frcP = modified projected forces (bond-length conserving)
c
c WARNING: This takes no account of periodicity
c          This will fail if bonds across cell boundaries
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  (zero=0.d0,two=2.d0)
c
      DIMENSION  islics(4,*)
      DIMENSION  ratm(3,*), frc(3,*)
      DIMENSION  ratmP(3,*), frcP(3,*)
c scratch
      DIMENSION  x12(3)
c Local diagnostics
      DATA  iprint / 2 /
c
c >>>> EXECUTABLE CODE:
c
      if( islics(1,islic).ne.4 ) RETURN
c
      iatm1 = islics(2,islic)
      iatm2 = islics(3,islic)
      if( iprint.gt.1 )
     $write(IWR,*) 'BOND-select: fixed length bond atoms=',iatm1,iatm2
c
c     Put R-CM into R-iatm1, and F-CM into F-iatm1
c     Put R-rot into R-iatm2, and F-rot into F-iatm22
      do  id=1,3
        ratmP(id,iatm1) = (ratm(id,iatm2)+ratm(id,iatm1))/two
        ratmP(id,iatm2) = (ratm(id,iatm2)-ratm(id,iatm1))/two
        frcP(id,iatm1) = (frc(id,iatm2)+frc(id,iatm1))/two
        frcP(id,iatm2) = (frc(id,iatm2)-frc(id,iatm1))/two
      enddo
c
c     Project out the interatomic-force along the bond vector
      r12 = SQRT(ratmP(1,iatm2)**2+ratmP(2,iatm2)**2+ratmP(3,iatm2)**2)
      do  id=1,3
        x12(id) = ratmP(id,iatm2)/r12
      enddo
      f12 = x12(1)*frcP(1,iatm2) + x12(2)*frcP(2,iatm2)
     $    + x12(3)*frcP(3,iatm2)
      do id=1,3
        frcP(id,iatm2) = frcP(id,iatm2) - f12*x12(id)
      enddo
c
      if( iprint.gt.1 ) write(IWR,914) 'Bond length(au)=',two*r12
      if( r12 .gt. 5.0 )then
        write(IWR,*) '***** WARNING: very long BOND (bohr)=',two*r12
      endif
      if( iprint.gt.3 )then
        write(IWR,*) 'BONDTO recap: iatm, R(iatm), Rpr(iatm);'
  913   format(3x,i5,5x,3f11.4,5x,3f11.4)
  914   format(3x,a,3f12.6)
        write(IWR,914) 'BOND axis=',x12
        write(IWR,*)'  SLIC/BOND recap: iatm, R(iatm), Rpr(iatm);'
        write(IWR,913)iatm1,(ratm(i,iatm1),i=1,3),(ratmP(i,iatm1),i=1,3)
        write(IWR,913)iatm2,(ratm(i,iatm2),i=1,3),(ratmP(i,iatm2),i=1,3)
        write(IWR,*)'  SLIN/BOND recap: iatm, F(iatm), Fpr(iatm);'
        write(IWR,913)iatm1,(frc(i,iatm1),i=1,3),(frcP(i,iatm1),i=1,3)
        write(IWR,913)iatm2,(frc(i,iatm2),i=1,3),(frcP(i,iatm2),i=1,3)
      endif
c  
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLBONDFR
c
c
      subroutine SLBONDFR( IWR, islic, islics,
     $ ratm0, ratmx, ratmP )
c---------------------------------------------------------------
c Purpose: transform back from bond coordinates
c Written:  1-Aug-2017, for v2.68
c Revisions: none
c---------------------------------------------------------------
c
c ratm0 = original cartesian coordinates
c ratmP = original projected coordinates (CM/12)
c ratmx = new project coordinates
c
c We use a simpler approach here
c  Rathe than doing a radial update (formally along sphere),
c  we do a cartesian update, and then contract the bond to fixed length
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  (zero=0.d0,two=2.d0)
c
      DIMENSION  islics(4,*)
      DIMENSION  ratm0(3,*)
      DIMENSION  ratmx(3,*)
      DIMENSION  ratmP(3,*)
c Local declarations
      DIMENSION  r1(3),r2(3)
c Local diagnostic
      DATA  iprint / 2 /
c
c >>>> EXECUTABLE CODE:
c
      if( islics(1,islic).ne.4 ) RETURN
c
      iatm1 = islics(2,islic)
      iatm2 = islics(3,islic)
      if( iprint.gt.1 )
     $ write(IWR,*)'BONDS-update: bonded atoms=',iatm1,iatm2
  913 format(3x,i5,5x,3f11.4,5x,3f11.4)
      if( iprint.gt.4 )then
       write(IWR,*)'  BOND recap: iatm, Rpr(iatm), Rpr-new(iatm);'
       write(IWR,913)iatm1,(ratmP(i,iatm1),i=1,3),(ratmx(i,iatm1),i=1,3)
       write(IWR,913)iatm2,(ratmP(i,iatm2),i=1,3),(ratmx(i,iatm2),i=1,3)
      endif
c  
      d12old = SQRT( ratmP(1,iatm2)**2 + ratmP(2,iatm2)**2
     $             + ratmP(3,iatm2)**2 ) * two
      d12new = SQRT( ratmx(1,iatm2)**2 + ratmx(2,iatm2)**2
     $             + ratmx(3,iatm2)**2 ) * two
      if( iprint.gt.2 )
     $ write(IWR,'(a,2f10.6)') 'Old,new bond=',d12old,d12new
c
c     Back transorm into native carteeian units
c     conserving bond distance
      fac = d12old/d12new
      do  id=1,3
        r1(id) = (ratmx(id,iatm1)-fac*ratmx(id,iatm2))
        r2(id) = (ratmx(id,iatm1)+fac*ratmx(id,iatm2))
      enddo
      do  id=1,3
        ratmx(id,iatm1) = r1(id)
        ratmx(id,iatm2) = r2(id)
      enddo
      if( iprint.gt.4 )then
       write(IWR,*)'  BOND recap: iatm, R0(iatm), Rxyz-new(iatm);'
       write(IWR,913)iatm1,(ratm0(i,iatm1),i=1,3),(ratmx(i,iatm1),i=1,3)
       write(IWR,913)iatm2,(ratm0(i,iatm2),i=1,3),(ratmx(i,iatm2),i=1,3)
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
c End
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLINKSEL
c
c
      subroutine SLINKSEL( IWR, islic, islics,islicats,vslics,
     $ ratm,frc, ratmP,frcP )
c---------------------------------------------------------------
c Purpose: transform into slinky/layers (slabs)
c Written: 2-Aug-2017, for v2.68
c Revisions: none
c---------------------------------------------------------------
c
c ratm = original cartesian coordinates
c ratmP = original projected coordinates (CM/12)
c frc = original cartesian forces
c frcP = modified projected forces (bond-length conserving)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  (zero=0.d0,two=2.d0)
c
      DIMENSION  islics(4,*), islicats(*), vslics(3,*)
c
      DIMENSION  ratm(3,*), frc(3,*)
      DIMENSION  ratmP(3,*), frcP(3,*)
c Local declarations
      DIMENSION  ftmp(3), uvw(3,3)
c Local diagnostic
      DATA  iprint / 2 /
c
c >>>> EXECUTABLE CODE:
c
      if( iprint.gt.1 ) write(IWR,*) 'SLINK-select (d-layer) start'
c
      if( islics(1,islic).ne.9 ) RETURN
c
      nslink = islics(2,islic)
      if( iprint.gt.2 ) write(IWR,*) 'SLINK number atoms=',nslink
      if( nslink.le.0 ) RETURN
c
      isl1 = islics(3,islic)
      isln = isl1 + nslink - 1
      islicv = islics(4,islic)
      iatm1 = islicats(isl1)
c
c     Create coordinate system around slinky axis
      call UVWCOORD( uvw, vslics(1,islicv) )
c
      if( iprint.gt.4 )then
        write(IWR,*) 'SLINK vector and projected coordinates'
        write(IWR,'(a,3f10.6)') 'slinky =',(vslics(i,islicv),i=1,3)
        write(IWR,'(a,3f10.6)') '  u =',(uvw(i,1),i=1,3)
        write(IWR,'(a,3f10.6)') '  v =',(uvw(i,2),i=1,3)
        write(IWR,'(a,3f10.6)') '  w =',(uvw(i,3),i=1,3)
      endif
c
c     Do the transformation in layer separation in new coordinates
c     Note that the final (basement) atom remains cartesian coordinates
      call MKZERO( 3, ftmp )
      iatmlast = 0
      do 200 islink=isl1,isln
        iatm = islicats(islink)
        do  iuvw=1,3
          ratmP(iuvw,iatm) = DDOT( 3, uvw(1,iuvw),1, ratm(1,iatm),1 )
          frcP(iuvw,iatm)  = DDOT( 3, uvw(1,iuvw),1,  frc(1,iatm),1 )
        enddo
        if( iatmlast.gt.0 )then
c         Install layer separation into previous slinky atom ...
          ratmP(3,iatmlast) = ratmP(3,iatmlast) - ratmP(3,iatm)
c         Accumulate projected forces into current slinky atom..
          frcP(3,iatm) = frcP(3,iatm) + frcP(3,iatmlast)
        endif
        iatmlast = iatm
  200 continue
c     And then convert the final (base) atom back to cartesian
      wfrc = frcP(3,iatm)
      do  id=1,3
c       Coordinates had not been shift
        frcP(id,iatmlast) = frc(id,iatmlast) + wfrc*uvw(id,3)
      enddo
c
cdbgDEV
      if( iprint.gt.4 )then
        write(IWR,*) 'SLINK-SELECT summary:'
        write(IWR,'(a,11x,a,20x,a)')'  Atom','Raw force','New force'
        do  islink=isl1,isln
          iatm = islicats(islink)
          write(IWR,'(i5,10x,3f11.6,10x,3f11.6)')  iatm,
     $        (frc(i,iatm),i=1,3),(frcP(i,iatm),i=1,3)
        enddo
        write(IWR,'(a,11x,a,20x,a)')'  Atom','Raw coord','projected'
        do  islink=isl1,isln
          iatm = islicats(islink)
          write(IWR,'(i5,10x,3f11.6,10x,3f11.6)')  iatm,
     $        (ratm(i,iatm),i=1,3),(ratmP(i,iatm),i=1,3)
        enddo
      endif
c  
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SLINKUPD
c
c
      subroutine SLINKUPD( IWR, islic, islics,islicats,vslics,
     $ ratm0, ratmx, ratmP )
c---------------------------------------------------------------
c Purpose: transform back from bond coordinates
c Written:  1-Aug-2017, for v2.68
c Revisions: none
c---------------------------------------------------------------
c
c ratm0 = original cartesian coordinates
c ratmP = original projected coordinates (CM/12)
c ratmx = new project coordinates
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  (zero=0.d0,two=2.d0)
c
      DIMENSION  islics(4,*),islicats(*),vslics(3,*)
      DIMENSION  ratm0(3,*)
      DIMENSION  ratmx(3,*)
      DIMENSION  ratmP(3,*)
c Local declarations
      DIMENSION  uvw(3,3), rx(3)
c Loca diagnostic
      DATA  iprint / 2 /
c
c >>>> EXECUTABLE CODE:
c
      if( iprint.gt.1 ) write(IWR,*)'SLINK-update starts'
c
      if( islics(1,islic).ne.9 ) RETURN
c
      nslink = islics(2,islic)
      if( iprint.gt.2 ) write(IWR,*) 'SLINK number atoms=',nslink
      if( nslink.le.0 ) RETURN
c
      isl1 = islics(3,islic)
      isln = isl1 + nslink - 1
      islicv = islics(4,islic)
      iatm1 = islicats(isl1)
      iatmn = islicats(isln)
c
c     Create coordinate system around slinky axis
      call UVWCOORD( uvw, vslics(1,islicv) )
c
      if( iprint.gt.4 )then
        write(IWR,*) 'SLINK-UPDATE projected summary:'
        write(IWR,9010) 'Atom','Orig  uvw','post-uvw'
 9010   format(1x,a,13x,a,24x,a)
 9011   format(1x,i5,4x,3f11.6,4x,3f11.6)
        do  islink=isl1,isln
          iatm = islicats(islink)
          write(IWR,9011)  iatm,
     $        (ratmP(i,iatm),i=1,3),(ratmx(i,iatm),i=1,3)
        enddo
      endif
c
      wbase = DDOT( 3, uvw(3,1),1, ratmx(1,iatmn),1 )
      do 200 islink=isln-1,isl1,-1
        iatm = islicats(islink)  
        wdelta = ratmx(3,iatm)
        wbase = wbase + wdelta
        ratmx(3,iatm) = wbase
c       Transform back to xyz
        do  ixyz=1,3
          rx(ixyz) = ratmx(1,iatm)*uvw(ixyz,1)
     $             + ratmx(2,iatm)*uvw(ixyz,2)
     $             + ratmx(3,iatm)*uvw(ixyz,3)
        enddo
        call DCOPY( 3, rx,1, ratmx(1,iatm),1 )
  200 continue
c
      if( iprint.gt.4 )then
        write(IWR,9010) 'SLINK-UPDATE summary:'
        write(IWR,9010) 'Atom','Start-xyz','post-xyz'
        do  islink=isl1,isln
          iatm = islicats(islink)
          write(IWR,9011)  iatm,
     $        (ratm0(i,iatm),i=1,3),(ratmx(i,iatm),i=1,3)
        enddo
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> UVWCOORD
c
c
      subroutine UVWCOORD( uvw, axisvec  )
c---------------------------------------------------------------
c Purpose: create new ortho coord system around given axis
c Written:  3-Aug-2017, for v2.68
c Revisions: none
c---------------------------------------------------------------
c
c  The axis will usually be the z-axis (i.e., slab)
c  We make this our w-vec.
c  We sweep to find x,y,z-axis with smallest projection on axis.
c    and then project out w from that cartesian axis to get u.
c  Which them typically would make u=x-axis, and v=y-axis,
c    so that u,v,w maps onto z,y,z in standard cases.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      PARAMETER  ( zero=0.d0, one=1.d0 )
c
      DIMENSION  uvw(3,3), axisvec(3)
c
c >>>> EXECUTABLE CODE
c
c Set our w-coord to be along our axis
      call MKZERO( 9, uvw )
      do  id=1,3
        uvw(id,3) = axisvec(id)
      enddo
c
c Find the x,y,z axis with smallest projection along axis
      iwmin = 1
      wmin = ABS( uvw(iwmin,3) )
      do  id=2,3
        if( ABS( uvw(id,3) ) .lt. wmin )then
          iwmin = id
          wmin = ABS( uvw(iwmin,3) )
        endif
      enddo
c     And project w from that to get u-axis
      uvw(iwmin,1) = one
      udotw = uvw(iwmin,3)
      umag = zero
      do  id=1,3
        uvw(id,1) = uvw(id,1) - udotw*uvw(id,3)
        umag = umag + uvw(id,1)**2
      enddo
      umag = SQRT( umag ) 
      do  i=1,3
        uvw(id,1) = uvw(id,1)/umag
      enddo
c
c And do cross product of w,u to get v-axis
      call CROSS( uvw(1,3),uvw(1,1), uvw(1,2) )
c
c    That's all, Folks!
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
