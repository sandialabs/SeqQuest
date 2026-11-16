c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> UNSPHERE
c
c
      subroutine UNSPHERE( IWR,natm,iatomp,slicopts,ratm0,ratmx,ratmP )
c---------------------------------------------------------------
c Purpose:  this routine de-projects forces and coordinates from spherical to cartesian
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
      INTEGER    slicopts(natm,3)
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
c--------
c
c >>>> EXECUTABLE CODE:
c
      write(IWR,*)'entering spherical deprojection routine'
c set reference W direction, based  on anchor atoms
      wsq = zero
      ia1 = slicopts(iatomp,2)
      write(IWR,*)'first anchor atom (sphere center) is # ',ia1
c      stop
      ia2 = slicopts(iatomp,3)
      write(IWR,*)'second anchor atom is # ',ia2
c      stop
c W direction is from a2 to a1
      do index = 1,3
         Wvec(index) = ratm0(index,ia1)-ratm0(index,ia2)
         wsq = wsq + (Wvec(index))**2
      enddo
      wmag = sqrt(wsq)
      write(IWR,*)'wmag (length of reference a1-a2 bond)' ,wmag
c      stop
      do index = 1,3
         Wnorm(index) = Wvec(index)/wmag
      enddo
c write out Wprime vector as check
      write(IWR,*)'Wvec = ', Wmag,' * ' ,Wnorm(1),Wnorm(2),Wnorm(3)
c      stop
c
c-------
c set rvec, based on ratom(iatomp) and anchor1
c first, get bondvec = rvec
c get Rcoord = |bondvec|=bondmag
c NOTE for nested scheme, the bond vector will have to
c be kept from the original coordinate transformation
c rather than re-calculated here
      bondsq = zero
      do index = 1,3
         bondvec(index)=ratm0(index,iatomp)-ratm0(index,ia1)
         bondsq = bondsq +(bondvec(index))**2
      enddo
      bondmag = sqrt(bondsq)
      Rcoord0 = bondmag
      write(IWR,*)'Rcoord0 = ',Rcoord0
      write(IWR,*)'bondvec',bondvec(1),bondvec(2),bondvec(3)
c also useful, the r unit vector
      do index = 1,3
         rnorm(index) = bondvec(index)/Rcoord0
      enddo
c
c      
c     get Wcoord = bond dot Wnorm
      Wcoord0 = zero
      do index = 1,3
         Wcoord0 = Wcoord0 + Wnorm(index)*bondvec(index)
      enddo
c
c     now get xvec =normaliWed component of bondvec perp to Wvec
      Uvecsq = zero
      do index = 1,3
         Uvec(index) = bondvec(index)-Wcoord0*Wnorm(index)
         Uvecsq = Uvecsq + (Uvec(index))**2
      enddo
      Ucoord0 = sqrt(Uvecsq)
c
      do index = 1,3
         Unorm(index) =  Uvec(index)/Ucoord0
      enddo
c
c---- 
c     as a check ,also get yprimecoord of iatomproj (should be zero!)
c use cross product  Wcrossr
      Vvec(1)=Wnorm(2)*Unorm(3)-Wnorm(3)*Unorm(2)
      Vvec(2)=Wnorm(3)*Unorm(1)-Wnorm(1)*Unorm(3)
      Vvec(3)=Wnorm(1)*Unorm(2)-Wnorm(2)*Unorm(1)
      write(IWR,*)'Vvec ',Vvec(1),Vvec(2),Vvec(3)
      Vmagsq = zero
      do index = 1,3
         Vmagsq = Vmagsq +(Vvec(index))**2
      enddo
      Vmag = sqrt(Vmagsq)
      write(IWR,*)'Vmag = ',Vmag
      do index = 1,3
         Vnorm(index)=Vvec(index)/Vmag
      enddo
c get Vcoord = bond dot Vnorm (should be zero!)
      Vcoord0 = zero
      do index = 1,3
         Vcoord0 = Vcoord0 + bondvec(index)*Vnorm(index)
      enddo
      write(IWR,*)'Vcoord0 =', Vcoord0
c
c check dot products
      VdotW = zero
      VdotU = zero
      WdotU = zero
      do index = 1,3
         VdotW = VdotW + Vvec(index)*Wvec(index)
         VdotU = VdotU + Vvec(index)*Uvec(index)
         WdotU = WdotU + Wvec(index)*Uvec(index)
      enddo
      write(IWR,*)'VdotW = ',VdotW
      write(IWR,*)'VdotU =', VdotU
      write(IWR,*)'WdotU =', WdotU
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
      Rcoord=ratmx(1,iatomp)
      phi = ratmx(2,iatomp)
      theta = ratmx(3,iatomp)
      write(IWR,*)'spherical coordinates, R,phi,theta'
      write(IWR,*)Rcoord,phi,theta
      cartnew(1) = Rcoord*sin(phi)*cos(theta)
      cartnew(2) = Rcoord*sin(phi)*sin(theta)
      cartnew(3) =Rcoord*cos(phi)
      write(IWR,*)'new coordinates',cartnew(1),cartnew(2),cartnew(3)
      write(IWR,*)'original shifted coordinates'
      write(IWR,*)Ucoord0,Vcoord0,Wcoord0
      write(IWR,*)'original unshifted coordinates'
      write(IWR,*)Ratm0(1,iatomp),Ratm0(2,iatomp),Ratm0(3,iatomp)
c rotate to original reference frame
      unrot(1) = cartnew(1)*Unorm(1)  +cartnew(2)*Vnorm(1)  
     &     +cartnew(3)*Wnorm(1)
      unrot(2) = cartnew(1)* Unorm(2) +cartnew(2)*Vnorm(2)  
     &     +cartnew(3)*Wnorm(2)
      unrot(3) = cartnew(1)*Unorm(3)  +cartnew(2)*Vnorm(3)  
     &     +cartnew(3)*Wnorm(3)
c check rotation
      write(IWR,*)'PROJ: unrotated coords', unrot(1),unrot(2),unrot(3)
c
c reset origin
      do i = 1,3
         unshift(i)=unrot(i)+Ratm0(i,ia1)
         write(IWR,*)ia1,Ratm0(i,ia1)
      enddo
      write(IWR,*)'PROJ: new coords',unshift(1),unshift(2),unshift(3)
c     stop
c
c return cartesian coordinates to ratmx array
      do index = 1,3
         ratmx(index,iatomp) = unshift(index)
      enddo
c
      RETURN
      END
