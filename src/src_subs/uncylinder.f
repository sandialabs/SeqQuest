c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> UNCYLINDER
c
c
      subroutine UNCYLINDER(IWR,natm,iatomp,slicopts,ratm0,ratmx,ratmP)
c---------------------------------------------------------------
c Purpose:  project curvelinear coodrdinates back to cartesian
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
      DIMENSION  ratm0(3,natm), ratmx(3,natm),ratmP(3,natm)
      INTEGER    slicopts(natm,3)
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
c >>>> EXECUTABLE CODE:
c
      write(IWR,*)'now entering uncylinder de-projection subroutine'
      Wsq = zero
c      iatomp=iatm
      ia1 = slicopts(iatomp,2)
      write(IWR,*)'first anchor atom is # ',ia1
      ia2 = slicopts(iatomp,3)
      write(IWR,*)'second anchor atom is # ',ia2
      do index = 1,3
         Wvec(index) = ratm0(index,ia2)-ratm0(index,ia1)
         Wsq = wsq + (Wvec(index))**2
      enddo
      Wmag = sqrt(wsq)
      write(IWR,*)'wmag' ,wmag
      do index = 1,3
         Wnorm(index) = Wvec(index)/wmag
      enddo
c write out zprime vector as check
      write(IWR,*)'Wvec = ', wmag,' * ' ,Wnorm(1),Wnorm(2),Wnorm(3)
c-------
c set rvec, based on ratom(iatomp) and new W
c first, get bondvec
      bondsq = zero
      do index = 1,3
         bondvec(index)=ratm0(index,iatomp)-ratm0(index,ia1)
         bondsq = bondsq +(bondvec(index))**2
      enddo
      bondmag = sqrt(bondsq)
      write(IWR,*)'bondmag',bondmag
      write(IWR,*)'bondvec',bondvec(1),bondvec(2),bondvec(3)
c     get Wcoord = bond dot znorm
      Wcoord = zero
      do index = 1,3
         Wcoord = Wcoord + Wnorm(index)*bondvec(index)
      enddo
      write(IWR,*)'Wcoord before step',Wcoord
      write(IWR,*)'check coordinate transfer'
      write(IWR,*)'Wcoord from previous projection',ratmP(3,iatomp)
c      stop
c
      rvecsq = zero
c     now get rvec and Rcoord
      do index = 1,3
         rvec(index) = bondvec(index)-Wcoord*Wnorm(index)
         rvecsq = rvecsq + (rvec(index))**2
      enddo
c     rvecmag = Rcoord = xprime
      Rcoord = sqrt (rvecsq)
      do index = 1,3
         rnorm(index) = rvec(index)/Rcoord
      enddo
      write(IWR,*)'Rcoord  ',Rcoord
      write(IWR,*)'rnorm ',rnorm(1), rnorm(2), rnorm(3)
      theta = zero
c---- 
c     as a check ,also get yprimecoord of iatomproj (should be zero!)
c use cross product  zcrossr
      Vvec(1)=Wcoord*Wnorm(2)*rvec(3)-Wcoord*Wnorm(3)*rvec(2)
      Vvec(2)=Wcoord*Wnorm(3)*rvec(1)-Wcoord*Wnorm(1)*rvec(3)
      Vvec(3)=Wcoord*Wnorm(1)*rvec(2)-Wcoord*Wnorm(2)*rvec(1)
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
c get ycoord = bond dot yprnorm (should be zero!)
      Vcoord = zero
      do index = 1,3
         Vcoord = Vcoord + bondvec(index)*Vnorm(index)
      enddo
      write(IWR,*)'Vcoord =', Vcoord
      theta = asin(Vcoord/rcoord)
      write(IWR,*)"THETA should be zero",theta
c
c check dot products
      VdotW = zero
      VdotU = zero
      WdotU = zero
      do index = 1,3
         VdotW = VdotW + Vvec(index)*Wvec(index)
         VdotU = VdotU + Vvec(index)*rvec(index)
         WdotU = WdotU + Wvec(index)*rvec(index)
      enddo
      write(IWR,*)'VdotW = ',VdotW
      write(IWR,*)'VdotU =', VdotU
      write(IWR,*)'WdotU =', WdotU
c-----
c     update polar coordinates of iatomproj
      write(IWR,*)'wcoord before update',Wcoord
c     
c       theta1 = theta
       Rcoord = ratmx(1,iatomp)
       theta = ratmx(2,iatomp)
       Wcoord = ratmx(3,iatomp)
       write(IWR,*)'Rcoord,theta,wcoord',rcoord,theta,wcoord
c       stop
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
c
c
c---- 
c     Backtransformation of coordinates
c     to cartesian
      write(IWR,*)'polar coordinates, R,theta,w',Rcoord,theta,Wcoord
      cartnew(1) = Rcoord*cos(theta)
      cartnew(2) = Rcoord*sin(theta)
c      cartnew(1) = Rcoord
c      cartnew(2) = theta
      cartnew(3) =Wcoord
      write(IWR,*)'new coordinates',cartnew(1),cartnew(2),cartnew(3)
      write(IWR,*)'original shifted coordinates'
      write(IWR,*)Rcoord,Vcoord,Wcoord
      write(IWR,*)'original unshifted coordinates'
      write(IWR,*)Ratm0(1,iatomp),Ratm0(2,iatomp),Ratm0(3,iatomp)
c rotate to original reference frame
      unrot(1) = cartnew(1)*rnorm(1)  +cartnew(2)*Vnorm(1)  
     &     +cartnew(3)*Wnorm(1)
      unrot(2) = cartnew(1)* rnorm(2) +cartnew(2)*Vnorm(2)  
     &     +cartnew(3)*Wnorm(2)
      unrot(3) = cartnew(1)*rnorm(3)  +cartnew(2)*Vnorm(3)  
     &     +cartnew(3)*Wnorm(3)
c check rotation
      write(IWR,*)'unrotated coordinates'
      write(IWR,*)unrot(1),unrot(2),unrot(3)
c
c reset origin
      do i = 1,3
         unshift(i)=unrot(i)+Ratm0(i,ia1)
c         write(IWR,*)ia1,Ratom(ia1,i)
      enddo
      write(IWR,*)'new coordinates',unshift(1),unshift(2),unshift(3)
c      stop
c
c return cartesian coordinates to ratmx array
      do index = 1,3
         ratmx(index,iatomp) = unshift(index)
      enddo
c
      RETURN
      END
