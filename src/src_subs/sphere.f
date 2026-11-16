c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SPHERE
c
c
      subroutine SPHERE(IWR,natm,iatomp,slicopts,ratm,frc,ratmP,frcP)
c---------------------------------------------------------------
c Purpose: projects the coordinates and forces to spherical coordinates
c
c Written: Renee M. Van Ginhoven
c
c Revision histoty:
c  28Sep08-PAS/2.62: cosmetic standardization
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
c
cRMV  inputs from gselect
      DIMENSION  ratm(3,natm), frc(3,natm), ratmP(3,natm),frcP(3,natm)
      INTEGER    slicopts(natm,3)
c local
c 
      DIMENSION  Uvec(3),Vvec(3),Wvec(3), rvec(3),bondvec(3)
      DIMENSION  Unorm(3),Vnorm(3),Wnorm(3),rnorm(3)
      DIMENSION  Fuvec(3),Fvvec(3),Fwvec(3),cartnew(3),Fuwperp(3)
      DIMENSION  Fuwperpnorm(3)
c    
      DATA  Pi / 3.141592654 /
      DATA  degperrad / 57.29577951 /
      DATA  zero / 0.d0 /
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
c--------
c
c >>>> EXECUTABLE CODE:
c
c set new z->W direction, based on anchor atoms
      Wsq = zero
      ia1 = slicopts(iatomp,2)
      write(IWR,*)'first anchor atom (sphere center) is # ',ia1
c      stop
      ia2 = slicopts(iatomp,3)
      write(IWR,*)'second anchor atom is # ',ia2
c z direction is from a2 to a1
      do index = 1,3
         Wvec(index) = ratm(index,ia1)-ratm(index,ia2)
         Wsq = Wsq + (Wvec(index))**2
      enddo
      Wmag = sqrt(Wsq)
      write(IWR,*)'Wmag (length of reference a1-a2 bond)' ,Wmag
      do index = 1,3
         Wnorm(index) = Wvec(index)/Wmag
      enddo
c
c-------
c set rvec, based on ratom(iatomp) and anchor1
c first, get bondvec = rvec
c get Rcoord = |bondvec|=bondmag
      bondsq = zero
      do index = 1,3
         bondvec(index)=ratm(index,iatomp)-ratm(index,ia1)
         bondsq = bondsq +(bondvec(index))**2
      enddo
      bondmag = sqrt(bondsq)
      Rcoord = bondmag
c
c also useful, the r unit vector
      do index = 1,3
         Rnorm(index) = bondvec(index)/Rcoord
      enddo
c     
c     get Wcoord = bond dot Wnorm
      Wcoord = zero
      do index = 1,3
         Wcoord = Wcoord + Wnorm(index)*bondvec(index)
      enddo
      write(IWR,*)'Wcoord before step',Wcoord
c
c     now get uvec =normalized component of bondvec perp to Wvec
      Uvecsq = zero
      do index = 1,3
         Uvec(index) = bondvec(index)-Wcoord*Wnorm(index)
         Uvecsq = Uvecsq + (Uvec(index))**2
      enddo
      Ucoord = sqrt(Uvecsq)
c      write(IWR,*)'Ucoord = ',Ucoord
      do index = 1,3
         Unorm(index) =  Uvec(index)/Ucoord
      enddo
c---- 
c     as a check ,also get Vcoord of iatomproj (should be zero!)
c use cross product  Wcrossr
      Vvec(1)=Wnorm(2)*Unorm(3)-Wnorm(3)*Unorm(2)
      Vvec(2)=Wnorm(3)*Unorm(1)-Wnorm(1)*Unorm(3)
      Vvec(3)=Wnorm(1)*Unorm(2)-Wnorm(2)*Unorm(1)
      write(IWR,*)'Vvec ',Vvec(1),Vvec(2),Vvec(3)
      Vprmagsq = zero
      do index = 1,3
         Vprmagsq = Vprmagsq +(Vvec(index))**2
      enddo
      Vprmag = sqrt(Vprmagsq)
      write(IWR,*)'Vprmag = ',Vprmag
      do index = 1,3
         Vnorm(index)=Vvec(index)/Vprmag
      enddo
c get Vcoord = bond dot Vnorm (should be zero!)
      Vcoord = zero
      do index = 1,3
         Vcoord = Vcoord + bondvec(index)*Vnorm(index)
      enddo
      write(IWR,*)'Vcoord(should be zero) =', Vcoord
c and for later we want the "Q" direction, which is just the 
c propjection of R into the UV plane
c  as a check, this initial Q should be the same as U
c
         Qcoordsq = Vcoord**2 + Ucoord**2
         Qcoord = sqrt(Qcoordsq)
         write(IWR,*)'PROJ: Qcoord = ',Qcoord,'  Ucoord = ',Ucoord
c
c-----
c
c  finish spherical projection of atom coordinates
c
      phi = acos(Wcoord/Rcoord)
      write(IWR,*)'phi  in radians from acos(z/r)= ',phi
c
c test that the initial theta is zero
      sinphi = sin(phi)
      write(IWR,*)'sinphi = ',sinphi
      ratio = Vcoord/Rcoord
c     write(IWR,*)'ratio = ',ratio
      theta = asin(ratio/sinphi)
c     write(IWR,*)'theta should be zero!  theta = ',theta
c-------
c     project forces
c     first onto rotated cartesians
      Fu = zero
      Fv = zero
      Fw = zero
      do index = 1,3
         Fu = Fu + frc(index,iatomp)*Unorm(index)
         Fv = Fv + frc(index,iatomp)*Vnorm(index)
         Fw =Fw + frc(index,iatomp)*Wnorm(index)
      enddo
c
c
cFr is the easiest
      Fr = zero
      do index = 1,3
         Fr =Fr + frc(index,iatomp)*rnorm(index)
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
      do index = 1,3
         Fuwperpmagsq = Fuwperpmagsq +Fuwperp(index)**2
      enddo
      Fuwperpmag = sqrt(Fuwperpmagsq)
      write(IWR,*)' Perp to R force',Fuwperpmag
c   
c normalize Fuwperp
      do index=1,3
         Fuwperpnorm(index)= Fuwperp(index)/Fuwperpmag
      enddo
c check math...
c
c find sign of phi force Fuwperp 
c Fuwperp dot Wvec
       FperpdotW = zero
       do index=1,3
          FperpdotW = FperpdotW + Fuwperp(index)*Wvec(index)
       enddo
c set direction of force
       sign = 1.0
       if (FperpdotW.gt.zero)then
          sign = -1.0
       endif
       write(IWR,*)'sign of force',sign
c set Phi force
      Fphi = asin(Fuwperpmag*sign/Rcoord)
      write(IWR,*)'Fphi from asin(Fuwcomp/Rcoord)',Fphi
c      stop
c if phi too close to 180, revert to cartesian
c      if (abs(phi).gt.3.0)then
c         Fphi = Fuwperpmag
c      endif
      sinphi = sin(phi)
      ratio = Fv/(Qcoord)
      write(IWR,*)'PROJ:Fv/Qcoor  ',ratio, Fv, Qcoord
c      ratiocheck = Fv/(Rcoord*sinphi)
      Ftheta = asin(ratio)
c      Fthetacheck = asin(ratiocheck)
      write(IWR,*)'PROJ: Fr,Fphi,Ftheta  ',Fr,Fphi,Ftheta
c,Fthetacheck
c    
c assign spherical forces to force vector
c return projected coordinates and forces in arrays ratmP frcP
      ratmP(1,iatomp) = Rcoord
      ratmP(2,iatomp) = phi
      ratmP(3,iatomp) = theta
c forces
c add scaling factors
      scalephi = 1.5
      scaletheta = 1.5
      frcP(1,iatomp) = Fr
      frcP(2,iatomp) = scalephi*Fphi
      frcP(3,iatomp) = scaletheta*Ftheta
c
      write(IWR,*)'check forces for r,phi,theta',Fr,Fphi,Ftheta
c 
      RETURN
      END
