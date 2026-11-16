c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CYLINDER
c
c
      subroutine CYLINDER(IWR,natm,iatomp,slicopts,ratm,frc,ratmP,frcP)
c---------------------------------------------------------------
c Purpose: projects the coordinates and forces to cylindrical coordinates
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
      INTEGER    slicopts(natm,3)
c
c slicopts contains the projection type and anchor atoms
c slicopts(natm,1) is the projection type
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
c--------
c set new z direction, based  on anchor atoms
      Wsq = zero
c      iatomp=iatm
      ia1 = slicopts(iatomp,2)
      write(IWR,*)'first anchor atom is # ',ia1
      ia2 = slicopts(iatomp,3)
      write(IWR,*)'second anchor atom is # ',ia2
      do index = 1,3
         Wvec(index) = ratm(index,ia2)-ratm(index,ia1)
         Wsq = wsq + (Wvec(index))**2
      enddo
      Wmag = sqrt(wsq)
      write(IWR,*)'wmag' ,Wmag
      do index = 1,3
         Wnorm(index) = Wvec(index)/Wmag
      enddo
c write out W vector as check
      write(IWR,*)'Wvec = ', wmag,' * ' ,Wnorm(1),Wnorm(2),Wnorm(3)
c-------
c set rvec, based on ratom(iatomp) and new w
c first, get bondvec
      bondsq = zero
      do index = 1,3
         bondvec(index)=ratm(index,iatomp)-ratm(index,ia1)
         bondsq = bondsq +(bondvec(index))**2
      enddo
      bondmag = sqrt(bondsq)
      write(IWR,*)'bondmag',bondmag
      write(IWR,*)'bondvec',bondvec(1),bondvec(2),bondvec(3)
c     get Wcoord = bond dot wnorm
      Wcoord = zero
      do index = 1,3
         Wcoord = Wcoord + Wnorm(index)*bondvec(index)
      enddo
      write(IWR,*)'Wcoord before step',Wcoord
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
      Vvec(1)=Wcoord*Wnorm(2)*rvec(3)-Wcoord*wnorm(3)*rvec(2)
      Vvec(2)=Wcoord*Wnorm(3)*rvec(1)-Wcoord*wnorm(1)*rvec(3)
      Vvec(3)=Wcoord*wnorm(1)*rvec(2)-Wcoord*wnorm(2)*rvec(1)
      write(IWR,*)'Vvec ',Vvec(1),Vvec(2),Vvec(3)
      Vmagsq = zero
      do index = 1,3
         Vmagsq = Vmagsq +(Vvec(index))**2
      enddo
      vmag = sqrt(vmagsq)
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
      write(IWR,*)"THETA",theta
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
c-------
c     project forces
c     
      Fw = zero
      do index = 1,3
         Fw =Fw + Frc(index,iatomp)*Wnorm(index)
      enddo
      do index=1,3
         Fwvec(index) =  Fw*wnorm(index) 
      enddo
c      
      Fr = zero
      do index = 1,3
         Fr =Fr + Frc(index,iatomp)*rnorm(index)
      enddo      
      do index=1,3
         Frvec(index) =  Fr*rnorm(index) 
      enddo
c     
c     
      Fv = zero
      Ftheta = zero
      do index = 1,3
        Fv = Fv + Frc(index,iatomp)*Vnorm(index)
      enddo
c 
c
      FthetaCH = asin(Fv/Rcoord)
      write(IWR,*)'check force magnitudes'
      cartfsq = zero
      do i=1,3
         cartfsq=cartfsq+frc(i,iatomp)*frc(i,iatomp)
      enddo
      cartF = sqrt(cartfsq)
      polfsq = zero
      rotcartfsq = Fr*Fr +Fv*Fv+ Fw*Fw
      rotcartf = sqrt(rotcartfsq)
      polfsq = Fr*Fr + FthetaCH*FthetaCH+Fw*Fw
      polf = sqrt(polfsq)
      write(IWR,*)'cart and pol forces',cartf,rotcartf,polf
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
c      Ftheta = asin(ratio)   
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
      Ftheta = 2.0*asin(Fv/Rcoord)
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
c       Ftheta = asin(Fy/rfixed)  
c-----------------------------------------------------
c set forces
c
      write(IWR,*)'check cyl forces',Fr,Ftheta,Fw
      write(IWR,*)'check Fxyz',frc(1,iatomp),frc(2,iatomp),frc(3,iatomp)
c set atom forces to projected forces
      frcP(1,iatomp)=Fr
      frcP(2,iatomp)=Ftheta
      frcP(3,iatomp)=Fw
C set atom coordinates to projected coordinates
      ratmP(1,iatomp)=Rcoord
      ratmP(2,iatomp)=theta
      ratmP(3,iatomp)=Wcoord
c
      RETURN
      END
