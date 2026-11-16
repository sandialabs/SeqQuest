c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FORCE_ELASTIC
c
c
      subroutine FORCE_ELASTIC( IWR,
     $     natom,ri,fi,kspring,fpt,fpn,fet,
     $     ti,driminus,driplus,energy_lo,energy_image,energy_hi,
     $     Loptall,Lclimb,fdottan )
c---------------------------------------------------------------
c     This subroutine takes the potential forces as input, and
c     computes the forces for the nudged elastic band simulation.
c     It does this by estimating the local MEP tangent vector for
c     each image, subtracting the forcefield component parallel to this
c     and adding on the elastic band force component parallel to this.
c
c Written: Aidan P. Thompson, Fall 2001
c
c Revision history:
c  14Jun05-PAS/2.59: File units; 27jan06 cosmetics
c  11May05-RMV/2.59: add climb, and post-processing output
c   7Oct02-PAS/2.54: adapted for Quest 2.54
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
      INTEGER  IWR
c
c     --- variables passed to and from subroutine
      real*8 ri(3,*),fi(3,*),ti(3,*)
      real*8 driminus(3,*),driplus(3,*)
c
      real*8 fpt,fpn,fet,kspring,fdottan
      real*8 energy_lo,energy_image,energy_hi
      integer natom
c
      LOGICAL  Loptall,Lclimb
c
c     --- local variables
      real*8 sumfpnsq,fti1,fti2,fti3,fni1,fni2,fni3
      real*8 deluplus,deluminus
      integer iatom
      real*8 zero
c
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE
c
c Calculate projection of potential force onto MEP
c     This is -delV(Ri).tau in Jonsson (1998)
cRMV original system force along path tangent
cRMV    Fsystem(atoms) dot Tangent(atoms)
cRMV  Tangent is already normalized
c
      fpt = zero
      fdottan= zero
      do  iatom=1,natom
         fpt = fpt + fi(1,iatom) * ti(1,iatom)
         fpt = fpt + fi(2,iatom) * ti(2,iatom)
         fpt = fpt + fi(3,iatom) * ti(3,iatom)
      enddo
cRMV output this force for useful post-processing
      fdottan = fpt
c
c Calculate projection of elastic force on MEP
c     This is Fis.tau in Jonsson (1998)
cRMV  elastic force along path tangent
cRMV    k *  Tangent(atoms) dot  DistancePrevNext(atoms)
c
      fet = zero
      do  iatom=1,natom
         fet = fet + kspring * ti(1,iatom) *
     $         ( driplus(1,iatom) - driminus(1,iatom) )
         fet = fet + kspring * ti(2,iatom) *
     $         ( driplus(2,iatom) - driminus(2,iatom) )
         fet = fet + kspring * ti(3,iatom) *
     $         ( driplus(3,iatom) - driminus(3,iatom) )
      enddo
c
c Calculate total force from normal and tangential components
c
      sumfpnsq = zero
cRMV begin loop over atoms
c
      do  iatom=1,natom
c
c        Calculate potential force normal to MEP
c            This is -delV(Ri)|perp in Jonnson (1998)
cRMV don't actually need to do this,
cRMV  useful for checking path convergence
c
         fni1 = fi(1,iatom) - fpt*ti(1,iatom)
         fni2 = fi(2,iatom) - fpt*ti(2,iatom)
         fni3 = fi(3,iatom) - fpt*ti(3,iatom)
         sumfpnsq = sumfpnsq + fni1**2 + fni2**2 + fni3**2
c
c        Calculate elastic force parallel to MEP
c            This is fis.tau tau in Jonnson (1998)
cRMV springforce * Tangent
c
         fti1 = fet*ti(1,iatom)
         fti2 = fet*ti(2,iatom)
         fti3 = fet*ti(3,iatom)
c
cRMV     climbing image method, to get clos(er)
cRMV  to the nearest saddle point 
cRMV  finds uphill extrema ONLY
c
         if( Lclimb )then  
c
            deluminus = energy_image - energy_lo
            deluplus = energy_hi - energy_image
c     
            if( deluminus .gt. zero .and. deluplus .lt. zero )then
c              Local maximum - climb to saddle point.
cRMV debugging
c               write(IWR,*)'climbing image invoked, this one climbs'
c              Reverse potential force parallel to MEP ...
               fi(1,iatom) = fni1 - fpt*ti(1,iatom)
               fi(2,iatom) = fni2 - fpt*ti(2,iatom)
               fi(3,iatom) = fni3 - fpt*ti(3,iatom)
c               ... and do not apply elastic force
               fet = zero
            else
c              Not a maxima, use normal NEB formula
c               write(IWR,*)'ci invoked, this images does not climb'
               fi(1,iatom) = fti1 + fni1
               fi(2,iatom) = fti2 + fni2
               fi(3,iatom) = fti3 + fni3
            endif
c     
         elseif( Loptall )then
c           Drive images to all local extrema stationary points
c     
            deluminus = energy_image - energy_lo
            deluplus = energy_hi - energy_image
c     
            if( deluminus*deluplus .ge. zero )then
c              Not an extremum - use normal NEB formula.
               fi(1,iatom) = fti1 + fni1
               fi(2,iatom) = fti2 + fni2
               fi(3,iatom) = fti3 + fni3
            elseif( deluminus .gt. zero )then
c              Local maximum - climb to saddle point.
c              Reverse potential force parallel to MEP ...
               fi(1,iatom) = fni1 - fpt*ti(1,iatom)
               fi(2,iatom) = fni2 - fpt*ti(2,iatom)
               fi(3,iatom) = fni3 - fpt*ti(3,iatom)
c               ... and do not apply elastic force
               fet = zero
            else
c              Local minimum - descend to minimum.
c              Use full (downhill) potential force ...
               fi(1,iatom) = fi(1,iatom)
               fi(2,iatom) = fi(2,iatom)
               fi(3,iatom) = fi(3,iatom)
c               ... and do not apply elastic force
               fet = zero
            endif
c     
         else
c           Normal nudged elastic band formula
c     
            fi(1,iatom) = fti1 + fni1
            fi(2,iatom) = fti2 + fni2
            fi(3,iatom) = fti3 + fni3
         endif
c
c        Close loop over atoms
      enddo
c
      fpn = SQRT( sumfpnsq )
c
c    That's all Folks!
c
      RETURN
      END
