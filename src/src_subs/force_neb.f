c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FORCE_NEB
c
c
      subroutine FORCE_NEB( IWRNEB,
     $     natom,riminus,ri,riplus,fi,
     $     ti,driminus,driplus,kspring,
     $     energy_lo,energy_image,energy_hi,
     $     drsqplus,drsqminus,fpt,fpn,fet,
     $     Lantikink,Loptall,Lclimb,
     $     fdottan,dprev,dnext )
c---------------------------------------------------------------
c Purpose: generate band forces for individual NEB image
c
c Written: Aidan Thompson, Fall 2001
c
c Revision history:
c  14Jun05-PAS/2.59: File units
c  11May05-RMV/2.59: add climb, pass dprev and dnext, pass fdottan
c  14Feb02-PAS/2.54: adapted to SeqQuest - cosmetics
c---------------------------------------------------------------
c
c  fdottan = projection of system force along tangent
c
      IMPLICIT NONE
c
c     --- variables passed to and from subroutine
      real*8 ri(3,*),riminus(3,*),riplus(3,*),fi(3,*)
      real*8 ti(3,*),driminus(3,*),driplus(3,*)
c
      real*8 energy_lo,energy_image,energy_hi
      real*8 drsqplus,drsqminus,fpt,fpn,fet
      real*8 kspring,fdottan,dprev,dnext
      integer natom,IWRNEB
c
      LOGICAL  Lantikink,Loptall,Lclimb
c
c >>>> EXECUTABLE CODE:
c
c Calculate local tangents and segment lengths
c
      call TANGENT( IWRNEB,
     $     natom,riminus,ri,
     $     riplus,ti,
     $     driminus,driplus,
     $     drsqplus,drsqminus,
     $     Lantikink,
     $     energy_lo,energy_image,energy_hi,
     $     dprev,dnext)
c
c Modify computed forces with band forces:
c
      call FORCE_ELASTIC( IWRNEB,
     $     natom,ri,fi,kspring,
     $     fpt,fpn,fet,
     $     ti,driminus,driplus,energy_lo,energy_image,energy_hi,
     $     Loptall,Lclimb,fdottan)
c
c    That's all Folks!
c
      RETURN
      END
