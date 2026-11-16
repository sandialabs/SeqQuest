c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DAT0NEB
c
c
      subroutine DAT0NEB( nimg_neb, nimg_max,
     $  ea_neb,eb_neb, spring_neb,
     $  Lantikink_neb,Loptall_neb, 
     $  Lclimb_neb, Lnoninteract, Lrestart_neb,
     $  imgstat, natm, ratm, ratma )
c---------------------------------------------------------------
c Purpose: initialize parameters for nudged elastic band calc
c
c Written: Peter A. Schultz, 10-Jan-2002, for v2.54
c
c Revision history:
c  11May05-RMV/2.59: new flags for climb, noninteracting images
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      LOGICAL    Lantikink_neb, Loptall_neb, Lrestart_neb
      LOGICAL    Lclimb_neb, Lnoninteract
      DIMENSION  imgstat(*)
      DIMENSION  ratm(3,natm),ratma(3,natm)
c
      DATA       zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Set total number of neb images to one (need as non-0 dim):
      nimg_neb = 1
c
c Total energy of initial and final (endpoint) geometries:
      ea_neb = zero
      eb_neb = zero
c
c Spring constant (artificial) between different images:
      spring_neb = 1.0d0
c
c NEB anti-kink (uphill-only tangents) machinery is: ON
      Lantikink_neb = .true.
c
c NEB optall (min/max images head down-/up-hill) is: OFF
      Loptall_neb = .false.
c
c NEB climbing image is (uphill hill opt only): OFF
      Lclimb_neb = .false.
c
c NEB noninteracting images is: OFF
      Lnoninteract = .false.
c
c NEB restart calculation is: OFF
      Lrestart_neb = .false.
c
c Initialize status to all empty images:
      do  img=1,nimg_max
        imgstat(img) = 0
      enddo
c
c Default the initial image coordinates to the setdata coords:
      do  iatm=1,natm
        ratma(1,iatm) = ratm(1,iatm)
        ratma(2,iatm) = ratm(2,iatm)
        ratma(3,iatm) = ratm(3,iatm)
      enddo
c
c    That's all Folks!
c
      RETURN
      END
