c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TANGENT
c
c
      subroutine TANGENT( IWRNEB,
     $     natom,riminus,ri,riplus,ti,driminus,driplus,
     $     sumdrisqplus,sumdrisqminus,Lantikink,
     $     energy_lo,energy_image,energy_hi,
     $      dprev, dnext )
c---------------------------------------------------------------
c
c     This subroutine calculates tangents and segment lengths
c     This subroutine calculate displacemnt vectors and
c     unit tangent vectors for each image.
c
c Written: Aidan P. Thompson, Fall 2001
c
c Revision history:
c   7Oct02-PAS/2.54: adapted for SeqQuest 2.54(NEB)
c  11May05-RMV/2.59 add pass of dprev, dnext
c---------------------------------------------------------------
c
      IMPLICIT  NONE
c
c     --- variables passed to and from subroutine
      real*8 ti(3,*),riminus(3,*),ri(3,*),riplus(3,*)
      real*8 driminus(3,*),driplus(3,*)
c
      real*8 sumdrisqplus,sumdrisqminus, dprev,dnext
      real*8 energy_lo,energy_image,energy_hi
      integer  natom
      INTEGER  IWRNEB
c
      LOGICAL  Lantikink
c
c     --- local variables
      real*8 drinorminvplus,drinorminvminus
      real*8 sumtsq,tnorminv,dotprod
      real*8 delumin,delumax,deluminus,deluplus
      integer iatom, i
      real*8  zero,one
      real*8  rsmall
c
      DATA  zero,one / 0.d0,1.d0 /
      DATA  rsmall  / 1.d-4 /
c
c >>>> EXECUTABLE CODE:
c
c Calculate displacement vectors:
c     driplus is defined as the vector from iimg to iimg+1,
c         ie R(iimg+1) - Ri(iimg)
c     driminus is defined as the vector from iimg-1 to iimg,
c         ie R(iimg) - Ri(iimg-1)
c
      sumdrisqplus = zero
      sumdrisqminus = zero
      do  iatom=1,natom
         driminus(1,iatom) = ri(1,iatom) - riminus(1,iatom)
         driminus(2,iatom) = ri(2,iatom) - riminus(2,iatom)
         driminus(3,iatom) = ri(3,iatom) - riminus(3,iatom)
         driplus(1,iatom) = riplus(1,iatom) - ri(1,iatom)
         driplus(2,iatom) = riplus(2,iatom) - ri(2,iatom)
         driplus(3,iatom) = riplus(3,iatom) - ri(3,iatom)
         sumdrisqminus = sumdrisqminus +
     $          driminus(1,iatom)*driminus(1,iatom) +
     $          driminus(2,iatom)*driminus(2,iatom) +
     $          driminus(3,iatom)*driminus(3,iatom)
         sumdrisqplus = sumdrisqplus +
     $          driplus(1,iatom)*driplus(1,iatom) +
     $          driplus(2,iatom)*driplus(2,iatom) +
     $          driplus(3,iatom)*driplus(3,iatom)
      enddo
c
c Check that neighboring images are not identical
c
  913 format(5x,3f15.6)
      if( sumdrisqplus .lt. rsmall )then
         write(IWRNEB,*)'***** ERROR: R(image+1) = R(image) *****'
         write(IWRNEB,*) 'SUM |R(I+1) - R(I)|^2 = ',sumdrisqplus
         do  iatom=1,natom
            write(IWRNEB,*)  'Atom/R(I)/R(I+1) = ',iatom
            write(IWRNEB,913)  ( ri    (i,iatom), i=1,3 )
            write(IWRNEB,913)  ( riplus(i,iatom), i=1,3 )
         enddo
         call STOPXERR( 'NEB TANGENT: R(image+1) = R(image)' )
      endif
      if( sumdrisqminus .lt. rsmall )then
         write(IWRNEB,*)'***** ERROR: R(image-1) = R(image) *****'
         write(IWRNEB,*) 'SUM |R(I-1) - R(I)|^2 = ',sumdrisqminus
         do  iatom=1,natom
            write(IWRNEB,*)  'Atom/R(I)/R(I-1) = ',iatom
            write(IWRNEB,913)  ( ri     (i,iatom), i=1,3 )
            write(IWRNEB,913)  ( riminus(i,iatom), i=1,3 )
         enddo
         call STOPXERR( 'NEB TANGENT: R(image-1) = R(image)' )
      endif
c
      dprev = SQRT( sumdrisqminus )
      dnext = SQRT( sumdrisqplus )
      drinorminvminus = one / SQRT( sumdrisqminus )
      drinorminvplus = one / SQRT( sumdrisqplus )
c
c Check for sharp angles
c
      dotprod = zero
      do  iatom=1,natom
         dotprod = dotprod + driplus(1,iatom)*driminus(1,iatom) +
     $                       driplus(2,iatom)*driminus(2,iatom) +
     $                       driplus(3,iatom)*driminus(3,iatom)
      enddo
c
      if( dotprod .le. zero )then
         write(IWRNEB,*)
         write(IWRNEB,*) '*************************************'
         write(IWRNEB,*) 'WARNING: NEB (I+1):(I):(I-1) Sharp angle'
         write(IWRNEB,*) 'Dotprod = ',dotprod
         do  iatom=1,natom
            write(IWRNEB,*)  'Atom/R(I-1)/R(I)/R(I+1) = ',iatom
            write(IWRNEB,913)  ( riminus(i,iatom), i=1,3 )
            write(IWRNEB,913)  ( ri     (i,iatom), i=1,3 )
            write(IWRNEB,913)  ( riplus (i,iatom), i=1,3 )
         enddo
         write(IWRNEB,*) '*************************************'
         write(IWRNEB,*)
      endif
c
c Calculate the unit tangent vector for each image.
c This is tau in Jonnson (1998)
c
      if( .not. Lantikink )then
c
c        This is the "bisector" method.
         do  iatom=1,natom
            ti(1,iatom) = driplus(1,iatom)*drinorminvplus +
     $                    driminus(1,iatom)*drinorminvminus
            ti(2,iatom) = driplus(2,iatom)*drinorminvplus +
     $                    driminus(2,iatom)*drinorminvminus
            ti(3,iatom) = driplus(3,iatom)*drinorminvplus +
     $                    driminus(3,iatom)*drinorminvminus
         enddo
c
      elseif( Lantikink )then
c
c        Use energy-biased tangent estimate of Henkelman&Jonsson (2000)
c        Vector Tminus refers to R(iimg) - R(iimg-1).
c        Vector Tplus refers to R(iimg+1) - R(iimg).
c
         deluminus = energy_image - energy_lo
         deluplus = energy_hi - energy_image
c
         if( deluminus.ge.zero .and. deluplus.ge.zero )then
c
c           We are on the uphill side, use Tplus as tangent
            do  iatom=1,natom
               ti(1,iatom) = driplus(1,iatom)
               ti(2,iatom) = driplus(2,iatom)
               ti(3,iatom) = driplus(3,iatom)
            enddo
c
         elseif( deluminus.le.zero .and. deluplus.le.zero )then
c
c           We are on the downhill side, use Tminus as tangent
            do  iatom=1,natom
               ti(1,iatom) = driminus(1,iatom)
               ti(2,iatom) = driminus(2,iatom)
               ti(3,iatom) = driminus(3,iatom)
            enddo
c
         else
c
c           We have a local maximum or local minimum,
c           use weighted average of Tminus and Tplus as tangent.
c
            delumax = MAX( ABS(deluminus), ABS(deluplus) )
            delumin = MIN( ABS(deluminus), ABS(deluplus) )
c
            if( energy_hi .ge. energy_lo )then
c              Weight Tplus higher
               do  iatom=1,natom
                  ti(1,iatom) = delumax*driplus(1,iatom) +
     $                          delumin*driminus(1,iatom)
                  ti(2,iatom) = delumax*driplus(2,iatom) +
     $                          delumin*driminus(2,iatom)
                  ti(3,iatom) = delumax*driplus(3,iatom) +
     $                          delumin*driminus(3,iatom)
               enddo
            else
c              Weight Tminus higher
               do  iatom=1,natom
                  ti(1,iatom) = delumin*driplus(1,iatom) +
     $                          delumax*driminus(1,iatom)
                  ti(2,iatom) = delumin*driplus(2,iatom) +
     $                          delumax*driminus(2,iatom)
                  ti(3,iatom) = delumin*driplus(3,iatom) +
     $                          delumax*driminus(3,iatom)
               enddo
            endif
c
         endif
c
c        End Lantikink branch
      endif
c
c Finally, we need to normalize the tangent.
c
      sumtsq = zero
      do  iatom=1,natom
         sumtsq = sumtsq + ti(1,iatom)**2
     $                   + ti(2,iatom)**2
     $                   + ti(3,iatom)**2
      enddo
c
      if( sumtsq .eq. zero )then
         call STOPXERR( 'NEB TANGENT: R(i+1) - R(i-1) = 0' )
      endif
c
      tnorminv = one / SQRT( sumtsq )
c
      do  iatom=1,natom
         ti(1,iatom) = ti(1,iatom)*tnorminv
         ti(2,iatom) = ti(2,iatom)*tnorminv
         ti(3,iatom) = ti(3,iatom)*tnorminv
      enddo
c
c    That's all Folks!
c
      RETURN
      END
