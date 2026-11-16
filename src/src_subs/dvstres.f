c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DVSTRES
c
c
      subroutine DVSTRES( IWR, ndim, vstress,
     $ n1r,n2r,n3r,nptr,weight, gprim,ws1,ws2,ws3,ftarray,
     $ rhoslo,  espot,rhofft )
c---------------------------------------------------------------
c Purpose: assemble hartree stress-potential for volume term
c
c Written: P.J. Feibelman, 10/12/90
c
c Revision history:
c  27Sep99-PAS/2.40: full three-dimensional stress tensor
c   9Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
c Notes:
c  Only need the diagonal terms.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /CONST/ pi,rtpi,cutexp
c output array:
      DIMENSION  vstress(3,3)
c input arrays:
      DIMENSION  rhoslo(nptr)
      DIMENSION  gprim(3,3)
c scratch arrays:
      DIMENSION  ws1(*),ws2(*),ws3(*), ftarray(*)
      DIMENSION  espot(2,nptr),rhofft(2,nptr)
c local declarations:
      DIMENSION  ave(2)
      DIMENSION  g(3),g2(3),g3(3)
      DATA  extra / 1.d-15 /
      DATA  zero, eight / 0.d0, 8.d0 /
      DATA  lstout / 0 /
      DATA  IWRT / 0 /
c
c >>>> EXECUTABLE CODE:
c
      eightpi = eight*pi
c
      if( lstout.gt.2 ) IWRT = IWR
      if( lstout.gt.1 )
     $ call GRIDAVE(IWR,n1r,n2r,n3r,1,rhoslo,ave,ftarray,'slow density')
c
      call MKZERO( 2*nptr, espot )
      call DCOPY( nptr, rhoslo,1, espot,2 )
c
c  Fourier transform density:
c
      call CFFT3D( espot, n1r,n2r,n3r, ws1,ws2,ws3, ftarray )
c
c  Set rho(g=0), ie. net charge, to zero:
c
      indexg0 = (n1r+2)/2 + n1r*(n2r/2) + n1r*n2r*(n3r/2)
      defect = weight*DBLE(nptr)*espot(1,indexg0)
      if( lstout.gt.2 )
     $ write(IWR,900) indexg0, espot(1,indexg0),espot(2,indexg0)
      if( lstout.gt.0 )
     $ write(IWR,910) 'charge defect (electrons)=',defect
  900 format(/1x,'indexg0=',i8,'; rho(g=0)=',1p2d22.12)
  910 format(1x,a,1p2d22.12)
c
      espot(1,indexg0) = zero
      espot(2,indexg0) = zero
      if( ndim.gt.1 )  call DCOPY( 2*nptr, espot,1, rhofft,1 )
c
      f1 = DBLE((n1r+2)/2)
      f2 = DBLE((n2r+2)/2)
      f3 = DBLE((n3r+2)/2)
c
c  Compute 0.5*Hartree stress potential in g-space.
c
      do 1000 idim=1,ndim
c
        if( idim.gt.1 )  call DCOPY( 2*nptr, rhofft,1, espot,1 )
c
        do  i=1,3
          g3(i) = - f1*gprim(i,1) - f2*gprim(i,2) - f3*gprim(i,3)
        enddo
c
        nr = 0
        do 53 j3=1,n3r
          do  i=1,3
            g3(i) = g3(i) + gprim(i,3)
            g2(i) = g3(i)
          enddo
          do 52 j2=1,n2r
            do  i=1,3
              g2(i) = g2(i) + gprim(i,2)
              g(i)  = g2(i)
            enddo
            do 51 j1=1,n1r
              g(1) = g(1) + gprim(1,1)
              g(2) = g(2) + gprim(2,1)
              g(3) = g(3) + gprim(3,1)
              gsq = g(1)*g(1) + g(2)*g(2) + g(3)*g(3)
              gsqsq = gsq*gsq + extra
              gii = g(idim)*g(idim)
              gfac = eightpi*gii/gsqsq
              nr = nr + 1
              espot(1,nr) = gfac*espot(1,nr)
              espot(2,nr) = gfac*espot(2,nr)
   51       continue
   52     continue
   53   continue
c
c  Fourier transform it back to coordinate space
c
        call CFFT3DB( espot, n1r,n2r,n3r, ws1,ws2,ws3, ftarray )
c
        call GRIDAVE( IWRT, n1r,n2r,n3r,2,espot,ave,ftarray, 'slow Ves')
c
c  Compute contribution to electrostatic stress:
c
c PAS: why is a vacuum potential not extracted here?
        esstres = zero
        do  nr=1,nptr
          esstres = esstres + rhoslo(nr)*espot(1,nr)
        enddo
        vstress(idim,idim) = esstres*weight
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
