c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEARBY
c
c
      subroutine NEARBY( nearest, natm, itypa, norba,
     $  ratm, rprim,hh,n1r,n2r,n3r,nptr, mngr3,mxgr3,ngr,
     $  r,  rsq )
c---------------------------------------------------------------
c Purpose: find mesh point nearest each atom
c
c Written: Peter A. Schultz, 30-August-1993
c
c Revision history:
c  15Feb08-PAS/2.60f: fix nearby indexing error; skip brute force
c  14Jun05-PAS/2.59: File units
c  21Jun01-PAS/2.48: replace print statements
c  21Dec99-PAS/2.42: comsetic clean-up
c  10Dec98-PAS/2.29: promote source to explicit double precision
c   1May97-PAS/2.17: naming conventions, new grid coord rtn
c  15Dec95-PAS/2.14: shell-by-shell grid ranges
c---------------------------------------------------------------
c
c  Original pjf scheme to find the grid point nearest an atom is
c  brute force search of all grid points within the primary unit
c  cell.  I translate each atom to put it into the primary unit
c  cell to guarantee accuracy (in all but pathological cases).
c
c  In new algorithm to find grid point nearest an atom, begin by
c  finding the coordinates of atom in terms of the grid-defining
c  vectors (solving a 3X3 system of equations), round to integers,
c  perform a very local search, and bring the resulting nearby
c  points (if necessary) back into the primary unit cell by a
c  translation.  This works only with "well-behaved" sets of mesh
c  definition vectors.  Sufficient but not necessary conditions:
c    1. vec3 normal to plane of vec1/vec2 and and projections of
c       vec1 and vec2 onto each other does not exceed the other
c OR 2. all vecs project onto one another < 1/2 and projection of
c       each vector onto the plane defined by the other two falls
c       the primary parallelpiped.
c  These conditions could be satisfied by transforming input grid
c  vectors into a more suitable set, something I do not attempt
c  here.  Instead, I just check simple stringent test, and default
c  to pjf brute force scheme if vectors fail (they usually pass).
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input arrays:
      DIMENSION  itypa(*),norba(*)
      DIMENSION  ratm(3,natm),hh(3,3),rprim(3,3)
      DIMENSION  mngr3(*),mxgr3(*),ngr(*)
c Output array:
      DIMENSION  nearest(3,natm)
c Scratch arrays:
      DIMENSION  r(3,nptr),rsq(nptr)
c
c Local declarations:
      DIMENSION  r0(3),nprim(3),v(3),rat(3),rnr(3),rv(3),gg(3,3)
      DATA  zero,half,one / 0.0d0,0.5d0,1.0d0 /
      DATA  del / 0.49999999d0 /
      DATA  lstout / 0 /
c
c >>>> EXECUTABLE CODE:
c
c First, check if unprepared projection scheme will work
c Actually, we will do "smart" scheme, do basic check, and then
c  revert to brute force (nrmeth=0) only if simple check actually
c  other than just assume smart scheme fails at the outset.
      h11 = hh(1,1)*hh(1,1) + hh(2,1)*hh(2,1) + hh(3,1)*hh(3,1)
      h22 = hh(1,2)*hh(1,2) + hh(2,2)*hh(2,2) + hh(3,2)*hh(3,2)
      h33 = hh(1,3)*hh(1,3) + hh(2,3)*hh(2,3) + hh(3,3)*hh(3,3)
      h12 = hh(1,1)*hh(1,2) + hh(2,1)*hh(2,2) + hh(3,1)*hh(3,2)
      h23 = hh(1,2)*hh(1,3) + hh(2,2)*hh(2,3) + hh(3,2)*hh(3,3)
      h13 = hh(1,3)*hh(1,1) + hh(2,3)*hh(2,1) + hh(3,3)*hh(3,1)
      nrmeth = 1
      if( h13.eq.zero .and. h23.eq.zero )then
        if( ABS(h12).gt.h11 .or. ABS(h12).gt.h22 )then
          nrmeth = 0
        endif
      else
        if( ABS(h12).gt.half*h11 .or. ABS(h12).gt.half*h22 )then
          nrmeth = 0
        endif
c       There are other conditions necessary to guarantee
        nrmeth = 0
      endif
c
c Set up ranges for test to see if we got it right:
      call CROSS( hh(1,1),hh(1,2), gg(1,3) )
      call CROSS( hh(1,2),hh(1,3), gg(1,1) )
      call CROSS( hh(1,3),hh(1,1), gg(1,2) )
      g1mag = SQRT( gg(1,1)**2 + gg(2,1)**2 + gg(3,1)**2 )
      g2mag = SQRT( gg(1,2)**2 + gg(2,2)**2 + gg(3,2)**2 )
      g3mag = SQRT( gg(1,3)**2 + gg(2,3)**2 + gg(3,3)**2 )
      do  ir=1,3
        gg(ir,1) = gg(ir,1)/g1mag
        gg(ir,2) = gg(ir,2)/g2mag
        gg(ir,3) = gg(ir,3)/g3mag
      enddo
      r1eff = ABS( hh(1,1)*gg(1,1) + hh(2,1)*gg(2,1) + hh(3,1)*gg(3,1) )
      r2eff = ABS( hh(1,2)*gg(1,2) + hh(2,2)*gg(2,2) + hh(3,2)*gg(3,2) )
      r3eff = ABS( hh(1,3)*gg(1,3) + hh(2,3)*gg(2,3) + hh(3,3)*gg(3,3) )
c     Set up length along long diagonal from center of box to corner
c      as approximate maximum distance from atom to nearest grid point:
      ratnrmax = SQRT( r1eff**2+r2eff**2+r3eff**2 ) / 1.99d0
c
      nprim(1) = n1r
      nprim(2) = n2r
      nprim(3) = n3r
c
      iorb = 0
      do 1000 iatm=1,natm
c
c       Project atom position onto grid vectors
        call PROJLV( ratm(1,iatm), hh(1,1),hh(1,2),hh(1,3), v )
        do  iprim=1,3
          iv = 0
          if( v(iprim).gt.-del )then
            iv = (v(iprim) + del)
          else
            iv = -ABS(v(iprim)-half)
          endif
          nearest(iprim,iatm) = iv
        enddo
c
        do  i=1,3
          rat(i) = ratm(i,iatm)
        enddo
        do  iprim=1,3
          nv = 0
          if( nearest(iprim,iatm).lt.0 )then
            nv = -nearest(iprim,iatm) + nprim(iprim) - 1
            nv = nv/nprim(iprim)
          elseif( nearest(iprim,iatm).ge.nprim(iprim) )then
            nv = -( nearest(iprim,iatm)/nprim(iprim) )
          endif
          nearest(iprim,iatm) = nearest(iprim,iatm) + nv*nprim(iprim) +1
c         Translate atom coordinates into primary unit cell:
          do  i=1,3
            rat(i) = rat(i) + DBLE(nv)*rprim(i,iprim)
          enddo
        enddo
c
        rsqmin = zero
        do  i=1,3
          rnr(i) = zero
          do  iprim=1,3
            rnr(i) = rnr(i) + DBLE(nearest(iprim,iatm)-1)*hh(i,iprim)
          enddo
          rv(i) = rnr(i) - rat(i)
          rsqmin = rsqmin + rv(i)*rv(i)
        enddo
c       Perform a local search about this point:
        nr1 = nearest(1,iatm)
        nr2 = nearest(2,iatm)
        nr3 = nearest(3,iatm)
        do 163 idif3= -1,1
          do 162 idif2= -1,1
            do 161 idif1= -1,1
              rsqdif = 0
              do  i=1,3
                v(i) = DBLE(idif1)*hh(i,1) + DBLE(idif2)*hh(i,2)
     $               + DBLE(idif3)*hh(i,3)
                rsqdif = rsqdif + (rv(i)+v(i))**2
              enddo
              if( rsqdif.lt.rsqmin )then
                rsqmin = rsqdif
                nearest(1,iatm) = nr1 + idif1
                nearest(2,iatm) = nr2 + idif2
                nearest(3,iatm) = nr3 + idif3
                rnr(1) = rat(1) + rv(1) + v(1)
                rnr(2) = rat(2) + rv(2) + v(2)
                rnr(3) = rat(3) + rv(3) + v(3)
              endif
  161       continue
  162     continue
  163   continue
        if( nearest(1,iatm).ne.nr1 .or. nearest(2,iatm).ne.nr2 .or.
     $      nearest(3,iatm).ne.nr3 )then
c         Check to see that have not inadvertantly left unit cell
          do  iprim=1,3
            if( nearest(iprim,iatm).lt.1 )then
              nearest(iprim,iatm) = nearest(iprim,iatm) + nprim(iprim)
            elseif( nearest(iprim,iatm).gt.nprim(iprim) )then
              nearest(iprim,iatm) = nearest(iprim,iatm) - nprim(iprim)
            endif
          enddo
        endif
c
        if( lstout.gt.0 )then
          call FLGETIWR( IWR )
          write(IWR,*) ' '
          write(IWR,*) 'NEARBY: results'
          write(IWR,*) 'Atom Posn ',iatm,-1,-1,-1,(ratm(i,iatm),i=1,3)
          write(IWR,*) 'Atom shft ',iatm,-2,-2,-2,(rat(i),i=1,3)
          write(IWR,*) 'Atom near ',iatm,(nearest(i,iatm),i=1,3),
     $                  (rnr(i),i=1,3)
        endif
c
        ratnr = (rat(1)-rnr(1))**2 + (rat(2)-rnr(2))**2
     $        + (rat(3)-rnr(3))**2
        ratnr = SQRT( ratnr )
c
        if( ratnr .gt. ratnrmax )then
C        if( nrmeth.eq.0 )then
c
c          PJF BRUTE FORCE SCHEME, MUCH REVISED:
c
          r0(1) = zero
          r0(2) = zero
          r0(3) = zero
          call MKRBOX( hh,r0, n1r,n2r,n3r, r )
c
          imin3 = mngr3(iorb+1)
          imax3 = mxgr3(iorb+1)
          do  iorbat=1,norba(itypa(iatm))
            iorb = iorb + 1
            if( mngr3(iorb).gt.imin3 ) imin3 = mngr3(iorb)
            if( mxgr3(iorb).lt.imax3 ) imax3 = mxgr3(iorb)
          enddo
          nrii = n1r*n2r*(imax3-imin3+1)
          nr0 = n1r*n2r*(imin3-1)
c
c         Compute r**2 at each point within range of atom
          do  nr=nr0+1,nr0+nrii
            rsq(nr) = (r(1,nr)-rat(1))**2 + (r(2,nr)-rat(2))**2 +
     $                (r(3,nr)-rat(3))**2
          enddo
c
          nearpt = n1r*n2r*(nearest(3,iatm)-1) +
     $             n1r*(nearest(2,iatm)-1) + nearest(1,iatm)
          rsqmin = rsq(nearpt)
          do  nr=nr0+1,nr0+nrii
            if( rsq(nr).lt.rsqmin )then
              rsqmin = rsq(nr)
              nearpt = nr
            endif
          enddo
c
          nearm1 = nearpt - 1
          n3a =  nearm1 / (n1r*n2r)
          n2a = (nearm1 - n1r*n2r*n3a) / n1r
          n1a =  nearm1 - n1r*n2r*n3a - n1r*n2a
c
          nearest(1,iatm) = n1a + 1
          nearest(2,iatm) = n2a + 1
          nearest(3,iatm) = n3a + 1
c
          if( lstout.gt.0 )then
           do  i=1,3
             rnr(i) = zero
             do  iprim=1,3
               rnr(i) = rnr(i) + DBLE(nearest(iprim,iatm)-1)*hh(i,iprim)
             enddo
           enddo
C           print*,'Atmnrpjf',iatm,(nearest(i,iatm),i=1,3),(rnr(i),i=1,3)
          endif
c
        endif
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
