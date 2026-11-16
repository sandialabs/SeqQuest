c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> POTNEAR
c
c
      subroutine POTNEAR( nearopt, natm,
     $ ratm,hh,rprim, mxwk,n1r,n2r,n3r,nptr, nearatm,
     $ vslnear,vesnear,excnear,
     $ vslo,vesslo,excslo, wtnr )
c---------------------------------------------------------------
c Purpose: construct values of grid potentials "nearby" atoms
c
c Written: Peter A. Schultz, 27-August-1993
c
c Revision history:
c  10Feb12-PAS/2.63: cut unused args, dignostic write, renamed nearest
c  14Jun05-PAS/2.59: File units
c  21Jul01-PAS/2.49: spin polarized calculations
c  28Jul98-PAS/2.37: altered diagnostics
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  28Apr97-PAS/2.17: naming conventions
c   8Jun94-PAS
c---------------------------------------------------------------
c Notes:
c   Code to do weighted sampling is ugly because need to account
c   for cases where nearby region overlaps more than one uc, in
c   arbitrary directions, and nearatm point is near an atom image
c   rather than atom.  Because of this, gaussian weighting is not
c   expensive.  Except for pathological mesh vectors, scheme should
c   not underflow, while giving a good sampling of nearby region.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c input arrays:
      DIMENSION  nearatm(3,natm), ratm(3,natm),hh(3,3),rprim(3,3)
      DIMENSION  vslo(n1r,n2r,n3r),vesslo(n1r,n2r,n3r)
      DIMENSION  excslo(n1r,n2r,n3r)
c output arrays:
      DIMENSION  vslnear(natm),vesnear(*),excnear(*)
c scratch arrays:
      DIMENSION  wtnr(*)
c local declarations:
      DIMENSION  nprim(3), nrng(3),nnr(3),nruc(3),nrst(3),nrnd(3)
      DIMENSION  rnr(3),rtmp(3),rat(3), gg(3,3)
      DATA  lstout / 0 /
      DATA  zero / 0.d0 /
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
c
      nprim(1) = n1r
      nprim(2) = n2r
      nprim(3) = n3r
      if( lstout.gt.0 )then
        write(IWR,900)
        if( lstout .gt. 2 ) write(IWR,901)
      endif
  900 format(/1x,'POTNEAR: nearby values of grid potentials')
  901 format(2x,'NEAR',' iatm',  8x,'vtot',11x,'ves',12x,'exc')
  902 format(2x,'NEAR',i4,       4(1x,f14.8) )
c
c  Set up ranges for nearby stuff:
c
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
c
c     DEFAULT IS SET ALL NEARBY VALUES TO ZERO:
      do  iatm=1,natm
        vslnear(iatm) = zero
        vesnear(iatm) = zero
        excnear(iatm) = zero
      enddo
c
      if( nearopt.gt.0 )then
c
c       TAKE VALUES AT POINT NEAREST TO EACH ATOM:
        do  iatm=1,natm
          n1a = nearatm(1,iatm)
          n2a = nearatm(2,iatm)
          n3a = nearatm(3,iatm)
          vslnear(iatm) = vslo(n1a,n2a,n3a)
          vesnear(iatm) = vesslo(n1a,n2a,n3a)
          excnear(iatm) = excslo(n1a,n2a,n3a)
        enddo
c
        if( lstout.gt.2 )then
          write(IWR,*) 'values at nearby points:'
          do 30 iatm=1,natm
            write(IWR,902)  iatm, vslnear(iatm),vesnear(iatm),
     $       excnear(iatm)
   30     continue
        endif
c
        if( nearopt.eq.2 )then
c
          if( lstout.gt.0 )
     $     write(IWR,*) 'Using weighted sampling for NEARby values'
c        [NB: code in "grid0xc" and "potnear" must match]
c
c  nearby sampling function parameters:
          alfsmpl = 1.2d0
          rsampl  = 1.75d0
c
c  Determine range in grid point dimensions for nearby sampling:
          rng1 = (rsampl/r1eff + 0.5d0)
          rng2 = (rsampl/r2eff + 0.5d0)
          rng3 = (rsampl/r3eff + 0.5d0)
          nrng(1) = rng1
          nrng(2) = rng2
          nrng(3) = rng3
          nnr(1) = 2*nrng(1) + 1
          nnr(2) = 2*nrng(2) + 1
          nnr(3) = 2*nrng(3) + 1
          nsampl = nnr(1)*nnr(2)*nnr(3)
c
          if( nsampl.gt.mxwk .or. lstout.gt.1 )then
c           Have a problem, let us diagnose it
            write(IWR,'(a,1x,3f8.2,1x,3i5,1x,2i12)') 'NEARATM RNG=',
     $         rng1,rng2,rng3, (nnr(i),i=1,3), nsampl,mxwk
          endif
c
c
          if( nsampl.le.mxwk )then
c
            do 500 iatm=1,natm
              rsqnr = zero
              do  i=1,3
                rat(i) = ratm(i,iatm)
                rnr(i) = zero
                do  iprim=1,3
                  rnr(i) = rnr(i) + (nearatm(iprim,iatm)-1)*hh(i,iprim)
                enddo
                rsqnr = rsqnr + (rnr(i)-rat(i))*(rnr(i)-rat(i))
              enddo
c check that have atom image nearatm "nearatm" point
              do 123 iuc3=-1,1
                do 122 iuc2=-1,1
                  do 121 iuc1=-1,1
                    do  i=1,3
                      rtmp(i) = rnr(i) - ratm(i,iatm) - iuc1*rprim(i,1)
     $                              - iuc2*rprim(i,2) - iuc3*rprim(i,3)
                    enddo
                    rsqnr2 = rtmp(1)**2 + rtmp(2)**2 + rtmp(3)**2
                    if( rsqnr2.lt.rsqnr )then
                      do  i=1,3
                        rat(i) = rnr(i) - rtmp(i)
                      enddo
                      rsqnr = rsqnr2
                    endif
  121             continue
  122           continue
  123         continue
              rn = SQRT(rsqnr)
c
c Check if point is far enough away from atom to maybe be not nearatm
              if( rn.gt.SQRT(r1eff**2+r2eff**2+r3eff**2)/1.99d0 )then
                write(IWR,*) '**** WARNING: POTNEAR weighting suspect'
                write(IWR,*)   'weighted nearby, atom#',iatm
                write(IWR,913) 'r(iatm)   = ',(ratm(i,iatm),i=1,3)
                write(IWR,913) 'r,shifted = ',(rat(i),i=1,3)
                write(IWR,913) 'nearby pt = ',(rnr(i),i=1,3)
  913           format(1x,a,2x,3f12.6)
              endif
c
              do  i=1,3
                rnr(i) = rnr(i) - rat(i)
                do  iprim=1,3
                  rnr(i) = rnr(i) - DBLE(nrng(iprim))*hh(i,iprim)
                enddo
              enddo
              wtsum = zero
              do  ir=1,nsampl
                wtnr(ir) = zero
              enddo
              call MKNRWT( alfsmpl,wtsum,rnr,hh,nnr,nnr(2),nnr(3),wtnr )
              do  iprim=1,3
                nrst(iprim) = nearatm(iprim,iatm) - nrng(iprim)
                if( nrst(iprim).lt.1 )then
                  nv = -nrst(iprim) + nprim(iprim)
                  nv = nv/nprim(iprim)
                  nrst(iprim) = nrst(iprim) + nv*nprim(iprim)
                elseif( nrst(iprim).gt.nprim(iprim) )then
                  nv = (nrst(iprim) + 1)/nprim(iprim)
                  nrst(iprim) = nrst(iprim) - nv*nprim(iprim)
                endif
                nlast = nrst(iprim) + nnr(iprim) - 1
                nruc(iprim) = (nlast-1)/nprim(iprim) + 1
                nrnd(iprim) = nlast - (nruc(iprim)-1)*nprim(iprim)
              enddo
              wts = zero
              vslnr = zero
              vesnr = zero
              excnr = zero
              j3st = 0
              i3st = nrst(3)
              i3nd = n3r
              jtot = 0
              do 333 iuc3=1,nruc(3)
                if( iuc3.eq.nruc(3) ) i3nd = nrnd(3)
                i2st = nrst(2)
                i2nd = n2r
                j2st = j3st
                do 332 iuc2=1,nruc(2)
                  if( iuc2.eq.nruc(2) ) i2nd = nrnd(2)
                  i1st = nrst(1)
                  i1nd = n1r
                  j1st = j2st
                  do 331 iuc1=1,nruc(1)
                    if( iuc1.eq.nruc(1) ) i1nd = nrnd(1)
c
                    jr00 = j1st
                    do 233 i3=i3st,i3nd
                      jr0 = jr00
                      do 232 i2=i2st,i2nd
                        jr = jr0
                        do 231 i1=i1st,i1nd
                          jtot = jtot + 1
                          jr = jr + 1
                          wts = wts + wtnr(jr)
                          vslnr = vslnr + wtnr(jr)*vslo(i1,i2,i3)
                          vesnr = vesnr + wtnr(jr)*vesslo(i1,i2,i3)
                          excnr = excnr + wtnr(jr)*excslo(i1,i2,i3)
  231                   continue
                        jr0 = jr0 + nnr(1)
  232                 continue
                      jr00 = jr00 + nnr(1)*nnr(2)
  233               continue
c
                    j1st = j1st + (i1nd-i1st+1)
                    i1st = 1
  331             continue
                  j2st = j2st + (i2nd-i2st+1)*nnr(1)
                  i2st = 1
  332           continue
                j3st = j3st + (i3nd-i3st+1)*nnr(1)*nnr(2)
                i3st = 1
  333         continue
c
              if( lstout.gt.3 )write(IWR,*)' wtsum,wts=',wtsum,jtot,wts
              vslnear(iatm) = vslnr/wtsum
              vesnear(iatm) = vesnr/wtsum
              excnear(iatm) = excnr/wtsum
c
  500       continue
c           Close nearopt=2, weighting option
          else
            write(IWR,*) '***** ERROR: weighted sampling not done'
            write(IWR,*) '>>>>> nearby option changed from 2 to 1'
            nearopt = 1
          endif
        endif
c       Close non-zero nearby section
      endif
c
      if( lstout.gt.2 .and. nearopt.ne.1 )then
        do  iatm=1,natm
          write(IWR,902)  iatm, vslnear(iatm),vesnear(iatm),
     $     excnear(iatm)
        enddo
      endif
c
c    That's all Folks!
c
      RETURN
      END
c         
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:        
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
