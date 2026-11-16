c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MATSTR
c
c
      subroutine MATTOTSTR( matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
c---------------------------------------------------------------
c Purpose: compute how many elements in striped-matrix
c
c Written: Peter A. Schultz,  2-October-2008, for 2.62
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c Output:
      INTEGER  matstr
c Input arrays:
      INTEGER  norb,nk,ncplx, natm,ntyp
      INTEGER  itypa, norba
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
c
c Local declatations:
      INTEGER  jatm, jtyp
      INTEGER  j0, norbatj, njstripe, matjstr, nmatstr
c
c >>>> EXECUTABLE CODE:
c
c Count up how much to clear:
      matstr = 0
      j0 = 0
      do  jatm=1,natm
        jtyp = itypa(jatm)
        norbatj = norba(jtyp)
        njstripe = norb - j0
        matjstr = norbatj*njstripe
        matstr = matstr + matjstr
        j0 = j0 + norbatj
      enddo
      nmatstr = ncplx*nk*matstr
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MATLOCSTR
c
c
      subroutine MATLOCSTR( iproc,nproc, nprocstr,
     $ myatm1, myatm2, myj0, maxmatpr, mymatpr,
     $ matstr,
     $ norb,nk,ncplx, natm,ntyp,itypa, norba )
c---------------------------------------------------------------
c Purpose: apportion stripes of matrix over processors
c
c Written: Peter A. Schultz,  2-October-2008, for 2.62
c
c Revision history:
c  none
c---------------------------------------------------------------
c  This scheme approtions stripes of a block diagoanl matrix
c  over processors; each processor will do a single sequence
c  of atom-blocked stripes, delineated by first (mayatm1) and
c  last (mayatm2) atom of that sequence, inclusive.
c
      IMPLICIT NONE
c
c Output:
      INTEGER  nprocstr
      INTEGER  myatm1,myatm2,myj0
      INTEGER  maxmatpr, mymatpr
c Input arrays:
      INTEGER  iproc,nproc
      INTEGER  matstr
      INTEGER  norb,nk,ncplx, natm,ntyp
      INTEGER  itypa, norba
      DIMENSION  itypa(natm)
      DIMENSION  norba(ntyp)
c
c Local declatations:
      INTEGER  jatm, jtyp
      INTEGER  j0,j0next,j0start, norbatj, njstripe, matjstr
      INTEGER  matavepr
      INTEGER  jproc, matpr, jatm1,jatm2
c
      INTEGER  IWR
      INTEGER  npdbg, npr
      PARAMETER  ( npdbg = 64 )
      INTEGER  matsum
      INTEGER  myats(2,npdbg), mymats(npdbg), myj0s(npdbg)
c
      INTEGER  lvldbg
      DATA     lvldbg / 0 /
c
c >>>> EXECUTABLE CODE:
c
      if( nproc.eq. 1 )then
c       Single-processr does it all ...
        myatm1 = 1
        myatm2 = natm
        myj0 = 0
        maxmatpr = matstr
        mymatpr = matstr
        nprocstr = 1
        goto 900
      endif
c
c Set defaults to do no atoms with no storage
c   (in case I end up with an idle processor)
      myatm1 = natm + 1
      myatm2 = 0
      mymatpr = 0
c
      do  jproc=1,npdbg
        myats(1,jproc) = 0
        myats(2,jproc) = 0
        mymats(jproc) = 0
      enddo
c
c
c Compute the target storage: average the total storage over procs:
      matavepr = ( matstr + nproc - 1 ) / nproc
      matsum = 0
      maxmatpr = 0
      matpr = 0
c
c Count up how much to clear:
      jproc = 0
      j0next = 0
      jatm1 = 1
      j0start = j0next
      do  jatm=1,natm
        jtyp = itypa(jatm)
        norbatj = norba(jtyp)
        j0 = j0next
        j0next = j0 + norbatj
        njstripe = norb - j0
c
        matjstr = norbatj*njstripe
        matpr = matpr + matjstr
        if( matpr.gt.matavepr .or. jatm.eq.natm )then
c         Another processor is full, load up its data
          jatm2 = jatm
          if( jproc.eq.iproc )then
            myatm1 = jatm1
            myatm2 = jatm2
            myj0 = j0start
            mymatpr = matpr
          endif
          if( matpr.gt.maxmatpr ) maxmatpr = matpr
          matsum = matsum + matpr
c
          if( lvldbg.gt.0 .and. (jproc+1).le.npdbg )then
            myats(1,jproc+1) = jatm1
            myats(2,jproc+1) = jatm2
            myj0s(jproc+1) = j0start
            mymats(jproc+1) = matpr
          endif
c
          nprocstr = jproc + 1
c
          jproc = jproc + 1
          jatm1 = jatm + 1
          j0start = j0next
          matpr = 0
c         Finished with this processor allotment
        endif
c
      enddo
c
  900 continue
c
      if( lvldbg.gt.0 )then
        call FLGETIWR( IWR )
        write(IWR,*) 'MATLOC: stripes by procs:'
        npr = MIN( npdbg,nproc )
        do  jproc=1,npr
          write(IWR,'(1x,a,1x,i4,1x,2i5,2x,i6,2x,i12)')
     $     'proc,a1,a2,i0,matpr=',
     $     jproc-1, myats(1,jproc),myats(2,jproc),
     $     myj0s(jproc), mymats(jproc)
        enddo
        write(IWR,*) 'MATLOC: average over processors=',matavepr
        write(IWR,*) 'MATLOC: largest over processors=',maxmatpr
        write(IWR,*) 'MATLOC: size of square matrix   =',norb*norb
        write(IWR,*) 'MATLOC:  ... striped matrix, chk=',matstr,matsum
      endif
c
c    That's all, Folks!
c
      RETURN
      END
