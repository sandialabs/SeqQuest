c***********************************************************************
c    BANDS - Quest band structure module
c    Written: Art Edwards, 2008-2012, for version 2.64 (kpp+b.s.)
c    Aided by: Andy Pineda
c
c***********************************************************************
c Revision history:
c  30Jan14-ACP/2.65b: patch to make b.s. symbols to work out right
c  31Jan13-PAS/2.64: installation into production code
c  2008-2012: ACP: merging into main code
c***********************************************************************
c
c    Special points from:
c     Setyawan and Curtarolo, Comp. Mat. Sci. 49,(2010), 299, 312.
c
c    Routines setting up bands for different crystals
c     BANDSBCC
c     BANDSCUB
c     BANDSFCC2
c     BANDSFCC
c     BANDSHEX
c     BANDSMCL
c     BANDSMCL2
c     BANDSORC
c     BANDSRHL
c     BANDSSC
c     BANDSSTC
c     BANDSTC
c     BANDSTET
c     BANDSTRI
c
c    Utility function:
c     INTEGER FUNCTION LOOKUP_SYMPT( sympt,symlabels,nsyms )
c
c***********************************************************************
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSBCC
c
c
c REVISION HISTORY:
c   4Sep12-ACP: Forked from BANDSHEX
c   3May12-ACP: Merged into 2.64
c  20Dec10-ACP: Fully merged into 2.62dev.
c  07Oct08-AHE: Written for 2.61b
c
      SUBROUTINE BANDSBCC(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (numpts=4)
c
c     BCC Bravais Lattice
c     One of a suite of routines that calculate a set of branches for
c     band structure calculations.
c
c     Conventional Lattice Setting
c      a1=(a,0,0)
c      a2=(0,a,0)
c      a3=(0,0,a)
c     Primitive Lattice Setting
c      a1=(-a/2,a/2,a/2)
c      a2=(a/2,-a/2,a/2)
c      a3=(a/2,a/2,-a/2)
c     
c     First attempt to parse the bsymb routine into constituent
c     characters, and check against the set of special points in BCC.
c
      CHARACTER*10 thisroutine
      CHARACTER*4  latname
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER symlabels(numpts),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),gprim(3,3)
      DIMENSION gspec(3,numpts),dx(3),dk(nknscfd),nbr(nbranchd)
      DIMENSION symcoeffs(3,numpts)
      DATA thisroutine /'BANDSBCC'/
      DATA latname     /'BCC'/
      DATA symlabels/'G','H','P','N'/
      DATA symcoeffs / 0.0d0,  0.0d0,  0.0d0,
     $                 0.5d0, -0.5d0,  0.5d0,
     $                 0.25d0, 0.25d0, 0.25d0,
     $                 0.0d0,  0.0d0,  0.5d0 /
      np1=30
c
c     Use the primitive lattice vectors to make special points
c
      WRITE(IWR,*)thisroutine//': primitive reciprical lattice vectors'
      DO ivec=1,3
         WRITE(IWR,*)(gprim(icart,ivec),icart=1,3)
      ENDDO
      write(IWR,*)
      WRITE(IWR,*)thisroutine//': symcoeffs'
      DO ispec=1,numpts
        WRITE(IWR,*)symlabels(ispec),':',
     $    (symcoeffs(ivec,ispec),ivec=1,3)
      END DO
      WRITE(IWR,*)
c     Special points from Setyawan and Curtarolo, Comp. Mat. Sci. 49,
c     (2010), 299, 312.
      DO ispec=1,numpts
        DO icart = 1,3
          gspec(icart,ispec)=0.0d0
          DO ivec = 1,3
            gspec(icart,ispec) = gspec(icart,ispec)
     $        + gprim(icart,ivec)*symcoeffs(ivec,ispec)
          END DO
        END DO
      END DO
c
c     Now write out the values of the special points
c
      DO ispec=1,numpts
         WRITE(IWR,*) '  ',symlabels(ispec),' ',
     $    (gspec(icart,ispec),icart=1,3)
         WRITE(IWR,*)' '
      ENDDO      
c Use default path if nbranch<=0
      IF (nbranch.LE.0) THEN
        WRITE(IWR,*)'Default '//latname//
     $    ' bandstructure option selected'
        bsymb(1)='GH'
        bsymb(2)='HN'
        bsymb(3)='NG'
        bsymb(4)='GP'
        bsymb(5)='PH'
        bsymb(6)='PN'
        nbranch=6
      ENDIF
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         ibegin=LOOKUP_SYMPT(symb, symlabels, numpts)
         IF (ibegin.LT.1) THEN
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
         iend=LOOKUP_SYMPT(syme, symlabels, numpts)
         IF (iend.LT.1) THEN
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ. LOOKUP_SYMPT('G', symlabels, numpts)) THEN 
            IF(iend.EQ. LOOKUP_SYMPT('H', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('P', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('N', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('H', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('P', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('N', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('P', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('H', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('N', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('N', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('H', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('P', symlabels, numpts )) GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)thisroutine//': d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/dfloat(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)thisroutine//': nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSCUB
c
c
c REVISION HISTORY:
c   4Sep12-ACP: Forked from BANDSHEX
c   3May12-ACP: Merged into 2.64
c  20Dec10-ACP: Fully merged into 2.62dev.
c  07Oct08-AHE: Written for 2.61b
c
      SUBROUTINE BANDSCUB(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (numpts=4)
c
c     CUB Bravais Lattice
c     One of a suite of routines that calculate a set of branches for
c     band structure calculations.
c     
c     Conventional Lattice Setting
c      a1=(a,0,0)
c      a2=(0,a,0)
c      a3=(0,0,a)
c
c     First attempt to parse the bsymb routine into constituent
c     characters, and check against the set of special points in CUB.
c
      CHARACTER*10 thisroutine
      CHARACTER*4  latname
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER symlabels(numpts),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),gprim(3,3)
      DIMENSION gspec(3,numpts),dx(3),dk(nknscfd),nbr(nbranchd)
      DIMENSION symcoeffs(3,numpts)
      DATA thisroutine /'BANDSCUB'/
      DATA latname     /'CUB'/
      DATA symlabels/'G','M','R','X'/
      DATA symcoeffs / 0.0d0, 0.0d0, 0.0d0,
     $                 0.5d0, 0.5d0, 0.0d0,
     $                 0.5d0, 0.5d0, 0.5d0,
     $                 0.0d0, 0.5d0, 0.0d0 /
      np1=30
c
c     Use the primitive lattice vectors to make special points
c
      WRITE(IWR,*)thisroutine//': primitive reciprical lattice vectors'
      DO ivec=1,3
         WRITE(IWR,*)(gprim(icart,ivec),icart=1,3)
      ENDDO
      WRITE(IWR,*)
      WRITE(IWR,*)thisroutine//': symcoeffs'
      DO ispec=1,numpts
        WRITE(IWR,*)symlabels(ispec),':',
     $    (symcoeffs(ivec,ispec),ivec=1,3)
      END DO
      WRITE(IWR,*)
c     Special points from Setyawan and Curtarolo, Comp. Mat. Sci. 49,
c     (2010), 299, 312.
      DO ispec=1,numpts
        DO icart = 1,3
          gspec(icart,ispec)=0.0d0
          DO ivec = 1,3
            gspec(icart,ispec) = gspec(icart,ispec)
     $        + gprim(icart,ivec)*symcoeffs(ivec,ispec)
          END DO
        END DO
      END DO
c
c     Now write out the values of the special points
c
      DO ispec=1,numpts
         WRITE(IWR,*) '  ',symlabels(ispec),' ',
     $    (gspec(icart,ispec),icart=1,3)
         WRITE(IWR,*)' '
      ENDDO      
c Use default path if nbranch<=0
      IF (nbranch.LE.0) THEN
        WRITE(IWR,*)'Default '//latname//
     $    ' bandstructure option selected'
        bsymb(1)='GX'
        bsymb(2)='XM'
        bsymb(3)='MG'
        bsymb(4)='GR'
        bsymb(5)='RX'
        bsymb(6)='MR'
        nbranch=6
      ENDIF
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         ibegin=LOOKUP_SYMPT(symb, symlabels, numpts)
         IF (ibegin.LT.1) THEN
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
         iend=LOOKUP_SYMPT(syme, symlabels, numpts)
         IF (iend.LT.1) THEN
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ. LOOKUP_SYMPT('G', symlabels, numpts)) THEN 
            IF(iend.EQ. LOOKUP_SYMPT('R', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('M', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('R', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('R', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('X', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('R', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)thisroutine//': d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/dfloat(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)thisroutine//': nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSFCC2
c
c
c REVISION HISTORY:
c   4Sep12-ACP: Forked from BANDSHEX
c   3May12-ACP: Merged into 2.64
c  20Dec10-ACP: Fully merged into 2.62dev.
c  07Oct08-AHE: Written for 2.61b
c
      SUBROUTINE BANDSFCC2(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (numpts=6)
c
c     FCC Bravais Lattice
c     One of a suite of routines that calculate a set of branches for
c     band structure calculations.
c     Conventional Lattice Setting
c      a1=(a,0,0)
c      a2=(0,a,0)
c      a3=(0,0,a)
c     Primitive Lattice Setting
c      a1=(0,a/2,a/2)
c      a2=(a/2,0,a/2)
c      a3=(a/2,a/2,0)
c     First attempt to parse the bsymb routine into constituent
c     characters, and check against the set of special points in FCC.
c
      CHARACTER*10 thisroutine
      CHARACTER*4  latname
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER symlabels(numpts),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),gprim(3,3)
      DIMENSION gspec(3,numpts),dx(3),dk(nknscfd),nbr(nbranchd)
      DIMENSION symcoeffs(3,numpts)
      DATA thisroutine /'BANDSFCC2'/
      DATA latname     /'FCC'/
      DATA symlabels/'G','K','L','U','W','X'/
      DATA symcoeffs / 0.000d0, 0.000d0, 0.000d0,
     $                 0.375d0, 0.375d0, 0.750d0,
     $                 0.500d0, 0.500d0, 0.500d0,
     $                 0.625d0, 0.250d0, 0.625d0,
     $                 0.500d0, 0.250d0, 0.750d0,
     $                 0.500d0, 0.000d0, 0.500d0 /
      np1=30
c
c     Use the primitive lattice vectors to make special points
c
      WRITE(IWR,*)thisroutine//': primitive reciprical lattice vectors'
      DO ivec=1,3
         WRITE(IWR,*)(gprim(icart,ivec),icart=1,3)
      ENDDO
      WRITE(IWR,*)
      WRITE(IWR,*)thisroutine//': symcoeffs'
      DO ispec=1,numpts
        WRITE(IWR,*)symlabels(ispec),':',
     $    (symcoeffs(ivec,ispec),ivec=1,3)
      END DO
      WRITE(IWR,*)
c     Special points from Setyawan and Curtarolo, Comp. Mat. Sci. 49,
c     (2010), 299, 312.
      DO ispec=1,numpts
        DO icart = 1,3
          gspec(icart,ispec)=0.0d0
          DO ivec = 1,3
            gspec(icart,ispec) = gspec(icart,ispec)
     $        + gprim(icart,ivec)*symcoeffs(ivec,ispec)
          END DO
        END DO
      END DO
c
c     Now write out the values of the special points
c
      DO ispec=1,numpts
         WRITE(IWR,*) '  ',symlabels(ispec),' ',
     $    (gspec(icart,ispec),icart=1,3)
         WRITE(IWR,*)' '
      ENDDO      
c Use default path if nbranch<=0
      IF (nbranch.LE.0) THEN
        WRITE(IWR,*)'Default '//latname//
     $    ' bandstructure option selected'
        bsymb(1)='GX'
        bsymb(2)='XW'
        bsymb(3)='WK'
        bsymb(4)='KG'
        bsymb(5)='GL'
        bsymb(6)='LU'
        bsymb(7)='UW'
        bsymb(8)='WL'
        bsymb(9)='LK'
        bsymb(10)='UX'
        nbranch=10
      ENDIF
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         ibegin=LOOKUP_SYMPT(symb, symlabels, numpts)
         IF (ibegin.LT.1) THEN
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
         iend=LOOKUP_SYMPT(syme, symlabels, numpts)
         IF (iend.LT.1) THEN
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ. LOOKUP_SYMPT('G', symlabels, numpts)) THEN 
            IF(iend.EQ. LOOKUP_SYMPT('L', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('K', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('K', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('L', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('W', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('L', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('K', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('W', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('U', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('U', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('L', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('W', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('W', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('L', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('U', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('K', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('X', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('U', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('W', symlabels, numpts )) GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)thisroutine//': d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/dfloat(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)thisroutine//': nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSFCC
c
c
c REVISION HISTORY:
c   3May12-ACP: Merged into 2.64
c  20Dec10-ACP: Fully merged into 2.62dev.
c  07Oct08-AHE: Written for 2.61b.
c
      SUBROUTINE BANDSFCC(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
cAHE 100708
c
c     One of a suite of routines that calculate a set of branches for band structure 
c     calculations.
c     
c     First attempt to parse the bsymb routine into constituent characters, and check
c     against the set of special points in FCC
c
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER FCC(7),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),gprim(3,3)
      DIMENSION gspec(3,7),dx(3),dk(nknscfd),nbr(nbranchd)
      DIMENSION xnorm(3),wp(3)
      DATA FCC/'G','X','L','K','W','U','A'/
      np1=30
c
c     Use the primitive lattice vectors to make special points
c
      DO i = 1,3
c     G (gamma)
         gspec(i,1)=0.0d0
cc     K
c         gspec(i,4)=(gprim(i,1)+gprim(i,2)+2.0*gprim(i,3))*0.5
c     X
         gspec(i,2)=(gprim(i,2)+gprim(i,3))*0.5
c     L
         gspec(i,3)=(gprim(i,1)+gprim(i,2)+gprim(i,3))*0.5
c     A another gamma point
         gspec(i,7)=gprim(i,1)+gprim(i,2)+gprim(i,3)
      ENDDO
c
c     These are the edge and corner points that are derived from other special points
c     We have to test whether the new points are related correctly to the old points
c
c     K
      IF(gspec(3,2).GT.1.0d-5) THEN
         gspec(1,4)=1.5*gspec(1,3)
         gspec(3,4)=gspec(1,4)
         gspec(2,4)=0.0d0
      ELSE
         gspec(1,4)=1.5*gspec(1,3)
         gspec(2,4)=gspec(1,4)
         gspec(3,4)=0.0d0
      ENDIF
c
c     W and U
c
      ddX=0.0d0
      DO icart=1,3
         ddX=ddX+gspec(icart,2)*gspec(icart,2)
      ENDDO
      ddX=SQRT(ddX)
      WRITE(IWR,*)'ddX = ', ddx
      DO icart = 1,3
         IF(gspec(icart,2).GT.1.0d-5) THEN 
            gspec(icart,5)=gspec(icart,2)
            gspec(icart,6)=gspec(icart,2)
         ELSEIF(gspec(icart,4).LT.1.0d-5) THEN
            gspec(icart,5)=0.0
            gspec(icart,6)=0.5*(3*gspec(icart,3)-ddX)
         ELSE
            gspec(icart,5)=3*gspec(icart,3)-ddX
            gspec(icart,6)=0.5*(3*gspec(icart,3)-ddX)
         ENDIF
      ENDDO
c
c     Now write out all of the special points
c
c
c     Now write out the values of the special points
c
      DO ispec=1,7
         WRITE(IWR,*) '  ',FCC(ispec),' ',(gspec(icart,ispec),icart=1,3)
         WRITE(IWR,*)' '
      ENDDO
c Use default path if nbranch<=0
      IF (nbranch.LE.0) THEN
        WRITE(IWR,*)'Default FCC'//
     $    ' bandstructure option selected'
        bsymb(1)='GX'
        bsymb(2)='XW'
        bsymb(3)='WK'
        bsymb(4)='KG'
        bsymb(5)='GL'
        bsymb(6)='LU'
        bsymb(7)='UW'
        bsymb(8)='WL'
        bsymb(9)='LK'
        bsymb(10)='UX'
        nbranch=10
      ENDIF
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         write(IWR,*)symb,'  ',syme
         IF(symb.EQ.'G') THEN
            ibegin=1
         ELSEIF(symb.EQ.'X') THEN
            ibegin=2
         ELSEIF(symb.EQ.'L')THEN
            ibegin=3
         ELSEIF(symb.EQ.'K')THEN
            ibegin=4
         ELSEIF(symb.EQ.'W') THEN
            ibegin=5
         ELSEIF(symb.EQ.'U')THEN
            ibegin=6
         ELSEIF(symb.EQ.'A')THEN
            ibegin=7
         ELSE
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'not in the FCC set of special points.  Bye!'
            STOP
         ENDIF
         IF(syme.EQ.'G') THEN
            iend=1
         ELSEIF(syme.EQ.'X') THEN
            iend=2
         ELSEIF(syme.EQ.'L')THEN
            iend=3
         ELSEIF(syme.EQ.'K')THEN
            iend=4
         ELSEIF(syme.EQ.'W')THEN
            iend=5
         ELSEIF(syme.EQ.'U')THEN
            iend=6
         ELSEIF(syme.EQ.'A')THEN
            iend=7
         ELSE
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'not in the FCC set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ.1.AND. iend.NE.1) THEN 
            GO TO 10
         ELSEIF(ibegin.EQ.2) THEN
            IF(iend.NE.4.AND. iend.NE.3) GO TO 10
         ELSEIF(ibegin.EQ.3) THEN
            IF(iend.EQ.4.OR.iend.EQ.1.OR.iend.EQ.5.OR.iend.EQ.6)  
     $           GO TO 10
         ELSEIF(ibegin.EQ.4) THEN
            if(iend.EQ.1.OR.iend.EQ.3) GO TO 10
         ELSEIF(ibegin.EQ.5.AND.iend.NE.5) THEN
            GO TO 10
         ELSEIF(ibegin.EQ.6)THEN
            if(iend.NE.6.AND.iend.NE.4)GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)'In BANDSFCC, d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/DFLOAT(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)'in BANDSFCC, nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSHEX
c
c
c REVISION HISTORY:
c  21Feb13-ACP:  Forked from BANDSFCC2 (to fix bug)
c   3May12-ACP: Merged into 2.64
c  20Dec10-ACP: Fully merged into 2.62dev.
c  07Oct08-AHE: Written for 2.61b
c
      SUBROUTINE BANDSHEX(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (numpts=6)
c
c     HEX Bravais Lattice
c     One of a suite of routines that calculate a set of branches for
c     band structure calculations.
c     Lattice Setting
c      a1=(a/2, -a*sqrt(3)/2, 0)
c      a2=(a/2, a*sqrt(3)/2,  0)
c      a3=(0,   0,            c)
c     First attempt to parse the bsymb routine into constituent
c     characters, and check against the set of special points in FCC.
c
      CHARACTER*10 thisroutine
      CHARACTER*4  latname
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER symlabels(numpts),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),gprim(3,3)
      DIMENSION gspec(3,numpts),dx(3),dk(nknscfd),nbr(nbranchd)
      DIMENSION symcoeffs(3,numpts)
      DATA thisroutine /'BANDSHEX'/
      DATA latname     /'HEX'/
      DATA symlabels/'G','A','H','K','L','M'/
      DATA symcoeffs / 0.000d0, 0.000d0, 0.000d0,
     $                 0.000d0, 0.000d0, 0.500d0,
     $          0.333333333333333333d0, 0.333333333333333333d0, 0.500d0,
     $          0.333333333333333333d0, 0.333333333333333333d0, 0.000d0,
     $                 0.500d0, 0.000d0, 0.500d0,
     $                 0.500d0, 0.000d0, 0.000d0 /
      np1=30
c
c     Use the primitive lattice vectors to make special points
c
      WRITE(IWR,*)thisroutine//': primitive reciprical lattice vectors'
      DO ivec=1,3
         WRITE(IWR,*)(gprim(icart,ivec),icart=1,3)
      ENDDO
      WRITE(IWR,*)
      WRITE(IWR,*)thisroutine//': symcoeffs'
      DO ispec=1,numpts
        WRITE(IWR,*)symlabels(ispec),':',
     $    (symcoeffs(ivec,ispec),ivec=1,3)
      END DO
      WRITE(IWR,*)
c     Special points from Setyawan and Curtarolo, Comp. Mat. Sci. 49,
c     (2010), 299, 312.
      DO ispec=1,numpts
        DO icart = 1,3
          gspec(icart,ispec)=0.0d0
          DO ivec = 1,3
            gspec(icart,ispec) = gspec(icart,ispec)
     $        + gprim(icart,ivec)*symcoeffs(ivec,ispec)
          END DO
        END DO
      END DO
c
c     Now write out the values of the special points
c
      DO ispec=1,numpts
         WRITE(IWR,*) '  ',symlabels(ispec),' ',
     $    (gspec(icart,ispec),icart=1,3)
         WRITE(IWR,*)' '
      ENDDO      
c Use default path if nbranch<=0
      IF (nbranch.LE.0) THEN
        WRITE(IWR,*)'Default '//latname//
     $    ' bandstructure option selected'
        bsymb(1)='GM'
        bsymb(2)='MK'
        bsymb(3)='KG'
        bsymb(4)='GA'
        bsymb(5)='AL'
        bsymb(6)='LH'
        bsymb(7)='HA'
        bsymb(8)='LM'
        bsymb(9)='KH'
        nbranch=9
      ENDIF
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         ibegin=LOOKUP_SYMPT(symb, symlabels, numpts)
         IF (ibegin.LT.1) THEN
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
         iend=LOOKUP_SYMPT(syme, symlabels, numpts)
         IF (iend.LT.1) THEN
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ. LOOKUP_SYMPT('G', symlabels, numpts)) THEN 
            IF(iend.EQ. LOOKUP_SYMPT('A', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('K', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('A', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('L', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('H', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('L', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('A', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('H', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('M', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('K', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('L', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('H', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('A', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('L', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('K', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('K', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('H', symlabels, numpts )) GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)thisroutine//': d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/dfloat(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)thisroutine//': nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSMCL
c
c
c REVISION HISTORY:
c   4Sep12-ACP: Forked from BANDSHEX
c   3May12-ACP: Merged into 2.64
c  20Dec10-ACP: Fully merged into 2.62dev.
c  07Oct08-AHE: Written for 2.61b
c
      SUBROUTINE BANDSMCL(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, rprim, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (numpts=16)
c
c     Simple Monoclinic Bravais Lattice
c      Ordering of the lattice a,b<=c; alpha< 90 deg; beta=gamma=90 deg
c     Lattice Setting
c      a1=(a, 0, 0)
c      a2=(0, b, 0)
c      a3=(0, c*cos(alpha), c*sin(alpha))
c     
      CHARACTER*10 thisroutine
      CHARACTER*4  latname
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER symlabels(numpts),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),rprim(3,3),gprim(3,3)
      DIMENSION gspec(3,numpts),dx(3),dk(nknscfd),nbr(nbranchd)
      DIMENSION symcoeffs(3,numpts)
      DATA thisroutine /'BANDSMCL'/
      DATA latname     /'MCL'/
c
c     symlabels:  GAMMA-A-C-D-D1-E-H-H1-H2-M-M1-M2-X-Y-Y1-Z
c
      DATA symlabels/'G','A','C',
     $               'D','d','E',
     $               'H','h','i',
     $               'M','m','n',
     $               'X','Y','y','Z'/
c
      DATA symcoeffs / 0.0d0,  0.0d0,  0.0d0,
     $                 0.5d0,  0.5d0,  0.0d0,
     $                 0.0d0,  0.5d0,  0.5d0,
     $                 0.5d0,  0.0d0,  0.5d0,
     $                 0.5d0,  0.0d0, -0.5d0,
     $                 0.5d0,  0.5d0,  0.5d0,
     $                 0.0d0,  100.0d0,  100.0d0,
     $                 0.0d0,  100.0d0,  100.0d0,
     $                 0.0d0,  100.0d0,  100.0d0,
     $                 0.5d0,  100.0d0,  100.0d0,
     $                 0.5d0,  100.0d0,  100.0d0,
     $                 0.5d0,  100.0d0,  100.0d0,
     $                 0.0d0,  0.5d0,  0.0d0,
     $                 0.0d0,  0.0d0,  0.5d0,
     $                 0.0d0,  0.0d0, -0.5d0,
     $                 0.5d0,  0.0d0,  0.0d0 /
c
      np1=30
c
c Compute special points that depend on rprim.
c
      a=DDOT(3,rprim(1,1),1,rprim(1,1),1)
      a=SQRT(a)
      b=DDOT(3,rprim(1,2),1,rprim(1,2),1)
      b=SQRT(b)
      c=DDOT(3,rprim(1,3),1,rprim(1,3),1)
      c=SQRT(c)
c
      cosalpha=rprim(2,3)/c
      IF (cosalpha.LE.0.0d0) THEN
        CALL BANDSMCL2(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $    nknscfd, rprim, gprim,dk,nbr,IWR)
        RETURN
      ENDIF
c
      sinalpha=rprim(3,3)/c
      IF (sinalpha.LT.0.0d0) THEN
        STOP 'BANDSMCL: sin(alpha) < 0.0d0'
      ENDIF

      should_be_zero=MAX(ABS(rprim(2,1)),ABS(rprim(3,1)),
     $                   ABS(rprim(1,2)),ABS(rprim(3,3)),
     $                   ABS(rprim(1,3)),
     $                   ABS(1.0d0-cosalpha**2-sinalpha**2))
c
      IF(should_be_zero.GT.1.d-8) THEN
        STOP 'BANDSMCL: error in rprim'
      ENDIF
c
      eta=0.5d0*(1.0d0-b*cosalpha/c)/sinalpha/sinalpha
      gnu=0.5d0-eta*c*cosalpha/b
c  H
      symcoeffs(1,7)=0.0d0
      symcoeffs(2,7)=eta
      symcoeffs(3,7)=1.0d0-gnu
c  H1=h
      symcoeffs(1,8)=0.0d0
      symcoeffs(2,8)=1.0d0-eta
      symcoeffs(3,8)=gnu
c  H2=i
      symcoeffs(1,9)=0.0d0
      symcoeffs(2,9)=eta
      symcoeffs(3,9)=-gnu
c  M
      symcoeffs(1,10)=0.5d0
      symcoeffs(2,10)=eta
      symcoeffs(3,10)=1.0d0-gnu
c  M1=m
      symcoeffs(1,11)=0.5d0
      symcoeffs(2,11)=1.0-eta
      symcoeffs(3,11)=gnu
c  M2=n
      symcoeffs(1,12)=0.5d0
      symcoeffs(2,12)=eta
      symcoeffs(3,12)=-gnu
c
c     Use the primitive lattice vectors to make special points
c
      WRITE(IWR,*)thisroutine//': primitive reciprical lattice vectors'
      DO ivec=1,3
         WRITE(IWR,*)(gprim(icart,ivec),icart=1,3)
      ENDDO
c
      WRITE(IWR,*)
      WRITE(IWR,*)thisroutine//': symcoeffs'
      DO ispec=1,numpts
        WRITE(IWR,*)symlabels(ispec),':',
     $    (symcoeffs(ivec,ispec),ivec=1,3)
      END DO
      WRITE(IWR,*)
c     Special points from Setyawan and Curtarolo, Comp. Mat. Sci. 49,
c     (2010), 299, 312.
      DO ispec=1,numpts
        DO icart = 1,3
          gspec(icart,ispec)=0.0d0
          DO ivec = 1,3
            gspec(icart,ispec) = gspec(icart,ispec)
     $        + gprim(icart,ivec)*symcoeffs(ivec,ispec)
          END DO
        END DO
      END DO
c
c     Now write out the values of the special points
c
      DO ispec=1,numpts
         WRITE(IWR,*) '  ',symlabels(ispec),' ',
     $    (gspec(icart,ispec),icart=1,3)
         WRITE(IWR,*)' '
      ENDDO      
c Use default path if nbranch<=0
      IF (nbranch.LE.0) THEN
        WRITE(IWR,*)'Default '//latname//
     $    ' bandstructure option selected'
        bsymb(1)='GY'
        bsymb(2)='YH'
        bsymb(3)='HC'
        bsymb(4)='CE'
        bsymb(5)='Em'
        bsymb(6)='mA'
        bsymb(7)='AX'
        bsymb(8)='Xh'
        bsymb(9)='MD'
        bsymb(10)='DZ'
        bsymb(11)='YD'
        nbranch=11
      ENDIF
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         ibegin=LOOKUP_SYMPT(symb, symlabels, numpts)
         IF (ibegin.LT.1) THEN
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
         iend=LOOKUP_SYMPT(syme, symlabels, numpts)
         IF (iend.LT.1) THEN
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ. LOOKUP_SYMPT('G', symlabels, numpts)) THEN 
            IF(iend.EQ. LOOKUP_SYMPT('Y', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('y', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('A', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('m', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('n', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('C', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('H', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('h', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('E', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('D', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('Z', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Y', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('d', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('Z', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('n', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('y', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('E', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('m', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('C', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('H', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('C', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Y', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('h', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('C', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('m', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('i', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('y', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('n', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('M', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('H', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('D', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('E', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('m', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('A', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('E', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('h', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('n', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('A', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('d', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('i', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('X', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('A', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('h', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('i', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('Y', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('H', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('D', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('y', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('i', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('d', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('Z', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('D', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('d', symlabels, numpts )) GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)thisroutine//': d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/dfloat(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)thisroutine//': nnscf = ',nnscf
      RETURN
      END

      SUBROUTINE BANDSMCL2(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, rprim, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (numpts=16)
c
c     Simple Monoclinic Bravais Lattice
c      Ordering of the lattice a,b<=c; alpha > 90 deg; beta=gamma=90 deg
c     Lattice Setting
c      a1=(a, 0, 0)
c      a2=(0, b, 0)
c      a3=(0, c*cos(alpha), c*sin(alpha))
c     
      CHARACTER*10 thisroutine
      CHARACTER*4  latname
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER symlabels(numpts),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),rprim(3,3),gprim(3,3)
      DIMENSION gspec(3,numpts),dx(3),dk(nknscfd),nbr(nbranchd)
      DIMENSION symcoeffs(3,numpts)
      DATA thisroutine /'BANDSMCL2'/
      DATA latname     /'MCL'/
c
c     symlabels:  GAMMA-A-C-D-D1-E-H-H1-H2-M-M1-M2-X-Y-Y1-Z
c
      DATA symlabels/'G','A','C',
     $               'D','d','E',
     $               'H','h','i',
     $               'M','m','n',
     $               'X','Y','y','Z'/
c
      DATA symcoeffs / 0.0d0,  0.0d0,  0.0d0,
     $                 0.5d0,  0.5d0,  0.0d0,
     $                 0.0d0,  0.5d0, -0.5d0,
     $                 0.5d0,  0.0d0,  0.5d0,
     $                 0.5d0,  0.0d0, -0.5d0,
     $                 0.5d0,  0.5d0, -0.5d0,
     $                 0.0d0,  100.0d0,  100.0d0,
     $                 0.0d0,  100.0d0,  100.0d0,
     $                 0.0d0,  100.0d0,  100.0d0,
     $                 0.5d0,  100.0d0,  100.0d0,
     $                 0.5d0,  100.0d0,  100.0d0,
     $                 0.5d0,  100.0d0,  100.0d0,
     $                 0.0d0,  0.5d0,  0.0d0,
     $                 0.0d0,  0.0d0,  0.5d0,
     $                 0.0d0,  0.0d0, -0.5d0,
     $                 0.5d0,  0.0d0,  0.0d0 /
c
      np1=30
c
c Compute special points that depend on rprim.
c
      a=DDOT(3,rprim(1,1),1,rprim(1,1),1)
      a=SQRT(a)
      b=DDOT(3,rprim(1,2),1,rprim(1,2),1)
      b=SQRT(b)
      c=DDOT(3,rprim(1,3),1,rprim(1,3),1)
      c=SQRT(c)
c
      cosalpha=rprim(2,3)/c
      IF (cosalpha.GT.0.0d0) THEN
        STOP 'BANDSMCL2: cos(alpha) > 0.0d0'
      ENDIF
c
      sinalpha=rprim(3,3)/c
      should_be_zero=MAX(ABS(rprim(2,1)),ABS(rprim(3,1)),
     $                   ABS(rprim(1,2)),ABS(rprim(3,3)),
     $                   ABS(rprim(1,3)),
     $                   ABS(1.0d0-cosalpha**2-sinalpha**2))
c
      IF(should_be_zero.GT.1.d-8) THEN
        STOP 'BANDSMCL2: error in rprim'
      ENDIF
c
      eta=0.5d0*(1.0d0+b*cosalpha/c)/sinalpha/sinalpha
      gnu=0.5d0-eta*c*cosalpha/b
c  H
      symcoeffs(1,7)=0.0d0
      symcoeffs(2,7)=eta
      symcoeffs(3,7)=1.0d0-gnu
c  H1=h
      symcoeffs(1,8)=0.0d0
      symcoeffs(2,8)=1.0d0-eta
      symcoeffs(3,8)=-(1.0d0-gnu)
c  H2=i
      symcoeffs(1,9)=0.0d0
      symcoeffs(2,9)=eta
      symcoeffs(3,9)=-gnu
c  M
      symcoeffs(1,10)=0.5d0
      symcoeffs(2,10)=eta
      symcoeffs(3,10)=1.0d0-gnu
c  M1=m
      symcoeffs(1,11)=0.5d0
      symcoeffs(2,11)=1.0-eta
      symcoeffs(3,11)=-(1.0d0-gnu)
c  M2=n
      symcoeffs(1,12)=0.5d0
      symcoeffs(2,12)=eta
      symcoeffs(3,12)=-gnu
c
c     Use the primitive lattice vectors to make special points
c
      WRITE(IWR,*)thisroutine//': primitive reciprical lattice vectors'
      DO ivec=1,3
         WRITE(IWR,*)(gprim(icart,ivec),icart=1,3)
      ENDDO
c
      WRITE(IWR,*)
      WRITE(IWR,*)thisroutine//': symcoeffs'
      DO ispec=1,numpts
        WRITE(IWR,*)symlabels(ispec),':',
     $    (symcoeffs(ivec,ispec),ivec=1,3)
      END DO
      WRITE(IWR,*)
c     Special points from Setyawan and Curtarolo, Comp. Mat. Sci. 49,
c     (2010), 299, 312.
      DO ispec=1,numpts
        DO icart = 1,3
          gspec(icart,ispec)=0.0d0
          DO ivec = 1,3
            gspec(icart,ispec) = gspec(icart,ispec)
     $        + gprim(icart,ivec)*symcoeffs(ivec,ispec)
          END DO
        END DO
      END DO
c
c     Now write out the values of the special points
c
      DO ispec=1,numpts
         WRITE(IWR,*) '  ',symlabels(ispec),' ',
     $    (gspec(icart,ispec),icart=1,3)
         WRITE(IWR,*)' '
      ENDDO      
c Use default path if nbranch<=0
      IF (nbranch.LE.0) THEN
        WRITE(IWR,*)'Default '//latname//
     $    ' bandstructure option selected'
        bsymb(1)='GY'
        bsymb(2)='YH'
        bsymb(3)='HX'
        bsymb(4)='XA'
        bsymb(5)='Am'
        bsymb(6)='mE'
        bsymb(7)='EC'
        bsymb(8)='Ch'
        bsymb(9)='MD'
        bsymb(10)='DZ'
        bsymb(11)='YD'
        nbranch=11
      ENDIF
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         ibegin=LOOKUP_SYMPT(symb, symlabels, numpts)
         IF (ibegin.LT.1) THEN
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
         iend=LOOKUP_SYMPT(syme, symlabels, numpts)
         IF (iend.LT.1) THEN
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ. LOOKUP_SYMPT('G', symlabels, numpts)) THEN 
c               GAMMA - Z, Y, Y1
            IF(iend.EQ. LOOKUP_SYMPT('Z', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Y', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('y', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('A', symlabels, numpts)) THEN
c               A  - X, M, M1
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('m', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('C', symlabels, numpts)) THEN
c               C  - E, H1, H2
            IF(iend.EQ. LOOKUP_SYMPT('E', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('h', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('i', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('D', symlabels, numpts)) THEN
c               D  - Y, Z, M
            IF(iend.EQ. LOOKUP_SYMPT('Y', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Z', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('d', symlabels, numpts)) THEN
c                D1 - Y1, Z, M2
            IF(iend.EQ. LOOKUP_SYMPT('y', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Z', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('n', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('E', symlabels, numpts)) THEN
c                E  - C, M1, M2
            IF(iend.EQ. LOOKUP_SYMPT('C', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('m', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('n', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('H', symlabels, numpts)) THEN
c                H  - M, X, Y
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Y', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('h', symlabels, numpts)) THEN
c                H1 - M1, X, C
            IF(iend.EQ. LOOKUP_SYMPT('m', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('C', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('i', symlabels, numpts)) THEN
c                H2 - M2, C, Y1
            IF(iend.EQ. LOOKUP_SYMPT('n', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('C', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('y', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('M', symlabels, numpts)) THEN
c                M  - H, D, A
            IF(iend.EQ. LOOKUP_SYMPT('H', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('D', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('A', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('m', symlabels, numpts)) THEN
c                M1 - H1, A, E
            IF(iend.EQ. LOOKUP_SYMPT('h', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('A', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('E', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('n', symlabels, numpts)) THEN
c                M2 -  H2, E, D1
            IF(iend.EQ. LOOKUP_SYMPT('i', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('E', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('d', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('X', symlabels, numpts)) THEN
c                X  -  A, H, H1
            IF(iend.EQ. LOOKUP_SYMPT('A', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('H', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('h', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('Y', symlabels, numpts)) THEN
c                Y  -  D, GAMMA, H
            IF(iend.EQ. LOOKUP_SYMPT('D', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('H', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('y', symlabels, numpts)) THEN
c                Y1 - D1, GAMMA, H2
            IF(iend.EQ. LOOKUP_SYMPT('d', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('i', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('Z', symlabels, numpts)) THEN
c                Z - GAMMA, D, D1
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('D', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('d', symlabels, numpts )) GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)thisroutine//': d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/dfloat(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)thisroutine//': nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSORC
c
c
c REVISION HISTORY:
c   4Sep12-ACP: Forked from BANDSHEX
c   3May12-ACP: Merged into 2.64
c  20Dec10-ACP: Fully merged into 2.62dev.
c  07Oct08-AHE: Written for 2.61b
c
      SUBROUTINE BANDSORC(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (numpts=8)
c
c     Orthrhombic Bravais Lattice
c     Assumed Lattice ordering:   a<b<c
c        a1=(a,0,0)
c        a2=(0,b,0)
c        a3=(0,0,c)
c     One of a suite of routines that calculate a set of branches for
c     band structure calculations.
c     
c     First attempt to parse the bsymb routine into constituent
c     characters, and check against the set of special points in ORC.
c
      CHARACTER*10 thisroutine
      CHARACTER*4  latname
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER symlabels(numpts),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),gprim(3,3)
      DIMENSION gspec(3,numpts),dx(3),dk(nknscfd),nbr(nbranchd)
      DIMENSION symcoeffs(3,numpts)
      DATA thisroutine /'BANDSORC'/
      DATA latname     /'ORC'/
      DATA symlabels/'G','R','S','T','U','X','Y','Z'/
      DATA symcoeffs / 0.0d0, 0.0d0, 0.0d0,
     $                 0.5d0, 0.5d0, 0.5d0,
     $                 0.5d0, 0.5d0, 0.0d0,
     $                 0.0d0, 0.5d0, 0.5d0,
     $                 0.5d0, 0.0d0, 0.5d0,
     $                 0.5d0, 0.0d0, 0.0d0,
     $                 0.0d0, 0.5d0, 0.0d0,
     $                 0.0d0, 0.0d0, 0.5d0 /
      np1=30
c
c     Use the primitive lattice vectors to make special points
c
      WRITE(IWR,*)thisroutine//': primitive reciprical lattice vectors'
      DO ivec=1,3
         WRITE(IWR,*)(gprim(icart,ivec),icart=1,3)
      ENDDO
      WRITE(IWR,*)
      WRITE(IWR,*)thisroutine//': symcoeffs'
      DO ispec=1,numpts
        WRITE(IWR,*)symlabels(ispec),':',
     $    (symcoeffs(ivec,ispec),ivec=1,3)
      END DO
      WRITE(IWR,*)
c     Special points from Setyawan and Curtarolo, Comp. Mat. Sci. 49,
c     (2010), 299, 312.
      DO ispec=1,numpts
        DO icart = 1,3
          gspec(icart,ispec)=0.0d0
          DO ivec = 1,3
            gspec(icart,ispec) = gspec(icart,ispec)
     $        + gprim(icart,ivec)*symcoeffs(ivec,ispec)
          END DO
        END DO
      END DO
c
c     Now write out the values of the special points
c
      DO ispec=1,numpts
         WRITE(IWR,*) '  ',symlabels(ispec),' ',
     $    (gspec(icart,ispec),icart=1,3)
         WRITE(IWR,*)' '
      ENDDO      
c Use default path if nbranch<=0
      IF (nbranch.LE.0) THEN
        WRITE(IWR,*)'Default '//latname//
     $    ' bandstructure option selected'
        bsymb(1)='GX'
        bsymb(2)='XS'
        bsymb(3)='SY'
        bsymb(4)='YG'
        bsymb(5)='GZ'
        bsymb(6)='ZU'
        bsymb(7)='UR'
        bsymb(8)='RT'
        bsymb(9)='TZ'
        bsymb(10)='YT'
        bsymb(11)='UX'
        bsymb(12)='SR'
        nbranch=12
      ENDIF
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         ibegin=LOOKUP_SYMPT(symb, symlabels, numpts)
         IF(ibegin.LT.1) THEN
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'not in the ORC set of special points.  Bye!'
            STOP
         ENDIF
         iend=LOOKUP_SYMPT(syme, symlabels, numpts)
         IF(iend.LT.1) THEN
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'not in the ORC set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ. LOOKUP_SYMPT('G', symlabels, numpts)) THEN 
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Y', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Z', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('R', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('S', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('T', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('U', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('S', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('R', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Y', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('T', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('R', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Y', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Z', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('U', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('R', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Z', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('X', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('U', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('S', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('Y', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('S', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('T', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('Z', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('U', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('T', symlabels, numpts )) GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)thisroutine//': d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/DFLOAT(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)thisroutine//': nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSRHL
c
c
c REVISION HISTORY:
c   4Sep12-ACP: Forked from BANDSHEX
c   3May12-ACP: Merged into 2.64
c  20Dec10-ACP: Fully merged into 2.62dev.
c  07Oct08-AHE: Written for 2.61b
c
      SUBROUTINE BANDSRHL(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, rprim, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
c   Variation 1 has 12 symmetry points
      PARAMETER (numpts1=12)
      PARAMETER (numpath1=9)
c   Variation 2 has 8 symmetry points
      PARAMETER (numpts2=8)
      PARAMETER (numpath2=9)
c 
      PARAMETER (numptsmax=12)
c
c     Rhombohedral Lattice
c     Lattice Setting
c      a1=(a*cos(alpha/2), -a*sin(alpha/2), 0)
c      a2=(a*cos(alpha/2),  a*sin(alpha/2), 0)
c      a3=(a*cos(alpha)/cos(alpha/2), 0, a*sqrt(1-cos^2(alpha)/cos^2(alpha/2)))
c     
c     There are two variations alpha<90 (cos(alpha)>0) and alpha>90 (cos(alpha)<0).
      CHARACTER*10 thisroutine
      CHARACTER*4  latname
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER*2 RHL1PATH(numpath1), RHL2PATH(numpath2)
      CHARACTER symlabels(numptsmax), symb, syme, symeold
      CHARACTER RHL1LABELS(numpts1)
      CHARACTER RHL2LABELS(numpts2)
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),rprim(3,3),gprim(3,3)
      DIMENSION gspec(3,numptsmax),dx(3),dk(nknscfd),nbr(nbranchd)
      DIMENSION symcoeffs(3,numptsmax)
      DATA thisroutine /'BANDSRHL'/
      DATA latname     /'RHL'/
c     Special points from Setyawan and Curtarolo, Comp. Mat. Sci. 49,
c     (2010), 299, 312.
c     cos(alpha)>0
c     RHL1LABELS:  GAMMA-B-B1-F-L-L1-P-P1-P2-Q-X-Z
      DATA RHL1LABELS/'G','B','b',
     $                'F','L','l',
     $                'P','p','o',
     $                'Q','X','Z'/
      DATA RHL1PATH /'GL', 'Lb', 'BZ',
     $               'ZG', 'GX', 'QF',
     $               'Fp', 'pZ', 'LP' /
c     cos(alpha)<0
c     RHL2LABELS: GAMMA-F-L-P-P1-Q-Q1-Z
      DATA RHL2LABELS /'G','F','L','P',
     $                 'p','Q','q','Z' /
      DATA RHL2PATH /'GP', 'PZ', 'ZQ',
     $               'QG', 'GF', 'Fp',
     $               'pq', 'qL', 'LZ' /
c
      np1=30
      should_be_zero=MAX(ABS(rprim(3,1)),
     $                   ABS(rprim(3,2)),
     $                   ABS(rprim(2,3)))
      IF(should_be_zero.GT.1.d-8) THEN
        STOP 'BANDSRHL: error in rprim'
      ENDIF
c
c Figure out which variation to use.
c
      a=DDOT(3,rprim(1,1),1,rprim(1,1),1)
      a=SQRT(a)
      a1a3=DDOT(3,rprim(1,1),1,rprim(1,3),1)
      cosalpha=a1a3/a**2
      tanhalfalpha=rprim(2,2)/rprim(1,2)
      IF (cosalpha.GT.0.0d0) THEN
c  Variation RHL1
        eta=(1.0d0+4.0d0*cosalpha)/(2.0d0+4.0d0*cosalpha)
        gnu=0.75d0-eta*0.5d0
        numpts=numpts1
        DO ispec=1,numpts
          symlabels(ispec)=RHL1LABELS(ispec)
        END DO
c
        ispec=LOOKUP_SYMPT('G', symlabels, numpts)
        symcoeffs(1,ispec)=0.0d0
        symcoeffs(2,ispec)=0.0d0
        symcoeffs(3,ispec)=0.0d0
c
        ispec=LOOKUP_SYMPT('B', symlabels, numpts)
        symcoeffs(1,ispec)=eta
        symcoeffs(2,ispec)=0.5d0
        symcoeffs(3,ispec)=1-eta
c
        ispec=LOOKUP_SYMPT('b', symlabels, numpts)
        symcoeffs(1,ispec)=0.5d0
        symcoeffs(2,ispec)=1.d0-eta
        symcoeffs(3,ispec)=eta-1.0d0
c
        ispec=LOOKUP_SYMPT('F', symlabels, numpts)
        symcoeffs(1,ispec)=0.5d0
        symcoeffs(2,ispec)=0.5d0
        symcoeffs(3,ispec)=0.0d0
c
        ispec=LOOKUP_SYMPT('L', symlabels, numpts)
        symcoeffs(1,ispec)=0.5d0
        symcoeffs(2,ispec)=0.0d0
        symcoeffs(3,ispec)=0.0d0
c
        ispec=LOOKUP_SYMPT('l', symlabels, numpts)
        symcoeffs(1,ispec)=0.0d0
        symcoeffs(2,ispec)=0.0d0
        symcoeffs(3,ispec)=-0.5d0
c
        ispec=LOOKUP_SYMPT('P', symlabels, numpts)
        symcoeffs(1,ispec)=eta
        symcoeffs(2,ispec)=gnu
        symcoeffs(3,ispec)=gnu
c
        ispec=LOOKUP_SYMPT('p', symlabels, numpts)
        symcoeffs(1,ispec)=1.0d0-gnu
        symcoeffs(2,ispec)=1.0d0-gnu
        symcoeffs(3,ispec)=1.0d0-eta
c
        ispec=LOOKUP_SYMPT('o', symlabels, numpts)
        symcoeffs(1,ispec)=gnu
        symcoeffs(2,ispec)=gnu
        symcoeffs(3,ispec)=eta-1.0d0
c
        ispec=LOOKUP_SYMPT('Q', symlabels, numpts)
        symcoeffs(1,ispec)=1.0d0-gnu
        symcoeffs(2,ispec)=gnu
        symcoeffs(3,ispec)=0.0d0
c
        ispec=LOOKUP_SYMPT('X', symlabels, numpts)
        symcoeffs(1,ispec)=gnu
        symcoeffs(2,ispec)=0.0d0
        symcoeffs(3,ispec)=-gnu
c
        ispec=LOOKUP_SYMPT('Z', symlabels, numpts)
        symcoeffs(1,ispec)=0.5d0
        symcoeffs(2,ispec)=0.5d0
        symcoeffs(3,ispec)=0.5d0
c
c Use default path if nbranch<=0
        IF (nbranch.LE.0) THEN
          WRITE(IWR,*)'Default '//latname//
     $      '1 bandstructure option selected'
          nbranch=numpath1
          DO i=1,nbranch
            bsymb(i)=RHL1PATH(i)
          END DO
        ENDIF
      ELSEIF (cosalpha.LT.0.0d0) THEN
c  Variation RHL2
        eta=1.d0/(2.d0*tanhalfalpha**2)
        gnu=0.75d0-eta*0.5d0
        numpts=numpts2
        DO ispec=1,numpts
          symlabels(ispec)=RHL2LABELS(ispec)
        END DO
c
        ispec=LOOKUP_SYMPT('G', symlabels, numpts)
        symcoeffs(1,ispec)=0.0d0
        symcoeffs(2,ispec)=0.0d0
        symcoeffs(3,ispec)=0.0d0
c
        ispec=LOOKUP_SYMPT('F', symlabels, numpts)
        symcoeffs(1,ispec)= 0.5d0
        symcoeffs(2,ispec)=-0.5d0
        symcoeffs(3,ispec)= 0.0d0
c
        ispec=LOOKUP_SYMPT('L', symlabels, numpts)
        symcoeffs(1,ispec)=0.5d0
        symcoeffs(2,ispec)=0.0d0
        symcoeffs(3,ispec)=0.0d0
c
        ispec=LOOKUP_SYMPT('P', symlabels, numpts)
        symcoeffs(1,ispec)=1.0d0-gnu
        symcoeffs(2,ispec)=-gnu
        symcoeffs(3,ispec)=1.0d0-gnu
c
        ispec=LOOKUP_SYMPT('p', symlabels, numpts)
        symcoeffs(1,ispec)=gnu
        symcoeffs(2,ispec)=gnu-1.0d0
        symcoeffs(3,ispec)=gnu-1.0d0
c
        ispec=LOOKUP_SYMPT('Q', symlabels, numpts)
        symcoeffs(1,ispec)=eta
        symcoeffs(2,ispec)=eta
        symcoeffs(3,ispec)=eta
c
        ispec=LOOKUP_SYMPT('q', symlabels, numpts)
        symcoeffs(1,ispec)=1.0d0-eta
        symcoeffs(2,ispec)=-eta
        symcoeffs(3,ispec)=-eta
c
        ispec=LOOKUP_SYMPT('Z', symlabels, numpts)
        symcoeffs(1,ispec)= 0.5d0
        symcoeffs(2,ispec)=-0.5d0
        symcoeffs(3,ispec)= 0.5d0
c
c Use default path if nbranch<=0
        IF (nbranch.LE.0) THEN
          WRITE(IWR,*)'Default '//latname//
     $      '2 bandstructure option selected'
          nbranch=numpath2
          DO i=1,nbranch
            bsymb(i)=RHL2PATH(i)
          END DO
        ENDIF
      ELSE
        STOP "Wrong lattice vectors"
      ENDIF
c

c
c     Use the primitive lattice vectors to make special points
c
      WRITE(IWR,*)thisroutine//': primitive reciprical lattice vectors'
      DO ivec=1,3
         WRITE(IWR,*)(gprim(icart,ivec),icart=1,3)
      ENDDO
c
      WRITE(IWR,*)
      WRITE(IWR,*)thisroutine//': symcoeffs'
      DO ispec=1,numpts
        WRITE(IWR,*)symlabels(ispec),':',
     $    (symcoeffs(ivec,ispec),ivec=1,3)
      END DO
c
      WRITE(IWR,*)
      DO ispec=1,numpts
        DO icart = 1,3
          gspec(icart,ispec)=0.0d0
          DO ivec = 1,3
            gspec(icart,ispec) = gspec(icart,ispec)
     $        + gprim(icart,ivec)*symcoeffs(ivec,ispec)
          END DO
        END DO
      END DO
c
c     Now write out the values of the special points
c
      DO ispec=1,numpts
         WRITE(IWR,*) '  ',symlabels(ispec),' ',
     $    (gspec(icart,ispec),icart=1,3)
         WRITE(IWR,*)' '
      ENDDO      
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         ibegin=LOOKUP_SYMPT(symb, symlabels, numpts)
         IF (ibegin.LT.1) THEN
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
         iend=LOOKUP_SYMPT(syme, symlabels, numpts)
         IF (iend.LT.1) THEN
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path (SKIPPED)
c
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)thisroutine//': d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/dfloat(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)thisroutine//': nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSSC
c
c
c REVISION HISTORY:
c   3May12-ACP: Merged into 2.64
c  14Oct08-AHE: Written for 2.61b
c  20Dec10-ACP: Fully merged into 2.62dev.
c
      SUBROUTINE BANDSSC(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
cAHE 101408
c
c     One of a suite of routines that calculate a set of branches for band structure 
c     calculations.
c     
c     First attempt to parse the bsymb routine into constituent characters, and check
c     against the set of special points in SC
c
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER SC(4),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),gprim(3,3)
      DIMENSION gspec(3,6),dx(3),dk(nknscfd),nbr(nbranchd)
      DATA SC/'G','X','M','R'/
      np1=60
c
c     Use the primitive lattice vectors to make special points
c
      DO i = 1,3
c     G (gamma)
         gspec(i,1)=0.0d0
c     R
         gspec(i,4)=(gprim(i,1)+gprim(i,2)+gprim(i,3))*0.5
c     X
         gspec(i,2)=(gprim(i,1))*0.5
c     M
         gspec(i,3)=(gprim(i,1)+gprim(i,2))*0.5
c     
      ENDDO
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         IF(symb.EQ.'G') THEN
            ibegin=1
         ELSEIF(symb.EQ.'X') THEN
            ibegin=2
         ELSEIF(symb.EQ.'M')THEN
            ibegin=3
         ELSEIF(symb.EQ.'R')THEN
            ibegin=4
         ELSE
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'not in the SC set of special points.  Bye!'
            STOP
         ENDIF
         IF(syme.EQ.'G') THEN
            iend=1
         ELSEIF(syme.EQ.'X') THEN
            iend=2
         ELSEIF(syme.EQ.'M')THEN
            iend=3
         ELSEIF(syme.EQ.'R')THEN
            iend=4
         ELSE
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'not in the SC set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ.1.AND. iend.NE.1) THEN 
            GO TO 10
         ELSEIF(ibegin.EQ.2.AND.iend.NE.2) THEN
            GO TO 10
         ELSEIF(ibegin.EQ.3.AND.iend.NE.3) THEN
            GO TO 10
         ELSEIF(ibegin.EQ.4.AND.iend.NE.4) THEN
            GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)'BANDSSC: d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/DFLOAT(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)'in BANDSC, nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSSTC
c
c
c REVISION HISTORY:
c   3May12-ACP: Merged into 2.64
c  04Jan09-AHE: Written for 2.61b.
c  20Dec10-ACP: Fully merged into 2.62dev.
c
      SUBROUTINE BANDSSTC(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
cAHE 010409
c
c     One of a suite of routines that calculate a set of branches for band structure 
c     calculations.
c     
c     This is the simple triclinic system, not particularly attractive. We start with 
c     one line from the gamma point to the midpoint along the
c     shortest reciprocal lattice vector. We will label thie GX
c     First attempt to parse the 
c
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER TC(4),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),gprim(3,3)
      DIMENSION gspec(3,6),dx(3),dk(nknscfd),nbr(nbranchd)
      DIMENSION gdis(3)
      DATA TC/'O','G','B','F'/
      np1=60
c
c     Use the primitive lattice vectors to make special points
c
c     Order the directions by size of the reciprocal vectors.
c     
      gmin=0.0d0
      DO  j=1,3
         gdis(j)=0.0d0
         DO i = 1,3
            gdis(i)=gdis(i)+gprim(i,j)*gprim(i,j)
         ENDDO
         gdis(j)=sqrt(gdis(j))
      ENDDO
      jmin=1
      jmax=1
      DO j = 2,3
         IF(gdis(j).LT.gdis(jmin))jmin=j
         IF(gdis(j).GT.gdis(jmax))jmax=j
      ENDDO
      DO j=1,3
         IF(j.NE.jmin.AND.j.NE.jmax) jmid=j
      ENDDO
      DO i = 1,3
c     O (gamma)
         gspec(i,1)=0.0d0
c     G (minimum length)
         gspec(i,2)=gprim(i,jmin)*0.5
c     B (middle length)
         gspec(i,3)=gprim(i,jmid)*0.5
c     F (maximum length)
         gspec(i,4)=gprim(i,jmax)*0.5
      ENDDO
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.eq.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         IF(symb.EQ.'O') THEN
            ibegin=1
         ELSEIF(symb.EQ.'G') THEN
            ibegin=2
         ELSEIF(symb.EQ.'B') THEN
            ibegin=3
         ELSEIF(symb.EQ.'F') THEN
            ibegin=4
         ELSE
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'not in the STC set of special points.  Bye!'
            STOP
         ENDIF
         IF(syme.EQ.'O') THEN
            iend=1
         ELSEIF(syme.EQ.'G') THEN
            iend=2
         ELSEIF(syme.EQ.'B') THEN
            iend=3
         ELSEIF(syme.EQ.'F') THEN
            iend=4
         ELSE
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'not in the STC set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ.1.AND. iend.NE.1) THEN 
            GO TO 10
         ELSEIF(ibegin.EQ.2.AND.iend.NE.2) THEN
            GO TO 10
         ELSEIF(ibegin.EQ.3.AND.iend.NE.3) THEN
            GO TO 10
         ELSEIF(ibegin.EQ.4.AND.iend.NE.4) THEN
            GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)'BANDSSTC: d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/DFLOAT(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)'in BANDSTC, nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSTC
c
c
c REVISION HISTORY:
c   3May12-ACP: Merged into 2.64
c  04Jan09-AHE: Written for 2.61b.
c  20Dec10-ACP: Fully merged into 2.62dev.
c
      SUBROUTINE BANDSTC(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
cAHE 010409
c
c     One of a suite of routines that calculate a set of branches for band structure 
c     calculations.
c     
c     This is the triclinic system, not particularly attractive. We start with 
c     one line from the gamma point to the midpoint along the
c     shortest reciprocal lattice vector. We will label thie GX
c     First attempt to parse the 
c
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER*1 TC(2)
      CHARACTER SC(2),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),gprim(3,3)
      DIMENSION gspec(3,6),dx(3),dk(nknscfd),nbr(nbranchd)
      DIMENSION gdis(3)
      DATA TC/'G','X'/
      np1=60
c
c     Use the primitive lattice vectors to make special points
c
c     Find the shortest of the gprims
c     
      gmin=0.0
      DO  j=1,3
         gdis(j)=0.0d0
         DO i = 1,3
            gdis(i)=gdis(i)+gprim(i,j)*gprim(i,j)
         ENDDO
         gdis(j)=SQRT(gdis(j))
      ENDDO
      jmin=1
      DO j = 2,3
         IF(gdis(j).lt.gdis(jmin))jmin=j
      ENDDO
      DO i = 1,3
c     G (gamma)
         gspec(i,1)=0.0d0
c     X
         gspec(i,2)=(gprim(i,jmin))*0.5
      ENDDO
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         IF(symb.EQ.'G') THEN
            ibegin=1
         ELSEIF(symb.EQ.'X') THEN
            ibegin=2
         ELSE
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'not in the TC set of special points.  Bye!'
            STOP
         ENDIF
         IF(syme.EQ.'G') THEN
            iend=1
         ELSEIF(syme.EQ.'X') THEN
            iend=2
         ELSE
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'not in the SC set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ.1.AND. iend.NE.1) THEN 
            GO TO 10
         ELSEIF(ibegin.EQ.2.AND.iend.NE.2) THEN
            GO TO 10
         ELSEIF(ibegin.EQ.3.AND.iend.NE.3) THEN
            GO TO 10
         ELSEIF(ibegin.EQ.4.AND.iend.NE.4) THEN
            GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)'BANDSTC: d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/DFLOAT(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)'in BANDTC, nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSTET
c
c
c REVISION HISTORY:
c   4Sep12-ACP: Forked from BANDSHEX
c   3May12-ACP: Merged into 2.64
c  20Dec10-ACP: Fully merged into 2.62dev.
c  07Oct08-AHE: Written for 2.61b
c
      SUBROUTINE BANDSTET(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (numpts=6)
c
c     Tetragonal Bravais Lattice
c     One of a suite of routines that calculate a set of branches for
c     band structure calculations.
c     
c     Conventional Lattice Setting
c      a1=(a,0,0)
c      a2=(0,a,0)
c      a3=(0,0,c)
c
c     First attempt to parse the bsymb routine into constituent
c     characters, and check against the set of special points in TET.
c
      CHARACTER*10 thisroutine
      CHARACTER*4  latname
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER symlabels(numpts),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),gprim(3,3)
      DIMENSION gspec(3,numpts),dx(3),dk(nknscfd),nbr(nbranchd)
      DIMENSION symcoeffs(3,numpts)
      DATA thisroutine /'BANDSTET'/
      DATA latname     /'TET'/
      DATA symlabels/'G','A','M','R','X','Z'/
      DATA symcoeffs / 0.0d0, 0.0d0, 0.0d0,
     $                 0.5d0, 0.5d0, 0.5d0,
     $                 0.5d0, 0.5d0, 0.0d0,
     $                 0.0d0, 0.5d0, 0.5d0,
     $                 0.0d0, 0.5d0, 0.0d0,
     $                 0.0d0, 0.0d0, 0.5d0 /
      np1=30
c
c     Use the primitive lattice vectors to make special points
c
      WRITE(IWR,*)thisroutine//': primitive reciprical lattice vectors'
      DO ivec=1,3
         WRITE(IWR,*)(gprim(icart,ivec),icart=1,3)
      ENDDO
      WRITE(IWR,*)
      WRITE(IWR,*)thisroutine//': symcoeffs'
      DO ispec=1,numpts
        WRITE(IWR,*)symlabels(ispec),':',
     $    (symcoeffs(ivec,ispec),ivec=1,3)
      END DO
      WRITE(IWR,*)
c     Special points from Setyawan and Curtarolo, Comp. Mat. Sci. 49,
c     (2010), 299, 312.
      DO ispec=1,numpts
        DO icart = 1,3
          gspec(icart,ispec)=0.0d0
          DO ivec = 1,3
            gspec(icart,ispec) = gspec(icart,ispec)
     $        + gprim(icart,ivec)*symcoeffs(ivec,ispec)
          END DO
        END DO
      END DO
c
c     Now write out the values of the special points
c
      DO ispec=1,numpts
         WRITE(IWR,*) '  ',symlabels(ispec),' ',
     $    (gspec(icart,ispec),icart=1,3)
         WRITE(IWR,*)' '
      ENDDO      
c Use default path if nbranch<=0
      IF (nbranch.LE.0) THEN
        WRITE(IWR,*)'Default '//latname//
     $    ' bandstructure option selected'
        bsymb(1)='GX'
        bsymb(2)='XM'
        bsymb(3)='MG'
        bsymb(4)='GZ'
        bsymb(5)='ZR'
        bsymb(6)='RA'
        bsymb(7)='AZ'
        bsymb(8)='XR'
        bsymb(9)='MA'
        nbranch=9
      ENDIF
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         ibegin=LOOKUP_SYMPT(symb, symlabels, numpts)
         IF (ibegin.LT.1) THEN
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
         iend=LOOKUP_SYMPT(syme, symlabels, numpts)
         IF (iend.LT.1) THEN
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ. LOOKUP_SYMPT('G', symlabels, numpts)) THEN 
            IF(iend.EQ. LOOKUP_SYMPT('Z', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('Z', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('A', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('R', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('A', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Z', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('R', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('R', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('Z', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('A', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('M', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('A', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.EQ.LOOKUP_SYMPT('X', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('R', symlabels, numpts )) GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)thisroutine//': d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c Write out last band label
      WRITE(IBLABEL,*)syme,dk(nnscf)
      CLOSE(IBLABEL,status='keep')
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/DFLOAT(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)thisroutine//': nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BANDSTRI
c
c
c REVISION HISTORY:
c   4Sep12-ACP: Forked from BANDSHEX
c   3May12-ACP: Merged into 2.64
c  20Dec10-ACP: Fully merged into 2.62dev.
c  07Oct08-AHE: Written for 2.61b
c
      SUBROUTINE BANDSTRI(nbranch,nbranchd,bsymb,veckn,wtkn,nnscf,
     $     nknscfd, gprim,dk,nbr,IWR)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (numpts=8)
c
c     Triclinic Bravais Lattice
c     One of a suite of routines that calculate a set of branches for
c     band structure calculations.
c
c     Conventional Triclinic Lattice Setting (a<b<c)
c
c      a1=(a,0,0)
c      a2=(b*cos(gamma),b*sin(gamma),0)
c      a3=(c*cos(beta),c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma),
c        (c/sin(gamma))*sqrt(sin(gamma)**2-cos(alpha)**2-cos(beta)**2+2*cos(alpha)*cos(beta)*cos(gamma)))
c
c      a1.a2=a*b*cos(gamma)
c      a1.a3=a*c*cos(beta)
c      a2.a3=b*c*cos(alpha)
c     
c     There are 4 variations that depend on the angles between the
c     reciprocal lattice vectors. The variations break into two sets of 2
c     variations that share the same coefficients for symmetry points.
c
      CHARACTER*10 thisroutine
      CHARACTER*4  latname
      CHARACTER*2 bsymb(nbranchd),symd
      CHARACTER symlabels(numpts),symb,syme,symeold
      DIMENSION veckn(3,nknscfd),wtkn(nknscfd),gprim(3,3)
      DIMENSION gspec(3,numpts),dx(3),dk(nknscfd),nbr(nbranchd)
c               symcoeffs(icart,ispec,ivariation)
      DIMENSION symcoeffs(3,numpts,2)
      DATA thisroutine /'BANDSTRI'/
      DATA latname     /'TRI'/
      DATA symlabels/'G','L','M','N','R','X','Y','Z'/
      DATA symcoeffs /
     $  0.0d0,  0.0d0,  0.0d0,
     $  0.5d0,  0.5d0,  0.0d0,
     $  0.0d0,  0.5d0,  0.5d0,
     $  0.5d0,  0.0d0,  0.5d0,
     $  0.5d0,  0.5d0,  0.5d0,
     $  0.5d0,  0.0d0,  0.0d0,
     $  0.0d0,  0.5d0,  0.0d0,
     $  0.0d0,  0.0d0,  0.5d0,
     $  0.0d0,  0.0d0,  0.0d0,
     $  0.5d0, -0.5d0,  0.0d0,
     $  0.0d0,  0.0d0,  0.5d0,
     $ -0.5d0, -0.5d0,  0.5d0,
     $  0.0d0, -0.5d0,  0.5d0,
     $  0.0d0, -0.5d0,  0.0d0,
     $  0.5d0,  0.0d0,  0.0d0,
     $ -0.5d0,  0.0d0,  0.5d0
     $ /
      np1=30
c
c     Use the primitive lattice vectors to make special points
c
      WRITE(IWR,*)thisroutine//': primitive reciprical lattice vectors'
      DO ivec=1,3
         WRITE(IWR,*)(gprim(icart,ivec),icart=1,3)
      ENDDO
      WRITE(IWR,*)
c
c     Decide which variation of the triclinic system we are.
c
      coskalpha=DDOT(3,gprim(1,2),1,gprim(1,3),1)
      coskbeta=DDOT(3,gprim(1,1),1,gprim(1,3),1)
      coskgamma=DDOT(3,gprim(1,1),1,gprim(1,2),1)
      smallabs=MIN(ABS(coskalpha), ABS(coskbeta), ABS(coskgamma))
      smallest=MIN(coskalpha, coskbeta, coskgamma)
      biggest=MAX(coskalpha, coskbeta, coskgamma)
      cutoff=1.0d-8
      IF(ABS(smallabs).GT.cutoff) THEN
        IF(biggest.GT.cutoff) THEN
c        TRI2b
          ivar=2
        ELSEIF (smallest.LT.-cutoff) THEN
c        TRI2a
          ivar=1
        ELSE
          STOP "Unrecognized variation"
        ENDIF
      ELSEIF(coskalpha.LT.0.0d0.AND.coskbeta.LT.0.0d0
     $    .AND.coskgamma.LT.0.0d0) THEN
c        TRI1a
        ivar=1
      ELSEIF(coskalpha.GT.0.0d0.AND.coskbeta.GT.0.0d0
     $    .AND.coskgamma.GT.0.0d0) THEN
c        TRI1b
        ivar=2
      ELSE
        STOP "Unrecognized variation"
      ENDIF

      WRITE(IWR,*) thisroutine//
     $  ": dot products of reciprocal lattice vectors" 
      WRITE(IWR,*)'cos(kalpha)=',coskalpha
      WRITE(IWR,*)'cos(kbeta)=',coskbeta
      WRITE(IWR,*)'cos(kgamma)=',coskgamma
      WRITE(IWR,*)
      WRITE(IWR,*)thisroutine//': symcoeffs(variation=',ivar,')'
      DO ispec=1,numpts
        WRITE(IWR,*)symlabels(ispec),':',
     $    (symcoeffs(ivec,ispec,ivar),ivec=1,3)
      END DO
      WRITE(IWR,*)
c     Special points from Setyawan and Curtarolo, Comp. Mat. Sci. 49,
c     (2010), 299, 312.
      
      DO ispec=1,numpts
        DO icart = 1,3
          gspec(icart,ispec)=0.0d0
          DO ivec = 1,3
            gspec(icart,ispec) = gspec(icart,ispec)
     $        + gprim(icart,ivec)*symcoeffs(ivec,ispec,ivar)
          END DO
        END DO
      END DO
c
c     Now write out the values of the special points
c
      DO ispec=1,numpts
         WRITE(IWR,*) '  ',symlabels(ispec),' ',
     $    (gspec(icart,ispec),icart=1,3)
         WRITE(IWR,*)' '
      ENDDO      
c Use default path if nbranch<=0
c The default path is the same for all variations.
      IF (nbranch.LE.0) THEN
        WRITE(IWR,*)'Default '//latname//
     $    ' bandstructure option selected'
        bsymb(1)='XG'
        bsymb(2)='GY'
        bsymb(3)='LG'
        bsymb(4)='GZ'
        bsymb(5)='NG'
        bsymb(6)='GM'
        bsymb(7)='RG'
        nbranch=7
      ENDIF
      nnscf=1
      dk(1)=0.0d0
      CALL FLOPENA(IBLABEL, 'bandlbls')
      symeold=' '
      DO ibr=1, nbranch
         symd=bsymb(ibr)
         symb=symd(1:1)
         iskip=0
         IF(ibr.GT.1) THEN
            IF(symb.EQ.syme) THEN
c     reduce number of points by one
               iskip=1
            ENDIF
         ENDIF
         syme=symd(2:2)
         WRITE(IWR,*)symb,'  ',syme
         ibegin=LOOKUP_SYMPT(symb, symlabels, numpts)
         IF (ibegin.lt.1) THEN
            WRITE(IWR,*)'symb = ', symb
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
         iend=LOOKUP_SYMPT(syme, symlabels, numpts)
         IF (iend.LT.1) THEN
            WRITE(IWR,*)'syme = ', syme
            WRITE(IWR,*)'Not in the '//latname//
     $        ' set of special points.  Bye!'
            STOP
         ENDIF
c
c     Now check that the endpoints are proper. Do they trace a well described path
c
         IF(ibegin.EQ. LOOKUP_SYMPT('G', symlabels, numpts)) THEN 
            IF(iend.EQ. LOOKUP_SYMPT('L', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('M', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('N', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('R', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('X', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Y', symlabels, numpts )) GO TO 10
            IF(iend.EQ. LOOKUP_SYMPT('Z', symlabels, numpts )) GO TO 10
         ELSEIF(ibegin.NE.LOOKUP_SYMPT('G', symlabels, numpts)) THEN
            IF(iend.EQ. LOOKUP_SYMPT('G', symlabels, numpts )) GO TO 10
         ELSE
            WRITE(IWR,*)symb,' and ',syme,' are not on a proper line.'
            WRITE(IWR,*)'ibegin = ', ibegin,'  iend = ',iend
            WRITE(IWR,*)'bye!'
            STOP
         ENDIF
 10      CONTINUE
c
c     Now build the points for the branch, and add them to list of NSCF points
c
         d1=0.0d0
         DO icart=1,3
            dx(icart)=(gspec(icart,iend)-gspec(icart,ibegin))/
     $           FLOAT(np1-1)
            d1=d1+dx(icart)*dx(icart)
         ENDDO
         d1=SQRT(d1)
         WRITE(IWR,*)thisroutine//': d1 = ', d1
c         write(IWR,*)'dx = ',(dx(icart),icart=1,3)
         IF(iskip.EQ.0) THEN
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)
            ENDDO
            IF(nnscf.EQ.1)THEN
               dk(nnscf)=0.0d0
            ELSE
               dk(nnscf)=dk(nnscf-1)
            ENDIF
         ELSE
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
         ENDIF
c Write out band label for this branch.
         IF(nnscf.eq.1) THEN
c First point in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf)
           ENDIF
         ELSE
c Intermediate points in the band structure
           IF((symeold.EQ.' ').OR.(symb.EQ.symeold)) THEN
             WRITE(IBLABEL,*)symb,dk(nnscf-1)
           ELSE
             WRITE(IBLABEL,*)symeold//symb,dk(nnscf-1)
           ENDIF
         ENDIF
         symeold=syme
c Write out the branch
         DO ik=2,np1-iskip
            nnscf=nnscf+1
            DO icart=1,3
               veckn(icart,nnscf)=gspec(icart,ibegin)+
     $              (ik-1+iskip)*dx(icart)
            ENDDO
            dk(nnscf)=dk(nnscf-1)+d1
c            write(IWR,*)'knscf(',nnscf,') = ',(veckn(icart,nnscf),
c     $           icart=1,3)
         ENDDO
      ENDDO
c
c     Now set the artificial weights so that the occupation is correct
c
      wt1=1.0d0/DFLOAT(nnscf)
      DO innscf=1, nnscf
         wtkn(innscf)=wt1
      ENDDO
      WRITE(IWR,*)thisroutine//': nnscf = ',nnscf
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LOOKUP_SYMPT
c
c
      INTEGER FUNCTION LOOKUP_SYMPT( sympt,symlabels,nsyms )
c
c     LOOKUP_SYMPT looks up the index of a special point symbol in a
c     character array.
c
      INTEGER nsyms
      CHARACTER sympt
      CHARACTER symlabels(*)
      INTEGER ispec, iresult
c
      iresult = -1
      DO ispec=1,nsyms
        IF(sympt.EQ.symlabels(ispec)) THEN
          iresult=ispec
        ENDIF
      END DO
c
      LOOKUP_SYMPT=iresult
c
      RETURN
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000
characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
