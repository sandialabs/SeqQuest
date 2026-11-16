      PROGRAM RINGS
c***************************************************************
c
c Program: RINGS
c
c Purpose: compute ring statistics for amorphous cell
c
c Written: Peter A. Schultz, 20-June-1996
c
c Last revised: 21Oct96-PAS / alter a0, fix r(2nn) breakdown
c
c***************************************************************
c
c Input files:  <geomfil> queried of user
c  1) filename: <geomfil>
c     unit number: IDAT=7
c     contents: atom coordinates - see routine GETGEOM
c Output files:  <outfil> name queried of user
c  1) filename: <outfil>.range
c     unit number: IRNG=8
c     contents: nearest neighbor bond list
c  2) filename: <outfil>.rings
c     unit number: IRNG=8 (same as above, explicitly open/close)
c     contents: detailed ring/range analysis
c  3) filename: <outfil>.ij4
c     unit number: IDIJ=9
c     contents: for 64 atom cell, puts out matrix of nn range
c
c Notes:
c
c 1Jul96-PAS:
c This program computes ring statistics for amorphous unit cells.
c Criterion defining ring: there is no shorter path between two
c atoms in the ring than a path along the ring.  This is exactly
c the Franzblau criterion [D.S. Franzblau, PRB 44, 4925 (1991)].
c I assume the shortest path connecting any atom to its image in
c another unit cell is 6nn.  If less to an image, the routines
c making 4nn sets and 6-member plus rings fail (this is not hard
c to fix).  Program will bail out if it finds this condition.
c In the dimensioning I assume that there will be at most a
c single 5-fold atom, or at least that 5-fold atoms are not very
c close to one another.  Program checks that its dimensioning
c (of the later nn* arrays) is adequate in the MKNN* routines
c and let you know if there is a problem.
c
      PARAMETER  (maxatm=225)
c
      DIMENSION  ratm(3,maxatm)
c
      DIMENSION  ijdist(maxatm*maxatm)
c
c  I will not allow more than 5nn for a full rings analysis:
      PARAMETER  (maxnnx=5)
      DIMENSION  num1(maxatm), nn1(5*maxatm), nn1lat(3*5*maxatm)
      DIMENSION  num2(maxatm), nn2(3*10*maxatm)
      DIMENSION  num3(maxatm), nn3(7*36*maxatm)
      DIMENSION  num4(maxatm), nn4(9*108*maxatm)
c
      CHARACTER*40  geomfil,outname,rangfl,ringfl,dijfil
      DATA  IRD,IWR,IDAT,IRNG,IDIJ / 5,6,7,8,9 /
c
c >>>> EXECUTABLE CODE:
c
      write(IWR,9000) maxatm
 9000 format(1x,72('=')
     $     /15x,'Program RINGS: ring statistics program'
     $     /20x,'Written by P.A. Schultz - 20Jun96'
     $     /20x,'Last Revised: 12Oct96-PAS'
     $     / 1x,72('=')
     $    // 5x,'Maximum number of atoms = ',i4/ )
c
      write(IWR,*) 'Enter geometry file name:'
      read(IRD,8000)  geomfil
 8000 format(a)
      if( geomfil.eq.' ' ) geomfil = 'atoms.in'
c
c Read in coordinates of atoms in system:
c
      OPEN(IDAT,file=geomfil,form='formatted',status='old',err=1301)
c
      call GETGEOM( IWR,IDAT,ierr, maxatm,natom,rcube,ratm )
c
      CLOSE(unit=IDAT)
      if( ierr.ne.0 ) goto 1302
      write(IWR,*) 'Geometry successfully read'
      write(IWR,*)
c
c Take care of some bookkeeping for output: filenames, etc.
c
   10 write(IWR,*) 'Enter output name:'
      read(IRD,8000)  outname
      if( outname.eq.' ' ) outname = 'rings'
      do 20 i=1,36
        lename = i - 1
        if( outname(i:i) .eq. ' ' ) goto 21
   20 continue
   21 continue
      if( lename.eq.0 ) goto 10 
      rangfl = ' '
      rangfl(1:lename) = outname(1:lename)
      rangfl(lename+1:lename+6) = '.range'
      ringfl = ' '
      ringfl(1:lename) = outname(1:lename)
      ringfl(lename+1:lename+6) = '.rings'
      dijfil = ' '
      dijfil(1:lename) = outname(1:lename)
      dijfil(lename+1:lename+6) = '.ij4'
c
c Ask for distance criterion for bond
c
      write(IWR,*) 'Enter bond range cutoff (1.80A):'
      read(IRD,*) bond
      if( bond .lt. 0.10 .or. bond .gt. 4.00 ) bond = 1.80
c
c Make nearest neighbor list:
c
      OPEN(unit=IRNG,file=rangfl,form='formatted',status='unknown')
c
      call MK1NN( IWR,IRNG, bond, rcube, natom, ratm, ijdist,
     $ maxnnx, mxnn, num1,nn1, nn1lat )
c
      CLOSE(unit=IRNG)
      if( mxnn.gt.maxnnx ) STOP 'mxnn.gt.5'
c
c Open output files:
c
      OPEN(unit=IRNG,file=ringfl,form='formatted',status='unknown')
      OPEN(unit=IDIJ,file=dijfil,form='formatted',status='unknown')
c
c Find 3-member rings:
c
      call COUNT3( IWR,IRNG, n3ring, natom, ijdist, num1,nn1 )
c
      write(IWR ,*) 'Number of 3-member rings found = ',n3ring
      write(IRNG,*) 'Number of 3-member rings found = ',n3ring
c
c Make second nearest neighbor list:
c
      call MK2NN( IWR,IRNG, natom, ijdist, num1,nn1, num2,nn2 )
c
c Find 4-member rings:
c
      call COUNT4( IWR,IRNG, n4ring, natom, ijdist,
     $ num1,nn1, num2,nn2 )
c
      write(IWR ,*) 'Number of 4-member rings found = ',n4ring
      write(IRNG,*) 'Number of 4-member rings found = ',n4ring
c
c Find 5-member rings:
c
      call COUNT5( IWR,IRNG, n5ring, n5chain, natom, ijdist,
     $ num1,nn1,nn1lat, num2,nn2)
c
      write(IWR ,*) 'Number of 5-member rings found = ',n5ring
      write(IRNG,*) 'Number of 5-member rings found = ',n5ring
c
      if( n5chain.gt.0 )then
c       Have a problem: 5nn connects an atom to its image in
c       another unit cell.  Later routines depend on this.
        write(IWR ,*) 'Number of 5-member 1-D chains = ',n5chain
        write(IRNG,*) 'Number of 5-member 1-D chains = ',n5chain
        STOP '5nn1d'
      endif
c
c Make 3rd nn list:
c
      call MKNN3( IWR,IRNG, numnn3, natom, ijdist, 
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3 )
c
C      write(IWR ,*) 'Number of 3nn sets =',numnn3
      write(IRNG,*) 'Number of 3nn sets =',numnn3
c
c Find 6-member rings:
c
      call COUNT6( IWR,IRNG, n6ring, n6chain, natom, ijdist, 
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3 )
c
      write(IWR ,*) 'Number of 6-member rings found = ',n6ring
      write(IRNG,*) 'Number of 6-member rings found = ',n6ring
c
C      write(IWR ,*) 'Number of 6-member 1-D chains = ',n6chain
      write(IRNG,*) 'Number of 6-member 1-D chains = ',n6chain
c
c Find 7-member rings:
c
      call COUNT7( IWR,IRNG, n7ring, n7chain, natom, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3 )
c
      write(IWR ,*) 'Number of 7-member rings found = ',n7ring
      write(IRNG,*) 'Number of 7-member rings found = ',n7ring
c
C      write(IWR ,*) 'Number of 7-member 1-D chains = ',n7chain
      write(IRNG,*) 'Number of 7-member 1-D chains = ',n7chain
CCc
CCc Find 8-member rings (from 3nn): (this routine replaced)
CCc
CC      call MK8NN3( IWR,IRNG, n8ring, n8chain, natom, ijdist,
CC     $ num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
CCc
CC      write(IWR ,*) 'Number of 8-member rings found = ',n8ring
CC      write(IRNG,*) 'Number of 8-member rings found = ',n8ring
CCc
CCC      write(IWR ,*) 'Number of 8-member 1-D chains = ',n8chain
CC      write(IRNG,*) 'Number of 8-member 1-D chains = ',n8chain
c
c Make 4nn sets:
c
      call MKNN4( IWR,IRNG,IDIJ, numnn4, natom, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
c
C      write(IWR ,*) 'Number of 4nn sets =',numnn4
      write(IRNG,*) 'Number of 4nn sets =',numnn4
c
c Find 8-member rings:
c
      call COUNT8( IWR,IRNG, n8ring, n8chain, natom, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
c
      write(IWR ,*) 'Number of 8-member rings found = ',n8ring
      write(IRNG,*) 'Number of 8-member rings found = ',n8ring
c
C      write(IWR ,*) 'Number of 8-member 1-D chains = ',n8chain
      write(IRNG,*) 'Number of 8-member 1-D chains = ',n8chain
c
c Find 9-member rings:
c
      call COUNT9( IWR,IRNG, n9ring, n9chain, natom, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
c
      write(IWR ,*) 'Number of 9-member rings found = ',n9ring
      write(IRNG,*) 'Number of 9-member rings found = ',n9ring
c
C      write(IWR ,*) 'Number of 9-member 1-D chains = ',n9chain
      write(IRNG,*) 'Number of 9-member 1-D chains = ',n9chain
c
c Do some post-analysis:
c
      call AUTOPSY( IWR,IRNG, rcube, natom, ratm, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3 )
c
      CLOSE(unit=IRNG)
      CLOSE(unit=IDIJ)
c
      STOP 'all done'
c
 1301 write(IWR,*)  '***** ERROR: no such file'
      STOP 'tryagain'
 1302 write(IWR,*)  '***** ERROR: reading geometry file'
      STOP 'tryagain'
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GETGEOM
c
c
      SUBROUTINE GETGEOM( IWR,IGIN,ierr, maxatm, natom, rcube,ratm )
c
      DIMENSION  ratm(3,*)
      CHARACTER  wantlbl*6,readlbl*6,line*72
      EQUIVALENCE  (readlbl,line)
      DATA  Aperau  / 0.5291771d0 /
c parameters for atC cells:
C      DATA  r64 / 14.3052324648 /,r216 / 10.874762 /
c parameters for Si cells:
      DATA  r64 / 20.40d0 /,r216 / 30.60d0 /
c
 8000 format(a)
      ierr = 0
      read(IGIN,8000)  line
      write(IWR,8000)  line
      if( readlbl .ne. 'iterat' ) goto 1310
      read(IGIN,*) iter
c
      iscale = 0
      scale = 1.d0
      read(IGIN,8000)  line
      write(IWR,8000)  line
      wantlbl = 'number'
      if( readlbl .ne. wantlbl )then
        wantlbl = 'center'
        if( readlbl .ne. wantlbl ) goto 1310
c        I have a QUEST format geometry
        iform = 1
        read(IGIN,8000)  line
        read(line,*,err=20)   natom,lenlbl,rcube
        goto 30
   20   read(line,*)  natom,lenlbl
        if( natom.eq.64 )then
          rcube = r64
        elseif( natom.eq.216 )then
          rcube = r216
        else
          ierr = 3
          write(IWR,*) '***** ERROR: need a lattice parameter'
          RETURN
        endif
   30   continue
        write(IWR,*)  'quest: number of atoms=',natom
        read(IGIN,8000)  line
        wantlbl = 'scale'
        if( readlbl .eq. wantlbl )then
          iscale = 1
          read( line(6:72) ,*) scale
        else
          BACKSPACE(unit=IGIN)
        endif
      else
c       FEIBELMAN format
        iform = 2
        read(IGIN,'(i2)') natom
        lenlbl = 8
        read(IGIN,8000)  line
        write(IWR,*) 'serial: number of atoms=',natom
      endif
c
c get cube dim and scale, in Angstrom
      rcube = rcube*Aperau
      scale = scale*Aperau
c
      if( natom.gt.maxatm )then
        write(IWR,*) '***** ERROR: too many atoms'
        write(IWR,*) '>>>>> have ',natom,', maxatm=',maxatm
        STOP 'max atom'
      endif
c
      write(*,'(a,2f12.8)') 'rcube,scale=',rcube,scale
      lrd = lenlbl + 1
      do 100 iatm=1,natom
        read(IGIN,8000)  line
        read(line(lrd:72),*)  (ratm(i,iatm),i=1,3)
        do 90 i=1,3
          ratm(i,iatm) = scale*ratm(i,iatm)
   90   continue
  100 continue
c
      RETURN
c
 1310 write(IWR,*)  '***** ERROR: looking for label:',wantlbl
      write(IWR,*)  '                         found:',readlbl
      ierr = 1
      RETURN
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MK1NN
c
c
      SUBROUTINE MK1NN( IWR,ILST, bond,rcube, natom, ratm, ijdist,
     $ maxnnx, mxnn, num1,nn1, nn1lat )
c---------------------------------------------------------------
c Purpose: this routine sets up everything that follows.  It
c takes the input coordinates and sets up the connectivity,
c taking into account whether the bond crosses unit cell bounds.
c
c Input:
c   bond  - bond distance criterion, smaller: yes, larger: no
c   rcube - unit cell assumed a cube, this is length of cube edge
c   natom - number of atoms
c   ratm(3,natom) - coordinates of each atom
c Output:
c   num1(natom)      - number of n.n. for each atom
c   nn1(maxnn,natom) - list of n.n. atoms
c   nn1lat(3,maxnn,natom) - tracks crossing of unit cell bounds
c     for bond - needed to distinguish rings from 1-D chains
c     within the algorithm I use to create rings
c---------------------------------------------------------------
c 
      PARAMETER  (maxnn=12)
      DIMENSION  ratm(3,natom)
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(maxnnx,natom), nn1lat(3,maxnnx,natom)
c local:
      DIMENSION  nfold(maxnn), jnn(maxnn),rnn(maxnn),ll(3),lnn(3,maxnn)
      DATA  zero / 0.d0 /
c
c initialize nn lists
c
      ierr = 0
      mxnn = 0
      do 10 iatm=1,natom
        do 8 i=1,5
          nn1(i,iatm) = 0
          nn1lat(1,i,iatm) = 0
          nn1lat(2,i,iatm) = 0
          nn1lat(3,i,iatm) = 0
    8   continue
        do 9 jatm=1,natom
          ijdist(iatm,jatm) = natom
    9   continue
        ijdist(iatm,iatm) = 0
   10 continue
c
      do 20 i=1,maxnn
        nfold(i) = 0
   20 continue
c
      bndave = zero
      bnd2ave = zero
      bnd3ave = zero
      bnd4ave = zero
      bndsqr = zero
      bndrms = zero
      bnderr = zero
c
      write(ILST,*)  ' Atom  bonds   nn.  Rnn'
c
      nbond = 0
      bave = zero
      brms = zero
c
      do 200 iatm=1,natom
        nn = 0
        do 100 jatm=1,natom
          if( jatm.eq.iatm ) goto 100
c
          call RNG2AT( rng,ll, rcube, ratm(1,iatm),ratm(1,jatm) )
          if( rng .le. bond )then
            nbond = nbond + 1
            bave = bave + rng
            brms = brms + rng*rng
            nn = nn + 1
            if( nn.gt.maxnn ) goto 100
            jnn(nn) = jatm
            rnn(nn) = rng
            lnn(1,nn) = ll(1)
            lnn(2,nn) = ll(2)
            lnn(3,nn) = ll(3)
          endif
  100   continue
c
c load up nn info
c
        if( nn.gt.mxnn ) mxnn = nn
        if( nn.le.maxnn )then
          nfold(nn) = nfold(nn) + 1
        else
          ierr = MAX( ierr, 2 )
        endif
        num1(iatm) = nn
        if( nn.gt.maxnnx )then
          ierr = MAX( ierr, 1 )
        else
          do 120 i=1,nn
            inn = jnn(i)
            nn1(i,iatm) = inn
            nn1lat(1,i,iatm) = lnn(1,i)
            nn1lat(2,i,iatm) = lnn(2,i)
            nn1lat(3,i,iatm) = lnn(3,i)
            ijdist(iatm,inn) = 1
            ijdist(inn,iatm) = 1
  120     continue
        endif
c
c compute trivial bond order analysis (carbon-based)
c
        bndord = 0.d0
        bnd2ord = 0.d0
        bnd3ord = 0.d0
        bnd4ord = 0.d0
        if( nn.gt.maxnn ) goto 180
        do 150 i=1,nn
          if( rnn(i) .lt. 1.75 )then
            bndord = bndord + (1.75-rnn(i))/.20
            if( rnn(i) .lt. 1.74 )then
              bnd2ord = bnd2ord + (1.74-rnn(i))/.20
              if( rnn(i) .lt. 1.73 )then
                bnd3ord = bnd3ord + (1.73-rnn(i))/.19
                bnd4ord = bnd4ord + (1.73-rnn(i))/.20
              endif
            endif
          endif
  150   continue
c
        bndave = bndave + bndord
        bnd2ave = bnd2ave + bnd2ord
        bnd3ave = bnd3ave + bnd3ord
        bnd4ave = bnd4ave + bnd4ord
        bndsqr = bndsqr + bndord**2
        bnderr = bnderr + ABS(bndord-4.d0)
        bndrms = bndrms + (bndord-4.d0)**2
c
  180   write(ILST,9100)  iatm,bndord,(jnn(i),rnn(i),i=1,nn)
 9100   format(1x,i3,f7.2, 6(i6,f7.3) )
c
  200 continue
c
      bndave = bndave/DBLE(natom)
      bnd2ave = bnd2ave/DBLE(natom)
      bnd3ave = bnd3ave/DBLE(natom)
      bnd4ave = bnd4ave/DBLE(natom)
      bndsqr = bndsqr/DBLE(natom)
      bnddev = SQRT( bndsqr - bndave**2)
      bnderr = bnderr/DBLE(natom)
      bndrms = SQRT(bndrms/DBLE(natom))
c
      bave = bave/DBLE(nbond)
      brms = brms/DBLE(nbond)
      bdev = SQRT( brms - bave*bave )
      brms = SQRT( brms )
      do 400 i=1,maxnn
        if( nfold(i).ne.0 )then
          write(IWR ,'(6(i3,"-fold atoms:",i4))') i,nfold(i)
          write(ILST,'(6(i3,"-fold atoms:",i4))') i,nfold(i)
        endif
  400 continue
c
      write(IWR ,9500) bndave,bnddev, bnderr,bndrms
      write(ILST,9500) bndave,bnddev, bnderr,bndrms
      write(IWR ,9501) bnd2ave,bnd3ave,bnd4ave
      write(ILST,9501) bnd2ave,bnd3ave,bnd4ave
 9501 format(1x,'ord[1.54,1.34]=', f7.4, ', ord[1.54,1.35]=',f7.4
     $        ', ord[1.53,1.33]=', f7.4 )
 9500 format(1x,'bond ord=',f7.4,' +/-',f7.4,' err,rms=',2f7.4)
CC      write(IWR ,9600) nbond/2, bave, bdev, brms
      write(ILST,9600) nbond/2, bave, bdev, brms
 9600 format(1x,'Number of bonds =',i5,
     $      /1x,'Ave bond length=',f7.4,'+/-',f6.4,', rms=',f7.4)
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RNG2AT
c
c
      SUBROUTINE RNG2AT( rij,ll, rcube, ri,rj )
c
      DIMENSION  ll(3),ri(3),rj(3), rt(3)
      DATA  half / 0.5d0 /
c
      do 10 i=1,3
        rt(i) = ri(i) - rj(i)
   10 continue
      do 20 i=1,3
        ll(i) = 0
        if( rt(i) .gt. half*rcube )then
          ll(i) = -1
          rt(i) = rt(i) - rcube
        elseif( rt(i) .lt. -half*rcube )then
          ll(i) =  1
          rt(i) = rt(i) + rcube
        endif
   20 continue
      rij = SQRT( rt(1)**2 + rt(2)**2 + rt(3)**2 )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> AUTOPSY
c
c
      SUBROUTINE AUTOPSY( IWR,ILST, rcube, natom, ratm, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3 )
c
      DIMENSION  ratm(3,natom)
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom), nn1lat(3,5,natom)
      DIMENSION  num2(natom), nn2(3,10,natom)
      DIMENSION  num3(natom), nn3(7,36,natom)
c local:
      DIMENSION  nxnn(10),r1sum(10),r2sum(10)
      DIMENSION  rmin(10),rmax(10),ijmin(2,10),ijmax(2,10)
      DIMENSION  ll(3)
      DATA  zero / 0.d0 /
c
      call RNG2AT( rng,ll, rcube, ratm(1,iatm),ratm(1,jatm) )
c
      n3f = 0
      n4f = 0
      nxf = 0
      n33 = 0
      n34 = 0
      n44 = 0
      r33 = zero
      r34 = zero
      r44 = zero
      r33sq = zero
      r34sq = zero
      r44sq = zero
      do 10 nndist=1,10
        nxnn(nndist) = 0
        r1sum(nndist) = zero
        r2sum(nndist) = zero
        rmin(nndist) = 99.d0
        rmax(nndist) = zero
   10 continue
c do analysis thru nndmax'th n.n.
      nndmax = 3
c
c accumulate nn distances
c
      do 100 jatm=1,natom
        nf = num1(jatm)
        if( nf.eq.3 )then
          n3f = n3f + 1
        elseif( nf.eq.4 )then
          n4f = n4f + 1
        else
          nxf = nxf + 1
        endif
        do 90 iatm=1,natom
          if( iatm.eq.jatm ) goto 90
          call RNG2AT( rng,ll, rcube, ratm(1,iatm),ratm(1,jatm) )
          ijd = ijdist(iatm,jatm)
          if( ijd.lt.1 .or. ijd.gt.nndmax )then
            nndist = nndmax + 1
            nxnn(nndist) = nxnn(nndist) + 1
            if( rng.lt.rmin(nndist) )then
              rmin(nndist) = rng
              ijmin(1,nndist) = iatm
              ijmin(2,nndist) = jatm
            elseif( rng.gt.rmax(nndist) )then
              rmax(nndist) = rng
              ijmax(1,nndist) = iatm
              ijmax(2,nndist) = jatm
            endif
            goto 90
          endif
          do 80 nndist=1,nndmax
            if( ijd.eq.nndist )then
              nxnn(nndist) = nxnn(nndist) + 1
              r1sum(nndist) = r1sum(nndist) + rng
              r2sum(nndist) = r2sum(nndist) + rng*rng
              if( nndist.eq.1 )then
                if( num1(iatm)+num1(jatm) .eq. 8 )then
                  n44 = n44 + 1
                  r44 = r44 + rng
                  r44sq = r44sq + rng*rng
                elseif( num1(iatm)+num1(jatm) .eq. 7 )then
                  n34 = n34 + 1
                  r34 = r34 + rng
                  r34sq = r34sq + rng*rng
                elseif( num1(iatm).eq.3 .and. num1(jatm).eq.3 )then
                  n33 = n33 + 1
                  r33 = r33 + rng
                  r33sq = r33sq + rng*rng
                endif
              endif
              if( rng.lt.rmin(nndist) )then
                rmin(nndist) = rng
                ijmin(1,nndist) = iatm
                ijmin(2,nndist) = jatm
              endif
              if( rng.gt.rmax(nndist) )then
                rmax(nndist) = rng
                ijmax(1,nndist) = iatm
                ijmax(2,nndist) = jatm
              endif
            endif
   80     continue
   90   continue
  100 continue
c
c Do breakdown
c
      dat = DBLE(natom)
      do 200 nndist=1,nndmax
        numxnn = nxnn(nndist)
        dnum = DBLE(numxnn)
        dave = r1sum(nndist)/dnum
        d2av = r2sum(nndist)/dnum
        ddev = SQRT( d2av - dave*dave )
        drms = SQRT( d2av )
        numxnn = numxnn/2
        dnum = dnum/dat
        dmin = rmin(nndist)
        dmax = rmax(nndist)
        write(IWR ,9200) nndist,numxnn,dnum,dave,ddev,drms,dmin,dmax
        write(ILST,9200) nndist,numxnn,dnum,dave,ddev,drms,dmin,dmax
 9200   format(1x,'nn=',i2,': n=',i4,', n/at=',f7.4,
     $   ', r=',f7.4,' +/-',f7.4,', rms=',f7.4,', [',f6.4,',',f6.4,']')
 9201   format(1x,' else' ,': n=',i4,'       ', 7x ,
     $   '    ', 7x ,'    ', 7x ,'      ', 7x ,', [',f6.4,',',f6.4,']')
        imn = ijmin(1,nndist)
        jmn = ijmin(2,nndist)
        imx = ijmax(1,nndist)
        jmx = ijmax(2,nndist)
 9211   format(5x,'imin=',i3,1x,3f8.4,5x,'imax=',i3,1x,3f8.4)
 9212   format(5x,'jmin=',i3,1x,3f8.4,5x,'jmax=',i3,1x,3f8.4)
CC        write(IWR ,9211) imn,(ratm(i,imn),i=1,3),imx,(ratm(i,imx),i=1,3) 
CC        write(IWR ,9212) jmn,(ratm(i,jmn),i=1,3),jmx,(ratm(i,jmx),i=1,3) 
        write(ILST,9211) imn,(ratm(i,imn),i=1,3),imx,(ratm(i,imx),i=1,3) 
        write(ILST,9212) jmn,(ratm(i,jmn),i=1,3),jmx,(ratm(i,jmx),i=1,3) 
  200 continue
c
      nndist = nndmax + 1
      if( nxnn(nndist).gt.0 )then
        numxnn = nxnn(nndist)
        dmin = rmin(nndist)
        dmax = rmax(nndist)
        write(IWR ,9201)        numxnn,                    dmin,dmax
        imn = ijmin(1,nndist)
        jmn = ijmin(2,nndist)
        imx = ijmax(1,nndist)
        jmx = ijmax(2,nndist)
CC        write(IWR ,9211) imn,(ratm(i,imn),i=1,3),imx,(ratm(i,imx),i=1,3) 
CC        write(IWR ,9212) jmn,(ratm(i,jmn),i=1,3),jmx,(ratm(i,jmx),i=1,3) 
        write(ILST,9211) imn,(ratm(i,imn),i=1,3),imx,(ratm(i,imx),i=1,3) 
        write(ILST,9212) jmn,(ratm(i,jmn),i=1,3),jmx,(ratm(i,jmx),i=1,3) 
      endif
      write(IWR ,*) 'Note that 3nn is not clean if have 6nn 1d chains'
      write(ILST,*) 'Note that 3nn is not clean if have 6nn 1d chains'
c
c split 1nn analysis - 3-3 bond, 3-4 bond, 4-4 bond
c
      write(IWR ,*) 'Number of 3-fold, 4-fold, else atoms=',n3f,n4f,nxf
      write(IWR ,*) 'nn breakdown by types of atoms bonded:'
      write(ILST,*) 'Number of 3-fold, 4-fold, else atoms=',n3f,n4f,nxf
      write(ILST,*) 'nn breakdown by types of atoms bonded:'
 9300 format(5x,a,i4,', N1/at=',f7.4,', Rave=',f7.4,' +/-',f7.4)
      if( n33.gt.0 )then
        dnum = DBLE(n33)
        dave = r33/dnum
        d2av = r33sq/dnum
        ddev = SQRT( d2av - dave*dave )
        drms = SQRT( d2av )
        dnum = dnum/dat
        write(IWR ,9300) '[3-3]',n33/2,dnum,dave,ddev
        write(ILST,9300) '[3-3]',n33/2,dnum,dave,ddev
      endif
      if( n34.gt.0 )then
        dnum = DBLE(n34)
        dave = r34/dnum
        d2av = r34sq/dnum
        ddev = SQRT( d2av - dave*dave )
        drms = SQRT( d2av )
        dnum = dnum/dat
        write(IWR ,9300) '[3-4]',n34/2,dnum,dave,ddev
        write(ILST,9300) '[3-4]',n34/2,dnum,dave,ddev
      endif
      if( n44.gt.0 )then
        dnum = DBLE(n44)
        dave = r44/dnum
        d2av = r44sq/dnum
        ddev = SQRT( d2av - dave*dave )
        drms = SQRT( d2av )
        dnum = dnum/dat
        write(IWR ,9300) '[4-4]',n44/2,dnum,dave,ddev
        write(ILST,9300) '[4-4]',n44/2,dnum,dave,ddev
      endif
      if( (n34+n44).gt.0 )then
        dnum = DBLE(n44+n34)
        dave = (r34+r44)/dnum
        d2av = (r34sq+r44sq)/dnum
        ddev = SQRT( d2av - dave*dave )
        drms = SQRT( d2av )
        dnum = dnum/dat
        write(IWR ,9300) '[3,4-4]',(n34+n44)/2,dnum,dave,ddev
        write(ILST,9300) '[3,4-4]',(n34+n44)/2,dnum,dave,ddev
      endif
c
c split 2nn analysis - in 4-ring/5-ring/else
c
      do i=1,10
        nxnn(i) = 0
        r1sum(i) = zero
        r2sum(i) = zero
      enddo
c
      n6pl = 0
      r1pl = zero
      r2pl = zero
      numnn2 = 0
c
      do 400 i2=1,natom
        do 380 i2nn=1,num2(i2)
          numnn2 = numnn2 + 1
          myring = nn2(2,i2nn,i2)
          if( myring.eq.3 )then
c           Skip
            n3use = n3use + 1
            goto 380
          elseif( myring.eq.0 )then
            myring = 10
            notuse = notuse + 1
          endif
          iatm = nn2(1,i2nn,i2)
          jatm = nn2(3,i2nn,i2)
          call RNG2AT( rng,ll, rcube, ratm(1,iatm),ratm(1,jatm) )
          nxnn(myring) = nxnn(myring) + 1
          r1sum(myring) = r1sum(myring) + rng
          r2sum(myring) = r2sum(myring) + rng*rng
          if( myring.ge.6 )then
            n6pl = n6pl + 1
            r1pl = r1pl + rng
            r2pl = r2pl + rng*rng
          endif
  380   continue
  400 continue
c
c break down 2nn analysis
c
      write(IWR ,*)  '2nn analysis by smallest ring it is in:'
      write(ILST,*)  '2nn analysis by smallest ring it is in:'
      do 500 myring=4,10
        nx = nxnn(myring)
        if( nx.eq.0 ) goto 500
        dnum = DBLE(nx)
        dave = r1sum(myring)/dnum
        d2av = r2sum(myring)/dnum
        ddev = SQRT( d2av - dave*dave )
        drms = SQRT( d2av )
        dnum = dnum/dat
c fix count (4-rings double-count 2nn)
        if( myring.eq.4 )then
          nx = nx/2
        else
          dnum = dnum + dnum
        endif
        write(IWR ,9500)  myring,nx,dnum,dave,ddev,drms
        write(ILST,9500)  myring,nx,dnum,dave,ddev,drms
 9500   format(1x,'ring=',i2,': n=',i4,', n/at=',f7.4,
     $   ', r=',f7.4,' +/-',f7.4,', rms=',f7.4)
  500 continue
c
c combine 2nn in 6+ rings
c
      nx = n6pl
      dnum = DBLE(nx)
      dave = r1pl/dnum
      d2av = r2pl/dnum
      ddev = SQRT( d2av - dave*dave )
      drms = SQRT( d2av )
      dnum = (dnum+dnum)/dat
      write(IWR ,9506)  nx,dnum,dave,ddev,drms
      write(ILST,9506)  nx,dnum,dave,ddev,drms
 9506 format(1x,'6+ rings: n=',i4,', n/at=',f7.4,
     $ ', r=',f7.4,' +/-',f7.4,', rms=',f7.4)
c
      write(IWR ,*) 'Number of unused nn2 sets=',notuse,' of',numnn2
      write(ILST,*) 'Number of unused nn2 sets=',notuse,' of',numnn2
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKZERI
c
c
      SUBROUTINE MKZERI( n, iarr )
      DIMENSION  iarr(n)
      do 10 i=1,n
        iarr(i) = 0
   10 continue
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LOC2NN
c
c
      SUBROUTINE LOC2NN( i2nn, i1,i2,i3, natom,num2,nn2 )
c---------------------------------------------------------------
c Purpose: find index of 2nn string in list
c---------------------------------------------------------------
      DIMENSION  num2(natom), nn2(3,10,natom)
c
c orient string correctly
      if( i1.lt.i3 )then
        j1 = i1
        j3 = i3
      else
        j1 = i3
        j3 = i1
      endif
c
c search through 2nn list
      do 10 j2nn=1,num2(i2)
        if( nn2(1,j2nn,i2).eq.j1 .and.
     $      nn2(3,j2nn,i2).eq.j3 ) then
          i2nn = j2nn
          RETURN
        endif
   10 continue
c
      i2nn = -1
      STOP 'bad 2nn'
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MK2NN
c
c
      SUBROUTINE MK2NN( IWR,ILST, natom, ijdist, num1,nn1, num2,nn2 )
c---------------------------------------------------------------
c Purpose: create list of sets of non-duplicate 2nn sets
c---------------------------------------------------------------
c Note: I record "2nn" sets that result in 1nn (i.e. 3-rings)
c       I track them by recording ring-index value of 3
c
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom)
      DIMENSION  num2(natom), nn2(3,10,natom)
c
c >>>> EXECUTABLE CODE:
c
      nn2tot = 0
      do 100 i2=1,natom
c
c clear out 2nn lists
        do 10 i2nn=1,10
          do 10 i=1,3
            nn2(i,i2nn,i2) = 0
   10   continue
c
        nn = num1(i2)
        num2(i2) = nn*(nn-1)/2
c
c load 2nn lists [low-end,ring-order,hi-end]
        i2ii = 0
        do 20 inn=1,nn-1
          do 20 jnn=inn+1,nn
            i1 = nn1(inn,i2)
            i3 = nn1(jnn,i2)
            if( i3.lt.i1 )then
              ix = i1
              i1 = i3
              i3 = ix
            endif
            i2ii = i2ii + 1
            nn2tot = nn2tot + 1
            nn2(1,i2ii,i2) = i1
            nn2(3,i2ii,i2) = i3
            if( ijdist(i1,i3).eq.1 )then
              nn2(2,i2ii,i2) = 3
            else
              ijdist(i1,i3) = 2
              ijdist(i3,i1) = 2
            endif
   20   continue
c
  100 continue
c
C      write(IWR ,*) 'Number of 2nn sets =',nn2tot
      write(ILST,*) 'Number of 2nn sets =',nn2tot
      n1nn = 0
      n2nn = 0
      numnn2 = 0
      do 200 iatm=1,natom
        n1nn = n1nn + num1(iatm)
        do 120 i2=1,num2(iatm)
          if( nn2(2,i2,iatm).ne.3 ) numnn2 = numnn2 + 1
  120   continue
        do 140 jatm=iatm,natom
          if( ijdist(iatm,jatm).eq.2 ) n2nn = n2nn + 1
  140   continue
  200 continue
      write(ILST,*) 'Number of clean 2nn sets/atoms =',numnn2,n2nn
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKNN3
c
c
      SUBROUTINE MKNN3( IWR,IRNG, numnn3, natom, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3 )
c---------------------------------------------------------------
c Purpose: create list of sets of non-duplicate true 3nn sets
c---------------------------------------------------------------
c
c Note: the string is guaranteed 3nn, but note that endpoints
c       may be less along different path crossing u.c.
c
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom), nn1lat(3,5,natom)
      DIMENSION  num2(natom), nn2(3,10,natom)
      DIMENSION  num3(natom), nn3(7,36,natom)
c local:
      LOGICAL  badrng
c
c >>>> EXECUTABLE CODE:
c
      maxnn3 = 36
      numnn3 = 0
      call MKZERI( natom, num3 )
      call MKZERI( 7*maxnn3*natom, nn3 )
c
      do 300 i2=1,natom
c
        do 200 i2ii=1,num2(i2)
c
          if( nn2(2,i2ii,i2).eq.3 ) goto 200
c
          iend1 = nn2(1,i2ii,i2)
          iend3 = nn2(3,i2ii,i2)
          i1 = iend1
          i3 = iend3
   10     continue
c
          if( i2.lt.i3 )then
            do 100 nni3=1,num1(i3)
              i4 = nn1(nni3,i3)
              if( i4.eq.i2 .or. ijdist(i1,i4).lt.2 
     $                     .or. ijdist(i2,i4).ne.2 ) goto 100
c
c just check in case 5 laps cell
              if( ijdist(i1,i4).eq.2 )then
                do 20 i4nn=1,num1(i4)
                  i5 = nn1(i4nn,i4)
                  if( ijdist(i5,i1).eq.1 )then
                    call CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $                           5, i1,i2,i3,i4,i5,i6,i7,i8,i9 )
                    if( badrng )then
                      goto 21
                    else
                      goto 100
                    endif
                  endif
   20           continue
                STOP '3NNchain'
   21           continue
              else
                if( ijdist(i1,i4).lt.3 ) STOP '3NN-3rng'
                ijdist(i1,i4) = 3
                ijdist(i4,i1) = 3
              endif
c
c These guys are true 3rd nn, [i1:i2:i3:i4], load it up
c
              call LOC2NN( i3ii, i2,i3,i4, natom, num2,nn2 )
              numnn3 = numnn3 + 1
              if( i4.lt.i1 )then
                num3(i4) = num3(i4) + 1
                nj = num3(i4)
                nn3(1,nj,i4) = i4
                nn3(2,nj,i4) = i3
                nn3(3,nj,i4) = i2
                nn3(4,nj,i4) = i1
                nn3(5,nj,i4) = 0
                nn3(6,nj,i4) = i3ii
                nn3(7,nj,i4) = i2ii
              else
                num3(i1) = num3(i1) + 1
                nj = num3(i1)
                nn3(1,nj,i1) = i1
                nn3(2,nj,i1) = i2
                nn3(3,nj,i1) = i3
                nn3(4,nj,i1) = i4
                nn3(5,nj,i1) = 0
                nn3(6,nj,i1) = i2ii
                nn3(7,nj,i1) = i3ii
              endif
  100       continue
          endif
c
          if( i1.ne.iend3 )then
            i1 = iend3
            i3 = iend1
            goto 10
          endif
c
  200   continue
c
  300 continue
c
      mxnum3 = num3(1)
      do 310 iatm=1,natom
        if( num3(iatm).gt.mxnum3 ) mxnum3 = num3(iatm)
  310 continue
      if( mxnum3.gt.maxnn3 )then
        write(IWR ,*) 'Too many nn3 =',mxnum3,', MAX =',maxnn3
        write(IRNG,*) 'Too many nn3 =',mxnum3,', MAX =',maxnn3
        STOP 'MKNN3'
      endif
CCcpas:start
CCc WRITE OUT VIABLE 3NN SETS
CC      OPEN(unit=13,file='tmp.nn3',form='formatted',status='unknown')
CC      write(13,*) 'nn3 sets='
CC      do 400 i1=1,natom
CC        n3s = num3(i1)
CC        do 400 i3s=1,n3s
CC          write(13,413) i3s,(nn3(i,i3s,i1),i=1,5)
CC  413     format(1x,"set#",i3," = ",4i3,i6)
CC  400 continue
CC      CLOSE(unit=13)
CCcpas:end
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKNN4
c
c
      SUBROUTINE MKNN4( IWR,IRNG,IDIJ, numnn4, natom, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
c---------------------------------------------------------------
c Purpose: create list of sets of non-duplicate true 4nn sets
c---------------------------------------------------------------
c
c Note: the string is guaranteed 4nn, but note that endpoints
c       may be less along different path crossing u.c.
c
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom), nn1lat(3,5,natom)
      DIMENSION  num2(natom), nn2(3,10,natom)
      DIMENSION  num3(natom), nn3(7,36,natom)
      DIMENSION  num4(natom), nn4(9,108,natom)
c local:
      LOGICAL    badrng
c
c >>>> EXECUTABLE CODE:
c
      maxnn4 = 108
      numnn4 = 0
      call MKZERI( natom, num4 )
      call MKZERI( 9*maxnn4*natom, nn4 )
c
      do 300 i1atm=1,natom
c
        do 200 i3s=1,num3(i1atm)
          i1 = i1atm
          i2 = nn3(2,i3s,i1atm)
          i3 = nn3(3,i3s,i1atm)
          i4 = nn3(4,i3s,i1atm)
c
   10     continue
c
          do 100 i4nn=1,num1(i4)
            i5 = nn1(i4nn,i4)
c will catch this from other 3-ring
            if( i5.lt.i1 ) goto 100
c
c check if backtrack/assume need at least 6nn to find an image
            if( ijdist(i5,i3).ne.2 .or. ijdist(i5,i2).ne.3 ) goto 100
c
            if( ijdist(i5,i1).eq.2 )then
              do 20 i5nn=1,num1(i5)
                j6 = nn1(i5nn,i5)
                if( ijdist(j6,i1).eq.1 )then
                  call CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $                         6, i1,i2,i3,i4,i5,j6,i7,i8,i9 )
                  if( .not. badrng ) goto 100
                  goto 21
                endif
   20         continue
              STOP 'MK4NN-1'
   21         continue
            elseif( ijdist(i5,i1).eq.3 )then
              do 45 i5nn=1,num1(i5)
                j6 = nn1(i5nn,i5)
                do 40 i1nn=1,num1(i1)
                  j7 = nn1(i1nn,i1)
                  if( ijdist(j6,j7).eq.1 )then
                    call CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $                           7, i1,i2,i3,i4,i5,j6,j7,i8,i9 )
                    if( .not.badrng ) goto 100
                  endif
   40           continue
   45         continue
            else
              ijdist(i1,i5) = 4
              ijdist(i5,i1) = 4
            endif
c
            numnn4 = numnn4 + 1
            num4(i1) = num4(i1) + 1
            nj = num4(i1)
            nn4(1,nj,i1) = i1
            nn4(2,nj,i1) = i2
            nn4(3,nj,i1) = i3
            nn4(4,nj,i1) = i4
            nn4(5,nj,i1) = i5
  100     continue
c
          if( i1.eq.i1atm )then
c           Reverse the ring and try extending again
            i1 = nn3(4,i3s,i1atm)
            i2 = nn3(3,i3s,i1atm)
            i3 = nn3(2,i3s,i1atm)
            i4 = i1atm
            goto 10
          endif
              
  200   continue
c
  300 continue
c
      mxnum4 = num4(1)
      do 310 iatm=1,natom
        if( num4(iatm).gt.mxnum4 ) mxnum4 = num4(iatm)
  310 continue
      if( mxnum4.gt.maxnn4 )then
        write(IWR ,*) 'Too many nn4 =',mxnum4,', MAX =',maxnn4
        write(IRNG,*) 'Too many nn4 =',mxnum4,', MAX =',maxnn4
        STOP 'MKNN4'
      endif
c
      if( natom.eq.64 )then
c       Write out ijdist matrix
        write(IDIJ,*) 'ijdist matrix'
        write(IDIJ,1390)  ((i,i=1,9),j=1,6),(i,i=1,4)
 1390   format(5x,9(' '),'1 ',9('1'),'2 ',9('2'),'3 ',9('3'),'4 ',
     $            9('4'),'5 ',9('5'),'6 ',4('6')
     $        /4x,7(1x,9i1,'0') )
        write(IDIJ,*)
        do 1300 j=1,64
          write(IDIJ,1391)  j,(ijdist(i,j),i=1,64)
 1391     format(1x,i2,1x,7(1x,10i1))
          if( (j/10)*10.eq.j ) write(IDIJ,*)
 1300   continue
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CHKLAT
c
c
      SUBROUTINE CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $             nring, i1,i2,i3,i4,i5,i6,i7,i8,i9 )
c---------------------------------------------------------------
c Purpose: check if links return to same atom giving a ring, or
c          if links lead to atom image in different unit cell
c          giving a 1-D periodic chain.
c---------------------------------------------------------------
c
      DIMENSION  num1(natom), nn1(5,natom), nn1lat(3,5,natom)
c local:
      PARAMETER  (maxrng=9)
      DIMENSION  iring(maxrng+1),ll(3)
      LOGICAL    badrng
c
      ichain = 0
      if( nring.lt.3 ) STOP 'CHKLAT-d'
      if( nring.gt.maxrng ) STOP 'CHKLAT-D'
c
      iring(1) = i1
      iring(2) = i2
      iring(3) = i3
      if( nring.ge.4 ) iring(4) = i4
      if( nring.ge.5 ) iring(5) = i5
      if( nring.ge.6 ) iring(6) = i6
      if( nring.ge.7 ) iring(7) = i7
      if( nring.ge.8 ) iring(8) = i8
      if( nring.ge.9 ) iring(9) = i9
      if( nring.gt.9 ) STOP 'CHKLAT-d'
      iring(nring+1) = i1
      ll(1) = 0
      ll(2) = 0
      ll(3) = 0
      do 100 ir=1,nring
        iatm = iring(ir)
        jatm = iring(ir+1)
        do 20 inn=1,num1(iatm)
          if( nn1(inn,iatm).eq.jatm )then
            ll(1) = ll(1) + nn1lat(1,inn,iatm)
            ll(2) = ll(2) + nn1lat(2,inn,iatm)
            ll(3) = ll(3) + nn1lat(3,inn,iatm)
            goto 100
          endif
   20   continue
c we did not find iatm:jatm to be nn!!
        write(*,*) 'CHKLAT: bad nn=(',iatm,jatm,') in',
     $   (iring(is),is=1,nring)
        STOP 'CHKLATnn'
  100 continue
c
      badrng = .false.
      if( ll(1).ne.0 .or. ll(2).ne.0 .or. ll(3).ne.0 ) badrng = .true.
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COUNT3
c
c
      SUBROUTINE COUNT3( IWR,IRNG, n3ring, natom, ijdist, num1,nn1 )
c---------------------------------------------------------------
c Purpose: find and count 3-member rings
c---------------------------------------------------------------
c
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom)
c
 9030 format(1x,'3-Ring     = ',9i4)
 9031 format(1x,'3-Ring#',i4,' = ',9i4)
c
      n3ring = 0
c
      do 100 i1=1,natom
        nni = num1(i1)
        do 80 inn=1,nni-1
          i2 = nn1(inn,i1)
          if( i1.lt.i2 )then
            do 70 jnn=inn+1,nni
              i3 = nn1(jnn,i1)
              if( i1.lt.i3 )then
                if( ijdist(i2,i3).eq.1 )then
c                 Have a 3-Ring
                  n3ring = n3ring + 1
                  if( i2.lt.i3 )then
                    write(IRNG,9030)  i1,i2,i3
CNUM                    write(IRNG,9031)  n3ring, i1,i2,i3
                  else
                    write(IRNG,9030)  i1,i2,i3
CNUM                    write(IRNG,9030)  n3ring, i1,i3,i2
                  endif
                endif
              endif
   70       continue
          endif
   80   continue
c
  100 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COUNT4
c
c
      SUBROUTINE COUNT4( IWR,IRNG, n4ring, natom, ijdist,
     $ num1,nn1, num2,nn2 )
c---------------------------------------------------------------
c Purpose: find and count 4-member rings
c---------------------------------------------------------------
c
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom)
      DIMENSION  num2(natom), nn2(3,10,natom)
c
 9040 format(1x,'4-Ring     = ',9i4)
 9044 format(1x,'4-Ring#    ',i4,' = ',9i4)
 9045 format(1x,'4-Ringlap3#',i4,' = ',9i4)
c
      n4ring = 0
      n4lap3 = 0
c
      do 200 i2=1,natom
c loop over 2nn sets for this atom
        do 100 i2ii=1,num2(i2)
c get index of atoms at ends of this 2nn
          i1 = nn2(1,i2ii,i2)
          i3 = nn2(3,i2ii,i2)
          if( i1.lt.i2 .or. i3.lt.i2 .or. ijdist(i1,i3).lt.2 ) goto 100
c
          nni1 = num1(i1)
          nni3 = num1(i3)
          do 20 i1nn=1,nni1
            i4 = nn1(i1nn,i1)
            if( i4.le.i2 ) goto 20
            do 10 i3nn=1,nni3
              if( nn1(i3nn,i3).eq.i4 )then
                if( ijdist(i2,i4).eq.2 )then
                  n4ring = n4ring + 1
                  call LOC2NN( i1ii, i2,i1,i4, natom, num2,nn2 )
                  call LOC2NN( i3ii, i2,i3,i4, natom, num2,nn2 )
                  call LOC2NN( i4ii, i1,i4,i3, natom, num2,nn2 )
                  if( i3.lt.i1 )then
                    write(IRNG,9040)  i2,i3,i4,i1
CNUM                    write(IRNG,9044)  n4ring, i2,i3,i4,i1
                  else
                    write(IRNG,9040)  i2,i1,i4,i3
CNUM                    write(IRNG,9044)  n4ring, i2,i1,i4,i3
                  endif
                  nn2(2,i1ii,i1) = 4
                  nn2(2,i2ii,i2) = 4
                  nn2(2,i3ii,i3) = 4
                  nn2(2,i4ii,i4) = 4
                elseif( ijdist(i2,i4).eq.1 )then
                  n4lap3 = n4lap3 + 1
                  call LOC2NN( i4ii, i3,i4,i1, natom,num2,nn2 )
                  if( nn2(2,i2ii,i2).eq.0 ) nn2(2,i2ii,i2) = 4
                  if( nn2(2,i4ii,i4).eq.0 ) nn2(2,i4ii,i4) = 4
                  write(IRNG,9045)  n4lap3, i1,i2,i3,i4
                endif
              endif
   10       continue
   20     continue
c
  100   continue
c
  200 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COUNT5
c
c
      SUBROUTINE COUNT5( IWR,IRNG, n5ring, n5chain, natom, ijdist, 
     $ num1,nn1,nn1lat, num2,nn2)
c---------------------------------------------------------------
c Purpose: find and count 5-member rings
c---------------------------------------------------------------
c
c This routine finds 5-rings by starting at fixed atom and finding
c different 2nn sets whose endpoints are n.n.  This routine will
c check if 5-"rings" actually are 1-D periodic chains connecting
c an atom and its image in another unit cell.
c
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom), nn1lat(3,5,natom)
      DIMENSION  num2(natom), nn2(3,10,natom)
c local:
      LOGICAL    badrng
c
c >>>> EXECUTABLE CODE:
c
 9050 format(1x,'5-Ring     = ',9i4)
 9051 format(1x,'5-Ringlap4 = ',9i4)
 9150 format(1x,'5-Chain    ',i4,' = ',9i4)
 9055 format(1x,'5-Ring#    ',i4,' = ',9i4)
 9056 format(1x,'5-Ringlap4#',i4,' = ',9i4)
c
      n5ring = 0
c
      do 300 i2=1,natom
c
        do 200 i2ii=1,num2(i2)
c
          i1 = nn2(1,i2ii,i2)
          i3 = nn2(3,i2ii,i2)
          if( i1.lt.i2 .or. i3.lt.i2 .or. ijdist(i1,i3).ne.2 ) goto 200
c
          do 100 nni3=1,num1(i3)
            i4 = nn1(nni3,i3)
            if( i4.le.i2 .or. ijdist(i1,i4).eq.1 
     $                   .or. ijdist(i2,i4).ne.2 ) goto 100
c
            if( ijdist(i1,i4).ne.2 ) goto 100
c
c           We have a five - find missing atom(s)
c
            call LOC2NN( i3ii, i2,i3,i4, natom, num2,nn2)
            if( nn2(2,i2ii,i2).eq.0 ) nn2(2,i2ii,i2) = 5
            if( nn2(2,i3ii,i3).eq.0 ) nn2(2,i3ii,i3) = 5
c
            do 65 nni4=1,num1(i4)
              i5 = nn1(nni4,i4)
              if( i5.lt.i2 .or. ijdist(i5,i3).ne.2 
     $                     .or. ijdist(i5,i2).ne.2 ) goto 65
c
              do 60 nni1=1,num1(i1)
                if( nn1(nni1,i1).ne.i5 .or. ijdist(i5,i2).ne.2 ) goto 60
c
c Trace steps to be sure have true ring and not periodic chain
                call CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $                       5, i1,i2,i3,i4,i5,i6,i7,i8,i9 )
                if( badrng )then
                  n5chain = n5chain + 1
                  write(IRNG,9150)  n5chain, i1,i2,i3,i4,i5
                  goto 60
                endif
c
                n5ring = n5ring + 1
                call LOC2NN( i4ii, i3,i4,i5, natom, num2,nn2)
                call LOC2NN( i5ii, i4,i5,i1, natom, num2,nn2)
                call LOC2NN( i1ii, i5,i1,i2, natom, num2,nn2)
                if( nn2(2,i4ii,i4).eq.0 ) nn2(2,i4ii,i4) = 5
                if( nn2(2,i5ii,i5).eq.0 ) nn2(2,i5ii,i5) = 5
                if( nn2(2,i1ii,i1).eq.0 ) nn2(2,i1ii,i1) = 5
                irngs = nn2(2,i1ii,i1) + nn2(2,i2ii,i2) + nn2(2,i3ii,i3)
     $                + nn2(2,i4ii,i4) + nn2(2,i5ii,i5)
                if( (irngs/5)*5.eq.irngs )then
c                 Have a clean 5-ring
                  write(IRNG,9050)  i2,i1,i5,i4,i3
CNUM                  write(IRNG,9055)  n5ring, i2,i1,i5,i4,i3
                else
c                 Have a 2nn already used in a 4-ring in the set
                  write(IRNG,9051)  i2,i1,i5,i4,i3
CNUM                  write(IRNG,9056)  n5ring, i2,i1,i5,i4,i3
                endif
c
   60         continue
c
   65       continue
c
  100     continue
c
  200   continue
c
  300 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COUNT6
c
c
      SUBROUTINE COUNT6( IWR,IRNG, n6ring, n6chain, natom, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3 )
c---------------------------------------------------------------
c Purpose: find and count 6-member rings
c---------------------------------------------------------------
c
c This routine finds 6-rings by starting at fixed atom and finding
c different 3nn sets whose endpoints match.  This routine will
c check if 6-"rings" actually are 1-D periodic chains connecting
c an atom and its image in another unit cell.  Also checks that
c all diagonals are no less than three atoms apart, taking into
c account possibility of linking between rings across u.c. bound.
c
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom), nn1lat(3,5,natom)
      DIMENSION  num2(natom), nn2(3,10,natom)
      DIMENSION  num3(natom), nn3(7,36,natom)
c local:
      LOGICAL    badrng
c
c >>>> EXECUTABLE CODE:
c
 9060 format(1x,'6-Ring     = ',9i4)
 9061 format(1x,'6-Ringlaps = ',9i4)
 9160 format(1x,'6-Chain#   ',i4,' = ',9i4)
 9066 format(1x,'6-Ring#    ',i4,' = ',9i4)
 9067 format(1x,'6-Ringlaps#',i4,' = ',9i4)
c
      n6ring = 0
      n6chain = 0
c
      do 300 i1=1,natom
        n3s = num3(i1)
        if( n3s.le.1 ) goto 300
c
        do 200 j3s=1,n3s-1
          i2 = nn3(2,j3s,i1)
          i3 = nn3(3,j3s,i1)
          i4 = nn3(4,j3s,i1)
          if( i2.lt.i1 .or. i3.lt.i1 .or. i4.lt.i1 ) goto 200
c
          if( ijdist(i1,i4).ne.3 )then
            write(IWR ,*) 'ERROR: 3nn.ne.3',i1,i2,i3,i4
            write(IWR ,*) 'ERROR: 3nn.ne.3',i1,i2,i3,i4
            STOP 'nn3.ne.3'
          endif
c
          do 100 k3s=j3s+1,n3s
c
c check that endpoints match
            if( nn3(4,k3s,i1).ne.i4 ) goto 100
c
c check indices greater than kickoff point (prevents duplicates)
            i5 = nn3(3,k3s,i1)
            i6 = nn3(2,k3s,i1)
            if( i5.lt.i1 .or. i6.lt.i1 ) goto 100
c
c check that have legit ring - no shorter paths.
            if( ijdist(i2,i5).ne.3 .or. ijdist(i3,i6).ne.3 ) goto 100
c
c Trace steps to be sure have true ring and not periodic chain
            call CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $                   6, i1,i2,i3,i4,i5,i6,i7,i8,i9 )
            if( badrng )then
              n6chain = n6chain + 1
C              write(IRNG,9160)  n6chain, i1,i2,i3,i4,i5,i6
              goto 100
            endif
c
c we have a 6-ring, go ahead and log it
            n6ring = n6ring + 1
            call LOC2NN( i1ii, i6,i1,i2, natom,num2,nn2 )
            i2ii = nn3(6,j3s,i1)
            i3ii = nn3(7,j3s,i1)
            call LOC2NN( i4ii, i3,i4,i5, natom,num2,nn2 )
            i5ii = nn3(7,k3s,i1)
            i6ii = nn3(6,k3s,i1)
            if( nn2(2,i1ii,i1).eq.0 ) nn2(2,i1ii,i1) = 6
            if( nn2(2,i2ii,i2).eq.0 ) nn2(2,i2ii,i2) = 6
            if( nn2(2,i3ii,i3).eq.0 ) nn2(2,i3ii,i3) = 6
            if( nn2(2,i4ii,i4).eq.0 ) nn2(2,i4ii,i4) = 6
            if( nn2(2,i5ii,i5).eq.0 ) nn2(2,i5ii,i5) = 6
            if( nn2(2,i6ii,i6).eq.0 ) nn2(2,i6ii,i6) = 6
            irngs = nn2(2,i1ii,i1) + nn2(2,i2ii,i2) + nn2(2,i3ii,i3)
     $            + nn2(2,i4ii,i4) + nn2(2,i5ii,i5) + nn2(2,i6ii,i6)
            if( irngs.eq.36 )then
c             Clean 6-ring
              write(IRNG,9060)  i1,i2,i3,i4,i5,i6
CNUM              write(IRNG,9066)  n6ring, i1,i2,i3,i4,i5,i6
            else
c             Part of it has been used before in a smaller ring
              write(IRNG,9061)  i1,i2,i3,i4,i5,i6
CNUM              write(IRNG,9067)  n6ring, i1,i2,i3,i4,i5,i6
            endif
c
  100     continue
c
  200   continue
c
  300 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COUNT7
c
c
      SUBROUTINE COUNT7( IWR,IRNG, n7ring, n7chain, natom, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3 )
c---------------------------------------------------------------
c Purpose: find and count 7-member rings
c---------------------------------------------------------------
c
c This routine finds 7-rings by starting at fixed atom and finding
c different 3nn sets whose endpoints are n.n.  This routine will
c check if 7-"rings" actually are 1-D periodic chains connecting
c an atom and its image in another unit cell.  Also checks that
c all diagonals are no less than three atoms apart, taking into
c account possibility of linking between rings across u.c. bound.
c
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom), nn1lat(3,5,natom)
      DIMENSION  num2(natom), nn2(3,10,natom)
      DIMENSION  num3(natom), nn3(7,36,natom)
c local:
      LOGICAL    badrng
c
c >>>> EXECUTABLE CODE:
c
c
 9070 format(1x,'7-Ring     = ',9i4)
 9071 format(1x,'7-Ringlaps = ',9i4)
 9170 format(1x,'7-Chain#   ',i4,' = ',9i4)
 9077 format(1x,'7-Ring#    ',i4,' = ',9i4)
 9078 format(1x,'7-Ringlaps#',i4,' = ',9i4)
c
      n7ring = 0
      n7chain = 0
c
      do 300 i1=1,natom
c
        do 200 j3s=1,num3(i1)-1
          i2 = nn3(2,j3s,i1)
          i3 = nn3(3,j3s,i1)
          i4 = nn3(4,j3s,i1)
          if( i2.lt.i1 .or. i3.lt.i1 .or. i4.lt.i1 ) goto 200
c
          do 100 k3s=j3s+1,num3(i1)
            i5 = nn3(4,k3s,i1)
            if( ijdist(i5,i4).ne.1 ) goto 100
c
            i6 = nn3(3,k3s,i1)
            i7 = nn3(2,k3s,i1)
            if( i5.lt.i1 .or. i6.lt.i1 .or. i7.lt.i1 ) goto 100
c
c could be a valid seven-ring - check paths.
            if( ijdist(i2,i5).ne.3 .or. ijdist(i2,i6).ne.3 .or.
     $          ijdist(i7,i4).ne.3 .or. ijdist(i7,i3).ne.3 .or.
     $          ijdist(i3,i6).ne.3 ) goto 100
c
c Trace steps to be sure have true ring and not periodic chain
            call CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $                   7, i1,i2,i3,i4,i5,i6,i7,i8,i9 )
            if( badrng )then
              n7chain = n7chain + 1
C              write(IRNG,9170)  n7chain, i1,i2,i3,i4,i5,i6,i7
              goto 100
            endif
c
c Have a valid 7-member ring
            n7ring = n7ring + 1
c
            call LOC2NN(i1ii, i7,i1,i2, natom,num2,nn2)
            i2ii = nn3(6,j3s,i1)
            i3ii = nn3(7,j3s,i1)
            call LOC2NN(i4ii, i3,i4,i5, natom,num2,nn2)
            call LOC2NN(i5ii, i4,i5,i6, natom,num2,nn2)
            i6ii = nn3(7,k3s,i1)
            i7ii = nn3(6,k3s,i1)
            if( nn2(2,i1ii,i1).eq.0 ) nn2(2,i1ii,i1) = 7
            if( nn2(2,i2ii,i2).eq.0 ) nn2(2,i2ii,i2) = 7
            if( nn2(2,i3ii,i3).eq.0 ) nn2(2,i3ii,i3) = 7
            if( nn2(2,i4ii,i4).eq.0 ) nn2(2,i4ii,i4) = 7
            if( nn2(2,i5ii,i5).eq.0 ) nn2(2,i5ii,i5) = 7
            if( nn2(2,i6ii,i6).eq.0 ) nn2(2,i6ii,i6) = 7
            if( nn2(2,i7ii,i7).eq.0 ) nn2(2,i7ii,i7) = 7
            irngs = nn2(2,i1ii,i1) + nn2(2,i2ii,i2) + nn2(2,i3ii,i3)
     $            + nn2(2,i4ii,i4) + nn2(2,i5ii,i5) + nn2(2,i6ii,i6)
     $            + nn2(2,i7ii,i7)
            if( irngs.eq.49 )then
c             Have a virgin ring
              write(IRNG,9070)  i1,i2,i3,i4,i5,i6,i7
CNUM              write(IRNG,9077)  n7ring, i1,i2,i3,i4,i5,i6,i7
            else
c             Part of it has been used before
              write(IRNG,9071)  i1,i2,i3,i4,i5,i6,i7
CNUM              write(IRNG,9078)  n7ring, i1,i2,i3,i4,i5,i6,i7
            endif
c       
  100     continue
c
  200   continue
c
  300 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MK8NN3
c
c
      SUBROUTINE MK8NN3( IWR,IRNG, n8ring, n8chain, natom, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
c---------------------------------------------------------------
c Purpose: find and count 8-member rings
c---------------------------------------------------------------
c
c This routine finds 8-rings by starting at fixed atom and finding
c different 3nn sets whose endpoints are 2nn.  This routine will
c check if 8-"rings" actually are 1-D periodic chains connecting
c an atom and its image in another unit cell.  Also checks that
c all diagonals are no less than four atoms apart, taking into
c account possibility of linking between rings across u.c. bound.
c Note: this routine replaced by 4nn-generating routine.
c
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom), nn1lat(3,5,natom)
      DIMENSION  num2(natom), nn2(3,10,natom)
      DIMENSION  num3(natom), nn3(7,36,natom)
      DIMENSION  num4(natom), nn4(9,108,natom)
c local:
      LOGICAL    badrng
      COMMON  /ringat/ i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
      DIMENSION  iring(10)
      EQUIVALENCE  (iring,i1)
c
c >>>> EXECUTABLE CODE:
c
c
 9080 format(1x,'8-Ring     = ',9i4)
 9081 format(1x,'8-Ringlaps = ',9i4)
 9180 format(1x,'8-Chain#   ',i4,' = ',9i4)
 9088 format(1x,'8-Ring#    ',i4,' = ',9i4)
 9089 format(1x,'8-Ringlaps#',i4,' = ',9i4)
c
      n8ring = 0
      n8chain = 0
c
      do 300 i1=1,natom
c
        do 200 j3s=1,num3(i1)-1
          i2 = nn3(2,j3s,i1)
          i3 = nn3(3,j3s,i1)
          i4 = nn3(4,j3s,i1)
          if( i2.lt.i1 .or. i3.lt.i1 .or. i4.lt.i1 ) goto 200
c
          do 100 k3s=j3s+1,num3(i1)
            i6 = nn3(4,k3s,i1)
            if( ijdist(i6,i4).ne.2 ) goto 100
c
            i7 = nn3(3,k3s,i1)
            i8 = nn3(2,k3s,i1)
            if( i6.lt.i1 .or. i7.lt.i1 .or. i8.lt.i1 ) goto 100
c
c find missing atom link(s)
            do 80 i4nn=1,num1(i4)
              i5 = nn1(i4nn,i4)
              if( ijdist(i5,i6).ne.1 .or. i5.lt.i1 ) goto 80
c
c Trace steps to be sure have true ring and not periodic chain
              call CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $                     8, i1,i2,i3,i4,i5,i6,i7,i8,i9 )
              if( badrng )then
                n8chain = n8chain + 1
                goto 80
              endif
c
c could be a valid 8-ring - need to check paths ...
              if( ijdist(i1,i5).lt.4 )then
                call CROSS4( badrng, i1,i2,i3,i4,i5, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 80
              endif
              if( ijdist(i2,i6).lt.4 )then
                call CROSS4( badrng, i2,i3,i4,i5,i6, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 80
              endif
              if( ijdist(i3,i7).lt.4 )then
                call CROSS4( badrng, i3,i4,i5,i6,i7, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 80
              endif
              if( ijdist(i4,i8).lt.4 )then
                call CROSS4( badrng, i4,i5,i6,i7,i8, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 80
              endif
c
c Have a valid 8-member ring
              n8ring = n8ring + 1
c
              call LOC2NN(i1ii, i8,i1,i2, natom,num2,nn2)
              i2ii = nn3(6,j3s,i1)
              i3ii = nn3(7,j3s,i1)
              call LOC2NN(i4ii, i3,i4,i5, natom,num2,nn2)
              call LOC2NN(i5ii, i4,i5,i6, natom,num2,nn2)
              call LOC2NN(i6ii, i5,i6,i7, natom,num2,nn2)
              i7ii = nn3(7,k3s,i1)
              i8ii = nn3(6,k3s,i1)
              if( nn2(2,i1ii,i1).eq.0 ) nn2(2,i1ii,i1) = 8
              if( nn2(2,i2ii,i2).eq.0 ) nn2(2,i2ii,i2) = 8
              if( nn2(2,i3ii,i3).eq.0 ) nn2(2,i3ii,i3) = 8
              if( nn2(2,i4ii,i4).eq.0 ) nn2(2,i4ii,i4) = 8
              if( nn2(2,i5ii,i5).eq.0 ) nn2(2,i5ii,i5) = 8
              if( nn2(2,i6ii,i6).eq.0 ) nn2(2,i6ii,i6) = 8
              if( nn2(2,i7ii,i7).eq.0 ) nn2(2,i7ii,i7) = 8
              if( nn2(2,i8ii,i8).eq.0 ) nn2(2,i8ii,i8) = 8
              irngs = nn2(2,i1ii,i1) + nn2(2,i2ii,i2) + nn2(2,i3ii,i3)
     $              + nn2(2,i4ii,i4) + nn2(2,i5ii,i5) + nn2(2,i6ii,i6)
     $              + nn2(2,i7ii,i7) + nn2(2,i8ii,i8)
              if( irngs.eq.64 )then
c               Have a virgin ring
                write(IRNG,9080)  i1,i2,i3,i4,i5,i6,i7,i8
CNUM                write(IRNG,9088)  n8ring, i1,i2,i3,i4,i5,i6,i7,i8
              else
c               Part of it has been used before
                write(IRNG,9081)  i1,i2,i3,i4,i5,i6,i7,i8
CNUM                write(IRNG,9089)  n8ring, i1,i2,i3,i4,i5,i6,i7,i8
              endif
   80       continue
c       
  100     continue
c
  200   continue
c
  300 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COUNT8
c
c
      SUBROUTINE COUNT8( IWR,IRNG, n8ring, n8chain, natom, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
c---------------------------------------------------------------
c Purpose: find and count 8-member rings
c---------------------------------------------------------------
c
c This routine finds 8-rings by starting at fixed atom and finding
c different 4nn sets whose endpoints match.  This routine will
c check if 8-"rings" actually are 1-D periodic chains connecting
c an atom and its image in another unit cell.  Also checks that
c all diagonals are no less than four atoms apart, taking into
c account possibility of linking between rings across u.c. bound.
c
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom), nn1lat(3,5,natom)
      DIMENSION  num2(natom), nn2(3,10,natom)
      DIMENSION  num3(natom), nn3(7,36,natom)
      DIMENSION  num4(natom), nn4(9,108,natom)
c local:
      LOGICAL    badrng
c
c >>>> EXECUTABLE CODE:
c
c
 9080 format(1x,'8-Ring     = ',9i4)
 9081 format(1x,'8-Ringlaps = ',9i4)
 9180 format(1x,'8-Chain#   ',i4,' = ',9i4)
 9088 format(1x,'8-Ring#    ',i4,' = ',9i4)
 9089 format(1x,'8-Ringlaps#',i4,' = ',9i4)
c
      n8ring = 0
      n8chain = 0
c
      do 300 i1=1,natom
        n4s = num4(i1)
        if( n4s.le.1 ) goto 300
c
        do 200 i4s=1,n4s-1
          i2 = nn4(2,i4s,i1)
          i3 = nn4(3,i4s,i1)
          i4 = nn4(4,i4s,i1)
          i5 = nn4(5,i4s,i1)
          if( i2.lt.i1 .or. i3.lt.i1 .or. i4.lt.i1 ) goto 200
c
          do 100 j4s=i4s+1,n4s
            if( nn4(5,j4s,i1).eq.i5 )then
c
              i8 = nn4(2,j4s,i1)
              i7 = nn4(3,j4s,i1)
              i6 = nn4(4,j4s,i1)
              if( i8.lt.i1 .or. i7.lt.i1 .or. i6.lt.i1 ) goto 100
c
c Have an viable unique candidate
c   Trace steps to be sure have true ring and not periodic chain
              call CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $                     8, i1,i2,i3,i4,i5,i6,i7,i8,i9 )
              if( badrng )then
                n8chain = n8chain + 1
                goto 100
              endif
c
c could be a valid 8-ring - need to check paths ...
              if( ijdist(i1,i5).lt.4 )then
                call CROSS4( badrng, i1,i2,i3,i4,i5, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
              if( ijdist(i2,i6).lt.4 )then
                call CROSS4( badrng, i2,i3,i4,i5,i6, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
              if( ijdist(i3,i7).lt.4 )then
                call CROSS4( badrng, i3,i4,i5,i6,i7, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
              if( ijdist(i4,i8).lt.4 )then
                call CROSS4( badrng, i4,i5,i6,i7,i8, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
c
c Have a valid 8-member ring
              n8ring = n8ring + 1
c
              call LOC2NN(i1ii, i8,i1,i2, natom,num2,nn2)
              call LOC2NN(i2ii, i1,i2,i3, natom,num2,nn2)
              call LOC2NN(i3ii, i2,i3,i4, natom,num2,nn2)
              call LOC2NN(i4ii, i3,i4,i5, natom,num2,nn2)
              call LOC2NN(i5ii, i4,i5,i6, natom,num2,nn2)
              call LOC2NN(i6ii, i5,i6,i7, natom,num2,nn2)
              call LOC2NN(i7ii, i6,i7,i8, natom,num2,nn2)
              call LOC2NN(i8ii, i7,i8,i1, natom,num2,nn2)
              if( nn2(2,i1ii,i1).eq.0 ) nn2(2,i1ii,i1) = 8
              if( nn2(2,i2ii,i2).eq.0 ) nn2(2,i2ii,i2) = 8
              if( nn2(2,i3ii,i3).eq.0 ) nn2(2,i3ii,i3) = 8
              if( nn2(2,i4ii,i4).eq.0 ) nn2(2,i4ii,i4) = 8
              if( nn2(2,i5ii,i5).eq.0 ) nn2(2,i5ii,i5) = 8
              if( nn2(2,i6ii,i6).eq.0 ) nn2(2,i6ii,i6) = 8
              if( nn2(2,i7ii,i7).eq.0 ) nn2(2,i7ii,i7) = 8
              if( nn2(2,i8ii,i8).eq.0 ) nn2(2,i8ii,i8) = 8
              irngs = nn2(2,i1ii,i1) + nn2(2,i2ii,i2) + nn2(2,i3ii,i3)
     $              + nn2(2,i4ii,i4) + nn2(2,i5ii,i5) + nn2(2,i6ii,i6)
     $              + nn2(2,i7ii,i7) + nn2(2,i8ii,i8)
              if( irngs.eq.64 )then
c               Have a virgin ring
                write(IRNG,9080)  i1,i2,i3,i4,i5,i6,i7,i8
CNUM                write(IRNG,9088)  n8ring, i1,i2,i3,i4,i5,i6,i7,i8
              else
c               Part of it has been used before
                write(IRNG,9081)  i1,i2,i3,i4,i5,i6,i7,i8
CNUM                write(IRNG,9089)  n8ring, i1,i2,i3,i4,i5,i6,i7,i8
              endif
c
            endif
c
  100     continue
c
  200   continue
c
  300 continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CROSS4
c
c
       SUBROUTINE CROSS4( linked, i1,i2,i3,i4,i5, natom,
     $   ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
c---------------------------------------------------------------
c Purpose: I have a closed 8(9)-member string.  I want to check
c that path i1-i5 across the ring is not less than four.  This
c is tricky because I have to check that shorter path does not
c link between unit cells rather than across the ring. I assume
c 4-rings are guaranteed and that larger "rings" are suspect.
c---------------------------------------------------------------
c
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom), nn1lat(3,5,natom)
      DIMENSION  num2(natom), nn2(3,10,natom)
      DIMENSION  num3(natom), nn3(7,36,natom)
      DIMENSION  num4(natom), nn4(9,108,natom)
      LOGICAL    linked
c local:
      LOGICAL    badrng
c
      linked = .false.
c
      if( ijdist(i1,i5).ge.4 )then
c if not already recorded as closer than four, I am home
        if( ijdist(i1,i5).gt.4 )then
          write(*,*) 'CROSS4: missing 4nn=',i1,i5
          ijdist(i1,i5) = 4
          ijdist(i5,i1) = 4
        endif
c
      elseif( ijdist(i1,i5).eq.0 )then
c       Linked squares
        linked = .true.
c
      elseif( ijdist(i1,i5).eq.1 )then
c       Linked i1-i5 giving 5-ring?
        call CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $               5, i1,i2,i3,i4,i5,i6,i7,i8,i9 )
        if( .not. badrng ) linked = .true.
c
      elseif( ijdist(i1,i5).eq.2 )then
c       Linked i5-j6-i1 giving 6-ring?  Find j6.
        do 60 i5nn=1,num1(i5)
          j6 = nn1(i5nn,i5)
          if( ijdist(j6,i1).eq.1 )then
            call CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $                   6, i1,i2,i3,i4,i5,j6,i7,i8,i9 )
            if( .not. badrng ) linked = .true.
            goto 61
          endif
   60   continue
        STOP  'CROSS4-1'
   61   continue
c
      elseif( ijdist(i1,i5).eq.3 )then
c       Linked i5-j6-j7-i1 giving 7-ring?  Find connecting triple.
c       Actually indexing hard, try nn-nn match as straightforward
        do 75 i5nn=1,num1(i5)
          j6 = nn1(i5nn,i5)
          do 70 i1nn=1,num1(i1)
            j7 = nn1(i1nn,i1)
            if( ijdist(j6,j7).eq.1 )then
              call CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $                     7, i1,i2,i3,i4,i5,j6,j7,i8,i9 )
              if( .not. badrng )then
                linked = .true.
                goto 76
              endif
            endif
   70     continue
   75   continue
   76   continue
c
      endif
c       
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COUNT9
c
c
      SUBROUTINE COUNT9( IWR,IRNG, n9ring, n9chain, natom, ijdist,
     $ num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
c---------------------------------------------------------------
c Purpose: find and count 9-member rings
c---------------------------------------------------------------
c
c This routine finds 9-rings by starting at fixed atom and finding
c different 4nn sets whose endpoints are n.n.  This routine will
c check if 9-"rings" actually are 1-D periodic chains connecting
c an atom and its image in another unit cell.  Also checks that
c all diagonals are no less than four atoms apart, taking into
c account possibility of linking between rings across u.c. bound.
c
      DIMENSION  ijdist(natom,natom)
      DIMENSION  num1(natom), nn1(5,natom), nn1lat(3,5,natom)
      DIMENSION  num2(natom), nn2(3,10,natom)
      DIMENSION  num3(natom), nn3(7,36,natom)
      DIMENSION  num4(natom), nn4(9,108,natom)
c local:
      LOGICAL    badrng
c
c >>>> EXECUTABLE CODE:
c
c
 9090 format(1x,'9-Ring     = ',9i4)
 9091 format(1x,'9-Ringlaps = ',9i4)
 9190 format(1x,'9-Chain#   ',i4,' = ',9i4)
 9099 format(1x,'9-Ring#    ',i4,' = ',9i4)
 9098 format(1x,'9-Ringlaps#',i4,' = ',9i4)
c
      n9ring = 0
      n9chain = 0
c
      do 300 i1=1,natom
        n4s = num4(i1)
        if( n4s.le.1 ) goto 300
c
        do 200 i4s=1,n4s-1
          i2 = nn4(2,i4s,i1)
          i3 = nn4(3,i4s,i1)
          i4 = nn4(4,i4s,i1)
          i5 = nn4(5,i4s,i1)
          if( i2.lt.i1 .or. i3.lt.i1 .or. i4.lt.i1 ) goto 200
c
          do 100 j4s=i4s+1,n4s
            i6 = nn4(5,j4s,i1)
            if( ijdist(i5,i6).eq.1 .and. i1.lt.i6 )then
c
              i9 = nn4(2,j4s,i1)
              i8 = nn4(3,j4s,i1)
              i7 = nn4(4,j4s,i1)
              if( i8.lt.i1 .or. i7.lt.i1 .or. i6.lt.i1 ) goto 100
c
c Have an viable unique candidate
c   Trace steps to be sure have true ring and not periodic chain
              call CHKLAT( badrng,  natom, num1,nn1,nn1lat,
     $                     9, i1,i2,i3,i4,i5,i6,i7,i8,i9 )
              if( badrng )then
                n9chain = n9chain + 1
                goto 100
              endif
c
c could be a valid 9-ring - need to check paths ...
              if( ijdist(i1,i5).lt.4 )then
                call CROSS4( badrng, i1,i2,i3,i4,i5, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
              if( ijdist(i2,i6).lt.4 )then
                call CROSS4( badrng, i2,i3,i4,i5,i6, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
              if( ijdist(i3,i7).lt.4 )then
                call CROSS4( badrng, i3,i4,i5,i6,i7, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
              if( ijdist(i4,i8).lt.4 )then
                call CROSS4( badrng, i4,i5,i6,i7,i8, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
              if( ijdist(i5,i9).lt.4 )then
                call CROSS4( badrng, i5,i6,i7,i8,i9, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
              if( ijdist(i6,i1).lt.4 )then
                call CROSS4( badrng, i6,i7,i8,i9,i1, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
              if( ijdist(i7,i2).lt.4 )then
                call CROSS4( badrng, i7,i8,i9,i1,i2, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
              if( ijdist(i8,i3).lt.4 )then
                call CROSS4( badrng, i8,i9,i1,i2,i3, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
              if( ijdist(i9,i4).lt.4 )then
                call CROSS4( badrng, i9,i1,i2,i3,i4, natom,
     $          ijdist, num1,nn1,nn1lat, num2,nn2, num3,nn3, num4,nn4 )
                if( badrng ) goto 100
              endif
c
c Have a valid 9-member ring
              n9ring = n9ring + 1
c
              call LOC2NN(i1ii, i9,i1,i2, natom,num2,nn2)
              call LOC2NN(i2ii, i1,i2,i3, natom,num2,nn2)
              call LOC2NN(i3ii, i2,i3,i4, natom,num2,nn2)
              call LOC2NN(i4ii, i3,i4,i5, natom,num2,nn2)
              call LOC2NN(i5ii, i4,i5,i6, natom,num2,nn2)
              call LOC2NN(i6ii, i5,i6,i7, natom,num2,nn2)
              call LOC2NN(i7ii, i6,i7,i8, natom,num2,nn2)
              call LOC2NN(i8ii, i7,i8,i9, natom,num2,nn2)
              call LOC2NN(i9ii, i8,i9,i1, natom,num2,nn2)
              if( nn2(2,i1ii,i1).eq.0 ) nn2(2,i1ii,i1) = 9
              if( nn2(2,i2ii,i2).eq.0 ) nn2(2,i2ii,i2) = 9
              if( nn2(2,i3ii,i3).eq.0 ) nn2(2,i3ii,i3) = 9
              if( nn2(2,i4ii,i4).eq.0 ) nn2(2,i4ii,i4) = 9
              if( nn2(2,i5ii,i5).eq.0 ) nn2(2,i5ii,i5) = 9
              if( nn2(2,i6ii,i6).eq.0 ) nn2(2,i6ii,i6) = 9
              if( nn2(2,i7ii,i7).eq.0 ) nn2(2,i7ii,i7) = 9
              if( nn2(2,i8ii,i8).eq.0 ) nn2(2,i8ii,i8) = 9
              if( nn2(2,i9ii,i9).eq.0 ) nn2(2,i9ii,i9) = 9
              irngs = nn2(2,i1ii,i1) + nn2(2,i2ii,i2) + nn2(2,i3ii,i3)
     $              + nn2(2,i4ii,i4) + nn2(2,i5ii,i5) + nn2(2,i6ii,i6)
     $              + nn2(2,i7ii,i7) + nn2(2,i8ii,i8) + nn2(2,i9ii,i9)
              if( irngs.eq.81 )then
c               Have a virgin ring
                write(IRNG,9090)  i1,i2,i3,i4,i5,i6,i7,i8,i9
CNUM                write(IRNG,9099)  n8ring, i1,i2,i3,i4,i5,i6,i7,i8,i9
              else
c               Part of it has been used before
                write(IRNG,9091)  i1,i2,i3,i4,i5,i6,i7,i8,i9
CNUM                write(IRNG,9098)  n8ring, i1,i2,i3,i4,i5,i6,i7,i8,i9
              endif
c
            endif
c
  100     continue
c
  200   continue
c
  300 continue
c
      RETURN
      END
