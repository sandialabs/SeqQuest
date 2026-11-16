      PROGRAM  E X B I R F T
C#######################################################################
C     Fits energy versus volume data to the form
C
C     E(V) = Sum [ a(n) v^(-2n/3) , {n,0,N}]
C
C     and then calculates E(Vo), Vo, Ko, and Ko'
C
C     This version uses the Least Squares fitting routine of
C      the Numerical Recipes book.
C
C     We've added allowances for weighting of energy points
C                                           mjm  1 March 1990
C Revisions:
c  15Dec14-PAS/1.10: install better continuation options, clean version code
c   8Aug08-PAS: to work off volume rather than lattice
c  PAS - to integrate in to quest environment, 1992-2000
C#######################################################################
c
      IMPLICIT DOUBLE PRECISION (a-h,k,o-z)
c
      CHARACTER  version*5,         lastup*11
      DATA       version /'1.10 '/, lastup/'15Dec14-PAS'/
c
      PARAMETER (maxn=1000,maxm=10,maxm1=maxm+1)
      PARAMETER (maxcol=10)
C
      COMMON  v(maxn),e(maxn),a(0:maxm),p(0:maxm,maxn),
     1   earray(0:3),sig(maxn),cvm(0:maxm,0:maxm),asig(0:maxm),
     2   col(maxcol),lista(0:maxm)
C
      LOGICAL    details,brief
      LOGICAL    cona,disb,engr,minfound,lsig
C  Conversion factors:  V=CONV*(A**3)
      DIMENSION  conv(6)
C                         fcc   sc   bcc   Td
      DATA       conv / 0.25d0,1.d0,0.5d0,2.0d0,0.d0,1.d0 /
      CHARACTER*72 datafl,outfl
      CHARACTER*75 title
      CHARACTER*1  blank,ans,yesu,yesl
      DATA  blank/' '/,yesu/'y'/,yesl/'Y'/
C     New (1986) Values for 1 Bohr in  and 1 Rydberg in eV:
      DATA  auoan/5.29177249d-1/, ryoev/1.36056981d1/
      DATA  third/3.33333333333333333d-1/, two3/6.666666666666666667d-1/
      DATA  zero,one/0.d0,1.d0/
C set io units:
      DATA  IRD,IWR,IDAT /5,6,7/
C
C >>>>> ECECUTABLE CODE:
C
      write(IWR,9000)  version, lastup
 9000 format(1x,'Program BIRCH: fit EoS with Birch-Murnaghan equation',
     $  /1x,'Version:  ',a,'     Last revised: ',a,
     $  /1x,'Written by Mike Mehl (1990)',
     $  /1x,'Cite use with:',
     $  /5x,'M.J. Mehl, Phys. Rev. B 47, 2493 (1993)',
     $  /5x,'F. Birch, J. Geophys. Res. 83, 1257 (1978)',
     $  /5x,'F.D. Murnaghan, Am. J. Math. 49, 235 (1937)',
     $  / )
      details = .false.
      brief = .true.
      au3 = auoan*auoan*auoan
C  We assume all fitting parameters need to be fitted.
c  Note that LISTA requires 1 + our index:
      do 10 i=0,maxm
        lista(i) = i + 1
   10 continue
c
c     **********  Get file to be fitted  **********
c
  100 continue
c
      write(IWR,*) 'Enter name of data file:'
      read(IRD,8100,err=999,end=999)  datafl
 8100 format(a72)
      if( datafl .eq. blank ) STOP
C
  110 write(IWR,9110)  
 9110 format(1x,'Enter lattice type'
     $, ' (1-fcc;2-sc;3-bcc;4-Td;5-X;6=Vol):  ',$)
      read(IRD,*)  latopt
c
      if( latopt.eq.0 )then
c       toggle brief
        brief = .not. brief
        goto 110
      elseif( latopt.gt.4 .or. latopt.lt.0 )then
c       Non-standard cell volume input
        if( latopt .eq. 6 )then
c         We are inputing the volume directly, not lattice constant
          iwhich = 1
          goto 129
        else
          latopt = 5
          write(IWR,9111)
 9111     format(1x,'Enter volume factor q (such that V=q*a**3):  ',$)
          read(IRD,*)  conv(5)
          if( conv(5).le.zero ) goto 110
        endif
      endif
C
      iwhich = 2
      if( details )then
        write(IWR,9120)  
 9120   format(1x,'Enter 1 for volume, 2 for lattice constant:  ',$)
        read(IRD,*)  iwhich
      endif
  129 continue
      cona = iwhich.eq.2
      if( .not.cona ) iwhich = 1
C
      iaunit = 2
      if( details )then
        write(IWR,9140)  
 9140   format(1x,'Distances in (1) Angstroms, or (2) Bohrs?  ',$)
        read(IRD,*)  iaunit
      endif
      disb = iaunit.eq.2
      if( .not.disb ) iaunit = 1
      if( disb ) then
        aunit = one
      else
        aunit = auoan
      endif
      vunit = aunit*aunit*aunit
C
      ieunit = 2
      if( details )then
        write(IWR,9150)
 9150   format(1x,'Energy units in (1) eV, or (2) Rydbergs?  ',$)
        read(ird,*)  ieunit
      endif
      engr = ieunit.eq.2
      if(.not.engr) ieunit = 1
      if( engr ) then
        eunit = one
      else
        eunit = ryoev
      endif
C
  160 continue
      write(IWR,9160)
 9160 format(1x,'How many lines should we skip?  ',$)
      read(IRD,*)  lskip
      OPEN(unit=IDAT,file=datafl,status='old',form='formatted',err=100)
      do 170 ls=1,lskip
        read(IDAT,8170)  title
        write(IWR,8170)  title
 8170   format(a75)
  170 continue
C
  180 continue
      write(IWR,9180)  maxcol
 9180 format( ' Enter column: for size, for energy (MAX=',i3,'):  ',$)
      read(IRD,*)  icx,ice
      icol = MAX(icx,ice)
      if( icol.gt.maxcol ) goto 180
C
      lsig = .false.
      if( details )then
        write(IWR,9190)
 9190   format(1x,'Are there error estimates for the energies?  ',$)
        read(IRD,8190)  ans
 8190   format(a1)
        lsig = (ans.eq.yesu) .or. (ans.eq.yesl)
      endif
      if( lsig )then
        write(IWR,9191)
 9191   format(1x,'Enter column for the errors:  ',$)
        read(IRD,*) icsig
      endif
C
C  To find minimum energy, we need estimated starting value.
C  We might as well use volume of the lowest energy as estimated Vo
C
      emin = 1.d10
      vmin = zero
C
C  Read data in the form X,E (x=volume or lattice constant, E=energy)
C  where Angstroms and eV or atomic units (Bohrs and Rydbergs) are
C  used depending upon the setting of IxUNITs
C
      n = 0
  200 continue
c     read(IDAT,*,err=200,end=290) (col(i),i=1,icol)
      read(IDAT,*,err=290,end=290) (col(i),i=1,icol)
      n = n + 1
      x = col(icx)
      ee = col(ice)
      if( cona )then
        va = conv(latopt)*x*x*x
      else
        va = x
      endif
C  Volume in au**3:
      v(n) = va/vunit
C  Energy in Rydbergs:
      e(n) = ee/eunit
 9210 format(1x,i5,5f14.5)
      if( lsig )then
        sig(n) = col(icsig)/eunit
        if( .not. brief ) write(IWR,9210) n,va,ee,v(n),e(n),sig(n)
      else
C       Scale to an arbitrary constant value:
        sig(n) = one
        if( .not. brief ) write(IWR,9210)  n,va,ee,v(n),e(n)
      endif
      if(e(n).lt.emin) then
        emin = e(n)
        vmin = v(n)
      endif
      goto 200
C
  290 continue
      CLOSE(unit=IDAT)
      write(IWR,*) 'Number of valid data points found = ',n
      if( n.eq.0 ) STOP 'no point'
C   XSTART is the starting value for the minimum search:
      xstart = vmin**(-two3)
      write(IWR,9290)  emin,vmin,xstart
 9290 format(/' Minimum energy in file = ',F15.5,' at V = ',2F15.7/)
C
C     Set up the fit.  Note that the functions to be fitted are
C      v^(-2n/3),n=0,1,2,...MAXM
C
  300 continue
      write(IWR,9300)  maxm
 9300 format(1x,'Enter order of v^(-2/3) polynomial (.le.',i3,'):  ',$)
      read(IRD,*,err=999,end=999)  meos
      if( meos.eq.0 ) goto 160
      if( meos.lt.0 .or. meos.gt.maxm ) goto 999
C
C  Establish the basis functions:
      meos1 = meos + 1
      do 342 i=1,n
        pn = one
        x = v(i)**(-two3)
        do 341 j=0,meos
          p(j,i) = pn
          pn = pn * x
  341   continue
  342 continue
C
      call LFIT(e,sig,n,a,meos1,lista,meos1,cvm,maxm1,chisq,p)
C
C  Figure the rms "error" in each term:
C  CHI is the numerical recipies definition, unless N=2, in which
C  case we won't get a good fit anyway
C
      if( n.gt.2 )then
        chi = SQRT(chisq/(n-2))
      else
        chi = SQRT(chisq)
        write(IWR,9350)
 9350   format(/' Warning:  N<3, so delta A''s are not well defined')
      endif
      if( lsig )then
C       CVM's give direct estimates of the errors:
        do 360 i=0,meos
          asig(i) = SQRT(cvm(i,i))
  360   continue
      else
C       Weight errors by CHI:
        do 370 i=0,meos
          asig(i) = SQRT(cvm(i,i))*chi
  370   continue
      endif
      if( .not. brief ) write(IWR,9380)  (i,a(i),asig(i),i=0,meos)
 9380 format(/' Fitting coefficients:'/(1x,i5,1p,2d16.8))
C
C  Now for the error checking:
C
      erms = zero
      emax = zero
      do 390 i=1,n
        xo = v(i)**(-two3)
        call PLYEVL(m,zero,a,xo,0,earray)
        eck = earray(0)
        err = eck - e(i)
        if( ABS(err).gt.ABS(emax) ) emax = err
        erms = erms + err*err
        write(IWR,9390)  i,v(i),e(i),eck,err
 9390   format(1x,i5,f12.5,3f16.6)
  390 continue
      erms = SQRT(erms/n)
C  Convert the errors to eV
      ermsev = erms*ryoev
      emaxev = emax*ryoev
C
C  Now we must find the equilibrium volume, VO, if we can.  We'll use
C  Newton's method.  Note we can write the energy as E(v)=f(v^(-2/3)).  
C  Then every extremum of E is also an extremum of f.
C  Pick a trial value for XO.  Eventually, VO=XO**(-3/2).  Logical
C  starting value is VMIN**(-2/3), which we calculated above. After
C  first minimization, we should have good estimate of Vo, use that
C
      xo = xstart
C  Since we can calculate exact second derivatives of the fit for M>1,
C  we can use Newton's method to find the root.  Try max of 100 times
      difmax = 1.d-6
      do 410 it=1,100
        call PLYEVL(m,zero,a,xo,2,earray)
c       EARRAY contains f and its derivatives
        diff = earray(1)/earray(2)
        if( .not. brief )write(IWR,9410)  xo,earray(1),earray(2),diff
 9410   format(1x,1p4d16.7)
        xo = xo - diff
        if( ABS(diff).lt.difmax ) goto 420
  410 continue
      minfound = .false.
C  No roots found.  Ask if fit should be printed out anyway:
      write(IWR,9413)
 9413 FORMAT(/1x,'No minimum found after 100 iterations.')
      goto 500
C
  420 minfound = .true.
C  Use this new value of XO as the new starting value
      xstart = xo
C  XOI is the 2/3rd root of the volume:
      xoi = one/xo
      vo = xoi**1.5d0
C  Now for the other equilibrium constants:
      call PLYEVL(m,zero,a,xo,3,earray)
      if( .not. brief ) write(IWR,9420)  xo,(earray(i),i=0,3)
 9420 format(1x,1p,5d15.6)
      eo = earray(0)
C     To compare with "universal" equations of state, find the
C      predicted energy at infinity:
      einf = a(0)
      ecoh = einf - eo
C     Check that we've really found an extremum:
      po = -two3*earray(1)*(xo**2.5d0)
      if( .not. brief ) write(IWR,9430)  po
 9430 format(' "Equilibrium" pressure = ',1p,e15.6)
      ko = (1d1*earray(1)+4d0*xo*earray(2))*xo/(9d0*vo)
      kop = (2.5d1*earray(1)+4d0*xo*(6d0*earray(2)+xo*earray(3)))/
     1      (3d0*(5d0*earray(1)+2d0*xo*earray(2)))
C Convert to standard units:
C  Energies in eV:
      eov = eo*ryoev
      ecohv = ecoh*ryoev
C  Pressure in Mbar (1986 conversion factor)
      kov = 1.47105164d2*ko
      kog = 100.d0 * kov
C  Volume in A^3
      vov = vo*au3
C  Lattice constant in A
      if( cona )then
c       Working off lattice constant
        alat = (vo/conv(latopt))**third
        alatv = alat*auoan
      else
c       Working off volume, "assume" cubic lattice
        alat = vo**third
        alatv = alat*auoan
      endif
C
C Now use the ASIG to estimate error terms.  We'll ignore correlation
C between errors for now, and assume that all terms are additive:
C We'll only calculate DELVO, DELEO, and DELKO.  All polynomials in XO
      deleo = zero
      delvo = zero
      delko = zero
      do 440 i=meos,0,-1
        absag = ABS(asig(i))
        deleo = deleo*xo+absag
        tx = two3*i
        delvo = delvo*xo + absag*tx
        delko = delko*xo + absag*tx*(tx+one)
  440 continue
C DELVO and DELKO need corrections:
      delvo = delvo/ko
      delko = delko/vo + abs(ko*kop*delvo/vo)
      delao = third*delvo*alat/vo
C Now convert to standard units:
      deleov = deleo*ryoev
      delko = 1.47105164d2*delko
      delvov = delvo*au3
      delaov = delao*auoan
C
      if( cona )then
c       Lattice parameter optimized
        write(IWR,9491)
     $   meos, alat,alatv, eo,eov, kog,kov,kop, vo,vov, ermsev, emaxev
 9491   format( /1x,'Equilibrium lattice for Birch-Murnaghan',
     $   ' equation of order ',i2,':'
     $   /4x,'a0 = ' ,f10.5, ' Bohr = ', f9.5,' A',
     $   /4x,'Eo = ' ,f20.6, ' Ry = ',f17.6,' eV',
     $   /4x,'Ko = ' ,f10.3, ' GPa  = ', f9.5, ' Mbar',
     5   /4x,'Ko'' = ',     17x         ,f9.5, ' Mbar',
     $   /4x,'Vo = ' ,f20.5, ' Bohr**3 = ',f15.5,' A**3',
     $   /4x,'RMS error in energy fit = ',f13.5,' eV',
     $   /4x,'Maximum error in energy fit = ',f9.5,' eV')
      else
c       Volume-optimized
        write(IWR,9492)
     $   meos, alat,alatv, eo,eov, kog,kov,kop, vo,vov, ermsev, emaxev
 9492   format( /1x,'Equilibrium volume  for Birch-Murnaghan',
     $   ' equation of order ',i2,':'
     $   /4x,'V^(1/3) = ' ,f10.5, ' Bohr = ', f9.5,' A',
     $   /4x,'Eo = ' ,f15.6, ' Ry = ',f17.6,' eV',
     $   /4x,'Ko = ' ,f10.3, ' GPar = ',f9.5, ' Mbar',
     5   /4x,'Ko'' = ',    17x          ,f9.5, ' Mbar',
     $   /4x,'Vo = ' ,f15.5, ' Bohr**3 = ',f15.5,' A**3',
     $   /4x,'RMS error in energy fit = ',f13.5,' eV',
     $   /4x,'Maximum error in energy fit = ',f9.5,' eV')
      endif
C
C     *********  Done, what is next? **********
C
  500 continue
      write(IWR,*)
      write(IWR,*) 'What next(s=Save,p=Poly,f=new Fit,q=Quit)?'
      read(IRD,8190)  ans
c
      if(     ans .eq. 's' .or. ans .eq. 'S' )then
        goto 600
      elseif( ans .eq. 'p' .or. ans .eq. 'P' .or.
     $        ans .eq. 'n' .or. ans .eq. 'N' )then
        goto 300
      elseif( ans .eq. 'f' .or. ans .eq. 'F' )then
        goto 100
      elseif( ans .eq. 'q' .or. ans .eq. 'Q' )then
        goto 900
      else
        write(IWR,9513)  ans
 9513   format(1x,'**** ERROR: invalid option=',a1,
     $   /1x,'Valid options=',
     $   /4x,'  s = (S)ave fit to output file', 
     $   /4x,'  p = new (P)polynomial',
     $   /4x,'  f = (F)it new (F)ile',
     $   /4x,'  q = (Q)uit' )
      endif
      goto 500
c
c     Write
c
  600 continue
      write(IWR,*)  'Name of output file?'
      read(IRD,8100)  outfl
      OPEN(unit=IDAT,file=outfl,status='new',err=600)
c
c     Write out the fit parameters:
      write(IDAT,8520)  meos,(a(i),i=0,meos)
 8520 format(1x,'EoS fit parameters (order',i2,' equation:',
     $ /4x, (1p,4e20.12) )
c
C     Write out the estimated RMS variances:
      write(IDAT,8530)  (asig(i),i=0,meos)
 8530 format(1x,'RMS variances:'/ (1p,4e20.12) )
C     Add the equilibrium functions as a reference:
      if( minfound )then
        write(IDAT,9496) 
     $   meos, alat,alatv, eo,eov, kog,kov,kop, vo,vov
      endif
 9496 format( /1x,'Equilibrium parameters for Birch-Murnaghan',
     $ ' equation of order ',i2,':'
     $ /4x,'a0 = ' ,f10.5, ' Bohr = ', f9.5,' A',
     $ /4x,'Eo = ' ,f20.6, ' Ry =      ',f17.6,' eV',
     $ /4x,'Ko = ' ,f10.3, ' GPa  = ',f9.5, ' Mbar', 
     5 /4x,'Ko'' = ',     17x         ,f9.5, ' Mbar',
     $ /4x,'Vo = ' ,f20.5, ' Bohr**3 = ',f15.5,' A**3' )
C     Also print the errors in the calculation
      write(IDAT,9495)  erms,ermsev,emax,emaxev
 9495 format(' RMS error = ',1p,e11.3,' Ry = ',e11.3,' eV'/
     $ ' Max error = ',e11.3,' Ry = ',e11.3,' eV')
      CLOSE(unit=IDAT)
c
c  Get new option from user:
      goto 500
c
c     **********  All done  **********
  900 continue
  999 continue
c
c    That's all, Folks!
c
      STOP 'bye guy'
      END
C
C
C_______________________________________________________________________
C
C
      SUBROUTINE  PLYEVL (m,x0,a,x,n,p)
C
C     Evaluates the polynomial given by
C       p(x) = Sum [ a(i) (x-x0)^i ,{i,0,m}]
C       and its first N derivatives
C
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
C     Note the dummy indexing on A and P:
      DIMENSION  a(0:m),p(0:n)
C
      y = x - x0
C     Zeroth order term (in Y)
      iprod = 1
      do 10 j=0,n
        p(j) = iprod*a(m)
        iprod = iprod*(m-j)
   10 continue
      do 22 i=m-1,0,-1
        iprod = 1
        do 21 j=0,n
          if( iprod.gt.0 ) p(j) = p(j)*y + iprod*a(i)
          iprod = iprod*(i-j)
   21   continue
   22 continue
c
      RETURN
      END
C
C
C_______________________________________________________________________
C
C     Numerical Recipes routines.  See pages 513 ff of the Fortran edition
C
C_______________________________________________________________________
C
C
      SUBROUTINE LFIT(Y,SIG,NDATA,A,MA,LISTA,MFIT,COVAR,NCVM,CHISQ,P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXN=1000,MAXM=10,MAXM1=MAXM+1)
      PARAMETER (MMAX=MAXM1)
      DIMENSION Y(NDATA),SIG(NDATA),A(MA),LISTA(MFIT),
     *    COVAR(NCVM,NCVM),BETA(MMAX),P(MAXM1,MAXN)
      DATA ZERO,one/0.d0,1.d0/
c
      KK=MFIT+1
      DO 12 J=1,MA
       IHIT=0
       DO 11 K=1,MFIT
        IF (LISTA(K).EQ.J) IHIT=IHIT+1
11      CONTINUE
       IF (IHIT.EQ.0) THEN
        LISTA(KK)=J
        KK=KK+1
       ELSE IF (IHIT.GT.1) THEN
       STOP 'Improper set in LISTA'
       ENDIF
12     CONTINUE
      IF (KK.NE.(MA+1)) PAUSE 'Improper set in LISTA'
      DO 14 J=1,MFIT
       DO 13 K=1,MFIT
        COVAR(J,K)=ZERO
13      CONTINUE
       BETA(J)=ZERO
14    CONTINUE
      DO 18 I=1,NDATA
       YM=Y(I)
       IF(MFIT.LT.MA) THEN
cpas: this has to be a bug, A() is not preset here
        STOP 'bad egg'
        DO 15 J=MFIT+1,MA
         YM=YM-A(LISTA(J))*P(LISTA(J),I)
15       CONTINUE
       ENDIF
       SIG2I=one/SIG(I)**2
       DO 17 J=1,MFIT
        WT=P(LISTA(J),I)*SIG2I
        DO 16 K=1,J
         COVAR(J,K)=COVAR(J,K)+WT*P(LISTA(K),I)
16       CONTINUE
        BETA(J)=BETA(J)+YM*WT
17      CONTINUE
18     CONTINUE
      IF (MFIT.GT.1) THEN
       DO 21 J=2,MFIT
        DO 19 K=1,J-1
         COVAR(K,J)=COVAR(J,K)
19       CONTINUE
21      CONTINUE
      ENDIF
C
      CALL GAUSSJ(COVAR,MFIT,NCVM,BETA,1,1)
C
      DO 22 J=1,MFIT
       A(LISTA(J))=BETA(J)
22     CONTINUE
      CHISQ=ZERO
      DO 24 I=1,NDATA
       SUM=ZERO
       DO 23 J=1,MA
        SUM=SUM+A(J)*P(J,I)
23      CONTINUE
       CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
24     CONTINUE
C
      CALL COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
c
      RETURN
      END
C_______________________________________________________________________
C
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXN=1000,MAXM=10,MAXM1=MAXM+1)
      PARAMETER (NMAX=MAXM1)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DATA ZERO,ONE/0.D0,1.D0/
      DO 11 J=1,N
       IPIV(J)=0
11    CONTINUE
C
      DO 22 I=1,N
C
       BIG=ZERO
       DO 13 J=1,N
        IF(IPIV(J).NE.1)THEN
         DO 12 K=1,N
          IF (IPIV(K).EQ.0) THEN
           IF (ABS(A(J,K)).GE.BIG)THEN
            BIG=ABS(A(J,K))
            IROW=J
            ICOL=K
           ENDIF
          ELSE IF (IPIV(K).GT.1) THEN
           STOP 'Singular matrix'
          ENDIF
12        CONTINUE
        ENDIF
13      CONTINUE
C
       IPIV(ICOL)=IPIV(ICOL)+1
       IF (IROW.NE.ICOL) THEN
        DO 14 L=1,N
         DUM=A(IROW,L)
         A(IROW,L)=A(ICOL,L)
         A(ICOL,L)=DUM
14       CONTINUE
        DO 15 L=1,M
         DUM=B(IROW,L)
         B(IROW,L)=B(ICOL,L)
         B(ICOL,L)=DUM
15       CONTINUE
       ENDIF
       INDXR(I)=IROW
       INDXC(I)=ICOL
       IF (A(ICOL,ICOL).EQ.ZERO) STOP 'Singular matrix.'
C
       PIVINV=ONE/A(ICOL,ICOL)
       A(ICOL,ICOL)=ONE
       DO 16 L=1,N
        A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
       DO 17 L=1,M
        B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
C
       DO 21 LL=1,N
        IF(LL.NE.ICOL)THEN
         DUM=A(LL,ICOL)
         A(LL,ICOL)=ZERO
         DO 18 L=1,N
          A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18        CONTINUE
         DO 19 L=1,M
          B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19        CONTINUE
        ENDIF
21      CONTINUE
C
22     CONTINUE
C
      DO 24 L=N,1,-1
       IF(INDXR(L).NE.INDXC(L))THEN
        DO 23 K=1,N
         DUM=A(K,INDXR(L))
         A(K,INDXR(L))=A(K,INDXC(L))
         A(K,INDXC(L))=DUM
23       CONTINUE
       ENDIF
24    CONTINUE
C
      RETURN
      END
C_______________________________________________________________________
C
      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COVAR(NCVM,NCVM),LISTA(MFIT)
      DATA ZERO/0.D0/
      DO 12 J=1,MA-1
       DO 11 I=J+1,MA
        COVAR(I,J)=ZERO
11      CONTINUE
12     CONTINUE
      DO 14 I=1,MFIT-1
       DO 13 J=I+1,MFIT
        IF(LISTA(J).GT.LISTA(I)) THEN
         COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
        ELSE
         COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
        ENDIF
13      CONTINUE
14     CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
       COVAR(1,J)=COVAR(J,J)
       COVAR(J,J)=ZERO
15     CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
       COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16     CONTINUE
      DO 18 J=2,MA
       DO 17 I=1,J-1
        COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18     CONTINUE
      RETURN
      END
