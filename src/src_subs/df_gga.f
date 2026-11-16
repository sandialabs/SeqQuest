c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DF_GGA
c
c
c---------------------------------------------------------------
c Purpose: combine under one roof all gga functional code
c
c Adapted: Peter A. Schultz, 28-December-2000, for SeqQuest 2.46
c
c Original code from official version of PBE kept by Kieron Burke:
c http://crab.rutgers.edu/~dft/pubs/PBE.asc
c
c Revision history:
c  13Nov12-JEM/2.65: HSE semilocal term added
c  17Dec08-RJM/2.62: AM05 spin added
c  24May06-AEM/2.60: AM05 added
c  11Oct03-PAS/2.57: fix pw91x stress terms
c   8Oct03-AEM/2.57: stress corrected for spin. argument eliminated.
c   7Oct03-AEM/2.57: stress added to blypsp. Non spin lyp eliminated.
c   2Jul03-AEM/2.56: blypsp added
c  11Jun03-AEM/2.55: blyp added (with stress)
c   8Apr02-PAS/2.53: gga (unpolarized) stress
c  28Dec00-PAS/2.46: adapted for Quest - routine names
c---------------------------------------------------------------
c
c Notes:
c  I substituted calls to the "slimmed" version of routine "GCOR"
c  in the PW91 potential, as had been done for the PBE potential
c  in the official version of the PBE code.
c  I also made changes in routine names to standardize everything
c  All modifications are denoted by "cpas"
c   - 28Dec00-PAS
c
c  I have implemeted blyp directly from Becke PRA 38, 3098 (1988)
c  and Lee, Yang, and Parr PRB 37, 785 (1988). All errors are mine.
c   - 11Jun03-AEM
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
cpas      SUBROUTINE CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
      SUBROUTINE DF_LSDC(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
C  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
C  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
C  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
C     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
C  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(gam=0.5198421D0,fzz=1.709921D0)
      parameter(thrd=0.333333333333D0,thrd4=1.333333333333D0)
      F = ((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
cpas: added rtRS, and replaced calls to gcor with calls to df_gcor/28dec00
      rtRS = DSQRT( RS )
cpas: I disapprove of passing explicit constants, but ...
cpas      CALL GCOR(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
cpas     1    0.49294D0,1.00D0,RS,EU,EURS)
      CALL DF_GCOR(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     1    0.49294D0,   rtRS,EU,EURS)
cpas      CALL GCOR(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
cpas     1    0.62517D0,1.00D0,RS,EP,EPRS)
      CALL DF_GCOR(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     1    0.62517D0,   rtRS,EP,EPRS)
cpas      CALL GCOR(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
cpas     1    0.49671D0,1.00D0,RS,ALFM,ALFRSM)
      CALL DF_GCOR(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     1    0.49671D0,   rtRS,ALFM,ALFRSM)
C  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZET**4
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
C  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
cpas      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
      SUBROUTINE DF_GCOR(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
c slimmed down version of GCOR used in PW91 routines, to interpolate
c LSD correlation energy, as given by (10) of
c J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c K. Burke, May 11, 1996.
      IMPLICIT REAL*8 (A-H,O-Z)
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DF_GGA
c
c
cpas      SUBROUTINE EXCHPBE(rho,S,U,V,lgga,lpot,EX,VX)
      SUBROUTINE DF_PBEX( rho, S,U,V, lgga,lpot, EX,VX, dexs )
c----------------------------------------------------------------------
C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c  K Burke's modification of PW91 codes, May 14, 1996
c  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  INPUT rho : DENSITY
C  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
C  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
c   (for U,V, see PW86(24))
c  input lgga:  (=0=>don't put in gradient corrections, just LDA)
c  input lpot:  (=0=>don't get potential and don't need U and V)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (E).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c   	e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c	e_x[PBE]=e_x[unif]*FxPBE(s)
c	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c uk, ul defined after [a](13)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density
      exunif = AX*rho**THRD
      if(lgga.eq.0)then
        ex=exunif
        vx=ex*thrd4
        return
      endif
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct PBE enhancement factor
      S2 = S*S
      P0=1.d0+ul*S2
      FxPBE = 1d0+uk-uk/P0
      EX = exunif*FxPBE
      if(lpot.eq.0)return
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE POTENTIAL:
c  find first and second derivatives of Fx w.r.t s.
c  Fs=(1/s)*d FxPBE/ ds
c  Fss=d Fs/ds
      Fs=2.d0*uk*ul/(P0*P0)
      Fss=-4.d0*ul*S*Fs/P0
cpas: stress term dexs = (1/s) * d[Ex]/ds = exunif * (1/s) * dFxBPE/ds
      dexs = exunif*Fs
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c calculate potential from [b](24)
      VX = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
cpas      SUBROUTINE CORPBE(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn,
      SUBROUTINE DF_PBEC(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn,
     1                  H,DVCUP,DVCDN, dht )
c----------------------------------------------------------------------
c  Official PBE correlation code. K. Burke, May 14, 1996.
C  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
C       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
C       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
C       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
C       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
C       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
c       :  UU,VV,WW, only needed for PBE potential
c       : lgga=flag to do gga (0=>LSD only)
c       : lpot=flag to do potential (0=>energy only)
c  output: ec=lsd correlation energy from [a]
c        : vcup=lsd up correlation potential
c        : vcdn=lsd dn correlation potential
c        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
c        : dvcup=nonlocal correction to vcup
c        : dvcdn=nonlocal correction to vcdn
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof,
c     {\sl Generalized gradient approximation made simple}, sub.
c     to Phys. Rev.Lett. May 1996.
c [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
c     construction of a generalized gradient approximation:  The PW91
c     density functional}, submitted to Phys. Rev. B, Feb. 1996.
c [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
c thrd*=various multiples of 1/3
c numbers for use in LSD energy spin-interpolation formula, [c](9).
c      GAM= 2^(4/3)-2
c      FZZ=f''(0)= 8/(9*GAM)
c numbers for construction of PBE
c      gamma=(1-log(2))/pi^2
c      bet=coefficient in gradient expansion for correlation, [a](4).
c      eta=small number to stop d phi/ dzeta from blowing up at
c          |zeta|=1.
      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(GAM=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*GAM))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gamma)
      parameter(eta=1.d-12)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c find LSD energy contributions, using [c](10) and Table I[c].
c EU=unpolarized LSD correlation energy
c EURS=dEU/drs
c EP=fully polarized LSD correlation energy
c EPRS=dEP/drs
c ALFM=-spin stiffness, [c](3).
c ALFRSM=-dalpha/drs
c F=spin-scaling factor from [c](9).
c construct ec, using [c](8)
      rtrs=dsqrt(rs)
cpas: I still disapprove of passing explicit constants, but ...
cpas      CALL gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
      CALL DF_GCOR(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     1    0.49294D0,rtrs,EU,EURS)
cpas      CALL gcor2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
      CALL DF_GCOR(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     1    0.62517D0,rtRS,EP,EPRS)
cpas      CALL gcor2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
      CALL DF_GCOR(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     1    0.49671D0,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c LSD potential from [c](A1)
c ECRS = dEc/drs [c](A2)
c ECZET=dEc/dzeta [c](A3)
c FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      if(lgga.eq.0)return
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c PBE correlation energy
c G=phi(zeta), given after [a](3)
c DELT=bet/gamma
c B=A of [a](8)
      G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
cpas: to refactor, and add compute of stress term
cpas      H = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      hfac = G3*(BET/DELT)
      hargln = 1.D0+DELT*Q4*T2/Q5
c
c gradient correction to energy density:
      H = hfac *DLOG( hargln )
c
      if( lpot.eq.0 ) RETURN
c
c gradient derivative for stress:
      Q7 = 1.d0 + 2.d0*B*T2
c     Stress term: dht = (1/t) * dH[t]/dt
      dht = 2.d0*hfac*(DELT/hargln) * ( Q7 / Q5**2 )
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3.D0
      GZ=(((1.d0+zet)**2+eta)**sixthm-
     1((1.d0-zet)**2+eta)**sixthm)/3.d0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      hRST = RSTHRD*T2*hBT*BEC*ECRS
      hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2.d0*BET*G3*Q9/Q8
      hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE DF_HSEX_HELPER(S,A,F,X,DXDS,DXDA,DXDF,D2XDS2,D2XDA2,
     1                          D2XDF2,D2XDSA,D2XDSF,D2XDAF)
c----------------------------------------------------------------------
c  COMPUTES THE EXCHANGE HOLE INTEGRAL THAT IS CENTRAL TO HSE EXCHANGE
c----------------------------------------------------------------------
c  INPUT S: same as in DF_HSEX (S>=0)
c  INPUT A: exponent of gaussian envelope in the integral (A>0)
c  INPUT F: parameter that defines the exchange hole (F>=0)
c  OUTPUT X: exchange hole integral
c  OUTPUT DXDS: derivative of X with S
c  OUTPUT DXDA: derivative of X with A
c  OUTPUT DXDF: derivative of X with F
c  OUTPUT D2XDS2: 2nd derivative of X with S
c  OUTPUT D2XDA2: 2nd derivative of X with A
c  OUTPUT D2XDF2: 2nd derivative of X with F
c  OUTPUT D2XDSA: derivative of X with S & A
c  OUTPUT D2XDSF: derivative of X with S & F
c  OUTPUT D2XDAF: derivative of X with A & F
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(euler=0.5772156649015329d0,c1=0.4516064d0)
      parameter(c2=0.57786348d0,c3=0.16520372d0,c4=0.0068635965d0)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c rescale the S input to avoid the undefined large-S regime of HSE
      if(S.le.8.d0) then
        S2 = S
        DS2DS = 1.d0
        D2S2DS2 = 0.d0
      else
        S2 = S-DEXP(DLOG(S-8.d0)-3.d0/(S-8.d0))
        DS2DS = 1.d0-(1.d0/(S-8.d0)+3.d0/(S-8.d0)**2)
     1          *DEXP(DLOG(S-8.d0)-3.d0/(S-8.d0))
        D2S2DS2 = (1.d0/(S-8.d0)**2+6.d0/(S-8.d0)**3
     1            -(1.d0/(S-8.d0)+3.d0/(S-8.d0)**2)**2)
     2            *DEXP(DLOG(S-8.d0)-3.d0/(S-8.d0))
      endif
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c compute the set of redundant terms that appear in several places
      r1 = DSQRT(c1)*DEXP(F/c1)*DERFC(DSQRT(F/c1))
      Dr1DF = r1/c1-DSQRT(1.0/(pi*F))
      r2 = 1.d0/(c2+F+A)
      r3 = (c4+F*2.d0/45.d0+S2**2*4.d0/1215.d0)*(c2+F+5.d0*A)
     1      +(8.d0/15.d0)*(2.d0*c1*(c2+F)**3-c3*(c2+F)**2)
     2      +(16.d0*DSQRT(pi)/15.d0)*(2.d0/3.d0-r1)
     3       *(c2+F)**3.5d0
      Dr3DF = (2.d0/45.d0)*(c2+F+5.d0*A)+c4+F*2.d0/45.d0
     1        +S2**2*4.d0/1215.d0
     2        +(8.d0/15.d0)*(6.d0*c1*(c2+F)**2-2.d0*c3*(c2+F))
     3        +(56.d0*DSQRT(pi)/15.d0)*(2.d0/3.d0-r1)*(c2+F)**2.5d0
     4        -(16.d0*DSQRT(pi)/15.d0)*Dr1DF*(c2+F)**3.5d0
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c special function is ln(1+(c2/c1)/y) - exp(y)*E1(y)
      y = (A+F)/c1
      SP = 0.d0
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c power series expansion for y <= 2.0
      if(y.le.2.d0) then
        do i = 22, 1, -1
          SP = -(y/DBLE(i))*(1.d0/DBLE(i) + SP)
        enddo
        SP = (DEXP(y)-1.d0)*DLOG(y)+DLOG(y+c2/c1)+DEXP(y)*(euler+SP)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c continued fraction expansion for y > 2.0
      else
        do i = 45, 1, -1
          SP = DBLE(i*i)/(y+1.d0+2.d0*DBLE(i) - SP)
        enddo
        SP = DLOG(1.d0+(c2/c1)/y) - 1.d0/(y+1.d0-SP)
      endif
      dSP = SP+1.d0/(y+c2/c1)-DLOG(1.d0+(c2/c1)/y)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c value of the integral
      X = c1*SP + c3*r2 + r3*r2**3
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c first derivatives
      DXDS = S2*DS2DS*(8.d0/1215.d0)*(c2+F+5.d0*A)*r2**3
      DXDy = dSP-c3*r2**2-3.d0*r3*r2**4
      DXDA = DXDy+5.d0*(c4+F*2.d0/45.d0+S2**2*4.d0/1215.d0)*r2**3
      DXDF = DXDy+Dr3DF*r2**3
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c second derivatives
      D2XDS2 = (DS2DS**2+S2*D2S2DS2)*(8.d0/1215.d0)*(c2+F+5.d0*A)*r2**3
      D2XDSA = -3.d0*DXDS*r2+S2*DS2DS*(8.d0/243.d0)*r2**3
      D2XDSF = -3.d0*DXDS*r2+S2*DS2DS*(8.d0/1215.d0)*r2**3

      D2XDy2 = (dSP-1.d0/(y+c2/c1)**2+1.d0/(y+(c1/c2)*y**2))/c1
     1         +2.d0*c3*r2**3+12.d0*r3*r2**5
      D2XDA2 = D2XDy2-30.d0*(c4+F*2.d0/45.d0+S2**2*4.d0/1215.d0)*r2**4
      D2XDAF = D2XDy2-15.d0*(c4+F*2.d0/45.d0+S2**2*4.d0/1215.d0)*r2**4
     1         -3.d0*Dr3DF*r2**4+(2.d0/9.d0)*r2**3
      D2XDF2 = D2XDy2-6.d0*Dr3DF*r2**4
     1         +(4.d0/45.d0+(8.d0/15.d0)*(12.d0*c1*(c2+F)-2.d0*c3)
     2         +(28.d0*DSQRT(pi)/3.d0)*(2.d0/3.d0-r1)*(c2+F)**1.5d0
     3       +(56.d0/15.d0)*(F**(-0.5d0)-DSQRT(pi)*r1/c1)*(c2+F)**2.5d0
     4         -(8.d0/15.d0)*(c2+F)**3.5d0*F**(-1.5d0)
     5         -(16.d0*DSQRT(pi)/15.d0)*((c2+F)/c1+3.5d0)*(c2+F)**2.5d0
     6         *Dr1DF)*r2**3
c----------------------------------------------------------------------
      RETURN
      END
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      SUBROUTINE DF_HSEX( rho, S,U,V, omega, EX,VX, dexs )
c----------------------------------------------------------------------
C  HSE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c----------------------------------------------------------------------
C  INPUT rho : DENSITY
C  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
C  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
C  INPUT V:  (LAPLACIAN rho)/(rho*(2*KF)**2)
c  INPUT omega:  (range separation parameter, =0.11 in HSE06)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX), POTENTIAL (VX),
C           and stress term dexs
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a]J. E. Moussa and P. A. Schultz, J. Comp. Chem. ??, ??? (2012).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 w(26),a(26)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(um=0.2195149727645171d0,uk=0.804d0)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c coefficients of erfc(x) ~= w(i) exp(-a(i)*x**2) fit
      w(1) = 2.6444678d-1
      w(2) = 2.1518752d-1
      w(3) = 1.5602599d-1
      w(4) = 1.0941077d-1
      w(5) = 7.6737994d-2
      w(6) = 5.4113816d-2
      w(7) = 3.8233937d-2
      w(8) = 2.6951074d-2
      w(9) = 1.8891822d-2
      w(10) = 1.3138701d-2
      w(11) = 9.0507222d-3
      w(12) = 6.1667803d-3
      w(13) = 4.1504553d-3
      w(14) = 2.7553333d-3
      w(15) = 1.8013008d-3
      w(16) = 1.1574260d-3
      w(17) = 7.2924425d-4
      w(18) = 4.4921162d-4
      w(19) = 2.6952465d-4
      w(20) = 1.5673828d-4
      w(21) = 8.7754201d-5
      w(22) = 4.6852295d-5
      w(23) = 2.3512291d-5
      w(24) = 1.0830436d-5
      w(25) = 4.3824327d-6
      w(26) = 1.4089207d-6
      a(1) = 1.0461980d0
      a(2) = 1.4529966d0
      a(3) = 2.4816609d0
      a(4) = 4.6668302d0
      a(5) = 9.1647121d0
      a(6) = 1.8453344d1
      a(7) = 3.7911389d1
      a(8) = 7.9437595d1
      a(9) = 1.6994660d2
      a(10) = 3.7184035d2
      a(11) = 8.3377708d2
      a(12) = 1.9206304d3
      a(13) = 4.5580211d3
      a(14) = 1.1182009d4
      a(15) = 2.8473892d4
      a(16) = 7.5635070d4
      a(17) = 2.1088326d5
      a(18) = 6.2204592d5
      a(19) = 1.9611688d6
      a(20) = 6.7003372d6
      a(21) = 2.5288181d7
      a(22) = 1.0845899d8
      a(23) = 5.5279640d8
      a(24) = 3.6218918d9
      a(25) = 3.5946468d10
      a(26) = 8.7787680d11
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c use Newton's method to find the value of the HSE fitting parameter
      if(S .gt. 1d16) then
        S = 1d16
      endif
      if(S .lt. 0.0002d0) then
c regularization of S: this perturbs the result less than the
c single-precision accuracy target & avoids resolving 0/0 cancellations
        S =   0.0002D0
        F =   6.387488276512606D-018
        DFDS =   2.613912033810124D-013
        D2FDS2 =   6.655331272050626D-009
      else
        F = (0.00979681d0*S**4 + 0.0410834d0*S**6) / (1.d0
     1       + 0.187440d0*S**4 + 0.00120824d0*S**5 + 0.0347188d0*S**6)
        call DF_HSEX_HELPER(S,0.d0,F,X,
     1        DXDS,DXDA,DXDF,D2XDS2,D2XDA2,D2XDF2,D2XDSA,D2XDSF,D2XDAF)
        Rnew = 1.d0+uk-uk/(1.d0+(um/uk)*S**2)-X
        Rold = DABS(Rnew) + 1.d0
        do while(DABS(Rnew) .lt. DABS(Rold))
          Rold = Rnew
          dF = Rold/DXDF
cxxxhse       if(F+dF > 0.d0) then
          if( (F+dF) .gt. 0.d0 )then
            F = F+dF
          else
            F = 0.5d0*F
          endif
          call DF_HSEX_HELPER(S,0.d0,F,X,
     1        DXDS,DXDA,DXDF,D2XDS2,D2XDA2,D2XDF2,D2XDSA,D2XDSF,D2XDAF)
          Rnew = 1.d0+uk-uk/(1.d0+(um/uk)*S**2)-X
        enddo
        DX0DS = 2.d0*um*S/(1.d0+(um/uk)*S**2)**2
        D2X0DS2 = (2.d0*um-6.d0*(um*S)**2/uk)/(1.d0+(um/uk)*S**2)**3
        DFDS = (DX0DS-DXDS)/DXDF
        D2FDS2 = (D2X0DS2-D2XDS2-D2XDF2*DFDS**2-2.d0*D2XDSF*DFDS)/DXDF
      endif
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct HSE enhancement factor, X2, and derivatives
      X2 = 0.d0
      DX2DP = 0.d0
      DX2DS = 0.d0
      D2X2DPS = 0.d0
      D2X2DS2 = 0.d0
      fermiK = (3.d0*pi**2*rho)**(1.d0/3.d0)
      afact = (omega/fermiK)**2
      do i = 1, 26
        call DF_HSEX_HELPER(S,a(i)*afact,F,X,
     1        DXDS,DXDA,DXDF,D2XDS2,D2XDA2,D2XDF2,D2XDSA,D2XDSF,D2XDAF)
        X2 = X2 + w(i)*X
        DX2DP = DX2DP - w(i)*(2.d0/3.d0)*afact*a(i)*DXDA/rho
        DX2DS = DX2DS + w(i)*(DXDS+DXDF*DFDS)
        D2X2DPS = D2X2DPS - w(i)*(2.d0/3.d0)*afact*a(i)
     1                     *(D2XDSA+D2XDAF*DFDS)/rho
        D2X2DS2 = D2X2DS2 + w(i)*(D2XDS2+2.d0*D2XDSF*DFDS+DXDF*D2FDS2
     1                            +D2XDF2*DFDS**2)
      enddo
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA & HSE exchange energy density (and stress term)
      EXlda = -(0.75d0/pi)*fermiK
      dEXlda = EXlda/(3.d0*rho)
      EX = EXlda*X2
cpas: stress term dexs = (1/s) * d[Ex]/ds
      dexs = EXlda*DX2DS/S
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c convert derivatives to rho/sigma(=M) to simplify potential formula
      G1 = S*rho*2.d0*fermiK
      G2 = V*rho*(2.d0*fermiK)**2
      G3 = U*rho**2*(2.d0*fermiK)**3
      DSDP = -(4.d0/3.d0)*S/rho
      DSDG = 1.d0/(rho*2.d0*fermiK)
      D2SDPG = -(4.d0/3.d0)*DSDG/rho
      DEXDP = dEXlda*X2+EXlda*(DX2DP+DX2DS*DSDP)
      DEXDG = EXlda*DX2DS*DSDG
      D2EXDPG = dEXlda*DX2DS*DSDG+EXlda*(D2X2DPS*DSDG+D2X2DS2*DSDP*DSDG
     1         +DX2DS*D2SDPG)
      D2EXDG2 = EXlda*D2X2DS2*DSDG**2
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c potential formula from [PRB 54, 1568 (1996)]
      VX = EX+rho*DEXDP-G1*DEXDG-rho*G1*D2EXDPG-rho*DEXDG*G2/G1
     1    +(rho/G1)*(DEXDG/G1-D2EXDG2)*G3
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
cpas      SUBROUTINE EXCHPW91(D,S,U,V,EX,VX)
      SUBROUTINE DF_PW91X( D, S,U,V, EX,VX, dexs )
C  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
C  INPUT D : DENSITY
C  INPUT S:  ABS(GRAD D)/(2*KF*D)
C  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(a1=0.19645D0,a2=0.27430D0,a3=0.15084D0,a4=100.d0)
      parameter(ax=-0.7385588D0,a=7.7956D0,b1=0.004d0)
      parameter(thrd=0.333333333333D0,thrd4=1.33333333333D0)
c for Becke exchange, set a3=b1=0
      FAC = AX*D**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.D0/DSQRT(1.D0+A*A*S2)
      P1 = DLOG(A*S+1.D0/P0)
      P2 = DEXP(-A4*S2)
      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*F
c
caem: start of stress term additions
cpas: stress term dexs = (1/s) * d[Ex]/ds = FAC * (1/s) * dFx[s]/ds
      if( s .gt. 1.d-8 )then
cpas: There was a couple of sign errors here in the original code
        dfx1 = -( a1*P1 + 4.d0*b1*S3 + a1*a*s*P0 ) * F
        dfx2 = a1*P1 + a1*a*s*P0 + 2.d0*s*(a2-a3*P2) + 2.d0*a4*a3*S3*P2
        dfxs = P3 * ( dfx1 + dfx2 )
        dexs =  FAC * dfxs / s
      else
cpas:   dexs = a2 + a3 ! since I multiply by "s" later, can blank out:
        dexs = 0.d0
      endif
caem: end of stress term additions
c
C  LOCAL EXCHANGE OPTION
C     EX = FAC
C  ENERGY DONE. NOW THE POTENTIAL:
      P5 = B1*S2-(A2-A3*P2)
      P6 = A1*S*(P1+A*S*P0)
      P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
      FS = P3*(P3*P5*P6+P7)
      P8 = 2.D0*S*(B1-A3*A4*P2)
      P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
      P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
      P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
      FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
      VX = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)
C  LOCAL EXCHANGE OPTION:
C     VX = FAC*THRD4
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
cpas      SUBROUTINE CORpw91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H,
      SUBROUTINE DF_PW91C( RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H,
     1                   DVCUP,DVCDN, dht )
C  pw91 CORRELATION, modified by K. Burke to put all arguments
c  as variables in calling statement, rather than in common block
c  May, 1996.
C  INPUT RS: SEITZ RADIUS
C  INPUT ZET: RELATIVE SPIN POLARIZATION
C  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
C  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
C  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
C  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
C  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
C  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(xnu=15.75592D0,cc0=0.004235D0,cx=-0.001667212D0)
      parameter(alf=0.09D0)
      parameter(c1=0.002568D0,c2=0.023266D0,c3=7.389D-6,c4=8.723D0)
      parameter(c5=0.472D0,c6=7.389D-2,a4=100.D0)
      parameter(thrdm=-0.333333333333D0,thrd2=0.666666666667D0)
      BET = XNU*CC0
      DELT = 2.D0*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.D0+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = 0.663436444d0*rs
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.D0*CX/7.D0
      R2 = XNU*COEFF*G3
      R3 = DEXP(-R1*T2)
      H0 = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0 + H1
C  LOCAL CORRELATION OPTION:
C     H = 0.0D0
C  ENERGY DONE. NOW THE POTENTIAL:
      CCRS = (C2+2.*C3*RS)/Q7 - Q6*(C4+2.*C5*RS+3.*C6*RS2)/Q7**2
      RSTHRD = RS/3.D0
      R4 = RSTHRD*CCRS/COEFF
      GZ = ((1.D0+ZET)**THRDM - (1.D0-ZET)**THRDM)/3.D0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
cpasaem:start add terms for stress calculation
c stress terms: dht=(1/t)*dH[t]/dt
      dht = 2.d0*( G3*BET*Q9/Q8 + R2*R3*( 1.D0 - R1*T2 ) )
cpasaem:end
      H0B = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      H0RS = -RSTHRD*H0B*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      H0BT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
      H0Z = 3.D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
      H0T = 2.*BET*G3*Q9/Q8
      H0ZT = 3.D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      H0TT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      H1RS = R3*R2*T2*(-R4+R1*T2/3.D0)
      FACT4 = 2.D0-R1*T2
      H1RST = R3*R2*T2*(2.D0*R4*(1.D0-R1*T2)-THRD2*R1*T2*FACT4)
      H1Z = GZ*R3*R2*T2*(3.D0-4.D0*R1*T2)/G
      H1T = 2.D0*R3*R2*(1.D0-R1*T2)
      H1ZT = 2.D0*GZ*R3*R2*(3.D0-11.D0*R1*T2+4.D0*R1*R1*T4)/G
      H1TT = 4.D0*R3*R2*R1*T*(-2.D0+R1*T2)
      HRS = H0RS+H1RS
      HRST = H0RST+H1RST
      HT = H0T+H1T
      HTT = H0TT+H1TT
      HZ = H0Z+H1Z
      HZT = H0ZT+H1ZT
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
C  LOCAL CORRELATION OPTION:
C     DVCUP = 0.0D0
C     DVCDN = 0.0D0
      RETURN
      END
c
c
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      subroutine DF_B1X( rho, s, u, v, ex, vx, dexs )
c Becke exchange for a spin-unpolarized electron system
c A. D. Becke, PRA 38, 3098 (1988).
c Input quantities:
c   rho = total electron number density
c   s = ABS(GRAD D)/(2*KF*D)
c   u = (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
c   v = (LAPLACIAN D)/(D*(2*KF)**2)
c Output quantities:
c   ex = exchange energy per electron
c   vx = exchange potential (in rydbergs)
c   dexs = term needed in stress calculation. See GGASTRESS.
c
c Note: Becke's original form (with beta=0.0042 as parameter)
c   is translated into the PW91 form. This means:
c   c3 = 2 (6 pi^2)^(1/3), c2 = 6 beta c3, and
c   c1 = beta c3^2 / (3/4 (6/pi)^(1/3)). Note the similarity
c   to the PW91 parameters above (c2 to a1, c1 to a2, and c3 to a).
c   ex0 is the LDA prefactor: -3/4 (3/pi)^(1/3).
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER ( c1 = 0.2742931042d0, c2 = 0.1964479653d0 )
      PARAMETER ( c3 = 7.795554179d0, ex0 = -0.7385587664d0 )
      PARAMETER ( one = 1.d0, two = 2.d0, three = 3.d0, four = 4.d0 )
      PARAMETER ( five = 5.d0, thrd = one / three )
c
c Intermediate quantities formed.
c   fx = enhancement factor Fx(s)
c   vx0 = zeroth order prefactor in vx
c       = 3/4 Fx(s)
c   vxv = first order prefactor for v in vx
c       = -1/s dFx(s)/ds
c   vxu = first order prefactor for u in vx
c       = - d( 1/s dFx(s)/ds )/ds
c   vxs3 = first order prefactor for s^3 in vx
c        = 4/3 d( 1/s dFx(s)/ds )/ds
c   exlda = LDA exchange
c
      a = SQRT( one + (c3*s)**2 )
      b = LOG( c3*s + a )
      c = one + c2*s*b
      d = five + four*(c3*s)**2 - two*c2*c3*a*s**2
      e = c2*c3*s**2 + three*a**3
      fx = one + c1*s**2 / c
      vx0 = four*fx / three
      vxv = -c1*( two*a + c2*s*a*b - c2*c3*s**2 ) / ( a*c**2 )
      vxu = c1*c2*( c3*s*d + b*e + c2*s*a**3*b**2 ) / ( (a*c)**3 )
      vxs3 = -four*vxu / three
      exlda = ex0*rho**thrd
c
c Assemble ex, vx and dexs
c
      ex = exlda*fx
      vx = exlda*( vx0 + vxv*v + vxu*u + vxs3*s**3 )
      dexs = -exlda*vxv
c
c
c
c
      RETURN
      END
c
c
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      subroutine DF_LYP( rhup,rhdn, grho, grhup,grhdn, g2rhup,g2rhdn,
     $   ec, vcup, vcdn, dcgrho, dcgrup, dcgrdn, dcg2rup, dcg2rdn )
c LYP correlation for a spin-polarized electron system.
c Lee, Yang, and Parr, PRB 37, 785 (1988).
c Input quantities:
c   rhup   = spin-up electron density
c   rhdn   = spin-down electron density
c   grho   = | GRAD[rho] |
c   grhup  = | GRAD[rhup] |
c   grhdn  = | GRAD[rhdn] |
c   g2rhup = GRAD**2[rhup]
c   g2rhdn = GRAD**2[rhdn]
c Output quantities:
c   ec   = correlation energy per particle
c   vcup = spin-up correlation potential
c   vcdn = spin-down correlation potential
c   dcgrho = stress correction for gradients (grho)
c   dcgrup = spin-up stress correction for gradients (grhup)
c   dcgrdn = spin-down stress correction for gradients (grhdn)
c   dcg2rup = spin-up stress correction for laplacians
c   dcg2rdn = spin-down stress correction for laplacians
c
c Note: The constants a, b, c, d, have dimensions
c   Hartree, dimensionless, bohr and bohr
c   and are given in the original LYP paper.
c   cf = 3/10 (3 pi^2)^(2/3)
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER ( a = 0.04918d0, b = 0.132d0 )
      PARAMETER ( c = 0.2533d0, d = 0.349d0, cf = 2.871234 )
      PARAMETER ( one = 1.d0, two = 2.d0, three = 3.d0, four = 4.d0 )
      PARAMETER ( five = 5.d0, thrd = one / three, snth = -one / 9.d0 )
c
c Intermediate quantities are formed.
c      ec = ec0 + ecgg*grho**2 + ecggu*grhup**2 +
c           + ecggd*grhdn + eclu*g2rhup + ecld*g2rhdn
c      vca = vca0 + vcagg*grho**2 + vcaggu*grhup**2 +
c           + vcaggd*grhdn + vcalu*g2rhup + vcald*g2rhdn
c      where a is u(p) or d(own).
c      ec = fc/rho where fc is the integrand of Eq. (22) in the LYP paper.
c      Each of the three terms in the definition of vca,
c      vcup = d[fc]/d(rhup) - GRAD dot d[fc]/d(GRAD[rhup]) +
c             + GRAD**2 d[fc]/d(GRAD**2[rhup]),
c      contribute to the prefactors above.
c      In LYP the task of forming vcud and vcdn is somewhat simplified by simple
c      relations between all ec prefactors except ec0.
c      I use ec0 and ecgg as primary quantities.
c      Primary intermediate quantities that are formed are:
c      vcu0 and vcd0: spin up and spin down derivatives of rho*ec0
c      db0du and db0dd: spin up and spin down derivatives of rho*ecgg
c      d2b0du2, d2b0dd2, and d2b0dudd: second derivatives of rho*ecgg
c
c      Most quantities are spin independent, if not they are clearly marked
c      with u(p) or d(own) in the parameter name.
c
      rho = rhup + rhdn
      r13 = rho**thrd
      r23 = r13*r13
      r33 = rho
      r43 = rho*r13
      r53 = r43*r13
      ru83 = rhup**(8.d0*thrd)
      rd83 = rhdn**(8.d0*thrd)
      rurdr = ( ru83 + rd83 ) / r43**2
      paru = 12.d0*rhup + 9.d0*rhdn
      pard = 12.d0*rhdn + 9.d0*rhup
      aa = d + r13
      argexp = -c / r13
      if ( argexp .gt. -60.d0 ) then
        bb = a*b*exp( argexp )
      else
        bb = 0.d0
      endif
      a0 = two**( five*thrd ) * cf * bb
      b0 = bb / ( aa * r13 * rho**4 )
c
c Energy density prefactors formed
c     In lyp ecggu=ecggd=ecggs
c
      ec0 = -four*rhup*rhdn*( a0*rurdr + a ) / ( aa * r53 )
      ecgg = rhup*rhdn*b0
      ecggs = snth*ecgg
      eclu = ecggs*paru
      ecld = ecggs*pard
c
c vcu0 and vcd0 formed
c
      ceven = ru83 + rd83
      codd = rhdn*ru83 - rhup*rd83
      ccu = ( c*rhup + three*rhdn*r13 )*ceven + 8.d0*r13*codd
      ccd = ( c*rhdn + three*rhup*r13 )*ceven - 8.d0*r13*codd
      ddu = a0*( aa*ccu + d*rhup*r13*ceven ) +
     $             a*rho**3*(d*rhup + three*rhdn*aa)
      ddd = a0*( aa*ccd + d*rhdn*r13*ceven ) +
     $             a*rho**3*(d*rhdn + three*rhup*aa)
      vcu0 = -four*rhdn*ddu*r13 / ( three * aa**2 * rho**5 )
      vcd0 = -four*rhup*ddd*r13 / ( three * aa**2 * rho**5 )
c
c Mess, sorry no way to make it nice
c
      ee = rhup*rhdn*( c*aa / r13 - 7.d0*aa - r13 )
      db0du = ( ee + three * rhdn**2 * aa ) * b0 / ( aa * three )
      db0dd = ( ee + three * rhup**2 * aa ) * b0 / ( aa * three )
      ff0 = r13*aa*( 30.d0*aa + three*r13 )
      ff1 = r13*( 35.d0*aa**2 + 8.d0*aa*r13 + r23 )
      ff2 = c*aa*three*aa
      ff3 = c*aa*( 9.d0*aa + r13 )
      ff5 = two*r13*( ff0 - ff2 )
      ff6 = ( c*aa )**2 + two*r13*( ff1 - ff3 )
      ffu = rhdn*ff5 - rhup*ff6
      ffd = rhup*ff5 - rhdn*ff6
      d2b0du2 = snth*bb*ffu*rhdn / ( rho**6 * aa**3 )
      d2b0dd2 = snth*bb*ffd*rhup / ( rho**6 * aa**3 )
      sq1 = three*( rhup**2 - 6.d0*rhup*rhdn + rhdn**2 )
      ff7 = ( 19.d0 * aa**2 ) + d*r13 + c*aa
      ff8 = ( c*aa )**2 - two*r23*ff7
      ff9 = r13*aa*( 7.d0*r13*aa - c*aa + r23 )
      ffud = rhup*rhdn*ff8 - sq1*ff9
      d2b0dudd =-snth*bb*ffud / ( rho**6 * aa**3 )
c
c Remaining potential prefactors formed
c
      vcugg = snth*(d2b0dudd*paru + 9.d0*db0du + 20.d0*db0dd)
      vcdgg = snth*(d2b0dudd*pard + 9.d0*db0dd + 20.d0*db0du)
      vcuggu = snth*((d2b0du2-d2b0dudd)*paru+23.d0*db0du-20.d0*db0dd)
      vcdggu = snth*((d2b0du2-d2b0dudd)*pard-17.d0*db0dd+16.d0*db0du)
      vcuggd = snth*((d2b0dd2-d2b0dudd)*paru-17.d0*db0du+16.d0*db0dd)
      vcdggd = snth*((d2b0dd2-d2b0dudd)*pard+23.d0*db0dd-20.d0*db0du)
      vculu = snth*(two*db0du*paru + 40.d0*rho*ecgg)
      vcdlu = snth*(db0dd*paru + db0du*pard) - four*rho*ecgg
      vculd = snth*(db0du*pard + db0dd*paru) - four*rho*ecgg
      vcdld = snth*(two*db0dd*pard + 40.d0*rho*ecgg)
c
c Final quantities.
c
      ec = ec0 + ecgg * grho**2 + ecggs * ( grhup**2 + grhdn**2 ) +
     $     eclu * g2rhup + ecld * g2rhdn
      vcup = vcu0 + vcugg * grho**2 +
     $       vcuggu * grhup**2 + vcuggd * grhdn**2 +
     $       vculu * g2rhup + vculd * g2rhdn
      vcdn = vcd0 + vcdgg * grho**2 +
     $       vcdggu * grhup**2 + vcdggd * grhdn**2 +
     $       vcdlu * g2rhup + vcdld * g2rhdn
      dcgrho = two * ecgg * rho
      dcgrup = two * ecggs * rho
      dcgrdn = two * ecggs * rho
      dcg2rup = -two * eclu * rho
      dcg2rdn = -two * ecld * rho
c
c
c
      RETURN
      END
c
c
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      subroutine DF_AM05(n,s,u,v,ex,ec,vx,vc,dexs,decs)
c
c Armiento Mattsson am05 functional for exchange and correlation 
c "R. Armiento and A. E. Mattsson, PRB 72, 085108 (2005). " 
c input
c   n = electron density, D [bohr**(-3)]
c   s = ABS(GRAD D)/(2*KF*D) [dimless]
c   u = (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3) [dimless]
c   v = (LAPLACIAN D)/(D*(2*KF)**2) [dimless]
c output
c   ex       exchange energy per electron [hartree]
c   ec       correlation energy per electron [hartree]
c   vx       exchange potential [hartree]
c   vc       correlation potential [hartree]
c   dexs     exchange stress term. See GGASTRESS.
c   decs     correlation stress term. See GGASTRESS.
c
      implicit none
c
c     ** Input parameters
      real*8 n, s, u, v
c
c     ** Output parameters
      real*8 ex, ec, vx, vc, dexs, decs
c
c     ** Constants and shorthand notation
      real*8 pi, g, a, s2, c
      parameter (pi = 3.141592653589793238462643383279502884197d0)
      parameter (g = 0.8098d0, a = 2.804d0)
      parameter (c = 0.7168d0)
c
c     ** Local variables
      real*8 exlda, vxlda, eclda, vclda, X, Xsos, Xsossos
      real*8 Hx, Hxsos, Hxsossos, Hc, Hcsos, Hcsossos
      real*8 F, Fsos, s2Fsossos
      real*8 szsoz, mixder
      real*8 denom, denomsos, sdenomsoss
      real*8 zfac, zosn, w
c
c     ** Cutoff 
      if((n .le. 1.0d-16)) then
         ex = 0.0d0
         ec = 0.0d0
         vx = 0.0d0
         vc = 0.0d0
         return
      endif
c
      s2 = s**2
c
c     *******************
c       LDA correlation
c     *******************
      call am05_ldapwc(n,eclda,vclda)
c
c     *******************
c        LDA exchange
c     *******************
      vxlda = -(3.0d0*n/pi)**(1.0d0/3.0d0)
      exlda = (3.0d0/4.0d0*vxlda)
c      
c     ********************
c          Exchange
c     ********************
c     ** Interpolation index
      X = 1.0d0/(1.0d0 + a*s2)
c      
c     ** Airy LAA refinement function
      call am05_labertw(s**(3.0d0/2.0d0)/sqrt(24.0d0),w)
c
c     ** am05_lambertw give back argument if it is < 1.0e-20
c     ** (1.0e-14)^{3/2} = 1.0e-21 => give  low s limit for z/s
c     ** zosn = normalized z/s
      if (s .lt. 1.0e-14) then
              zosn = 1.0d0
      else
              zosn = 24.0d0**(1.0d0/3.0d0)*w**(2.0d0/3.0d0)/s
      end if
      zfac = s2*(zosn*27.0d0/32.0d0/pi**2)**2
c
c     ** denom = denominator of Airy LAA refinement function
      denom = 1.0d0 + c*s2*zosn*(1.0d0 + zfac)**(1.0d0/4.0d0)
      F = (c*s2 + 1.0d0)/denom
c      
c     ** Exchange refinement function
      Hx = X + (1.0d0 - X)*F
c
c     ** Exchange energy per particle, Ex = Integrate[n*ex]
      ex = exlda*Hx
c      
c     ********************
c         Correlation
c     ********************
c     ** Correlation refinement function 
      Hc = X + g*(1.0d0 - X)
c      
c     ** Correlation energy per particle, Ec = Integrate[n*ec]
      ec = eclda*Hc
c      
c     ***************************
c         Exchange potential
c     ***************************
c     ** Interpolation index derivatives: 1/s dX/ds, 1/s d/ds(1/s dX/ds)                   
      Xsos = -2.0d0*a*X**2
      Xsossos = 8.0d0*a**2*X**3
c
c     ** Airy LAA refinement function derivatives, 1/s dF/ds and s^2 1/s d/ds (1/s dF/ds) 
c     ** szsoz = s*(dz/ds)/z, mixder = szsoz + s^2*(d^2z/ds^2)/z
      szsoz = 1.0d0/(1.0d0 + w)
      mixder = (2.0d0 - w)/(2.0d0*(1.0d0 + w)**3)
c
c     ** denomsos = 1/s d(denom)/ds, sdenomsoss = s*d/ds(1/s d(denom)/ds))
      denomsos = c*zosn/(1.0d0 + zfac)**(3.0d0/4.0d0)*
     *     (1.0d0 + zfac + (1.0d0 + 3.0d0/2.0d0*zfac)*szsoz)
      sdenomsoss = c*zosn/(1.0d0 + zfac)**(7.0d0/4.0d0)*
     *     (-1.0d0 - zfac*(2.0d0 + zfac)
     *      + (1.0d0 + zfac/2.0d0*(5.0d0 + 3.0d0*zfac))*mixder
     *      + 3.0d0/2.0d0*zfac*(1.0d0 + zfac/2.0d0)*szsoz**2)

      Fsos = c/denom**2*(2.0d0 - zosn*
     *            ((1.0d0 - c*s2)*(1.0d0 + zfac)**(1.0d0/4.0d0) +
     *             (1.0d0 + c*s2)*(1.0d0 + 3.0d0/2.0d0*zfac)/
     *                      (1.0d0 + zfac)**(3.0d0/4.0d0)*szsoz))
      s2Fsossos = (-4.0d0*c*s2*denom*denomsos + (c*s2 + 1.0d0)*
     *                    (2.0d0*s2*denomsos**2 - denom*sdenomsoss))/
     *                      denom**3
c
c     ** Refinement function derivatives, 1/s dHx/ds, 1/s d/ds (1/s dHx/ds) 
c     ** We use that (1 - X) = a*X*s2
      Hxsos = (1.0d0 - X)*Fsos - (F - 1.0d0)*Xsos
      Hxsossos = - 2.0d0*Fsos*Xsos + a*X*s2Fsossos -
     *                  (F - 1.0d0)*Xsossos
c
c     ** vx formula for gradient dependent functional,
      vx = vxlda*(Hx - s2*Hxsos) +
     *     exlda*((4.0d0/3.0d0*s2-v)*Hxsos - 
     *     (u-4.0d0/3.0d0*s**3)*s*Hxsossos)
c     **exchange stress term
      dexs = exlda*Hxsos
c
c     *****************************
c        Correlation potential
c     *****************************
c     ** Correlation refinement function derivatives, dF/ds 
      Hcsos = Xsos*(1.0d0 - g)
      Hcsossos = Xsossos*(1.0d0 - g)
c      
c     ** vc formula for gradient dependent functional,
c     Generalized form of Eq. (24) in PRB 33, 8800 (1986) 
      vc = vclda*(Hc - s2*Hcsos) +
     *     eclda*((4.0d0/3.0d0*s2 - v)*Hcsos - 
     *     (u - 4.0d0/3.0d0*s**3)*s*Hcsossos)
c     ** correlation stress term      
      decs = eclda*Hcsos
c
      return
      end
c
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
c I am repeating the LDA PW routine here even if it is available
c in DF_GCOR above. Once I clean up the PBE and PW91 implementations
c I'll eliminate DF_GCOR and replace it with something nicer. Until then
c I keep the am05 specific name since only DF_AM05 uses it.
c
      subroutine am05_ldapwc(n,ec,vc)
c***********************************************
c Local density approximation correlation
c
c input
c   n        electron density [bohr**(-3)]
c
c output
c   ec       correlation energy per electron [hartree]
c   vc       correlation potential [hartree]
c
c As parameterized by Perdew Wang,
c Phys. Rev. B 45, 13244 (1992) 
c Based on monte carlo data by Ceperley Alder, 
c Phys. Rev. Lett. 45, 566 (1980)
c
c (Clean room implementation from paper)
c
c***********************************************
      implicit none
c     ** Input parameters
      real*8 n
c
c     ** Output parameters
      real*8 ec, vc
c
c     ** Constants
      real*8 pi
      real*8 A0,a01,b01,b02,b03,b04
      parameter (pi = 3.141592653589793238462643383279502884197d0)
      parameter (a01 = 0.21370d0)
      parameter (b01 = 7.5957d0)
      parameter (b02 = 3.5876d0)
      parameter (b03 = 1.6382d0)
      parameter (b04 = 0.49294d0)
c     ** Paper actually use this:
c      parameter (A0 = 0.031091d0)
c     ** But routines now "defacto standard" was distributed using:
      parameter (A0 = 0.0310907d0)
c
c     ** Local variables
      real*8 rsq
      real*8 Q0, Q1, Q1p, ecrs
c
      rsq = (3.0d0/(4.0d0*pi*n))**(1.0d0/6.0d0)
c
      ec = -2.0d0*A0*(1.0d0 + a01*rsq**2)* 
     *     log(1.0d0 + 1.0d0/
     *     (2.0d0*A0*rsq*(b01 + rsq*(b02 + rsq*(b03 + b04*rsq)))))
c
      Q0 = -2.0d0*A0*(1.0d0 + a01*rsq**2)
      Q1 = 2.0d0*A0*rsq*(b01 + rsq*(b02 + rsq*(b03 + b04*rsq)))
      Q1p = A0*(b01/rsq+2.0d0*b02+3.0d0*b03*rsq+4.0d0*b04*rsq**2)
      ecrs = -2.0d0*A0*a01*log(1.0d0 + 1.0d0/Q1)-Q0*Q1p/(Q1**2+Q1)
c
      vc=ec - rsq**2/3.0d0*ecrs
c      
      return
      end
c
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
c This routine is only used by DF_AM05, and it will never be used by 
c any other routine. That's why I do not rename it according to QUEST
c standard.
c 
c 
c
      subroutine am05_labertw(z,result)
c     ***********************************************
c       LambertW function. 
c
c       Corless, Gonnet, Hare, Jeffrey, and Knuth (1996), 
c         Adv. in Comp. Math. 5(4):329-359. 
c       Implementation approach loosely inspired by the 
c       GNU Octave version by N. N. Schraudolph, but this 
c       implementation is only for real values and 
c       principal branch.
c
c       Copyright (c) 2005, Rickard Armiento
c       All rights reserved.
c     ***********************************************      
      implicit none
c
c     input
      real*8 z
c     output
      real*8 result
c     local variables
      real*8 e,t,p    
      integer i
c
c     ** If z too low, go with the first term of the power expansion, z
      if( z .lt. 1.0d-20) then
         result = z
         return
      endif
c
      e = exp(1.0d0)
c
c     ** Inital guess
      if( abs(z + 1.0d0/e) .gt. 1.45d0 ) then
c        ** Asymptotic expansion at 0 and Inf
         result = log(z)
         result = result - log(result)
      else
c        ** Series expansion about -1/e to first order
         result = 1.0d0*sqrt(2.0d0*e*z + 2.0d0) - 1.0d0
      endif
c
c     ** Find result through iteration
      do i=1,10
         p = exp(result)
         t = result*p - z
         if( result .ne. -1.0d0 ) then
            t = t/(p*(result + 1.0d0) - 
     *           0.5d0*(result + 2.0d0)*t/(result + 1.0d0))
         else
            t = 0.0d0
         endif
         result = result - t
         if(abs(t) .lt. (2.48d0*1.0d-14)*(1.0d0 + abs(result))) then
            return
         endif
      enddo
c     ** This should never happen!;
      write(*,*) 'am05_labertw: iteration limit reached.' 
      write(*,*) 'Should never happen: execution aborted.'
      stop
      end
c
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      subroutine DF_AM05s1(n,nop,s,u,v,ex,vx,dexs)
c
c Armiento Mattsson am05 functional for exchange and correlation 
c "R. Armiento and A. E. Mattsson, PRB 72, 085108 (2005). " 
C 9/17/8:RJM spin dependence added 
C     Note: subroutine outputs ex 
C     and some convenient mixture of vx and vc 
C
c input
c   n = electron density, D [bohr**(-3)]
c   s = ABS(GRAD D)/(2*KF*D) [dimless]
c   u = (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3) [dimless]
c   v = (LAPLACIAN D)/(D*(2*KF)**2) [dimless]
c output
c   ex       exchange energy per electron [hartree]
c   ec       correlation energy per electron [hartree]
c   vx       exchange potential [hartree]
c   vc       correlation potential [hartree]
c   dexs     exchange stress term. See GGASTRESS.
c   decs     correlation stress term. See GGASTRESS.
c
      implicit none
c
c     ** Input parameters
      real*8 n,nop, s, u, v
c
c     ** Output parameters
      real*8 ex, vx, dexs
c     ** Auxilliary parameters
      real*8 ec, vc, decs
c
c     ** Constants and shorthand notation
      real*8 pi, g, a, s2, c,two,third
      parameter (pi = 3.141592653589793238462643383279502884197d0)
      parameter (g = 0.8098d0, a = 2.804d0)
      parameter (c = 0.7168d0)
      parameter (two = 2.d0, third = 1d0/3d0)
c
c     ** Local variables
      real*8 exlda, vxlda, eclda, vclda, X, Xsos, Xsossos
      real*8 Hx, Hxsos, Hxsossos, Hc, Hcsos, Hcsossos
      real*8 F, Fsos, s2Fsossos
      real*8 szsoz, mixder
      real*8 denom, denomsos, sdenomsoss
      real*8 zfac, zosn, w
	  real*8 rs,zet,a1,a2,a3,a4,ntot	
c
c     ** Cutoff 
      if((n .le. 1.0d-16)) then
         ex = 0.0d0
         ec = 0.0d0
         vx = 0.0d0
         vc = 0.0d0
         return
      endif
c
      s2 = s**2
c
c     *******************
c       LDA correlation
c     *******************
      ntot=n/two+nop
      rs=(3./(4.*pi*ntot))**third
      zet=(n/two-nop)/ntot
      call DF_LSDC(rs,zet,eclda,vclda,a1,a2,a3,a4)
c     *******************
c        LDA exchange
c     *******************
      vxlda = -(3.0d0*n/pi)**(1.0d0/3.0d0)
      exlda = (3.0d0/4.0d0*vxlda)
c      
c     ********************
c          Exchange
c     ********************
c     ** Interpolation index
      X = 1.0d0/(1.0d0 + a*s2)
c      
c     ** Airy LAA refinement function
      call am05_labertw(s**(3.0d0/2.0d0)/sqrt(24.0d0),w)
c
c     ** am05_lambertw give back argument if it is < 1.0e-20
c     ** (1.0e-14)^{3/2} = 1.0e-21 => give  low s limit for z/s
c     ** zosn = normalized z/s
      if (s .lt. 1.0e-14) then
              zosn = 1.0d0
      else
              zosn = 24.0d0**(1.0d0/3.0d0)*w**(2.0d0/3.0d0)/s
      end if
      zfac = s2*(zosn*27.0d0/32.0d0/pi**2)**2
c
c     ** denom = denominator of Airy LAA refinement function
      denom = 1.0d0 + c*s2*zosn*(1.0d0 + zfac)**(1.0d0/4.0d0)
      F = (c*s2 + 1.0d0)/denom
c      
c     ** Exchange refinement function
      Hx = X + (1.0d0 - X)*F
c
c     ** Exchange energy per particle, Ex = Integrate[n*ex]
      ex = exlda*Hx
c      
c     ********************
c         Correlation
c     ********************
c     ** Correlation refinement function 
      Hc = X + g*(1.0d0 - X)
c      
c     ***************************
c         Exchange potential
c     ***************************
c     ** Interpolation index derivatives: 1/s dX/ds, 1/s d/ds(1/s dX/ds)                   
      Xsos = -2.0d0*a*X**2
      Xsossos = 8.0d0*a**2*X**3
c
c     ** Airy LAA refinement function derivatives, 1/s dF/ds and s^2 1/s d/ds (1/s dF/ds) 
c     ** szsoz = s*(dz/ds)/z, mixder = szsoz + s^2*(d^2z/ds^2)/z
      szsoz = 1.0d0/(1.0d0 + w)
      mixder = (2.0d0 - w)/(2.0d0*(1.0d0 + w)**3)
c
c     ** denomsos = 1/s d(denom)/ds, sdenomsoss = s*d/ds(1/s d(denom)/ds))
      denomsos = c*zosn/(1.0d0 + zfac)**(3.0d0/4.0d0)*
     *     (1.0d0 + zfac + (1.0d0 + 3.0d0/2.0d0*zfac)*szsoz)
      sdenomsoss = c*zosn/(1.0d0 + zfac)**(7.0d0/4.0d0)*
     *     (-1.0d0 - zfac*(2.0d0 + zfac)
     *      + (1.0d0 + zfac/2.0d0*(5.0d0 + 3.0d0*zfac))*mixder
     *      + 3.0d0/2.0d0*zfac*(1.0d0 + zfac/2.0d0)*szsoz**2)

      Fsos = c/denom**2*(2.0d0 - zosn*
     *            ((1.0d0 - c*s2)*(1.0d0 + zfac)**(1.0d0/4.0d0) +
     *             (1.0d0 + c*s2)*(1.0d0 + 3.0d0/2.0d0*zfac)/
     *                      (1.0d0 + zfac)**(3.0d0/4.0d0)*szsoz))
      s2Fsossos = (-4.0d0*c*s2*denom*denomsos + (c*s2 + 1.0d0)*
     *                    (2.0d0*s2*denomsos**2 - denom*sdenomsoss))/
     *                      denom**3
c
c     ** Refinement function derivatives, 1/s dHx/ds, 1/s d/ds (1/s dHx/ds) 
c     ** We use that (1 - X) = a*X*s2
      Hxsos = (1.0d0 - X)*Fsos - (F - 1.0d0)*Xsos
      Hxsossos = - 2.0d0*Fsos*Xsos + a*X*s2Fsossos -
     *                  (F - 1.0d0)*Xsossos
c
c     ** vx formula for gradient dependent functional,
      vx = vxlda*(Hx - s2*Hxsos) +
     *     exlda*((4.0d0/3.0d0*s2-v)*Hxsos - 
     *     (u-4.0d0/3.0d0*s**3)*s*Hxsossos)
c     **exchange stress term
      dexs = exlda*Hxsos
c
c     *****************************
c        Correlation potential
c     *****************************
c     ** Correlation refinement function derivatives, dF/ds 
      Hcsos = Xsos*(1.0d0 - g)
      Hcsossos = Xsossos*(1.0d0 - g)
c      
c     ** vc formula for gradient dependent functional,
c     Generalized form of Eq. (24) in PRB 33, 8800 (1986) 
      vc = vclda*(Hc - s2*Hcsos) +
     *     eclda*((4.0d0/3.0d0*s2 - v)*Hcsos  
     *     +(4.0d0/3.0d0*s**3-u)*s*Hcsossos)
c
      vx=vx+vc
C
      return
      end
c
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      subroutine DF_AM05s2(nup,ndn,dup,ddn,dtot,ec,
     &                     vcup,vcdn,decsup,decsdn)
c
c Armiento Mattsson am05 functional for exchange and correlation 
c "R. Armiento and A. E. Mattsson, PRB 72, 085108 (2005). " 
C 9/17/8:RJM spin dependence added 
C Note: subroutine outputs ec and some fraction of vc 
C
c input
c   nup,ndn = electron density, D [bohr**(-3)]
c output
c   ec       correlation energy per electron [hartree]
c   vc       correlation potential [hartree]
c   decs     correlation stress term. See GGASTRESS.
c
      implicit none
c
c     ** Input parameters
      real*8 n,kf,ex,vx,ntot,s,v,u
      real*8 nup,ndn,dup,ddn,sup,sdn,dtot
      real*8 vcldaup,vcldadn,dudotdd
      real*8 decsup,decsdn
c
c     ** Auxilliary parameters
      real*8 ec, vc, decs
c
c     ** Constants and shorthand notation
      real*8 pi, g, a, s2, c, two, third
      parameter (pi = 3.141592653589793238462643383279502884197d0)
      parameter (g = 0.8098d0, a = 2.804d0)
      parameter (c = 0.7168d0)
      parameter (two = 2.d0, third =1.d0/3.d0 )
c
c     ** Local variables
      real*8 eclda, vclda
	  real*8 rs,zet,a1,a2,a3,a4	
          real*8 vcup,vcdn,vcuplda,vcdnlda
          real*8 xup,xdn,xsosup,xsosdn,hcup,hcdn
	Real*8 hcsosup,hcsosdn 
c
c     ** Cutoff 
      if((n .le. 1.0d-16)) then
         ex = 0.0d0
         ec = 0.0d0
         vx = 0.0d0
         vc = 0.0d0
         return
      endif
c
c     *******************
c       LDA correlation
c     *******************
           ntot=nup+ndn
C check factors of 2. before nup
           kf=(3.*pi**2.*two*nup)**third
           sup=dup/(2.*kf*nup)
           kf=(3.*pi**2.*two*ndn)**third
           sdn=ddn/(two*kf*ndn)
           rs=(3./(4.*pi*ntot))**third
           zet=(nup-ndn)/(ntot)
           call DF_LSDC(rs,zet,eclda,vcldaup,vcldadn,a2,a3,a4)
C note eclda is the energy per particle eps
           dudotdd=(dtot**two-ddn**two-dup**two)/two
c     ********************
c         Correlation Energy
c     ********************
      s2=sup**2.
c     ** Interpolation index
      Xup = 1.0d0/(1.0d0 + a*s2)
c     ** Correlation refinement function 
      Hcup = Xup + g*(1.0d0 - Xup)
      s2=sdn**2.
c     ** Interpolation index
      Xdn = 1.0d0/(1.0d0 + a*s2)
c     ** Correlation refinement function 
      Hcdn = Xdn + g*(1.0d0 - Xdn)
c      
c     ** Correlation energy per particle, Ec = Integrate[n*ec]
      ec = eclda*(Hcup*nup+Hcdn*ndn)/ntot
C
c     *****************************
c        Correlation potential
c     *****************************
c     ** Interpolation index derivatives: 1/s dX/ds
      Xsosup = -2.0d0*a*Xup**2
      Xsosdn = -2.0d0*a*Xdn**2
c     ** Correlation refinement function derivatives, dF/ds 
      Hcsosup = Xsosup*(1.0d0 - g)
      Hcsosdn = Xsosdn*(1.0d0 - g)
c     ** vc formula for gradient dependent functional,
c     Generalized form of Eq. (24) in PRB 33, 8800 (1986) 
      vcup = 
     & (eclda-vcuplda)*ndn/ntot*(Hcup-Hcdn-sup**2.*Hcsosup)
     & +(eclda-vcdnlda)*nup/ntot*(sup**2.*Hcsosup)*dudotdd/dup**two
      vcdn = 
     & (eclda-vcdnlda)*nup/ntot*(Hcdn-Hcup-sdn**2.*Hcsosdn)
     & +(eclda-vcuplda)*ndn/ntot*(sdn**2.*Hcsosdn)*dudotdd/ddn**two
c     ** correlation stress term      
      decsup = eclda*(Hcsosup)
      decsdn = eclda*(Hcsosdn)
      return
      end
