c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c Elements moduls
c Purpose: contains atomic (element data)
c Written: Peter A. Schultz, Febraury-2013, for 2.65
c
c No private data or private routines.
c Public routines:
c  ELEMSYMB2Z - returns Z for an input atomic symbol
c  ELEMZ2SYMB - returns atomic symbol for an input Z
c  ELEM_ULG12 - vdW parameters for given input Z, for Kim&WAG ULG corerction
c  ELEM_DFTD2 - vdW parameters for given input Z, for Grimme D2 correction
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ELEMSYMB2Z
c
c
      subroutine ELEMSYMB2Z( iatZ, atsymb )
c---------------------------------------------------------------
c Purpose: convert atom Z into atom symbol, and vice-versa
c Written: Peter A. Schultz, for 2.65 (vdW corrections)
c---------------------------------------------------------------
c
c  Returns integer "Z" for valid atomic symbol
c  Returns a "0' for a blank '  ' symbol input
c  Returns a '-1' for an unknown synbol
c
c Use 120 elements, to provide latitutde for artificial atom types
c  (e.g. floating orbitals, terminators)
c
      IMPLICIT NONE
c
      INTEGER  iatZ
      CHARACTER*2  atsymb
c
      INTEGER      isymb
      CHARACTER*2  symtable(120), symb
      DATA         symtable
     $ / 'H ','He','Li','Be','B ',  'C ','N ','O ','F ','Ne',
     $    'Na','Mg','Al','Si','P ',  'S ','Cl','Ar','K ','Ca',
     $     'Sc','Ti','V ','Cr','Mn',  'Fe','Co','Ni','Cu','Zn',
     $      'Ga','Ge','As','Se','Br',  'Kr','Rb','Sr','Y ','Zr',
     $       'Nb','Mo','Tc','Ru','Rh',  'Pd','Ag','Cd','In','Sn',
     $   'Sb','Te','I ','Xe','Cs',  'Ba','La','Ce','Pr','Nd',
     $    'Pm','Sm','Eu','Gd','Tb',  'Dy','Ho','Er','Tm','Yb',
     $     'Lu','Hf','Ta','W ','Re',  'Os','Ir','Pt','Au','Hg',
     $      'Tl','Pb','Bi','Po','At',  'Rn','Fr','Ra','Ac','Th',
     $       'Pa','U ','Np','Pu','Am',  'Cm','Bk','Cf','Es','Fm',
     $   'Md','No','Lr','A4','A5',  'A6','A7','A8','A9','T0',
     $    'T1','T2','T3','T4','T5',  'T6','T7','T8','T9','FO' /
c
      iatZ = 0
      if( atsymb .eq. '  ' ) RETURN
c
      symb = atsymb
c     Left-justify symbol, to ensure valid comparison to table:
      if( symb(1:1) .eq. ' ' )then
        symb(1:1) = symb(2:2)
        symb(2:2) = ' '
      endif
c
      do  isymb=1,120
c       Find the atom in the table:
        if( symb.eq.symtable(isymb) )then
          iatZ = isymb
          RETURN
        endif
      enddo
      iatZ = -1
c
      RETURN
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ELEMZ2SYMB
c
      ENTRY ELEMZ2SYMB( iatZ, atsymb )
c
      if( iatZ.gt.0 .and. iatZ.le.120 )then
        atsymb = symtable(iatZ)
      else
c       Not a defined Z
        atsymb = '??'
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ELEM_ULG12
c
c
      subroutine ELEM_ULG12( iatZ, c6, r0 )
c---------------------------------------------------------------
c Purpose: c6/r0 parameters for ULG vdW correction
c Written: Peter A. Schultz, for 2.65 (vdW corrections)
c---------------------------------------------------------------
c
c Following: H.Kim, J.-M.Choi, and WAG, JPhChL 3, 360 (2012).
c   vdW Data: D0 is in units of (kcal/mol)*A^6, R0 is in units of Angstrom
c   Output c6 is in units of Ry, r0 in units of bohr
c  
      IMPLICIT NONE
c
      REAL*8   zero,one,two
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0 )
c
c Input:
      INTEGER  iatZ
c Output:
      REAL*8   c6,r0
c
c Local:
      REAL*8   d0
c Data:
      INTEGER  nzulg
      PARAMETER  ( nzulg = 103 )
c        'H ',  'He',  'Li',  'Be',    'B ',  'C ',  'N ',  'O ',
c         'F ',  'Ne',  'Na',  'Mg',    'Al',  'Si',  'P ',  'S ',
c          'Cl',  'Ar',  'K ',  'Ca',    'Sc',  'Ti',  'V ',  'Cr',
c           'Mn',  'Fe',  'Co',  'Ni',    'Cu',  'Zn',  'Ga',  'Ge',
c        'As',  'Se',  'Br',  'Kr',    'Rb',  'Sr',  'Y ',  'Zr',
c         'Nb',  'Mo',  'Tc',  'Ru',    'Rh',  'Pd',  'Ag',  'Cd',
c          'In',  'Sn',  'Sb',  'Te',    'I ',  'Xe',  'Cs',  'Ba',
c           'La',  'Ce',  'Pr',  'Nd',    'Pm',  'Sm',  'Eu',  'Gd',
c        'Tb',  'Dy',  'Ho',  'Er',    'Tm',  'Yb',  'Lu',  'Hf',
c         'Ta',  'W ',  'Re',  'Os',    'Ir',  'Pt',  'Au',  'Hg',
c          'Tl',  'Pb',  'Bi',  'Po',    'At',  'Rn',  'Fr',  'Ra',
c           'Ac',  'Th',  'Pa',  'U ',    'Np',  'Pu',  'Am',  'Cm',
c        'Bk',  'Cf',  'Es',  'Fm',    'Md',  'No',  'Lr'
      REAL*8  xkcal2Ry, xA2bohr
      PARAMETER  ( xkcal2Ry=0.0031864d0, xA2bohr=1.8897265d0 )
      REAL*8  D0ulg, R0ulg
      DIMENSION  D0ulg(nzulg), R0ulg(nzulg)
      DATA    D0ulg
     $ / 0.0440,0.0560,0.0250,0.0850,  0.1800,0.1050,0.0690,0.0600,
     $    0.0500,0.0420,0.0300,0.1110,  0.5050,0.4020,0.3050,0.2740,
     $     0.2270,0.1850,0.0350,0.2380,  0.0190,0.0170,0.0160,0.0150,
     $      0.0130,0.0130,0.0140,0.0150,  0.0050,0.1240,0.4150,0.3790,
     $   0.3090,0.2910,0.2510,0.2200,  0.0400,0.2350,0.0720,0.0690,
     $    0.0590,0.0560,0.0480,0.0560,  0.0530,0.0480,0.0360,0.2280,
     $     0.5990,0.5670,0.4490,0.3980,  0.3390,0.3320,0.0450,0.3640,
     $      0.0170,0.0130,0.0100,0.0100,  0.0090,0.0080,0.0080,0.0090,
     $   0.0070,0.0070,0.0070,0.0070,  0.0060,0.2280,0.0410,0.0720,
     $    0.0810,0.0670,0.0660,0.0370,  0.0730,0.0800,0.0390,0.3850,
     $     0.6800,0.6630,0.5180,0.3250,  0.2840,0.2480,0.0500,0.4040,
     $      0.0330,0.0260,0.0220,0.0220,  0.0190,0.0160,0.0140,0.0130,
     $   0.0130,0.0130,0.0120,0.0120,  0.0110,0.0110,0.0110 / 
      DATA     R0ulg
     $ / 2.8860,2.3620,2.4510,2.7450,  4.0830,3.8510,3.6600,3.5000,
     $    3.3640,3.2430,2.9830,3.0210,  4.4990,4.2950,4.1470,4.0350,
     $     3.9470,3.8680,3.8120,3.3990,  3.2950,3.1750,3.1440,3.0230,
     $      2.9610,2.9120,2.8720,2.8340,  3.4950,2.7630,4.3830,4.2800,
     $   4.2300,4.2050,4.1890,4.1410,  4.1140,3.6410,3.3450,3.1240,
     $    3.1650,3.0520,2.9980,2.9630,  2.9290,2.8990,3.1480,2.8480,
     $     4.4630,4.3920,4.4200,4.4700,  4.5000,4.4040,4.5170,3.7030,
     $      3.5220,3.5560,3.6060,3.5750,  3.5470,3.5200,3.4930,3.3680,
     $   3.4510,3.4280,3.4090,3.3910,  3.3740,3.3550,3.6400,3.1410,
     $    3.1700,3.0690,2.9540,3.1200,  2.8400,2.7540,3.2930,2.7050,
     $     4.3470,4.2970,4.3700,4.7090,  4.7500,4.7650,4.9000,3.6770,
     $      3.4780,3.3960,3.4240,3.3950,  3.4240,3.4240,3.3810,3.3260,
     $   3.3390,3.3130,3.2990,3.2860,  3.2740,3.2480,3.2360 /
c
c >>>> EXECUTABLE CODE
c
      if( iatZ.gt.0 .and. iatZ.le.nzulg )then
        d0 = D0ulg(iatZ)
        r0 = R0ulg(iatZ)
        c6 = two*d0*(r0**6)
      else
        c6 = -1.
        r0 = one
      endif
c Convert units to Ry and bohr:
      c6 = c6 * xkcal2Ry * (xA2bohr**6)
      r0 = r0 * xA2bohr
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ELEM_DFTD2
c
c
      subroutine ELEM_DFTD2( iatZ, c6, r0 )
c---------------------------------------------------------------
c Purpose: c6/r0 parameters for Grimme06 DFT-D2 vdW correction
c Written: Peter A. Schultz, for 2.65 (vdW corrections)
c---------------------------------------------------------------
c
c Following: S.Grimme, J.Comput.Chem. 27, 1787 (2006).
c
c  vdW Data: C6 in units of J(nm^6)*mol^-1, R0 in units of Angstrom
c  Output c6 is in units of Ry*bohr^6, r0 in units of bohr
c  Conversion arithmetic:
c  1 J = 2.2937d+17 Ha = 4.5874d+17 Ry
c  1 nm = 10 A, 1 A = 1.8897265 bohr
c  Avogadro's = 6.022137d+23 mol^-1
c  1 Jnm6/mol = (4.5874d17 Ry/J) * (10A/nm)**6 * (1/6.022137d+23) = 0.76176 Ry*A^6
c
      IMPLICIT NONE
c
      REAL*8   zero,one,two
      PARAMETER  ( zero=0.d0, one=1.d0, two=2.d0 )
c
c Input:
      INTEGER  iatZ
c Output:
      REAL*8   c6,r0
c
c Local:
      REAL*8   xJtoRy, xA2bohr
      PARAMETER  ( xJtoRy=0.76176d0, xA2bohr=1.8897265d0 )
      INTEGER  nzd2
      PARAMETER  ( nzd2 = 54 )
c        'H ',  'He',  'Li',  'Be',    'B ',  'C ',  'N ',  'O ',
c         'F ',  'Ne',  'Na',  'Mg',    'Al',  'Si',  'P ',  'S ',
c          'Cl',  'Ar',  'K ',  'Ca',    'Sc',  'Ti',  'V ',  'Cr',
c           'Mn',  'Fe',  'Co',  'Ni',    'Cu',  'Zn',  'Ga',  'Ge',
c        'As',  'Se',  'Br',  'Kr',    'Rb',  'Sr',  'Y ',  'Zr',
c         'Nb',  'Mo',  'Tc',  'Ru',    'Rh',  'Pd',  'Ag',  'Cd',
c          'In',  'Sn',  'Sb',  'Te',    'I ',  'Xe'
      REAL*8     C6D2, R0D2
      DIMENSION  C6D2(nzd2), R0D2(nzd2)
      DATA    C6D2
     $ /  0.14,  0.08,  1.61,  1.61,    3.13,  1.75,  1.23,  0.70, 
     $     0.75,  0.63,  5.71,  5.71,   10.79,  9.23,  7.84,  5.57, 
     $      5.07,  4.61, 10.80, 10.80,   10.80, 10.80, 10.80, 10.80,
     $      10.80, 10.80, 10.80, 10.80,   10.80, 10.80, 16.99, 17.10, 
     $   16.37, 12.64, 12.47, 12.01,   24.67, 24.67, 24.67, 24.67,
     $    24.67, 24.67, 24.67, 24.67,   24.67, 24.67, 24.67, 24.67,
     $     37.32, 38.71, 38.44, 31.74,   31.50, 29.99 /
      DATA     R0D2
     $ / 1.001, 1.012, 0.825, 1.408,   1.485, 1.452, 1.397, 1.342,
     $    1.287, 1.243, 1.144, 1.364,   1.639, 1.716, 1.705, 1.683,
     $     1.639, 1.595, 1.485, 1.474,   1.562, 1.562, 1.562, 1.562,
     $      1.562, 1.562, 1.562, 1.562,   1.562, 1.562, 1.650, 1.727,
     $   1.760, 1.771, 1.749, 1.727,   1.628, 1.606, 1.639, 1.639,
     $    1.639, 1.639, 1.639, 1.639,   1.639, 1.639, 1.639, 1.639,
     $     1.672, 1.804, 1.881, 1.892,   1.892, 1.881 /
c
c >>>> EXECUTABLE CODE
c
      if( iatZ.gt.0 .and. iatZ.le.nzd2 )then
        c6 = C6D2(iatZ)
        r0 = R0D2(iatZ)
      else
c       If not in set, set to invalid values
        c6 = -1.
        r0 = one
      endif
c
c Convert units to Ry and bohr:
      c6 = c6 * xJtoRy * (xA2bohr**6)
      r0 = r0 * xA2bohr
c
      RETURN
      END
