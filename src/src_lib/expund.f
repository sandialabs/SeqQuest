c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EXPUND
c
c
      SUBROUTINE EXPUND( argmax )
c---------------------------------------------------------------
c Purpose: return value x for which exp(-x) underflows
c
c Written: 8Jan91-PAS
c
c Revision history:
c   1May01-PAS/2.47: eliminate machine dependence
c  10Dec98-PAS/2.29: promote source to explicit double precision
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
CCC      DATA  rtpi / 1.772453850 9055160d0 /
CCCc
CCC F1MACH referred to the correct type R1MACH/D1MACH routine
CCC      xunderf= SQRT(- LOG(rtpi*F1MACH(1)))
CCC      xunderf = xunderf - 0.5* LOG(xunderf)/xunderf - 0.01
CCC      argmax = xunderf*xunderf
c
c The standard result gave 704.016; e(-704)~=10**(-306)
c A value of 60 --> 10**(-26.1) is more than adequate in tests
c
      argmax = 60.d0
c
      RETURN
      END
