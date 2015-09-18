CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC ELIT.F from SciPy
CCC Inputs
CCC   HK = [Elliptic] Modulus, k = Eccentricity, e
CCC   PHI = Eccentric Anomaly, degrees
CCC Outputs
CCC   FE = F(k,phi), Incomplete Elliptic Integral of the First Kind
CCC   EE = F(k,phi), Incomplete Elliptic Integral of the Second Kind
CCC      = integral(sqrt(1 - (1 - HK**2) sin(t)**2)) for t = 0 to PHI
CCC      = arc length of ellipse from vertex i.e. tip of semi-MINOR axis
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        SUBROUTINE ELIT(HK,PHI,FE,EE)
C       ==================================================
C       Purpose: Compute complete and incomplete elliptic
C                integrals F(k,phi) and E(k,phi)
C       Input  : HK  --- Modulus k ( 0 ≤ k ≤ 1 )
C                Phi --- Argument ( in degrees )
C       Output : FE  --- F(k,phi)
C                EE  --- E(k,phi)
C       ==================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        G=0.0D0
        PI=3.14159265358979D0
        A0=1.0D0
        B0=DSQRT(1.0D0-HK*HK)
        D0=(PI/180.0D0)*PHI
        R=HK*HK
        IF (HK.EQ.1.0D0.AND.PHI.EQ.90.0D0) THEN
           FE=1.0D+300
           EE=1.0D0
        ELSE IF (HK.EQ.1.0D0) THEN
           FE=DLOG((1.0D0+DSIN(D0))/DCOS(D0))
           EE=DSIN(D0)
        ELSE
           FAC=1.0D0
           D=0.0D0
           DO 10 N=1,40
              A=(A0+B0)/2.0D0
              B=DSQRT(A0*B0)
              C=(A0-B0)/2.0D0
              FAC=2.0D0*FAC
              R=R+FAC*C*C
              IF (PHI.NE.90.0D0) THEN
                 D=D0+DATAN((B0/A0)*DTAN(D0))
                 G=G+C*DSIN(D)
                 D0=D+PI*INT(D/PI+.5D0)
              ENDIF
              A0=A
              B0=B
              IF (C.LT.1.0D-7) GO TO 15
10         CONTINUE
15         CK=PI/(2.0D0*A)
           CE=PI*(2.0D0-R)/(4.0D0*A)
           IF (PHI.EQ.90.0D0) THEN
              FE=CK
              EE=CE
           ELSE
              FE=D/(FAC*A)
              EE=FE*CE/CK+G
           ENDIF
        ENDIF
        RETURN
        END
