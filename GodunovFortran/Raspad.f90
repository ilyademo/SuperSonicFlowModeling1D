subroutine RASPAD(AK,P1,R1,U1,P2,R2,U2,PBIG,UBIG,RBIG1,RBIG2,DL1,DL2,DP1,DP2) bind(c, name = 'raspad') 
    use iso_c_binding    
    IMPLICIT NONE
    REAL(c_double) AK,P1,R1,U1,P2,R2,U2

      REAL(c_double) PBIG,UBIG       ! parameters on contact surface
      REAL(c_double) RBIG1,RBIG2     ! density right and left of the surface
      REAL(c_double) DL1,DL2,DP1,DP2 ! speed of characteristics

      LOGICAL(c_bool) FLAG
      REAL(c_double) AK1,AK2,C1,C2,ST1,ST2,A1,A2,C2Z,C1Z,PI,EPSP,EPS
      REAL(c_double) PBIGN

      REAL,EXTERNAL :: FGOD,DIVFGOD

      EPSP=1.0E-7
      EPS=1.0E-7

      AK1=AK  
      AK2=AK  
      C1=SQRT(AK1*P1/R1)
      C2=SQRT(AK2*P2/R2)

      IF ( U1-U2+2.0*(C1/(AK1-1.0)+C2/(AK2-1.0)) < 0.0 ) THEN
        PBIG=1.0E-05
        RBIG1=0.0
        RBIG2=0.0
        RETURN
      ENDIF  

      FLAG=.FALSE.
      PBIG=(P1*R2*C2+P2*R1*C1+(U1-U2)*R1*C1*R2*C2)/(R1*C1+R2*C2)
      IF (PBIG<0.0) PBIG=1.0E-05
1     PBIGN=PBIG-(FGOD(PBIG,P1,R1,AK1)+FGOD(PBIG,P2,R2,AK2)-(U1-U2))/(DIVFGOD(PBIG,P1,R1,AK1)+DIVFGOD(PBIG,P2,R2,AK2))
            
      IF (PBIGN<0.0) THEN
        !WRITE (*,*) PBIG,FLAG
        !PAUSE
        PBIG=1.0E-05
        RBIG1=0.0
        RBIG2=0.0
        IF (FLAG) RETURN
        FLAG=.TRUE.
        GO TO 1
      ENDIF
      IF (ABS(PBIGN/PBIG-1.0) > EPS) THEN
        PBIG=PBIGN
        GO TO 1
      ENDIF

      ST1=0.5*(AK1-1.0)/AK1
      ST2=0.5*(AK2-1.0)/AK2
      IF ( PBIG >= P1-EPS ) THEN
        A1=SQRT(R1*(0.5*(AK1+1.0)*PBIG+0.5*(AK1-1.0)*P1))
      ELSE
        PI=PBIG/P1
        A1=0.5*(AK1-1.0)/AK1*R1*C1*(1.0-PI)/(1.0-PI**ST1)
      ENDIF
      IF ( PBIG >= P2-EPS ) THEN
        A2=SQRT(R2*(0.5*(AK2+1.0)*PBIG+0.5*(AK2-1.0)*P2))
      ELSE
        PI=PBIG/P2
        A2=0.5*(AK2-1.0)/AK2*R2*C2*(1.0-PI)/(1.0-PI**ST2)
      ENDIF

      UBIG=(A1*U1+A2*U2+P1-P2)/(A1+A2)
      
      IF ( P1-EPS<PBIG ) THEN
        RBIG1=R1*A1/(A1-R1*(U1-UBIG))
        DL1=U1-A1/R1
        DL2=DL1
      ELSE
        C1Z=C1+0.5*(AK1-1.0)*(U1-UBIG)
        RBIG1=AK1*PBIG/(C1Z*C1Z)
        DL1=U1-C1
        DL2=UBIG-C1Z
      ENDIF
      IF ( P2-EPS<PBIG ) THEN
        RBIG2=R2*A2/(A2+R2*(U2-UBIG))
        DP1=U2+A2/R2
        DP2=DP1
      ELSE
        C2Z=C2-0.5*(AK2-1.0)*(U2-UBIG)
        RBIG2=AK2*PBIG/(C2Z*C2Z)
        DP2=U2+C2
        DP1=UBIG+C2Z
      ENDIF

      RETURN
      END subroutine

      REAL function FGOD(PBIG,P,R,AK) bind(c, name = 'fgod')
      use iso_c_binding  
      REAL(c_double) PBIG,P,R,AK,PI,ST

      PI=PBIG/P
      ST=0.5*(AK-1.0)/AK
      IF (PI >= 1.0) THEN
        FGOD=SQRT(P/R)*(PI-1.0)/SQRT(0.5*(AK+1.0)*PI+0.5*(AK-1.0))
      ELSE
        FGOD=2.0/(AK-1.0)*SQRT(AK*P/R)*(PI**ST-1.0)
      ENDIF
      RETURN
      END

      REAL function DIVFGOD(PBIG,P,R,AK) bind(c, name = 'divfgod')
      use iso_c_binding  
      REAL(c_double) PBIG,P,R,AK,PI,ST
      REAL(c_double):: power = 3
      PI=PBIG/P
      ST=0.5*(AK-1.0)/AK
      IF (PI.GE.1.0) THEN
        DIVFGOD=((AK+1.0)*PI+(3.0*AK-1.0))/(4.0*AK*R*SQRT(AK*P/R)*SQRT((0.5*(AK+1.0)/AK*PI+0.5*(AK-1.0)/AK)**power))
      ELSE
         DIVFGOD=1.0/AK/PBIG*SQRT(AK*P/R)*PI**ST
      ENDIF
      RETURN
      END