C
C   PROGRAM TO SOLVE FOURTH ORDER BVP WITH CUBIC B-SPLINES BY GALERKIN METHOD
C
C	
	DIMENSION X(-3:53),A(0:50,7),RHS(0:50),ALP(0:50),WT(4),GP(4)
	COMMON/NI/N
	COMMON/PTS/X
	COMMON/BC/AA0,AA1,AC0,AC1
	COMMON/LIMITS/LL,IUL

	N=10
	X(0)=-1
	X(N)=1
	H=(X(N)-X(0))/N
	DO 10 I=-3,N+3
	X(I)=X(0)+I*H
10	PRINT*,X(I)

	AA0=0
	AC0=0
	AA1=(SINH(2.)-SIN(2.))/(4.*COSH(2.)+4.*COS(2.))
	AC1=-(SINH(2.)-SIN(2.))/(4.*COSH(2.)+4.*COS(2.))

	WT(1)=0.3478549
	WT(2)=0.6521452
	WT(3)=0.6521452
	WT(4)=0.3478549
	GP(1)=-0.8611363
	GP(2)=-0.339981
	GP(3)=0.339981
	GP(4)=0.8611363

	LL=0
	IUL=N

	DO 11 I=LL,IUL
	RHS(I)=0
	ALP(I)=0
	DO 12 J=1,7
12	A(I,J)=0
11	CONTINUE

	DO 30 I=LL,IUL
	DO 40 J=I-3,I+3
	IF(J.LT.LL .OR. J.GT.IUL) GOTO 40
	NCOL=J-I+4
	DO 50 II=I-2,I+1
	IF(II.LT.0 .OR. II.GT.N-1) GOTO 50
	DO 60 JJ=J-2,J+1
	IF(II.NE.JJ) GOTO 60
	DO 70 KK=1,4
	PT=((X(II+1)-X(II))/2.)*GP(KK)+(X(II+1)+X(II))/2.
	A(I,NCOL)=A(I,NCOL)+((-D3A0(PT)*BSPT(I,PT)*DBSPT(J,PT))
     1	          +(-3*D2A0(PT)*DBSPT(I,PT)*DBSPT(J,PT))
     2 	          +(-3*DA0(PT)*D2BSPT(I,PT)*DBSPT(J,PT))
     3	          +(-A0(PT)*D3BSPT(I,PT)*DBSPT(J,PT))
     4	          +D2A1(PT)*BSPT(I,PT)*DBSPT(J,PT)
	5              +2*DA1(PT)*DBSPT(I,PT)*DBSPT(J,PT)
	6              +A1(PT)*D2BSPT(I,PT)*DBSPT(J,PT)
	7              +A2(PT)*D2BSPT(J,PT)*BSPT(I,PT)
	8              +A3(PT)*DBSPT(J,PT)*BSPT(I,PT)
	9              +A4(PT)*BSPT(I,PT)*BSPT(J,PT))
	1              *WT(KK)*((X(II+1)-X(II))/2.)
70    CONTINUE
60    CONTINUE
50    CONTINUE
	A(I,NCOL)=A(I,NCOL)+(-DA0(X(N))*BSPT(I,X(N))*D2BSPT(J,X(N)))
	1           +(-A0(X(N))*DBSPT(I,X(N))*D2BSPT(J,X(N)))
	2           +(DA0(X(0))*BSPT(I,X(0))*D2BSPT(J,X(0)))
	3           +(A0(X(0))*DBSPT(I,X(0))*D2BSPT(J,X(0)))
40    CONTINUE
      DO 80 II=I-2,I+1
	IF(II.LT.0 .OR. II.GT.N-1) GOTO 80
	DO 90 KK=1,4
	PT=((X(II+1)-X(II))/2.)*GP(KK)+(X(II+1)+X(II))/2.	
	RHS(I)=RHS(I)+(B(PT)*BSPT(I,PT)+D3A0(PT)*BSPT(I,PT)*DW(PT)
	1    +3*D2A0(PT)*DBSPT(I,PT)*DW(PT)+3*DA0(PT)*D2BSPT(I,PT)*DW(PT)
	2    +A0(PT)*D3BSPT(I,PT)*DW(PT)-D2A1(PT)*BSPT(I,PT)*DW(PT)
	3    -2*DA1(PT)*DBSPT(I,PT)*DW(PT)-A1(PT)*D2BSPT(I,PT)*DW(PT)
	4    -A2(PT)*BSPT(I,PT)*D2W(PT)-A3(PT)*BSPT(I,PT)*DW(PT)
	5    -A4(PT)*BSPT(I,PT)*W(PT))*WT(KK)*((X(II+1)-X(II))/2.)
90    CONTINUE
80    CONTINUE
      RHS(I)=RHS(I)+(DA0(X(N))*BSPT(I,X(N))*D2W(X(N))
	1    +A0(X(N))*DBSPT(I,X(N))*D2W(X(N))
     2    -DA0(X(0))*BSPT(I,X(0))*D2W(X(0))
	3    -A0(X(0))*DBSPT(I,X(0))*D2W(X(0))
	4    -D2A0(X(N))*BSPT(I,X(N))*AC1-2*DA0(X(N))*DBSPT(I,X(N))*AC1
	5    -A0(X(N))*D2BSPT(I,X(N))*AC1+D2A0(X(0))*BSPT(I,X(0))*AA1
     6    +2*DA0(X(0))*DBSPT(I,X(0))*AA1+A0(X(0))*D2BSPT(I,X(0))*AA1
	7    +DA1(X(N))*BSPT(I,X(N))*AC1+A1(X(N))*DBSPT(I,X(N))*AC1
	8    -DA1(X(0))*BSPT(I,X(0))*AA1-A1(X(0))*DBSPT(I,X(0))*AA1)
30	CONTINUE

      CALL BAND(7,A,RHS,ALP)

	PRINT*,'POINT     SOLUTION     EXACT SOLUTION     ABSOLUTE ERROR'
	PRINT*,'--------------------------------------------------------'
	DO 31 I=1,N-1
	PT=X(I)
	SUM=W(PT)
	DO 32 J=LL,IUL
32	SUM=SUM+ALP(J)*BSPT(J,PT)
	Y=SUM
	PRINT*,PT,Y,EXACT(PT),ABS(Y-EXACT(PT))
31	CONTINUE
	STOP
	END
C---------FUNCTION A0(PT)	COFFICIENT OF D4Y
	FUNCTION A0(PT)
	A0=1
	RETURN
	END

C---------FUNCTION A1(PT)	COFFICIENT OF D3Y
	FUNCTION A1(PT)
	A1=0
	RETURN
	END

C---------FUNCTION A2(PT) COFFICIENT OF D2Y
	FUNCTION A2(PT)
	A2=0
	RETURN
	END

C---------FUNCTION A3(PT) COFFICIENT OF DY
	FUNCTION A3(PT)
	A3=0
	RETURN
	END
C---------FUNCTION A4(PT)	COFFICIENT OF Y
	FUNCTION A4(PT)
	A4=4
	RETURN
	END
C---------FUNCTION B(PT) RHS FUNCTION
	FUNCTION B(PT)
	B=1
	RETURN
	END
C     -----------------------------------------------------------------
      FUNCTION DA0(PT)
	DA0=0
	RETURN
	END
C     -----------------------------------------------------------------
      FUNCTION D2A0(PT)
	D2A0=0
	RETURN
	END
C     -----------------------------------------------------------------
      FUNCTION D3A0(PT)
	D3A0=0
	RETURN
	END
C     -----------------------------------------------------------------
      FUNCTION DA1(PT)
	DA1=0
	RETURN
	END
C     -----------------------------------------------------------------
      FUNCTION D2A1(PT)
	D2A1=0
	RETURN
	END
C---------FUNCTION EXACT(PT)
	FUNCTION EXACT(PT)
	EXACT=(.25)*(1-(2.*(SINH(1.)*SIN(1.)*SINH(PT)*SIN(PT)+COSH(1.)*
	1       COS(1.)*COSH(PT)*COS(PT))/(COS(2.)+COSH(2.))))
	RETURN
	END
C-------------------------------------------------------------
C	BSP=FUNCTION FOR ORIGINAL BASIS FUNCTION 
      FUNCTION BSP(I,PT)
      COMMON/PTS/X(-3:53)
      IF (PT.LE.X(I-2).OR.PT.GE.X(I+2)) THEN
            BSP=0.
      ELSE IF (PT.GT.X(I-2).AND.PT.LT.X(I-1)) THEN
	T1=(X(I-1)-X(I-2))*(X(I-1)-X(I))*(X(I-1)-X(I+1))*
     1     (X(I-1)-X(I+2))
	T2=(X(I)-X(I-2))*(X(I)-X(I-1))*(X(I)-X(I+1))*
     1     (X(I)-X(I+2))
	T3=(X(I+1)-X(I-2))*(X(I+1)-X(I-1))*(X(I+1)-X(I))*
     1     (X(I+1)-X(I+2))
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*
     1     (X(I+2)-X(I+1))
      BSP=(X(I-1)-PT)**3/T1+(X(I)-PT)**3/T2          
     1          +(X(I+1)-PT)**3/T3+(X(I+2)-PT)**3/T4         
      ELSE IF (PT.GE.X(I-1).AND.PT.LT.X(I)) THEN
      T2=(X(I)-X(I-2))*(X(I)-X(I-1))*(X(I)-X(I+1))*
     1     (X(I)-X(I+2))
	T3=(X(I+1)-X(I-2))*(X(I+1)-X(I-1))*(X(I+1)-X(I))*
     1     (X(I+1)-X(I+2))
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*
     1     (X(I+2)-X(I+1))
      BSP=(X(I)-PT)**3/T2+(X(I+1)-PT)**3/T3          
     1      +(X(I+2)-PT)**3/T4               
      ELSE IF (PT.GE.X(I).AND.PT.LT.X(I+1)) THEN
	T3=(X(I+1)-X(I-2))*(X(I+1)-X(I-1))*(X(I+1)-X(I))*
     1     (X(I+1)-X(I+2))
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*
     1     (X(I+2)-X(I+1))
      BSP=(X(I+1)-PT)**3/T3+(X(I+2)-PT)**3/T4
      ELSE IF (PT.GE.X(I+1).AND.PT.LT.X(I+2)) THEN
      T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*
     1     (X(I+2)-X(I+1))
      BSP=(X(I+2)-PT)**3/T4      
	ENDIF
	RETURN
	END
C----------------------------------------------------------------------
C	DBSP=DERIVATIVE OF BSP 
      FUNCTION DBSP(I,PT)
      COMMON/PTS/X(-3:53)
      IF (PT.LE.X(I-2).OR.PT.GE.X(I+2)) THEN
            DBSP=0.
      ELSE IF (PT.GT.X(I-2).AND.PT.LT.X(I-1)) THEN
	T1=(X(I-1)-X(I-2))*(X(I-1)-X(I))*(X(I-1)-X(I+1))*
     1     (X(I-1)-X(I+2))
	T2=(X(I)-X(I-2))*(X(I)-X(I-1))*(X(I)-X(I+1))*
     1     (X(I)-X(I+2))
	T3=(X(I+1)-X(I-2))*(X(I+1)-X(I-1))*(X(I+1)-X(I))*
     1     (X(I+1)-X(I+2))
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*
     1     (X(I+2)-X(I+1))
      DBSP=-3.*((X(I-1)-PT)**2/T1+(X(I)-PT)**2/T2          
     1          +(X(I+1)-PT)**2/T3+(X(I+2)-PT)**2/T4)         
      ELSE IF (PT.GE.X(I-1).AND.PT.LT.X(I)) THEN
      T2=(X(I)-X(I-2))*(X(I)-X(I-1))*(X(I)-X(I+1))*
     1     (X(I)-X(I+2))
	T3=(X(I+1)-X(I-2))*(X(I+1)-X(I-1))*(X(I+1)-X(I))*
     1     (X(I+1)-X(I+2))
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*
     1     (X(I+2)-X(I+1))
      DBSP=-3.*((X(I)-PT)**2/T2+(X(I+1)-PT)**2/T3          
     1      +(X(I+2)-PT)**2/T4)               
      ELSE IF (PT.GE.X(I).AND.PT.LT.X(I+1)) THEN
	T3=(X(I+1)-X(I-2))*(X(I+1)-X(I-1))*(X(I+1)-X(I))*
     1     (X(I+1)-X(I+2))
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*
     1     (X(I+2)-X(I+1))
      DBSP=-3.*((X(I+1)-PT)**2/T3+(X(I+2)-PT)**2/T4)
      ELSE IF (PT.GE.X(I+1).AND.PT.LT.X(I+2)) THEN
      T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*
     1     (X(I+2)-X(I+1))
      DBSP=-3.*(X(I+2)-PT)**2/T4      
	ENDIF
	RETURN
	END
C----------------------------------------------------------------------
C	D2BSP=SECOND DERIVATIVE OF BSP 
      FUNCTION D2BSP(I,PT)
      COMMON/PTS/X(-3:53)
      IF (PT.LE.X(I-2).OR.PT.GE.X(I+2)) THEN
            D2BSP=0.
      ELSE IF (PT.GT.X(I-2).AND.PT.LT.X(I-1)) THEN
	T1=(X(I-1)-X(I-2))*(X(I-1)-X(I))*(X(I-1)-X(I+1))*
     1     (X(I-1)-X(I+2))
	T2=(X(I)-X(I-2))*(X(I)-X(I-1))*(X(I)-X(I+1))*
     1     (X(I)-X(I+2))
	T3=(X(I+1)-X(I-2))*(X(I+1)-X(I-1))*(X(I+1)-X(I))*
     1     (X(I+1)-X(I+2))
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*
     1     (X(I+2)-X(I+1))
      D2BSP=6.*((X(I-1)-PT)/T1+(X(I)-PT)/T2          
     1          +(X(I+1)-PT)/T3+(X(I+2)-PT)/T4)         
      ELSE IF (PT.GE.X(I-1).AND.PT.LT.X(I)) THEN
      T2=(X(I)-X(I-2))*(X(I)-X(I-1))*(X(I)-X(I+1))*
     1     (X(I)-X(I+2))
	T3=(X(I+1)-X(I-2))*(X(I+1)-X(I-1))*(X(I+1)-X(I))*
     1     (X(I+1)-X(I+2))
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*
     1     (X(I+2)-X(I+1))
      D2BSP=6.*((X(I)-PT)/T2+(X(I+1)-PT)/T3          
     1      +(X(I+2)-PT)/T4)               
      ELSE IF (PT.GE.X(I).AND.PT.LT.X(I+1)) THEN
	T3=(X(I+1)-X(I-2))*(X(I+1)-X(I-1))*(X(I+1)-X(I))*
     1     (X(I+1)-X(I+2))
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*
     1     (X(I+2)-X(I+1))
      D2BSP=6.*((X(I+1)-PT)/T3+(X(I+2)-PT)/T4)
      ELSE IF (PT.GE.X(I+1).AND.PT.LT.X(I+2)) THEN
      T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*
     1     (X(I+2)-X(I+1))
      D2BSP=6.*(X(I+2)-PT)/T4      
	ENDIF
	RETURN
	END
C     -----------------------------------------------------------------
C     D3BSP REPRESENTS DOUBLE DERIVATIVE OF ORIGINAL BASIS FUNCTIONS
      FUNCTION D3BSP(I,PT)
	COMMON/PTS/X(-3:53)
	IF(PT.LE.X(I-2).OR.PT.GE.X(I+2))THEN
	          D2BSP=0
	ELSE IF((X(I-2).LE.PT).AND.(PT.LE.X(I-1)))THEN
	T1=(X(I-1)-X(I-2))*(X(I-1)-X(I))*(X(I-1)-X(I+1))*(X(I-1)-X(I+2))
	T2=(X(I)-X(I-2))*(X(I)-X(I-1))*(X(I)-X(I+1))*(X(I)-X(I+2))
	T3=(X(I+1)-X(I-2))*(X(I+1)-X(I-1))*(X(I+1)-X(I))*(X(I+1)-X(I+2))
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*(X(I+2)-X(I+1))
	D3BSP=(-6/T1)+(-6/T2)+(-6/T3)+(-6/T4)
      ELSE IF((X(I-1).LE.PT).AND.(PT.LE.X(I)))THEN
	T2=(X(I)-X(I-2))*(X(I)-X(I-1))*(X(I)-X(I+1))*(X(I)-X(I+2))
	T3=(X(I+1)-X(I-2))*(X(I+1)-X(I-1))*(X(I+1)-X(I))*(X(I+1)-X(I+2))
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*(X(I+2)-X(I+1))
	D3BSP=(-6/T2)+(-6/T3)+(-6/T4)
      ELSE IF((X(I).LE.PT).AND.(PT.LE.X(I+1)))THEN
	T3=(X(I+1)-X(I-2))*(X(I+1)-X(I-1))*(X(I+1)-X(I))*(X(I+1)-X(I+2))
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*(X(I+2)-X(I+1))
	D3BSP=(-6/T3)+(-6/T4)
	ELSE IF((X(I+1).LE.PT).AND.(PT.LE.X(I+2)))THEN
	T4=(X(I+2)-X(I-2))*(X(I+2)-X(I-1))*(X(I+2)-X(I))*(X(I+2)-X(I+1))
	D3BSP=(-6/T4)
	ENDIF
	RETURN
	END
C     -----------------------------------------------------------------
C     BSPT REPRESENTS RE-ARRANGED BASIS FUNCTIONS
      FUNCTION BSPT(I,PT)
	COMMON/NI/N
	COMMON/PTS/X(-3:53)
	IF(I.EQ.0) BSPT=BSP(0,PT)-BSP(0,X(0))*BSP(-1,PT)/BSP(-1,X(0))
      IF(I.EQ.1) BSPT=BSP(1,PT)-BSP(1,X(0))*BSP(-1,PT)/BSP(-1,X(0))
	IF(2.LE.I.AND.I.LE.N-2) BSPT=BSP(I,PT)
	IF(I.EQ.N-1) THEN
	  NM1=N-1
	  BSPT=BSP(NM1,PT)-BSP(NM1,X(N))*BSP(N+1,PT)/BSP(N+1,X(N))
	ENDIF
	IF(I.EQ.N)THEN
	   BSPT=BSP(N,PT)-BSP(N,X(N))*BSP(N+1,PT)/BSP(N+1,X(N))
	END IF
	RETURN
	END
C     -----------------------------------------------------------------
C     DBSPT REPRESENTS DERIVATIVE OF RE-ARRANGED BASIS FUNCTION
      FUNCTION DBSPT(I,PT)
      COMMON/NI/N
	COMMON/PTS/X(-3:53)
      IF(I.EQ.0) DBSPT=DBSP(0,PT)-BSP(0,X(0))*DBSP(-1,PT)/BSP(-1,X(0))
      IF(I.EQ.1) DBSPT=DBSP(1,PT)-BSP(1,X(0))*DBSP(-1,PT)/BSP(-1,X(0))
	IF(2.LE.I.AND.I.LE.N-2) DBSPT=DBSP(I,PT)
	IF(I.EQ.N-1) THEN
	  NM1=N-1
	  DBSPT=DBSP(NM1,PT)-BSP(NM1,X(N))*DBSP(N+1,PT)/BSP(N+1,X(N))
	END IF
	IF(I.EQ.N) THEN
	  DBSPT=DBSP(N,PT)-BSP(N,X(N))*DBSP(N+1,PT)/BSP(N+1,X(N))
     	END IF
	RETURN
	END
C     -----------------------------------------------------------------
C     D2BSPT REPRESENT DOUBLE DERIVATIVE OF RE-ARRANGED BASIS FUNCTIONS
      FUNCTION D2BSPT(I,PT)
	COMMON/NI/N
	COMMON/PTS/X(-3:53)
	IF(I.EQ.0)D2BSPT=D2BSP(0,PT)-BSP(0,X(0))*
     1                             D2BSP(-1,PT)/BSP(-1,X(0))
	IF(I.EQ.1)D2BSPT=D2BSP(1,PT)-BSP(1,X(0))*
     1	                             D2BSP(-1,PT)/BSP(-1,X(0))
	IF(2.LE.I.AND.I.LE.N-2)D2BSPT=D2BSP(I,PT)
	IF(I.EQ.N-1)THEN
	   NM1=N-1
	   D2BSPT=D2BSP(NM1,PT)-BSP(NM1,X(N))*D2BSP(N+1,PT)/BSP(N+1,X(N))
	ENDIF
	IF(I.EQ.N)THEN
	   D2BSPT=D2BSP(N,PT)-BSP(N,X(N))*D2BSP(N+1,PT)/BSP(N+1,X(N))
	ENDIF
	RETURN
	END
C     -----------------------------------------------------------------
C     D3BSPT REPRESENTS THIRD DERIVATIVE OF RE-ARRANGED BASIS FUNCTIONS
	FUNCTION D3BSPT(I,PT)
	COMMON/NI/N
	COMMON/PTS/X(-3:53)
	IF(I.EQ.0)D3BSPT=D3BSP(0,PT)-BSP(0,X(0))*
	1                             D3BSP(-1,PT)/BSP(-1,X(0))
	IF(I.EQ.1)D3BSPT=D3BSP(1,PT)-BSP(1,X(0))*
     1	                             D3BSP(-1,PT)/BSP(-1,X(0))
	IF(2.LE.I.AND.I.LE.N-2)D3BSPT=D3BSP(I,PT)
	IF(I.EQ.N-1)THEN
	   NM1=N-1
	   D3BSPT=D3BSP(NM1,PT)-BSP(NM1,X(N))*D3BSP(N+1,PT)/BSP(N+1,X(N))
	ENDIF
	IF(I.EQ.N)THEN
	   D3BSPT=D3BSP(N,PT)-BSP(N,X(N))*D3BSP(N+1,PT)/BSP(N+1,X(N))
	ENDIF
	RETURN
	END
C    	-----------------------------------------------------------------
C     W(PT)=THE FUNCTION WHICH TAKES CARE OF GIVEN BOUNDARY CONDITIONS
      FUNCTION W(PT)
	COMMON/NI/N
	COMMON/PTS/X(-3:53)
	COMMON/BC/AA0,AA1,AC0,AC1
	W=AA0*BSP(-1,PT)/BSP(-1,X(0))+AC0*BSP(N+1,PT)/BSP(N+1,X(N))
	RETURN
	END
C     -----------------------------------------------------------------
C     DW=DERIVATIVE OF W(PT)
      FUNCTION DW(PT)
	COMMON/PTS/X(-3:53)
	COMMON/NI/N
	COMMON/BC/AA0,AA1,AC0,AC1
      DW=AA0*DBSP(-1,PT)/BSP(-1,X(0))+AC0*DBSP(N+1,PT)/BSP(N+1,X(N))
      RETURN
	END
C     -----------------------------------------------------------------
C     D2W=DOUBLE DERIVATIVE OF W(PT)
      FUNCTION D2W(PT)
	COMMON/PTS/X(-3:53)
	COMMON/NI/N
	COMMON/BC/AA0,AA1,AC0,AC1
      D2W=AA0*D2BSP(-1,PT)/BSP(-1,X(0))+AC0*D2BSP(N+1,PT)/BSP(N+1,X(N))
      RETURN
	END
