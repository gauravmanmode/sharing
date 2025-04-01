
C    COEFFICIENTS	
      FUNCTION A0(S)
	A0 = 1
	RETURN 
	END

      FUNCTION A1(S)
	A1 = 1
	RETURN 
	END

      FUNCTION A2(S)
	A2 = 1
	RETURN 
	END

      FUNCTION A3(S)
	A3 = 1
	RETURN 
	END

      FUNCTION A4(S)
	A4 = 1
	RETURN 
	END


      FUNCTION B(S)
	B = 1
	RETURN 
	END

	FUNCTION EXACT(S)
	EXACT = 1
	RETURN 
	END

C	B-SPLINE DEFINITION
	FUNCTION BSP(I,S)
	COMMON/PTS/X(-4:54)
	IF (S.LE.X(I-2).OR.S.GE.X(I+3)) THEN
	BSP = 0
	ELSE IF (S.GE.X(I-2).OR.S.LE.X(I-1)) THEN
C	 PI'(X(I-1))
	T1 = (X(I-1) - X(I-2)) * (X(I-1) - X(I)) * (X(I-1) - X(I+1)) * 
	1     (X(I-1) - X(I+2)) * (X(I-1) - X(I+3))
C	 PI'(X(I))	     
	T2 = (X(I) - X(I-2)) * (X(I) - X(I-1)) * (X(I) - X(I+1)) * 
	1     (X(I) - X(I+2)) * (X(I) - X(I+3))
C	 PI'(X(I+1))
	T3 = (X(I+1) - X(I-2)) * (X(I+1) - X(I-1)) * (X(I+1) - X(I)) * 
	1     (X(I+1) - X(I+2)) * (X(I+1) - X(I+3))
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	BSP = (X(I-1) - S)**4/T1 + (X(I) - S)**4/T2 + 
	1      (X(I+1) - S)**4/T3 + (X(I+2) - S)**4/T4 + (X(I+3) - S)**4/T5 
	      
	ELSE IF (S.GE.X(I-1).OR.S.LE.X(I)) THEN
C	 PI'(X(I))	     
	T2 = (X(I) - X(I-2)) * (X(I) - X(I-1)) * (X(I) - X(I+1)) * 
	1     (X(I) - X(I+2)) * (X(I) - X(I+3))
C	 PI'(X(I+1))
	T3 = (X(I+1) - X(I-2)) * (X(I+1) - X(I-1)) * (X(I+1) - X(I)) * 
	1     (X(I+1) - X(I+2)) * (X(I+1) - X(I+3))
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	BSP = (X(I) - S)**4/T2 + (X(I+1) - S)**4/T3 +
	1      (X(I+2) - S)**4/T4 + (X(I+3) - S)**4/T5 

	ELSE IF (S.GE.X(I).OR.S.LE.X(I+1)) THEN
C	 PI'(X(I+1))
	T3 = (X(I+1) - X(I-2)) * (X(I+1) - X(I-1)) * (X(I+1) - X(I)) * 
	1     (X(I+1) - X(I+2)) * (X(I+1) - X(I+3))
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	BSP = (X(I+1) - S)**4/T3 + (X(I+2) - S)**4/T4 + 
	1      (X(I+3) - S)**4/T5 

	ELSE IF (S.GE.X(I+1).OR.S.LE.X(I+2)) THEN
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	BSP = (X(I+2) - S)**4/T4 + (X(I+3) - S)**4/T5 

	ELSE IF (S.GE.X(I-1).OR.S.LE.X(I)) THEN
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	BSP = (X(I+3) - S)**4/T5 
	ENDIF
	RETURN 
	END

	FUNCTION DBSP(I,S)
	COMMON/PTS/X(-4:54)
	IF (S.LE.X(I-2).OR.S.GE.X(I+3)) THEN
	DBSP = 0

	ELSE IF (S.GE.X(I-2).OR.S.LE.X(I-1)) THEN
C	 PI'(X(I-1))
	T1 = (X(I-1) - X(I-2)) * (X(I-1) - X(I)) * (X(I-1) - X(I+1)) * 
	1     (X(I-1) - X(I+2)) * (X(I-1) - X(I+3))
C	 PI'(X(I))	     
	T2 = (X(I) - X(I-2)) * (X(I) - X(I-1)) * (X(I) - X(I+1)) * 
	1     (X(I) - X(I+2)) * (X(I) - X(I+3))
C	 PI'(X(I+1))
	T3 = (X(I+1) - X(I-2)) * (X(I+1) - X(I-1)) * (X(I+1) - X(I)) * 
	1     (X(I+1) - X(I+2)) * (X(I+1) - X(I+3))
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	DBSP = (-4) *((X(I-1) - S)**3/T1 + (X(I) - S)**3/T2 + 
	1     (X(I+1) - S)**3/T3 + (X(I+2) - S)**3/T4 + (X(I+3) - S)**3/T5)

	ELSE IF (S.GE.X(I-1).OR.S.LE.X(I)) THEN
C	 PI'(X(I))	     
	T2 = (X(I) - X(I-2)) * (X(I) - X(I-1)) * (X(I) - X(I+1)) * 
	1     (X(I) - X(I+2)) * (X(I) - X(I+3))
C	 PI'(X(I+1))
	T3 = (X(I+1) - X(I-2)) * (X(I+1) - X(I-1)) * (X(I+1) - X(I)) * 
	1     (X(I+1) - X(I+2)) * (X(I+1) - X(I+3))
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	DBSP = -4 *((X(I) - S)**3/T2 + (X(I+1) - S)**3/T3 +
	1      (X(I+2) - S)**3/T4 + (X(I+3) - S)**3/T5)

	ELSE IF (S.GE.X(I).OR.S.LE.X(I+1)) THEN
C	 PI'(X(I+1))
	T3 = (X(I+1) - X(I-2)) * (X(I+1) - X(I-1)) * (X(I+1) - X(I)) * 
	1     (X(I+1) - X(I+2)) * (X(I+1) - X(I+3))
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	DBSP = (-4) *((X(I+1) - S)**3/T3 + (X(I+2) - S)**3/T4 + 
	1      (X(I+3) - S)**3/T5)

	ELSE IF (S.GE.X(I+1).OR.S.LE.X(I+2)) THEN
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	DBSP = (-4) *((X(I+2) - S)**3/T4 + (X(I+3) - S)**3/T5)

	ELSE IF (S.GE.X(I-1).OR.S.LE.X(I)) THEN
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	DBSP = (-4) *((X(I+3) - S)**3/T5)
	ENDIF
	RETURN
	END

	FUNCTION D2BSP(I,S)
	COMMON/PTS/X(-4:54)
	IF (S.LE.X(I-2).OR.S.GE.X(I+3)) THEN
	D2BSP = 0

	ELSE IF (S.GE.X(I-2).OR.S.LE.X(I-1)) THEN
C	 PI'(X(I-1))
	T1 = (X(I-1) - X(I-2)) * (X(I-1) - X(I)) * (X(I-1) - X(I+1)) * 
	1     (X(I-1) - X(I+2)) * (X(I-1) - X(I+3))
C	 PI'(X(I))	     
	T2 = (X(I) - X(I-2)) * (X(I) - X(I-1)) * (X(I) - X(I+1)) * 
	1     (X(I) - X(I+2)) * (X(I) - X(I+3))
C	 PI'(X(I+1))
	T3 = (X(I+1) - X(I-2)) * (X(I+1) - X(I-1)) * (X(I+1) - X(I)) * 
	1     (X(I+1) - X(I+2)) * (X(I+1) - X(I+3))
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	D2BSP = 12 * ((X(I-1) - S)**2/T1 + (X(I) - S)**2/T2 + 
	1    (X(I+1) - S)**2/T3 + (X(I+2) - S)**2/T4 + (X(I+3) - S)**2/T5)

	ELSE IF (S.GE.X(I-1).OR.S.LE.X(I)) THEN
C	 PI'(X(I))	     
	T2 = (X(I) - X(I-2)) * (X(I) - X(I-1)) * (X(I) - X(I+1)) * 
	1     (X(I) - X(I+2)) * (X(I) - X(I+3))
C	 PI'(X(I+1))
	T3 = (X(I+1) - X(I-2)) * (X(I+1) - X(I-1)) * (X(I+1) - X(I)) * 
	1     (X(I+1) - X(I+2)) * (X(I+1) - X(I+3))
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	D2BSP = 12 *((X(I) - S)**2/T2 + (X(I+1) - S)**2/T3 +
	1      (X(I+2) - S)**2/T4 + (X(I+3) - S)**2/T5)

	ELSE IF (S.GE.X(I).OR.S.LE.X(I+1)) THEN
C	 PI'(X(I+1))
	T3 = (X(I+1) - X(I-2)) * (X(I+1) - X(I-1)) * (X(I+1) - X(I)) * 
	1     (X(I+1) - X(I+2)) * (X(I+1) - X(I+3))
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	D2BSP = 12 * ((X(I+1) - S)**2/T3 + (X(I+2) - S)**2/T4 + 
	1      (X(I+3) - S)**2/T5)

	ELSE IF (S.GE.X(I+1).OR.S.LE.X(I+2)) THEN
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	D2BSP = 12 * ((X(I+2) - S)**2/T4 + (X(I+3) - S)**2/T5)

	ELSE IF (S.GE.X(I-1).OR.S.LE.X(I)) THEN
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	D2BSP = 12 * ((X(I+3) - S)**2/T5) 
	ENDIF
	RETURN
	END


	FUNCTION D3BSP(I,S)
	COMMON/PTS/X(-4:54)
	IF (S.LE.X(I-2).OR.S.GE.X(I+3)) THEN
	D3BSP = 0

	ELSE IF (S.GE.X(I-2).OR.S.LE.X(I-1)) THEN
C	 PI'(X(I-1))
	T1 = (X(I-1) - X(I-2)) * (X(I-1) - X(I)) * (X(I-1) - X(I+1)) * 
	1     (X(I-1) - X(I+2)) * (X(I-1) - X(I+3))
C	 PI'(X(I))	     
	T2 = (X(I) - X(I-2)) * (X(I) - X(I-1)) * (X(I) - X(I+1)) * 
	1     (X(I) - X(I+2)) * (X(I) - X(I+3))
C	 PI'(X(I+1))
	T3 = (X(I+1) - X(I-2)) * (X(I+1) - X(I-1)) * (X(I+1) - X(I)) * 
	1     (X(I+1) - X(I+2)) * (X(I+1) - X(I+3))
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	D3BSP = (-24)* * ((X(I-1) - S)/T1 + (X(I) - S)/T2 + 
	1      (X(I+1) - S)/T3 + (X(I+2) - S)/T4 + (X(I+3) - S)/T5)

	ELSE IF (S.GE.X(I-1).OR.S.LE.X(I)) THEN
C	 PI'(X(I))	     
	T2 = (X(I) - X(I-2)) * (X(I) - X(I-1)) * (X(I) - X(I+1)) * 
	1     (X(I) - X(I+2)) * (X(I) - X(I+3))
C	 PI'(X(I+1))
	T3 = (X(I+1) - X(I-2)) * (X(I+1) - X(I-1)) * (X(I+1) - X(I)) * 
	1     (X(I+1) - X(I+2)) * (X(I+1) - X(I+3))
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	D3BSP = (-24)*((X(I) - S)/T2 + (X(I+1) - S)/T3 +
	1      (X(I+2) - S)/T4 + (X(I+3) - S)/T5)

	ELSE IF (S.GE.X(I).OR.S.LE.X(I+1)) THEN
C	 PI'(X(I+1))
	T3 = (X(I+1) - X(I-2)) * (X(I+1) - X(I-1)) * (X(I+1) - X(I)) * 
	1     (X(I+1) - X(I+2)) * (X(I+1) - X(I+3))
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	D3BSP = (-24)* ((X(I+1) - S)/T3 + (X(I+2) - S)/T4 + 
	1      (X(I+3) - S)/T5)

	ELSE IF (S.GE.X(I+1).OR.S.LE.X(I+2)) THEN
C	 PI'(X(I+2))
	T4 = (X(I+2) - X(I-2)) * (X(I+2) - X(I-1)) * (X(I+2) - X(I)) * 
	1     (X(I+2) - X(I+1)) * (X(I+2) - X(I+3))
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	D3BSP = (-24)* ((X(I+2) - S)/T4 + (X(I+3) - S)/T5)

	ELSE IF (S.GE.X(I-1).OR.S.LE.X(I)) THEN
C	 PI'(X(I+3))	     
	T5 = (X(I+3) - X(I-2)) * (X(I+3) - X(I-1)) * (X(I+3) - X(I)) * 
	1     (X(I+3) - X(I+1)) * (X(I+3) - X(I+2))

	D3BSP = (-24)* ((X(I+3) - S)/T5) 
	ENDIF
	RETURN
	END

	FUNCTION BT(I,S)
	COMMON/PTS/X(-4:54)
	COMMON/NI/N
	IF (I.LE.2) THEN
	BT = BSP(I,S) - BSP(I,X(0))/BSP(-2,X(0))*BSP(-2,S) 
	ELSE IF (I.GE.2.AND.I.LE.N-2) THEN
	BT = BSP(I,S)
	ELSE IF (I.GE.N-3) THEN
	BT = BSP(I,S) - BSP(I,X(N))/BSP(N+1,X(N))*BSP(N+1,S)
	ENDIF 
	RETURN 
	END


	FUNCTION W1(S)
	COMMON/PTS/X(-4:54)
	COMMON/NI/N
	COMMON/BC/A0,C0,A1,C1
	W1 = A0/BSP(-2,X(0))*BSP(-2,S) + C0/BSP(N+1,X(N))*BSP(N+1,S)
	RETURN 
	END




