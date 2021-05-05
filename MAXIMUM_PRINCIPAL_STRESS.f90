	SUBROUTINE MAXIMUM_PRINCIPAL_STRESS(KI,KII,KIII,THETAP)
!
	INTEGER::I,N,K
!
    REAL*8::KI,KII,KIII,THETAP
	REAL*8::F,X,Xl,Xu,TAL,X1,X2,Fl,Fu,F1,F2,PI
!
    PI=DACOS(-1.D0)
!
	Xl=-PI
	Xu=PI
	N=100
!
	X=Xl
	CALL PRINCIPAL_STRESS_FUNCTION(KI,KII,KIII,X,F)
	Fl=F
!
	X=Xu
	CALL PRINCIPAL_STRESS_FUNCTION(KI,KII,KIII,X,F)
	Fu=F
!
	TAL=0.381966D0
!
	X1=(1.D0-TAL)*Xl+TAL*Xu
	X=X1
	CALL PRINCIPAL_STRESS_FUNCTION(KI,KII,KIII,X,F)
	F1=F
!
	X2=TAL*Xl+(1.D0-TAL)*Xu	
	X=X2
	CALL PRINCIPAL_STRESS_FUNCTION(KI,KII,KIII,X,F)
	F2=F
!
	K=3
	DO I=K,N
		IF(F1.LT.F2) THEN
			Xl=X1
			Fl=F1
			X1=X2
			F1=F2
			X2=TAL*Xl+(1.D0-TAL)*Xu
			X=X2
			CALL PRINCIPAL_STRESS_FUNCTION(KI,KII,KIII,X,F)
			F2=F
		ELSE
			Xu=X2
			Fu=F2
			X2=X1
			F2=F1
			X1=(1.D0-TAL)*Xl+TAL*Xu
			X=X1
			CALL PRINCIPAL_STRESS_FUNCTION(KI,KII,KIII,X,F)
			F1=F
		ENDIF
	ENDDO
!
	THETAP=((X1+X2)/(2.D0))
!	
	END SUBROUTINE
!
	SUBROUTINE PRINCIPAL_STRESS_FUNCTION(KI,KII,KIII,X,F)
!
    REAL*8::KI,KII,KIII,X,F
!    
    F=KI*(3*COS(X/2.D0)+COS(3*X/2.D0))-KII*(3*SIN(X/2.D0)+3*SIN(3*X/2.D0))+DSQRT((KI*(3*COS(X/2.D0)+COS(3*X/2.D0))-KII*(3*SIN(X/2.D0)+3*SIN(3*X/2.D0)))**2+64*KIII**2*COS(X/2.D0)**2)
!    
	END SUBROUTINE    
	