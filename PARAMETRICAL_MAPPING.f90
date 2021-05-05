    SUBROUTINE PARAMETRICAL_MAPPING(PATCH,A,VALUES,PHI,QSI1,QSI2)
!
!	IN THIS SUBROUTINE, THE USER MUST INFORM THE NUMBER OF THE PARAMETRICAL TRIMMING CURVE (PATCH),THE PARAMETRIC COORDINATES QSI1 AND QSI2 AND 
!   THE COORDINATES OF THE CONTROL POINTS (VALUES). THE COORDINATES OF THE INTERPOLATED POINTS (X1,X2,X3) ARE THE OUTPUTS.  
!    
    USE NURBS_SURFACES
    USE TRIMMING_CURVES
!
    IMPLICIT NONE
!    
    INTEGER::PATCH,I,J,P,CP,Ni   
!   
    REAL*8::QSI,VALUES(N_BFCURVES(PATCH),2),A,QSI1,QSI2,PHI(N_BFCURVES(PATCH)),N(N_BFCURVES(PATCH)),Nplus1(N_BFCURVES(PATCH)+POLYNOMIAL_ORDER(PATCH)),&
            SOMA,Ai,Aiplus1,AiplusP,AiplusPplus1,DENOMINATOR1,DENOMINATOR2 
                      
!
    N=0.D0     
!
! UNIVARIATE BSPLINE BASIS FUNCTIONS 
!  
   !ZERO ORDER BSPLINES
    Ni=N_BFCURVES(PATCH)
    DO I=1,N_BFCURVES(PATCH)  
        Ai=KNOT_VECTOR(PATCH,Ni)
        Aiplus1=KNOT_VECTOR(PATCH,Ni+1)            
        IF(A.GE.Ai.AND.A.LT.Aiplus1)THEN
            N(I)=1.D0
        ELSE
            N(I)=0.D0      
        ENDIF                       
        Ni=Ni-1    
    ENDDO
   !ZERO ORDER BSPLINES(i+P)
    Ni=N_BFCURVES(PATCH)+POLYNOMIAL_ORDER(PATCH)
    DO I=1,N_BFCURVES(PATCH)+POLYNOMIAL_ORDER(PATCH)              
        Ai=KNOT_VECTOR(PATCH,Ni)
        Aiplus1=KNOT_VECTOR(PATCH,Ni+1)        
        IF(A.GE.Ai.AND.A.LT.Aiplus1)THEN
            Nplus1(I)=1.D0
        ELSE
            Nplus1(I)=0.D0       
        ENDIF                  
        Ni=Ni-1    
    ENDDO    
!    
!   COX-DE-BOOR RECURSION FORMULA  
!   
    DO P=1,POLYNOMIAL_ORDER(PATCH)
       !P ORDER BSPLINES 
        Ni=N_BFCURVES(PATCH)
        DO I=1,N_BFCURVES(PATCH)                               
            Ai=KNOT_VECTOR(PATCH,Ni)
            Aiplus1=KNOT_VECTOR(PATCH,Ni+1)
            AiplusP=KNOT_VECTOR(PATCH,Ni+P)
            AiplusPplus1=KNOT_VECTOR(PATCH,Ni+P+1)
            DENOMINATOR1=AiplusP-Ai
            DENOMINATOR2=AiplusPplus1-Aiplus1
            IF(DENOMINATOR1.NE.0.D0.AND.DENOMINATOR2.NE.0.D0)THEN
                N(I)=(A-Ai)*Nplus1(I+POLYNOMIAL_ORDER(PATCH)-P+1)/DENOMINATOR1+(AiplusPplus1-A)*Nplus1(I+POLYNOMIAL_ORDER(PATCH)-P)/DENOMINATOR2
            ENDIF    
            IF(DENOMINATOR1.NE.0.D0.AND.DENOMINATOR2.EQ.0.D0)THEN
                N(I)=(A-Ai)*Nplus1(I+POLYNOMIAL_ORDER(PATCH)-P+1)/DENOMINATOR1
            ENDIF 
            IF(DENOMINATOR1.EQ.0.D0.AND.DENOMINATOR2.NE.0.D0)THEN
                N(I)=(AiplusPplus1-A)*Nplus1(I+POLYNOMIAL_ORDER(PATCH)-P)/DENOMINATOR2
            ENDIF       
            IF(DENOMINATOR1.EQ.0.D0.AND.DENOMINATOR2.EQ.0.D0)THEN
                N(I)=0.D0
            ENDIF                                            
           Ni=Ni-1    
        ENDDO 
        IF(P.LT.POLYNOMIAL_ORDER(PATCH))THEN
           !P ORDER BSPLINES(i+P) 
            Ni=N_BFCURVES(PATCH)+POLYNOMIAL_ORDER(PATCH)-P  
            DO I=1,N_BFCURVES(PATCH)+POLYNOMIAL_ORDER(PATCH)-P                                       
                Ai=KNOT_VECTOR(PATCH,Ni)
                Aiplus1=KNOT_VECTOR(PATCH,Ni+1)
                AiplusP=KNOT_VECTOR(PATCH,Ni+P)
                AiplusPplus1=KNOT_VECTOR(PATCH,Ni+P+1)
                DENOMINATOR1=AiplusP-Ai
                DENOMINATOR2=AiplusPplus1-Aiplus1
                IF(DENOMINATOR1.NE.0.D0.AND.DENOMINATOR2.NE.0.D0)THEN
                    Nplus1(I)=(A-Ai)*Nplus1(I+1)/DENOMINATOR1+(AiplusPplus1-A)*Nplus1(I)/DENOMINATOR2
                ENDIF    
                IF(DENOMINATOR1.NE.0.D0.AND.DENOMINATOR2.EQ.0.D0)THEN
                    Nplus1(I)=(A-Ai)*Nplus1(I+1)/DENOMINATOR1
                ENDIF 
                IF(DENOMINATOR1.EQ.0.D0.AND.DENOMINATOR2.NE.0.D0)THEN
                    Nplus1(I)=(AiplusPplus1-A)*Nplus1(I)/DENOMINATOR2
                ENDIF       
                IF(DENOMINATOR1.EQ.0.D0.AND.DENOMINATOR2.EQ.0.D0)THEN
                    Nplus1(I)=0.D0
                ENDIF                                            
               Ni=Ni-1    
            ENDDO
        ENDIF             
    ENDDO
    IF(A.EQ.KNOT_VECTOR(PATCH,N_BFCURVES(PATCH)+POLYNOMIAL_ORDER(PATCH)+1))THEN
        N(1)=1.D0    
    ENDIF            
!
!   UNIVARIATE NURBS BASIS FUNCTIONS 
!
    PHI=0.D0
    SOMA=0.D0 
    DO I=1,N_BFCURVES(PATCH)
        CP=PARAMETRICCONTROLPOINTS_CONNECTIVITY(PATCH,I)
        SOMA=SOMA+N(I)*COORD_PARAMETRICCONTROLPOINTS(CP,3) 
    ENDDO
!
    DO I=1,N_BFCURVES(PATCH)
        CP=PARAMETRICCONTROLPOINTS_CONNECTIVITY(PATCH,I)    
        PHI(I)=N(I)*COORD_PARAMETRICCONTROLPOINTS(CP,3)/SOMA 
    ENDDO                  
! 
    QSI1=0.D0
    QSI2=0.D0	            
    DO I=1,N_BFCURVES(PATCH)
        QSI1=QSI1+PHI(I)*VALUES(I,1)
        QSI2=QSI2+PHI(I)*VALUES(I,2)
    ENDDO
!    
    END SUBROUTINE PARAMETRICAL_MAPPING       