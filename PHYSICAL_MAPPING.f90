    SUBROUTINE PHYSICAL_MAPPING(PATCH,QSI1,QSI2,VALUES,PHI,X1,X2,X3)
!
!	IN THIS SUBROUTINE, THE USER MUST INFORM THE NUMBER OF THE PARAMETRICAL SURFACE (PATCH),THE PARAMETRIC COORDINATES QSI1 AND QSI2 AND 
!   THE COORDINATES OF THE CONTROL POINTS (VALUES). THE COORDINATES OF THE INTERPOLATED POINTS (X1,X2,X3) ARE THE OUTPUTS.  
!    
    USE NURBS_SURFACES
!
    IMPLICIT NONE
!    
    INTEGER::PATCH,I,J,P,GLOB_NUM,CP   
!   
    REAL*8::QSI1,QSI2,VALUES(N_GBF(PATCH),3),X1,X2,X3,PHI(N_GBF(PATCH)),N(N_BF(PATCH,1)),M(N_BF(PATCH,2)),&
            Nplus1(N_BF(PATCH,1)+POLYNOMIAL_ORDERS(PATCH,1)),Mplus1(N_BF(PATCH,2)+POLYNOMIAL_ORDERS(PATCH,2)),SOMA,&
            Ni,QSIi,QSIiplus1,QSIiplusP,QSIiplusPplus1,DENOMINATOR1,DENOMINATOR2           
!
    X1=0.D0
    X2=0.D0
    X3=0.D0	
    N=0.D0
    M=0.D0   
    SOMA=0.D0    
!
! N UNIVARIATE BSPLINE BASIS FUNCTIONS 
!  
   !ZERO ORDER BSPLINES
    Ni=N_BF(PATCH,1)
    DO I=1,N_BF(PATCH,1)  
        QSIi=KNOT_VECTORS(PATCH,1,Ni)
        QSIiplus1=KNOT_VECTORS(PATCH,1,Ni+1)            
        IF(QSI1.GE.QSIi.AND.QSI1.LT.QSIiplus1)THEN
            N(I)=1.D0
        ELSE
            N(I)=0.D0      
        ENDIF                       
        Ni=Ni-1    
    ENDDO
   !ZERO ORDER BSPLINES(i+P)
    Ni=N_BF(PATCH,1)+POLYNOMIAL_ORDERS(PATCH,1)
    DO I=1,N_BF(PATCH,1)+POLYNOMIAL_ORDERS(PATCH,1)              
        QSIi=KNOT_VECTORS(PATCH,1,Ni)
        QSIiplus1=KNOT_VECTORS(PATCH,1,Ni+1)        
        IF(QSI1.GE.QSIi.AND.QSI1.LT.QSIiplus1)THEN
            Nplus1(I)=1.D0
        ELSE
            Nplus1(I)=0.D0       
        ENDIF                  
        Ni=Ni-1    
    ENDDO    
!    
!   COX-DE-BOOR RECURSION FORMULA  
!   
    DO P=1,POLYNOMIAL_ORDERS(PATCH,1)
       !P ORDER BSPLINES 
        Ni=N_BF(PATCH,1)
        DO I=1,N_BF(PATCH,1)                               
            QSIi=KNOT_VECTORS(PATCH,1,Ni)
            QSIiplus1=KNOT_VECTORS(PATCH,1,Ni+1)
            QSIiplusP=KNOT_VECTORS(PATCH,1,Ni+P)
            QSIiplusPplus1=KNOT_VECTORS(PATCH,1,Ni+P+1)
            DENOMINATOR1=QSIiplusP-QSIi
            DENOMINATOR2=QSIiplusPplus1-QSIiplus1
            IF(DENOMINATOR1.NE.0.D0.AND.DENOMINATOR2.NE.0.D0)THEN
                N(I)=(QSI1-QSIi)*Nplus1(I+POLYNOMIAL_ORDERS(PATCH,1)-P+1)/DENOMINATOR1+(QSIiplusPplus1-QSI1)*Nplus1(I+POLYNOMIAL_ORDERS(PATCH,1)-P)/DENOMINATOR2
            ENDIF    
            IF(DENOMINATOR1.NE.0.D0.AND.DENOMINATOR2.EQ.0.D0)THEN
                N(I)=(QSI1-QSIi)*Nplus1(I+POLYNOMIAL_ORDERS(PATCH,1)-P+1)/DENOMINATOR1
            ENDIF 
            IF(DENOMINATOR1.EQ.0.D0.AND.DENOMINATOR2.NE.0.D0)THEN
                N(I)=(QSIiplusPplus1-QSI1)*Nplus1(I+POLYNOMIAL_ORDERS(PATCH,1)-P)/DENOMINATOR2
            ENDIF       
            IF(DENOMINATOR1.EQ.0.D0.AND.DENOMINATOR2.EQ.0.D0)THEN
                N(I)=0.D0
            ENDIF                                            
           Ni=Ni-1    
        ENDDO 
        IF(P.LT.POLYNOMIAL_ORDERS(PATCH,1))THEN
           !P ORDER BSPLINES(i+P) 
            Ni=N_BF(PATCH,1)+POLYNOMIAL_ORDERS(PATCH,1)-P  
            DO I=1,N_BF(PATCH,1)+POLYNOMIAL_ORDERS(PATCH,1)-P                                       
                QSIi=KNOT_VECTORS(PATCH,1,Ni)
                QSIiplus1=KNOT_VECTORS(PATCH,1,Ni+1)
                QSIiplusP=KNOT_VECTORS(PATCH,1,Ni+P)
                QSIiplusPplus1=KNOT_VECTORS(PATCH,1,Ni+P+1)
                DENOMINATOR1=QSIiplusP-QSIi
                DENOMINATOR2=QSIiplusPplus1-QSIiplus1
                IF(DENOMINATOR1.NE.0.D0.AND.DENOMINATOR2.NE.0.D0)THEN
                    Nplus1(I)=(QSI1-QSIi)*Nplus1(I+1)/DENOMINATOR1+(QSIiplusPplus1-QSI1)*Nplus1(I)/DENOMINATOR2
                ENDIF    
                IF(DENOMINATOR1.NE.0.D0.AND.DENOMINATOR2.EQ.0.D0)THEN
                    Nplus1(I)=(QSI1-QSIi)*Nplus1(I+1)/DENOMINATOR1
                ENDIF 
                IF(DENOMINATOR1.EQ.0.D0.AND.DENOMINATOR2.NE.0.D0)THEN
                    Nplus1(I)=(QSIiplusPplus1-QSI1)*Nplus1(I)/DENOMINATOR2
                ENDIF       
                IF(DENOMINATOR1.EQ.0.D0.AND.DENOMINATOR2.EQ.0.D0)THEN
                    Nplus1(I)=0.D0
                ENDIF                                            
               Ni=Ni-1    
            ENDDO
        ENDIF             
    ENDDO
    IF(QSI1.EQ.KNOT_VECTORS(PATCH,1,N_BF(PATCH,1)+POLYNOMIAL_ORDERS(PATCH,1)+1))THEN
        N(1)=1.D0    
    ENDIF            
!
! M UNIVARIATE BSPLINE BASIS FUNCTIONS 
!  
   !ZERO ORDER BSPLINES
    Ni=N_BF(PATCH,2)
    DO I=1,N_BF(PATCH,2)      
        QSIi=KNOT_VECTORS(PATCH,2,Ni)
        QSIiplus1=KNOT_VECTORS(PATCH,2,Ni+1)             
        IF(QSI2.GE.QSIi.AND.QSI2.LT.QSIiplus1)THEN
            M(I)=1.D0
        ELSE
            M(I)=0.D0      
        ENDIF                       
        Ni=Ni-1    
    ENDDO
   !ZERO ORDER BSPLINES(i+P)
    Ni=N_BF(PATCH,2)+POLYNOMIAL_ORDERS(PATCH,2)
    DO I=1,N_BF(PATCH,2)+POLYNOMIAL_ORDERS(PATCH,2)             
        QSIi=KNOT_VECTORS(PATCH,2,Ni)
        QSIiplus1=KNOT_VECTORS(PATCH,2,Ni+1)        
        IF(QSI2.GE.QSIi.AND.QSI2.LT.QSIiplus1)THEN
            Mplus1(I)=1.D0
        ELSE
            Mplus1(I)=0.D0       
        ENDIF                  
        Ni=Ni-1    
    ENDDO    
   !COX-DE-BOOR RECURSION FORMULA    
    DO P=1,POLYNOMIAL_ORDERS(PATCH,2)
       !P ORDER BSPLINES 
        Ni=N_BF(PATCH,2)
        DO I=1,N_BF(PATCH,2)                                  
            QSIi=KNOT_VECTORS(PATCH,2,Ni)
            QSIiplus1=KNOT_VECTORS(PATCH,2,Ni+1)
            QSIiplusP=KNOT_VECTORS(PATCH,2,Ni+P)
            QSIiplusPplus1=KNOT_VECTORS(PATCH,2,Ni+P+1)
            DENOMINATOR1=QSIiplusP-QSIi
            DENOMINATOR2=QSIiplusPplus1-QSIiplus1
            IF(DENOMINATOR1.NE.0.D0.AND.DENOMINATOR2.NE.0.D0)THEN
                M(I)=(QSI2-QSIi)*Mplus1(I+POLYNOMIAL_ORDERS(PATCH,2)-P+1)/DENOMINATOR1+(QSIiplusPplus1-QSI2)*Mplus1(I+POLYNOMIAL_ORDERS(PATCH,2)-P)/DENOMINATOR2
            ENDIF    
            IF(DENOMINATOR1.NE.0.D0.AND.DENOMINATOR2.EQ.0.D0)THEN
                M(I)=(QSI2-QSIi)*Mplus1(I+POLYNOMIAL_ORDERS(PATCH,2)-P+1)/DENOMINATOR1
            ENDIF 
            IF(DENOMINATOR1.EQ.0.D0.AND.DENOMINATOR2.NE.0.D0)THEN
                M(I)=(QSIiplusPplus1-QSI2)*Mplus1(I+POLYNOMIAL_ORDERS(PATCH,2)-P)/DENOMINATOR2
            ENDIF       
            IF(DENOMINATOR1.EQ.0.D0.AND.DENOMINATOR2.EQ.0.D0)THEN
                M(I)=0.D0
            ENDIF                                            
           Ni=Ni-1    
        ENDDO 
        IF(P.LT.POLYNOMIAL_ORDERS(PATCH,2))THEN
           !P ORDER BSPLINES(i+P) 
            Ni=N_BF(PATCH,2)+POLYNOMIAL_ORDERS(PATCH,2)-P  
            DO I=1,N_BF(PATCH,2)+POLYNOMIAL_ORDERS(PATCH,2)-P                                     
                QSIi=KNOT_VECTORS(PATCH,2,Ni)
                QSIiplus1=KNOT_VECTORS(PATCH,2,Ni+1)
                QSIiplusP=KNOT_VECTORS(PATCH,2,Ni+P)
                QSIiplusPplus1=KNOT_VECTORS(PATCH,2,Ni+P+1) 
                DENOMINATOR1=QSIiplusP-QSIi
                DENOMINATOR2=QSIiplusPplus1-QSIiplus1
                IF(DENOMINATOR1.NE.0.D0.AND.DENOMINATOR2.NE.0.D0)THEN
                    Mplus1(I)=(QSI2-QSIi)*Mplus1(I+1)/DENOMINATOR1+(QSIiplusPplus1-QSI2)*Mplus1(I)/DENOMINATOR2
                ENDIF    
                IF(DENOMINATOR1.NE.0.D0.AND.DENOMINATOR2.EQ.0.D0)THEN
                    Mplus1(I)=(QSI2-QSIi)*Mplus1(I+1)/DENOMINATOR1
                ENDIF 
                IF(DENOMINATOR1.EQ.0.D0.AND.DENOMINATOR2.NE.0.D0)THEN
                    Mplus1(I)=(QSIiplusPplus1-QSI2)*Mplus1(I)/DENOMINATOR2
                ENDIF       
                IF(DENOMINATOR1.EQ.0.D0.AND.DENOMINATOR2.EQ.0.D0)THEN
                    Mplus1(I)=0.D0
                ENDIF                                            
               Ni=Ni-1    
            ENDDO
        ENDIF             
    ENDDO  
    IF(QSI2.EQ.KNOT_VECTORS(PATCH,2,N_BF(PATCH,2)+POLYNOMIAL_ORDERS(PATCH,2)+1))THEN
        M(1)=1.D0    
    ENDIF         
!
!   BIVARIATE BASIS FUNCTIONS 
!              
    PHI=0.D0 
    GLOB_NUM=0 
    DO J=0,N_BF(PATCH,2)-1
        DO I=0,N_BF(PATCH,1)-1
            GLOB_NUM=GLOB_NUM+1
            CP=CONTROLPOINTS_CONNECTIVITY(PATCH,GLOB_NUM)
            PHI(GLOB_NUM)=N(N_BF(PATCH,1)-I)*M(N_BF(PATCH,2)-J)*COORD_CONTROLPOINTS(CP,4)
            SOMA=SOMA+PHI(GLOB_NUM)
        ENDDO
    ENDDO
    DO I=1,N_GBF(PATCH) 
        PHI(I)=PHI(I)/SOMA
    ENDDO
!             
    DO I=1,N_GBF(PATCH)
        X1=X1+PHI(I)*VALUES(I,1)
        X2=X2+PHI(I)*VALUES(I,2)
        X3=X3+PHI(I)*VALUES(I,3) 
    ENDDO
!    
    END SUBROUTINE PHYSICAL_MAPPING       