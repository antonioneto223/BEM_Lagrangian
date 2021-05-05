SUBROUTINE FUNCOES_DE_FORMA_MEC_MEF (FLAG,A,B,ELEM,VALORES,X,Y,Z,PHI)
    
!
!   A = COORDENADA ADMENSIONAL
!   B = NUMBER OS NODES OF THE ELEMENT
!   ELEM = NUMBER OF ELEMENT IN THE REINFORCEMENT MESH
!    
!	IN THIS SUBROUTINE, THE USER MUST INFORM THE QSI, THE ORDER OF THE 
!	ELEMENT, THE VALUES TO BE INTERPOLATED (VALORES) AND THE RESULTS ARE STORED IN
!	X, Y, Z AND PHI VARIABLES
!
!	IF FLAG=1, A VALUE IS INTERPOLATED USING THE SHAPE FUNCTIONS
!	THE RESULTS ARE STORED IN X AND Y AND Z VARIABLES.
!	IF FLAG=2, THE VALUES OF THE SHAPE FUNCTIONS ARE CALCULATED FOR A 
!	GIVEN ADIMENSIONAL COORDINATE. THE RESULTS ARE STORED IN PHI VARIABLE
!	IF FLAG=3, THE DERIVATIVES OF THE SHAPE FUNCTIONS AT A GIVEN ADIMENSIONAL
!	COORDINATE ARE CALCULATED.THE RESULTS ARE STORED IN PHI VARIABLE
!
!	IF FLAG=4 THE VALUES OF THE SHAPE FUNCTIONS ARE CALCULATED FOR A 
!   GIVEN ADIMENSIONAL COORDINATE, CONSIDERING GEOMETRIC NODES. RESULTS SROTED IN PHI VARIABLE
!   IF FLAG=5, A VALUE IS INTERPOLATED USING THE SHAPE FUNCTIONS CONSIDERING GEOMETRIC NODES
!   THE RESULTS ARE STORED IN X AND Y AND Z VARIABLES.
!   IF FLAG=6, THE DERIVATIVES OF THE SHAPE FUNCTIONS AT A GIVEN ADMIENSIONAL COORDINATE ARE CALCULATED, 
!   CONSIDERING GEOMETRIC NODES. THE RESULTS ARE STORED IN PHI VARIABLE.
!
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,J,K,B,FLAG,APONTA,ELEM
    		 
	REAL*8::A,VALORES(B,3),X,Y,Z,AUX1,AUX2,AUX3,AUX4,EEI,EEJ,PHI(B),PASSO
    
    IF (FLAG .EQ. 1) THEN !_______________________________________________________________________________________
        
        X=0.0D0
        Y=0.0D0
        Z=0.0D0
        
        DO I=1,B
			AUX1=1.D0
			EEI=QSI_MEF(ELEM,I)
			DO J=1,B
				EEJ=QSI_MEF(ELEM,j)
				IF(I.NE.J) THEN
					AUX1=AUX1*((A-EEJ)/(EEI-EEJ))
				ENDIF
			ENDDO
			X=X+AUX1*VALORES(I,1)
			Y=Y+AUX1*VALORES(I,2)
            Z=Z+AUX1*VALORES(I,3)
		ENDDO
        
    ELSE IF (FLAG .EQ. 2) THEN !__________________________________________________________________________________
        
        DO I=1,B
			AUX1=1.D0
			EEI=QSI_MEF(ELEM,I)
			DO J=1,B
				EEJ=QSI_MEF(ELEM,J)
				IF(I.NE.J) THEN
					AUX1=AUX1*((A-EEJ)/(EEI-EEJ))
				ENDIF
			ENDDO
			PHI(I)=AUX1
        ENDDO
        
    ELSE IF (FLAG .EQ. 3) THEN !__________________________________________________________________________________
        
        PHI=0.D0

		DO K=1,B
			AUX3=0.D0
			APONTA=0
			EEI=QSI_MEF(ELEM,K)
			DO I=1,(B-1)
				AUX1=1.D0
				APONTA=APONTA+1
				IF(APONTA.EQ.K) THEN
					APONTA=APONTA+1
				ENDIF

				DO J=1,B
					EEJ=QSI_MEF(ELEM,J)
					IF(K.NE.J) THEN
						IF(J.EQ.APONTA) THEN
							AUX4=1.D0
						ELSE
							AUX4=A-EEJ
						ENDIF
						AUX2=((AUX4)/(EEI-EEJ))
						AUX1=AUX1*AUX2
					ENDIF
				ENDDO
				AUX3=AUX3+AUX1
			ENDDO
			PHI(K)=AUX3
		ENDDO
        
    ELSE IF (FLAG .EQ. 4) THEN !__________________________________________________________________________________
        
        PHI=0.D0
		PASSO=((2.D0)/(B-1))
        
        DO I=1,B
			AUX1=1.D0
			EEI=-1.D0+(PASSO*(I-1))
			DO J=1,B
				EEJ=-1.D0+(PASSO*(J-1))
				IF(I.NE.J) THEN
					AUX1=AUX1*((A-EEJ)/(EEI-EEJ))
				ENDIF
			ENDDO
			PHI(I)=AUX1
        ENDDO
        
    ELSE IF (FLAG .EQ. 5) THEN !__________________________________________________________________________________
        
        X=0.D0
		Y=0.D0
        Z=0.0D0
		PASSO=((2.D0)/(B-1))

		DO I=1,B
			AUX1=1.D0
			EEI=-1.D0+(PASSO*(I-1))
			DO J=1,B
				EEJ=-1.D0+(PASSO*(J-1))
				IF(I.NE.J) THEN
					AUX1=AUX1*((A-EEJ)/(EEI-EEJ))
				ENDIF
			ENDDO
			X=X+AUX1*VALORES(I,1)
			Y=Y+AUX1*VALORES(I,2)
            Z=Z+AUX1*VALORES(I,3)
        ENDDO
            
    ELSE IF (FLAG .EQ. 6) THEN !__________________________________________________________________________________
        				
        PHI=0.D0
		PASSO=((2.D0)/(B-1))

		DO K=1,B
			AUX3=0.D0
			APONTA=0
			EEI=-1.D0+(PASSO*(K-1))
			DO I=1,(B-1)
				AUX1=1.D0
				APONTA=APONTA+1
				IF(APONTA.EQ.K) THEN
					APONTA=APONTA+1
				ENDIF

				DO J=1,B
					EEJ=-1.D0+(PASSO*(J-1))
					IF(K.NE.J) THEN
						IF(J.EQ.APONTA) THEN
							AUX4=1.D0
						ELSE
							AUX4=A-EEJ
						ENDIF
						AUX2=((AUX4)/(EEI-EEJ))
						AUX1=AUX1*AUX2
					ENDIF
				ENDDO
				AUX3=AUX3+AUX1
			ENDDO
			PHI(K)=AUX3
        ENDDO
    
    ENDIF !_______________________________________________________________________________________________________
    
    
END SUBROUTINE ! _________________________________________________________________________________________________
    
    
    
    
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  
    

    
    
 SUBROUTINE FUNCOES_FORMA_MEF_P_INT(I,QSI,NUM_P_INT,QSI_P_INT,PHI)
     
     !
     !SUBROTINA QUE CALCULA OS VALORES DAS FUNCOES DE FORMA NO PONTO QSI PARA O ELEMENTO FINITO DEFINIDO PELOS PONTOS DE INTEGRACAO
     !USADA PARA INTERPOLACAO DOS VALORES DE TENSAO NO PROCESSO ELASTOPLASTICO
     !I: ELEMENTO EM QUESTAO
     !QSI: VALOR DA COORDENADA ADMENSIONAL DO PONTO ONDE VAI CALCULAR AS FUNCOES
     !PHI: VETOR VETOR COM OS VALORES DAS FUNCOES DE FORMA NO PONTO QSI (SAIDA)
     !NUM_P_INT: NUMERO DE PONTOS DE INTEGRACAO NO ELEMENTO
     !QSI_P_INT: VETOR COORDENADAS ADIMENSIONAIS DOS PONTOS DE INTEGRACAO
     !
 
	  USE REINFORCEMENTS
	  
	  IMPLICIT NONE
	  
      INTEGER,INTENT(IN)::I,NUM_P_INT
      REAL*8,INTENT(IN)::QSI,QSI_P_INT(NUM_P_INT)
      REAL*8,INTENT(OUT)::PHI(NUM_P_INT)
      
	  INTEGER::J,JJ
	  REAL*8::AUX1,EEI,EEJ
	  
        DO J=1,NUM_P_INT
	      AUX1=1.D0
			EEI=QSI_P_INT(J)
			DO JJ=1,NUM_P_INT
				EEJ=QSI_P_INT(JJ)
				IF(J.NE.JJ) THEN
					AUX1=AUX1*((QSI-EEJ)/(EEI-EEJ))
				ENDIF
			ENDDO
			PHI(J)=AUX1
		ENDDO
	
	  
END SUBROUTINE FUNCOES_FORMA_MEF_P_INT      
    
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  
   
    
    
    
SUBROUTINE D_FUNCOES_FORMA_MEF_P_INT(I,QSI,NUM_P_INT,QSI_P_INT,D_PHI)

    !SUBROTINA QUE CALCULA OS VALORES DAS DERIVADAS DAS FUNCOES DE FORMA NO PONTO QSI PARA O ELEMENTO FINITO DEFINIDO PELOS PONTOS DE INTEGRACAO
    !USADA PARA INTERPOLACAO DOS VALORES DE TENSAO NO PROCESSO ELASTOPLASTICO
    !I: ELEMENTO EM QUESTAO
    !QSI: VALOR DA COORDENADA ADMENSIONAL DO PONTO ONDE VAI CALCULAR AS FUNCOES
    !D_PHI: VETOR COM OS VALORES DAS DERIVADAS DAS FUNCOES DE FORMA NO PONTO QSI (SAIDA)
    !NUM_P_INT: NUMERO DE PONTOS DE INTEGRACAO NO ELEMENTO
    !QSI_P_INT: COORDENADAS ADIMENSIONAIS DOS PONTOS DE INTEGRACAO
	  
	    USE REINFORCEMENTS
	    
	    IMPLICIT NONE
	    
        INTEGER,INTENT(IN)::I,NUM_P_INT
        REAL*8,INTENT(IN)::QSI,QSI_P_INT(NUM_P_INT)
        REAL*8,INTENT(OUT)::D_PHI(NUM_P_INT)
        
	    INTEGER::J,JJ,APONTA,II,K
	    REAL*8::AUX1,EEI,EEJ,AUX3,AUX4,AUX2
	    
	    D_PHI=0.D0
        DO K=1,NUM_P_INT
	        AUX3=0.D0
	        APONTA=0
	        EEI=QSI_P_INT(K)
	        DO II=1,(NUM_P_INT-1)
		        AUX1=1.D0
		        APONTA=APONTA+1
		        IF(APONTA.EQ.K) THEN
			        APONTA=APONTA+1
		        ENDIF

		        DO J=1,NUM_P_INT
      		        EEJ=QSI_P_INT(J)
			        IF(K.NE.J) THEN
				        IF(J.EQ.APONTA) THEN
					        AUX4=1.D0
				        ELSE
					        AUX4=QSI-EEJ
				        ENDIF
				        AUX2=((AUX4)/(EEI-EEJ))
				        AUX1=AUX1*AUX2
			        ENDIF
		        ENDDO
		        AUX3=AUX3+AUX1
	        ENDDO
	        D_PHI(K)=AUX3
        ENDDO
	 
END SUBROUTINE
    
        
        
    
    
    
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  
    

    
 
SUBROUTINE UPDATE_VISCO_VARIABLES(TIME)
    !
    ! SUBROUTINE USED IN ANALYZES (VISCO CASE) TO UPDATE ALL VISCOELASTIC VARIABLES AT EACH TIME 
    !

    USE VISCO
    USE REINFORCEMENTS
    USE ANALYSIS
    
    IMPLICIT NONE
    
    INTEGER,INTENT(IN)::TIME
    INTEGER::MAX_ORDER,I
    
     IF (TIME .EQ. 1) THEN
            Fs=0.D0
            U_ANTERIOR=0.0D0
            T_ANTERIOR=0.0D0
            U_INT_ANTERIOR=0.0D0
            S_ANTERIOR=0.0D0
            S_EL_ANTERIOR=0.0D0
            
            IF(EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
                U_ANTERIOR_MEF = 0.0D0
                P_ANTERIOR_MEF = 0.0D0
                IF (VISCO_ANALYSIS_REINF) THEN
                    MAX_ORDER=0
                    DO I=1,N_ELEMENTOS_MEF
                        IF (ORDEM_MEF(I).GT.MAX_ORDER) MAX_ORDER=ORDEM_MEF(I)
                    END DO
                    ALLOCATE(E_NORMAL_ANTERIOR(N_ELEMENTOS_MEF,MAX_ORDER+1),&
                        E_NORMAL_EL_ANTERIOR(N_ELEMENTOS_MEF,MAX_ORDER+1))
                    E_NORMAL_ANTERIOR=0.0D0
                    E_NORMAL_EL_ANTERIOR=0.0D0
                ENDIF
            ENDIF
        ELSE
            U_ANTERIOR = U
            T_ANTERIOR = T
            U_INT_ANTERIOR=U_INT_VISCO
            S_ANTERIOR = S_INT_VISCO
            S_EL_ANTERIOR = S_EL
            IF(EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
                U_ANTERIOR_MEF = U_MEF
                P_ANTERIOR_MEF = P_MEF
                IF (VISCO_ANALYSIS_REINF) THEN
                    E_NORMAL_ANTERIOR=E_NORMAL
                    E_NORMAL_EL_ANTERIOR=E_NORMAL_EL
                ENDIF
            ENDIF
        ENDIF
        
END SUBROUTINE