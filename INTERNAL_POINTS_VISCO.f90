SUBROUTINE INTERNAL_POINTS_VISCO(TIME_STEP)
    
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
	USE REMESHING
	USE PROPAGATION
    USE REINFORCEMENTS
    USE VISCO
    
    IMPLICIT NONE
    
    INTERFACE 
        SUBROUTINE INTERNAL_POINTS_MATRIX(H_IC,G_IC,COORDS,PROP,DOMAIN)
            INTEGER::DOMAIN    
            REAL*8,DIMENSION(:,:)::H_IC,G_IC
            REAL*8::COORDS(3),PROP(5)
        END SUBROUTINE
    END INTERFACE
    INTERFACE 
        SUBROUTINE INTERNAL_POINTS_MATRIX_HP(H_IC,G_IC,COORDS,PROP,DOMAIN)
            INTEGER::DOMAIN    
            REAL*8,DIMENSION(:,:)::H_IC,G_IC
            REAL*8::COORDS(3),PROP(5)
        END SUBROUTINE
    END INTERFACE
    

    INTEGER,INTENT(IN)::TIME_STEP
    LOGICAL::INT_2D
    INTEGER::I,J,M,K,KK,II,JJ,DOMAIN,NUM_LEI,P_INTEGRACAO
    REAL*8::GAMA,Eve,Ee,K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,VET_COORD(3),VET_PROP(5),PI,CONST,Nu,Mu
    REAL*8,DIMENSION(:),ALLOCATABLE::U_PONTO,T_PONTO,P_PONTO_MEF
    REAL*8,DIMENSION(:,:),ALLOCATABLE::DG_CF,QSIW
            
    ! STARTING LOCAL VARIABLES
    ALLOCATE(U_PONTO(3*N_COLLOCPOINTS),T_PONTO(3*N_COLLOCPOINTS))
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') ALLOCATE(P_PONTO_MEF(3*N_NOS_MEF))
    
    IF (TIME_STEP.EQ.1) THEN    ! FIRST TIME STEP, CALCULATE ALL THE MATRIXES
        
        ! STARTING GLOBAL VARIABLES
        II=N_COLLOCPOINTS
        JJ=N_NODES_INT_VISCO
        ALLOCATE(HH1(3*JJ,3*II),GG1(3*JJ,3*II),HHH1(9*JJ,3*II),GGG1(9*JJ,3*II),HH2(3*JJ,3*II),&
            GG2(3*JJ,3*II),HHH2(9*JJ,3*II),GGG2(9*JJ,3*II))    
        IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
            II=N_NOS_MEF
            JJ=N_NODES_INT_VISCO
            ALLOCATE(GG1F(3*JJ,3*II),GG2F(3*JJ,3*II),GGG1F(9*JJ,3*II),GGG2F(9*JJ,3*II))
            GG1F=0.0D0
            GG2F=0.0D0
            GGG1F=0.0D0
            GGG2F=0.0D0
        ENDIF
        
        ! STARTING LOCAL INTEGRATION VARIABLES
        PI=DACOS(-1.D0)
        JJ=N_COLLOCPOINTS
        
        DO I=1,N_NODES_INT_VISCO
            
            ! Parametros viscoelasticos
            NUM_LEI = TIME_MODEL(NDI_VISCO(I))
            GAMA = EMP_VISCO(NDI_VISCO(I),1)
            Eve  = EMP_VISCO(NDI_VISCO(I),2)
            Ee   = EMP(NDI_VISCO(I),1)
            SELECT CASE(NUM_LEI)
            CASE(1)!KELVIN'S MODEL
                K1= -(1.0D0 + GAMA/DELTA_T_VISCO)
                K2=  1.0D0
                K3=  GAMA/DELTA_T_VISCO
                K4=  0.0D0
                K5=  1.0D0
                K6= -1.0D0
                K7=  0.0D0
                K8= -GAMA
                K9=  GAMA/DELTA_T_VISCO
                K10= 1.0D0/(1.0D0+GAMA/DELTA_T_VISCO)
                K11= 0.0D0
                K12= 1.0D0 
            CASE(2)!BOLTZMANN'S MODEL
	            K1= -(1.0D0 + GAMA/DELTA_T_VISCO)
                K2=  GAMA/DELTA_T_VISCO + (Ee + Eve)/Eve
                K3=  GAMA/DELTA_T_VISCO
                K4= -GAMA/DELTA_T_VISCO
                CONST= 1.D0/(1.0D0 + GAMA*Eve/(DELTA_T_VISCO*(Ee+Eve)))
                K5=  CONST
                K6= -CONST*Eve/(Ee+Eve)
                K7=  CONST*GAMA*Eve/(Ee+Eve)
                K8= -CONST*GAMA*Eve/(Ee+Eve)
                K9= GAMA/DELTA_T_VISCO
                K10= 1.0D0/(1.0D0+GAMA/DELTA_T_VISCO)
                K11= GAMA*Eve/(Ee+Eve)/DELTA_T_VISCO
                K12= CONST
            CASE(3)!MAXWELL'S MODEL
   	            K1= -(1.0D0)
                K2=  1.0D0 + DELTA_T_VISCO/GAMA
                K3=  1.0D0
                K4= -1.0D0
                K5=  DELTA_T_VISCO/GAMA
                K6=  0.0D0
                K7=  DELTA_T_VISCO
                K8= -DELTA_T_VISCO
                K9=  1.0D0
                K10= 1.0D0
                K11= 1.0D0
                K12= 1/(1.0D0 + DELTA_T_VISCO/GAMA) 
            CASE(4)!HOOKE'S MODEL
   	            K1= -(1.0D0)
                K2= 1.0D0
                K3= 0.0D0
                K4= 0.0D0
                K5= 1.0D0
                K6=-1.0D0
                K7= 0.0D0
                K8= 0.0D0
                K9= 0.0D0
                K10=1.0D0
                K11=0.0D0
                K12=1.0D0 
            ENDSELECT
            
            ! Parametros para integracao
            DOMAIN=NDI_VISCO(I)
            VET_PROP(1)=EMP(DOMAIN,1)/(2*(1+EMP(DOMAIN,2)))                 !Mu
	        VET_PROP(2)=EMP(DOMAIN,2)                                       !Nu
	        VET_PROP(3)=((1.D0)/(16.D0*PI*VET_PROP(1)*(1.D0-VET_PROP(2))))  !C1
	        VET_PROP(4)=((-1.D0)/(8.D0*PI*(1.D0-VET_PROP(2))))              !C2
	        VET_PROP(5)=VET_PROP(1)/(4.0D0*PI*(1.D0-VET_PROP(2)))           !C3
            VET_COORD(1)=COORDS_INT_VISCO(I,1)
            VET_COORD(2)=COORDS_INT_VISCO(I,2)
            VET_COORD(3)=COORDS_INT_VISCO(I,3)
    
            ! Realizando integracoes
            CALL INTERNAL_POINTS_MATRIX(HH1(3*I-2:3*I,:),GG1(3*I-2:3*I,:),VET_COORD,VET_PROP,DOMAIN)
            CALL INTERNAL_POINTS_MATRIX_HP(HHH1(9*I-8:9*I,:),GGG1(9*I-8:9*I,:),VET_COORD,VET_PROP,DOMAIN)
            
            ! Atualizando matrizes globais com coeficientes
            DO K=1,3
                !HH1(3*(I-1)+K,:) = HH_AUX(K,:)*K1
                !GG1(3*(I-1)+K,:) = GG_AUX(K,:)*K2
                !HH2(3*(I-1)+K,:) = HH_AUX(K,:)*K3
                !GG2(3*(I-1)+K,:) = GG_AUX(K,:)*K4
                HH2(3*(I-1)+K,:) = HH1(3*(I-1)+K,:)*K3
                HH1(3*(I-1)+K,:) = HH1(3*(I-1)+K,:)*K1
                GG2(3*(I-1)+K,:) = GG1(3*(I-1)+K,:)*K4
                GG1(3*(I-1)+K,:) = GG1(3*(I-1)+K,:)*K2
            ENDDO
            DO K=1,9
                !HHH1(9*(I-1)+K,:) = HHH_AUX(K,:)*K6
                !GGG1(9*(I-1)+K,:) = GGG_AUX(K,:)*K5
                !HHH2(9*(I-1)+K,:) = HHH_AUX(K,:)*K8
                !GGG2(9*(I-1)+K,:) = GGG_AUX(K,:)*K7
                HHH2(9*(I-1)+K,:) = HHH1(9*(I-1)+K,:)*K8
                HHH1(9*(I-1)+K,:) = HHH1(9*(I-1)+K,:)*K6
                GGG2(9*(I-1)+K,:) = GGG1(9*(I-1)+K,:)*K7
                GGG1(9*(I-1)+K,:) = GGG1(9*(I-1)+K,:)*K5
            ENDDO
            
            IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
                ! Integration points
                P_INTEGRACAO=INT(0.6*N_GAUSSPOINTS)
                IF (P_INTEGRACAO.LT.30) P_INTEGRACAO=30
                ALLOCATE(QSIW(P_INTEGRACAO,2))
                CALL GAUSS_POINTS(P_INTEGRACAO,QSIW)
                        
                DO J=1,N_ELEMENTOS_MEF
                    IF (ND_MEF(J).EQ.NDI_VISCO(I)) THEN
                        !Strating variables
                        ALLOCATE(DG_CF(3,3*(ORDEM_MEF(J)+1)))
                        Nu=EMP(NDI_VISCO(I),2)                          ! POISSON (Nu)
                        Mu=EMP(NDI_VISCO(I),1)/(2*(1+EMP(NDI_VISCO(I),2)))   ! G (Mu)
                                  
                        ! Integration as cylinder or as line for singular equation
                        CALL FIBER_INT_DECISION(INT_2D,J,VET_COORD)
                        IF (INT_2D) THEN                                        
                            CALL FIBER_SING_S(J,Nu,Mu,P_INTEGRACAO,QSIW,VET_COORD,DG_CF)
                        ELSE                                                    
                            CALL FIBER_NONSING_S(J,Nu,Mu,VET_COORD,DG_CF)
                        ENDIF
                        ! Contribuing in the global matrix for J element
                        DO K=1,ORDEM_MEF(J)+1
                            KK=CONECTI_MEF(J,K)
                            DO M=1,3
                                IF (DABS(DG_CF(M,3*K-2)).LT.TOLER) DG_CF(M,3*K-2)=0.0D0
                                IF (DABS(DG_CF(M,3*K-1)).LT.TOLER) DG_CF(M,3*K-1)=0.0D0
                                IF (DABS(DG_CF(M,3*K)).LT.TOLER) DG_CF(M,3*K)=0.0D0
                            
                                GG1F(3*(I-1)+M,3*KK-2)=GG1F(3*(I-1)+M,3*KK-2)+DG_CF(M,3*K-2)*K2
                                GG1F(3*(I-1)+M,3*KK-1)=GG1F(3*(I-1)+M,3*KK-1)+DG_CF(M,3*K-1)*K2
                                GG1F(3*(I-1)+M,3*KK)=GG1F(3*(I-1)+M,3*KK)+DG_CF(M,3*K)*K2
                            
                                GG2F(3*(I-1)+M,3*KK-2)=GG2F(3*(I-1)+M,3*KK-2)+DG_CF(M,3*K-2)*K4
                                GG2F(3*(I-1)+M,3*KK-1)=GG2F(3*(I-1)+M,3*KK-1)+DG_CF(M,3*K-1)*K4
                                GG2F(3*(I-1)+M,3*KK)=GG2F(3*(I-1)+M,3*KK)+DG_CF(M,3*K)*K4
                            ENDDO
                        ENDDO
                        DEALLOCATE(DG_CF)
                        
                        ! Integration as cylinder or as line for singular equation
                        ALLOCATE(DG_CF(9,3*(ORDEM_MEF(J)+1)))
                        IF (INT_2D) THEN                                        
                            !CALL FIBER_SING_HP(J,Nu,Mu,P_INTEGRACAO,QSIW,VET_COORD,DG_CF)
                            CALL FIBER_SING_HP(J,Nu,Mu,P_INTEGRACAO,QSIW,VET_COORD,DG_CF)
                        ELSE                                                    
                            CALL FIBER_NONSING_HP(J,Nu,Mu,VET_COORD,DG_CF)
                        ENDIF
                        ! Contribuing in the global matrix for J element
                        DO K=1,ORDEM_MEF(J)+1
                            KK=CONECTI_MEF(J,K)
                            DO M=1,9
                                IF (DABS(DG_CF(M,3*K-2)).LT.TOLER) DG_CF(M,3*K-2)=0.0D0
                                IF (DABS(DG_CF(M,3*K-1)).LT.TOLER) DG_CF(M,3*K-1)=0.0D0
                                IF (DABS(DG_CF(M,3*K)).LT.TOLER) DG_CF(M,3*K)=0.0D0
                            
                                GGG1F(9*(I-1)+M,3*KK-2)=GGG1F(9*(I-1)+M,3*KK-2)+DG_CF(M,3*K-2)*K5
                                GGG1F(9*(I-1)+M,3*KK-1)=GGG1F(9*(I-1)+M,3*KK-1)+DG_CF(M,3*K-1)*K5
                                GGG1F(9*(I-1)+M,3*KK)=GGG1F(9*(I-1)+M,3*KK)+DG_CF(M,3*K)*K5
                                      
                                GGG2F(9*(I-1)+M,3*KK-2)=GGG2F(9*(I-1)+M,3*KK-2)+DG_CF(M,3*K-2)*K7
                                GGG2F(9*(I-1)+M,3*KK-1)=GGG2F(9*(I-1)+M,3*KK-1)+DG_CF(M,3*K-1)*K7
                                GGG2F(9*(I-1)+M,3*KK)=GGG2F(9*(I-1)+M,3*KK)+DG_CF(M,3*K)*K7
                            ENDDO
                        ENDDO
                        DEALLOCATE(DG_CF)
                    ENDIF
                ENDDO
                
                DEALLOCATE(QSIW)
            ENDIF
                        
            ! ATUALIZANDO AS GRANDEZAS INTERNAS NO PASSO ANTERIOR            
	        U_INT_VISCO(3*I-2)=U_INT_ANTERIOR(3*I-2)*K9
            U_INT_VISCO(3*I-1)=U_INT_ANTERIOR(3*I-1)*K9
	        U_INT_VISCO(3*I)=U_INT_ANTERIOR(3*I)*K9
            DO K=1,9
                S_INT_VISCO(9*(I-1)+K)=S_ANTERIOR(9*(I-1)+K)*K11
            ENDDO      
            
        ENDDO
                
    ELSE                        ! NOT FIRST TIME STEP, MATRIXES ALREADY CALCULATED
        
        DO I=1,N_NODES_INT_VISCO
            
            ! Parametros viscoelasticos
            NUM_LEI = TIME_MODEL(NDI_VISCO(I))
            GAMA = EMP_VISCO(NDI_VISCO(I),1)
            Eve  = EMP_VISCO(NDI_VISCO(I),2)
            Ee   = EMP(NDI_VISCO(I),1)
            SELECT CASE(NUM_LEI)
            CASE(1)!KELVIN'S MODEL
                K9=  GAMA/DELTA_T_VISCO
                K10= 1/(1+GAMA/DELTA_T_VISCO)
                K11= 0.0D0 
                K12= 1.0D0
            CASE(2)!BOLTZMANN'S MODEL
                CONST= 1.D0/(1.0D0 + GAMA*Eve/(DELTA_T_VISCO*(Ee+Eve)))
                K9= GAMA/DELTA_T_VISCO
                K10= 1.0D0/(1.0D0+GAMA/DELTA_T_VISCO)
                K11= GAMA*Eve/(Ee+Eve)/DELTA_T_VISCO
                K12= CONST
            CASE(3)!MAXWELL'S MODEL
                K9=  1.0D0
                K10= 1.0D0
                K11= 1.0D0
                K12= 1/(1.0D0 + DELTA_T_VISCO/GAMA) 
            CASE(4)!HOOKE'S MODEL
                K9=  0.0D0
                K10= 1.0D0
                K11= 0.0D0
                K12= 1.0D0
            ENDSELECT
            
            ! ATUALIZANDO AS GRANDEZAS INTERNAS NO PASSO ANTERIOR            
	        U_INT_VISCO(3*I-2)=U_INT_ANTERIOR(3*I-2)*K9
            U_INT_VISCO(3*I-1)=U_INT_ANTERIOR(3*I-1)*K9
	        U_INT_VISCO(3*I)=U_INT_ANTERIOR(3*I)*K9
            DO K=1,9
                S_INT_VISCO(9*(I-1)+K)=S_ANTERIOR(9*(I-1)+K)*K11
            ENDDO   
            
        ENDDO
        
    ENDIF
    
    ! OBTAINING INTERNAL POINTS DISPLACEMENTS
    U_INT_VISCO = U_INT_VISCO + MATMUL(HH1,U)
    U_INT_VISCO = U_INT_VISCO + MATMUL(GG1,T)
    U_INT_VISCO = U_INT_VISCO + MATMUL(HH2,U_ANTERIOR)
    U_INT_VISCO = U_INT_VISCO + MATMUL(GG2,T_ANTERIOR)
    
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        U_INT_VISCO = U_INT_VISCO + MATMUL(GG1F,P_MEF)
        U_INT_VISCO = U_INT_VISCO + MATMUL(GG2F,P_ANTERIOR_MEF)
    ENDIF
        
    ! OBTAINING DISPLACEMENTS AND TRACTIONS TIME DERIVATIVES 
    T_PONTO = (T - T_ANTERIOR)/DELTA_T_VISCO
    U_PONTO = (U - U_ANTERIOR)/DELTA_T_VISCO
    
    ! OBTAINING INTERNAL POINTS STRESSES 
    S_INT_VISCO = S_INT_VISCO + MATMUL(GGG1,T)
    S_INT_VISCO = S_INT_VISCO + MATMUL(HHH1,U)
    S_INT_VISCO = S_INT_VISCO + MATMUL(GGG2,T_PONTO)
    S_INT_VISCO = S_INT_VISCO + MATMUL(HHH2,U_PONTO)
    
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        P_PONTO_MEF = (P_MEF - P_ANTERIOR_MEF)/DELTA_T_VISCO    
    
        S_INT_VISCO = S_INT_VISCO + MATMUL(GGG1F,P_MEF)
        S_INT_VISCO = S_INT_VISCO + MATMUL(GGG2F,P_PONTO_MEF)
    ENDIF
    
    DO I=1,N_NODES_INT_VISCO
        
        ! Parametros viscoelasticos
        NUM_LEI = TIME_MODEL(NDI_VISCO(I))
        GAMA = EMP_VISCO(NDI_VISCO(I),1)
        Eve  = EMP_VISCO(NDI_VISCO(I),2)
        Ee   = EMP(NDI_VISCO(I),1)
        SELECT CASE(NUM_LEI)
        CASE(1)!KELVIN'S MODEL
            K10= 1/(1+GAMA/DELTA_T_VISCO)
            K12= 1.0D0
        CASE(2)!BOLTZMANN'S MODEL
            CONST= 1.D0/(1.0D0 + GAMA*Eve/(DELTA_T_VISCO*(Ee+Eve)))
            K10= 1.0D0/(1.0D0+GAMA/DELTA_T_VISCO)
            K12= CONST
        CASE(3)!MAXWELL'S MODEL
            K10= 1.0D0
            K12= 1/(1.0D0 + DELTA_T_VISCO/GAMA) 
        CASE(4)!HOOKE'S MODEL
            K10= 1.0D0
            K12= 1.0D0
        ENDSELECT
        
        U_INT_VISCO(3*I-2)=U_INT_VISCO(3*I-2)*K10
        U_INT_VISCO(3*I-1)=U_INT_VISCO(3*I-1)*K10
	    U_INT_VISCO(3*I)=U_INT_VISCO(3*I)*K10
        DO K=1,9
            S_INT_VISCO(9*(I-1)+K)=S_INT_VISCO(9*(I-1)+K)*K12
        ENDDO  
        
    ENDDO
    
    ! OBTAINING ELASTIC AND VISCO PORTIONS OF STRESS
    DO I=1,N_NODES_INT_VISCO 

        ! Parametros viscoelasticos
        NUM_LEI = TIME_MODEL(NDI_VISCO(I))
        GAMA = EMP_VISCO(NDI_VISCO(I),1)
        Eve  = EMP_VISCO(NDI_VISCO(I),2)
        Ee   = EMP(NDI_VISCO(I),1)

        SELECT CASE(NUM_LEI)
        CASE(1)!KELVIN'S MODEL
            K1= GAMA/DELTA_T_VISCO
            K2= 1.0D0
        CASE(2)!BOLTZMANN'S MODEL
	        K1= GAMA/DELTA_T_VISCO
            K2= 1.0D0
        CASE(3)!MAXWELL'S MODEL
   	        K1= 0.0D0
            K2= 2.0D0
        CASE(4)!HOOKE'S MODEL
   	        K1= 0.0D0
            K2= 1.0D0
        ENDSELECT

        DO K=1,9
            S_EL(9*(I-1)+K) = (S_INT_VISCO(9*(I-1)+K)+K1*S_EL_ANTERIOR(9*(I-1)+K))/(1.0D0+K1)    
        ENDDO
        !S_EL(4*I-3)=(S_INT_VISCO(4*I-3)+K1*S_EL_ANTERIOR(4*I-3))/(1.0D0+K1)
        !S_EL(4*I-2)=(S_INT(4*I-2)+K1*S_EL_ANTERIOR(4*I-2))/(1.0D0+K1)
        !S_EL(4*I-1)=(S_INT(4*I-1)+K1*S_EL_ANTERIOR(4*I-1))/(1.0D0+K1)
        !S_EL(4*I  )=(S_INT(4*I  )+K1*S_EL_ANTERIOR(4*I  ))/(1.0D0+K1)
        
        DO K=1,9
            S_VISCOSA(9*(I-1)+K) = K2*S_INT_VISCO(9*(I-1)+K)-S_EL(9*(I-1)+K)
        ENDDO

    ENDDO

    DEALLOCATE(U_PONTO,T_PONTO)
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') DEALLOCATE(P_PONTO_MEF)
    
END SUBROUTINE