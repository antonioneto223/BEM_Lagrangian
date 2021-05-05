SUBROUTINE REINFORCEMENTS_EP_PROCESS    

    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE REINFORCEMENTS
    USE OMP_LIB
    
    IMPLICIT NONE
    
    INTEGER::I,J,K,II,JJ,MAX_INT_P_MEF,PASSO,ITER_EP,contelem_pnod[allocatable](:)
    
    REAL*8::ERRO,AUX1,AUX2,AUX3
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::CRIT,D_ALFA
    
    ! Finding max number of integration points in one FE
    ! Finding total number of local element's nodes (for the normal effort matrixes)
    MAX_INT_P_MEF=0
    TOTAL_ELEM_NODES_MEF=0
    DO I=1,N_ELEMENTOS_MEF
       TOTAL_ELEM_NODES_MEF=TOTAL_ELEM_NODES_MEF+ORDEM_MEF(I)+1
        IF (N_INT_P_MEF (I) .GT. MAX_INT_P_MEF) MAX_INT_P_MEF = N_INT_P_MEF(I) 
    ENDDO
    
    ! Starting variables ----------
    II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
    JJ=3*N_COLLOCPOINTS+6*N_NOS_MEF
    
    ALLOCATE(D_DEF_PL_INTP(N_ELEMENTOS_MEF,MAX_INT_P_MEF),DEF_INTP(N_ELEMENTOS_MEF,MAX_INT_P_MEF),&
    LAMBDA(N_ELEMENTOS_MEF,MAX_INT_P_MEF),ALFA(N_ELEMENTOS_MEF,MAX_INT_P_MEF),D_ALFA(N_ELEMENTOS_MEF,MAX_INT_P_MEF),P_ACUM_PL(3*N_COLLOCPOINTS),&
    U_MEF_ACUM(3*N_NOS_MEF),U_ACUM_PL(3*N_COLLOCPOINTS),SIGMA_INTP(N_ELEMENTOS_MEF,MAX_INT_P_MEF),D_SIGMA_INTP(N_ELEMENTOS_MEF,MAX_INT_P_MEF),&
    D_DEF_INTP(N_ELEMENTOS_MEF,MAX_INT_P_MEF),CRIT(N_ELEMENTOS_MEF,MAX_INT_P_MEF),P_MEF_AT(3*N_NOS_MEF),P_MEF_ACUM(3*N_NOS_MEF),&
    CARGA(3*N_COLLOCPOINTS,2))
    
    ALLOCATE(GRAND_ACUM_NOD(TOTAL_ELEM_NODES_MEF,2)) 
    
    ALLOCATE (SIGMA_NOD(N_ELEMENTOS_MEF,MAXORDER_MEF+1),DEF_PL_NOD(N_ELEMENTOS_MEF,MAXORDER_MEF+1),&
    F_DES_ANT(JJ),F_ALVO(JJ),F_DES(JJ),F_ACUM(JJ))  
    
    ALLOCATE(ALFA_CONECT(N_INTERFACE_MEF,2),SIGMA_INTP_CONECT(N_INTERFACE_MEF,2),D_DEF_PL_INTP_CONECT(N_INTERFACE_MEF,2),LAMBDA_CONECT(N_INTERFACE_MEF,2),&
    D_SIGMA_INTP_CONECT(N_INTERFACE_MEF,2),D_DEF_INTP_CONECT(N_INTERFACE_MEF,2))      
     
    DEF_INTP = 0.0D0
    SIGMA_INTP = 0.0D0
    SIGMA_INTP_CONECT = 0.0D0
    LAMBDA = 0.0D0
    LAMBDA_CONECT = 0.0D0
    ALFA = 0.0D0
    ALFA_CONECT = 0.0D0
    U_MEF_ACUM = 0.0D0
    U_ACUM_PL = 0.0D0
    P_ACUM_PL = 0.0D0
    P_MEF_ACUM = 0.0D0
    F_ALVO = 0.0D0
    GRAND_ACUM_NOD = 0.0D0
    CARGA(:,1) = U(:) / DFLOAT(N_PASSOS)
    CARGA(:,2) = T(:) / DFLOAT(N_PASSOS)
    OPEN(17,FILE='Output_data/REINFORCEMENT_ELASTOPLASTIC_DATA.TXT',STATUS='UNKNOWN')
    ! ---------------------------
    
    DO PASSO=1,N_PASSOS ! Load step loop ________________
        
        U(:) = CARGA(:,1)
        T(:) = CARGA(:,2)
        F_DES_ANT = 0.0D0
        F_DES = 0.0D0
        F_ACUM = 0.0D0
    
        ITER_EP = 0
        ERRO = 10.0D0
        WRITE(17,*)''
        WRITE(17,'(a,i3,a)')'____________________|| PASSO DE CARGA:',passo,'||____________________'
        WRITE(17,*)''
        
        DO WHILE (ERRO.GT.EP_TOL) ! Iterations loop _____
            
            ITER_EP = ITER_EP + 1
            
            ! Aplicando cargas
            CALL REINF_SOLVER_PLAST(PASSO,ITER_EP)
            
            ! Calculando tensoes e deformações nos pontos de gauss
            CALL INTERPOLA_GRAND(1)
            
            ! Verificacao e correcao do processo iterativo
            DO I=1,N_ELEMENTOS_MEF
                DO J=1,N_INT_P_MEF(I)
                
                    ! Calculando criterio
                    AUX1 = SIGMA_INTP(I,J) + D_SIGMA_INTP(I,J)
                    CRIT(I,J) = DABS(AUX1) - (MP_EP_MEF(I,1)+MP_EP_MEF(I,2)*ALFA(I,J))
                    IF (CRIT(I,J) .LT. 0.0D0) CRIT(I,J) = 0.0D0
                    
                    ! Evolucao das deformações plasticas
                    LAMBDA(i,j) = CRIT(i,j)/(MP_MEF(I,1)+MP_EP_MEF(I,2))
                    D_ALFA(i,j) = LAMBDA(i,j)
                    ALFA(i,j) = ALFA(i,j) + D_ALFA(i,j)
                    IF(AUX1 .GT. 0.0D0) THEN
                        D_DEF_PL_INTP(I,J) = LAMBDA(I,J)
                    ELSE
                        D_DEF_PL_INTP(I,J) = -LAMBDA(I,J)
                    ENDIF
                    
                    ! Atualizando tensao normal
                    D_SIGMA_INTP(I,J) = D_SIGMA_INTP(I,J) - MP_MEF(I,1)*D_DEF_PL_INTP(I,J)                    
                    
                ENDDO
            ENDDO
            
            ! Processo para os CONECTION ELEMENTS
            DO I=1,N_INTERFACE_MEF
                DO J=1,2
                    !CALCULANDO CRITERO
                    AUX1 = SIGMA_INTP_CONECT(I,J) + D_SIGMA_INTP_CONECT(I,J)
                    
                    DO k=1,N_ELEMENTOS_MEF !encontrando elemento que chega em no_i e no_j
                        IF((NI_MEF(I,1) .EQ. CONECTI_MEF(K,1)) .OR. (NI_MEF(I,1) .EQ. CONECTI_MEF(K,ORDEM_MEF(K)+1))) II=K
                        IF((NI_MEF(I,2) .EQ. CONECTI_MEF(K,1)) .OR. (NI_MEF(I,2) .EQ. CONECTI_MEF(K,ORDEM_MEF(K)+1))) JJ=K
                    ENDDO
                    AUX2 = (MP_EP_MEF(II,1)+MP_EP_MEF(JJ,1))/2.0D0
                    AUX3 = (MP_EP_MEF(II,2)+MP_EP_MEF(JJ,2))/2.0D0
                    
                    CRIT(I,J) = DABS(AUX1) - (AUX2+AUX3*ALFA_CONECT(I,J))
                    IF (CRIT(I,J) .LT. 0.0D0) CRIT(I,J) = 0.0D0
                    
                    AUX2 = (MP_MEF(II,1)+MP_MEF(JJ,1))/2.0D0
                    
                    ! EVOLUCAO DAS DEFORMACOES PLASTICAS
                    LAMBDA_CONECT(i,j) = CRIT(i,j)/(AUX2+AUX3)
!                    D_ALFA(i,j) = LAMBDA(i,j)
                    ALFA_CONECT(i,j) = ALFA_CONECT(i,j) + LAMBDA_CONECT(i,j)
                    
                    IF(AUX1 .GT. 0.0D0) THEN
                        D_DEF_PL_INTP_CONECT(I,J) = LAMBDA_CONECT(I,J)
                    ELSE
                        D_DEF_PL_INTP_CONECT(I,J) = -LAMBDA_CONECT(I,J)
                    ENDIF
                                        
                    !ATUALIZANDO TENSAO NORMAL
                    D_SIGMA_INTP_CONECT(I,J) = D_SIGMA_INTP_CONECT(I,J) - AUX2*D_DEF_PL_INTP_CONECT(I,J)     
                ENDDO
            ENDDO
                        
            ! Passando tensao dos ptos de integracao para forças nodais
            CALL INTERPOLA_GRAND(2)
            
            ! Acumulando grandezas
            SIGMA_INTP = SIGMA_INTP + D_SIGMA_INTP
            DEF_INTP = DEF_INTP + D_DEF_INTP
            SIGMA_INTP_CONECT = SIGMA_INTP_CONECT + D_SIGMA_INTP_CONECT
            
            ! Calculando vetor de forças desequilibradas
            F_DES_ANT(:) = F_DES(:)
            CALL CALCULA_F_DES
            
            ! Encontrando erro
            CALL CALCULA_ERRO(ERRO)
            
            ! Escreve resultados no txt
            WRITE(17,*)''            
            WRITE(17,'(a,i3,a)')'---------------ITERACAO:',ITER_EP,'-------------------'
            WRITE(17,*)''
            WRITE(17,*)'ERRO NA ITERACAO:',ERRO
            WRITE(17,*)''
            WRITE(17,*)'GRANDEZAS ATUAIS DA ITERACAO NOS NÓS DOS ELEMENTOS:'
            WRITE(17,'(A)')'NO_TOTAL   ELEM  NO_LOCAL    NORMAL             DEF    '
            K=0
            DO I=1,N_ELEMENTOS_MEF
                DO J=1,ORDEM_MEF(I)+1
                    K=K+1
                    WRITE(17,'(I8,I6,I6,E18.6,E18.6)')K,I,J,GRAND_ACUM_NOD(K,1),GRAND_ACUM_NOD(K,2)
                ENDDO
            ENDDO
            WRITE(17,*)''
            
            
        ENDDO ! End of iterations loop __________________
        
        IF (NODE_HIST.NE.'N') CALL NODE_HISTORIC_PRINT(PASSO)
         
    ENDDO ! End of load step loop _______________________
    
    
    ! SAVING MATRIXES AND VECTORS TO CONTINUE THE CODE
    ALLOCATE(PL_STRAIN_NOD(N_NOS_MEF),contelem_pnod(N_NOS_MEF))
    E_NORMAL_NOD=0.0D0
    contelem_pnod=0
    PL_STRAIN_NOD=0.0D0
    
    K=0
    DO I=1,N_ELEMENTOS_MEF
        DO J=1,ORDEM_MEF(I)+1
            K=K+1
            E_NORMAL(I,J)=GRAND_ACUM_NOD(K,1)
            
            E_NORMAL_NOD(CONECTI_MEF(I,J))=E_NORMAL_NOD(CONECTI_MEF(I,J))+E_NORMAL(I,J)
            PL_STRAIN_NOD(CONECTI_MEF(I,J))=PL_STRAIN_NOD(CONECTI_MEF(I,J))+GRAND_ACUM_NOD(K,2)
            contelem_pnod(CONECTI_MEF(I,J))=contelem_pnod(CONECTI_MEF(I,J))+1       
        ENDDO
    ENDDO
    
    DO I=1,N_NOS_MEF
        E_NORMAL_NOD(I)=E_NORMAL_NOD(I)/DFLOAT(contelem_pnod(I))
        PL_STRAIN_NOD(I)=PL_STRAIN_NOD(I)/DFLOAT(contelem_pnod(I))
    ENDDO
    DEALLOCATE(contelem_pnod)
    
    U_MEF(:) = U_MEF_ACUM(:)
    P_MEF(:) = P_MEF_ACUM(:)
    U(:) = U_ACUM_PL(:)
    T(:) = P_ACUM_PL(:)    
    
    CLOSE (17)
        
    
END SUBROUTINE REINFORCEMENTS_EP_PROCESS
        
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

SUBROUTINE INTERPOLA_GRAND(FLAG)  
    
    !
    ! SUBROUTINE TO INTERPOLATE STRESS AND STRAIN VALUES FOR THE ELASTOPLASTIC ANALYSIS 
    ! IF FLAG = 1: FIND VALUES IN THE INTEGRATION POINTS, USING NODAL VALUES 
    ! IF FLAG = 2: FIND VALUES IN THE NODAL POINTS(INCLUDING P_MEF_AT), USING THE VALUES IN INTEGRATION POINTS
    !

    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER,INTENT(IN)::FLAG
    
    INTEGER::I,J,K,II,JJ,KK,KK_G,CONT1,CONT2,N_NOS,P_INT_MEF,IPIV[ALLOCATABLE](:),INFO
    
    REAL*8::AUX1,AUX2,AUX3,AUX4,X,Y,Z,TANG(4),QSI_SOURCE,QSI_ELEMENT,JACOBIAN_SIDE,E1
    
    REAL*8,DIMENSION(:),ALLOCATABLE::SOLUTION,PHI_MEF,E_NORMAL_L,D_PHI_MEF,FFLOCAL,SOLUCAO3,SIGMA_L,FF,&
        GI_MEF,OME_MEF,SOLUCAO2,SOLUCAO
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::TEN_NOD_L,QSIW_MEF,VALUES,MR,MAT,MAT_A,MAT_B
    
    
    
    IF (FLAG .EQ. 1) THEN ! ______________ FLAG 1: FIND INTEGRATION POINTS VALUES ________________
        
        ALLOCATE(TEN_NOD_L(N_ELEMENTOS_MEF,MAXORDER_MEF+1))
                
        D_SIGMA_INTP=0.0D0
        D_DEF_INTP = 0.0D0
        
        DO I=1,N_ELEMENTOS_MEF
            
            ALLOCATE(VALUES(ORDEM_MEF(I)+1,3),PHI_MEF(ORDEM_MEF(I)+1),QSIW_MEF(N_INT_P_MEF(I),2))
            
            DO J=1,ORDEM_MEF(I)+1
                ! Encontrando valores nodais locais de tensao
                TEN_NOD_L(I,J)=E_NORMAL(I,J)/MP_MEF(I,2)
                
!                VALORES(JJ,1)=COORDNOS_MEF(CONECTI_MEF(I,JJ),1)
!                VALORES(JJ,2)=COORDNOS_MEF(CONECTI_MEF(I,JJ),2) !DUVIDA
                VALUES(J,1)=COORDPCOLOC_MEF(CONECTI_MEF(I,J),1)
                VALUES(J,2)=COORDPCOLOC_MEF(CONECTI_MEF(I,J),2)
                VALUES(J,3)=COORDPCOLOC_MEF(CONECTI_MEF(I,J),3)
            ENDDO
            
            ! Pontos de integracao
            CALL GAUSS_POINTS(N_INT_P_MEF(I),QSIW_MEF)
                
            DO J=1,N_INT_P_MEF(I)
                
                !OBS::: VER SE O ESFORÇO NORMAL TA SENDO CALCULADO NO PTO COLOCACAO OU NO GEOMETRICO DA FIBRA
                
                !CALL FUNCOES_DE_FORMA_MEC_MEF(2,QSIW_MEF(J,1),(ORDEM_MEF(I)+1),I,VALORES,X,Y,PHI_MEF)
                CALL FUNCOES_DE_FORMA_MEC_MEF(2,QSIW_MEF(J,1),(ORDEM_MEF(I)+1),I,VALUES,X,Y,Z,PHI_MEF)
                                
                AUX1=0.0D0
                DO JJ=1,ORDEM_MEF(I)+1
                    AUX1 = AUX1 + PHI_MEF(JJ)*TEN_NOD_L(I,JJ)
                ENDDO
                
                D_SIGMA_INTP(I,J) = AUX1
                D_DEF_INTP(I,J) = D_SIGMA_INTP(I,J)/MP_MEF(I,1)  
                
            ENDDO
            
            DEALLOCATE(VALUES,QSIW_MEF,PHI_MEF)
            
        ENDDO
             
        ! _________________________________ PARTE DO JOINT ELEMENT ____________________________________
        ALLOCATE(QSIW_MEF(2,2))
        CALL GAUSS_POINTS(2,QSIW_MEF)   
        
        DO I=1,N_INTERFACE_MEF
                         
            !encontrando elemento que chega em no_i e no_j
            DO K=1,N_ELEMENTOS_MEF 
                IF(NI_MEF(I,1) .EQ. CONECTI_MEF(K,1)) THEN
                    II=K
                    KK=1
                ELSE IF (NI_MEF(I,1) .EQ. CONECTI_MEF(K,ORDEM_MEF(K)+1)) THEN
                    II=K
                    KK=ORDEM_MEF(K)+1
                ENDIF
                IF(NI_MEF(I,2) .EQ. CONECTI_MEF(K,1)) THEN
                    JJ=K
                    KK_G=1
                ELSE IF (NI_MEF(I,2) .EQ. CONECTI_MEF(K,ORDEM_MEF(K)+1)) THEN
                    JJ=K
                    KK_G=ORDEM_MEF(K)+1
                ENDIF
            ENDDO     
            
            ! CALCULANDO TENSAO E DEFORMACAO       
            CONT1=0
            CONT2=0           
            DO J=1,II-1
                DO K=1,ORDEM_MEF(J)+1
                    CONT1=CONT1+1
                ENDDO
            ENDDO
            DO J=1,JJ-1
                DO K=1,ORDEM_MEF(J)+1
                    CONT2=CONT2+1
                ENDDO
            ENDDO
            
            DO J=1,2
                !AUX1 = TEN_NOD(CONT1+KK)
                !AUX2 = TEN_NOD(CONT2+KK_G)
                AUX1= TEN_NOD_L(II,KK)
                AUX2= TEN_NOD_L(JJ,KK_G)
                
                D_SIGMA_INTP_CONECT(I,J) = AUX1 + 0.50D0*(QSIW_MEF(J,1)+1.0D0)*(AUX2-AUX1)
                
                AUX1 = AUX1/MP_MEF(II,1)
                AUX2 = AUX2/MP_MEF(JJ,1)
                D_DEF_INTP_CONECT(I,J) =  AUX1 + 0.50D0*(QSIW_MEF(J,1)+1.0D0)*(AUX2-AUX1)           
            ENDDO                                         
        ENDDO 
        ! ___________________________ FIM DA PARTE DO JOINT ELEMENT _________________________________
        
        DEALLOCATE(QSIW_MEF)
        DEALLOCATE(TEN_NOD_L)
                    
    ELSE IF (FLAG .EQ. 2) THEN ! _________ FLAG 2: FIND NODES VALUES _____________________________
        
        P_MEF_AT = 0.0D0 
        SIGMA_NOD = 0.0D0
        DEF_PL_NOD = 0.0D0
        CONT1=0
        
        DO I=1,N_ELEMENTOS_MEF
            
            ! Pontos de integracao
            ALLOCATE(QSIW_MEF(N_INT_P_MEF(I),2))
            CALL GAUSS_POINTS(N_INT_P_MEF(I),QSIW_MEF)
            
            ! _________________ ENCONTRANDO ESFORÇOS NORMAIS NODAIS ______________________
            ALLOCATE(VALUES(ORDEM_MEF(I)+1,3),PHI_MEF(N_INT_P_MEF(I)),SOLUTION(N_INT_P_MEF(I)),E_NORMAL_L(2*(ORDEM_MEF(I)+1)))

            DO K=1,ORDEM_MEF(I)+1
!                VALORES(K,1)=COORDNOS_MEF(CONECTI_MEF(I,K),1)
!                VALORES(K,2)=COORDNOS_MEF(CONECTI_MEF(I,K),2)
                VALUES(K,1)=COORDPCOLOC_MEF(CONECTI_MEF(I,K),1)
                VALUES(K,2)=COORDPCOLOC_MEF(CONECTI_MEF(I,K),2)
                VALUES(K,3)=COORDPCOLOC_MEF(CONECTI_MEF(I,K),3)                
            ENDDO
                  
            E_NORMAL_L=0.0D0
            DO J=1,ORDEM_MEF(I)+1 
                
                SOLUTION(:)=QSIW_MEF(:,1)  ! Coordenadas admensionais dos pontos de integracao
                
                CALL FUNCOES_FORMA_MEF_P_INT(I,QSI_MEF(I,J),N_INT_P_MEF(I),SOLUTION,PHI_MEF)
                
                AUX1=0.0D0
                AUX2=0.0D0
                AUX3=0.0D0
                DO JJ=1,N_INT_P_MEF(I)
                    AUX1 = AUX1 + PHI_MEF(JJ)*(-MP_MEF(I,1)*D_DEF_PL_INTP(I,JJ))
                    AUX2 = AUX2 + PHI_MEF(JJ)*D_DEF_PL_INTP(I,JJ)
                    AUX3 = AUX3 + PHI_MEF(JJ)*(D_SIGMA_INTP(I,JJ))
                ENDDO
                
                SIGMA_NOD(I,J) = AUX1   ! VARIACAO DE TENSAO
                DEF_PL_NOD(I,J) = AUX2
                
                CONT1=CONT1+1
                GRAND_ACUM_NOD(CONT1,1)=GRAND_ACUM_NOD(CONT1,1)+AUX3*MP_MEF(I,2)        ! Esforço normal
                GRAND_ACUM_NOD(CONT1,2)=GRAND_ACUM_NOD(CONT1,2)+AUX2                    ! Deformacao plastica
                
                E_NORMAL_L(2*J-1)=SIGMA_NOD(I,J)*MP_MEF(I,2)
                E_NORMAL_L(2*J)=0.0D0 
                
                !E_NORMAL(CONT1,3) = E_NORMAL(CONT1,3)*E_NORMAL(CONT1,4) + SIGMA_NOD(I,J) 
                E_NORMAL(I,J) = E_NORMAL(I,J) + SIGMA_NOD(I,J) 
                
            ENDDO
            DEALLOCATE(SOLUTION,PHI_MEF)
            ! ______________ FIM DE ENCONTRAR ESFORÇOS NORMAIS NODAIS ____________________
            
            ! _____________ ENCONTRANDO FORÇAS NODAIS GLOBAIS (P_MEF_AT) _________________
            N_NOS=ORDEM_MEF(I)+1
            ALLOCATE(MAT(N_INT_P_MEF(I),N_NOS),SIGMA_L(N_INT_P_MEF(I)),FF(3*N_NOS),D_PHI_MEF(N_NOS),&
                FFLOCAL(N_NOS),MR(N_NOS,3*N_NOS),PHI_MEF(N_NOS))
            SIGMA_L=0.0D0
            MR=0.0D0
            MAT=0.0D0
            
            DO J=1,N_NOS       ! CALCULATING MR
                CALL FUNCOES_DE_FORMA_MEC_MEF(3,QSI_MEF(I,J),N_NOS,I,VALUES,X,Y,Z,D_PHI_MEF)
                 
                AUX1=0.0D0
                AUX2=0.0D0
                AUX3=0.0D0
                DO K=1,N_NOS
		            AUX1=AUX1+D_PHI_MEF(K)*COORDPCOLOC_MEF(CONECTI_MEF(I,K),1)
		            AUX2=AUX2+D_PHI_MEF(K)*COORDPCOLOC_MEF(CONECTI_MEF(I,K),2)
                    AUX3=AUX3+D_PHI_MEF(K)*COORDPCOLOC_MEF(CONECTI_MEF(I,K),3)
                ENDDO
                AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
                AUX1=AUX1/AUX4
                AUX2=AUX2/AUX4
                AUX3=AUX3/AUX4
            
                !TO AVOID NUMERICAL ERRORS:            
                IF (ABS(AUX1).LT.TOLER) AUX1=0.D0
                IF (ABS(AUX2).LT.TOLER) AUX2=0.D0
                IF (ABS(AUX3).LT.TOLER) AUX3=0.D0

                MR(J,3*J-2) = AUX1
                MR(J,3*J-1) = AUX2
                MR(J,3*J) = AUX3
            ENDDO
                
            ! GUARDANDO INFORMACOES DOS PONTS DE INTEGRACAO
            ALLOCATE(GI_MEF(N_INT_P_MEF(I)),OME_MEF(N_INT_P_MEF(I)))
            GI_MEF(:)=QSIW_MEF(:,1)
            OME_MEF(:)=QSIW_MEF(:,2)
            DEALLOCATE(QSIW_MEF)
            
            ! INTEGRACAO INTERNA
            P_INT_MEF=CEILING((ORDEM_MEF(I)+1)/2.0D0) 
            ALLOCATE(QSIW_MEF(P_INT_MEF,2))
            CALL GAUSS_POINTS(P_INT_MEF,QSIW_MEF)
            
            DO K=1,N_INT_P_MEF(I)
                
                ! MONTANDO VETOR A DIREITA (SOLUCAO1)
                SIGMA_L(K) = E_NORMAL_L(1) - (-MP_MEF(I,1)*D_DEF_PL_INTP(I,K))*MP_MEF(I,2)
                
                QSI_SOURCE=GI_MEF(K)
                JACOBIAN_SIDE=(qsi_source+1.0d0)/2.0d0
                
                DO II=1,N_NOS
                    
                    AUX1=0.0D0
                    
                    DO J=1,P_INT_MEF
                        QSI_ELEMENT=(QSIW_MEF(J,1)+1.0D0)*JACOBIAN_SIDE-1.0d0
                        
                        CALL FUNCOES_DE_FORMA_MEC_MEF(2,QSI_ELEMENT,N_NOS,I,VALUES,X,Y,Z,PHI_MEF)
                        CALL FUNCOES_DE_FORMA_MEC_MEF(3,QSI_ELEMENT,N_NOS,I,VALUES,X,Y,Z,D_PHI_MEF)
                        
                        TANG=0.D0
                        DO JJ=1,ORDEM_MEF(I)+1
!	                            TAN1=TAN1+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(I,JJ),1) 
!		                        TAN2=TAN2+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(I,JJ),2) 
	                        TANG(1)=TANG(1)+D_PHI_MEF(JJ)*COORDPCOLOC_MEF(CONECTI_MEF(I,JJ),1) !TANGENT VECTORS
		                    TANG(2)=TANG(2)+D_PHI_MEF(JJ)*COORDPCOLOC_MEF(CONECTI_MEF(I,JJ),2) !TANGENT VECTORS
                            TANG(3)=TANG(3)+D_PHI_MEF(JJ)*COORDPCOLOC_MEF(CONECTI_MEF(I,JJ),3) !TANGENT VECTORS
        	            ENDDO
                        TANG(4)=DSQRT(TANG(1)**2.0D0+TANG(2)**2.0D0+TANG(3)**2.0D0) !TANGENT VECTOR NORM, JACOBIAN OF THE POINT
                        
                        AUX1 = AUX1 + PHI_MEF(II)*JACOBIAN_SIDE*QSIW_MEF(J,2)*TANG(4)
                    ENDDO
                    
                    MAT(K,II) = AUX1
                    
                ENDDO ! END OF N_NOS LOOP
            ENDDO ! END OF N_INT_P_MEF(I) LOOP
            
            ! RESOLVENDO SISTEMA PARA ENCONTRAR FFLOCAL: 
            IF (N_INT_P_MEF(I) .NE. N_NOS) THEN
                II=N_INT_P_MEF(I) 
                JJ=ORDEM_MEF(I)+1                                              
                ALLOCATE(MAT_A(JJ,JJ),IPIV(JJ))
                FFLOCAL=0.0D0 
                CALL DMURRV(JJ,II,TRANSPOSE(MAT),JJ,II,SIGMA_L,1,JJ,FFLOCAL)
                CALL DMRRRR(JJ,II,TRANSPOSE(MAT),JJ,II,JJ,MAT,II,JJ,JJ,MAT_A,JJ)
                !CALL DLSLRG(JJ,MAT_A,JJ,SOLUCAO2,1,FFLOCAL)
                CALL DGESV(JJ,1,MAT_A,JJ,IPIV,FFLOCAL,JJ,INFO)
                DEALLOCATE(MAT_A,IPIV)  
            ELSE
                II=N_INT_P_MEF(I)
                ALLOCATE(IPIV(II))
                !CALL DLSLRG(II,MAT,II,SIGMA_L,1,FFLOCAL)
                FFLOCAL=SIGMA_L
                CALL DGESV(II,1,MAT,II,IPIV,FFLOCAL,II,INFO)
                DEALLOCATE(IPIV)
            ENDIF
            
            ! GIRANDO FFLOCAL NAS COORDENADS GLOBAIS (FF):
            FF=MATMUL(TRANSPOSE(MR),FFLOCAL) 
            
            DO J=1,N_NOS
                K=CONECTI_MEF(I,J)
                IF ((DABS(P_MEF_AT(3*K-2)).LE.TOLER).AND.(DABS(P_MEF_AT(3*K-1)).LE.TOLER).AND.(DABS(P_MEF_AT(3*K)).LE.TOLER)) THEN
                    P_MEF_AT(3*K-2) = -FF(3*J-2)
                    P_MEF_AT(3*K-1) = -FF(3*J-1)
                    P_MEF_AT(3*K) = -FF(3*J)
                ELSE
                    P_MEF_AT(3*K-2) = (P_MEF_AT(3*K-2) - FF(3*J-2))/2.0D0
                    P_MEF_AT(3*K-1) = (P_MEF_AT(3*K-1) - FF(3*J-1))/2.0D0
                    P_MEF_AT(3*K) = (P_MEF_AT(3*K) - FF(3*J))/2.0D0   
!                    P_MEF_AT(2*K-1) = (P_MEF_AT(2*K-1) - FF(2*J-1))/1.0D0
!                    P_MEF_AT(2*K) = (P_MEF_AT(2*K) - FF(2*J))/1.0D0 
                ENDIF
            ENDDO 
            
            DEALLOCATE(VALUES,QSIW_MEF,GI_MEF,OME_MEF,E_NORMAL_L,MR,SIGMA_L,FF,MAT,FFLOCAL,D_PHI_MEF,PHI_MEF)            
            
        ENDDO ! END OF FE LOOP
        ! ___________ FIM DE ENCONTRAR FORÇAS NODAIS GLOBAIS (P_MEF_AT) _______________

        ! _________________________ PARTE DO JOINT ELEMENT _______________________
        ALLOCATE(QSIW_MEF(2,2))
        CALL GAUSS_POINTS(2,QSIW_MEF)
	    
	    ALLOCATE(PHI_MEF(2),D_PHI_MEF(2),FFLOCAL(2),SOLUCAO3(6),MR(2,6),FF(6),SOLUCAO(2))
        
        DO I=1,N_INTERFACE_MEF
        
            ! ENCONTRANDO ELEMENTO QUE CHEGA EM no_i e no_j
            DO K=1,N_ELEMENTOS_MEF 
                IF(NI_MEF(I,1) .EQ. CONECTI_MEF(K,1)) THEN
                    CONT1=K
                    KK=1
                ELSE IF (NI_MEF(I,1) .EQ. CONECTI_MEF(K,ORDEM_MEF(K)+1)) THEN
                    CONT1=K
                    KK=ORDEM_MEF(K)+1
                ENDIF
                IF(NI_MEF(I,2) .EQ. CONECTI_MEF(K,1)) THEN
                    CONT2=K
                    KK_G=1
                ELSE IF (NI_MEF(I,2) .EQ. CONECTI_MEF(K,ORDEM_MEF(K)+1)) THEN
                    CONT2=K
                    KK_G=ORDEM_MEF(K)+1
                ENDIF
            ENDDO  
            
            SOLUCAO(1)=-1.0D0
            SOLUCAO(2)=1.0D0
            
            MR=0.0D0
            DO J=1,2       ! CALCULATING MR
                CALL D_FUNCOES_FORMA_MEF_P_INT(I,SOLUCAO(J),2,SOLUCAO,D_PHI_MEF)
                 
                AUX1=0.0D0
                AUX2=0.0D0
                AUX3=0.0D0
                DO K=1,2
		            AUX1=AUX1+D_PHI_MEF(K)*COORDPCOLOC_MEF(CONECTI_MEF(I,K),1)
		            AUX2=AUX2+D_PHI_MEF(K)*COORDPCOLOC_MEF(CONECTI_MEF(I,K),2)
                    AUX3=AUX3+D_PHI_MEF(K)*COORDPCOLOC_MEF(CONECTI_MEF(I,K),3)
                ENDDO
                AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
                AUX1=AUX1/AUX4
                AUX2=AUX2/AUX4
                AUX3=AUX3/AUX4
            
                !TO AVOID NUMERICAL ERRORS:            
                IF (ABS(AUX1).LT.TOLER) AUX1=0.D0
                IF (ABS(AUX2).LT.TOLER) AUX2=0.D0
                IF (ABS(AUX3).LT.TOLER) AUX3=0.D0

                MR(J,3*J-2) = AUX1
                MR(J,3*J-1) = AUX2
                MR(J,3*J) = AUX3
            ENDDO
            
            !ENCONTRANDO ESFORÇO NORMAL NA PRIMEIRA EXTREMIDADE (E_NORMAL_L(1) = E1)
            AUX1=(MP_MEF(CONT1,1)+MP_MEF(CONT2,1))/2.0D0
            AUX4=(MP_MEF(CONT1,2)+MP_MEF(CONT2,2))/2.0D0
            
            AUX2= -AUX1*D_DEF_PL_INTP_CONECT(I,1)*AUX4
            AUX3= -AUX1*D_DEF_PL_INTP_CONECT(I,2)*AUX4
            E1 = AUX2 + ((AUX3-AUX2)/(QSIW_MEF(2,1)-QSIW_MEF(1,1)))*(-1.0D0-QSIW_MEF(1,1))
            
            ALLOCATE(SIGMA_L(2),MAT(2,2),IPIV(2))
            
            DO K=1,2
                                
                ! MONTANDO VETOR A DIREITA (SOLUCAO1)
                SIGMA_L(K) = E1 - (-AUX1*D_DEF_PL_INTP_CONECT(I,K))*AUX4
                               
                QSI_SOURCE=QSIW_MEF(K,1)
                JACOBIAN_SIDE=(QSI_SOURCE+1.0d0)/2.0d0
                
                DO II=1,2
                    
                    AUX2=0.0D0
                    DO J=1,2
                        QSI_ELEMENT=(QSIW_MEF(J,1)+1.0D0)*JACOBIAN_SIDE-1.0d0

                        CALL FUNCOES_FORMA_MEF_P_INT(I,QSI_ELEMENT,2,SOLUCAO,PHI_MEF)
                        CALL D_FUNCOES_FORMA_MEF_P_INT(I,QSI_ELEMENT,2,SOLUCAO,D_PHI_MEF)
                            
                        TANG=0.0D0
                        DO JJ=1,3       ! TANGENT VECTORS
                            TANG(JJ)=TANG(JJ)+D_PHI_MEF(1)*COORDPCOLOC_MEF(CONECTI_MEF(CONT1,KK),JJ)
                            TANG(JJ)=TANG(JJ)+D_PHI_MEF(2)*COORDPCOLOC_MEF(CONECTI_MEF(CONT2,KK_G),JJ)
                        ENDDO      	                
        	            TANG(4)=DSQRT(TANG(1)**2.0D0+TANG(2)**2.0D0+TANG(3)**2.0D0) !TANGENT VECTOR NORM, JACOBIAN OF THE POINT
                        
                        AUX2 = AUX2 + PHI_MEF(II)*JACOBIAN_SIDE*QSIW_MEF(J,2)*TANG(4)
                    ENDDO
                    
                    MAT(K,II) = AUX2
                    
                ENDDO
            ENDDO
            
            II=2
            FFLOCAL=SIGMA_L
            !CALL DGESV(II,1,MAT,II,IPIV,FFLOCAL,II,INFO)
            CALL GESV(MAT,FFLOCAL,IPIV=IPIV,INFO=INFO)
            
            ! GIRANDO FFLOCAL NAS COORDENADS GLOBAIS (FF):
            FF=MATMUL(TRANSPOSE(MR),FFLOCAL)   
            
            DEALLOCATE(MAT,SIGMA_L,IPIV)
            
            DO J=1,2
!                K=CONECTI_MEF(I,J)
                IF (J.EQ.1) K=CONECTI_MEF(CONT1,KK)
                IF (J.EQ.2) K=CONECTI_MEF(CONT2,KK_G)

                IF ((DABS(P_MEF_AT(3*K-2)).LE.TOLER).AND.(DABS(P_MEF_AT(3*K-1)).LE.TOLER).AND.(DABS(P_MEF_AT(3*K)).LE.TOLER)) THEN
                    P_MEF_AT(3*K-2) = -FF(3*J-2)
                    P_MEF_AT(3*K-1) = -FF(3*J-1)
                    P_MEF_AT(3*K) = -FF(3*J)
                ELSE
                    P_MEF_AT(3*K-2) = (P_MEF_AT(3*K-2) - FF(3*J-2))/2.0D0
                    P_MEF_AT(3*K-1) = (P_MEF_AT(3*K-1) - FF(3*J-1))/2.0D0
                    P_MEF_AT(3*K) = (P_MEF_AT(3*K) - FF(3*J))/2.0D0   
!                    P_MEF_AT(2*K-1) = (P_MEF_AT(2*K-1) - FF(2*J-1))/1.0D0
!                    P_MEF_AT(2*K) = (P_MEF_AT(2*K) - FF(2*J))/1.0D0 
                ENDIF
            ENDDO 
            
        ENDDO
        
        DEALLOCATE(QSIW_MEF,PHI_MEF,D_PHI_MEF,FFLOCAL,SOLUCAO3,SOLUCAO,MR,FF)           
        ! _______________________ FIM DA PARTE DO JOINT ELEMENT __________________
        
        
    ELSE ! _____________________________ NO - FLAG _____________________________
        
        WRITE(*,*)'******* ERROR IN THE INTERPOLA_GRAND SUBROUTINE ******'
        READ(*,*)
        
    ENDIF
    
END SUBROUTINE INTERPOLA_GRAND
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

SUBROUTINE CALCULA_ERRO(ERR)

    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    REAL*8,INTENT(OUT)::ERR

    INTEGER::I,J,K,II,JJ
    
    REAL*8::AUX1,AUX2,SOLUCAO[ALLOCATABLE](:)
    
    II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
    JJ=3*N_COLLOCPOINTS+6*N_NOS_MEF
    
    ALLOCATE(SOLUCAO(JJ))
    SOLUCAO(:) = F_DES(:) - F_DES_ANT(:)
    
    AUX1=0.0D0
    AUX2=0.0D0
    DO I=1,JJ
        AUX1=AUX1+(F_ALVO(I))**2.0D0
        AUX2=AUX2+(SOLUCAO(I))**2.0D0
    ENDDO
    
    ERR = DSQRT(AUX2)/DSQRT(AUX1)
    
    DEALLOCATE(SOLUCAO)

END SUBROUTINE CALCULA_ERRO
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

SUBROUTINE CALCULA_F_DES

    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,J,K,II,JJ,CONT
        
    REAL*8,DIMENSION(:),ALLOCATABLE::SOLUCAO,F
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::AGLOBAL
    
    !STARTING VARIABLES
    II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
	JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
    ALLOCATE(SOLUCAO(JJ),F(JJ))
    
    ! CREATING SOLUCAO
    DO I=1,N_COLLOCPOINTS
	    II=0
	    JJ=0
	    IF(DUAL_BEM(I).EQ.'B')THEN  ! POINT AT BOUNDARY
		    DO J=1,N_INTERFACE_POINTS
			    IF(I.EQ.NI(J,1)) THEN
				    II=-1
				    JJ=NI(J,2)
			    ELSE
				    IF(I.EQ.NI(J,2)) THEN
					    II=NI(J,1)
					    JJ=-1
				    ENDIF
			    ENDIF
		    ENDDO	
    
		    IF(II.EQ.-1) THEN
			    SOLUCAO(3*I-2)=U(3*I-2)
			    SOLUCAO(3*I-1)=U(3*I-1)
			    SOLUCAO(3*I)=U(3*I)
    
			    SOLUCAO(3*I-2)=U(3*JJ-2)
			    SOLUCAO(3*I-1)=U(3*JJ-1)
			    SOLUCAO(3*I)=U(3*JJ)
		    ELSE
			    IF(JJ.EQ.-1) THEN
				    SOLUCAO(3*I-2)=-T(3*I-2)
				    SOLUCAO(3*I-1)=-T(3*I-1)
			        SOLUCAO(3*I)=-T(3*I)
    
				    SOLUCAO(3*I-2)=T(3*II-2)
				    SOLUCAO(3*I-1)=T(3*II-1)
				    SOLUCAO(3*I)=T(3*II)
			    ELSE
				    IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					    SOLUCAO(3*I-2)=T(3*I-2)
				    ELSE
					    SOLUCAO(3*I-2)=U(3*I-2)
				    ENDIF
    
				    IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					    SOLUCAO(3*I-1)=T(3*I-1)
				    ELSE
					    SOLUCAO(3*I-1)=U(3*I-1)
				    ENDIF
    
				    IF (B_CONDITIONS(3*I).EQ.0) THEN
					    SOLUCAO(3*I)=T(3*I)
				    ELSE
					    SOLUCAO(3*I)=U(3*I)
				    ENDIF
			    ENDIF
		    ENDIF
        ELSE   ! POINT AT CRACK CONTOURN
            !
            ! UNDER CONSTRUCTION
            !
        ENDIF
    ENDDO
    II=3*N_COLLOCPOINTS
    JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
    DO I=1,(3*N_NOS_MEF)
        SOLUCAO(II+I)=U_MEF(I)
        SOLUCAO(JJ+I)=P_MEF(I) + P_MEF_AT(I)
    ENDDO 
    
    P_MEF_ACUM(:) = P_MEF_ACUM(:) + P_MEF(:) - P_MEF_AT(:)
    
    !II=3*N_COLLOCPOINTS+3*N_NOS_MEF
    !JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
    !DO J=1,JJ
    !    DO I=1,3*N_NOS_MEF
    !        AGLOBAL(I+II,J)=0.0D0
    !    ENDDO
    !ENDDO
    
    JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
    CALL DGEMV('N',JJ,JJ,1.0D0,A_GLOBAL_ELASTIC,JJ,SOLUCAO,1,0.0D0,F,1)
    
    II=3*N_COLLOCPOINTS+3*N_NOS_MEF
    DO I=1,3*N_NOS_MEF
        F(I+II)=0.0D0
    ENDDO
    
    F_ACUM(:) = F_ACUM(:) + F(:)
       
	F_DES(:) =  F_ALVO(:) - F_ACUM(:)
    
    DEALLOCATE(F,SOLUCAO)    
    

END SUBROUTINE CALCULA_F_DES