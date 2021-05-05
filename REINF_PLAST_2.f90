SUBROUTINE REINFORCEMENTS_EP_PROCESS_2
        
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
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::CRIT,D_LAMBDA,D_ALFA
    
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
    F_DES_ANT(3*N_NOS_MEF),F_ALVO(JJ),F_DES(JJ),F_ACUM(JJ))  
    
    ALLOCATE(ALFA_CONECT(N_INTERFACE_MEF,2),SIGMA_INTP_CONECT(N_INTERFACE_MEF,2),D_DEF_PL_INTP_CONECT(N_INTERFACE_MEF,2),LAMBDA_CONECT(N_INTERFACE_MEF,2),&
    D_SIGMA_INTP_CONECT(N_INTERFACE_MEF,2),D_DEF_INTP_CONECT(N_INTERFACE_MEF,2),NI_SLIP_MEF(N_NOS_MEF))      
     
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
            CALL REINF_SOLVER_PLAST_2(PASSO,ITER_EP)
            
            write(*,*)'passou reinf solver plast'
            
            ! Calculando tensoes e deformações nos pontos de gauss
            CALL INTERPOLA_GRAND_2(1)
            
            write(*,*)'passou interpola grand 1'
            
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
            ! UNDER CONSTRUCTION
            
            write(*,*)'passou verificacao dos elementos'
            
            ! Passando tensao dos ptos de integracao para forças nodais
            NI_SLIP_MEF=0
            CALL INTERPOLA_GRAND_2(2)
            
            write(*,*)'passou interpola grand 2'
            
            ! Acumulando grandezas
            SIGMA_INTP = SIGMA_INTP + D_SIGMA_INTP
            DEF_INTP = DEF_INTP + D_DEF_INTP
            !SIGMA_INTP_FAKE = SIGMA_INTP_FAKE + D_SIGMA_INTP_FAKE
            
            WRITE(*,*)'PASSOU ACUMULA GRANDEZAS'
            
            ! Calculando vetor de forças desequilibradas
            !F_DES_ANT(:) = F_DES(:)
            CALL CALCULA_F_DES_2(F_DES)
            
            WRITE(*,*)'PASSOU CALCULA F DES'
            
            ! Encontrando erro
            CALL CALCULA_ERRO_2(ERRO)
            
            write(*,*)'passou calcula erro'
            
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
            
            WRITE(*,'(A,I5,A,F15.7)')'ITERACAO:',ITER_EP,'  ERRO:',ERRO
            
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
        
    
END SUBROUTINE REINFORCEMENTS_EP_PROCESS_2
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

SUBROUTINE INTERPOLA_GRAND_2(FLAG)  
    
    !
    ! SUBROUTINE TO INTERPOLATE STRESS AND STRAIN VALUES FOR THE ELASTOPLASTIC ANALYSIS 
    ! IF FLAG = 1: FIND VALUES IN THE INTEGRATION POINTS, USING NODAL VALUES 
    ! IF FLAG = 2: FIND VALUES IN THE NODAL POINTS(INCLUDING P_MEF_AT), USING THE VALUES IN INTEGRATION POINTS
    !

    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER,INTENT(IN)::FLAG
    
    INTEGER::I,J,K,II,JJ,KK,CONT1,N_NOS,P_INT_MEF,INFO
    
    INTEGER,DIMENSION(:),ALLOCATABLE::IPIV,contelem_pnod,c2
    
    REAL*8::AUX1,AUX2,AUX3,AUX4,X,Y,Z,TANG(4),QSI_SOURCE,QSI_ELEMENT,JACOBIAN_SIDE
    
    REAL*8,DIMENSION(:),ALLOCATABLE::SOLUTION,PHI_MEF,E_NORMAL_L,D_PHI_MEF,FFLOCAL,SOLUCAO3,SIGMA_L,FF,&
        GI_MEF,OME_MEF,SOLUCAO2
    
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
             
        ! PARTE DO JOINT ELEMENT
        ! UNDER CONSTRUCTION
        
        DEALLOCATE(TEN_NOD_L)
                    
    ELSE IF (FLAG .EQ. 2) THEN ! ______________ FLAG 2: FIND NODES VALUES ________________________
        
        ALLOCATE(contelem_pnod(3*N_NOS_MEF),c2(3*n_nos_mef))
        P_MEF_AT = 0.0D0 
        SIGMA_NOD = 0.0D0
        DEF_PL_NOD = 0.0D0
        CONT1=0
        contelem_pnod=0
        
        c2=0
        F_DES_ANT=0.0D0
        
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
                DEF_PL_NOD(I,J) = AUX2  ! VARIACAO DE DEFORMAÇÃO PLASTICA
                
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
            
            ! mudou
            N_NOS=ORDEM_MEF(I)+1
            ALLOCATE(MR(N_NOS,3*N_NOS),SOLUTION(N_NOS),SOLUCAO2(3*N_NOS),D_PHI_MEF(N_NOS))
            MR=0.0D0
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
                
                SOLUTION(J)=SIGMA_NOD(I,J)*MP_MEF(I,2)
            ENDDO
            SOLUCAO2=MATMUL(TRANSPOSE(MR),SOLUTION)
            DO J=1,N_NOS
                F_DES_ANT(3*CONECTI_MEF(I,J)-2)=F_DES_ANT(3*CONECTI_MEF(I,J)-2)+SOLUCAO2(3*J-2)
                F_DES_ANT(3*CONECTI_MEF(I,J)-1)=F_DES_ANT(3*CONECTI_MEF(I,J)-1)+SOLUCAO2(3*J-1)
                F_DES_ANT(3*CONECTI_MEF(I,J))=F_DES_ANT(3*CONECTI_MEF(I,J))+SOLUCAO2(3*J)
                
                IF (DABS(SOLUCAO2(3*J-2)).GT.TOLER) C2(3*CONECTI_MEF(I,J)-2)=C2(3*CONECTI_MEF(I,J)-2)+1
                IF (DABS(SOLUCAO2(3*J-1)).GT.TOLER) C2(3*CONECTI_MEF(I,J)-1)=C2(3*CONECTI_MEF(I,J)-1)+1
                IF (DABS(SOLUCAO2(3*J)).GT.TOLER) C2(3*CONECTI_MEF(I,J))=C2(3*CONECTI_MEF(I,J))+1
            ENDDO
                        
            DEALLOCATE(MR,SOLUTION,SOLUCAO2,D_PHI_MEF)
            ! ACABOU O MUDOU
            
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
                !IF ((DABS(P_MEF_AT(3*K-2)).LE.TOLER).AND.(DABS(P_MEF_AT(3*K-1)).LE.TOLER).AND.(DABS(P_MEF_AT(3*K)).LE.TOLER)) THEN
                !    P_MEF_AT(3*K-2) = -FF(3*J-2)
                !    P_MEF_AT(3*K-1) = -FF(3*J-1)
                !    P_MEF_AT(3*K) = -FF(3*J)
                !ELSE
                !    P_MEF_AT(3*K-2) = (P_MEF_AT(3*K-2) - FF(3*J-2))/2.0D0
                !    P_MEF_AT(3*K-1) = (P_MEF_AT(3*K-1) - FF(3*J-1))/2.0D0
                !    P_MEF_AT(3*K) = (P_MEF_AT(3*K) - FF(3*J))/2.0D0   
                !ENDIF
                
                ! MUDOU
                P_MEF_AT(3*K-2)=P_MEF_AT(3*K-2)-FF(3*J-2)
                P_MEF_AT(3*K-1)=P_MEF_AT(3*K-1)-FF(3*J-1)
                P_MEF_AT(3*K)=P_MEF_AT(3*K)-FF(3*J)

                IF (DABS(FF(3*J-2)).GT.TOLER) contelem_pnod(3*K-2)=contelem_pnod(3*K-2)+1   
                IF (DABS(FF(3*J-1)).GT.TOLER) contelem_pnod(3*K-1)=contelem_pnod(3*K-1)+1   
                IF (DABS(FF(3*J)).GT.TOLER) contelem_pnod(3*K)=contelem_pnod(3*K)+1   
            ENDDO 
            
            DEALLOCATE(VALUES,QSIW_MEF,GI_MEF,OME_MEF,E_NORMAL_L,MR,SIGMA_L,FF,MAT,FFLOCAL,D_PHI_MEF,PHI_MEF)            
            
        ENDDO ! END OF FE LOOP
        
        ! FAZENDO MEDIA EM NÓS QUE RECEBERAM MAIS DE UM VALOR
        DO I=1,3*N_NOS_MEF
            IF (contelem_pnod(I).EQ.0) contelem_pnod(I)=1
            P_MEF_AT(I)=P_MEF_AT(I)/DFLOAT(contelem_pnod(I))
        ENDDO
        DEALLOCATE(contelem_pnod)
        
        ! mudou
        DO I=1,3*N_NOS_MEF
            IF (c2(I).EQ.0) c2(I)=1
            f_des_ant(I)=f_des_ant(I)/DFLOAT(c2(I))
        ENDDO
        DEALLOCATE(c2)
        
        ! VERIFICANDO NÓS QUE TEM CARGA REAPLICADA
        DO I=1,N_NOS_MEF
            AUX1=DABS(P_MEF_AT(3*J-2))
            AUX2=DABS(P_MEF_AT(3*J-1))
            AUX3=DABS(P_MEF_AT(3*J))
            IF ((AUX1.GT.TOLER).OR.(AUX2.GT.TOLER).OR.(AUX3.GT.TOLER)) THEN
                NI_PLAST_MEF(I)=1
            ENDIF
        ENDDO
        
        ! _____________________ FIM DE ENCONTRAR FORÇAS NODAIS GLOBAIS (P_MEF_AT) _________________________

        ! _________________________ PARTE DO JOINT ELEMENT ________________________________________________
        !   UNDER CONSTRUCTION
        ! _______________________ FIM DA PARTE DO JOINT ELEMENT ___________________________________________
        
        
    ELSE ! _____________________________ NO - FLAG ________________________________________________________
        
        WRITE(*,*)'******* ERROR IN THE INTERPOLA_GRAND SUBROUTINE ******'
        READ(*,*)
        
    ENDIF
    
END SUBROUTINE INTERPOLA_GRAND_2
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

SUBROUTINE CALCULA_ERRO_2(ERR)

    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    REAL*8,INTENT(OUT)::ERR

    INTEGER::I,J,JJ
    
    REAL*8::AUX1,AUX2
    
    AUX1=0.0D0
    AUX2=0.0D0
    JJ=3*N_NOS_MEF
    
    DO I=1,JJ
        AUX1=AUX1+(P_MEF_AT(I))**2.0D0
        AUX2=AUX2+(P_MEF_ACUM(I))**2.0D0
    ENDDO
    
    ERR=DSQRT(AUX1)/DSQRT(AUX2)
    
END SUBROUTINE CALCULA_ERRO_2
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

SUBROUTINE CALCULA_F_DES_2(VETOR)

    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,J,K,II,JJ,KK
        
    REAL*8,INTENT(OUT)::VETOR(3*N_COLLOCPOINTS+6*N_NOS_MEF)
    
    REAL*8,DIMENSION(:),ALLOCATABLE::SOLUTION
        
    IF (.NOT.ALLOCATED(B_GLOBAL_EP_D)) THEN ! CRIA MATRIZ PARA MULTIPLICAR PELO VETOR DESEQUILIBRADO
        
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
        ALLOCATE(B_GLOBAL_EP_D(JJ,3*N_NOS_MEF))
        
        DO J=1,3*N_NOS_MEF
            DO I=1,3*N_COLLOCPOINTS
                B_GLOBAL_EP_D(I,J)=G_CF(I,J)
            ENDDO
            II=3*N_COLLOCPOINTS
            KK=3*N_NOS_MEF
            DO I=1,3*N_NOS_MEF
                B_GLOBAL_EP_D(II+I,J)=G_FF(I,J)
                B_GLOBAL_EP_D(II+KK+I,J)=-GGLOBAL(I,J)
            ENDDO
        ENDDO
        
    ENDIF
    
    II=3*N_NOS_MEF
    JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF    
    ! MULTIPLICANDO BGLOBAL_EP_D POR P_MEF_AT (VETOR DESEQUILIBRADO)
    !CALL DGEMV('N',JJ,II,1.0D0,B_GLOBAL_EP_D,JJ,-P_MEF_AT,1,0.0D0,SOLUTION,1)  ! DUVIDA NO SINAL DE P_MEF_AT
    VETOR=MATMUL(B_GLOBAL_EP_D,P_MEF_AT)
    
    ! MUDOU
    VETOR=0.0D0
    II=3*N_COLLOCPOINTS+3*N_NOS_MEF
    DO I=1,3*N_NOS_MEF
        VETOR(II+I)=f_des_ant(I)
    ENDDO
    ! FIM DE MUDOU
    
    
    WRITE(*,*)'PASSOU MULTIPLICACAO'
        
    WRITE(*,*)'PASSOU igualar f des'  
    
        WRITE(*,*)'PASSOU dealoca solution'
    
    ! ACUMULANDO P_MEF
    !P_MEF_ACUM(:) = P_MEF_ACUM(:) + P_MEF(:) - P_MEF_AT(:)
    P_MEF_ACUM = P_MEF_ACUM + P_MEF ! MUDOU
    
    WRITE(*,*)'PASSOU ACUMULA P MEF ACUM'
    

END SUBROUTINE CALCULA_F_DES_2