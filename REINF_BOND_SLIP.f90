SUBROUTINE REINFORCEMENTS_SLIP_PROCESS
            
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE REINFORCEMENTS
    USE OMP_LIB
    
    IMPLICIT NONE
    
    LOGICAL::NI_CHANGE
    
    INTEGER::I,J,K,KK_G,PASSO,ITER_SLIP,TIPO_PASS_RED
    
    INTEGER,DIMENSION(:),ALLOCATABLE::NI_SLIP_MEF_ANT,NI_SLIP_CONVERG
    
    REAL*8,DIMENSION(:),ALLOCATABLE::DELTA_F_C,D_P_LOCAL,DELTA_F_C_GLOBAL,DELTA_S,S_AUX,CARGA_NODAL
    
    REAL*8::ERRO,F_MODELO,AUX1,AUX2,F_MOD,PASS_RAT
    
    ! Finding total number of local element's nodes (for matrixes)
    TOTAL_ELEM_NODES_MEF=0
    DO I=1,N_ELEMENTOS_MEF
       TOTAL_ELEM_NODES_MEF=TOTAL_ELEM_NODES_MEF+ORDEM_MEF(I)+1
    ENDDO
    
    !TIPO_PASS_RED=1
    !TIPO_PASS_RED=2
    
    IF (SLIP_LAW.EQ.'LINEAR') TIPO_PASS_RED=2
    IF (SLIP_LAW.EQ.'HUANG2') TIPO_PASS_RED=1
    
    ! Initiating variables
    ALLOCATE(SLIP_MEF(N_NOS_MEF),CARGA(3*N_COLLOCPOINTS,2),DELTA_F_C(TOTAL_ELEM_NODES_MEF),P_MEF_ACUM_LOCAL(TOTAL_ELEM_NODES_MEF),&
        D_P_LOCAL(TOTAL_ELEM_NODES_MEF),DELTA_F_C_GLOBAL(3*N_NOS_MEF),S_ACUM(TOTAL_ELEM_NODES_MEF),P_MEF_ACUM(3*N_NOS_MEF),&
        U_MEF_ACUM(3*N_NOS_MEF),U_ACUM_PL(3*N_COLLOCPOINTS),P_ACUM_PL(3*N_COLLOCPOINTS),NI_SLIP_MEF(N_NOS_MEF),&
        F_APP_DISP_ACUM(3*N_NODAL_FORCES),NI_SLIP_MEF_ANT(N_NOS_MEF),CARGA_NODAL(N_NODAL_FORCES),NI_SLIP_CONVERG(N_NOS_MEF))
    SLIP_MEF=0.0D0
    DELTA_F_C=0.0D0
    DELTA_F_C_GLOBAL=0.0D0
    D_P_LOCAL=0.0D0
    U_ACUM_PL=0.0D0
    P_ACUM_PL=0.0D0
    P_MEF_ACUM=0.0D0
    U_MEF_ACUM=0.0D0
    F_APP_DISP_ACUM=0.0D0
    S_ACUM=0.0D0
    P_MEF_ACUM_LOCAL=0.0D0
    NI_SLIP_MEF=0  ! TODOS OS NÓS ESTAO COLADOS (UD=UE)
    NI_SLIP_MEF_ANT=0
    NI_SLIP_CONVERG=0
    !CARGA(:,1) = U(:) / DFLOAT(N_PASSOS)
    !CARGA(:,2) = T(:) / DFLOAT(N_PASSOS)
    !NODAL_FORCES_VALUE(:) = NODAL_FORCES_VALUE(:) / DFLOAT(N_PASSOS)
    CARGA(:,1) = U(:)
    CARGA(:,2) = T(:)
    CARGA_NODAL(:) = NODAL_FORCES_VALUE(:)
    IF (TIPO_PASS_RED.EQ.1) PASS_RAT=1.750D0/(DFLOAT(N_PASSOS)*(DFLOAT(N_PASSOS-1)))
    IF (TIPO_PASS_RED.EQ.2) PASS_RAT=1.90D0/(DFLOAT(N_PASSOS)*(DFLOAT(N_PASSOS-1)))
    ! --------------------
    
    ! ABRINDO TXT PARA ESCREVER RESUTADOS AO LONGO DO PROCESSO
    OPEN(17,FILE='Output_data/REINFORCEMENT_BONDSLIP_DATA.TXT',STATUS='UNKNOWN')
    
    DO PASSO=1,N_PASSOS ! _____________________ Load step loop ____________________
        
        ! ATUALIZANDO CARGAS
        IF (TIPO_PASS_RED.EQ.1) THEN
            U(:) = CARGA(:,1)*(1.0D0/(8.0D0*N_PASSOS)+DFLOAT(PASSO-1)*PASS_RAT)
            T(:) = CARGA(:,2)*(1.0D0/(8.0D0*N_PASSOS)+DFLOAT(PASSO-1)*PASS_RAT)
            NODAL_FORCES_VALUE(:) = CARGA_NODAL(:)*(1.0D0/(8.0D0*N_PASSOS)+DFLOAT(PASSO-1)*PASS_RAT)
        ELSE IF (TIPO_PASS_RED.EQ.2) THEN
            U(:) = CARGA(:,1)*(1.0D0/(20.0D0*N_PASSOS)+DFLOAT(PASSO-1)*PASS_RAT)
            T(:) = CARGA(:,2)*(1.0D0/(20.0D0*N_PASSOS)+DFLOAT(PASSO-1)*PASS_RAT)
            NODAL_FORCES_VALUE(:) = CARGA_NODAL(:)*(1.0D0/(20.0D0*N_PASSOS)+DFLOAT(PASSO-1)*PASS_RAT)
        ENDIF
        
        !F_DES_ANT = 0.0D0
        !F_DES = 0.0D0
        !F_ACUM = 0.0D0
    
        !NI_SLIP_MEF=0  ! DUVIA, ZERO OU NAO
        !NI_SLIP_MEF=NI_SLIP_CONVERG
        NI_CHANGE=.TRUE.
        ITER_SLIP = 0
        ERRO = 10.0D0
                
        DO WHILE (ERRO.GT.SLIP_TOL) ! ________________ Iterations loop ______________
            
            ITER_SLIP = ITER_SLIP + 1
            NI_SLIP_MEF_ANT=NI_SLIP_MEF
            
            ! Aplicando cargas
            CALL REINF_SOLVER_SLIP(PASSO,ITER_SLIP,-DELTA_F_C_GLOBAL,NI_CHANGE)
            
            ! Calculando força de contato em cada nó, eixo axial local
            CALL ROTATE_TO_LOCAL(-P_MEF,D_P_LOCAL)
            
            ! Verificacao do modelo de aderencia para cada nó
            KK_G=0
            !NI_SLIP_MEF=0  ! NAO ZERA
            !NI_SLIP_MEF=NI_SLIP_CONVERG  ! NAO VOLTA AO CONVERGIDO
            DO I=1,N_ELEMENTOS_MEF
                DO J=1,ORDEM_MEF(I)+1
                    
                    KK_G=KK_G+1
                    
                    !AUX1=-P_MEF_ACUM(CONECTI_MEF(I,J))+D_P_LOCAL(KK_G)     
                    AUX1=-P_MEF_ACUM_LOCAL(KK_G)+D_P_LOCAL(KK_G)
                    
                    IF (AUX1.GT.0.0D0) THEN
                        F_MOD=F_MODELO(S_ACUM(KK_G),AUX1)       
                    ELSE
                        F_MOD=-F_MODELO(S_ACUM(KK_G),AUX1)       
                    ENDIF
                    
                    !IF (SLIP_LAW.EQ.'HUANG2') THEN
                    AUX2=dabs(AUX1)-dabs(F_MOD)
                    IF (dabs(AUX2).LT.SLIP_TOL) THEN
                        DELTA_F_C(KK_G)=0.0D0
                    ELSE
                        DELTA_F_C(KK_G)=-(AUX1-F_MOD)
                        IF (NI_SLIP_MEF(CONECTI_MEF(I,J)).NE.1) NI_SLIP_MEF(CONECTI_MEF(I,J))=1
                    ENDIF
                    !ELSE
                    !    !AUX2=DABS(D_P_LOCAL(KK_G))-DABS(F_MODELO(DELTA_S(CONECTI_MEF(I,J))))
                    !    AUX2=DABS(AUX1)-DABS(F_MOD)
                    !    IF (AUX2.LT.SLIP_TOL) THEN
                    !        DELTA_F_C(KK_G)=0.0D0
                    !    ELSE
                    !        !DELTA_F_C(KK_G)=-(D_P_LOCAL(KK_G)-F_MODELO(DELTA_S(CONECTI_MEF(I,J))))
                    !        DELTA_F_C(KK_G)=-(AUX1-F_MOD)  
                    !        ! Guarda que esse nó descolou
                    !        IF (NI_SLIP_MEF(CONECTI_MEF(I,J)).NE.1) NI_SLIP_MEF(CONECTI_MEF(I,J))=1
                    !    ENDIF
                    !ENDIF
                    !AUX2=(AUX1)-(F_MOD)
                    !IF (DABS(AUX2).LT.SLIP_TOL) THEN
                    !    DELTA_F_C(KK_G)=0.0D0
                    !ELSE
                    !    DELTA_F_C(KK_G)=-AUX2
                    !    IF (NI_SLIP_MEF(CONECTI_MEF(I,J)).NE.1) NI_SLIP_MEF(CONECTI_MEF(I,J))=1
                    !ENDIF                
                ENDDO
            ENDDO
            
            ! Verificando mudança no ni_slip_mef
            NI_CHANGE=.FALSE.
            DO I=1,N_NOS_MEF
                IF (NI_SLIP_MEF(I).NE.NI_SLIP_MEF_ANT(I)) NI_CHANGE=.TRUE.
            ENDDO
            
            ! Calculando delta_f_c nos nós globais
            CALL ROTATE_TO_GLOBAL(DELTA_F_C,DELTA_F_C_GLOBAL)
            
            ! Acumulando força de aderencia na matriz
            P_MEF_ACUM_LOCAL(:)=P_MEF_ACUM_LOCAL(:)-D_P_LOCAL(:)-DELTA_F_C(:)
            P_MEF_ACUM(:)=P_MEF_ACUM(:)+P_MEF(:)-DELTA_F_C_GLOBAL(:)
            
            ! Calculando erro
            CALL CALCULA_ERRO_SLIP(DELTA_F_C,ERRO)  
                                    
        ENDDO ! ______________________ End of iterations loop _____________________
        
        NI_SLIP_CONVERG=NI_SLIP_MEF
        
        WRITE(*,*)''
        WRITE(*,'(A,I3,A,I5,A,F15.7)')'STEP:',PASSO,'  N_ITERATIONS:',ITER_SLIP,'  FINAL ERROR:',ERRO
        CALL TXT_OUTPUT_BOND_SLIP(PASSO,ITER_SLIP,ERRO)
        CALL NODE_HISTORIC_PRINT(PASSO)
        
    ENDDO ! ___________________________ End of load step loop _____________________
    
    
    ! SAVING MATRIXES AND VECTORS TO CONTINUE THE CODE
    U_MEF(:) = U_MEF_ACUM(:)
    P_MEF(:) = P_MEF_ACUM(:)
    U(:) = U_ACUM_PL(:)
    T(:) = P_ACUM_PL(:)    
    ALLOCATE(S_AUX(TOTAL_ELEM_NODES_MEF))
    S_AUX=S_ACUM
    DEALLOCATE(S_ACUM)
    ALLOCATE(S_ACUM(3*N_NOS_MEF))
    CALL ROTATE_TO_GLOBAL(S_AUX,S_ACUM)
    
    CALL REINFORCEMENTS_STRESS
    
    CLOSE (17)    
    
END SUBROUTINE
        
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

SUBROUTINE ROTATE_TO_LOCAL(P_GLOBAL,P_LOCAL_TOTAL)
    
    !
    ! CALCULATE CONTACT FORCE IN THE AXIAL DIRECTION OF THE ELEMENT FOR ALL FE NODES
    !
    
    USE REINFORCEMENTS
    
    IMPLICIT NONE
        
    INTEGER::I,J,K,N_NO,NODE_G
    
    REAL*8::X,Y,Z,AUX(4)
    
    REAL*8,INTENT(IN)::P_GLOBAL(3*N_NOS_MEF)
    REAL*8,INTENT(OUT)::P_LOCAL_TOTAL(TOTAL_ELEM_NODES_MEF)
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::VALORES,MR
    REAL*8,DIMENSION(:),ALLOCATABLE::P_LOCAL,D_PHI_MEF,P_LOCAL_G
    
    NODE_G=0
    DO I=1,N_ELEMENTOS_MEF
        
        N_NO=ORDEM_MEF(I)+1
        ALLOCATE(D_PHI_MEF(N_NO),VALORES(N_NO,3),MR(N_NO,3*N_NO),P_LOCAL(N_NO),P_LOCAL_G(3*N_NO))
        
        MR=0.0D0
        DO J=1,N_NO
            ! CALCULANDO MATRIZ DE ROTACAO MR
            CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSI_MEF(I,J),(ORDEM_MEF(I)+1),I,VALORES,X,Y,Z,D_PHI_MEF)
            AUX=0.0D0
            DO K=1,N_NO
		        AUX(1)=AUX(1)+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),1)
		        AUX(2)=AUX(2)+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),2)
                AUX(3)=AUX(3)+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),3)
            ENDDO
            AUX(4)=DSQRT(AUX(1)*AUX(1)+AUX(2)*AUX(2)+AUX(3)*AUX(3))
            AUX(1)=AUX(1)/AUX(4)
            AUX(2)=AUX(2)/AUX(4)
            AUX(3)=AUX(3)/AUX(4)
            !TO AVOID NUMERICAL ERRORS:            
            IF (DABS(AUX(1)).LT.TOLER) AUX(1)=0.D0
            IF (DABS(AUX(2)).LT.TOLER) AUX(2)=0.D0
            IF (DABS(AUX(3)).LT.TOLER) AUX(3)=0.D0
            MR(J,3*J-2) = AUX(1)
            MR(J,3*J-1) = AUX(2)
            MR(J,3*J) = AUX(3)
            
            ! PEGANDO VALORES DO GLOBAL
            K=CONECTI_MEF(I,J)
            P_LOCAL_G(3*J-2)=P_GLOBAL(3*K-2)
            P_LOCAL_G(3*J-1)=P_GLOBAL(3*K-1)
            P_LOCAL_G(3*J)=P_GLOBAL(3*K)        
        ENDDO
        
        ! ROTACIONANDO GLOBAL PARA LOCAL
        P_LOCAL=MATMUL(MR,P_LOCAL_G)
        
        ! GUARDANDO NA MATRIZ DE TODOS OS NOS
        DO J=1,N_NO
            NODE_G=NODE_G+1
            P_LOCAL_TOTAL(NODE_G)=P_LOCAL(J)
        ENDDO
        
        DEALLOCATE(D_PHI_MEF,MR,VALORES,P_LOCAL,P_LOCAL_G)
    ENDDO

END SUBROUTINE
                                               
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

SUBROUTINE ROTATE_TO_GLOBAL(F_C,F_C_GLOBAL)

    !
    !   CALCULATE FORCE TO BE REAPLIED (F_C_GLOBAL) GIVEN THE LOCAL DISEQUILIBRATED FORCES (F_C)
    !

    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,J,K,N_NO,KK_G,contelem_pnod[ALLOCATABLE](:)
    
    REAL*8,INTENT(IN)::F_C(TOTAL_ELEM_NODES_MEF)
    REAL*8,INTENT(OUT)::F_C_GLOBAL(3*N_NOS_MEF)
    
    REAL*8::AUX(4),X,Y,Z
    REAL*8,DIMENSION(:,:),ALLOCATABLE::VALORES,MR,MRT
    REAL*8,DIMENSION(:),ALLOCATABLE::D_PHI_MEF,F_C_L,F_C_G
    
    ALLOCATE(contelem_pnod(3*N_NOS_MEF))
    F_C_GLOBAL=0.0D0
    contelem_pnod=0
    
    KK_G=0
    DO I=1,N_ELEMENTOS_MEF
        
        N_NO=ORDEM_MEF(I)+1
        ALLOCATE(D_PHI_MEF(N_NO),VALORES(N_NO,3),MR(N_NO,3*N_NO),MRT(3*N_NO,N_NO),F_C_L(N_NO),F_C_G(3*N_NO))
        
        MR=0.0D0
        MRT=0.0D0
        F_C_L=0.0D0
        F_C_G=0.0D0
        
        DO J=1,N_NO
            ! CALCULANDO MATRIZ DE ROTACAO MR E TRANSPOSTA MRT
            CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSI_MEF(I,J),(ORDEM_MEF(I)+1),I,VALORES,X,Y,Z,D_PHI_MEF)
            AUX=0.0D0
            DO K=1,N_NO
		        AUX(1)=AUX(1)+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),1)
		        AUX(2)=AUX(2)+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),2)
                AUX(3)=AUX(3)+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),3)
            ENDDO
            AUX(4)=DSQRT(AUX(1)*AUX(1)+AUX(2)*AUX(2)+AUX(3)*AUX(3))
            AUX(1)=AUX(1)/AUX(4)
            AUX(2)=AUX(2)/AUX(4)
            AUX(3)=AUX(3)/AUX(4)
            !TO AVOID NUMERICAL ERRORS:            
            IF (DABS(AUX(1)).LT.TOLER) AUX(1)=0.D0
            IF (DABS(AUX(2)).LT.TOLER) AUX(2)=0.D0
            IF (DABS(AUX(3)).LT.TOLER) AUX(3)=0.D0
            MR(J,3*J-2) = AUX(1)
            MR(J,3*J-1) = AUX(2)
            MR(J,3*J) = AUX(3)     
            
            ! PEGANDO VALORES PARA ESTE ELEMENTO
            KK_G=KK_G+1
            F_C_L(J)=F_C(KK_G)
        ENDDO
        MRT=TRANSPOSE(MR)
        
        ! ROTACIONANDO PARA VALORES NA COORDENADA GLOBAL
        F_C_G=MATMUL(MRT,F_C_L)
        
        ! MANDANDO PARA NÓS GLOBAIS
        DO J=1,N_NO
            K=CONECTI_MEF(I,J)
            F_C_GLOBAL(3*K-2)=F_C_GLOBAL(3*K-2)+F_C_G(3*J-2)
            F_C_GLOBAL(3*K-1)=F_C_GLOBAL(3*K-1)+F_C_G(3*J-1)
            F_C_GLOBAL(3*K)=F_C_GLOBAL(3*K)+F_C_G(3*J)

            contelem_pnod(3*K-2)=contelem_pnod(3*K-2)+1   
            contelem_pnod(3*K-1)=contelem_pnod(3*K-1)+1   
            contelem_pnod(3*K)=contelem_pnod(3*K)+1   
        ENDDO
        
        DEALLOCATE(D_PHI_MEF,VALORES,MR,MRT,F_C_L,F_C_G)
    ENDDO

    ! FAZENDO MEDIA EM NÓS QUE RECEBERAM MAIS DE UM VALOR
    DO I=1,3*N_NOS_MEF
        F_C_GLOBAL(I)=F_C_GLOBAL(I)/DFLOAT(contelem_pnod(I))
    ENDDO
    DEALLOCATE(contelem_pnod)
        
    END SUBROUTINE
                
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

SUBROUTINE CALCULA_ERRO_SLIP(F_C,ERR)

    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    REAL*8,INTENT(IN)::F_C(TOTAL_ELEM_NODES_MEF)
    REAL*8,INTENT(OUT)::ERR
    
    INTEGER::I,J,K,KK
    REAL*8::AUX1,AUX2,AUX3,SOLUCAO[ALLOCATABLE](:)
        
    AUX1=0.0D0
    AUX2=0.0D0
    KK=0
    DO I=1,N_ELEMENTOS_MEF
        DO J=1,ORDEM_MEF(I)+1
            KK=KK+1
            AUX1=AUX1+(F_C(KK))**2.0D0
            K=CONECTI_MEF(I,J)
            AUX3=P_MEF_ACUM(3*K-2)**2.0D0+P_MEF_ACUM(3*K-1)**2.0D0+P_MEF_ACUM(3*K)**2.0D0
            AUX2=AUX2+AUX3
        ENDDO
    ENDDO
    
    IF (AUX2.GT.TOLER) THEN
        ERR = DSQRT(AUX1)/DSQRT(AUX2)
    ELSE
        ERR = DSQRT(AUX1)
    ENDIF
    

END SUBROUTINE    
                                                  
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

SUBROUTINE ROTATE_NI_TO_GLOBAL(NI_L,NI_G)

    !
    !   CALCULATE NI_SLIP_MEF IN THE GLOBAL SYSTEM (3*N) GIVEN NI_SLIP_MEF IN THE LOCAL SYSTEM (N)
    !

    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,J,K,N_NO,KK_G,contelem_pnod[ALLOCATABLE](:)
    
    INTEGER,INTENT(IN)::NI_L(N_NOS_MEF)
    INTEGER,INTENT(OUT)::NI_G(3*N_NOS_MEF)
    
    REAL*8::AUX(4),X,Y,Z
    REAL*8,DIMENSION(:,:),ALLOCATABLE::VALORES,MR,MRT
    REAL*8,DIMENSION(:),ALLOCATABLE::D_PHI_MEF,NI_G_TOTAL_AUX,NI_G_AUX
    
    INTEGER,DIMENSION(:),ALLOCATABLE::NI_L_AUX    
    
    ALLOCATE(NI_G_TOTAL_AUX(3*N_NOS_MEF))
    NI_G_TOTAL_AUX=0.0D0
    NI_G=0
    
    DO I=1,N_ELEMENTOS_MEF
        
        N_NO=ORDEM_MEF(I)+1
        ALLOCATE(D_PHI_MEF(N_NO),VALORES(N_NO,3),MR(N_NO,3*N_NO),MRT(3*N_NO,N_NO),NI_L_AUX(N_NO),NI_G_AUX(3*N_NO))
        
        MR=0.0D0
        MRT=0.0D0
        NI_L_AUX=0
        NI_G_AUX=0.0D0
        
        DO J=1,N_NO
            ! CALCULANDO MATRIZ DE ROTACAO MR E TRANSPOSTA MRT
            CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSI_MEF(I,J),(ORDEM_MEF(I)+1),I,VALORES,X,Y,Z,D_PHI_MEF)
            AUX=0.0D0
            DO K=1,N_NO
		        AUX(1)=AUX(1)+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),1)
		        AUX(2)=AUX(2)+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),2)
                AUX(3)=AUX(3)+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),3)
            ENDDO
            AUX(4)=DSQRT(AUX(1)*AUX(1)+AUX(2)*AUX(2)+AUX(3)*AUX(3))
            AUX(1)=AUX(1)/AUX(4)
            AUX(2)=AUX(2)/AUX(4)
            AUX(3)=AUX(3)/AUX(4)
            !TO AVOID NUMERICAL ERRORS:            
            IF (DABS(AUX(1)).LT.TOLER) AUX(1)=0.D0
            IF (DABS(AUX(2)).LT.TOLER) AUX(2)=0.D0
            IF (DABS(AUX(3)).LT.TOLER) AUX(3)=0.D0
            MR(J,3*J-2) = AUX(1)
            MR(J,3*J-1) = AUX(2)
            MR(J,3*J) = AUX(3)     
            
            ! PEGANDO VALORES PARA ESTE ELEMENTO
            NI_L_AUX(J)=NI_L(CONECTI_MEF(I,J))
        ENDDO
        MRT=TRANSPOSE(MR)
        
        ! ROTACIONANDO PARA VALORES NA COORDENADA GLOBAL
        NI_G_AUX=MATMUL(MRT,NI_L_AUX)
        
        ! MANDANDO PARA NÓS GLOBAIS
        DO J=1,N_NO
            K=CONECTI_MEF(I,J)
            NI_G_TOTAL_AUX(3*K-2)=NI_G_TOTAL_AUX(3*K-2)+NI_G_AUX(3*J-2)
            NI_G_TOTAL_AUX(3*K-1)=NI_G_TOTAL_AUX(3*K-1)+NI_G_AUX(3*J-1)
            NI_G_TOTAL_AUX(3*K)=NI_G_TOTAL_AUX(3*K)+NI_G_AUX(3*J)
        ENDDO
        
        DEALLOCATE(D_PHI_MEF,VALORES,MR,MRT,NI_L_AUX,NI_G_AUX)
    ENDDO

    ! REGULARIZNDO NI_G 
    DO I=1,3*N_NOS_MEF
        IF (DABS(NI_G_TOTAL_AUX(I)).GT.TOLER) THEN
            NI_G(I)=1
        ELSE
            NI_G(I)=0
        ENDIF
    ENDDO
        
END SUBROUTINE