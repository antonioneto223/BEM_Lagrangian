	SUBROUTINE SOLVE_COHESIVE_CRACK_CO
!	
    USE ISOPARAMETRIC_MESH
    USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS	
	USE COHESIVE   	 
    USE XBEM_SUPPORT_VARIABLES
!
    IMPLICIT NONE   
!
    LOGICAL::CHANGE_A
!
    INTEGER::I,J,K,II,KK,NEN,CONT,IT,INTERFACE_NODE,INTERFACE_NODE1,INTERFACE_NODE2,ELEM,LOCAL_NODE,POS,DIM_A,&
    NI_BACKUP[ALLOCATABLE](:,:),NI_DISC[ALLOCATABLE](:,:),NI_PRV[ALLOCATABLE](:,:),DISCONNECTED_PAIR,&
    INFO,N_INTERFACE_POINTS_BACKUP,CONT_NONCHS_NODE
!
    REAL*8::DX,DY,DZ,DR,DIF_F_EXC,DIF_F_EXC_PRV,DU_CR
    REAL*8,ALLOCATABLE,DIMENSION(:)::DU_TOTAL,U_TOTAL,T_TOTAL,U_IN,T_IN,U_EL,T_EL,T_INTERFACE,&
    T_REAL,T_EXC,VALUES_CON_SUPPS_BACKUP,XBEM_SUPP_FINAL,XBEM_SUPP_REACTIONS_EL
    REAL*8::VALUES[ALLOCATABLE](:,:),ETA(3),VJC(3),DVJC(3,2),JC,TOL1
    REAL*8::ROTATION_MATRIX_S(3,3),ROTATION_MATRIX_H(3,3),CHS_T_LOCAL(3),IPIV(3),F_T,CHS_T_NORMAL,TOL_DU
    REAL*8::U_N1,U_N2,FACTOR
    REAL*8::RES_GRAPH(N_COHESIVE_INCREMENTS,3),P_RES(3)
    REAL*8::CMOD_AB
!
    ALLOCATE(U_TOTAL(SIZE(U,1)),T_TOTAL(SIZE(U,1)),CHS_B_CONDITIONS(SIZE(U,1)),CHS_T(SIZE(U,1)),CHS_U(SIZE(U,1)),U_IN(SIZE(U,1)),T_IN(SIZE(U,1)),NI_BACKUP(SIZE(NI,1),SIZE(NI,2)),NI_DISC(SIZE(NI,1),SIZE(NI,2)),NI_PRV(SIZE(NI,1),SIZE(NI,2)))
    ALLOCATE(VALUES_CON_SUPPS_BACKUP(SIZE(VALUES_CON_SUPPS,1)),XBEM_SUPP_FINAL(SIZE(XBEM_SUPP_REACTIONS)))
    ALLOCATE(DU_TOTAL(2*N_INTERFACE_POINTS),T_INTERFACE(2*N_INTERFACE_POINTS),T_REAL(2*N_INTERFACE_POINTS),T_EXC(2*N_INTERFACE_POINTS))
    ALLOCATE(U_EL(SIZE(U,1)),T_EL(SIZE(U,1)),XBEM_SUPP_REACTIONS_EL(SIZE(XBEM_SUPP_REACTIONS)))
!
    TOL_DU=1.E-8
    TOL1=TOL_F_EXC*1.D-3
    DU_TOTAL=0.D0
    U_IN=U
    T_IN=T
    U_TOTAL=0.D0
    T_TOTAL=0.D0
    DISCONNECTED_PAIR=0
    IF(NUM_CON_SUPPS.GT.0)THEN
        XBEM_SUPP_FINAL=0.D0
        VALUES_CON_SUPPS_BACKUP=VALUES_CON_SUPPS
    ELSE
        VALUES_CON_SUPPS_BACKUP=0.D0
    ENDIF
!   STORING INTERFACE DATA
    NI_BACKUP=NI
    NI_DISC=NI
    NI_PRV=NI
!
    OPEN(9,file='Output_data\GRAPH.txt',status='unknown')
    !OPEN(10,file='Output_data\GRAPH_REACTIONS.txt',status='unknown')
    ! ELASTIC PREDICTION
    !WRITE(*,*)'ITERATION=',KK
    IT=0
    U=U_IN/DBLE(N_COHESIVE_INCREMENTS)
    T=T_IN/DBLE(N_COHESIVE_INCREMENTS)
    VALUES_CON_SUPPS=VALUES_CON_SUPPS_BACKUP/DBLE(N_COHESIVE_INCREMENTS)
    CHS_U=U
    CHS_T=T
    NI_DISC=NI
    N_INTERFACE_POINTS_BACKUP=N_INTERFACE_POINTS
    ! RELEASING DISCONNECTED INTERFACE POINTS
    ! OBS: TRAÇÃO PURA NÃO DÁ ZERO SE ESSA PARTE TIVER OPERANDO
!    DO I=1,N_INTERFACE_POINTS_BACKUP
!        IF(DU_TOTAL(I).GT.DU_CR)THEN
			!NI_DISC(I,1)=0
!            NI_DISC(I,2)=0       
!            N_INTERFACE_POINTS=N_INTERFACE_POINTS-1
!        ENDIF
!    ENDDO
!    NI=NI_DISC
!    N_INTERFACE_POINTS=N_INTERFACE_POINTS_BACKUP
!    NI=NI_BACKUP
!
    !WRITE(*,*)'ELASTIC PREDICTION'
    CALL SOLVE_BVP
    U_EL=U
    T_EL=T
    IF(NUM_CON_SUPPS.GT.0)THEN
        XBEM_SUPP_REACTIONS_EL=XBEM_SUPP_REACTIONS
    ENDIF
    
    CALL COLLOCATION_OUTPUT(-1)
    ! FIXING PRESCRIBED DISPLACEMENT AS ZERO
    CHS_B_CONDITIONS=B_CONDITIONS
    DO I=1,N_COLLOCPOINTS
        DO J=1,3
            IF(CHS_B_CONDITIONS(3*I-J+1).EQ.0)THEN
                CHS_U(3*I-J+1)=0.D0
                CHS_B_CONDITIONS(3*I-J+1)=0
            ENDIF
        ENDDO
        IF(NUM_CON_SUPPS.GT.0)THEN
            VALUES_CON_SUPPS=0.d0
        ENDIF
    ENDDO
!   COHESIVE INCREMENTAL ANALYSIS
    DO KK=1,N_COHESIVE_INCREMENTS
!       ACUMULATING ELASTIC PREDICTION TO FINAL ANSWER
        U_TOTAL=U_TOTAL+U_EL
        T_TOTAL=T_TOTAL+T_EL
        IF(NUM_CON_SUPPS.GT.0)THEN
            XBEM_SUPP_FINAL=XBEM_SUPP_FINAL+XBEM_SUPP_REACTIONS_EL
        ENDIF
        DIF_F_EXC=1.D10
        IT=0
        ! ESQUEMA PARA COMEÇAR A ANÁLISE JÁ PERTO DO SOFTENING
        FACTOR=0.D0
        !IF(KK.EQ.1)THEN
        !    U_TOTAL=U_TOTAL+U_EL*FACTOR
        !    T_TOTAL=T_TOTAL+T_EL*FACTOR
        !    IF(NUM_CON_SUPPS.GT.0)THEN
        !        XBEM_SUPP_FINAL=XBEM_SUPP_FINAL+XBEM_SUPP_REACTIONS_EL*FACTOR
        !    ENDIF            
        !ENDIF
        DO WHILE ((DABS(DIF_F_EXC).GT.TOL_F_EXC))
            CHS_U=0.D0
            CHS_T=0.D0
            T_REAL=0.D0
            T_EXC=0.D0
            CONT_NONCHS_NODE=0
            ! CALCULATING THE REAL TRACTION ON POINT
            DIF_F_EXC_PRV=DIF_F_EXC
!
!!$OMP PARALLEL PRIVATE(I,II,POS,INTERFACE_NODE,ROTATION_MATRIX_S,ROTATION_MATRIX_H,CHS_T_LOCAL,IPIV,INFO)
!!$OMP DO SCHEDULE(DYNAMIC)  
            !WRITE(*,*)'EXCEDENT TRACTION'
            DO I=1,N_INTERFACE_POINTS
                IF(TRIM(CM_TYPE(NI_INTERFACE_NUMBER(I))).NE.'NONE')THEN
                    DO II=1,2
                        POS=I+N_INTERFACE_POINTS*(II-1)
                        INTERFACE_NODE=NI(I,II)
!
                        ! ORGANIZING LOCAL ROTATION MATRIX                    
                        ROTATION_MATRIX_S=ROTATION_INTERFACE_S(:,:,POS)
                        ROTATION_MATRIX_H=ROTATION_INTERFACE_H(:,:,POS)
!
                        ! CALCULATING NORMAL TRACTON AT INTERFACE
                        T_INTERFACE(POS)=T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(1,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(1,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(1,3)
                        CALL COHESIVE_LAW(NI_INTERFACE_NUMBER(I),DU_TOTAL(POS),T_INTERFACE(POS),T_REAL(POS))
!
                        T_EXC(POS)=T_INTERFACE(POS)-T_REAL(POS)
                        CHS_T_LOCAL(1)=-T_EXC(POS)
                        ! COHESIVE LAW PARAMETER
                        DU_CR=CMP(NI_INTERFACE_NUMBER(I),2)
                        !
                        !IF(IT.EQ.0)THEN ! RETIRA CISALHAMENTO SE FISSURA FICTÍCIA ABRE
                        IF(DU_TOTAL(I).GT.DU_CR)THEN ! RETIRA CISALHAMENTO SE FISSURA REAL ABRE
                            CHS_T_LOCAL(2)=-(T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(2,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(2,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(2,3))
                            CHS_T_LOCAL(3)=-(T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(3,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(3,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(3,3))
                        ELSE
                            CHS_T_LOCAL(2)=0.D0
                            CHS_T_LOCAL(3)=0.D0                   
                        ENDIF
                        ! SOLVING P(GLOBAL)=R^-1*P(LOCAL)
                        CALL DGESV(3,1,ROTATION_MATRIX_S,3,IPIV,CHS_T_LOCAL,3,INFO)
                        !
                        ! APPLYING BOUNDARY CONDITIONS OF NONLINEAR PROBLEM
                        CHS_T(3*INTERFACE_NODE-2)=CHS_T_LOCAL(1)
                        CHS_T(3*INTERFACE_NODE-1)=CHS_T_LOCAL(2)
                        CHS_T(3*INTERFACE_NODE)=CHS_T_LOCAL(3)
                        !! VERIFYING IF NODE STILL IS COHESIVE OR ELASTIC
                        IF(II.EQ.1)THEN
                            IF((DABS(T_EXC(POS)).LT.TOL1).AND.DU_TOTAL(I).LT.TOL_DU)THEN
                                CONT_NONCHS_NODE=CONT_NONCHS_NODE+1
                            ENDIF
                        ENDIF
                    ENDDO
                ELSE
                    CONT_NONCHS_NODE=CONT_NONCHS_NODE+1
                ENDIF
            ENDDO
!!$OMP END DO
!!$OMP END PARALLEL  
!
            ! CALCULATING T_EXC NORM
            DIF_F_EXC=0.D0
            DO I=1,2*N_INTERFACE_POINTS
                DIF_F_EXC=DIF_F_EXC+T_EXC(I)*T_EXC(I)
            ENDDO
!
            DIF_F_EXC=DSQRT(DIF_F_EXC)
            ! CHECKING IF CHS_INV_A NEEDS TO BE CHANGED
            ! NI_PRV STORES LAST INTERFACE INFO TO CHECK IF IT CHANGED OR NOT
            NI_PRV=NI_DISC
            !NI_DISC=NI_BACKUP            
!
            ! NONLINEAR PROBLEM: APPLYING EXCEDENT FORCE IN SEPARATED PROBLEM
            IF((DABS(DIF_F_EXC).GT.TOL_F_EXC))THEN
                IF(ALLOCATED(CHS_INT)) DEALLOCATE(CHS_INT)
                ALLOCATE(CHS_INT(N_INTERFACE_POINTS-CONT_NONCHS_NODE,2))
                CONT=0
                DO I=1,N_COLLOCPOINTS
                    DO J=1,N_INTERFACE_POINTS
                        IF(I.EQ.NI(J,1))THEN
                            IF(TRIM(CM_TYPE(NI_INTERFACE_NUMBER(J))).NE.'NONE')THEN
                                IF(DABS(T_EXC(J)).GT.TOL1)THEN
				                    NI_DISC(J,1)=0
                                    NI_DISC(J,2)=0
                                ENDIF
                                ! COESIVO ROTACIONADO
                                !IF((DABS(T_EXC(J)).GE.TOL1).OR.(DU_TOTAL(J).GE.TOL_DU))THEN
                                !    IF(I.EQ.NI(J,1)) THEN
                                !        CONT=CONT+1
                                !        CHS_INT(CONT,1)=NI(J,1)
                                !        CHS_INT(CONT,2)=NI(J,2)
                                !    ENDIF
                                !ENDIF
                                ! ACABA COESIVO ROTACIONADO
                            ENDIF
                        ENDIF
                    ENDDO
                ENDDO 
!
                CHANGE_A=.FALSE.
                DISCONNECTED_PAIR=0
                DO I=1,N_INTERFACE_POINTS
                    IF(NI_DISC(I,1).NE.NI_PRV(I,1))THEN
                        CHANGE_A=.TRUE.
                    ENDIF
                    IF(NI_DISC(I,1).EQ.0)THEN
                        DISCONNECTED_PAIR=DISCONNECTED_PAIR+1
                    ENDIF
                ENDDO
!
                NI=NI_DISC
                IF(CHANGE_A)THEN
                    WRITE(*,*)'NUMBER OF DISCONNECTED COLLOC POINTS:',DISCONNECTED_PAIR
                    N_CHS_INTERFACE_POINTS=N_INTERFACE_POINTS-DISCONNECTED_PAIR
                    CALL COHESIVE_INFLUENCE_MATRICES(DISCONNECTED_PAIR)
                ELSEIF(IT.EQ.0.AND.KK.EQ.1)THEN
                    N_CHS_INTERFACE_POINTS=N_INTERFACE_POINTS
                    CALL COHESIVE_INFLUENCE_MATRICES(DISCONNECTED_PAIR)
                ENDIF
                IT=IT+1
                WRITE(*,*)'     ITERATION=',IT
                WRITE(*,*)'     DIF_F_EXC=',DIF_F_EXC
                CALL COHESIVE_SYSTEM
                !CALL COHESIVE_STEP_PLOT(IT)
                U_TOTAL=U_TOTAL+CHS_U
                T_TOTAL=T_TOTAL+CHS_T
                IF(NUM_CON_SUPPS.GT.0)THEN
                    XBEM_SUPP_FINAL=XBEM_SUPP_FINAL+XBEM_SUPP_REACTIONS
                ENDIF
                NI=NI_BACKUP
! 
            ENDIF
!!$OMP PARALLEL PRIVATE(I,II,POS,INTERFACE_NODE,ROTATION_MATRIX_S,ROTATION_MATRIX_H,U_N1,U_N2)
!!$OMP DO SCHEDULE(DYNAMIC)    
            DO I=1,N_INTERFACE_POINTS
                DO II=1,2
                    INTERFACE_NODE=NI(I,II)
                    POS=I+N_INTERFACE_POINTS*(II-1)
                    ! ROTATION MATRICES
                    ROTATION_MATRIX_S=ROTATION_INTERFACE_S(:,:,POS)
                    ROTATION_MATRIX_H=ROTATION_INTERFACE_H(:,:,POS)
                    IF (II.EQ.1)THEN
                        U_N1=U_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(1,1)+U_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(1,2)+U_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(1,3)
                    ELSE
                        U_N2=U_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(1,1)+U_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(1,2)+U_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(1,3)
                    ENDIF
                ENDDO
                DU_TOTAL(I)=-U_N1-U_N2
                DU_TOTAL(I+N_INTERFACE_POINTS)=-U_N1-U_N2
            ENDDO
!!$OMP END DO
!!$OMP END PARALLEL 
            NI=NI_BACKUP
        ENDDO
        U=U_TOTAL
        T=T_TOTAL
        ! CHANGING NUMBERS BELOW TOLERANCE TO ZERO - PLOT ISSUES
        DO I=1,SIZE(U,1)
            IF(DABS(U(I)).LT.1.D-16)THEN
                U(I)=0.D0
            ENDIF
            IF(DABS(T(I)).LT.1.D-16)THEN
                T(I)=0.D0
            ENDIF
        ENDDO
        
        CALL COLLOCATION_OUTPUT(KK)
        CALL INTERFACE_OUTPUT(KK,DU_TOTAL)
        !--------------------------------------------------
        ! EXEMPLO 1
        IF(TRIM(INPUT_FILE).EQ.'E2.TXT'.OR.TRIM(INPUT_FILE).EQ.'e2.txt')THEN
            RES_GRAPH(KK,1)=U_TOTAL(223)
            RES_GRAPH(KK,2)=-T_TOTAL(313)
            !CALL RESULTING_FORCE_PRESCRIBED_DISPL(1,CMOD_AB,'U')
            !--------------------------------------------------
            ! EXEMPLO 1 XBEM
            !RES_GRAPH(KK,1)=KK*VALUES_CON_SUPPS_BACKUP(1)/DBLE(N_COHESIVE_INCREMENTS)
            !RES_GRAPH(KK,2)=T_TOTAL(3*21-2)
            !DO I=1,16
            !    RES_GRAPH(KK,2)=RES_GRAPH(KK,2)+XBEM_SUPP_FINAL(I)
            !ENDDO  
        !!--------------------------------------------------
        !! EXEMPLOS 2 E 3
!        ELSEIF(TRIM(INPUT_FILE).EQ.'E1.TXT'.OR.TRIM(INPUT_FILE).EQ.'e1.txt')THEN
!            !RES_GRAPH(KK,1)=KK*VALUES_CON_SUPPS_BACKUP(1)/DBLE(N_COHESIVE_INCREMENTS)
!            !adaptação pra começar a análise mais perto da carga de pico
!            IF(NUM_CON_SUPPS.GT.0)THEN
!                IF(KK.EQ.1)THEN
!                    RES_GRAPH(KK,1)=(FACTOR+1.D0)*VALUES_CON_SUPPS_BACKUP(1)/DBLE(N_COHESIVE_INCREMENTS)    
!                ELSE
!                    RES_GRAPH(KK,1)=RES_GRAPH(KK-1,1)+VALUES_CON_SUPPS_BACKUP(1)/(DBLE(N_COHESIVE_INCREMENTS))
!                ENDIF
!!            
!                RES_GRAPH(KK,2)=0.D0
!                DO I=1,2*(NUM_CON_SUPPS-2)/5
!                    RES_GRAPH(KK,2)=RES_GRAPH(KK,2)+XBEM_SUPP_FINAL(I)
!                ENDDO
!                DO I=2*(NUM_CON_SUPPS-2)/5+1,4*(NUM_CON_SUPPS-2)/5
!                    RES_GRAPH(KK,3)=RES_GRAPH(KK,3)+XBEM_SUPP_FINAL(I)
!                ENDDO
!            ELSE
!                ! EXEMPLO 2 APOIOS DISTRIBUÍDOS
!                ! ESPESSURA UNITÁRIA
!                ! RES_GRAPH(KK,1)=KK*(-2.d-4)/DBLE(N_COHESIVE_INCREMENTS)
!                ! ESPESSURA SALEH E ALIABADI
!                RES_GRAPH(KK,1)=KK*(-1.d-3)/DBLE(N_COHESIVE_INCREMENTS)
!                CALL RESULTING_FORCE_PRESCRIBED_DISPL(P_RES)
!                RES_GRAPH(KK,2)=P_RES(2)
!            ENDIF
!        ENDIF
        !! EXEMPLO 4 - BROKENSHIRE
        ELSEIF(TRIM(INPUT_FILE).EQ.'E1.TXT'.OR.TRIM(INPUT_FILE).EQ.'e1.txt')THEN
            CALL CMOD_EX4(CMOD_AB)
            RES_GRAPH(KK,1)=CMOD_AB
            RES_GRAPH(KK,2)=XBEM_SUPP_FINAL(1)           
        ENDIF
        !   PULLOUT
        !ELSEIF(TRIM(INPUT_FILE).EQ.'E1.TXT'.OR.TRIM(INPUT_FILE).EQ.'e1.txt')THEN
        !    CALL RESULTING_FORCE_PRESCRIBED_DISPL(P_RES)
        !    RES_GRAPH(KK,1)=0.4*KK/N_COHESIVE_INCREMENTS
        !    RES_GRAPH(KK,2)=P_RES(3)
        !ENDIF
        !!--------------------------------------------------
        !!GENERAL
        WRITE(*,*)'---------------------------------------------'
        WRITE(*,*)'END OF INCREMENT ',KK
        WRITE(*,*)'U',RES_GRAPH(KK,1)
        WRITE(*,*)'LOAD',RES_GRAPH(KK,2)
        !WRITE(*,*)'LOAD_REACTION',RES_GRAPH(KK,3)
        WRITE(*,*)'---------------------------------------------'
        WRITE(9,100)RES_GRAPH(KK,1),RES_GRAPH(KK,2),IT
        !WRITE(10,100)RES_GRAPH(KK,1),RES_GRAPH(KK,3),IT
        !!--------------------------------------------------

    ENDDO
!
    CLOSE(9)
    !CLOSE(10)
!
100	FORMAT(16x,ES16.6,16x,ES16.6,8x,i4)
!    
    END SUBROUTINE SOLVE_COHESIVE_CRACK_CO
    
    SUBROUTINE CMOD_EX4(CMOD_AB)
!   
    USE ISOPARAMETRIC_MESH
    USE ANALYSIS
!
    INTEGER::ELEM_A,ELEM_B,NEN,I,J
    REAL*8::COORDS_A(3),COORDS_B(3),QSI_A(2),QSI_B(2),VALUES[ALLOCATABLE](:,:),AUX,U_A(3),U_B(3),CMOD_AB
    REAL*8::COEFFICIENTS[ALLOCATABLE](:,:),PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:)
!
!   POINT A
    COORDS_A(1)=0.1249D0
    COORDS_A(2)=0.05D0
    COORDS_A(3)=0.2089D0
    ! TESTE
    !COORDS_A(1)=0.D0
    !COORDS_A(2)=0.25D0
    !COORDS_A(3)=0.25D0 
!   
    CALL FIND_ELEMENT_WITH_PNT(COORDS_A,ELEM_A,QSI_A)
!
    NEN=ELEM_TYPE(ELEM_A)*ORDER_ELEM(ELEM_A)+(ELEM_TYPE(ELEM_A)-3)*(ORDER_ELEM(ELEM_A)-1)*POL_FAMILY(ELEM_A)
    ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2))
!
    DO J=1,NEN
        VALUES(J,:)=COORD_NODES(NODES_CONNECTIVITY(ELEM_A,J),:)
    ENDDO
!
    CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM_A,NEN,COEFFICIENTS)
    CALL SHAPE_FUNCTIONS(2,QSI_A(1),QSI_A(2),NEN,VALUES,COEFFICIENTS,AUX,AUX,AUX,PHI,DPHI,D2PHI)
!
    U_A=0.D0
    DO I=1,3
        DO J=1,NEN
            U_A(I)=U_A(I)+U(3*COLLOCPOINTS_CONNECTIVITY(ELEM_A,J)+I-3)*PHI(J)
        ENDDO
    ENDDO
!
!   POINT B   
    ! CORRETO
    COORDS_B(1)=0.1286D0
    COORDS_B(2)=0.05D0
    COORDS_B(3)=0.2126D0
    ! TESTE
    !COORDS_B(1)=4.D0
    !COORDS_B(2)=0.25D0
    !COORDS_B(3)=0.25D0   
!
    CALL FIND_ELEMENT_WITH_PNT(COORDS_B,ELEM_B,QSI_B)
!
    NEN=ELEM_TYPE(ELEM_B)*ORDER_ELEM(ELEM_B)+(ELEM_TYPE(ELEM_B)-3)*(ORDER_ELEM(ELEM_B)-1)*POL_FAMILY(ELEM_B)
    DEALLOCATE(VALUES,COEFFICIENTS,PHI,DPHI,D2PHI)
    ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2))
!
    DO J=1,NEN
        VALUES(J,:)=COORD_NODES(NODES_CONNECTIVITY(ELEM_B,J),:)
    ENDDO
!
    CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM_B,NEN,COEFFICIENTS)
    CALL SHAPE_FUNCTIONS(2,QSI_B(1),QSI_B(2),NEN,VALUES,COEFFICIENTS,AUX,AUX,AUX,PHI,DPHI,D2PHI)
!
    U_B=0.D0
    DO I=1,3
        DO J=1,NEN
            U_B(I)=U_B(I)+U(3*COLLOCPOINTS_CONNECTIVITY(ELEM_B,I)+I-3)*PHI(J)
        ENDDO
    ENDDO
!
    CMOD_AB=DSQRT(2.D0)/2.D0*(-U_A(1)-U_A(3)+U_B(1)*U_B(3))
!
    WRITE(*,*)'ELEM_A = ', ELEM_A
    WRITE(*,*)'ELEM_B = ', ELEM_B
    !WRITE(*,*)'U_A(1) = ',U_A(1)
    !WRITE(*,*)'U_A(2) = ',U_A(2)
    !WRITE(*,*)'U_A(3) = ',U_A(3)
    !WRITE(*,*)'U_B(1) = ',U_B(1)
    !WRITE(*,*)'U_B(2) = ',U_B(2)
    !WRITE(*,*)'U_B(3) = ',U_B(3)
!
    END SUBROUTINE