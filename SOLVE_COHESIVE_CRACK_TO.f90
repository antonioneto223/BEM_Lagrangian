	SUBROUTINE SOLVE_COHESIVE_CRACK_TO
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
    LOGICAL::CHANGE_A,CHANGE_TANGENT,FIRST_CHANGE
!
    INTEGER::I,J,K,II,KK,NEN,CONT,IT,INTERFACE_NODE,INTERFACE_NODE1,INTERFACE_NODE2,ELEM,LOCAL_NODE,POS,DIM_A,&
    NI_BACKUP[ALLOCATABLE](:,:),NI_DISC[ALLOCATABLE](:,:),NI_FINAL[ALLOCATABLE](:,:),NI_PRV[ALLOCATABLE](:,:),DISCONNECTED_PAIR,INFO,&
    CONT_NONCHS_NODE,N_CHS_POINTS,N_INTERFACE_POINTS_BACKUP
    INTEGER::CHS_B_CONDITIONS_INTERFACE(SIZE(U,1)),CHS_B_CONDITIONS_BACKUP(SIZE(U,1))
!
    REAL*8::DX,DY,DZ,DR,DIF_F_EXC,DIF_F_EXC_PRV,TOL1,FACTOR
    REAL*8::DU,DU_TOTAL(2*N_INTERFACE_POINTS),U_TOTAL[ALLOCATABLE](:),T_TOTAL[ALLOCATABLE](:),U_IN[ALLOCATABLE](:),&
    T_IN[ALLOCATABLE](:),XBEM_SUPP_REACTIONS_EL[ALLOCATABLE](:),U_EL[ALLOCATABLE](:),T_EL[ALLOCATABLE](:)
    REAL*8::T_INTERFACE(2*N_INTERFACE_POINTS),T_REAL(2*N_INTERFACE_POINTS),T_EXC(2*N_INTERFACE_POINTS),CHS_T_LOCAL(3),CHS_COD_BACKUP(N_COLLOCPOINTS),T_EXC_CORR(2*N_INTERFACE_POINTS)
    REAL*8::VALUES[ALLOCATABLE](:,:),ETA(3),VJC(3),DVJC(3,2),JC
    REAL*8::RES_GRAPH(N_COHESIVE_INCREMENTS,3)
    REAL*8::DU_CR,DU_2,F_T,F_T2,G_F,TOL_DU,T_INTERFACE_BILINEAR,T_REAL_BILINEAR,CHS_T_TO(SIZE(T,1))
    REAL*8::CHS_F1,CHS_F2,CHS_T_NORMAL
    REAL*8::ROTATION_MATRIX_S(3,3),ROTATION_MATRIX_H(3,3),IPIV(3),T_GLOBAL(3)
    REAL*8::VALUES_CON_SUPPS_BACKUP[ALLOCATABLE](:),XBEM_SUPP_FINAL[ALLOCATABLE](:)
    REAL*8::CHS_U_BACKUP[ALLOCATABLE](:),CHS_T_BACKUP[ALLOCATABLE](:)
    REAL*8::U_N1,U_N2,CMOD_AB,P_RES(3)
!
    ALLOCATE(U_TOTAL(SIZE(U,1)),T_TOTAL(SIZE(U,1)),CHS_B_CONDITIONS(SIZE(U,1)),CHS_T(SIZE(U,1)),CHS_U(SIZE(U,1)),U_IN(SIZE(U,1)),T_IN(SIZE(U,1)),NI_BACKUP(SIZE(NI,1),SIZE(NI,2)),NI_DISC(SIZE(NI,1),SIZE(NI,2)),CHS_U_BACKUP(SIZE(U,1)),CHS_T_BACKUP(SIZE(U,1)))
    ALLOCATE(VALUES_CON_SUPPS_BACKUP(SIZE(VALUES_CON_SUPPS,1)),XBEM_SUPP_FINAL(SIZE(XBEM_SUPP_REACTIONS)))
    ALLOCATE(CHS_COD(N_COLLOCPOINTS))
    ALLOCATE(U_EL(SIZE(U,1)),T_EL(SIZE(U,1)),XBEM_SUPP_REACTIONS_EL(SIZE(XBEM_SUPP_REACTIONS)))
!
    DIF_F_EXC_PRV=0.D0
    CONT_NONCHS_NODE=0
    !TOL_F_EXC=1.0D-6
    TOL1=TOL_F_EXC*1.D-3
    TOL_DU=1.D-8
    DU=0.D0
    DU_TOTAL=0.D0
    CHS_U=0.D0
    CHS_T=0.D0
    CHS_COD=0.D0
    T_EXC=0.D0
    U_IN=U
    T_IN=T
    U_TOTAL=0.D0
    T_TOTAL=0.D0
    CHS_U_BACKUP=0.D0
    CHS_T_BACKUP=0.D0
    DISCONNECTED_PAIR=0
    FACTOR=0.D0
!
    IF(NUM_CON_SUPPS.GT.0)THEN
        XBEM_SUPP_FINAL=0.D0
        VALUES_CON_SUPPS_BACKUP=VALUES_CON_SUPPS
    ELSE
        VALUES_CON_SUPPS_BACKUP=0.D0
    ENDIF
!   STORING INTERFACE DATA
    NI_BACKUP=NI
    NI_DISC=NI
!
    OPEN(9,file='Output_data\GRAPH.txt',status='unknown')
    ! ELASTIC PREDICTION
    WRITE(*,*)'    ELASTIC PREDICTION'
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
    ENDDO
    IF(NUM_CON_SUPPS.GT.0)THEN
        VALUES_CON_SUPPS=0.D0
    ENDIF    
    DO KK=1,N_COHESIVE_INCREMENTS
        IT=0
        CHS_T_TO=0.D0
        DIF_F_EXC=1.D10
 !      ACUMULATING ELASTIC PREDICTION TO FINAL ANSWER
        U_TOTAL=U_TOTAL+U_EL
        T_TOTAL=T_TOTAL+T_EL
        IF(NUM_CON_SUPPS.GT.0)THEN
            XBEM_SUPP_FINAL=XBEM_SUPP_FINAL+XBEM_SUPP_REACTIONS_EL
        ENDIF
       
        DO WHILE ((DABS(DIF_F_EXC).GT.TOL_F_EXC))
            CONT_NONCHS_NODE=0
            T_REAL=0.D0
            ! CALCULATING THE REAL TRACTION ON POINT
            DIF_F_EXC_PRV=DIF_F_EXC
            DIF_F_EXC=0.D0
            CHS_T=0.D0
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
!                       CALCULATING NORMAL TRACTON AT INTERFACE
                        T_INTERFACE(POS)=T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(1,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(1,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(1,3)
                        CALL COHESIVE_LAW(NI_INTERFACE_NUMBER(I),DU_TOTAL(POS),T_INTERFACE(POS),T_REAL(POS))
!
                        T_EXC(POS)=T_INTERFACE(POS)-T_REAL(POS)
                        CHS_T_LOCAL(1)=-T_EXC(POS)
!
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
                        !
                        !CALL DGESV(3,1,ROTATION_MATRIX_S,3,IPIV,CHS_T_LOCAL,3,INFO)
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
                        DIF_F_EXC=DIF_F_EXC+T_EXC(POS)*T_EXC(POS)
                    ENDDO
                ELSE
                    CONT_NONCHS_NODE=CONT_NONCHS_NODE+1
                ENDIF
            ENDDO
            DIF_F_EXC=DSQRT(DIF_F_EXC)
!         
!           NONLINEAR PROBLEM: APPLYING EXCEDENT FORCE IN TANGENT OPERATOR
            IF((DABS(DIF_F_EXC).GT.TOL_F_EXC))THEN
                ALLOCATE(CHS_INT(N_INTERFACE_POINTS-CONT_NONCHS_NODE,2),TANGENT_BILINEAR(N_INTERFACE_POINTS-CONT_NONCHS_NODE))
                TANGENT_BILINEAR=1
                CONT=0
! 
                DO I=1,N_COLLOCPOINTS
                    DO J=1,N_INTERFACE_POINTS
                        IF(I.EQ.NI(J,1))THEN
                            IF(TRIM(CM_TYPE(NI_INTERFACE_NUMBER(J))).NE.'NONE')THEN
                                IF(DABS(T_EXC(J)).GT.TOL1)THEN
				                    NI_DISC(J,1)=0
                                    NI_DISC(J,2)=0
                                ENDIF
                                IF((DABS(T_EXC(J)).GE.TOL1).OR.(DU_TOTAL(J).GE.TOL_DU))THEN
                                    IF(I.EQ.NI(J,1)) THEN
                                        CONT=CONT+1
                                        CHS_INT(CONT,1)=NI(J,1)
                                        CHS_INT(CONT,2)=NI(J,2)
                                    ENDIF
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDDO
                ENDDO  
!
                ! ORGANIZING POS_CHS
                ALLOCATE(POS_CHS(SIZE(CHS_INT,1),2))
                DO I=1,SIZE(CHS_INT,1)
                    DO II=1,2
                        INTERFACE_NODE=CHS_INT(I,II)
                        DO J=1,N_INTERFACE_POINTS
                            IF(NI_BACKUP(J,II).EQ.INTERFACE_NODE)THEN
                                POS_CHS(I,II)=J+N_INTERFACE_POINTS*(II-1)
                            ENDIF
                        ENDDO
                    ENDDO
                ENDDO   
!
                !ORGANIZING TANGENT_BILINEAR
                SELECT CASE(TRIM(CM_TYPE(1)))
                CASE('LINEAR')
                    DO I=1,SIZE(CHS_INT,1)
                        INTERFACE_NODE=CHS_INT(I,1)
                        DU_CR=CMP(NI_INTERFACE_NUMBER(POS_CHS(I,1)),2)
                        IF(CHS_COD(INTERFACE_NODE).GT.DU_CR)THEN
                            TANGENT_BILINEAR(I)=2
                        ENDIF
                    ENDDO
                CASE('BILINEAR')
                    DO I=1,SIZE(CHS_INT,1)
                        INTERFACE_NODE=CHS_INT(I,1)
                        F_T=CMP(NI_INTERFACE_NUMBER(POS_CHS(I,1)),1)
                        G_F=CMP(NI_INTERFACE_NUMBER(POS_CHS(I,1)),3)
                        DU_2=0.8D0*G_F/F_T
                        DU_CR=3.6D0*G_F/F_T
                        IF(CHS_COD(INTERFACE_NODE).GT.DU_2)THEN
                            IF(CHS_COD(INTERFACE_NODE).LT.DU_CR)THEN
                                TANGENT_BILINEAR(I)=2
                            ELSE
                                TANGENT_BILINEAR(I)=3
                            ENDIF
                        ENDIF
                    ENDDO            
                END SELECT
                NI=NI_DISC
                IT=IT+1
                WRITE(*,*)'     ITERATION=',IT
                WRITE(*,*)'     DIF_F_EXC=',DIF_F_EXC 
                ! SAVING U AND COD INFORMATIONS
                CHS_COD_BACKUP=CHS_COD
                CHS_U_BACKUP=CHS_U
                CHS_T_BACKUP=CHS_T
                CALL COHESIVE_TANGENT_OPERATOR
                !CALL COHESIVE_STEP_PLOT(IT)              
!            
                ! LINEAR COHESIVE LAW: CORRECTING IF TANGENT BECOME 0 IN STEP
                SELECT CASE(TRIM(CM_TYPE(1)))
                CASE('LINEAR')
                    T_EXC_CORR=0.D0
                    CHANGE_TANGENT=.FALSE.
                    FIRST_CHANGE=.TRUE.
                    DO I=1,SIZE(CHS_INT,1)
                        F_T=CMP(NI_INTERFACE_NUMBER(POS_CHS(I,1)),1)
                        DU_CR=CMP(NI_INTERFACE_NUMBER(POS_CHS(I,1)),2)
                        IF((CHS_COD(CHS_INT(I,1)).GT.DU_CR).AND.(TANGENT_BILINEAR(I).EQ.1))THEN
                            IF(FIRST_CHANGE)THEN
                                FIRST_CHANGE=.FALSE.
                                CHS_T=CHS_T_BACKUP
                            ENDIF
                            CHANGE_TANGENT=.TRUE.
                            TANGENT_BILINEAR(I)=2
                            ! APPLYING EXCEDENT TRACTION ON TANGENT OPERATOR
                            DO II=1,2
                                INTERFACE_NODE=CHS_INT(I,II)
                                POS=POS_CHS(I,II)
                                ! ORGANIZING LOCAL ROTATION MATRIX                    
                                ROTATION_MATRIX_S=ROTATION_INTERFACE_S(:,:,POS)
                                ROTATION_MATRIX_H=ROTATION_INTERFACE_H(:,:,POS)                                
                                ! EXISTENT TRACTION VALUE
                                CHS_F1=T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(1,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(1,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(1,3)
                                ! TRACTION VALUE WITH BILINEAR LAW SECOND PART
                                ! TRACTION VALUE EQUAL TO 0
                                CHS_F2=0.D0
                                ! TRACTION EXCEDENT THAT MUST BE REAPPLIED
                                CHS_T_LOCAL(1)=-(CHS_F1-CHS_F2)
!
                                !IF(IT.EQ.1)THEN ! RETIRA CISALHAMENTO SE FISSURA FICTÍCIA ABRE
                                IF(CHS_COD(CHS_INT(I,1)).GT.DU_CR)THEN ! RETIRA CISALHAMENTO SE FISSURA REAL ABRE
                                    CHS_T_LOCAL(2)=-(T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(2,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(2,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(2,3))
                                    CHS_T_LOCAL(3)=-(T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(3,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(3,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(3,3))
                                ELSE
                                    CHS_T_LOCAL(2)=0.D0
                                    CHS_T_LOCAL(3)=0.D0                   
                                ENDIF
!                                
                                ! debug
                                T_EXC_CORR(POS)=CHS_T_LOCAL(1)
                                ! SOLVING P(GLOBAL)=R^-1*P(LOCAL)
                                ! CALL DGESV(3,1,ROTATION_MATRIX_S,3,IPIV,CHS_T_LOCAL,3,INFO)
                                ! APPLYING BOUNDARY CONDITIONS OF NONLINEAR PROBLEM
                                CHS_T(3*INTERFACE_NODE-2)=CHS_T_LOCAL(1)!+CHS_T(3*INTERFACE_NODE-2)
                                CHS_T(3*INTERFACE_NODE-1)=CHS_T_LOCAL(2)!+!CHS_T(3*INTERFACE_NODE-1)
                                CHS_T(3*INTERFACE_NODE)=CHS_T_LOCAL(3)!+CHS_T(3*INTERFACE_NODE)
                            ENDDO
!
                                      
                        ELSEIF((CHS_COD(CHS_INT(I,1)).LE.DU_CR).AND.(TANGENT_BILINEAR(I).EQ.2))THEN
                            IF(FIRST_CHANGE)THEN
                                FIRST_CHANGE=.FALSE.
                                CHS_T=CHS_T_BACKUP
                            ENDIF
                            CHANGE_TANGENT=.TRUE.
                            TANGENT_BILINEAR(I)=1
                            ! APPLYING EXCEDENT TRACTION ON TANGENT OPERATOR
                            DO II=1,2
                                INTERFACE_NODE=CHS_INT(I,II)
                                POS=POS_CHS(I,II)
                                ! ORGANIZING LOCAL ROTATION MATRIX                    
                                ROTATION_MATRIX_S=ROTATION_INTERFACE_S(:,:,POS)
                                ROTATION_MATRIX_H=ROTATION_INTERFACE_H(:,:,POS)     
                                ! EXISTENT TRACTION VALUE
                                CHS_F1=0.D0!T_TOTAL(3*INTERFACE_NODE-2)*ETA(1)+T_TOTAL(3*INTERFACE_NODE-1)*ETA(2)+T_TOTAL(3*INTERFACE_NODE)*ETA(3)
                                ! NEGATIVE TRACTION VALUE ASSOCIATED TO COD PROJECTION ON LINEAR COHESIVE LAW
                                CHS_F2=F_T*(1.D0-CHS_COD_BACKUP(INTERFACE_NODE)/DU_CR)
                                ! TRACTION EXCEDENT THAT MUST BE REAPPLIED
                                CHS_T_LOCAL(1)=-(CHS_F1-CHS_F2)
!
                                !IF(IT.EQ.1)THEN ! RETIRA CISALHAMENTO SE FISSURA FICTÍCIA ABRE
                                IF(CHS_COD(CHS_INT(I,1)).GT.DU_CR)THEN ! RETIRA CISALHAMENTO SE FISSURA REAL ABRE
                                    CHS_T_LOCAL(2)=-(T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(2,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(2,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(2,3))
                                    CHS_T_LOCAL(3)=-(T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(3,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(3,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(3,3))
                                ELSE
                                    CHS_T_LOCAL(2)=0.D0
                                    CHS_T_LOCAL(3)=0.D0                   
                                ENDIF
!                                
                                ! debug
                                T_EXC_CORR(POS)=CHS_T_LOCAL(1)
                                ! SOLVING P(GLOBAL)=R^-1*P(LOCAL)
                                !CALL DGESV(3,1,ROTATION_MATRIX_S,3,IPIV,CHS_T_LOCAL,3,INFO)
                                ! APPLYING BOUNDARY CONDITIONS OF NONLINEAR PROBLEM
                                CHS_T(3*INTERFACE_NODE-2)=CHS_T_LOCAL(1)!+CHS_T(3*INTERFACE_NODE-2)
                                CHS_T(3*INTERFACE_NODE-1)=CHS_T_LOCAL(2)!+!CHS_T(3*INTERFACE_NODE-1)
                                CHS_T(3*INTERFACE_NODE)=CHS_T_LOCAL(3)!+CHS_T(3*INTERFACE_NODE)
                            ENDDO
                        ENDIF
                    ENDDO
!                   TANGENT OPERATOR
                    IF(CHANGE_TANGENT)THEN
                        CHS_COD=CHS_COD_BACKUP
                        CHS_U=CHS_U_BACKUP
                        CALL COHESIVE_TANGENT_OPERATOR
                        !CALL COHESIVE_STEP_PLOT(IT)
                    ENDIF
                CASE('BILINEAR')
                ! BILINEAR COHESIVE LAW
                    CHANGE_TANGENT=.FALSE.
                    FIRST_CHANGE=.TRUE.
                    DO I=1,SIZE(CHS_INT,1)
                        F_T=CMP(NI_INTERFACE_NUMBER(POS_CHS(I,1)),1)
                        G_F=CMP(NI_INTERFACE_NUMBER(POS_CHS(I,1)),3)
                        DU_2=0.8D0*G_F/F_T
                        DU_CR=3.6D0*G_F/F_T
                        F_T2=F_T/3.D0
!
                        IF((CHS_COD(CHS_INT(I,1)).GT.DU_2).AND.(TANGENT_BILINEAR(I).EQ.1))THEN
                            ! APPLYING TRACTION EXCEDENT ON TANGENT OPERATOR
                            IF(FIRST_CHANGE)THEN
                                FIRST_CHANGE=.FALSE.
                                CHS_T=CHS_T_BACKUP
                            ENDIF
                            CHANGE_TANGENT=.TRUE.
                            DO II=1,2
                                INTERFACE_NODE=CHS_INT(I,II)
                                POS=POS_CHS(I,II)
                                ! ORGANIZING LOCAL ROTATION MATRIX                    
                                ROTATION_MATRIX_S=ROTATION_INTERFACE_S(:,:,POS)
                                ROTATION_MATRIX_H=ROTATION_INTERFACE_H(:,:,POS)     
!
                                ! EXISTENT TRACTION VALUE
                                CHS_F1=T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(1,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(1,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(1,3)
! 
                                IF(CHS_COD(CHS_INT(I,1)).GT.DU_CR)THEN
                                    TANGENT_BILINEAR(I)=3
                                    CHS_F2=0.D0
                                ELSE
                                    TANGENT_BILINEAR(I)=2
                                    CHS_F2=F_T2*CHS_COD_BACKUP(INTERFACE_NODE)/(DU_2-DU_CR)+F_T2*(1.D0-DU_2/(DU_2-DU_CR))
                                ENDIF
!                           
                                ! TRACTION EXCEDENT THAT MUST BE REAPPLIED
                                CHS_T_LOCAL(1)=-(CHS_F1-CHS_F2)
                                !IF(IT.EQ.1)THEN ! RETIRA CISALHAMENTO SE FISSURA FICTÍCIA ABRE
                                IF(CHS_COD(CHS_INT(I,1)).GT.DU_CR)THEN ! RETIRA CISALHAMENTO SE FISSURA REAL ABRE
                                    CHS_T_LOCAL(2)=-(T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(2,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(2,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(2,3))
                                    CHS_T_LOCAL(3)=-(T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(3,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(3,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(3,3))
                                ELSE
                                    CHS_T_LOCAL(2)=0.D0
                                    CHS_T_LOCAL(3)=0.D0                   
                                ENDIF                 
!                                
                                ! SOLVING P(GLOBAL)=R^-1*P(LOCAL)
                                !CALL DGESV(3,1,ROTATION_MATRIX_S,3,IPIV,CHS_T_LOCAL,3,INFO)
                                ! APPLYING BOUNDARY CONDITIONS OF NONLINEAR PROBLEM
                                CHS_T(3*INTERFACE_NODE-2)=CHS_T_LOCAL(1)!+CHS_T(3*INTERFACE_NODE-2)
                                CHS_T(3*INTERFACE_NODE-1)=CHS_T_LOCAL(2)!+!CHS_T(3*INTERFACE_NODE-1)
                                CHS_T(3*INTERFACE_NODE)=CHS_T_LOCAL(3)!+CHS_T(3*INTERFACE_NODE)
                            ENDDO
                        ELSEIF(((CHS_COD(CHS_INT(I,1)).GT.DU_CR).OR.(CHS_COD(CHS_INT(I,1)).LT.DU_2)).AND.(TANGENT_BILINEAR(I).EQ.2))THEN
                            IF(FIRST_CHANGE)THEN
                                FIRST_CHANGE=.FALSE.
                                CHS_T=CHS_T_BACKUP
                            ENDIF
                            ! APPLYING TRACTION EXCEDENT ON TANGENT OPERATOR
                            CHANGE_TANGENT=.TRUE.
                            DO II=1,2
                                INTERFACE_NODE=CHS_INT(I,II)
                                POS=POS_CHS(I,II)
                                ! ORGANIZING LOCAL ROTATION MATRIX                    
                                ROTATION_MATRIX_S=ROTATION_INTERFACE_S(:,:,POS)
                                ROTATION_MATRIX_H=ROTATION_INTERFACE_H(:,:,POS)     
!
                                ! EXISTENT TRACTION VALUE
                                CHS_F1=T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(1,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(1,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(1,3)
                                ! TRACTION VALUE WITH CORRESPONDING COHESIVE LAW
                                IF(CHS_COD(CHS_INT(I,1)).GT.DU_CR)THEN
                                    TANGENT_BILINEAR(I)=3
                                    CHS_F2=0.D0
                                ELSEIF(CHS_COD(CHS_INT(I,1)).LT.DU_2)THEN
                                    TANGENT_BILINEAR(I)=1
                                    CHS_F2=F_T-((F_T-F_T2)/DU_2)*CHS_COD_BACKUP(INTERFACE_NODE)
                                ENDIF
!                                
                                ! TRACTION EXCEDENT THAT MUST BE REAPPLIED
                                CHS_T_LOCAL(1)=-(CHS_F1-CHS_F2)
                                !IF(IT.EQ.1)THEN ! RETIRA CISALHAMENTO SE FISSURA FICTÍCIA ABRE
                                IF(CHS_COD(CHS_INT(I,1)).GT.DU_CR)THEN ! RETIRA CISALHAMENTO SE FISSURA REAL ABRE
                                    CHS_T_LOCAL(2)=-(T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(2,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(2,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(2,3))
                                    CHS_T_LOCAL(3)=-(T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(3,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(3,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(3,3))
                                ELSE
                                    CHS_T_LOCAL(2)=0.D0
                                    CHS_T_LOCAL(3)=0.D0                   
                                ENDIF                  
!                                
                                ! SOLVING P(GLOBAL)=R^-1*P(LOCAL)
                                !CALL DGESV(3,1,ROTATION_MATRIX_S,3,IPIV,CHS_T_LOCAL,3,INFO)
                                ! APPLYING BOUNDARY CONDITIONS OF NONLINEAR PROBLEM
                                CHS_T(3*INTERFACE_NODE-2)=CHS_T_LOCAL(1)!+CHS_T(3*INTERFACE_NODE-2)
                                CHS_T(3*INTERFACE_NODE-1)=CHS_T_LOCAL(2)!+!CHS_T(3*INTERFACE_NODE-1)
                                CHS_T(3*INTERFACE_NODE)=CHS_T_LOCAL(3)!+CHS_T(3*INTERFACE_NODE)
                            ENDDO
                        ELSEIF(CHS_COD(CHS_INT(I,1)).LT.DU_CR.AND.(TANGENT_BILINEAR(I).EQ.3))THEN
                            IF(FIRST_CHANGE)THEN
                                FIRST_CHANGE=.FALSE.
                                CHS_T=CHS_T_BACKUP
                            ENDIF
                            CHANGE_TANGENT=.TRUE.
                            !APPLYING TRACTION EXCEDENT ON TANGENT OPERATOR
                            DO II=1,2
                                INTERFACE_NODE=CHS_INT(I,II)
                                POS=POS_CHS(I,II)
                                ! ORGANIZING LOCAL ROTATION MATRIX                    
                                ROTATION_MATRIX_S=ROTATION_INTERFACE_S(:,:,POS)
                                ROTATION_MATRIX_H=ROTATION_INTERFACE_H(:,:,POS)     
!
                                ! EXISTENT TRACTION VALUE
                                CHS_F1=0.D0!T_TOTAL(3*INTERFACE_NODE-2)*ETA(1)+T_TOTAL(3*INTERFACE_NODE-1)*ETA(2)+T_TOTAL(3*INTERFACE_NODE)*ETA(3)
                                ! TRACTION VALUE WITH CORRESPONDING COHESIVE LAW
                                IF(CHS_COD(CHS_INT(I,1)).LT.DU_2)THEN
                                    TANGENT_BILINEAR(I)=1
                                    CHS_F2=F_T-((F_T-F_T2)/DU_2)*CHS_COD_BACKUP(INTERFACE_NODE)
                                ELSE
                                    TANGENT_BILINEAR(I)=2
                                    CHS_F2=F_T2*CHS_COD_BACKUP(INTERFACE_NODE)/(DU_2-DU_CR)+F_T2*(1.D0-DU_2/(DU_2-DU_CR))
                                ENDIF
!                                
                                ! TRACTION EXCEDENT THAT MUST BE REAPPLIED
                                CHS_T_LOCAL(1)=-(CHS_F1-CHS_F2)
                                !IF(IT.EQ.1)THEN ! RETIRA CISALHAMENTO SE FISSURA FICTÍCIA ABRE
                                IF(CHS_COD(CHS_INT(I,1)).GT.DU_CR)THEN ! RETIRA CISALHAMENTO SE FISSURA REAL ABRE
                                    CHS_T_LOCAL(2)=-(T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(2,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(2,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(2,3))
                                    CHS_T_LOCAL(3)=-(T_TOTAL(3*INTERFACE_NODE-2)*ROTATION_MATRIX_S(3,1)+T_TOTAL(3*INTERFACE_NODE-1)*ROTATION_MATRIX_S(3,2)+T_TOTAL(3*INTERFACE_NODE)*ROTATION_MATRIX_S(3,3))
                                ELSE
                                    CHS_T_LOCAL(2)=0.D0
                                    CHS_T_LOCAL(3)=0.D0                   
                                ENDIF                 
!                                
                                ! SOLVING P(GLOBAL)=R^-1*P(LOCAL)
                                !CALL DGESV(3,1,ROTATION_MATRIX_S,3,IPIV,CHS_T_LOCAL,3,INFO)
                                ! APPLYING BOUNDARY CONDITIONS OF NONLINEAR PROBLEM
                                CHS_T(3*INTERFACE_NODE-2)=CHS_T_LOCAL(1)!+CHS_T(3*INTERFACE_NODE-2)
                                CHS_T(3*INTERFACE_NODE-1)=CHS_T_LOCAL(2)!+!CHS_T(3*INTERFACE_NODE-1)
                                CHS_T(3*INTERFACE_NODE)=CHS_T_LOCAL(3)!+CHS_T(3*INTERFACE_NODE)
                            ENDDO
                        ENDIF              
                    ENDDO
!                   TANGENT OPERATOR
                    IF(CHANGE_TANGENT)THEN
                        CHS_COD=CHS_COD_BACKUP
                        CHS_U=CHS_U_BACKUP
                        CALL COHESIVE_TANGENT_OPERATOR
                        !CALL COHESIVE_STEP_PLOT(IT)
                    ENDIF
                END SELECT
                DEALLOCATE(POS_CHS)              
!               REAPPLYING CODS AS PRESCRIBED DISPLACEMENT TO FIND P(COD)
                CHS_B_CONDITIONS_BACKUP=CHS_B_CONDITIONS
                CHS_B_CONDITIONS_INTERFACE=CHS_B_CONDITIONS
                DO I=1,N_COLLOCPOINTS
                    DO J=1,2
                        DO K=1,SIZE(CHS_INT,1)
                            IF(CHS_INT(K,J).EQ.I)THEN
                                CHS_B_CONDITIONS_INTERFACE(3*I-2)=0
                                CHS_B_CONDITIONS_INTERFACE(3*I-1)=0
                                CHS_B_CONDITIONS_INTERFACE(3*I)=0
                            ENDIF
                        ENDDO
                    ENDDO
                ENDDO
                CHS_B_CONDITIONS=CHS_B_CONDITIONS_INTERFACE
                CALL COHESIVE_BVP
                !CALL COHESIVE_STEP_PLOT(IT)
                CHS_B_CONDITIONS=CHS_B_CONDITIONS_BACKUP
!               
                U_TOTAL=U_TOTAL+CHS_U
                T_TOTAL=T_TOTAL+CHS_T
                IF(NUM_CON_SUPPS.GT.0)THEN
                    XBEM_SUPP_FINAL=XBEM_SUPP_REACTIONS+XBEM_SUPP_FINAL
                ENDIF
                CALL COLLOCATION_OUTPUT(KK)
                NI=NI_BACKUP
            ENDIF
            IF(ALLOCATED(CHS_INT))  DEALLOCATE(CHS_INT,TANGENT_BILINEAR)

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
                !RES_GRAPH(KK,1)=U_TOTAL(223)
                !RES_GRAPH(KK,2)=-T_TOTAL(313)
                ! XBEM
                RES_GRAPH(KK,1)=KK*VALUES_CON_SUPPS_BACKUP(1)/DBLE(N_COHESIVE_INCREMENTS)
                RES_GRAPH(KK,2)=T_TOTAL(3*21-2)
                ! ROTACIONADO
                !RES_GRAPH(KK,1)=DSQRT(U_TOTAL(223)**2+U_TOTAL(224)**2)
                !RES_GRAPH(KK,2)=DSQRT(T_TOTAL(313)**2+T_TOTAL(314)**2) 
        !--------------------------------------------------
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
        ! EXEMPLO GALVÉZ
!        ELSEIF(TRIM(INPUT_FILE).EQ.'E1.TXT'.OR.TRIM(INPUT_FILE).EQ.'e1.txt')THEN
!            RES_GRAPH(KK,1)=KK*VALUES_CON_SUPPS_BACKUP(1)/DBLE(N_COHESIVE_INCREMENTS)
!            RES_GRAPH(KK,2)=0.D0
!!
!            DO I=1,(NUM_CON_SUPPS-6)/5
!                RES_GRAPH(KK,2)=RES_GRAPH(KK,2)+XBEM_SUPP_FINAL(I)
!            ENDDO
!        ENDIF
        ! EXEMPLO 4
        ELSEIF(TRIM(INPUT_FILE).EQ.'E1.TXT'.OR.TRIM(INPUT_FILE).EQ.'e1.txt')THEN
            CALL CMOD_EX4(CMOD_AB)
            RES_GRAPH(KK,1)=CMOD_AB
            RES_GRAPH(KK,2)=XBEM_SUPP_FINAL(1)           
        ENDIF
        !--------------------------------------------------
        !GENERAL
        WRITE(*,*)'---------------------------------------------'
        WRITE(*,*)'END OF INCREMENT ',KK
        WRITE(*,*)'PRESCRIBED U',RES_GRAPH(KK,1)
        WRITE(*,*)'LOAD',RES_GRAPH(KK,2)
        WRITE(*,*)'---------------------------------------------'
        WRITE(9,100)RES_GRAPH(KK,1),RES_GRAPH(KK,2),IT
        !WRITE(10,100)RES_GRAPH(KK,1),RES_GRAPH(KK,3)
        !--------------------------------------------------
    ENDDO
!
    CLOSE(9)
!
100	FORMAT(16x,ES16.6,16x,ES16.6,8x,i4)
!    
    END SUBROUTINE SOLVE_COHESIVE_CRACK_TO