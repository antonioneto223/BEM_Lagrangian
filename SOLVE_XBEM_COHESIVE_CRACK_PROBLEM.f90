	SUBROUTINE SOLVE_XBEM_COHESIVE_CRACK_PROBLEM
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
    INTEGER::I,J,K,II,KK,NEN,CONT,IT,INTERFACE_NODE,INTERFACE_NODE1,INTERFACE_NODE2,ELEM,LOCAL_NODE,POS,TOL_IT,&
    NI_BACKUP[ALLOCATABLE](:,:),CHECK_DISCONNECT
!
    REAL*8::DX,DY,DZ,DR,TOL_F_EXC,DIF_F_EXC,DIF_F_EXC_PRV,DELTA_DIF_F
    REAL*8::DU,DU_TOTAL(2*N_INTERFACE_POINTS),U_TOTAL[ALLOCATABLE](:),T_TOTAL[ALLOCATABLE](:),U_IN[ALLOCATABLE](:),T_IN[ALLOCATABLE](:)
    REAL*8::T_INTERFACE(2*N_INTERFACE_POINTS),T_REAL(2*N_INTERFACE_POINTS),T_EXC(2*N_INTERFACE_POINTS)
    REAL*8::VALUES[ALLOCATABLE](:,:),ETA(3),VJC(3),DVJC(3,2),JC,RES_GRAPH(N_COHESIVE_INCREMENTS,2)
    REAL*8::VALUES_CON_SUPPS_BACKUP[ALLOCATABLE](:),XBEM_SUPP_FINAL[ALLOCATABLE](:),TOL1,TOL_STEP_F
!
    ALLOCATE(U_TOTAL(SIZE(U,1)),T_TOTAL(SIZE(U,1)),CHS_B_CONDITIONS(SIZE(U,1)),CHS_T(SIZE(U,1)),CHS_U(SIZE(U,1)),U_IN(SIZE(U,1)),T_IN(SIZE(U,1)),NI_BACKUP(SIZE(NI,1),SIZE(NI,2)))
    ALLOCATE(VALUES_CON_SUPPS_BACKUP(SIZE(VALUES_CON_SUPPS,1)),XBEM_SUPP_REACTIONS(SIZE(VALUES_CON_SUPPS,1)),XBEM_SUPP_FINAL(SIZE(VALUES_CON_SUPPS,1)))
!
    TOL_F_EXC=5.D0   
    TOL_STEP_F=1.D-1
    TOL1=1.D-2
    DU=0.D0
    DU_TOTAL=0.D0
    CHS_U=0.D0
    CHS_T=0.D0
    U_IN=U
    T_IN=T
    U_TOTAL=0.D0
    T_TOTAL=0.D0

!
    OPEN(9,file='GRAPH.txt',status='unknown')
    DO KK=1,N_COHESIVE_INCREMENTS
        WRITE(*,*)'INCREMENT    ',KK
        IT=0
        U=U_IN/DBLE(N_COHESIVE_INCREMENTS)
        T=T_IN/DBLE(N_COHESIVE_INCREMENTS)
        VALUES_CON_SUPPS=VALUES_CON_SUPPS_BACKUP/DBLE(N_COHESIVE_INCREMENTS)
        DIF_F_EXC=1.D10
        DELTA_DIF_F=1.D10
        DO WHILE ((DABS(DIF_F_EXC).GT.TOL_F_EXC).AND.(DELTA_DIF_F.GT.TOL_STEP_F))!.AND.IT.LE.TOL_IT)
        !   STEP 0
            IF(IT.EQ.0)THEN
                CHS_U=U
                CHS_T=T
                CALL SOLVE_BVP
                CALL COLLOCATION_OUTPUT(100)
                U_TOTAL=U_TOTAL+U
                T_TOTAL=T_TOTAL+T
                XBEM_SUPP_FINAL=XBEM_SUPP_FINAL+XBEM_SUPP_REACTIONS
                CHS_B_CONDITIONS=B_CONDITIONS
                DO I=1,N_COLLOCPOINTS
                    DO J=1,3
                        IF(CHS_B_CONDITIONS(3*I-J+1).EQ.0)THEN
                            CHS_U(3*I-J+1)=0.D0
                            CHS_B_CONDITIONS(3*I-J+1)=0
                        ENDIF
                    ENDDO
                ENDDO
            ENDIF
            ! CALCULATING THE REAL TRACTION ON POINT
            DIF_F_EXC_PRV=DIF_F_EXC
            DIF_F_EXC=0.D0
            DO II=1,2
                DO I=1,N_INTERFACE_POINTS
                    POS=I+N_INTERFACE_POINTS*(II-1)
                    INTERFACE_NODE=NI(I,II)
                    ! ELEMENT THAT CONTAINS INTERFACE NODE
                    DO J=1,N_ELEM
                        NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                        DO K=1,NEN
                            IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.INTERFACE_NODE)THEN
                                ELEM=J
                                LOCAL_NODE=K
                            ENDIF
                        ENDDO
                    ENDDO
                    NEN=ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM)
                    ALLOCATE(VALUES(NEN,3))
                    DO J=1,NEN
		                VALUES(J,1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,J),1)
		                VALUES(J,2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,J),2)
		                VALUES(J,3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,J),3)
                    ENDDO
	                CALL NORMAL_OUTWARD(NEN,VALUES,QSI(ELEM,LOCAL_NODE,1),QSI(ELEM,LOCAL_NODE,2),ETA,VJC,DVJC,JC)
                    DEALLOCATE(VALUES)
                    T_INTERFACE(POS)=T_TOTAL(3*INTERFACE_NODE-2)*ETA(1)+T_TOTAL(3*INTERFACE_NODE-1)*ETA(2)+T_TOTAL(3*INTERFACE_NODE)*ETA(3)
                    CALL COHESIVE_LAW(1,DU_TOTAL(POS),T_INTERFACE(POS),T_REAL(POS))
                    ! APPLYING BOUNDARY CONDITIONS OF NONLINEAR PROBLEM
                    T_EXC(POS)=T_INTERFACE(POS)-T_REAL(POS)
                    ! ONLY FOR INTERFACE WITH ETA=(+-1,0,0)
                    CHS_T(3*INTERFACE_NODE-2)=-T_EXC(POS)*ETA(1)
                    IF(IT.EQ.0)THEN
                        CHS_T(3*INTERFACE_NODE-1)=-T(3*INTERFACE_NODE-1)
                        CHS_T(3*INTERFACE_NODE)=-T(3*INTERFACE_NODE)
                    ELSE
                        CHS_T(3*INTERFACE_NODE-1)=0.D0
                        CHS_T(3*INTERFACE_NODE)=0.D0                       
                    ENDIF
                    ! CALCULATING T_EXC NORM
                    DIF_F_EXC=DIF_F_EXC+T_EXC(POS)*T_EXC(POS)
                ENDDO
            ENDDO
            DIF_F_EXC=DSQRT(DIF_F_EXC)
            DELTA_DIF_F=DABS(DIF_F_EXC-DIF_F_EXC_PRV)
            IF(IT.EQ.0.AND.KK.EQ.1)THEN   
                !CALL COHESIVE_INFLUENCE_MATRICES
            ENDIF
            ! NONLINEAR PROBLEM: APPLYING EXCEDENT FORCE IN SEPARATED PROBLEM
            WRITE(*,*)'ITERACTION    ',IT
            WRITE(*,*)'DIF_F_EXC=    ',DIF_F_EXC
            WRITE(*,*)'DELTA_DIF_F=    ',DELTA_DIF_F
            IF((DABS(DIF_F_EXC).GT.TOL_F_EXC).AND.(DELTA_DIF_F.GT.TOL_STEP_F))THEN
                IT=IT+1
!               RELEASING COLLOC POINTS IN INTERFACE
                NI_BACKUP=NI
                DO I=1,N_COLLOCPOINTS
                    IF(DABS(CHS_T(3*I-2)).GT.TOL1)THEN
                        DO J=1,N_INTERFACE_POINTS
			                IF(I.EQ.NI(J,1)) THEN
				                NI(J,1)=0
                                NI(J,2)=0
			                ELSEIF(I.EQ.NI(J,2))THEN
					            NI(J,1)=0
                                NI(J,2)=0
                            ENDIF
			            ENDDO        
                    ENDIF
                ENDDO
                !CALL COHESIVE_SYSTEM
                CALL COHESIVE_BVP
                !CALL COHESIVE_STEP_PLOT(IT)
                NI=NI_BACKUP
                U_TOTAL=U_TOTAL+CHS_U
                T_TOTAL=T_TOTAL+CHS_T
                XBEM_SUPP_FINAL=XBEM_SUPP_FINAL+XBEM_SUPP_REACTIONS
                DO I=1,N_INTERFACE_POINTS
                    INTERFACE_NODE1=NI(I,1)
                    INTERFACE_NODE2=NI(I,2)
                    DU=DABS(CHS_U(3*INTERFACE_NODE1-2)-CHS_U(3*INTERFACE_NODE2-2))
                    DU_TOTAL(I)=DU_TOTAL(I)+DU
                    DU_TOTAL(I+N_INTERFACE_POINTS)=DU_TOTAL(I+N_INTERFACE_POINTS)+DU
                ENDDO
            ELSE
                WRITE(*,*)'END OF INCREMENT    ',KK
            ENDIF
        ENDDO
        U=U_TOTAL
        T=T_TOTAL
        CALL COLLOCATION_OUTPUT(KK)
        ! CRIAR 'GET' PARA PEGAR O NO DO PONTO QUE QUER FAZER O GRAFICO
        !RES_GRAPH(KK,1)=U_TOTAL(3*188-2)
        !RES_GRAPH(KK,2)=T_TOTAL(3*264-2)/1000
        ! EX 2 GRAPH
        RES_GRAPH(KK,1)=KK*VALUES_CON_SUPPS_BACKUP(1)/DBLE(N_COHESIVE_INCREMENTS)
        RES_GRAPH(KK,2)=0.D0
        DO I=1,NUM_CON_SUPPS*2/7
            RES_GRAPH(KK,2)=RES_GRAPH(KK,2)+XBEM_SUPP_FINAL(I)
        ENDDO
        WRITE(*,*)'---------------------------------------------'
        WRITE(*,*)'END OF INCREMENT ',KK
        WRITE(*,*)'PRESCRIBED U',RES_GRAPH(KK,1)
        WRITE(*,*)'LOAD',RES_GRAPH(KK,2)
        WRITE(*,*)'---------------------------------------------'
        WRITE(9,100)RES_GRAPH(KK,1),RES_GRAPH(KK,2)
    ENDDO
    CLOSE(9)
!
    OPEN(2,FILE="XBEM_SUPP_REACTIONS_FINAL.TXT",STATUS='UNKNOWN')
    DO I=1,NUM_CON_SUPPS
        WRITE(2,*)XBEM_SUPP_FINAL(I)
    ENDDO
    CLOSE(2,STATUS='KEEP')

!
100	FORMAT(16x,F14.9,16x,F14.9)
!    
    END SUBROUTINE SOLVE_XBEM_COHESIVE_CRACK_PROBLEM