    SUBROUTINE XBEM_SUPPORT_CHS(A,F,DIM_XBEM_SUPP,CHS_ANALYSIS)
!
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE XBEM_SUPPORT_VARIABLES
!   
    IMPLICIT NONE
!    
    CHARACTER*1::CHS_ANALYSIS
!
    INTEGER::I,J,K,L,M,II,JJ,KK,NEN,ELEMS,NODES,DIM_XBEM_SUPP,NUM_COLLOCPOINT,ELEM_NODE,&
    NUM_CON_SUPPS_CHS,IT_I
!
    REAL*8::DX,DY,DZ,DR(3),U_AST(3,3),D_AST(9,3),R,&
    Mu,Nu,C1,C2,AUX1,KR(3,3),PI,VALUES[ALLOCATABLE](:,:),&
    ETAS(3),MAT_ETAS(3,9),VJC(3),DVJC(3,2),JC,AUX_D(3,3),&
    QSI_PNT(2),A(DIM_XBEM_SUPP,DIM_XBEM_SUPP),F(DIM_XBEM_SUPP),&
    COEFFICIENTS[ALLOCATABLE](:,:),AUX_DBL,VALUES_ETA[ALLOCATABLE](:,:),&
    PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),TOLDIST
!
    PI=DACOS(-1.D0)
    KR(1,1)=1.D0
	KR(1,2)=0.D0
    KR(1,3)=0.D0
	KR(2,1)=0.D0
	KR(2,2)=1.D0
	KR(2,3)=0.D0
	KR(3,1)=0.D0
	KR(3,2)=0.D0
	KR(3,3)=1.D0
    TOLDIST=1.0E-5
!
    DO I=1,NUM_CON_SUPPS
        IF(XBEM_CHS_BC(I).EQ."Y")THEN
            NUM_CON_SUPPS_CHS=NUM_CON_SUPPS_CHS+1
        ENDIF
    ENDDO  
!    
    IT_I=0
    IF(TRIM(CHS_ANALYSIS).EQ."N")THEN
        ALLOCATE(INT_DIR_CON_SUPPS(NUM_CON_SUPPS-NUM_CON_SUPPS_CHS))
        DO I=1,NUM_CON_SUPPS
            IF(TRIM(XBEM_CHS_BC(I)).EQ."N")THEN
                IT_I=IT_I+1
!               DEFINING DIRECTION NUMBER OF CONCENTRATED SUPPORT
                IF(DIR_CON_SUPPS(I).EQ."X")THEN
                    INT_DIR_CON_SUPPS(IT_I)=1
                ELSEIF(DIR_CON_SUPPS(I).EQ."Y")THEN
                    INT_DIR_CON_SUPPS(IT_I)=2 
                ELSEIF(DIR_CON_SUPPS(I).EQ."Z")THEN
                    INT_DIR_CON_SUPPS(IT_I)=3
                ENDIF
!               IDENTIFYING ELEMENT CONTAINING CONCENTRATED SUPPORT AND CALCULATING ADIMENSIONAL PARAMETERS OF CONCENTRATED SUPPORT POINT
                ELEM_NODE=0
                CALL FIND_ELEMENT_WITH_PNT(COORDS_CON_SUPPS(I,:),ELEM_NODE,QSI_PNT)
!
                NEN=ELEM_TYPE(ELEM_NODE)*ORDER_ELEM(ELEM_NODE)+(ELEM_TYPE(ELEM_NODE)-3)*(ORDER_ELEM(ELEM_NODE)-1)*POL_FAMILY(ELEM_NODE)
                ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2))
                DO J=1,NEN
                    VALUES(J,:)=COORD_COLLOCPOINTS(NODES_CONNECTIVITY(ELEM_NODE,J),:)
                ENDDO
!
                DO II=1,N_COLLOCPOINTS
!                   FINDING ELEMENT THAT CONTAINS COLLOCATION POINT
!                   ELEMS = ELEMENT CONTAINING COLLOCATION POINT  
                    DO K=1,N_ELEM
		                NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			            DO KK=1,NEN
				            IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.II) THEN
					            ELEMS=K
				            ENDIF
			            ENDDO
                    ENDDO
!                   CHECKING IF SUPPORT IS AT SAME DOMAIN THAT COLLOCATION POINT            
                    IF (ND(ELEM_NODE).EQ.ND(ELEMS)) THEN
                        Mu=EMP(ND(ELEM_NODE),1)/(2*(1+EMP(ND(ELEM_NODE),2)))
	                    Nu=EMP(ND(ELEM_NODE),2)
                        C1=((1.D0)/(16.D0*PI*Mu*(1.D0-Nu)))
                        C2=((1.D0)/(8.D0*PI*(1.D0-Nu)))
    !                   CALCULATING U*
        	            DX=COORDS_CON_SUPPS(I,1)-COORD_COLLOCPOINTS(II,1)
		                DY=COORDS_CON_SUPPS(I,2)-COORD_COLLOCPOINTS(II,2)
		                DZ=COORDS_CON_SUPPS(I,3)-COORD_COLLOCPOINTS(II,3)
		                R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                        IF(R.LT.TOLDIST)THEN
                            WRITE(*,*)"CONCENTRATED SUPPORT AND COLLOCATION POINT AT SAME POSITION!! CHANGE MESH!"
                            READ(*,*)
                        ENDIF
		                DR(1)=DX/R
		                DR(2)=DY/R
		                DR(3)=DZ/R
                        IF (EQ_TYPE(II).EQ."S") THEN
                            DO K=1,3
                                DO L=1,3
                                    AUX1=(3.0D0-4*Nu)*KR(K,L)
                                    U_AST(K,L)=C1*(AUX1+DR(K)*DR(L))/R
                                ENDDO
                            ENDDO
                        ELSE
                        !---------------------------------------------------------
                        !       SOURCE POINT NORMAL OUTWARD MATRIX
                        !---------------------------------------------------------
!
                            NEN=ELEM_TYPE(ELEMS)*ORDER_ELEM(ELEMS)+(ELEM_TYPE(ELEMS)-3)*(ORDER_ELEM(ELEMS)-1)*POL_FAMILY(ELEMS)
	                        ALLOCATE(VALUES_ETA(NEN,3))
	                        NODES=0
!                           NODES = LOCAL NODE
                            DO JJ=1,NEN
		                        IF(COLLOCPOINTS_CONNECTIVITY(ELEMS,JJ).EQ.II) THEN
			                        NODES=JJ
		                        ENDIF
                            ENDDO
!                       VALUES_ETA = COORDINATES OF ELEMENT NODES FOR ETA DETERMINATION	
	                        DO JJ=1,NEN
		                        VALUES_ETA(JJ,1)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,JJ),1)
		                        VALUES_ETA(JJ,2)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,JJ),2)
		                        VALUES_ETA(JJ,3)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,JJ),3)
                            ENDDO
!
	                        CALL NORMAL_OUTWARD(NEN,VALUES_ETA,QSI(ELEMS,NODES,1),QSI(ELEMS,NODES,2),ETAS,VJC,DVJC,JC)
!
                            MAT_ETAS=0.D0
	                        MAT_ETAS(1,1)=ETAS(1)
	                        MAT_ETAS(1,4)=ETAS(2)
                            MAT_ETAS(1,7)=ETAS(3)
	                        MAT_ETAS(2,2)=ETAS(1)
	                        MAT_ETAS(2,5)=ETAS(2)
	                        MAT_ETAS(2,8)=ETAS(3)
	                        MAT_ETAS(3,3)=ETAS(1)
                            MAT_ETAS(3,6)=ETAS(2)
                            MAT_ETAS(3,9)=ETAS(3)
                            DEALLOCATE(VALUES_ETA)
                            DO K=1,3
			                    DO L=1,3
			                        DO M=1,3
					                    AUX1=0.D0
					                    AUX1=KR(K,M)*DR(L)+KR(L,M)*DR(K)    
					                    AUX1=AUX1-KR(K,L)*DR(M)
					                    AUX1=AUX1*(1.D0-2.D0*Nu)
					                    AUX1=AUX1+3.D0*DR(K)*DR(L)*DR(M)
					                    D_AST((L+3*K-3),M)=(1.D0/R**2)*C2*AUX1			        
				                    ENDDO
			                    ENDDO
                            ENDDO
                            ! IN THIS CONTEXT, U_AST WAS USED JUST TO NOT CHANGE AFTER ENDIF
                            U_AST=MATMUL(MAT_ETAS,D_AST)
                        ENDIF
                    ELSE
                        U_AST=0.D0
                    ENDIF
                    A(3*II-2,3*N_COLLOCPOINTS+IT_I)=-U_AST(1,INT_DIR_CON_SUPPS(IT_I))
                    A(3*II-1,3*N_COLLOCPOINTS+IT_I)=-U_AST(2,INT_DIR_CON_SUPPS(IT_I))
                    A(3*II,3*N_COLLOCPOINTS+IT_I)=-U_AST(3,INT_DIR_CON_SUPPS(IT_I))
                ENDDO
!               CALCULATING A_X AND F_X  
                CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM_NODE,NEN,COEFFICIENTS)
                CALL SHAPE_FUNCTIONS(2,QSI_PNT(1),QSI_PNT(2),NEN,VALUES,COEFFICIENTS,AUX_DBL,AUX_DBL,AUX_DBL,PHI,DPHI,D2PHI)
                DO JJ=1,NEN
                    NUM_COLLOCPOINT=NODES_CONNECTIVITY(ELEM_NODE,JJ)
                    IF(INT_DIR_CON_SUPPS(IT_I).EQ.1)THEN
                        IF(B_CONDITIONS(3*NUM_COLLOCPOINT-2).EQ.0)THEN
                            F(3*N_COLLOCPOINTS+IT_I)=-PHI(JJ)*U(3*NUM_COLLOCPOINT-2)+VALUES_CON_SUPPS(I)
                        ELSE
                            A(3*N_COLLOCPOINTS+IT_I,3*NUM_COLLOCPOINT-2)=PHI(JJ)
                            F(3*N_COLLOCPOINTS+IT_I)=VALUES_CON_SUPPS(I)
                        ENDIF
                    ELSEIF(INT_DIR_CON_SUPPS(IT_I).EQ.2)THEN
                        IF(B_CONDITIONS(3*NUM_COLLOCPOINT-1).EQ.0)THEN
                            F(3*N_COLLOCPOINTS+IT_I)=-PHI(JJ)*U(3*NUM_COLLOCPOINT-1)+VALUES_CON_SUPPS(I)
                        ELSE
                            A(3*N_COLLOCPOINTS+IT_I,3*NUM_COLLOCPOINT-1)=PHI(JJ)
                            F(3*N_COLLOCPOINTS+IT_I)=VALUES_CON_SUPPS(I)
                        ENDIF
                    ELSE
                        IF(B_CONDITIONS(3*NUM_COLLOCPOINT).EQ.0)THEN
                            F(3*N_COLLOCPOINTS+IT_I)=-PHI(JJ)*U(3*NUM_COLLOCPOINT)+VALUES_CON_SUPPS(I)
                        ELSE
                            A(3*N_COLLOCPOINTS+IT_I,3*NUM_COLLOCPOINT)=PHI(JJ)
                            F(3*N_COLLOCPOINTS+IT_I)=VALUES_CON_SUPPS(I)
                        ENDIF
                    ENDIF
                ENDDO
            !
                DEALLOCATE(VALUES,COEFFICIENTS,PHI,DPHI,D2PHI)    
            ENDIF
        ENDDO
    ELSEIF(TRIM(CHS_ANALYSIS).EQ."Y")THEN
        ALLOCATE(INT_DIR_CON_SUPPS(NUM_CON_SUPPS))
        DO I=1,NUM_CON_SUPPS
!           DEFINING DIRECTION NUMBER OF CONCENTRATED SUPPORT
            IF(DIR_CON_SUPPS(I).EQ."X")THEN
                INT_DIR_CON_SUPPS(I)=1
            ELSEIF(DIR_CON_SUPPS(I).EQ."Y")THEN
                INT_DIR_CON_SUPPS(I)=2 
            ELSEIF(DIR_CON_SUPPS(I).EQ."Z")THEN
                INT_DIR_CON_SUPPS(I)=3
            ENDIF
!           IDENTIFYING ELEMENT CONTAINING CONCENTRATED SUPPORT AND CALCULATING ADIMENSIONAL PARAMETERS OF CONCENTRATED SUPPORT POINT
            ELEM_NODE=0
            CALL FIND_ELEMENT_WITH_PNT(COORDS_CON_SUPPS(I,:),ELEM_NODE,QSI_PNT)
!
            NEN=ELEM_TYPE(ELEM_NODE)*ORDER_ELEM(ELEM_NODE)+(ELEM_TYPE(ELEM_NODE)-3)*(ORDER_ELEM(ELEM_NODE)-1)*POL_FAMILY(ELEM_NODE)
            ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2))
            DO J=1,NEN
                VALUES(J,:)=COORD_COLLOCPOINTS(NODES_CONNECTIVITY(ELEM_NODE,J),:)
            ENDDO
!
            DO II=1,N_COLLOCPOINTS
!               FINDING ELEMENT THAT CONTAINS COLLOCATION POINT
!               ELEMS = ELEMENT CONTAINING COLLOCATION POINT  
                DO K=1,N_ELEM
		            NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			        DO KK=1,NEN
				        IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.II) THEN
					        ELEMS=K
				        ENDIF
			        ENDDO
                ENDDO
!               CHECKING IF SUPPORT IS AT SAME DOMAIN THAT COLLOCATION POINT            
                IF (ND(ELEM_NODE).EQ.ND(ELEMS)) THEN
                    Mu=EMP(ND(ELEM_NODE),1)/(2*(1+EMP(ND(ELEM_NODE),2)))
	                Nu=EMP(ND(ELEM_NODE),2)
                    C1=((1.D0)/(16.D0*PI*Mu*(1.D0-Nu)))
                    C2=((1.D0)/(8.D0*PI*(1.D0-Nu)))
!                   CALCULATING U*
        	        DX=COORDS_CON_SUPPS(I,1)-COORD_COLLOCPOINTS(II,1)
		            DY=COORDS_CON_SUPPS(I,2)-COORD_COLLOCPOINTS(II,2)
		            DZ=COORDS_CON_SUPPS(I,3)-COORD_COLLOCPOINTS(II,3)
		            R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                    IF(R.LT.TOLDIST)THEN
                        WRITE(*,*)"CONCENTRATED SUPPORT AND COLLOCATION POINT AT SAME POSITION!! CHANGE MESH!"
                        READ(*,*)
                    ENDIF
		            DR(1)=DX/R
		            DR(2)=DY/R
		            DR(3)=DZ/R
                    IF (EQ_TYPE(II).EQ."S") THEN
                        DO K=1,3
                            DO L=1,3
                                AUX1=(3.0D0-4*Nu)*KR(K,L)
                                U_AST(K,L)=C1*(AUX1+DR(K)*DR(L))/R
                            ENDDO
                        ENDDO
                    ELSE
                    !---------------------------------------------------------
                    !       SOURCE POINT NORMAL OUTWARD MATRIX
                    !---------------------------------------------------------
!
                        NEN=ELEM_TYPE(ELEMS)*ORDER_ELEM(ELEMS)+(ELEM_TYPE(ELEMS)-3)*(ORDER_ELEM(ELEMS)-1)*POL_FAMILY(ELEMS)
	                    ALLOCATE(VALUES_ETA(NEN,3))
	                    NODES=0
!                   NODES = LOCAL NODE
                        DO JJ=1,NEN
		                    IF(COLLOCPOINTS_CONNECTIVITY(ELEMS,JJ).EQ.II) THEN
			                    NODES=JJ
		                    ENDIF
                        ENDDO
    !                   VALUES_ETA = COORDINATES OF ELEMENT NODES FOR ETA DETERMINATION	
	                    DO JJ=1,NEN
		                    VALUES_ETA(JJ,1)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,JJ),1)
		                    VALUES_ETA(JJ,2)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,JJ),2)
		                    VALUES_ETA(JJ,3)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,JJ),3)
	                    ENDDO
!
	                    CALL NORMAL_OUTWARD(NEN,VALUES_ETA,QSI(ELEMS,NODES,1),QSI(ELEMS,NODES,2),ETAS,VJC,DVJC,JC)
!
                        MAT_ETAS=0.D0
	                    MAT_ETAS(1,1)=ETAS(1)
	                    MAT_ETAS(1,4)=ETAS(2)
                        MAT_ETAS(1,7)=ETAS(3)
	                    MAT_ETAS(2,2)=ETAS(1)
	                    MAT_ETAS(2,5)=ETAS(2)
	                    MAT_ETAS(2,8)=ETAS(3)
	                    MAT_ETAS(3,3)=ETAS(1)
                        MAT_ETAS(3,6)=ETAS(2)
                        MAT_ETAS(3,9)=ETAS(3)
                        DEALLOCATE(VALUES_ETA)
                        DO K=1,3
			                DO L=1,3
			                    DO M=1,3
					                AUX1=0.D0
					                AUX1=KR(K,M)*DR(L)+KR(L,M)*DR(K)    
					                AUX1=AUX1-KR(K,L)*DR(M)
					                AUX1=AUX1*(1.D0-2.D0*Nu)
					                AUX1=AUX1+3.D0*DR(K)*DR(L)*DR(M)
					                D_AST((L+3*K-3),M)=(1.D0/R**2)*C2*AUX1			        
				                ENDDO
			                ENDDO
                        ENDDO
                        ! IN THIS CONTEXT, U_AST WAS USED JUST TO NOT CHANGE AFTER ENDIF
                        U_AST=MATMUL(MAT_ETAS,D_AST)
                    ENDIF
                ELSE
                    U_AST=0.D0
                ENDIF
                A(3*II-2,3*N_COLLOCPOINTS+I)=-U_AST(1,INT_DIR_CON_SUPPS(I))
                A(3*II-1,3*N_COLLOCPOINTS+I)=-U_AST(2,INT_DIR_CON_SUPPS(I))
                A(3*II,3*N_COLLOCPOINTS+I)=-U_AST(3,INT_DIR_CON_SUPPS(I))
            ENDDO
!           CALCULATING A_X AND F_X  
            CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM_NODE,NEN,COEFFICIENTS)
            CALL SHAPE_FUNCTIONS(2,QSI_PNT(1),QSI_PNT(2),NEN,VALUES,COEFFICIENTS,AUX_DBL,AUX_DBL,AUX_DBL,PHI,DPHI,D2PHI)
            DO JJ=1,NEN
                NUM_COLLOCPOINT=NODES_CONNECTIVITY(ELEM_NODE,JJ)
                IF(INT_DIR_CON_SUPPS(I).EQ.1)THEN
                    IF(B_CONDITIONS(3*NUM_COLLOCPOINT-2).EQ.0)THEN
                        F(3*N_COLLOCPOINTS+I)=-PHI(JJ)*U(3*NUM_COLLOCPOINT-2)+VALUES_CON_SUPPS(I)
                    ELSE
                        A(3*N_COLLOCPOINTS+I,3*NUM_COLLOCPOINT-2)=PHI(JJ)
                        F(3*N_COLLOCPOINTS+I)=VALUES_CON_SUPPS(I)
                    ENDIF
                ELSEIF(INT_DIR_CON_SUPPS(I).EQ.2)THEN
                    IF(B_CONDITIONS(3*NUM_COLLOCPOINT-1).EQ.0)THEN
                        F(3*N_COLLOCPOINTS+I)=-PHI(JJ)*U(3*NUM_COLLOCPOINT-1)+VALUES_CON_SUPPS(I)
                    ELSE
                        A(3*N_COLLOCPOINTS+I,3*NUM_COLLOCPOINT-1)=PHI(JJ)
                        F(3*N_COLLOCPOINTS+I)=VALUES_CON_SUPPS(I)
                    ENDIF
                ELSE
                    IF(B_CONDITIONS(3*NUM_COLLOCPOINT).EQ.0)THEN
                        F(3*N_COLLOCPOINTS+I)=-PHI(JJ)*U(3*NUM_COLLOCPOINT)+VALUES_CON_SUPPS(I)
                    ELSE
                        A(3*N_COLLOCPOINTS+I,3*NUM_COLLOCPOINT)=PHI(JJ)
                        F(3*N_COLLOCPOINTS+I)=VALUES_CON_SUPPS(I)
                    ENDIF
                ENDIF
            ENDDO
        !
            DEALLOCATE(VALUES,COEFFICIENTS,PHI,DPHI,D2PHI)    
        ENDDO
    ENDIF
    DEALLOCATE(INT_DIR_CON_SUPPS)
    END SUBROUTINE XBEM_SUPPORT_CHS