    SUBROUTINE XBEM_SUPPORT
!
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE XBEM_SUPPORT_VARIABLES
    USE XBEM_FORCE_VARIABLES
!   
    IMPLICIT NONE
!
    LOGICAL::SINGULAR_INT  
!
    INTEGER::I,J,K,L,M,II,JJ,KK,LL,NEN,ELEMS,NODES,NUM_COLLOCPOINT,ELEM_NODE,POS_COL,POS_LIN,N_GAUSS,&
    N_DIV,CONT
!
    REAL*8::DX,DY,DZ,DR(3),U_AST(3,3),D_AST(9,3),R,&
    Mu,Nu,C1,C2,AUX1,KR(3,3),PI,VALUES[ALLOCATABLE](:,:),&
    ETAS(3),MAT_ETAS(3,9),VJC(3),DVJC(3,2),JC,AUX_D(3,3),&
    QSI_PNT(2),ETA(3),VALUES_DIST_AREA(4,3),&
    COEFFICIENTS[ALLOCATABLE](:,:),AUX_DBL,VALUES_ETA[ALLOCATABLE](:,:),&
    PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),TOLDIST
    REAL*8::THETA,THETAI[ALLOCATABLE](:),THETAF[ALLOCATABLE](:),&
    RHO_BOUND,RHO_BOUND_AUX1,RHO_BOUND_AUX2(4),QSI1,QSI2,RHO,AUX_CTE1,&
    QSIW_GAUSS[ALLOCATABLE](:,:),N_P(3),TOL1,AUX_CTE2,U_AST_N(3,4)
    REAL*8::X,Y,Z,X_TEST,Y_TEST,Z_TEST,COEFFICIENTS_AREA[ALLOCATABLE](:,:),QSI_COLLOC(2)
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
    N_GAUSS=N_GAUSSPOINTS
    ALLOCATE(QSIW_GAUSS(N_GAUSS,2),COEFFICIENTS_AREA(4,4))
    CALL GAUSS_POINTS(N_GAUSS,QSIW_GAUSS)
!   
    IF (ALLOCATED(AS_COL)) DEALLOCATE(AS_COL,AS_LIN,XBEM_SUPP_REACTIONS)
    ALLOCATE(XBEM_SUPP_REACTIONS(NUM_CON_SUPPS+NODES_DIST_SUPPS))
    ALLOCATE(AS_COL(3*N_COLLOCPOINTS,NUM_CON_SUPPS+NODES_DIST_SUPPS),AS_LIN(NUM_CON_SUPPS+NODES_DIST_SUPPS,3*N_COLLOCPOINTS))
    AS_COL=0.D0
    AS_LIN=0.D0
!   CONCENTRATED SUPPORTS
!!$OMP PARALLEL PRIVATE(I,J,K,L,M,II,JJ,KK,NEN,DX,DY,DZ,DR,R,ELEM_NODE,QSI_PNT,VALUES,COEFFICIENTS,PHI,DPHI,D2PHI,ELEMS,Mu,Nu,C1,C2,AUX1,U_AST,NODES,VALUES_ETA,ETAS,VJC,DVJC,JC,MAT_ETAS,D_AST,NUM_COLLOCPOINT)
!!$OMP DO SCHEDULE(DYNAMIC)	
    DO I=1,NUM_CON_SUPPS
!
        ELEM_NODE=ELEM_CON_SUPP(I)
        QSI_PNT=QSI_CON_SUPP(I,:)
!
        NEN=ELEM_TYPE(ELEM_NODE)*ORDER_ELEM(ELEM_NODE)+(ELEM_TYPE(ELEM_NODE)-3)*(ORDER_ELEM(ELEM_NODE)-1)*POL_FAMILY(ELEM_NODE)
        ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2))
        DO J=1,NEN
            VALUES(J,:)=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM_NODE,J),:)
        ENDDO
!
        DO II=1,N_COLLOCPOINTS
!           FINDING ELEMENT THAT CONTAINS COLLOCATION POINT
!           ELEMS = ELEMENT CONTAINING COLLOCATION POINT  
            DO K=1,N_ELEM
		        NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			    DO KK=1,NEN
				    IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.II) THEN
					    ELEMS=K
				    ENDIF
			    ENDDO
            ENDDO
!           CHECKING IF SUPPORT IS AT SAME DOMAIN THAT COLLOCATION POINT            
            IF (ND(ELEM_NODE).EQ.ND(ELEMS)) THEN
                Mu=EMP(ND(ELEM_NODE),1)/(2*(1+EMP(ND(ELEM_NODE),2)))
	            Nu=EMP(ND(ELEM_NODE),2)
                C1=((1.D0)/(16.D0*PI*Mu*(1.D0-Nu)))
                C2=((1.D0)/(8.D0*PI*(1.D0-Nu)))
!               CALCULATING U*
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
		                VALUES_ETA(JJ,1)=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEMS,JJ),1)
		                VALUES_ETA(JJ,2)=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEMS,JJ),2)
		                VALUES_ETA(JJ,3)=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEMS,JJ),3)
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
            AS_COL(3*II-2,I)=-U_AST(1,INT_DIR_CON_SUPPS(I))
            AS_COL(3*II-1,I)=-U_AST(2,INT_DIR_CON_SUPPS(I))
            AS_COL(3*II,I)=-U_AST(3,INT_DIR_CON_SUPPS(I))
        ENDDO
!       CALCULATING A_X AND F_X  
        CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM_NODE,NEN,COEFFICIENTS)
        CALL SHAPE_FUNCTIONS(2,QSI_PNT(1),QSI_PNT(2),NEN,VALUES,COEFFICIENTS,AUX_DBL,AUX_DBL,AUX_DBL,PHI,DPHI,D2PHI)
        DO JJ=1,NEN
            NUM_COLLOCPOINT=COLLOCPOINTS_CONNECTIVITY(ELEM_NODE,JJ)
            AS_LIN(I,3*NUM_COLLOCPOINT-3+INT_DIR_CON_SUPPS(I))=PHI(JJ)
        ENDDO
    !
        DEALLOCATE(VALUES,COEFFICIENTS,PHI,DPHI,D2PHI)    
    ENDDO
!!$OMP END DO
!!$OMP END PARALLEL  
!   ********************************************
!               DISTRIBUTED SUPPORTS
!   ********************************************
    IF(.NOT. ALLOCATED(QSI_DIST)) THEN
        ALLOCATE(QSI_DIST(NUM_DIST_SUPPS,4,2))
!!   QSI_DIST STORES ADIMENSIONAL COORDINATES OF ENRICHED AREA
        QSI_DIST(:,1,1)=-1.D0
        QSI_DIST(:,1,2)=-1.D0
        QSI_DIST(:,2,1)=1.D0
        QSI_DIST(:,2,2)=-1.D0
        QSI_DIST(:,3,1)=1.D0
        QSI_DIST(:,3,2)=1.D0
        QSI_DIST(:,4,1)=-1.D0
        QSI_DIST(:,4,2)=1.D0
    ENDIF
!
    POS_COL=NUM_CON_SUPPS
    DO II=1,NUM_DIST_SUPPS_X
        !POS_COL=NUM_CON_SUPPS!+4*(II-1)
!       IDENTIFYING ELEMENT CONTAINING DISTRIBUTED LOAD AND CALCULATING ADIMENSIONAL 
!       PARAMETERS OF CONCENTRATED SUPPORT POINT
!       THE DISTRIBUTED FORCE MUST BE AT ONE DOMAIN ONLY
        ELEM_NODE=0
        CALL FIND_ELEMENT_WITH_PNT(COORDS_DIST_SUPPS_X(CONEC_DIST_SUPPS_X(II,1),:),ELEM_NODE,QSI_PNT)
!
        DO JJ=1,N_COLLOCPOINTS
            SINGULAR_INT=.FALSE.
!           FINDING ELEMENT THAT CONTAINS COLLOCATION POINT
            DO K=1,N_ELEM
		        NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			    DO KK=1,NEN
				    IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.JJ) THEN
					    ELEMS=K
                        EXIT
				    ENDIF
			    ENDDO
            ENDDO
!
            NEN=ELEM_TYPE(ELEMS)*ORDER_ELEM(ELEMS)+(ELEM_TYPE(ELEMS)-3)*(ORDER_ELEM(ELEMS)-1)*POL_FAMILY(ELEMS)
	        IF(ALLOCATED(VALUES)) DEALLOCATE(VALUES,PHI,DPHI)
            ALLOCATE(VALUES(NEN,3),PHI(NEN),DPHI(NEN,2))
!
            IF (ND(ELEM_NODE).EQ.ND(ELEMS)) THEN
!
!               VERIFYING IF COLLOCATION POINT IS INSIDE DISTRIBUTED AREA
                DO K=1,4
                    VALUES_DIST_AREA(K,1)=COORDS_DIST_SUPPS_X(CONEC_DIST_SUPPS_X(II,K),1)
                    VALUES_DIST_AREA(K,2)=COORDS_DIST_SUPPS_X(CONEC_DIST_SUPPS_X(II,K),2)
                    VALUES_DIST_AREA(K,3)=COORDS_DIST_SUPPS_X(CONEC_DIST_SUPPS_X(II,K),3)
                ENDDO
                CALL NEWTON_RAPHSON_ADIMENSIONAL_COORDINATES(4,VALUES_DIST_AREA,COORD_COLLOCPOINTS(JJ,:),QSI_COLLOC)
                CALL SHAPE_FUNCTIONS(1,QSI_COLLOC(1),QSI_COLLOC(2),4,VALUES_DIST_AREA,COEFFICIENTS_AREA,X_TEST,Y_TEST,Z_TEST,PHI,DPHI,D2PHI)
                DX=X_TEST-COORD_COLLOCPOINTS(JJ,1)
                DY=Y_TEST-COORD_COLLOCPOINTS(JJ,2)
                DZ=Z_TEST-COORD_COLLOCPOINTS(JJ,3)
                R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                IF(R.LT.TOLDIST)THEN
                    IF (QSI_COLLOC(1).LT.(1.D0-TOL1).AND.QSI_COLLOC(1).GT.(-1.D0+TOL1).AND.QSI_COLLOC(2).LT.(1.D0-TOL1).AND.QSI_COLLOC(2).GT.(-1.D0+TOL1)) SINGULAR_INT=.TRUE.
                ENDIF
!
                CALL SHAPE_FUNCTIONS_COEFFICIENTS_XBEM_DIST(II,COEFFICIENTS_AREA)
                U_AST_N=0.D0
                IF(SINGULAR_INT)THEN
                    ! CHANGING THE PRESCRIBED BOUNDARY CONDITION FROM TRACTION TO DISPLACEMENT
                    B_CONDITIONS(3*JJ-3+1)=0
                    CALL SHAPE_FUNCTIONS(2,QSI_COLLOC(1),QSI_COLLOC(2),4,VALUES_DIST_AREA,COEFFICIENTS_AREA,X_TEST,Y_TEST,Z_TEST,PHI,DPHI,D2PHI)
                    U(3*JJ-3+1)=0.D0
                    DO LL=1,4
                        U(3*JJ-3+1)=U(3*JJ-3+1)+VALUES_DIST_SUPPS_X(CONEC_DIST_SUPPS_X(II,LL))*PHI(LL)
                    ENDDO
                    ! POLAR INTEGRAL SUBDIVISION
                    N_DIV=4
                    IF (ALLOCATED(THETAI)) DEALLOCATE(THETAI,THETAF)
                    ALLOCATE(THETAI(N_DIV),THETAF(N_DIV))
                    THETAI(1)=-PI+DATAN((-1.D0-QSI_COLLOC(2))/(-1.D0-QSI_COLLOC(1)))
                    THETAF(1)=DATAN((-1.D0-QSI_COLLOC(2))/(1.D0-QSI_COLLOC(1)))
                    THETAI(2)=DATAN((-1.D0-QSI_COLLOC(2))/(1.D0-QSI_COLLOC(1)))
                    THETAF(2)=DATAN((1.D0-QSI_COLLOC(2))/(1.D0-QSI_COLLOC(1)))
                    THETAI(3)=DATAN((1.D0-QSI_COLLOC(2))/(1.D0-QSI_COLLOC(1)))
                    THETAF(3)=PI+DATAN((1.D0-QSI_COLLOC(2))/(-1.D0-QSI_COLLOC(1)))            
                    THETAI(4)=PI+DATAN((1.D0-QSI_COLLOC(2))/(-1.D0-QSI_COLLOC(1)))
                    THETAF(4)=PI+DATAN((-1.D0-QSI_COLLOC(2))/(-1.D0-QSI_COLLOC(1)))  
                    DO LL=1,N_DIV
                        DO I=1,N_GAUSSPOINTS        
                            THETA=(THETAF(LL)-THETAI(LL))/2*QSIW_GAUSS(I,1)+(THETAF(LL)+THETAI(LL))/2
!***************************************************************************************************************************  
!                           COMPUTING RHO_BOUND(THETA)     
!***************************************************************************************************************************    
                            !
                            RHO_BOUND_AUX2(1)=(1-QSI_COLLOC(1))/DCOS(THETA)
                            RHO_BOUND_AUX2(2)=(-1-QSI_COLLOC(1))/DCOS(THETA)
                            RHO_BOUND_AUX2(3)=(1-QSI_COLLOC(2))/DSIN(THETA)
                            RHO_BOUND_AUX2(4)=(-1-QSI_COLLOC(2))/DSIN(THETA)
    !                        
                            RHO_BOUND=DSQRT(8.D0)
                            DO K=1,4
                                IF(RHO_BOUND_AUX2(K).GT.1.D-10)THEN
                                    RHO_BOUND_AUX1=RHO_BOUND_AUX2(K)
                                    IF(RHO_BOUND_AUX1.LE.RHO_BOUND)THEN
                                        RHO_BOUND=RHO_BOUND_AUX1
                                    ENDIF
                                ENDIF
                            ENDDO                
    !***************************************************************************************************************************                        
                            DO J=1,N_GAUSSPOINTS
                                RHO=(RHO_BOUND/2.D0)*QSIW_GAUSS(J,1)+(RHO_BOUND/2.D0)            
		                        DR=0.d0
		                        ETA=0.D0
		                        PHI=0.D0
		                        QSI1=RHO*DCOS(THETA)+QSI_COLLOC(1)
		                        QSI2=RHO*DSIN(THETA)+QSI_COLLOC(2)
!
                                CALL SHAPE_FUNCTIONS(1,QSI1,QSI2,4,VALUES_DIST_AREA,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
!                           
		                        DX=X-COORD_COLLOCPOINTS(JJ,1)
		                        DY=Y-COORD_COLLOCPOINTS(JJ,2)
		                        DZ=Z-COORD_COLLOCPOINTS(JJ,3)
		                        R=DSQRT(DX*DX+DY*DY+DZ*DZ)
		                        DR(1)=DX/R
		                        DR(2)=DY/R
		                        DR(3)=DZ/R
!
		                        CALL NORMAL_OUTWARD(4,VALUES_DIST_AREA,QSI1,QSI2,ETA,VJC,DVJC,JC)
!
                                AUX_CTE1=(1.D0/R)*C1*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*RHO*JC*((THETAF(LL)-THETAI(LL))/2)*(RHO_BOUND/2.D0)
                                DO K=1,3
			                        DO L=1,3
				                        AUX1=(3.D0-4.D0*Nu)*KR(K,L)
				                        AUX1=AUX1+DR(K)*DR(L)
				                        U_AST(K,L)=AUX1*AUX_CTE1
			                        ENDDO
		                        ENDDO
!
                                CALL SHAPE_FUNCTIONS(2,QSI1,QSI2,4,VALUES,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                                DO K=1,3
                                    DO L=1,4
                                        U_AST_N(K,L)=U_AST_N(K,L)+U_AST(K,1)*PHI(L)
                                    ENDDO
                                ENDDO
                            ENDDO
                        ENDDO
                    ENDDO
                ELSE
!                   NON-SINGULAR INTEGRAL IN AREA    
                    Mu=EMP(ND(ELEM_NODE),1)/(2*(1+EMP(ND(ELEM_NODE),2)))
	                Nu=EMP(ND(ELEM_NODE),2)
                    C1=((1.D0)/(16.D0*PI*Mu*(1.D0-Nu)))
                    C2=((1.D0)/(8.D0*PI*(1.D0-Nu)))
                    DO I=1,N_GAUSS
                        DO J=1,N_GAUSS
                            CALL SHAPE_FUNCTIONS(1,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),NEN,VALUES_DIST_AREA,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                            CALL NORMAL_OUTWARD(4,VALUES_DIST_AREA,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),ETA,VJC,DVJC,JC)
		                    DX=X-COORD_COLLOCPOINTS(JJ,1)
		                    DY=Y-COORD_COLLOCPOINTS(JJ,2)
		                    DZ=Z-COORD_COLLOCPOINTS(JJ,3)
		                    R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                            IF(R.LT.TOLDIST)THEN
                                WRITE(*,*)"DISTRIBUTED SUPPORT AND COLLOCATION POINT AT SAME POSITION! CHANGE MESH!"
                                READ(*,*)
                            ENDIF
		                    DR(1)=DX/R
		                    DR(2)=DY/R
		                    DR(3)=DZ/R
                            AUX_CTE1=C1/R*JC*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)
                            IF (EQ_TYPE(JJ).EQ."S") THEN
                                DO K=1,3
                                    DO L=1,3
                                        AUX1=(3.0D0-4*Nu)*KR(K,L)
                                        U_AST(K,L)=AUX_CTE1*(AUX1+DR(K)*DR(L))
                                    ENDDO
                                ENDDO
!
                                CALL SHAPE_FUNCTIONS(2,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),4,VALUES,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                                DO K=1,3
                                    DO L=1,4
                                        U_AST_N(K,L)=U_AST_N(K,L)+U_AST(K,1)*PHI(L)
                                    ENDDO
                                ENDDO
                            ELSE
                                ! AINDA NÃO TESTEI
!                                !---------------------------------------------------------
!                                !       SOURCE POINT NORMAL OUTWARD MATRIX
!                                !---------------------------------------------------------
!!
                                NEN=ELEM_TYPE(ELEMS)*ORDER_ELEM(ELEMS)+(ELEM_TYPE(ELEMS)-3)*(ORDER_ELEM(ELEMS)-1)*POL_FAMILY(ELEMS)
            	                IF(ALLOCATED(VALUES)) DEALLOCATE(VALUES)
                                ALLOCATE(VALUES(NEN,3))
            	                NODES=0
                                !NODES = LOCAL NODE
                                DO KK=1,NEN
            		                IF(COLLOCPOINTS_CONNECTIVITY(ELEMS,KK).EQ.JJ) THEN
            			                NODES=KK
            		                ENDIF
                                ENDDO
                                !VALUES = COORDINATES OF ELEMENT NODES		
            	                DO KK=1,NEN
            		                VALUES(KK,1)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,KK),1)
            		                VALUES(KK,2)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,KK),2)
            		                VALUES(KK,3)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,KK),3)
            	                ENDDO
            
            	                CALL NORMAL_OUTWARD(NEN,VALUES,QSI(ELEMS,NODES,1),QSI(ELEMS,NODES,2),ETAS,VJC,DVJC,JC)
            
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
                                CALL NORMAL_OUTWARD(4,VALUES_DIST_AREA,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),ETA,VJC,DVJC,JC)
                                AUX_CTE2=C2/(R*R)*JC*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)
                                DO K=1,3
            			            DO L=1,3
            			                DO M=1,3
            					            AUX1=0.D0
            					            AUX1=KR(K,M)*DR(L)+KR(L,M)*DR(K)    
            					            AUX1=AUX1-KR(K,L)*DR(M)
            					            AUX1=AUX1*(1.D0-2.D0*Nu)
            					            AUX1=AUX1+3.D0*DR(K)*DR(L)*DR(M)
            					            D_AST((L+3*K-3),M)=AUX_CTE2*AUX1			        
            				            ENDDO
                                    ENDDO
                                ENDDO
                                AUX_D=MATMUL(MAT_ETAS,D_AST)
!                                
                                CALL SHAPE_FUNCTIONS(2,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),4,VALUES,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                                DO K=1,3
                                    DO L=1,4
                                        U_AST_N(K,L)=U_AST_N(K,L)+AUX_D(K,1)*PHI(L)
                                    ENDDO
                                ENDDO
                            ENDIF
                        ENDDO
                    ENDDO
                ENDIF
            ELSE
                U_AST_N=0.D0
            ENDIF
            DO L=1,4
                AS_COL(3*JJ-2:3*JJ,POS_COL+CONEC_DIST_SUPPS_X(II,L))=AS_COL(3*JJ-2:3*JJ,POS_COL+CONEC_DIST_SUPPS_X(II,L))-U_AST_N(:,L)
            ENDDO
        ENDDO
    ENDDO
!
    POS_COL=NUM_CON_SUPPS+NODES_DIST_SUPPS_X
    DO II=1,NUM_DIST_SUPPS_Y
!       IDENTIFYING ELEMENT CONTAINING DISTRIBUTED LOAD AND CALCULATING ADIMENSIONAL 
!       PARAMETERS OF CONCENTRATED SUPPORT POINT
!       THE DISTRIBUTED FORCE MUST BE AT ONE DOMAIN ONLY
        ELEM_NODE=0
        CALL FIND_ELEMENT_WITH_PNT(COORDS_DIST_SUPPS_Y(CONEC_DIST_SUPPS_Y(II,1),:),ELEM_NODE,QSI_PNT)
!
        DO JJ=1,N_COLLOCPOINTS
            SINGULAR_INT=.FALSE.
!           FINDING ELEMENT THAT CONTAINS COLLOCATION POINT
            DO K=1,N_ELEM
		        NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			    DO KK=1,NEN
				    IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.JJ) THEN
					    ELEMS=K
                        EXIT
				    ENDIF
			    ENDDO
            ENDDO
!
            NEN=ELEM_TYPE(ELEMS)*ORDER_ELEM(ELEMS)+(ELEM_TYPE(ELEMS)-3)*(ORDER_ELEM(ELEMS)-1)*POL_FAMILY(ELEMS)
	        IF(ALLOCATED(VALUES)) DEALLOCATE(VALUES,PHI,DPHI)
            ALLOCATE(VALUES(NEN,3),PHI(NEN),DPHI(NEN,2))
!
            IF (ND(ELEM_NODE).EQ.ND(ELEMS)) THEN
!
!               VERIFYING IF COLLOCATION POINT IS INSIDE DISTRIBUTED AREA
                DO K=1,4
                    VALUES_DIST_AREA(K,1)=COORDS_DIST_SUPPS_Y(CONEC_DIST_SUPPS_Y(II,K),1)
                    VALUES_DIST_AREA(K,2)=COORDS_DIST_SUPPS_Y(CONEC_DIST_SUPPS_Y(II,K),2)
                    VALUES_DIST_AREA(K,3)=COORDS_DIST_SUPPS_Y(CONEC_DIST_SUPPS_Y(II,K),3)
                ENDDO
                CALL NEWTON_RAPHSON_ADIMENSIONAL_COORDINATES(4,VALUES_DIST_AREA,COORD_COLLOCPOINTS(JJ,:),QSI_COLLOC)
                CALL SHAPE_FUNCTIONS(1,QSI_COLLOC(1),QSI_COLLOC(2),4,VALUES_DIST_AREA,COEFFICIENTS_AREA,X_TEST,Y_TEST,Z_TEST,PHI,DPHI,D2PHI)
                DX=X_TEST-COORD_COLLOCPOINTS(JJ,1)
                DY=Y_TEST-COORD_COLLOCPOINTS(JJ,2)
                DZ=Z_TEST-COORD_COLLOCPOINTS(JJ,3)
                R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                IF(R.LT.TOLDIST)THEN
                    IF (QSI_COLLOC(1).LT.(1.D0-TOL1).AND.QSI_COLLOC(1).GT.(-1.D0+TOL1).AND.QSI_COLLOC(2).LT.(1.D0-TOL1).AND.QSI_COLLOC(2).GT.(-1.D0+TOL1)) SINGULAR_INT=.TRUE.
                ENDIF
!
                CALL SHAPE_FUNCTIONS_COEFFICIENTS_XBEM_DIST(II,COEFFICIENTS_AREA)
                U_AST_N=0.D0
                IF(SINGULAR_INT)THEN
                    ! CHANGING THE PRESCRIBED BOUNDARY CONDITION FROM TRACTION TO DISPLACEMENT
                    B_CONDITIONS(3*JJ-3+2)=0
                    CALL SHAPE_FUNCTIONS(2,QSI_COLLOC(1),QSI_COLLOC(2),4,VALUES_DIST_AREA,COEFFICIENTS_AREA,X_TEST,Y_TEST,Z_TEST,PHI,DPHI,D2PHI)
                    U(3*JJ-3+2)=0.D0
                    DO LL=1,4
                        U(3*JJ-3+2)=U(3*JJ-3+2)+VALUES_DIST_SUPPS_Y(CONEC_DIST_SUPPS_Y(II,LL))*PHI(LL)
                    ENDDO
                    ! POLAR INTEGRAL SUBDIVISION
                    N_DIV=4
                    IF (ALLOCATED(THETAI)) DEALLOCATE(THETAI,THETAF)
                    ALLOCATE(THETAI(N_DIV),THETAF(N_DIV))
                    THETAI(1)=-PI+DATAN((-1.D0-QSI_COLLOC(2))/(-1.D0-QSI_COLLOC(1)))
                    THETAF(1)=DATAN((-1.D0-QSI_COLLOC(2))/(1.D0-QSI_COLLOC(1)))
                    THETAI(2)=DATAN((-1.D0-QSI_COLLOC(2))/(1.D0-QSI_COLLOC(1)))
                    THETAF(2)=DATAN((1.D0-QSI_COLLOC(2))/(1.D0-QSI_COLLOC(1)))
                    THETAI(3)=DATAN((1.D0-QSI_COLLOC(2))/(1.D0-QSI_COLLOC(1)))
                    THETAF(3)=PI+DATAN((1.D0-QSI_COLLOC(2))/(-1.D0-QSI_COLLOC(1)))            
                    THETAI(4)=PI+DATAN((1.D0-QSI_COLLOC(2))/(-1.D0-QSI_COLLOC(1)))
                    THETAF(4)=PI+DATAN((-1.D0-QSI_COLLOC(2))/(-1.D0-QSI_COLLOC(1)))  
                    DO LL=1,N_DIV
                        DO I=1,N_GAUSSPOINTS        
                            THETA=(THETAF(LL)-THETAI(LL))/2*QSIW_GAUSS(I,1)+(THETAF(LL)+THETAI(LL))/2
!***************************************************************************************************************************  
!                           COMPUTING RHO_BOUND(THETA)     
!***************************************************************************************************************************    
                            !
                            RHO_BOUND_AUX2(1)=(1-QSI_COLLOC(1))/DCOS(THETA)
                            RHO_BOUND_AUX2(2)=(-1-QSI_COLLOC(1))/DCOS(THETA)
                            RHO_BOUND_AUX2(3)=(1-QSI_COLLOC(2))/DSIN(THETA)
                            RHO_BOUND_AUX2(4)=(-1-QSI_COLLOC(2))/DSIN(THETA)
    !                        
                            RHO_BOUND=DSQRT(8.D0)
                            DO K=1,4
                                IF(RHO_BOUND_AUX2(K).GT.1.D-10)THEN
                                    RHO_BOUND_AUX1=RHO_BOUND_AUX2(K)
                                    IF(RHO_BOUND_AUX1.LE.RHO_BOUND)THEN
                                        RHO_BOUND=RHO_BOUND_AUX1
                                    ENDIF
                                ENDIF
                            ENDDO                
    !***************************************************************************************************************************                        
                            DO J=1,N_GAUSSPOINTS
                                RHO=(RHO_BOUND/2.D0)*QSIW_GAUSS(J,1)+(RHO_BOUND/2.D0)            
		                        DR=0.d0
		                        ETA=0.D0
		                        PHI=0.D0
		                        QSI1=RHO*DCOS(THETA)+QSI_COLLOC(1)
		                        QSI2=RHO*DSIN(THETA)+QSI_COLLOC(2)
!
                                CALL SHAPE_FUNCTIONS(1,QSI1,QSI2,4,VALUES_DIST_AREA,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
!                           
		                        DX=X-COORD_COLLOCPOINTS(JJ,1)
		                        DY=Y-COORD_COLLOCPOINTS(JJ,2)
		                        DZ=Z-COORD_COLLOCPOINTS(JJ,3)
		                        R=DSQRT(DX*DX+DY*DY+DZ*DZ)
		                        DR(1)=DX/R
		                        DR(2)=DY/R
		                        DR(3)=DZ/R
!
		                        CALL NORMAL_OUTWARD(4,VALUES_DIST_AREA,QSI1,QSI2,ETA,VJC,DVJC,JC)
!
                                AUX_CTE1=(1.D0/R)*C1*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*RHO*JC*((THETAF(LL)-THETAI(LL))/2)*(RHO_BOUND/2.D0)
                                DO K=1,3
			                        DO L=1,3
				                        AUX1=(3.D0-4.D0*Nu)*KR(K,L)
				                        AUX1=AUX1+DR(K)*DR(L)
				                        U_AST(K,L)=AUX1*AUX_CTE1
			                        ENDDO
		                        ENDDO
!
                                CALL SHAPE_FUNCTIONS(2,QSI1,QSI2,4,VALUES,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                                DO K=1,3
                                    DO L=1,4
                                        U_AST_N(K,L)=U_AST_N(K,L)+U_AST(K,2)*PHI(L)
                                    ENDDO
                                ENDDO
                            ENDDO
                        ENDDO
                    ENDDO
                ELSE
!                   NON-SINGULAR INTEGRAL IN AREA    
                    Mu=EMP(ND(ELEM_NODE),1)/(2*(1+EMP(ND(ELEM_NODE),2)))
	                Nu=EMP(ND(ELEM_NODE),2)
                    C1=((1.D0)/(16.D0*PI*Mu*(1.D0-Nu)))
                    C2=((1.D0)/(8.D0*PI*(1.D0-Nu)))
                    DO I=1,N_GAUSS
                        DO J=1,N_GAUSS
                            CALL SHAPE_FUNCTIONS(1,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),NEN,VALUES_DIST_AREA,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                            CALL NORMAL_OUTWARD(4,VALUES_DIST_AREA,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),ETA,VJC,DVJC,JC)
		                    DX=X-COORD_COLLOCPOINTS(JJ,1)
		                    DY=Y-COORD_COLLOCPOINTS(JJ,2)
		                    DZ=Z-COORD_COLLOCPOINTS(JJ,3)
		                    R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                            IF(R.LT.TOLDIST)THEN
                                WRITE(*,*)"DISTRIBUTED SUPPORT AND COLLOCATION POINT AT SAME POSITION! CHANGE MESH!"
                                READ(*,*)
                            ENDIF
		                    DR(1)=DX/R
		                    DR(2)=DY/R
		                    DR(3)=DZ/R
                            AUX_CTE1=C1/R*JC*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)
                            IF (EQ_TYPE(JJ).EQ."S") THEN
                                DO K=1,3
                                    DO L=1,3
                                        AUX1=(3.0D0-4*Nu)*KR(K,L)
                                        U_AST(K,L)=AUX_CTE1*(AUX1+DR(K)*DR(L))
                                    ENDDO
                                ENDDO
!
                                CALL SHAPE_FUNCTIONS(2,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),4,VALUES,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                                DO K=1,3
                                    DO L=1,4
                                        U_AST_N(K,L)=U_AST_N(K,L)+U_AST(K,2)*PHI(L)
                                    ENDDO
                                ENDDO
                            ELSE
                                ! AINDA NÃO TESTEI
!                                !---------------------------------------------------------
!                                !       SOURCE POINT NORMAL OUTWARD MATRIX
!                                !---------------------------------------------------------
!!
                                NEN=ELEM_TYPE(ELEMS)*ORDER_ELEM(ELEMS)+(ELEM_TYPE(ELEMS)-3)*(ORDER_ELEM(ELEMS)-1)*POL_FAMILY(ELEMS)
            	                IF(ALLOCATED(VALUES)) DEALLOCATE(VALUES)
                                ALLOCATE(VALUES(NEN,3))
            	                NODES=0
                                !NODES = LOCAL NODE
                                DO KK=1,NEN
            		                IF(COLLOCPOINTS_CONNECTIVITY(ELEMS,KK).EQ.JJ) THEN
            			                NODES=KK
            		                ENDIF
                                ENDDO
                                !VALUES = COORDINATES OF ELEMENT NODES		
            	                DO KK=1,NEN
            		                VALUES(KK,1)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,KK),1)
            		                VALUES(KK,2)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,KK),2)
            		                VALUES(KK,3)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,KK),3)
            	                ENDDO
            
            	                CALL NORMAL_OUTWARD(NEN,VALUES,QSI(ELEMS,NODES,1),QSI(ELEMS,NODES,2),ETAS,VJC,DVJC,JC)
            
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
                                CALL NORMAL_OUTWARD(4,VALUES_DIST_AREA,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),ETA,VJC,DVJC,JC)
                                AUX_CTE2=C2/(R*R)*JC*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)
                                DO K=1,3
            			            DO L=1,3
            			                DO M=1,3
            					            AUX1=0.D0
            					            AUX1=KR(K,M)*DR(L)+KR(L,M)*DR(K)    
            					            AUX1=AUX1-KR(K,L)*DR(M)
            					            AUX1=AUX1*(1.D0-2.D0*Nu)
            					            AUX1=AUX1+3.D0*DR(K)*DR(L)*DR(M)
            					            D_AST((L+3*K-3),M)=AUX_CTE2*AUX1			        
            				            ENDDO
                                    ENDDO
                                ENDDO
                                AUX_D=MATMUL(MAT_ETAS,D_AST)
!                                
                                CALL SHAPE_FUNCTIONS(2,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),4,VALUES,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                                DO K=1,3
                                    DO L=1,4
                                        U_AST_N(K,L)=U_AST_N(K,L)+AUX_D(K,2)*PHI(L)
                                    ENDDO
                                ENDDO
                            ENDIF
                        ENDDO
                    ENDDO
                ENDIF
            ELSE
                U_AST_N=0.D0
            ENDIF
            DO L=1,4
                AS_COL(3*JJ-2:3*JJ,POS_COL+CONEC_DIST_SUPPS_Y(II,L))=AS_COL(3*JJ-2:3*JJ,POS_COL+CONEC_DIST_SUPPS_Y(II,L))-U_AST_N(:,L)
            ENDDO
        ENDDO
    ENDDO
!
    POS_COL=NUM_CON_SUPPS+NODES_DIST_SUPPS_X+NODES_DIST_SUPPS_Y
    DO II=1,NUM_DIST_SUPPS_Z
!       IDENTIFYING ELEMENT CONTAINING DISTRIBUTED LOAD AND CALCULATING ADIMENSIONAL 
!       PARAMETERS OF CONCENTRATED SUPPORT POINT
!       THE DISTRIBUTED FORCE MUST BE AT ONE DOMAIN ONLY
        ELEM_NODE=0
        CALL FIND_ELEMENT_WITH_PNT(COORDS_DIST_SUPPS_Z(CONEC_DIST_SUPPS_Z(II,1),:),ELEM_NODE,QSI_PNT)
!
        DO JJ=1,N_COLLOCPOINTS
            SINGULAR_INT=.FALSE.
!           FINDING ELEMENT THAT CONTAINS COLLOCATION POINT
            DO K=1,N_ELEM
		        NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			    DO KK=1,NEN
				    IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.JJ) THEN
					    ELEMS=K
                        EXIT
				    ENDIF
			    ENDDO
            ENDDO
!
            NEN=ELEM_TYPE(ELEMS)*ORDER_ELEM(ELEMS)+(ELEM_TYPE(ELEMS)-3)*(ORDER_ELEM(ELEMS)-1)*POL_FAMILY(ELEMS)
	        IF(ALLOCATED(VALUES)) DEALLOCATE(VALUES,PHI,DPHI)
            ALLOCATE(VALUES(NEN,3),PHI(NEN),DPHI(NEN,2))
!
            IF (ND(ELEM_NODE).EQ.ND(ELEMS)) THEN
!
!               VERIFYING IF COLLOCATION POINT IS INSIDE DISTRIBUTED AREA
                DO K=1,4
                    VALUES_DIST_AREA(K,1)=COORDS_DIST_SUPPS_Z(CONEC_DIST_SUPPS_Z(II,K),1)
                    VALUES_DIST_AREA(K,2)=COORDS_DIST_SUPPS_Z(CONEC_DIST_SUPPS_Z(II,K),2)
                    VALUES_DIST_AREA(K,3)=COORDS_DIST_SUPPS_Z(CONEC_DIST_SUPPS_Z(II,K),3)
                ENDDO
                CALL NEWTON_RAPHSON_ADIMENSIONAL_COORDINATES(4,VALUES_DIST_AREA,COORD_COLLOCPOINTS(JJ,:),QSI_COLLOC)
                CALL SHAPE_FUNCTIONS(1,QSI_COLLOC(1),QSI_COLLOC(2),4,VALUES_DIST_AREA,COEFFICIENTS_AREA,X_TEST,Y_TEST,Z_TEST,PHI,DPHI,D2PHI)
                DX=X_TEST-COORD_COLLOCPOINTS(JJ,1)
                DY=Y_TEST-COORD_COLLOCPOINTS(JJ,2)
                DZ=Z_TEST-COORD_COLLOCPOINTS(JJ,3)
                R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                IF(R.LT.TOLDIST)THEN
                    IF (QSI_COLLOC(1).LT.(1.D0-TOL1).AND.QSI_COLLOC(1).GT.(-1.D0+TOL1).AND.QSI_COLLOC(2).LT.(1.D0-TOL1).AND.QSI_COLLOC(2).GT.(-1.D0+TOL1)) SINGULAR_INT=.TRUE.
                ENDIF
!
                CALL SHAPE_FUNCTIONS_COEFFICIENTS_XBEM_DIST(II,COEFFICIENTS_AREA)
                U_AST_N=0.D0
                IF(SINGULAR_INT)THEN
                    ! CHANGING THE PRESCRIBED BOUNDARY CONDITION FROM TRACTION TO DISPLACEMENT
                    B_CONDITIONS(3*JJ-3+3)=0
                    CALL SHAPE_FUNCTIONS(2,QSI_COLLOC(1),QSI_COLLOC(2),4,VALUES_DIST_AREA,COEFFICIENTS_AREA,X_TEST,Y_TEST,Z_TEST,PHI,DPHI,D2PHI)
                    U(3*JJ-3+3)=0.D0
                    DO LL=1,4
                        U(3*JJ-3+3)=U(3*JJ-3+3)+VALUES_DIST_SUPPS_Z(CONEC_DIST_SUPPS_Z(II,LL))*PHI(LL)
                    ENDDO
                    ! POLAR INTEGRAL SUBDIVISION
                    N_DIV=4
                    IF (ALLOCATED(THETAI)) DEALLOCATE(THETAI,THETAF)
                    ALLOCATE(THETAI(N_DIV),THETAF(N_DIV))
                    THETAI(1)=-PI+DATAN((-1.D0-QSI_COLLOC(2))/(-1.D0-QSI_COLLOC(1)))
                    THETAF(1)=DATAN((-1.D0-QSI_COLLOC(2))/(1.D0-QSI_COLLOC(1)))
                    THETAI(2)=DATAN((-1.D0-QSI_COLLOC(2))/(1.D0-QSI_COLLOC(1)))
                    THETAF(2)=DATAN((1.D0-QSI_COLLOC(2))/(1.D0-QSI_COLLOC(1)))
                    THETAI(3)=DATAN((1.D0-QSI_COLLOC(2))/(1.D0-QSI_COLLOC(1)))
                    THETAF(3)=PI+DATAN((1.D0-QSI_COLLOC(2))/(-1.D0-QSI_COLLOC(1)))            
                    THETAI(4)=PI+DATAN((1.D0-QSI_COLLOC(2))/(-1.D0-QSI_COLLOC(1)))
                    THETAF(4)=PI+DATAN((-1.D0-QSI_COLLOC(2))/(-1.D0-QSI_COLLOC(1)))  
                    DO LL=1,N_DIV
                        DO I=1,N_GAUSSPOINTS        
                            THETA=(THETAF(LL)-THETAI(LL))/2*QSIW_GAUSS(I,1)+(THETAF(LL)+THETAI(LL))/2
!***************************************************************************************************************************  
!                           COMPUTING RHO_BOUND(THETA)     
!***************************************************************************************************************************    
                            !
                            RHO_BOUND_AUX2(1)=(1-QSI_COLLOC(1))/DCOS(THETA)
                            RHO_BOUND_AUX2(2)=(-1-QSI_COLLOC(1))/DCOS(THETA)
                            RHO_BOUND_AUX2(3)=(1-QSI_COLLOC(2))/DSIN(THETA)
                            RHO_BOUND_AUX2(4)=(-1-QSI_COLLOC(2))/DSIN(THETA)
    !                        
                            RHO_BOUND=DSQRT(8.D0)
                            DO K=1,4
                                IF(RHO_BOUND_AUX2(K).GT.1.D-10)THEN
                                    RHO_BOUND_AUX1=RHO_BOUND_AUX2(K)
                                    IF(RHO_BOUND_AUX1.LE.RHO_BOUND)THEN
                                        RHO_BOUND=RHO_BOUND_AUX1
                                    ENDIF
                                ENDIF
                            ENDDO                
    !***************************************************************************************************************************                        
                            DO J=1,N_GAUSSPOINTS
                                RHO=(RHO_BOUND/2.D0)*QSIW_GAUSS(J,1)+(RHO_BOUND/2.D0)            
		                        DR=0.d0
		                        ETA=0.D0
		                        PHI=0.D0
		                        QSI1=RHO*DCOS(THETA)+QSI_COLLOC(1)
		                        QSI2=RHO*DSIN(THETA)+QSI_COLLOC(2)
!
                                CALL SHAPE_FUNCTIONS(1,QSI1,QSI2,4,VALUES_DIST_AREA,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
!                           
		                        DX=X-COORD_COLLOCPOINTS(JJ,1)
		                        DY=Y-COORD_COLLOCPOINTS(JJ,2)
		                        DZ=Z-COORD_COLLOCPOINTS(JJ,3)
		                        R=DSQRT(DX*DX+DY*DY+DZ*DZ)
		                        DR(1)=DX/R
		                        DR(2)=DY/R
		                        DR(3)=DZ/R
!
		                        CALL NORMAL_OUTWARD(4,VALUES_DIST_AREA,QSI1,QSI2,ETA,VJC,DVJC,JC)
!
                                AUX_CTE1=(1.D0/R)*C1*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*RHO*JC*((THETAF(LL)-THETAI(LL))/2)*(RHO_BOUND/2.D0)
                                DO K=1,3
			                        DO L=1,3
				                        AUX1=(3.D0-4.D0*Nu)*KR(K,L)
				                        AUX1=AUX1+DR(K)*DR(L)
				                        U_AST(K,L)=AUX1*AUX_CTE1
			                        ENDDO
		                        ENDDO
!
                                CALL SHAPE_FUNCTIONS(2,QSI1,QSI2,4,VALUES,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                                DO K=1,3
                                    DO L=1,4
                                        U_AST_N(K,L)=U_AST_N(K,L)+U_AST(K,3)*PHI(L)
                                    ENDDO
                                ENDDO
                            ENDDO
                        ENDDO
                    ENDDO
                ELSE
!                   NON-SINGULAR INTEGRAL IN AREA    
                    Mu=EMP(ND(ELEM_NODE),1)/(2*(1+EMP(ND(ELEM_NODE),2)))
	                Nu=EMP(ND(ELEM_NODE),2)
                    C1=((1.D0)/(16.D0*PI*Mu*(1.D0-Nu)))
                    C2=((1.D0)/(8.D0*PI*(1.D0-Nu)))
                    DO I=1,N_GAUSS
                        DO J=1,N_GAUSS
                            CALL SHAPE_FUNCTIONS(1,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),NEN,VALUES_DIST_AREA,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                            CALL NORMAL_OUTWARD(4,VALUES_DIST_AREA,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),ETA,VJC,DVJC,JC)
		                    DX=X-COORD_COLLOCPOINTS(JJ,1)
		                    DY=Y-COORD_COLLOCPOINTS(JJ,2)
		                    DZ=Z-COORD_COLLOCPOINTS(JJ,3)
		                    R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                            IF(R.LT.TOLDIST)THEN
                                WRITE(*,*)"DISTRIBUTED SUPPORT AND COLLOCATION POINT AT SAME POSITION! CHANGE MESH!"
                                READ(*,*)
                            ENDIF
		                    DR(1)=DX/R
		                    DR(2)=DY/R
		                    DR(3)=DZ/R
                            AUX_CTE1=C1/R*JC*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)
                            IF (EQ_TYPE(JJ).EQ."S") THEN
                                DO K=1,3
                                    DO L=1,3
                                        AUX1=(3.0D0-4*Nu)*KR(K,L)
                                        U_AST(K,L)=AUX_CTE1*(AUX1+DR(K)*DR(L))
                                    ENDDO
                                ENDDO
!
                                CALL SHAPE_FUNCTIONS(2,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),4,VALUES,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                                DO K=1,3
                                    DO L=1,4
                                        U_AST_N(K,L)=U_AST_N(K,L)+U_AST(K,3)*PHI(L)
                                    ENDDO
                                ENDDO
                            ELSE
                                ! AINDA NÃO TESTEI
!                                !---------------------------------------------------------
!                                !       SOURCE POINT NORMAL OUTWARD MATRIX
!                                !---------------------------------------------------------
!!
                                NEN=ELEM_TYPE(ELEMS)*ORDER_ELEM(ELEMS)+(ELEM_TYPE(ELEMS)-3)*(ORDER_ELEM(ELEMS)-1)*POL_FAMILY(ELEMS)
            	                IF(ALLOCATED(VALUES)) DEALLOCATE(VALUES)
                                ALLOCATE(VALUES(NEN,3))
            	                NODES=0
                                !NODES = LOCAL NODE
                                DO KK=1,NEN
            		                IF(COLLOCPOINTS_CONNECTIVITY(ELEMS,KK).EQ.JJ) THEN
            			                NODES=KK
            		                ENDIF
                                ENDDO
                                !VALUES = COORDINATES OF ELEMENT NODES		
            	                DO KK=1,NEN
            		                VALUES(KK,1)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,KK),1)
            		                VALUES(KK,2)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,KK),2)
            		                VALUES(KK,3)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,KK),3)
            	                ENDDO
            
            	                CALL NORMAL_OUTWARD(NEN,VALUES,QSI(ELEMS,NODES,1),QSI(ELEMS,NODES,2),ETAS,VJC,DVJC,JC)
            
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
                                CALL NORMAL_OUTWARD(4,VALUES_DIST_AREA,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),ETA,VJC,DVJC,JC)
                                AUX_CTE2=C2/(R*R)*JC*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)
                                DO K=1,3
            			            DO L=1,3
            			                DO M=1,3
            					            AUX1=0.D0
            					            AUX1=KR(K,M)*DR(L)+KR(L,M)*DR(K)    
            					            AUX1=AUX1-KR(K,L)*DR(M)
            					            AUX1=AUX1*(1.D0-2.D0*Nu)
            					            AUX1=AUX1+3.D0*DR(K)*DR(L)*DR(M)
            					            D_AST((L+3*K-3),M)=AUX_CTE2*AUX1			        
            				            ENDDO
                                    ENDDO
                                ENDDO
                                AUX_D=MATMUL(MAT_ETAS,D_AST)
!                                
                                CALL SHAPE_FUNCTIONS(2,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),4,VALUES,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                                DO K=1,3
                                    DO L=1,4
                                        U_AST_N(K,L)=U_AST_N(K,L)+AUX_D(K,3)*PHI(L)
                                    ENDDO
                                ENDDO
                            ENDIF
                        ENDDO
                    ENDDO
                ENDIF
            ELSE
                U_AST_N=0.D0
            ENDIF
            DO L=1,4
                AS_COL(3*JJ-2:3*JJ,POS_COL+CONEC_DIST_SUPPS_Z(II,L))=AS_COL(3*JJ-2:3*JJ,POS_COL+CONEC_DIST_SUPPS_Z(II,L))-U_AST_N(:,L)
            ENDDO
        ENDDO
    ENDDO
!
    POS_LIN=NUM_CON_SUPPS
    ! COMPATIBILITY: CALCULATING A_X AND F_X  
    DO II=1,NODES_DIST_SUPPS_X
        ELEM_NODE=ELEM_DIST_SUPP(II)
        QSI_PNT=QSI_DIST_SUPP(II,:)
!
        NEN=ELEM_TYPE(ELEM_NODE)*ORDER_ELEM(ELEM_NODE)+(ELEM_TYPE(ELEM_NODE)-3)*(ORDER_ELEM(ELEM_NODE)-1)*POL_FAMILY(ELEM_NODE)
        IF(ALLOCATED(VALUES)) DEALLOCATE(VALUES)
        IF(ALLOCATED(COEFFICIENTS)) DEALLOCATE(COEFFICIENTS)
        IF(ALLOCATED(PHI)) DEALLOCATE(PHI)
        IF(ALLOCATED(DPHI)) DEALLOCATE(DPHI)
        IF(ALLOCATED(D2PHI)) DEALLOCATE(D2PHI)
        ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2))
        DO J=1,NEN
            VALUES(J,:)=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM_NODE,J),:)
        ENDDO
!
        CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM_NODE,NEN,COEFFICIENTS)
        CALL SHAPE_FUNCTIONS(2,QSI_PNT(1),QSI_PNT(2),NEN,VALUES,COEFFICIENTS,AUX_DBL,AUX_DBL,AUX_DBL,PHI,DPHI,D2PHI)
        DO JJ=1,NEN
            NUM_COLLOCPOINT=COLLOCPOINTS_CONNECTIVITY(ELEM_NODE,JJ)
            AS_LIN(POS_LIN+II,3*NUM_COLLOCPOINT-3+1)=PHI(JJ)
        ENDDO
!
        DEALLOCATE(VALUES,COEFFICIENTS,PHI,DPHI,D2PHI)   
    ENDDO
!
    POS_LIN=NUM_CON_SUPPS+NODES_DIST_SUPPS_X
    CONT=NODES_DIST_SUPPS_X
    DO II=1,NODES_DIST_SUPPS_Y
        ELEM_NODE=ELEM_DIST_SUPP(CONT+II)
        QSI_PNT=QSI_DIST_SUPP(CONT+II,:)
!
        NEN=ELEM_TYPE(ELEM_NODE)*ORDER_ELEM(ELEM_NODE)+(ELEM_TYPE(ELEM_NODE)-3)*(ORDER_ELEM(ELEM_NODE)-1)*POL_FAMILY(ELEM_NODE)
        IF(ALLOCATED(VALUES)) DEALLOCATE(VALUES)
        IF(ALLOCATED(COEFFICIENTS)) DEALLOCATE(COEFFICIENTS)
        IF(ALLOCATED(PHI)) DEALLOCATE(PHI)
        IF(ALLOCATED(DPHI)) DEALLOCATE(DPHI)
        IF(ALLOCATED(D2PHI)) DEALLOCATE(D2PHI)
        ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2))
        DO J=1,NEN
            VALUES(J,:)=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM_NODE,J),:)
        ENDDO
!
        CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM_NODE,NEN,COEFFICIENTS)
        CALL SHAPE_FUNCTIONS(2,QSI_PNT(1),QSI_PNT(2),NEN,VALUES,COEFFICIENTS,AUX_DBL,AUX_DBL,AUX_DBL,PHI,DPHI,D2PHI)
        DO JJ=1,NEN
            NUM_COLLOCPOINT=COLLOCPOINTS_CONNECTIVITY(ELEM_NODE,JJ)
            AS_LIN(POS_LIN+II,3*NUM_COLLOCPOINT-3+2)=PHI(JJ)
        ENDDO
!
        DEALLOCATE(VALUES,COEFFICIENTS,PHI,DPHI,D2PHI)   
    ENDDO

    POS_LIN=NUM_CON_SUPPS+NODES_DIST_SUPPS_X+NODES_DIST_SUPPS_Y
    CONT=NODES_DIST_SUPPS_X+NODES_DIST_SUPPS_Y
    DO II=1,NODES_DIST_SUPPS_Z
        ELEM_NODE=ELEM_DIST_SUPP(CONT+II)
        QSI_PNT=QSI_DIST_SUPP(CONT+II,:)
!
        NEN=ELEM_TYPE(ELEM_NODE)*ORDER_ELEM(ELEM_NODE)+(ELEM_TYPE(ELEM_NODE)-3)*(ORDER_ELEM(ELEM_NODE)-1)*POL_FAMILY(ELEM_NODE)
        IF(ALLOCATED(VALUES)) DEALLOCATE(VALUES)
        IF(ALLOCATED(COEFFICIENTS)) DEALLOCATE(COEFFICIENTS)
        IF(ALLOCATED(PHI)) DEALLOCATE(PHI)
        IF(ALLOCATED(DPHI)) DEALLOCATE(DPHI)
        IF(ALLOCATED(D2PHI)) DEALLOCATE(D2PHI)
        ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2))
        DO J=1,NEN
            VALUES(J,:)=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM_NODE,J),:)
        ENDDO
!
        CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM_NODE,NEN,COEFFICIENTS)
        CALL SHAPE_FUNCTIONS(2,QSI_PNT(1),QSI_PNT(2),NEN,VALUES,COEFFICIENTS,AUX_DBL,AUX_DBL,AUX_DBL,PHI,DPHI,D2PHI)
        DO JJ=1,NEN
            NUM_COLLOCPOINT=COLLOCPOINTS_CONNECTIVITY(ELEM_NODE,JJ)
            AS_LIN(POS_LIN+II,3*NUM_COLLOCPOINT-3+3)=PHI(JJ)
        ENDDO
!
        DEALLOCATE(VALUES,COEFFICIENTS,PHI,DPHI,D2PHI)   
    ENDDO
!   ORGANIZING VALUES_DIST_SUPPS
!    OPEN(2,file='AS_COL.txt',status='unknown')
!    DO I=1,SIZE(AS_COL,1)
!        WRITE(2,'(10000F)')(AS_COL(I,K),K=1,SIZE(AS_COL,2))
!    ENDDO
!    CLOSE(2)
!!!
!    OPEN(2,file='AS_LIN.txt',status='unknown')
!    DO I=1,SIZE(AS_LIN,1)
!        WRITE(2,'(10000F)')(AS_LIN(I,K),K=1,SIZE(AS_LIN,2))
!    ENDDO
!    CLOSE(2)
!
    END SUBROUTINE XBEM_SUPPORT