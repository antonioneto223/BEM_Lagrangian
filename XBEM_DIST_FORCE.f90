    SUBROUTINE XBEM_DIST_FORCE
!
!
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE XBEM_FORCE_VARIABLES
!    
    IMPLICIT NONE
!    
    LOGICAL::SINGULAR_INT   
!    
    INTEGER::I,J,K,L,M,II,JJ,KK,LL,NEN,ELEMS,NODES,ELEM_NODE,&
    POS_COLLOCPOINT,N_GAUSS,N_DIV
!
    REAL*8::DX,DY,DZ,DR(3),U_AST(3,3), D_AST(9,3),R,&
    Mu,Nu,C1,C2,AUX1,KR(3,3),PI,&
    ETA(3),ETAS(3),MAT_ETAS(3,9),VJC(3),DVJC(3,2),JC,AUX_D(3,3),&
    QSI_PNT(2),TOLDIST
    REAL*8::THETA,THETAI[ALLOCATABLE](:),THETAF[ALLOCATABLE](:),&
    RHO_BOUND,RHO_BOUND_AUX1,RHO_BOUND_AUX2(4),QSI1,QSI2,RHO,AUX_CTE1,&
    QSIW_GAUSS[ALLOCATABLE](:,:),N_P(3),TOL1,AUX_CTE2
    REAL*8::X,Y,Z,X_TEST,Y_TEST,Z_TEST,COEFFICIENTS_AREA[ALLOCATABLE](:,:),VALUES[ALLOCATABLE](:,:),&
    PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),QSI_COLLOC(2)
! 
    N_GAUSS=N_GAUSSPOINTS
    !N_GAUSS=2
    ALLOCATE(QSIW_GAUSS(N_GAUSS,2),COEFFICIENTS_AREA(4,4))
    CALL GAUSS_POINTS(N_GAUSS,QSIW_GAUSS)
    PI=DACOS(-1.D0)
    !F_ENRICHED=0.0D0
    TOLDIST=1.0E-5
    TOL1=1.D-3
!
    KR(1,1)=1.D0
	KR(1,2)=0.D0
    KR(1,3)=0.D0
	KR(2,1)=0.D0
	KR(2,2)=1.D0
	KR(2,3)=0.D0
	KR(3,1)=0.D0
	KR(3,2)=0.D0
	KR(3,3)=1.D0
!
    IF(ALLOCATED(QSI_DIST)) DEALLOCATE(QSI_DIST)
    ALLOCATE(QSI_DIST(NUM_DIST_LOADS,4,2))
!   QSI_DIST STORES ADIMENSIONAL COORDINATES OF ENRICHED AREA
    QSI_DIST(:,1,1)=-1.D0
    QSI_DIST(:,1,2)=-1.D0
    QSI_DIST(:,2,1)=1.D0
    QSI_DIST(:,2,2)=-1.D0
    QSI_DIST(:,3,1)=1.D0
    QSI_DIST(:,3,2)=1.D0
    QSI_DIST(:,4,1)=-1.D0
    QSI_DIST(:,4,2)=1.D0
!
    DO II=1,NUM_DIST_LOADS
!       IDENTIFYING ELEMENT CONTAINING DISTRIBUTED LOAD AND CALCULATING ADIMENSIONAL PARAMETERS OF CONCENTRATED SUPPORT POINT
!       THE DISTRIBUTED FORCE MUST BE AT ONE DOMAIN ONLY
        ELEM_NODE=0
        CALL FIND_ELEMENT_WITH_PNT(COORDS_DIST_LOADS(II,1,:),ELEM_NODE,QSI_PNT)
!
        DO JJ=1,N_COLLOCPOINTS
            F_ENRICHED_LOCAL=0.0D0
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
                CALL NEWTON_RAPHSON_ADIMENSIONAL_COORDINATES(4,COORDS_DIST_LOADS(II,:,:),COORD_COLLOCPOINTS(JJ,:),QSI_COLLOC)
                CALL SHAPE_FUNCTIONS(1,QSI_COLLOC(1),QSI_COLLOC(2),4,COORDS_DIST_LOADS(II,:,:),COEFFICIENTS_AREA,X_TEST,Y_TEST,Z_TEST,PHI,DPHI,D2PHI)
                DX=X_TEST-COORD_COLLOCPOINTS(JJ,1)
                DY=Y_TEST-COORD_COLLOCPOINTS(JJ,2)
                DZ=Z_TEST-COORD_COLLOCPOINTS(JJ,3)
                R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                IF(R.LT.TOLDIST)THEN
                    IF (QSI_COLLOC(1).LT.(1.D0-TOL1).AND.QSI_COLLOC(1).GT.(-1.D0+TOL1).AND.QSI_COLLOC(2).LT.(1.D0-TOL1).AND.QSI_COLLOC(2).GT.(-1.D0+TOL1)) SINGULAR_INT=.TRUE.
                ENDIF
!
                CALL SHAPE_FUNCTIONS_COEFFICIENTS_XBEM_DIST(II,COEFFICIENTS_AREA)
                IF(SINGULAR_INT)THEN
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
                                CALL SHAPE_FUNCTIONS(1,QSI1,QSI2,4,COORDS_DIST_LOADS(II,:,:),COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
    !                           
		                        DX=X-COORD_COLLOCPOINTS(JJ,1)
		                        DY=Y-COORD_COLLOCPOINTS(JJ,2)
		                        DZ=Z-COORD_COLLOCPOINTS(JJ,3)
		                        R=DSQRT(DX*DX+DY*DY+DZ*DZ)
		                        DR(1)=DX/R
		                        DR(2)=DY/R
		                        DR(3)=DZ/R
!
		                        CALL NORMAL_OUTWARD(4,COORDS_DIST_LOADS(II,:,:),QSI1,QSI2,ETA,VJC,DVJC,JC)
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
!
                                N_P=0.D0
                                DO K=1,4
                                    N_P=N_P+PHI(K)*VALUES_DIST_LOADS(II,:,K)
                                ENDDO
                                F_ENRICHED_LOCAL=F_ENRICHED_LOCAL+MATMUL(U_AST,N_P)
!                                
                            ENDDO
                        ENDDO
                    ENDDO
                ELSE
!                   REGULAR INTEGRAL IN AREA    
                    DO I=1,N_GAUSS
                        DO J=1,N_GAUSS
                            CALL SHAPE_FUNCTIONS(1,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),NEN,COORDS_DIST_LOADS(II,:,:),COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                            CALL NORMAL_OUTWARD(4,COORDS_DIST_LOADS(II,:,:),QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),ETA,VJC,DVJC,JC)
!                        
                            Mu=EMP(ND(ELEM_NODE),1)/(2*(1+EMP(ND(ELEM_NODE),2)))
	                        Nu=EMP(ND(ELEM_NODE),2)
                            C1=((1.D0)/(16.D0*PI*Mu*(1.D0-Nu)))
                            C2=((1.D0)/(8.D0*PI*(1.D0-Nu)))
		                    DX=X-COORD_COLLOCPOINTS(JJ,1)
		                    DY=Y-COORD_COLLOCPOINTS(JJ,2)
		                    DZ=Z-COORD_COLLOCPOINTS(JJ,3)
		                    R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                            IF(R.LT.TOLDIST)THEN
                                WRITE(*,*)"DISTRIBUTED FORCE AND COLLOCATION POINT AT SAME POSITION! CHANGE MESH!"
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
                                ! CHANGE LATER FROM MATMUL TO DMURRV
                                CALL SHAPE_FUNCTIONS(2,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),4,VALUES,COEFFICIENTS_AREA,X,Y,Z,PHI,DPHI,D2PHI)
                                N_P=0.D0
                                DO K=1,4
                                    N_P=N_P+PHI(K)*VALUES_DIST_LOADS(II,:,K)
                                ENDDO
                                F_ENRICHED_LOCAL=F_ENRICHED_LOCAL+MATMUL(U_AST,N_P)
                            ELSE
                                !---------------------------------------------------------
                                !       SOURCE POINT NORMAL OUTWARD MATRIX
                                !---------------------------------------------------------
!
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
                                CALL NORMAL_OUTWARD(4,COORDS_DIST_LOADS(II,:,:),QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),ETA,VJC,DVJC,JC)
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
                                N_P=0.D0
                                DO K=1,4
                                    N_P=N_P+PHI(K)*VALUES_DIST_LOADS(II,:,K)
                                ENDDO
!                                
                                F_ENRICHED_LOCAL=F_ENRICHED_LOCAL+MATMUL(AUX_D,N_P)
                            ENDIF
                        ENDDO
                    ENDDO
                ENDIF
            ELSE
                F_ENRICHED_LOCAL=0.D0
            ENDIF
            F_ENRICHED(3*JJ-2)=F_ENRICHED(3*JJ-2)+F_ENRICHED_LOCAL(1)
            F_ENRICHED(3*JJ-1)=F_ENRICHED(3*JJ-1)+F_ENRICHED_LOCAL(2)
            F_ENRICHED(3*JJ)=F_ENRICHED(3*JJ)+F_ENRICHED_LOCAL(3)
        ENDDO
    ENDDO
!
!
END SUBROUTINE