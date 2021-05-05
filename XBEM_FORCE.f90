    SUBROUTINE XBEM_FORCE
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
    INTEGER::I,J,K,L,M,II,KK,NEN,ELEMS,NODES,ELEM_NODE,&
    POS_COLLOCPOINT
!
    REAL*8::DX,DY,DZ,DR(3),U_AST(3,3), D_AST(9,3),R,&
    Mu,Nu,C1,C2,AUX1,KR(3,3),PI,VALUES[ALLOCATABLE](:,:),&
    ETAS(3),MAT_ETAS(3,9),VJC(3),DVJC(3,2),JC,AUX_D(3,3),&
    QSI_PNT(2),TOLDIST
! 
    PI=DACOS(-1.D0)
    F_ENRICHED_LOCAL=0.0D0
    !F_ENRICHED=0.0D0
    TOLDIST=1.0E-5
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
    DO I=1,NUM_CON_LOADS
!       IDENTIFYING ELEMENT CONTAINING CONCENTRATED SUPPORT AND CALCULATING ADIMENSIONAL PARAMETERS OF CONCENTRATED SUPPORT POINT
        ELEM_NODE=0
        CALL FIND_ELEMENT_WITH_PNT(COORDS_CON_LOADS(I,:),ELEM_NODE,QSI_PNT)
!
        DO J=1,N_COLLOCPOINTS
!           FINDING ELEMENT THAT CONTAINS COLLOCATION POINT
            DO K=1,N_ELEM
		        NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			    DO KK=1,NEN
				    IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.J) THEN
					    ELEMS=K
				    ENDIF
			    ENDDO
            ENDDO
            IF (ND(ELEM_NODE).EQ.ND(ELEMS)) THEN
                Mu=EMP(ND(ELEM_NODE),1)/(2*(1+EMP(ND(ELEM_NODE),2)))
	            Nu=EMP(ND(ELEM_NODE),2)
                C1=((1.D0)/(16.D0*PI*Mu*(1.D0-Nu)))
                C2=((1.D0)/(8.D0*PI*(1.D0-Nu)))
        	    DX=COORDS_CON_LOADS(I,1)-COORD_COLLOCPOINTS(J,1)
		        DY=COORDS_CON_LOADS(I,2)-COORD_COLLOCPOINTS(J,2)
		        DZ=COORDS_CON_LOADS(I,3)-COORD_COLLOCPOINTS(J,3)
		        R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                IF(R.LT.TOLDIST)THEN
                    WRITE(*,*)"CONCENTRATED FORCE AND COLLOCATION POINT AT SAME POSITION! CHANGE MESH!"
                    READ(*,*)
                ENDIF
		        DR(1)=DX/R
		        DR(2)=DY/R
		        DR(3)=DZ/R
                IF (EQ_TYPE(J).EQ."S") THEN
                    DO K=1,3
                        DO L=1,3
                            AUX1=(3.0D0-4*Nu)*KR(K,L)
                            U_AST(K,L)=C1*(AUX1+DR(K)*DR(L))/R
                        ENDDO
                    ENDDO
                    ! CHANGE LATER FROM MATMUL TO DMURRV
                    !CALL DMRRRR(3,3,U_AST,3,1,VALUES_CON_LOADS(I,:),3,1,F_ENRICHED_LOCAL)
                    F_ENRICHED_LOCAL=MATMUL(U_AST,VALUES_CON_LOADS(I,:))
                ELSE
                    !---------------------------------------------------------
                    !       SOURCE POINT NORMAL OUTWARD MATRIX
                    !---------------------------------------------------------
!
                    NEN=ELEM_TYPE(ELEMS)*ORDER_ELEM(ELEMS)+(ELEM_TYPE(ELEMS)-3)*(ORDER_ELEM(ELEMS)-1)*POL_FAMILY(ELEMS)
	                ALLOCATE(VALUES(NEN,3))
	                NODES=0
!                   NODES = LOCAL NODE
                    DO II=1,NEN
		                IF(COLLOCPOINTS_CONNECTIVITY(ELEMS,II).EQ.J) THEN
			                NODES=II
		                ENDIF
                    ENDDO
!                   VALUES = COORDINATES OF ELEMENT NODES		
	                DO II=1,NEN
		                VALUES(II,1)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,II),1)
		                VALUES(II,2)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,II),2)
		                VALUES(II,3)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,II),3)
	                ENDDO
!
	                CALL NORMAL_OUTWARD(NEN,VALUES,QSI(ELEMS,NODES,1),QSI(ELEMS,NODES,2),ETAS,VJC,DVJC,JC)
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
                    DEALLOCATE(VALUES)
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
                AUX_D=MATMUL(MAT_ETAS,D_AST)
                F_ENRICHED_LOCAL=MATMUL(AUX_D,VALUES_CON_LOADS(I,:))
                ENDIF
            ELSE
                F_ENRICHED_LOCAL=0.D0
            ENDIF    
            F_ENRICHED(3*J-2)=F_ENRICHED(3*J-2)+F_ENRICHED_LOCAL(1)
            F_ENRICHED(3*J-1)=F_ENRICHED(3*J-1)+F_ENRICHED_LOCAL(2)
            F_ENRICHED(3*J)=F_ENRICHED(3*J)+F_ENRICHED_LOCAL(3)
        ENDDO
    ENDDO
!
    END SUBROUTINE XBEM_FORCE