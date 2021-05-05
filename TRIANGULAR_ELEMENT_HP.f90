	SUBROUTINE TRIANGULAR_ELEMENT_HP(NODE,ELEM,ELEMS,DOMAIN,DG,DH,N)
!
	USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE ANALYSIS
!
	IMPLICIT NONE
!
	INTEGER::NODE,LOCAL_NODE,ELEM,ELEMS,DOMAIN,I,J,K,L,M,NEN,N,N_GAUSS
!
	REAL*8::PI,DG(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),&
	DH(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),X,Y,Z,DX,DY,DZ,R,DRDN,VJC(3),DVJC(3,2),JC,Mu,Nu,C1,C2,AUX1,AUX2,AUX3,AUX4,DR(3),&
	ETA(3),ETAS(3),KR(3,3),G_LOCAL(9,3),H_LOCAL(9,3),MAT_ETAS(3,9),PHI[ALLOCATABLE](:),DGAUX1(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),&
	DHAUX1(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),DGAUX2(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),&
	DHAUX2(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),&
    DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),MATPHI[ALLOCATABLE](:,:),VALUES[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:)
    REAL*8::QSIW_GAUSS[ALLOCATABLE](:,:)
	
    REAL*8::THETA,RHO,THETAI,THETAF,RHO_BOUND,RHO_BOUND_AUX1,RHO_BOUND_AUX2(3),QSI1,QSI2,QSI0(2)

	PI=DACOS(-1.D0)
	Mu=EMP(DOMAIN,1)/(2*(1+EMP(DOMAIN,2)))
	Nu=EMP(DOMAIN,2)
	C1=((1.D0)/(8.D0*PI*(1.D0-Nu)))
	C2=((Mu)/(4.D0*PI*(1.D0-Nu)))
	KR(1,1)=1.D0
	KR(1,2)=0.D0
    KR(1,3)=0.D0
	KR(2,1)=0.D0
	KR(2,2)=1.D0
	KR(2,3)=0.D0
	KR(3,1)=0.D0
	KR(3,2)=0.D0
	KR(3,3)=1.D0
	THETAI=3*PI/4
    THETAF=PI
    QSI0(1)=1.D0
    QSI0(2)=0.D0
    N_GAUSS=30*N
    ALLOCATE(QSIW_GAUSS(N_GAUSS,2))
    CALL GAUSS_POINTS(N_GAUSS,QSIW_GAUSS)
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   SOURCE POINT NORMAL OUTWARD MATRIX
!------------------------------------------------------------------------------------------------------------------------------------------------------  
    NEN=ELEM_TYPE(ELEMS)*ORDER_ELEM(ELEMS)+(ELEM_TYPE(ELEMS)-3)*(ORDER_ELEM(ELEMS)-1)*POL_FAMILY(ELEMS)
	ALLOCATE(VALUES(NEN,3))
	LOCAL_NODE=0
	DO I=1,NEN
		IF(COLLOCPOINTS_CONNECTIVITY(ELEMS,I).EQ.NODE) THEN
			LOCAL_NODE=I
		ENDIF
	ENDDO	
!	
	DO I=1,NEN
		VALUES(I,1)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,I),1)
		VALUES(I,2)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,I),2)
		VALUES(I,3)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,I),3)
	ENDDO
!
	CALL NORMAL_OUTWARD(NEN,VALUES,QSI(ELEMS,LOCAL_NODE,1),QSI(ELEMS,LOCAL_NODE,2),ETAS,VJC,DVJC,JC)
!
	MAT_ETAS(1,1)=ETAS(1)
	MAT_ETAS(1,2)=0.D0
	MAT_ETAS(1,3)=0.D0
	MAT_ETAS(1,4)=ETAS(2)
	MAT_ETAS(1,5)=0.D0
    MAT_ETAS(1,6)=0.D0
    MAT_ETAS(1,7)=ETAS(3)
	MAT_ETAS(1,8)=0.D0
    MAT_ETAS(1,9)=0.D0
	MAT_ETAS(2,1)=0.D0
	MAT_ETAS(2,2)=ETAS(1)
	MAT_ETAS(2,3)=0.D0
	MAT_ETAS(2,4)=0.D0
	MAT_ETAS(2,5)=ETAS(2)
    MAT_ETAS(2,6)=0.D0
    MAT_ETAS(2,7)=0.D0
	MAT_ETAS(2,8)=ETAS(3)
    MAT_ETAS(2,9)=0.D0
    MAT_ETAS(3,1)=0.D0
	MAT_ETAS(3,2)=0.D0
	MAT_ETAS(3,3)=ETAS(1)
	MAT_ETAS(3,4)=0.D0
	MAT_ETAS(3,5)=0.D0
    MAT_ETAS(3,6)=ETAS(2)
    MAT_ETAS(3,7)=0.D0
	MAT_ETAS(3,8)=0.D0
    MAT_ETAS(3,9)=ETAS(3)
!    
    DEALLOCATE(VALUES)	
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   SOURCE POINT NORMAL OUTWARD MATRIX
!------------------------------------------------------------------------------------------------------------------------------------------------------ 
    NEN=ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM)
	ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),MATPHI(3,3*NEN))
!	
	CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM,NEN,COEFFICIENTS)
!	
	DO I=1,NEN
		VALUES(I,1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,I),1)
		VALUES(I,2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,I),2)
		VALUES(I,3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,I),3)
	ENDDO
!
	DG=0.D0
	DH=0.D0
!
    DO I=1,N_GAUSS    
        THETA=(THETAF-THETAI)/2*QSIW_GAUSS(I,1)+(THETAF+THETAI)/2
!***************************************************************************************************************************  
!       COMPUTING RHO_BOUND(THETA)     
!***************************************************************************************************************************             
        RHO_BOUND_AUX2(1)=(-QSI0(1))/DCOS(THETA)
        RHO_BOUND_AUX2(2)=(-QSI0(2))/DSIN(THETA)
        RHO_BOUND_AUX2(3)=(1-QSI0(1)-QSI0(2))/(DCOS(THETA)+DSIN(THETA))
                        
        RHO_BOUND=DSQRT(8.D0)
        DO K=1,3
            IF(RHO_BOUND_AUX2(K).GT.0.D0)THEN
                RHO_BOUND_AUX1=RHO_BOUND_AUX2(K)
                IF(RHO_BOUND_AUX1.LE.RHO_BOUND)THEN
                    RHO_BOUND=RHO_BOUND_AUX1
                ENDIF
            ENDIF
        ENDDO    
!***************************************************************************************************************************                        
        DO J=1,N_GAUSS
            RHO=(RHO_BOUND/2.D0)*QSIW_GAUSS(J,1)+(RHO_BOUND/2.D0)            
		    DGAUX1=0.D0
		    DHAUX1=0.D0
		    DR=0.d0
		    DRDN=0.D0
		    ETA=0.D0
		    PHI=0.D0
		    MATPHI=0.D0
		    G_LOCAL=0.D0
		    H_LOCAL=0.D0
		    QSI1=RHO*DCOS(THETA)+QSI0(1)
		    QSI2=RHO*DSIN(THETA)+QSI0(2)
!
            CALL SHAPE_FUNCTIONS(1,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!        
		    DX=X-COORD_COLLOCPOINTS(NODE,1)
		    DY=Y-COORD_COLLOCPOINTS(NODE,2)
		    DZ=Z-COORD_COLLOCPOINTS(NODE,3)
		    R=DSQRT(DX*DX+DY*DY+DZ*DZ)
		    DR(1)=DX/R
		    DR(2)=DY/R
		    DR(3)=DZ/R
!
		    CALL NORMAL_OUTWARD(NEN,VALUES,QSI1,QSI2,ETA,VJC,DVJC,JC)
!
		    DRDN=DR(1)*ETA(1)+DR(2)*ETA(2)+DR(3)*ETA(3)
!
		    DO K=1,3
			    DO L=1,3
			        DO M=1,3
					    AUX1=0.D0
					    AUX1=KR(K,M)*DR(L)+KR(L,M)*DR(K)    
					    AUX1=AUX1-KR(K,L)*DR(M)
					    AUX1=AUX1*(1.D0-2.D0*Nu)
					    AUX1=AUX1+3.D0*DR(K)*DR(L)*DR(M)
					    G_LOCAL((L+3*K-3),M)=(1.D0/R**2)*C1*AUX1*RHO*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*JC*((THETAF-THETAI)/2)*(RHO_BOUND/2.D0)
!
					    AUX2=0.D0
					    AUX2=(1.D0-2.D0*Nu)*KR(K,L)*DR(M)
					    AUX3=KR(K,M)*DR(L)+KR(L,M)*DR(K)
					    AUX3=AUX3*Nu
					    AUX2=AUX2+AUX3-5.D0*DR(K)*DR(L)*DR(M)
					    AUX2=AUX2*3.D0*DRDN
!
					    AUX3=ETA(K)*DR(L)*DR(M)
					    AUX3=AUX3+ETA(L)*DR(K)*DR(M)
					    AUX3=AUX3*3.D0*Nu
					    AUX2=AUX2+AUX3
!
					    AUX3=3.D0*ETA(M)*DR(K)*DR(L)
					    AUX3=AUX3+ETA(L)*KR(K,M)
					    AUX3=AUX3+ETA(K)*KR(L,M)
					    AUX3=AUX3*(1.D0-2.D0*Nu)
					    AUX2=AUX2+AUX3
!
					    AUX3=-(1.D0-4.D0*Nu)*ETA(M)*KR(K,L)
					    AUX2=AUX2+AUX3
					    H_LOCAL((L+3*K-3),M)=(1.D0/R**3)*C2*AUX2*RHO*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*JC*((THETAF-THETAI)/2)*(RHO_BOUND/2.D0)
			        ENDDO
			    ENDDO
		    ENDDO
!		        
!		    INCLUDING THE SOMATION WITH ETA (II) SOURCE POINT --------------------------------------------------------------------------------
		    CALL DMRRRR (3,9,MAT_ETAS,3,9,3,G_LOCAL,9,3,3,DGAUX2,3)
		    CALL DMRRRR (3,9,MAT_ETAS,3,9,3,H_LOCAL,9,3,3,DHAUX2,3)
!		        
!           MULTIPLYING BY SHAPE FUNCTIONS ---------------------------------------------------------------------------------------------------		
            CALL SHAPE_FUNCTIONS(2,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!	
		    DO K=1,NEN
			    MATPHI(1,(3*K-2))=PHI(K)
			    MATPHI(2,(3*K-1))=PHI(K)
			    MATPHI(3,(3*K))=PHI(K)
		    ENDDO
!
		    CALL DMRRRR (3,3,DGAUX2,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DGAUX1,3)
		    CALL DMRRRR (3,3,DHAUX2,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHAUX1,3)
!
		    DG=DG+DGAUX1
		    DH=DH+DHAUX1
		ENDDO    
	ENDDO	        
!	
	END SUBROUTINE