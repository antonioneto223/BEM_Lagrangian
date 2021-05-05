	SUBROUTINE SINGULAR_ELEMENT_HP(NODE,ELEM,DOMAIN,DG,DH)
!
	USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE ANALYSIS
	USE CRACKS_DUAL_BEM
!
	IMPLICIT NONE
!
	INTEGER::NODE,ELEM,DOMAIN,III,I,J,K,II,JJ,KK,L,M,NEN,LOCAL_NODE,N_DIV,KODE_MECDUAL
!
	REAL*8::PI,DG(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),DH(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),&
	X,Y,Z,DX,DY,DZ,R,DRDN,VJC(3),DVJC(3,2),JC,Mu,Nu,C1,C2,DR(3),ETA(3),ETAS(3),KR(3,3),G_LOCAL(9,3),H_LOCAL(9,3),MAT_ETAS(3,9),PHI[ALLOCATABLE](:),DGAUX1(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),&
	DHAUX1(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),DGAUX2(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),&
	DHAUX2(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),MATPHI[ALLOCATABLE](:,:),VALUES[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:)
	REAL*8::QSIW_GAUSS[ALLOCATABLE](:,:),AUX_CTE1,AUX_CTE2
!	
	REAL*8::A,Beta,Gama,VaVb,Va(3),Vb(3),RHO,COS,SEN,DXDQSI(2),DYDQSI(2),DZDQSI(2),D2XD2QSI(2,2),D2YD2QSI(2,2),D2ZD2QSI(2,2),VJCS(3),DVJCS(3,2),JCS,THETA,THETAI[ALLOCATABLE](:),THETAF[ALLOCATABLE](:),RHO_BOUND,RHO_BOUND_AUX1,RHO_BOUND_AUX2(4),QSI1,QSI2,&
	JC0(3),JC1(3),AUX1,AUX2,AUX3,AUX4,AUX5,AUX6,AUX7,AUX8,AUX9,S_3,S_2,CPV,HFP
!
	PI=DACOS(-1.D0)
    KODE_MECDUAL=1
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
	ALLOCATE(QSIW_GAUSS(N_GAUSSPOINTS,2))
	CALL GAUSS_POINTS(N_GAUSSPOINTS,QSIW_GAUSS)
!	
	NEN=ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM)
	LOCAL_NODE=0.D0
	DO K=1,NEN
        IF(COLLOCPOINTS_CONNECTIVITY(ELEM,K).EQ.NODE)THEN
		    LOCAL_NODE=K
		ENDIF
	ENDDO
!
	IF(LOCAL_NODE.EQ.0) THEN
		DO K=1,NEN
			AUX1=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,K),1)-COORD_COLLOCPOINTS(NODE,1)
			AUX2=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,K),2)-COORD_COLLOCPOINTS(NODE,2)
	        AUX3=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,K),3)-COORD_COLLOCPOINTS(NODE,3)
			AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
			IF(AUX4.LE.1.E-12) THEN
				LOCAL_NODE=K
				IF(DUAL_BEM(NODE).EQ."H")THEN
				    KODE_MECDUAL=-1
				ENDIF
			ENDIF
		ENDDO
	ENDIF
!
	ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2),MATPHI(3,3*NEN))
!
    CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM,NEN,COEFFICIENTS)
!	
	DO I=1,NEN
		VALUES(I,1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,I),1)
		VALUES(I,2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,I),2)
		VALUES(I,3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,I),3)
	ENDDO
!
!------------------------------------------------------------------------------------------------------------------------------------------------------
!       SOURCE POINT NORMAL OUTWARD MATRIX
!------------------------------------------------------------------------------------------------------------------------------------------------------  
!
	CALL NORMAL_OUTWARD(NEN,VALUES,QSI(ELEM,LOCAL_NODE,1),QSI(ELEM,LOCAL_NODE,2),ETAS,VJCS,DVJCS,JCS)
!
	MAT_ETAS(1,1)=ETAS(1)*KODE_MECDUAL
	MAT_ETAS(1,2)=0.D0
	MAT_ETAS(1,3)=0.D0
	MAT_ETAS(1,4)=ETAS(2)*KODE_MECDUAL
	MAT_ETAS(1,5)=0.D0
    MAT_ETAS(1,6)=0.D0
    MAT_ETAS(1,7)=ETAS(3)*KODE_MECDUAL
	MAT_ETAS(1,8)=0.D0
    MAT_ETAS(1,9)=0.D0
	MAT_ETAS(2,1)=0.D0
	MAT_ETAS(2,2)=ETAS(1)*KODE_MECDUAL
	MAT_ETAS(2,3)=0.D0
	MAT_ETAS(2,4)=0.D0
	MAT_ETAS(2,5)=ETAS(2)*KODE_MECDUAL
    MAT_ETAS(2,6)=0.D0
    MAT_ETAS(2,7)=0.D0
	MAT_ETAS(2,8)=ETAS(3)*KODE_MECDUAL
    MAT_ETAS(2,9)=0.D0
    MAT_ETAS(3,1)=0.D0
	MAT_ETAS(3,2)=0.D0
	MAT_ETAS(3,3)=ETAS(1)*KODE_MECDUAL
	MAT_ETAS(3,4)=0.D0
	MAT_ETAS(3,5)=0.D0
    MAT_ETAS(3,6)=ETAS(2)*KODE_MECDUAL
    MAT_ETAS(3,7)=0.D0
	MAT_ETAS(3,8)=0.D0
    MAT_ETAS(3,9)=ETAS(3)*KODE_MECDUAL
!------------------------------------------------------------------------------------------------------------------------------------------------------
!       SUBDIVISION OF THETA INTEGRAL
!------------------------------------------------------------------------------------------------------------------------------------------------------ 
	DG=0.D0
	DH=0.D0
!    
    SELECTCASE(ELEM_TYPE(ELEM))	
    CASE(3)
        IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,LOCAL_NODE)).EQ."DISCONTINUOUS")THEN
            N_DIV=3
            ALLOCATE(THETAI(N_DIV),THETAF(N_DIV))
            THETAI(1)=-PI+DATAN((0.D0-QSI(ELEM,LOCAL_NODE,2))/(0.D0-QSI(ELEM,LOCAL_NODE,1)))
            THETAF(1)=DATAN((0.D0-QSI(ELEM,LOCAL_NODE,2))/(1.D0-QSI(ELEM,LOCAL_NODE,1)))
            THETAI(2)=DATAN((0.D0-QSI(ELEM,LOCAL_NODE,2))/(1.D0-QSI(ELEM,LOCAL_NODE,1)))
            THETAF(2)=PI+DATAN((1.D0-QSI(ELEM,LOCAL_NODE,2))/(0.D0-QSI(ELEM,LOCAL_NODE,1)))
            THETAI(3)=PI+DATAN((1.D0-QSI(ELEM,LOCAL_NODE,2))/(0.D0-QSI(ELEM,LOCAL_NODE,1))) 
            THETAF(3)=PI+DATAN((0.D0-QSI(ELEM,LOCAL_NODE,2))/(0.D0-QSI(ELEM,LOCAL_NODE,1)))
        ELSE IF (DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,LOCAL_NODE)).EQ."CONTINUOUS")THEN
            N_DIV=1
            ALLOCATE(THETAI(N_DIV),THETAF(N_DIV))
            SELECTCASE(LOCAL_NODE)
            CASE(1)
                THETAI(1)=3*PI/4
                THETAF(1)=PI
            CASE(2)
                THETAI(1)=3*PI/2
                THETAF(1)=7*PI/4            
            CASE(3)
                THETAI(1)=0.D0
                THETAF(1)=PI/2             
            CASE(4)
                THETAI(1)=3*PI/4
                THETAF(1)=7*PI/4           
            CASE(5)
                THETAI(1)=-PI/2
                THETAF(1)=PI/2            
            CASE(6)
                THETAI(1)=0.D0
                THETAF(1)=PI 
            ENDSELECT
        ELSE IF (DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,LOCAL_NODE)).EQ."EDGEDISCONTINUOUS")THEN         
            N_DIV=1
            ALLOCATE(THETAI(N_DIV),THETAF(N_DIV))
            SELECT CASE(EDGE_DISCONTINUITY(ELEM,LOCAL_NODE))
            CASE(1)
                THETAI(1)=3*PI/4
                THETAF(1)=7*PI/4             
            CASE(2)
                THETAI(1)=-PI/2
                THETAF(1)=PI/2            
            CASE(3)
                THETAI(1)=0.D0
                THETAF(1)=PI            
            ENDSELECT                
        ENDIF
    CASE(4)
        IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,LOCAL_NODE)).EQ."DISCONTINUOUS")THEN
            N_DIV=4
            ALLOCATE(THETAI(N_DIV),THETAF(N_DIV))
            THETAI(1)=-PI+DATAN((-1.D0-QSI(ELEM,LOCAL_NODE,2))/(-1.D0-QSI(ELEM,LOCAL_NODE,1)))
            THETAF(1)=DATAN((-1.D0-QSI(ELEM,LOCAL_NODE,2))/(1.D0-QSI(ELEM,LOCAL_NODE,1)))
            THETAI(2)=DATAN((-1.D0-QSI(ELEM,LOCAL_NODE,2))/(1.D0-QSI(ELEM,LOCAL_NODE,1)))
            THETAF(2)=DATAN((1.D0-QSI(ELEM,LOCAL_NODE,2))/(1.D0-QSI(ELEM,LOCAL_NODE,1)))
            THETAI(3)=DATAN((1.D0-QSI(ELEM,LOCAL_NODE,2))/(1.D0-QSI(ELEM,LOCAL_NODE,1)))
            THETAF(3)=PI+DATAN((1.D0-QSI(ELEM,LOCAL_NODE,2))/(-1.D0-QSI(ELEM,LOCAL_NODE,1)))            
            THETAI(4)=PI+DATAN((1.D0-QSI(ELEM,LOCAL_NODE,2))/(-1.D0-QSI(ELEM,LOCAL_NODE,1)))
            THETAF(4)=PI+DATAN((-1.D0-QSI(ELEM,LOCAL_NODE,2))/(-1.D0-QSI(ELEM,LOCAL_NODE,1)))         
        ELSE IF (DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,LOCAL_NODE)).EQ."CONTINUOUS")THEN 
            N_DIV=1
            ALLOCATE(THETAI(N_DIV),THETAF(N_DIV))
            SELECTCASE(LOCAL_NODE)
            CASE(1)
                THETAI(1)=0
                THETAF(1)=PI/2
            CASE(2)
                THETAI(1)=PI/2
                THETAF(1)=PI            
            CASE(3)
                THETAI=PI
                THETAF=3*PI/2             
            CASE(4)
                THETAI(1)=3*PI/2
                THETAF(1)=2*PI           
            CASE(5)
                THETAI(1)=0.D0
                THETAF(1)=PI            
            CASE(6)
                THETAI(1)=PI/2
                THETAF(1)=3*PI/2
            CASE(7)
                THETAI(1)=PI
                THETAF(1)=2*PI
            CASE(8)
                THETAI(1)=-PI/2
                THETAF(1)=PI/2
            CASE(9)
                DEALLOCATE(THETAI,THETAF)
                N_DIV=4
                ALLOCATE(THETAI(N_DIV),THETAF(N_DIV))
                THETAI(1)=-3*PI/4
                THETAF(1)=-PI/4
                THETAI(2)=-PI/4
                THETAF(2)=PI/4  
                THETAI(3)=PI/4
                THETAF(3)=3*PI/4
                THETAI(4)=3*PI/4
                THETAF(4)=5*PI/4                        
            ENDSELECT
        ELSE IF (DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,LOCAL_NODE)).EQ."EDGEDISCONTINUOUS")THEN     
            N_DIV=1
            ALLOCATE(THETAI(N_DIV),THETAF(N_DIV))
            SELECT CASE(EDGE_DISCONTINUITY(ELEM,LOCAL_NODE))
            CASE(1)
                THETAI(1)=0.D0
                THETAF(1)=PI              
            CASE(2)
                THETAI(1)=PI/2
                THETAF(1)=3*PI/2          
            CASE(3)
                THETAI(1)=PI
                THETAF(1)=2*PI 
            CASE(4)
                THETAI(1)=-PI/2
                THETAF(1)=PI/2                          
            ENDSELECT            
        ENDIF    
    ENDSELECT    
!       
    DO III=1,N_DIV       
        DO I=1,N_GAUSSPOINTS        
            THETA=(THETAF(III)-THETAI(III))/2*QSIW_GAUSS(I,1)+(THETAF(III)+THETAI(III))/2
!***************************************************************************************************************************  
!           COMPUTING RHO_BOUND(THETA)     
!***************************************************************************************************************************    
            SELECTCASE(ELEM_TYPE(ELEM))
            CASE(3)          
                RHO_BOUND_AUX2(1)=(-QSI(ELEM,LOCAL_NODE,1))/DCOS(THETA)
                RHO_BOUND_AUX2(2)=(-QSI(ELEM,LOCAL_NODE,2))/DSIN(THETA)
                RHO_BOUND_AUX2(3)=(1-QSI(ELEM,LOCAL_NODE,1)-QSI(ELEM,LOCAL_NODE,2))/(DCOS(THETA)+DSIN(THETA))
!                        
                RHO_BOUND=DSQRT(8.D0)
                DO K=1,3
                    IF(RHO_BOUND_AUX2(K).GT.1.D-10)THEN
                        RHO_BOUND_AUX1=RHO_BOUND_AUX2(K)
                        IF(RHO_BOUND_AUX1.LE.RHO_BOUND)THEN
                            RHO_BOUND=RHO_BOUND_AUX1
                        ENDIF
                    ENDIF
                ENDDO
            CASE(4)
                RHO_BOUND_AUX2(1)=(1-QSI(ELEM,LOCAL_NODE,1))/DCOS(THETA)
                RHO_BOUND_AUX2(2)=(-1-QSI(ELEM,LOCAL_NODE,1))/DCOS(THETA)
                RHO_BOUND_AUX2(3)=(1-QSI(ELEM,LOCAL_NODE,2))/DSIN(THETA)
                RHO_BOUND_AUX2(4)=(-1-QSI(ELEM,LOCAL_NODE,2))/DSIN(THETA)
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
            ENDSELECT    
!***************************************************************************************************************************                        
            DO J=1,N_GAUSSPOINTS
                RHO=(RHO_BOUND/2.D0)*QSIW_GAUSS(J,1)+(RHO_BOUND/2.D0)            
		        DGAUX1=0.D0
		        DHAUX1=0.D0
		        DR=0.D0
		        DRDN=0.D0
		        ETA=0.D0
		        PHI=0.D0
		        MATPHI=0.D0
		        G_LOCAL=0.D0
		        H_LOCAL=0.D0
		        QSI1=RHO*DCOS(THETA)+QSI(ELEM,LOCAL_NODE,1)
		        QSI2=RHO*DSIN(THETA)+QSI(ELEM,LOCAL_NODE,2)
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
                AUX_CTE1=(1.D0/(R*R))*C1*RHO*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*JC*((THETAF(III)-THETAI(III))/2)*(RHO_BOUND/2.D0)
                AUX_CTE2=(1.D0/(R*R*R))*C2*RHO*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*JC*((THETAF(III)-THETAI(III))/2)*(RHO_BOUND/2.D0)
		        DO II=1,3
			        DO JJ=1,3
			            DO KK=1,3
					        AUX1=0.D0
					        AUX1=KR(II,KK)*DR(JJ)+KR(JJ,KK)*DR(II)    
					        AUX1=AUX1-KR(II,JJ)*DR(KK)
					        AUX1=AUX1*(1.D0-2.D0*Nu)
					        AUX1=AUX1+3.D0*DR(II)*DR(JJ)*DR(KK)
					        G_LOCAL((JJ+3*II-3),KK)=AUX1*AUX_CTE1
!
					        AUX2=0.D0
					        AUX2=(1.D0-2.D0*Nu)*KR(II,JJ)*DR(KK)
					        AUX3=KR(II,KK)*DR(JJ)+KR(JJ,KK)*DR(II)
					        AUX3=AUX3*Nu
					        AUX2=AUX2+AUX3-5.D0*DR(II)*DR(JJ)*DR(KK)
					        AUX2=AUX2*3.D0*DRDN
!
					        AUX3=ETA(II)*DR(JJ)*DR(KK)
					        AUX3=AUX3+ETA(JJ)*DR(II)*DR(KK)
					        AUX3=AUX3*3.D0*Nu
					        AUX2=AUX2+AUX3
!
					        AUX3=3.D0*ETA(KK)*DR(II)*DR(JJ)
					        AUX3=AUX3+ETA(JJ)*KR(II,KK)
					        AUX3=AUX3+ETA(II)*KR(JJ,KK)
					        AUX3=AUX3*(1.D0-2.D0*Nu)
					        AUX2=AUX2+AUX3
!
					        AUX3=-(1.D0-4.D0*Nu)*ETA(KK)*KR(II,JJ)
					        AUX2=AUX2+AUX3
					        H_LOCAL((JJ+3*II-3),KK)=AUX_CTE2*AUX2
			            ENDDO
			        ENDDO
		        ENDDO
!		        
!		        INCLUDING THE SOMATION WITH ETA (II) SOURCE POINT --------------------------------------------------------------------------------
		        CALL DMRRRR (3,9,MAT_ETAS,3,9,3,G_LOCAL,9,3,3,DGAUX2,3)
		        CALL DMRRRR (3,9,MAT_ETAS,3,9,3,H_LOCAL,9,3,3,DHAUX2,3)
!		        
!        		MULTIPLYING BY SHAPE FUNCTIONS ---------------------------------------------------------------------------------------------------		
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
!***************************************************************************************************************************
!               COMPUTING THE SST TERMS
!***************************************************************************************************************************	        
		        G_LOCAL=0.D0
		        H_LOCAL=0.D0
		        DXDQSI=0.D0
		        DYDQSI=0.D0
		        DZDQSI=0.D0	
		        D2XD2QSI=0.D0
		        D2YD2QSI=0.D0
		        D2ZD2QSI=0.D0			        	        		        		        		                                
!		        
		        CALL SHAPE_FUNCTIONS(3,QSI(ELEM,LOCAL_NODE,1),QSI(ELEM,LOCAL_NODE,2),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!		        
		        DO K=1,NEN
		            DXDQSI(1)=DXDQSI(1)+DPHI(K,1)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),1)
		            DXDQSI(2)=DXDQSI(2)+DPHI(K,2)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),1)
!		            
		            DYDQSI(1)=DYDQSI(1)+DPHI(K,1)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),2)
		            DYDQSI(2)=DYDQSI(2)+DPHI(K,2)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),2)
!		            
		            DZDQSI(1)=DZDQSI(1)+DPHI(K,1)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),3)
		            DZDQSI(2)=DZDQSI(2)+DPHI(K,2)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),3)
!
		        ENDDO
!
                Va(1)=DXDQSI(1)*DCOS(THETA)+DXDQSI(2)*DSIN(THETA)
                Va(2)=DYDQSI(1)*DCOS(THETA)+DYDQSI(2)*DSIN(THETA)
                Va(3)=DZDQSI(1)*DCOS(THETA)+DZDQSI(2)*DSIN(THETA)
!               
                A=DSQRT(Va(1)**2+Va(2)**2+Va(3)**2)
!		        
		        CALL SHAPE_FUNCTIONS(5,QSI(ELEM,LOCAL_NODE,1),QSI(ELEM,LOCAL_NODE,2),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!		        
		        DO K=1,NEN
		            D2XD2QSI(1,1)=D2XD2QSI(1,1)+D2PHI(K,1,1)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),1)
		            D2XD2QSI(1,2)=D2XD2QSI(1,2)+D2PHI(K,1,2)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),1)
		            D2XD2QSI(2,1)=D2XD2QSI(2,1)+D2PHI(K,2,1)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),1)
		            D2XD2QSI(2,2)=D2XD2QSI(2,2)+D2PHI(K,2,2)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),1)
!		            
		            D2YD2QSI(1,1)=D2YD2QSI(1,1)+D2PHI(K,1,1)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),2)
		            D2YD2QSI(1,2)=D2YD2QSI(1,2)+D2PHI(K,1,2)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),2)
		            D2YD2QSI(2,1)=D2YD2QSI(2,1)+D2PHI(K,2,1)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),2)
		            D2YD2QSI(2,2)=D2YD2QSI(2,2)+D2PHI(K,2,2)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),2)
!		            
		            D2ZD2QSI(1,1)=D2ZD2QSI(1,1)+D2PHI(K,1,1)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),3)
		            D2ZD2QSI(1,2)=D2ZD2QSI(1,2)+D2PHI(K,1,2)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),3)
		            D2ZD2QSI(2,1)=D2ZD2QSI(2,1)+D2PHI(K,2,1)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),3)
		            D2ZD2QSI(2,2)=D2ZD2QSI(2,2)+D2PHI(K,2,2)*COORD_NODES(NODES_CONNECTIVITY(ELEM,K),3)
!	                
		        ENDDO
!
                Vb(1)=D2XD2QSI(1,2)*DSIN(THETA)*DCOS(THETA)+0.5D0*D2XD2QSI(1,1)*DCOS(THETA)**2+0.5D0*D2XD2QSI(2,2)*DSIN(THETA)**2
                Vb(2)=D2YD2QSI(1,2)*DSIN(THETA)*DCOS(THETA)+0.5D0*D2YD2QSI(1,1)*DCOS(THETA)**2+0.5D0*D2YD2QSI(2,2)*DSIN(THETA)**2
                Vb(3)=D2ZD2QSI(1,2)*DSIN(THETA)*DCOS(THETA)+0.5D0*D2ZD2QSI(1,1)*DCOS(THETA)**2+0.5D0*D2ZD2QSI(2,2)*DSIN(THETA)**2
!               
                VaVb=Va(1)*Vb(1)+Va(2)*Vb(2)+Va(3)*Vb(3)                
!                
                S_3=1/(A**3)
                S_2=-3*VaVb/(A**5)
!
                JC0=VJCS
                JC1(:)=DVJCS(:,1)*DCOS(THETA)+DVJCS(:,2)*DSIN(THETA)            
!                                 
                AUX_CTE1=(1.D0/(RHO*A)*(RHO*A))*C1*RHO*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*JCS*((THETAF(III)-THETAI(III))/2)*(RHO_BOUND/2.D0)
                AUX_CTE2=C2*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*((THETAF(III)-THETAI(III))/2)*(RHO_BOUND/2.D0)
                DO II=1,3  
			        DO JJ=1,3 
			            DO KK=1,3 
			            	AUX1=0.D0
					        AUX1=KR(II,KK)*(Va(JJ)/A)+KR(JJ,KK)*(Va(II)/A)    
					        AUX1=AUX1-KR(II,JJ)*(Va(KK)/A)
					        AUX1=AUX1*(1.D0-2.D0*Nu)
					        AUX1=AUX1+3.D0*(Va(II)/A)*(Va(JJ)/A)*(Va(KK)/A)
					        G_LOCAL((JJ+3*II-3),KK)=AUX1*AUX_CTE1
!
                            AUX1=0.D0
                            AUX1=3.D0*(1-2*Nu)*KR(II,JJ)*Va(KK)/A**2
                            AUX1=AUX1+3*Nu*(KR(II,KK)*Va(JJ)/A**2+KR(JJ,KK)*Va(II)/A**2)
                            AUX1=AUX1-15.D0*Va(II)*Va(JJ)*Va(KK)/A**4
                            AUX1=AUX1*(Vb(1)*JC0(1)+Vb(2)*JC0(2)+Vb(3)*JC0(3)+Va(1)*JC1(1)+Va(2)*JC1(2)+Va(3)*JC1(3))
!                            
                            AUX2=0.D0
                            AUX2=3*Nu*(JC0(II)*Va(JJ)*Va(KK)/A**2+JC0(JJ)*Va(II)*Va(KK)/A**2)
!
                            AUX3=0.D0
                            AUX3=JC0(II)*((Va(KK)/A)*(Vb(JJ)/A-Va(JJ)*VaVb/A**3)+(Va(JJ)/A)*(Vb(KK)/A-Va(KK)*VaVb/A**3))+JC1(II)*Va(JJ)*Va(KK)/A**2
                            AUX3=AUX3+JC0(JJ)*((Va(KK)/A)*(Vb(II)/A-Va(II)*VaVb/A**3)+(Va(II)/A)*(Vb(KK)/A-Va(KK)*VaVb/A**3))+JC1(JJ)*Va(II)*Va(KK)/A**2
                            AUX3=3*Nu*AUX3
!
                            AUX4=0.D0
                            AUX4=3*JC0(KK)*Va(II)*Va(JJ)/A**2
                            AUX4=AUX4+JC0(JJ)*KR(II,KK)+JC0(II)*KR(JJ,KK)
                            AUX4=(1-2*Nu)*AUX4
!                            
                            AUX5=0.D0
                            AUX5=JC0(KK)*((Va(JJ)/A)*(Vb(II)/A-Va(II)*VaVb/A**3)+(Va(II)/A)*(Vb(JJ)/A-Va(JJ)*VaVb/A**3))+JC1(KK)*Va(II)*Va(JJ)/A**2
                            AUX5=3*AUX5
                            AUX5=AUX5+JC1(JJ)*KR(II,KK)+JC1(II)*KR(JJ,KK)
                            AUX5=(1-2*Nu)*AUX5
!                            
                            AUX6=0.D0
                            AUX7=0.D0
                            AUX6=-(1-4*Nu)*JC0(KK)*KR(II,JJ)
                            AUX7=-(1-4*Nu)*JC1(KK)*KR(II,JJ)
!                            
                            AUX8=AUX2+AUX4+AUX6
                            AUX9=AUX1+AUX3+AUX5+AUX7					       
					        H_LOCAL((JJ+3*II-3),KK)=(S_3*AUX8/(RHO*RHO)+(S_2*AUX8+S_3*AUX9)/RHO)*AUX_CTE2
			            ENDDO
			        ENDDO
		        ENDDO
!		        
!		        INCLUDING THE SOMATION WITH ETA (II) SOURCE POINT --------------------------------------------------------------------------------
		        CALL DMRRRR (3,9,MAT_ETAS,3,9,3,G_LOCAL,9,3,3,DGAUX2,3)
		        CALL DMRRRR (3,9,MAT_ETAS,3,9,3,H_LOCAL,9,3,3,DHAUX2,3)
!		        
!        		MULTIPLYING BY SHAPE FUNCTIONS ---------------------------------------------------------------------------------------------------
                CALL SHAPE_FUNCTIONS(2,QSI(ELEM,LOCAL_NODE,1),QSI(ELEM,LOCAL_NODE,2),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
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
		        DG=DG-DGAUX1
		        DH=DH-DHAUX1 
!
		        H_LOCAL=0.D0
                AUX_CTE1=(S_3/RHO)*C2*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*((THETAF(III)-THETAI(III))/2)*(RHO_BOUND/2.D0)
                DO II=1,3  
			        DO JJ=1,3 
			            DO KK=1,3                             
                            AUX2=0.D0
                            AUX2=3*Nu*(JC0(II)*Va(JJ)*Va(KK)/A**2+JC0(JJ)*Va(II)*Va(KK)/A**2)
!
                            AUX4=0.D0
                            AUX4=3*JC0(KK)*Va(II)*Va(JJ)/A**2
                            AUX4=AUX4+JC0(JJ)*KR(II,KK)+JC0(II)*KR(JJ,KK)
                            AUX4=(1-2*Nu)*AUX4
!                            
                            AUX6=0.D0
                            AUX6=-(1-4*Nu)*JC0(KK)*KR(II,JJ)
!                            
                            AUX8=AUX2+AUX4+AUX6					       
					        H_LOCAL((JJ+3*II-3),KK)=AUX8*AUX_CTE1
			            ENDDO
			        ENDDO
		        ENDDO
!		        
!		        INCLUDING THE SOMATION WITH ETA (II) SOURCE POINT --------------------------------------------------------------------------------
		        CALL DMRRRR (3,9,MAT_ETAS,3,9,3,H_LOCAL,9,3,3,DHAUX2,3)	
!		        
!        		MULTIPLYING BY SHAPE FUNCTIONS DERIVATIVE ----------------------------------------------------------------------------------------	        
                CALL SHAPE_FUNCTIONS(4,QSI(ELEM,LOCAL_NODE,1),QSI(ELEM,LOCAL_NODE,2),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!	
                MATPHI=0.D0
		        DO K=1,NEN
			        MATPHI(1,(3*K-2))=DPHI(K,1)*DCOS(THETA)+DPHI(K,2)*DSIN(THETA)
			        MATPHI(2,(3*K-1))=DPHI(K,1)*DCOS(THETA)+DPHI(K,2)*DSIN(THETA)
			        MATPHI(3,(3*K))=DPHI(K,1)*DCOS(THETA)+DPHI(K,2)*DSIN(THETA)
		        ENDDO		        
!		        
		        CALL DMRRRR (3,3,DHAUX2,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHAUX1,3)
!
		        DH=DH-DHAUX1 		        		          	             		        		        
	        ENDDO
!***************************************************************************************************************************
!           COMPUTING SEMIANALYTICALLY THE VPC AND HFP
!***************************************************************************************************************************
		    G_LOCAL=0.D0
		    H_LOCAL=0.D0		   
!
            Va(1)=DXDQSI(1)*DCOS(THETA)+DXDQSI(2)*DSIN(THETA)
            Va(2)=DYDQSI(1)*DCOS(THETA)+DYDQSI(2)*DSIN(THETA)
            Va(3)=DZDQSI(1)*DCOS(THETA)+DZDQSI(2)*DSIN(THETA)
!               
            A=DSQRT(Va(1)**2+Va(2)**2+Va(3)**2)
            Beta=1.D0/A           
!
            Vb(1)=D2XD2QSI(1,2)*DSIN(THETA)*DCOS(THETA)+0.5D0*D2XD2QSI(1,1)*DCOS(THETA)**2+0.5D0*D2XD2QSI(2,2)*DSIN(THETA)**2
            Vb(2)=D2YD2QSI(1,2)*DSIN(THETA)*DCOS(THETA)+0.5D0*D2YD2QSI(1,1)*DCOS(THETA)**2+0.5D0*D2YD2QSI(2,2)*DSIN(THETA)**2
            Vb(3)=D2ZD2QSI(1,2)*DSIN(THETA)*DCOS(THETA)+0.5D0*D2ZD2QSI(1,1)*DCOS(THETA)**2+0.5D0*D2ZD2QSI(2,2)*DSIN(THETA)**2
!               
            VaVb=Va(1)*Vb(1)+Va(2)*Vb(2)+Va(3)*Vb(3)
            Gama=-VaVb/(A**4.D0)
!                
            S_3=1.D0/(A**3.D0)
            S_2=-3.D0*VaVb/(A**5.D0)
!
            JC0=VJCS
            JC1(:)=DVJCS(:,1)*DCOS(THETA)+DVJCS(:,2)*DSIN(THETA) 
! 
            CPV=DLOG(ABS(RHO_BOUND/Beta))    
!
            HFP=-(Gama/(Beta**2.D0)+1.D0/RHO_BOUND)
!
            G_LOCAL=0.D0
		    H_LOCAL=0.D0    
            AUX_CTE1=(1.D0/(A*A))*C1*CPV*QSIW_GAUSS(I,2)*JCS*((THETAF(III)-THETAI(III))/2)
            AUX_CTE2=C2*QSIW_GAUSS(I,2)*((THETAF(III)-THETAI(III))/2)
            DO II=1,3
			    DO JJ=1,3
			        DO KK=1,3
			        	AUX1=0.D0
					    AUX1=KR(II,KK)*Va(JJ)/A+KR(JJ,KK)*Va(II)/A    
					    AUX1=AUX1-KR(II,JJ)*Va(KK)/A
					    AUX1=AUX1*(1.D0-2.D0*Nu)
					    AUX1=AUX1+3.D0*(Va(II)/A)*(Va(JJ)/A)*(Va(KK)/A)
					    G_LOCAL((JJ+3*II-3),KK)=AUX1*AUX_CTE1
!
                        AUX1=0.D0
                        AUX1=3.D0*(1-2*Nu)*KR(II,JJ)*Va(KK)/A**2
                        AUX1=AUX1+3*Nu*(KR(II,KK)*Va(JJ)/A**2+KR(JJ,KK)*Va(II)/A**2)
                        AUX1=AUX1-15.D0*Va(II)*Va(JJ)*Va(KK)/A**4
                        AUX1=AUX1*(Vb(1)*JC0(1)+Vb(2)*JC0(2)+Vb(3)*JC0(3)+Va(1)*JC1(1)+Va(2)*JC1(2)+Va(3)*JC1(3))
!                            
                        AUX2=0.D0
                        AUX2=3*Nu*(JC0(II)*Va(JJ)*Va(KK)/A**2+JC0(JJ)*Va(II)*Va(KK)/A**2)
!
                        AUX3=0.D0
                        AUX3=JC0(II)*((Va(KK)/A)*(Vb(JJ)/A-Va(JJ)*VaVb/A**3)+(Va(JJ)/A)*(Vb(KK)/A-Va(KK)*VaVb/A**3))+JC1(II)*Va(JJ)*Va(KK)/A**2
                        AUX3=AUX3+JC0(JJ)*((Va(KK)/A)*(Vb(II)/A-Va(II)*VaVb/A**3)+(Va(II)/A)*(Vb(KK)/A-Va(KK)*VaVb/A**3))+JC1(JJ)*Va(II)*Va(KK)/A**2
                        AUX3=3*Nu*AUX3
!
                        AUX4=0.D0
                        AUX4=3*JC0(KK)*Va(II)*Va(JJ)/A**2
                        AUX4=AUX4+JC0(JJ)*KR(II,KK)+JC0(II)*KR(JJ,KK)
                        AUX4=(1-2*Nu)*AUX4
!                            
                        AUX5=0.D0
                        AUX5=JC0(KK)*((Va(JJ)/A)*(Vb(II)/A-Va(II)*VaVb/A**3)+(Va(II)/A)*(Vb(JJ)/A-Va(JJ)*VaVb/A**3))+JC1(KK)*Va(II)*Va(JJ)/A**2
                        AUX5=3*AUX5
                        AUX5=AUX5+JC1(JJ)*KR(II,KK)+JC1(II)*KR(JJ,KK)
                        AUX5=(1-2*Nu)*AUX5
!                            
                        AUX6=0.D0
                        AUX7=0.D0
                        AUX6=-(1-4*Nu)*JC0(KK)*KR(II,JJ)
                        AUX7=-(1-4*Nu)*JC1(KK)*KR(II,JJ)
!                            
                        AUX8=AUX2+AUX4+AUX6
                        AUX9=AUX1+AUX3+AUX5+AUX7						       
					    H_LOCAL((JJ+3*II-3),KK)=(S_3*AUX8*HFP+(S_2*AUX8+S_3*AUX9)*CPV)*AUX_CTE2
			        ENDDO
			    ENDDO
		    ENDDO    
!		        
!		    INCLUDING THE SOMATION WITH ETA (II) SOURCE POINT --------------------------------------------------------------------------------
		    CALL DMRRRR (3,9,MAT_ETAS,3,9,3,G_LOCAL,9,3,3,DGAUX2,3)
		    CALL DMRRRR (3,9,MAT_ETAS,3,9,3,H_LOCAL,9,3,3,DHAUX2,3)
!		        
!           MULTIPLYING BY SHAPE FUNCTIONS ---------------------------------------------------------------------------------------------------
            CALL SHAPE_FUNCTIONS(2,QSI(ELEM,LOCAL_NODE,1),QSI(ELEM,LOCAL_NODE,2),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!	
            MATPHI=0.D0
		    DO K=1,NEN
			    MATPHI(1,(3*K-2))=PHI(K)
			    MATPHI(2,(3*K-1))=PHI(K)
			    MATPHI(3,(3*K))=PHI(K)
		    ENDDO
		    CALL DMRRRR (3,3,DGAUX2,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DGAUX1,3)
		    CALL DMRRRR (3,3,DHAUX2,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHAUX1,3)
!
		    DG=DG+DGAUX1
		    DH=DH+DHAUX1
!
            H_LOCAL=0.D0
            AUX_CTE1=(S_3*CPV)*C2*QSIW_GAUSS(I,2)*((THETAF(III)-THETAI(III))/2)
            DO II=1,3
			    DO JJ=1,3
			        DO KK=1,3                           
                        AUX2=0.D0
                        AUX2=3*Nu*(JC0(II)*Va(JJ)*Va(KK)/A**2+JC0(JJ)*Va(II)*Va(KK)/A**2)
!
                        AUX4=0.D0
                        AUX4=3*JC0(KK)*Va(II)*Va(JJ)/A**2
                        AUX4=AUX4+JC0(JJ)*KR(II,KK)+JC0(II)*KR(JJ,KK)
                        AUX4=(1-2*Nu)*AUX4
!                            
                        AUX6=0.D0
                        AUX6=-(1-4*Nu)*JC0(KK)*KR(II,JJ)
!                            
                        AUX8=AUX2+AUX4+AUX6				       
					    H_LOCAL((JJ+3*II-3),KK)=AUX_CTE1*AUX8
			        ENDDO
			    ENDDO
		    ENDDO    
!		        
!		    INCLUDING THE SOMATION WITH ETA (II) SOURCE POINT --------------------------------------------------------------------------------
		    CALL DMRRRR (3,9,MAT_ETAS,3,9,3,H_LOCAL,9,3,3,DHAUX2,3)
!		        
!           MULTIPLYING BY SHAPE FUNCTIONS ---------------------------------------------------------------------------------------------------
           CALL SHAPE_FUNCTIONS(4,QSI(ELEM,LOCAL_NODE,1),QSI(ELEM,LOCAL_NODE,2),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!	
            MATPHI=0.D0
		    DO K=1,NEN
			    MATPHI(1,(3*K-2))=DPHI(K,1)*DCOS(THETA)+DPHI(K,2)*DSIN(THETA)
			    MATPHI(2,(3*K-1))=DPHI(K,1)*DCOS(THETA)+DPHI(K,2)*DSIN(THETA)
			    MATPHI(3,(3*K))=DPHI(K,1)*DCOS(THETA)+DPHI(K,2)*DSIN(THETA)
		    ENDDO	
!		        
		    CALL DMRRRR (3,3,DHAUX2,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHAUX1,3)
!
		    DH=DH+DHAUX1   		       		                   
	    ENDDO        
	ENDDO        
!	
	END SUBROUTINE