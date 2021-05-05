	SUBROUTINE SINGULAR_ELEMENT_S(NODE,ELEM,DOMAIN,DG,DH)
!
	USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE ANALYSIS
!
	IMPLICIT NONE
!
	INTEGER::NODE,ELEM,DOMAIN,II,I,J,K,L,NEN,LOCAL_NODE,N_DIV
!
	REAL*8::PI,DG(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),DH(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),&
	X,Y,Z,DX,DY,DZ,R,DRDN,VJC(3),DVJC(3,2),JC,Mu,Nu,C1,C2,AUX1,AUX2,AUX3,AUX4,DR(3),ETA(3),ETAS(3),KR(3,3),G_LOCAL(3,3),H_LOCAL(3,3),PHI[ALLOCATABLE](:),DGAUX1(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),&
	DHAUX1(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),DGAUX2(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),&
	DHAUX2(3,3*(ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM))),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),MATPHI[ALLOCATABLE](:,:),VALUES[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:)
	REAL*8::QSIW_GAUSS[ALLOCATABLE](:,:),AUX_CTE1,AUX_CTE2
!	
	REAL*8::A,B,Va(3),RHO,COS,SEN,DXDQSI(2),DYDQSI(2),DZDQSI(2),VJCS(3),DVJCS(3,2),JCS,THETA,THETAI[ALLOCATABLE](:),THETAF[ALLOCATABLE](:),RHO_BOUND,RHO_BOUND_AUX1,RHO_BOUND_AUX2(4),QSI1,QSI2
!
	PI=DACOS(-1.D0)
	Mu=EMP(DOMAIN,1)/(2*(1+EMP(DOMAIN,2)))
	Nu=EMP(DOMAIN,2)
	C1=((1.D0)/(16.D0*PI*Mu*(1.D0-Nu)))
	C2=((-1.D0)/(8.D0*PI*(1.D0-Nu)))
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
			ENDIF
		ENDDO
	ENDIF
!
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
        ELSE IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,LOCAL_NODE)).EQ."CONTINUOUS")THEN
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
        ELSE IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,LOCAL_NODE)).EQ."EDGEDISCONTINUOUS")THEN   
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
        ELSE IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,LOCAL_NODE)).EQ."CONTINUOUS")THEN
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
        ELSE IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEM,LOCAL_NODE)).EQ."EDGEDISCONTINUOUS")THEN    
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
    DO II=1,N_DIV       
        DO I=1,N_GAUSSPOINTS        
            THETA=(THETAF(II)-THETAI(II))/2*QSIW_GAUSS(I,1)+(THETAF(II)+THETAI(II))/2
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
		        DR=0.d0
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
                AUX_CTE1=(1.D0/R)*C1*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*RHO*JC*((THETAF(II)-THETAI(II))/2)*(RHO_BOUND/2.D0)
		        AUX_CTE2=(1.D0/(R*R))*C2*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*RHO*JC*((THETAF(II)-THETAI(II))/2)*(RHO_BOUND/2.D0)
                DO K=1,3
			        DO L=1,3
				        AUX1=0.D0
				        AUX1=(3.D0-4.D0*Nu)*KR(K,L)
				        AUX1=AUX1+DR(K)*DR(L)
				        G_LOCAL(K,L)=AUX1*AUX_CTE1
!
				        AUX2=0.D0
				        AUX2=(1.D0-2.D0*Nu)*KR(K,L)+3.D0*DR(K)*DR(L)
				        AUX2=AUX2*DRDN
			            AUX2=AUX2-(1.D0-2.D0*Nu)*(DR(K)*ETA(L)-DR(L)*ETA(K))
				        H_LOCAL(K,L)=AUX2*AUX_CTE2
			        ENDDO
		        ENDDO
!
                CALL SHAPE_FUNCTIONS(2,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!	
		        DO K=1,NEN
			        MATPHI(1,(3*K-2))=PHI(K)
			        MATPHI(2,(3*K-1))=PHI(K)
			        MATPHI(3,(3*K))=PHI(K)
		        ENDDO
!
		        CALL DMRRRR (3,3,G_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DGAUX1,3)
		        CALL DMRRRR (3,3,H_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHAUX1,3)
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
		        CALL NORMAL_OUTWARD(NEN,VALUES,QSI(ELEM,LOCAL_NODE,1),QSI(ELEM,LOCAL_NODE,2),ETAS,VJCS,DVJCS,JCS)
!              
                AUX_CTE1=(1.D0/(RHO*A))*C1*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*RHO*JCS*((THETAF(II)-THETAI(II))/2)*(RHO_BOUND/2.D0)
                !AUX_CTE2=(1.D0/(RHO*A)*(RHO*A))*C2*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*RHO*JCS*((THETAF(II)-THETAI(II))/2)*(RHO_BOUND/2.D0)
                AUX_CTE2=(1.D0/((RHO*A)*(RHO*A)))*C2*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*RHO*JCS*((THETAF(II)-THETAI(II))/2)*(RHO_BOUND/2.D0)
                DO K=1,3
			        DO L=1,3
				        AUX1=0.D0
				        AUX1=(3.D0-4.D0*Nu)*KR(K,L)
				        AUX1=AUX1+(Va(K)/A)*(Va(L)/A)
				        G_LOCAL(K,L)=AUX1*AUX_CTE1
!
				        AUX2=0.D0
				        AUX2=(1.D0-2.D0*Nu)*KR(K,L)+3.D0*(Va(K)/A)*(Va(L)/A)
				        AUX2=AUX2*((Va(1)/A)*ETAS(1)+(Va(2)/A)*ETAS(2)+(Va(3)/A)*ETAS(3))
			            AUX2=AUX2-(1.D0-2.D0*Nu)*((Va(K)/A)*ETAS(L)-(Va(L)/A)*ETAS(K))
				        H_LOCAL(K,L)=AUX2*AUX_CTE2
			        ENDDO
		        ENDDO
!
                CALL SHAPE_FUNCTIONS(2,QSI(ELEM,LOCAL_NODE,1),QSI(ELEM,LOCAL_NODE,2),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!	
		        DO K=1,NEN
			        MATPHI(1,(3*K-2))=PHI(K)
			        MATPHI(2,(3*K-1))=PHI(K)
			        MATPHI(3,(3*K))=PHI(K)
		        ENDDO
!
		        CALL DMRRRR (3,3,G_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DGAUX1,3)
		        CALL DMRRRR (3,3,H_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHAUX1,3)
!
		        DG=DG-DGAUX1
		        DH=DH-DHAUX1        		        		        
	        ENDDO
!***************************************************************************************************************************
!           COMPUTING SEMIANALYTICALLY THE VPC
!***************************************************************************************************************************
		    G_LOCAL=0.D0
		    H_LOCAL=0.D0		   
!
            Va(1)=DXDQSI(1)*DCOS(THETA)+DXDQSI(2)*DSIN(THETA)
            Va(2)=DYDQSI(1)*DCOS(THETA)+DYDQSI(2)*DSIN(THETA)
            Va(3)=DZDQSI(1)*DCOS(THETA)+DZDQSI(2)*DSIN(THETA)
!               
            A=DSQRT(Va(1)**2+Va(2)**2+Va(3)**2)
            B=1/A            
!   
            AUX_CTE1=(1.D0/A)*C1*(RHO_BOUND)*QSIW_GAUSS(I,2)*JCS*((THETAF(II)-THETAI(II))/2)
            AUX_CTE2=(1.D0/(A*A))*C2*DLOG(RHO_BOUND/B)*QSIW_GAUSS(I,2)*JCS*((THETAF(II)-THETAI(II))/2)
            DO K=1,3
			    DO L=1,3
				    AUX1=0.D0
				    AUX1=(3.D0-4.D0*Nu)*KR(K,L)
				    AUX1=AUX1+(Va(K)/A)*(Va(L)/A)
				    G_LOCAL(K,L)=AUX1*AUX_CTE1
!
				    AUX2=0.D0
				    AUX2=(1.D0-2.D0*Nu)*KR(K,L)+3.D0*(Va(K)/A)*(Va(L)/A)
				    AUX2=AUX2*((Va(1)/A)*ETAS(1)+(Va(2)/A)*ETAS(2)+(Va(3)/A)*ETAS(3))
			        AUX2=AUX2-(1.D0-2.D0*Nu)*((Va(K)/A)*ETAS(L)-(Va(L)/A)*ETAS(K))
				    H_LOCAL(K,L)=AUX2*AUX_CTE2
			    ENDDO
		    ENDDO    
!
		    CALL DMRRRR (3,3,G_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DGAUX1,3)
		    CALL DMRRRR (3,3,H_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHAUX1,3)
!
		    DG=DG+DGAUX1
		    DH=DH+DHAUX1                
	    ENDDO        
	ENDDO        
!	
	END SUBROUTINE