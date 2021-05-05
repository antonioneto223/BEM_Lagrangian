	SUBROUTINE XBEM_K_SINGULAR_S(NODE,ELEM,DOMAIN,P_ENR,N)
!
	USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE ANALYSIS
    USE CRACKS_DUAL_BEM
    USE XBEM_CRACKFRONT_VARIABLES
    USE PROPAGATION	   
	USE FATIGUE
!
	IMPLICIT NONE
!
	INTEGER::NODE,ELEM,DOMAIN,II,I,J,K,L,NEN,LOCAL_NODE,N_DIV,N,N_CRACKTIP_PAIR,POS_CT_PAIR,QSI_REF,N_GAUSS,&
    POS_CT(N_K),POS_ELEM,ELEM_CRACKTIP,NEN_CT_ELEM,ELEM_CT_S,NEN_CT_S,NEN_TEST
!
	REAL*8::PI,X,Y,Z,DX,DY,DZ,R,DRDN,VJC(3),DVJC(3,2),JC,Mu,Nu,C1,C2,AUX1,AUX2,AUX3,AUX4,DR(3),ETA(3),ETAS(3),KR(3,3),PHI[ALLOCATABLE](:),&
	DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),VALUES[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:)
	REAL*8::QSIW_GAUSS[ALLOCATABLE](:,:),DP_ENR(3,3),P_ENR_AUX(3,3),PHI1D_CT(N_K),PHI1D_AL(N_K),VALUES_K_ALIGN(N_K,3),COEFFICIENTS_K_ALIGN(N_K),&
    QSI1D_CT,QSI1D_AL,QSI1D_AL_CP(N_K),X_TIP,Y_TIP,Z_TIP,DR_DQSI1_CT(3),DR_DQSI2_CT(3),COORDS_CP_PROJECTION(3),QSI_PNT(2),SIGNAL,&
    VALUES_ELEM_CRACKTIP[ALLOCATABLE](:,:),VALUES_CT_S[ALLOCATABLE](:,:),VALUES_TEST[ALLOCATABLE](:,:),X_CT_S,Y_CT_S,Z_CT_S,X_TEST,Y_TEST,Z_TEST
!
	REAL*8::A,B,Va(3),RHO,COS,SEN,DXDQSI(2),DYDQSI(2),DZDQSI(2),VJCS(3),DVJCS(3,2),JCS,THETA,THETAI[ALLOCATABLE](:),THETAF[ALLOCATABLE](:),RHO_BOUND,RHO_BOUND_AUX1,RHO_BOUND_AUX2(4),QSI1,QSI2
    REAL*8::R_TIP,THETA_TIP,PSI_TIP(3,3),Z_FUN(3,3),Z_AUX(3,3),COORDS_TIP(3),R_TIP_COLLOC,THETA_TIP_COLLOC,&
    PSI_TIP_COLLOCPNT(3,3,N_K),ROTATION_MATRIX(3,3),P_ENR(3*N_K,3),DIST_AUX(N_K)  
!
    !OPEN(6,file='DP_SING.m',status='unknown')
    !WRITE(6,*)"dp_sing=["
	PI=DACOS(-1.D0)
	Mu=EMP(DOMAIN,1)/(2*(1+EMP(DOMAIN,2)))
	Nu=EMP(DOMAIN,2)
	C1=((1.D0)/(16.D0*PI*Mu*(1.D0-Nu)))
	C2=((-1.D0)/(8.D0*PI*(1.D0-Nu)))
    ROTATION_MATRIX=0.D0
    ROTATION_MATRIX(1,1)=1.D0
    ROTATION_MATRIX(2,2)=1.D0
    ROTATION_MATRIX(3,3)=1.D0
	KR(1,1)=1.D0
	KR(1,2)=0.D0
    KR(1,3)=0.D0
	KR(2,1)=0.D0
	KR(2,2)=1.D0
	KR(2,3)=0.D0
	KR(3,1)=0.D0
	KR(3,2)=0.D0
	KR(3,3)=1.D0
    N_GAUSS=N
    P_ENR=0.D0
    DP_ENR=0.D0
	ALLOCATE(QSIW_GAUSS(N_GAUSS,2))
	CALL GAUSS_POINTS(N_GAUSS,QSIW_GAUSS)
!
!   DEFINING IF ELEMENT IS UPPER CRACK SURFACE OR LOWER CRACK SURFACE BY EQUATION TYPE
    IF(DUAL_BEM(COLLOCPOINTS_CONNECTIVITY(ELEM,1)).EQ."S")THEN
	    SIGNAL=1.D0
    ELSE
	    SIGNAL=-1.D0
    ENDIF
!   ------------------------------------------------------------------------------------------------------
!   FINDING SOURCE POINT ON ELEMENT
	NEN=ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM)
    LOCAL_NODE=0.D0
	DO K=1,NEN
        IF(COLLOCPOINTS_CONNECTIVITY(ELEM,K).EQ.NODE)THEN
		    LOCAL_NODE=K
		ENDIF
	ENDDO
!   IF SOURCE POINT DOES NOT BELONG TO ELEMENT, FINDING NODE AT SAME POSITION
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
!   ------------------------------------------------------------------------------------------------------    
!   ELEMENT NODES COORDINATES
	ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2))
!
    CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM,NEN,COEFFICIENTS)
!	
	DO I=1,NEN
		VALUES(I,1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,I),1)
		VALUES(I,2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,I),2)
		VALUES(I,3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,I),3)
    ENDDO
!   ------------------------------------------------------------------------------------------------------
!   FINDING CRACK TIP PAIR ASSOCIATED WITH ELEMENT
    N_CRACKTIP_PAIR=SIZE(CRACK_TIP_PAIR,1)
    DO I=1,SIZE(ENRICHED_ELEMS,1)
        IF(ELEM.EQ.ENRICHED_ELEMS(I,1))THEN
            POS_ELEM=I
        ENDIF
    ENDDO
!
    DO I=1,N_K
        DO K=1,N_CRACKTIPS
            IF(ENRICHED_ELEMS(POS_ELEM,I+1).EQ.CRACKTIP_COLLOCPOINTS(K))THEN
                POS_CT(I)=K
            ENDIF
        ENDDO
    ENDDO
!
    DO I=1,N_CRACKTIP_PAIR
        DO J=1,N_K
            IF(CRACK_TIP_PAIR(I,J+2).EQ.ENRICHED_ELEMS(POS_ELEM,J+1))THEN
                IF(DUAL_BEM(COLLOCPOINTS_CONNECTIVITY(ELEM,1)).EQ."S")THEN
                    ELEM_CRACKTIP=CRACK_TIP_PAIR(I,1)
                ELSE
                    ELEM_CT_S=CRACK_TIP_PAIR(I,1)
                    NEN_CT_S=ELEM_TYPE(ELEM_CT_S)*ORDER_ELEM(ELEM_CT_S)+(ELEM_TYPE(ELEM_CT_S)-3)*(ORDER_ELEM(ELEM_CT_S)-1)*POL_FAMILY(ELEM_CT_S)
                    ALLOCATE(VALUES_CT_S(NEN_CT_S,3))
                    DO II=1,NEN_CT_S
                        VALUES_CT_S(II,1)=COORD_NODES(NODES_CONNECTIVITY(ELEM_CT_S,II),1)
                        VALUES_CT_S(II,2)=COORD_NODES(NODES_CONNECTIVITY(ELEM_CT_S,II),2)
                        VALUES_CT_S(II,3)=COORD_NODES(NODES_CONNECTIVITY(ELEM_CT_S,II),3)
                    ENDDO
                    CALL SHAPE_FUNCTIONS(1,0.D0,0.D0,NEN_CT_S,VALUES_CT_S,COEFFICIENTS,X_CT_S,Y_CT_S,Z_CT_S,PHI,DPHI,D2PHI)
                    DO K=1,N_ELEM
                        IF(DUAL_BEM(COLLOCPOINTS_CONNECTIVITY(K,1)).EQ."H")THEN
                            NEN_TEST=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
                            ALLOCATE(VALUES_TEST(NEN_TEST,3))
                            DO II=1,NEN_TEST
                                VALUES_TEST(II,1)=COORD_NODES(NODES_CONNECTIVITY(K,II),1)
                                VALUES_TEST(II,2)=COORD_NODES(NODES_CONNECTIVITY(K,II),2)
                                VALUES_TEST(II,3)=COORD_NODES(NODES_CONNECTIVITY(K,II),3)
                            ENDDO
                            CALL SHAPE_FUNCTIONS(1,0.D0,0.D0,NEN_TEST,VALUES_TEST,COEFFICIENTS,X_TEST,Y_TEST,Z_TEST,PHI,DPHI,D2PHI)
                            AUX1=X_CT_S-X_TEST
			                AUX2=Y_CT_S-Y_TEST
	                        AUX3=Z_CT_S-Z_TEST
			                AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
                            IF(AUX4.LT.1E-12)THEN
                                ELEM_CRACKTIP=K
                            ENDIF
                            DEALLOCATE(VALUES_TEST)
                        ENDIF
                    ENDDO
                    DEALLOCATE(VALUES_CT_S)
                ENDIF
                POS_CT_PAIR=I
                QSI_REF=CRACK_TIP_PAIR(I,2)
            ENDIF
        ENDDO
    ENDDO   
!   ------------------------------------------------------------------------------------------------------
!   VALUES_K_ALIGN RECEIVES COORDINATES OF EACH CRACK TIP
    DO K=1,N_K
        VALUES_K_ALIGN(K,1)=COORD_CRACKTIP(POS_CT(K),1)
        VALUES_K_ALIGN(K,2)=COORD_CRACKTIP(POS_CT(K),2)
        VALUES_K_ALIGN(K,3)=COORD_CRACKTIP(POS_CT(K),3)
    ENDDO
!   ------------------------------------------------------------------------------------------------------
!   COEFFICIENTS_K_ALIGN
    NEN_CT_ELEM=ELEM_TYPE(ELEM_CRACKTIP)*ORDER_ELEM(ELEM_CRACKTIP)+(ELEM_TYPE(ELEM_CRACKTIP)-3)*(ORDER_ELEM(ELEM_CRACKTIP)-1)*POL_FAMILY(ELEM_CRACKTIP)
    ALLOCATE(VALUES_ELEM_CRACKTIP(NEN_CT_ELEM,3))
    DO I=1,NEN_CT_ELEM
        VALUES_ELEM_CRACKTIP(I,1)=COORD_NODES(NODES_CONNECTIVITY(ELEM_CRACKTIP,I),1)
		VALUES_ELEM_CRACKTIP(I,2)=COORD_NODES(NODES_CONNECTIVITY(ELEM_CRACKTIP,I),2)
		VALUES_ELEM_CRACKTIP(I,3)=COORD_NODES(NODES_CONNECTIVITY(ELEM_CRACKTIP,I),3)
    ENDDO
!
    DO K=1,N_K
        QSI_PNT=0.D0
        CALL NEWTON_RAPHSON_ADIMENSIONAL_COORDINATES(NEN,VALUES_ELEM_CRACKTIP,VALUES_K_ALIGN(K,:),QSI_PNT)
        IF(QSI_REF.EQ.1)THEN
            COEFFICIENTS_K_ALIGN(K)=QSI_PNT(2)
        ELSE
            COEFFICIENTS_K_ALIGN(K)=QSI_PNT(1)
        ENDIF
    ENDDO
    DEALLOCATE(VALUES_ELEM_CRACKTIP)
!   ------------------------------------------------------------------------------------------------------
!   QSI1D_AL_CP STORES QSI 1D OF COLLOCATION POINTS, 1 IS FOR FURTHER CP AND 2 IS FOR CLOSER CP, LINEAR ONLY
    QSI1D_AL_CP(1)=-.8D0
    QSI1D_AL_CP(2)=.8D0
    DO I=1,N_K
        IF(QSI_REF.EQ.1)THEN
            CALL SHAPE_FUNCTIONS(1,QSI1D_AL_CP(I),0.8D0,NEN,VALUES,COEFFICIENTS,X_TEST,Y_TEST,Z_TEST,PHI,DPHI,D2PHI)
        ELSE
            CALL SHAPE_FUNCTIONS(1,0.8D0,QSI1D_AL_CP(I),NEN,VALUES,COEFFICIENTS,X_TEST,Y_TEST,Z_TEST,PHI,DPHI,D2PHI)
        ENDIF
        AUX1=VALUES_K_ALIGN(I,1)-X_TEST
		AUX2=VALUES_K_ALIGN(I,2)-Y_TEST
	    AUX3=VALUES_K_ALIGN(I,3)-Z_TEST
		DIST_AUX(I)=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
    ENDDO
!   
    IF(DIST_AUX(1).LT.DIST_AUX(2))THEN
        QSI1D_AL_CP(1)=.8D0
        QSI1D_AL_CP(2)=-.8D0    
    ENDIF
!   ------------------------------------------------------------------------------------------------------
!   ELEMENT DIVISION FOR INTEGRATION
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
!           COMPUTING RHO_BOUND(THETA)     
            SELECTCASE(ELEM_TYPE(ELEM))
!            CASE(3)          
!                RHO_BOUND_AUX2(1)=(-QSI(ELEM,LOCAL_NODE,1))/DCOS(THETA)
!                RHO_BOUND_AUX2(2)=(-QSI(ELEM,LOCAL_NODE,2))/DSIN(THETA)
!                RHO_BOUND_AUX2(3)=(1-QSI(ELEM,LOCAL_NODE,1)-QSI(ELEM,LOCAL_NODE,2))/(DCOS(THETA)+DSIN(THETA))
!!                        
!                RHO_BOUND=DSQRT(8.D0)
!                DO K=1,3
!                    IF(RHO_BOUND_AUX2(K).GT.1.D-10)THEN
!                        RHO_BOUND_AUX1=RHO_BOUND_AUX2(K)
!                        IF(RHO_BOUND_AUX1.LE.RHO_BOUND)THEN
!                            RHO_BOUND=RHO_BOUND_AUX1
!                        ENDIF
!                    ENDIF
!                ENDDO
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
!           ------------------------------------------------------------------------------------------------------
!           NUMERICAL INTEGRATION
            DO J=1,N_GAUSSPOINTS
                RHO=(RHO_BOUND/2.D0)*QSIW_GAUSS(J,1)+(RHO_BOUND/2.D0)            
		        DR=0.D0
		        DRDN=0.D0
		        ETA=0.D0
		        PHI=0.D0
		        QSI1=RHO*DCOS(THETA)+QSI(ELEM,LOCAL_NODE,1)
		        QSI2=RHO*DSIN(THETA)+QSI(ELEM,LOCAL_NODE,2)
!
                CALL SHAPE_FUNCTIONS(1,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!
!               FUNDAMENTAL SOLUTION P*
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
                DXDQSI=0.D0
		        DYDQSI=0.D0
		        DZDQSI=0.D0		        		        		        		                                
!		        
		        CALL SHAPE_FUNCTIONS(3,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
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
                AUX1=(1.D0/(RHO*A)**2)*C2*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*RHO*JC*((THETAF(II)-THETAI(II))/2)*(RHO_BOUND/2.D0)
                DO K=1,3
			        DO L=1,3
				        AUX2=0.D0
				        AUX2=(1.D0-2.D0*Nu)*KR(K,L)+3.D0*(Va(K)/A)*(Va(L)/A)
				        AUX2=AUX2*((Va(1)/A)*ETA(1)+(Va(2)/A)*ETA(2)+(Va(3)/A)*ETA(3))
			            AUX2=AUX2-(1.D0-2.D0*Nu)*((Va(K)/A)*ETA(L)-(Va(L)/A)*ETA(K))
				        P_ENR_AUX(K,L)=AUX1*AUX2
			        ENDDO
                ENDDO
!               ------------------------------------------------------------------------------------------------------
!               ENRICHED FUNCTION Z_FUN 
!               CALCULATING CRACK TIP COORDINATES ALIGNED TO GAUSS POINT
                IF(QSI_REF.EQ.1)THEN
                    QSI1D_CT=QSI2
                    QSI1D_AL=QSI1
                ELSE
                    QSI1D_CT=QSI1
                    QSI1D_AL=QSI2
                ENDIF
!           
                CALL SHAPE_FUNCTIONS_1D(2,QSI1D_CT,N_K-1,VALUES_K_ALIGN,COEFFICIENTS_K_ALIGN,X_TIP,Y_TIP,Z_TIP,PHI1D_CT,DPHI)           
!               CRACK TIP COORDINATES ALIGNED TO FIELD POINT
                COORDS_TIP=0.D0
                DO K=1,N_K
                    COORDS_TIP(1)=COORDS_TIP(1)+PHI1D_CT(K)*COORD_CRACKTIP(POS_CT(K),1)
                    COORDS_TIP(2)=COORDS_TIP(2)+PHI1D_CT(K)*COORD_CRACKTIP(POS_CT(K),2)
                    COORDS_TIP(3)=COORDS_TIP(3)+PHI1D_CT(K)*COORD_CRACKTIP(POS_CT(K),3)
                ENDDO
!               ------------------------------------------------------------------------------------------------------
!               CRACK TIP ROTATION MATRIX
!                IF(QSI_REF.EQ.1)THEN
!	                CALL SHAPE_FUNCTIONS(3,CT_PAIR_QSI(POS_CT_PAIR),QSI1D_CT,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!	                CALL NORMAL_OUTWARD(NEN,VALUES,CT_PAIR_QSI(POS_CT_PAIR),QSI1D_CT,ETA,VJC,DVJC,JC)
!                ELSE
!	                CALL SHAPE_FUNCTIONS(3,QSI1D_CT,CT_PAIR_QSI(POS_CT_PAIR),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!	                CALL NORMAL_OUTWARD(NEN,VALUES,QSI1D_CT,CT_PAIR_QSI(POS_CT_PAIR),ETA,VJC,DVJC,JC)
!                ENDIF
!                DR_DQSI1_CT=0.D0
!                DR_DQSI2_CT=0.D0
!                DO K=1,3
!                    DO L=1,NEN
!                        DR_DQSI1_CT(K)=DR_DQSI1_CT(K)+DPHI(L,1)*VALUES(L,K)
!                        DR_DQSI2_CT(K)=DR_DQSI2_CT(K)+DPHI(L,2)*VALUES(L,K)
!                    ENDDO
!                ENDDO
!                DR_DQSI1_CT=DR_DQSI1_CT/(DSQRT(DR_DQSI1_CT(1)**2+DR_DQSI1_CT(2)**2+DR_DQSI1_CT(3)**2))
!                DR_DQSI2_CT=DR_DQSI2_CT/(DSQRT(DR_DQSI2_CT(1)**2+DR_DQSI2_CT(2)**2+DR_DQSI2_CT(3)**2))
!!               CHANGE DIRECTION OF NORMAL OUTWARD BECAUSE DISPLACEMENT EQUATION
!                IF(DUAL_BEM(COLLOCPOINTS_CONNECTIVITY(ELEM,1)).EQ."S")THEN
!                    ETA=-ETA
!                    !DR_DQSI1_CT=-DR_DQSI1_CT
!                    !DR_DQSI2_CT=-DR_DQSI2_CT
!                ENDIF
!                IF(QSI_REF.EQ.1)THEN
!                    ROTATION_MATRIX(:,1)=DR_DQSI1_CT
!                    ROTATION_MATRIX(:,2)=ETA
!                    ROTATION_MATRIX(:,3)=DR_DQSI2_CT
!                ELSE
!                    ROTATION_MATRIX(:,1)=DR_DQSI2_CT
!                    ROTATION_MATRIX(:,2)=ETA
!                    ROTATION_MATRIX(:,3)=DR_DQSI1_CT
!                ENDIF                     
!               ------------------------------------------------------------------------------------------------------
!               WILLIAMS DISPL FUNCTIONS OF EACH COLLOCATION POINT ALIGNED TO FIELD POINT
                DO K=1,N_K
                    IF(QSI_REF.EQ.1)THEN
                        CALL SHAPE_FUNCTIONS(1,QSI1D_AL_CP(K),QSI1D_CT,NEN,VALUES,COEFFICIENTS,COORDS_CP_PROJECTION(1),COORDS_CP_PROJECTION(2),COORDS_CP_PROJECTION(3),PHI,DPHI,D2PHI)
                    ELSE
                        CALL SHAPE_FUNCTIONS(1,QSI1D_CT,QSI1D_AL_CP(K),NEN,VALUES,COEFFICIENTS,COORDS_CP_PROJECTION(1),COORDS_CP_PROJECTION(2),COORDS_CP_PROJECTION(3),PHI,DPHI,D2PHI)
                    ENDIF
                    DX=COORDS_CP_PROJECTION(1)-COORDS_TIP(1)
                    DY=COORDS_CP_PROJECTION(2)-COORDS_TIP(2)
                    DZ=COORDS_CP_PROJECTION(3)-COORDS_TIP(3)
                    R_TIP_COLLOC=DSQRT(DX*DX+DY*DY+DZ*DZ)
                    THETA_TIP_COLLOC=SIGNAL*ACOS((ROTATION_MATRIX(1,1)*DX+ROTATION_MATRIX(2,1)*DY+ROTATION_MATRIX(3,1)*DZ)/R_TIP_COLLOC)     
                    CALL WILLIAMS_DISPL_FUNCTIONS(R_TIP_COLLOC,THETA_TIP_COLLOC,Mu,Nu,PSI_TIP_COLLOCPNT(:,:,K))
                ENDDO
!               ------------------------------------------------------------------------------------------------------
!               WILLIAMS FUNCTIONS OF FIELD POINT AND CRACK TIP
		        DX=X-COORDS_TIP(1)
		        DY=Y-COORDS_TIP(2)
		        DZ=Z-COORDS_TIP(3)
                R_TIP=DSQRT(DX*DX+DY*DY+DZ*DZ)
                THETA_TIP=SIGNAL*ACOS((ROTATION_MATRIX(1,1)*DX+ROTATION_MATRIX(2,1)*DY+ROTATION_MATRIX(3,1)*DZ)/R_TIP)
                CALL WILLIAMS_DISPL_FUNCTIONS(R_TIP,THETA_TIP,Mu,Nu,PSI_TIP)
!               ------------------------------------------------------------------------------------------------------
!               ENRICHED Z FUNCTION
                Z_AUX=0.D0
!               VALUES_K_ALIGN IS USELESS IN FUNCTION BELOW
                CALL SHAPE_FUNCTIONS_1D(2,QSI1D_AL,N_K-1,VALUES_K_ALIGN,QSI1D_AL_CP,X_TIP,Y_TIP,Z_TIP,PHI1D_AL,DPHI)
                DO K=1,N_K
                    Z_AUX=Z_AUX+PHI1D_AL(K)*(PSI_TIP-PSI_TIP_COLLOCPNT(:,:,K))
                ENDDO   
!
                CALL DMRRRR (3,3,ROTATION_MATRIX,3,3,3,Z_AUX,3,3,3,Z_FUN,3)
!               ------------------------------------------------------------------------------------------------------
!               P*Z AND INTERPOLATION TO CRACK TIPS CONSIDERED
		        CALL DMRRRR (3,3,P_ENR_AUX,3,3,3,Z_FUN,3,3,3,DP_ENR,3)
                !WRITE(6,*)DP_ENR(1,1),QSI1D_AL,QSI1D_CT
!               INTERPOLATION TO CRACK TIPS CONSIDERED
                DO K=1,N_K
                    DO L=1,3
                        P_ENR(3*K-2,L)=P_ENR(3*K-2,L)+PHI1D_CT(K)*DP_ENR(1,L)
                        P_ENR(3*K-1,L)=P_ENR(3*K-1,L)+PHI1D_CT(K)*DP_ENR(2,L)
                        P_ENR(3*K,L)=P_ENR(3*K,L)+PHI1D_CT(K)*DP_ENR(3,L)
                    ENDDO
                ENDDO
!               ------------------------------------------------------------------------------------------------------
	        ENDDO
	    ENDDO        
    ENDDO   
!
    DEALLOCATE(VALUES)!,QSIW_GAUSS,COEFFICIENTS,PHI,DPHI,D2PHI)
!	
    !WRITE(6,*)"];"
    !CLOSE(6)
	END SUBROUTINE XBEM_K_SINGULAR_S