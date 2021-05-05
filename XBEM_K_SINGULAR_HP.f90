	SUBROUTINE XBEM_K_SINGULAR_HP(NODE,ELEM,DOMAIN,S_ENR,N)
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
    CHARACTER(40)::NOME	
!
	INTEGER::NODE,ELEM,ELEMS,DOMAIN,III,I,J,K,II,JJ,KK,L,M,N,NEN,LOCAL_NODE,N_DIV,KODE_MECDUAL,N_GAUSS,&
    N_CRACKTIP_PAIR,POS_CT_PAIR,QSI_REF,POS_CT(N_K),POS_ELEM,ELEM_CRACKTIP,NEN_CT_ELEM,NEN_TEST,ELEM_CT_S,NEN_CT_S
!
	REAL*8::PI,X,Y,Z,DX,DY,DZ,R,DRDN,VJC(3),DVJC(3,2),JC,Mu,Nu,C1,C2,DR(3),ETA(3),ETAS(3),KR(3,3),&
    PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),VALUES[ALLOCATABLE](:,:),&
    COEFFICIENTS[ALLOCATABLE](:,:),PHI_AL(N_K),DPHI_AL(N_K),X_CP,Y_CP,Z_CP,COORDS_CP_PROJECTION(3),&
    VALUES_ELEM_CRACKTIP[ALLOCATABLE](:,:),VALUES_CT_S[ALLOCATABLE](:,:),VALUES_TEST[ALLOCATABLE](:,:)
	REAL*8::QSIW_GAUSS[ALLOCATABLE](:,:)
    REAL*8::VALUES_K_ALIGN(N_K,3),COEFFICIENTS_K_ALIGN(N_K),QSI1D_CT,QSI1D_AL,X_TIP,Y_TIP,Z_TIP,DR_DQSI1_CT(3),DR_DQSI2_CT(3),&
    QSI1D_AL_CP(2),Z_AUX_DAL(3,3),Z_AUX_DAL1(3,3),SIGNAL,QSI_PNT(2),SIGNAL_QSI_AL_CT
    REAL*8::X_CT_S,Y_CT_S,Z_CT_S,X_TEST,Y_TEST,Z_TEST
!	
	REAL*8::A,Beta,Gama,VaVb,Va(3),Vb(3),RHO,COS,SEN,DXDQSI(2),DYDQSI(2),DZDQSI(2),D2XD2QSI(2,2),D2YD2QSI(2,2),D2ZD2QSI(2,2),VJCS(3),&
    DVJCS(3,2),JCS,THETA,THETAI[ALLOCATABLE](:),THETAF[ALLOCATABLE](:),RHO_BOUND,RHO_BOUND_AUX1,RHO_BOUND_AUX2(4),QSI1,QSI2,&
	JC0(3),JC1(3),AUX1,AUX2,AUX3,AUX4,AUX5,AUX6,AUX7,AUX8,AUX9,S_3,S_2,CPV,AUX_CTE
!
    REAL*8::R_TIP,THETA_TIP,PSI_TIP(3,3),Z_FUN(3,3),COORDS_TIP(3),R_TIP_COLLOC,THETA_TIP_COLLOC,&
    PSI_TIP_COLLOCPNT(3,3,N_K),ROTATION_MATRIX(3,3),S_ENR(3*N_K,3),S_ENR_LOCAL(3,3),S_ENR_0(3,3),S_ENR_1(3,3),S_ENR_AUX1(3,3,3),S_ENR_AUX2(3,3),&
    Z_FUN_DAL(3,3),PSI_TIP_DERIV(3,3),Z_AUX(3,3),Z_AUX_D1(3,3),Z_AUX_D2(3,3),Z_FUN_D(3,3),AUX_DBL,PHI1D_CT(N_K),PHI1DS_CT(N_K),PHI1D_AL(N_K),PSI_SP(3,3),DIST_AUX(N_K)
!   TESTES
    REAL*8::S_ENR_1_TESTE(3*N_K,3),S_ENR_0_TESTE(3*N_K,3),S_ENR_VPC_TESTE(3*N_K,3),Z_FUN_TEST(3,3)
!    
	PI=DACOS(-1.D0)
    KODE_MECDUAL=1
	Mu=EMP(DOMAIN,1)/(2*(1+EMP(DOMAIN,2)))
	Nu=EMP(DOMAIN,2)
	C1=((1.D0)/(8.D0*PI*(1.D0-Nu)))
	C2=((Mu)/(4.D0*PI*(1.D0-Nu)))
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
    N_GAUSS=10*N
    S_ENR=0.D0
    S_ENR_0_TESTE=0.D0
    S_ENR_1_TESTE=0.D0
    S_ENR_VPC_TESTE=0.D0
	ALLOCATE(QSIW_GAUSS(N_GAUSSPOINTS,2))
	CALL GAUSS_POINTS(N_GAUSSPOINTS,QSIW_GAUSS)
    !OPEN(6,file='PSI_SING_HP_0.m',status='unknown')
    !OPEN(7,file='PSI_SING_HP_1.m',status='unknown')
    !OPEN(8,file='Z_SING_HP_0.m',status='unknown')
    !OPEN(9,file='Z_SING_HP_1.m',status='unknown')
    !OPEN(10,file='SZ_SING_HP_0.m',status='unknown')
    !OPEN(11,file='SZ_SING_HP_1.m',status='unknown')
    !WRITE(6,*)"psi_sing_hp_0=["
    !WRITE(7,*)"psi_sing_hp_1=["
    !WRITE(8,*)"z_sing_hp_0=["
    !WRITE(9,*)"z_sing_hp_1=["
    !WRITE(10,*)"s_z_sing_hp_0=["
    !WRITE(11,*)"s_z_sing_hp_1=["
!
!   DEFINING IF ELEMENT IS UPPER CRACK SURFACE OR LOWER CRACK SURFACE BY EQUATION TYPE    
    SIGNAL=0.D0
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
				IF(DUAL_BEM(NODE).EQ."H")THEN
				    KODE_MECDUAL=-1
				ENDIF
			ENDIF
		ENDDO
    ENDIF
!   ------------------------------------------------------------------------------------------------------    
!   ELEMENT NODES COORDINATES
    ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2))
!
    CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM,NEN,COEFFICIENTS)
!	
	DO I=1,NEN
		VALUES(I,1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,I),1)
		VALUES(I,2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,I),2)
		VALUES(I,3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,I),3)
    ENDDO   
!   SOURCE POINT NORMAL OUTWARD MATRIX
    CALL NORMAL_OUTWARD(NEN,VALUES,QSI(ELEM,LOCAL_NODE,1),QSI(ELEM,LOCAL_NODE,2),ETAS,VJCS,DVJCS,JCS)
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
    QSI1D_AL_CP(1)=-.8D0
    QSI1D_AL_CP(2)=.8D0
    SIGNAL_QSI_AL_CT=-1.D0
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
        SIGNAL_QSI_AL_CT=1.D0
    ENDIF
!   ----------------------------------------------------------------------------------------------- 	    
!   WILLIAMS DISPLACEMENT DERIVATIVE AT ADIMENSIONAL DIRECTIONS EVALUATED AT SOURCE POINT
!   CALCULATING COORDINATES OF CRACK TIP ASSOCIATED WITH SOURCE POINT
    IF(QSI_REF.EQ.1)THEN
        QSI1D_CT=QSI(ELEM,LOCAL_NODE,2)
        QSI1D_AL=QSI(ELEM,LOCAL_NODE,1)
    ELSE
        QSI1D_CT=QSI(ELEM,LOCAL_NODE,1)
        QSI1D_AL=QSI(ELEM,LOCAL_NODE,2)
    ENDIF
!   PHI 1D OF CRACK TIP IN SOURCE POINT
    CALL SHAPE_FUNCTIONS_1D(2,QSI1D_CT,N_K-1,VALUES_K_ALIGN,COEFFICIENTS_K_ALIGN,X_TIP,Y_TIP,Z_TIP,PHI1DS_CT,DPHI)         
!
    COORDS_TIP=0.D0
    DO K=1,N_K
        COORDS_TIP(1)=COORDS_TIP(1)+PHI1DS_CT(K)*COORD_CRACKTIP(POS_CT(K),1)
        COORDS_TIP(2)=COORDS_TIP(2)+PHI1DS_CT(K)*COORD_CRACKTIP(POS_CT(K),2)
        COORDS_TIP(3)=COORDS_TIP(3)+PHI1DS_CT(K)*COORD_CRACKTIP(POS_CT(K),3)
    ENDDO
!   CRACK TIP ROTATION MATRIX
!    IF(QSI_REF.EQ.1)THEN
!	    CALL SHAPE_FUNCTIONS(3,CT_PAIR_QSI(POS_CT_PAIR),QSI1D_CT,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!	    CALL NORMAL_OUTWARD(NEN,VALUES,CT_PAIR_QSI(POS_CT_PAIR),QSI1D_CT,ETA,VJC,DVJC,JC)
!    ELSE
!	    CALL SHAPE_FUNCTIONS(3,QSI1D_CT,CT_PAIR_QSI(POS_CT_PAIR),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!	    CALL NORMAL_OUTWARD(NEN,VALUES,QSI1D_CT,CT_PAIR_QSI(POS_CT_PAIR),ETA,VJC,DVJC,JC)
!    ENDIF
!    DR_DQSI1_CT=0.D0
!    DR_DQSI2_CT=0.D0
!    DO K=1,3
!        DO L=1,NEN
!            DR_DQSI1_CT(K)=DR_DQSI1_CT(K)+DPHI(L,1)*VALUES(L,K)
!            DR_DQSI2_CT(K)=DR_DQSI2_CT(K)+DPHI(L,2)*VALUES(L,K)
!        ENDDO
!    ENDDO
!    DR_DQSI1_CT=DR_DQSI1_CT/(DSQRT(DR_DQSI1_CT(1)**2+DR_DQSI1_CT(2)**2+DR_DQSI1_CT(3)**2))
!    DR_DQSI2_CT=DR_DQSI2_CT/(DSQRT(DR_DQSI2_CT(1)**2+DR_DQSI2_CT(2)**2+DR_DQSI2_CT(3)**2))
!!   CHANGE DIRECTION OF NORMAL OUTWARD BECAUSE DISPLACEMENT EQUATION
!    IF(DUAL_BEM(COLLOCPOINTS_CONNECTIVITY(ELEM,1)).EQ."S")THEN
!        ETA=-ETA
!        !DR_DQSI1_CT=-DR_DQSI1_CT
!        !DR_DQSI2_CT=-DR_DQSI2_CT
!    ENDIF
!    IF(QSI_REF.EQ.1)THEN
!        ROTATION_MATRIX(:,1)=DR_DQSI1_CT
!        ROTATION_MATRIX(:,2)=ETA
!        ROTATION_MATRIX(:,3)=DR_DQSI2_CT
!    ELSE
!        ROTATION_MATRIX(:,1)=DR_DQSI2_CT
!        ROTATION_MATRIX(:,2)=ETA
!        ROTATION_MATRIX(:,3)=DR_DQSI1_CT
!    ENDIF          
!   ----------------------------------------------------------------------------------------------- 
!   WILLIAMS DISPLACEMENT DERIVATIVE
    CALL WILLIAMS_DISPL_FUNCTIONS_DERIV(QSI(ELEM,LOCAL_NODE,1),QSI(ELEM,LOCAL_NODE,2),QSI_REF,COORDS_TIP,ROTATION_MATRIX,VALUES,NEN,ELEM,Mu,Nu,SIGNAL,SIGNAL_QSI_AL_CT,PSI_TIP_DERIV)
    CALL SHAPE_FUNCTIONS_1D(2,QSI1D_AL,N_K-1,VALUES,QSI1D_AL_CP,X,Y,Z,PHI_AL,DPHI_AL)
    CALL SHAPE_FUNCTIONS_1D(3,QSI1D_AL,N_K-1,VALUES,QSI1D_AL_CP,X,Y,Z,PHI_AL,DPHI_AL)
    Z_AUX_DAL=0.D0
!   WILLIAMS DISPLACEMENT FUNCTIONS AT SOURCE POINT
    IF(QSI_REF.EQ.1)THEN
        CALL SHAPE_FUNCTIONS(1,QSI1D_AL,QSI1D_CT,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
    ELSE
        CALL SHAPE_FUNCTIONS(1,QSI1D_CT,QSI1D_AL,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
    ENDIF    
    DX=X-COORDS_TIP(1)
    DY=Y-COORDS_TIP(2)
    DZ=Z-COORDS_TIP(3)
    R_TIP=DSQRT(DX*DX+DY*DY+DZ*DZ)
    THETA_TIP=SIGNAL*ACOS((ROTATION_MATRIX(1,1)*DX+ROTATION_MATRIX(2,1)*DY+ROTATION_MATRIX(3,1)*DZ)/R_TIP)
    CALL WILLIAMS_DISPL_FUNCTIONS(R_TIP,THETA_TIP,Mu,Nu,PSI_SP)
!
    DO K=1,N_K
        IF(QSI_REF.EQ.1)THEN
            CALL SHAPE_FUNCTIONS(1,QSI1D_AL_CP(K),QSI1D_CT,NEN,VALUES,COEFFICIENTS,X_CP,Y_CP,Z_CP,PHI,DPHI,D2PHI)
        ELSE
            CALL SHAPE_FUNCTIONS(1,QSI1D_CT,QSI1D_AL_CP(K),NEN,VALUES,COEFFICIENTS,X_CP,Y_CP,Z_CP,PHI,DPHI,D2PHI)
        ENDIF
!       R_TIP AND THETA_TIP OF COLLOCATION POINT 
        DX=X_CP-COORDS_TIP(1)
        DY=Y_CP-COORDS_TIP(2)
        DZ=Z_CP-COORDS_TIP(3)
        R_TIP_COLLOC=DSQRT(DX*DX+DY*DY+DZ*DZ)
        THETA_TIP_COLLOC=SIGNAL*ACOS((ROTATION_MATRIX(1,1)*DX+ROTATION_MATRIX(2,1)*DY+ROTATION_MATRIX(3,1)*DZ)/R_TIP_COLLOC)
        CALL WILLIAMS_DISPL_FUNCTIONS(R_TIP_COLLOC,THETA_TIP_COLLOC,Mu,Nu,PSI_TIP_COLLOCPNT(:,:,K))
!
        Z_AUX_DAL=Z_AUX_DAL+DPHI_AL(K)*(PSI_SP-PSI_TIP_COLLOCPNT(:,:,K))+PHI_AL(K)*PSI_TIP_DERIV
    ENDDO  
!
    CALL DMRRRR(3,3,ROTATION_MATRIX,3,3,3,Z_AUX_DAL,3,3,3,Z_AUX_DAL1,3) 
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
!           COMPUTING RHO_BOUND(THETA)     
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
                CALL SHAPE_FUNCTIONS(1,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
                CALL SHAPE_FUNCTIONS(2,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
                CALL SHAPE_FUNCTIONS(4,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,AUX_DBL,AUX_DBL,AUX_DBL,PHI,DPHI,D2PHI)
!
!               FUNDAMENTAL SOLUTION S*
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
                S_ENR_AUX1=0.D0
                AUX_CTE=0.D0
                AUX_CTE=(1.D0/(R*R*R))*C2*RHO*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*JC*((THETAF(III)-THETAI(III))/2)*(RHO_BOUND/2.D0)
                DO KK=1,3
			        DO II=1,3
			            DO JJ=1,3
                            AUX2=0.D0
                            AUX3=0.D0
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
					        S_ENR_AUX1(KK,II,JJ)=AUX_CTE*AUX2
			            ENDDO
			        ENDDO
                ENDDO
!		        
!		        INCLUDING THE SOMATION WITH ETA (II) SOURCE POINT 
                S_ENR_AUX2=0.D0
                DO KK=1,3
                    DO JJ=1,3
                        DO II=1,3
                            S_ENR_AUX2(KK,JJ)=S_ENR_AUX2(KK,JJ)+S_ENR_AUX1(KK,II,JJ)*ETAS(II)
                        ENDDO
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
!                    CALL SHAPE_FUNCTIONS(3,CT_PAIR_QSI(POS_CT_PAIR),QSI1D_CT,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!                    CALL NORMAL_OUTWARD(NEN,VALUES,CT_PAIR_QSI(POS_CT_PAIR),QSI1D_CT,ETA,VJC,DVJC,JC)
!                ELSE
!                    CALL SHAPE_FUNCTIONS(3,QSI1D_CT,CT_PAIR_QSI(POS_CT_PAIR),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!                    CALL NORMAL_OUTWARD(NEN,VALUES,QSI1D_CT,CT_PAIR_QSI(POS_CT_PAIR),ETA,VJC,DVJC,JC)
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
                R_TIP_COLLOC=0.D0
                THETA_TIP_COLLOC=0.D0
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
                DX=0.D0
                DY=0.D0
                DZ=0.D0
                DX=X-COORDS_TIP(1)
		        DY=Y-COORDS_TIP(2)
		        DZ=Z-COORDS_TIP(3)
                R_TIP=DSQRT(DX*DX+DY*DY+DZ*DZ)
                THETA_TIP=SIGNAL*ACOS((ROTATION_MATRIX(1,1)*DX+ROTATION_MATRIX(2,1)*DY+ROTATION_MATRIX(3,1)*DZ)/R_TIP)
                CALL WILLIAMS_DISPL_FUNCTIONS(R_TIP,THETA_TIP,Mu,Nu,PSI_TIP)
                !WRITE(6,*)PSI_TIP(2,1),QSI1D_AL,QSI1D_CT,";"
!               ------------------------------------------------------------------------------------------------------
!               ENRICHED Z FUNCTION
                Z_AUX=0.D0
!               VALUES_K_ALIGN IS USELESS IN FUNCTION BELOW
                CALL SHAPE_FUNCTIONS_1D(2,QSI1D_AL,N_K-1,VALUES_K_ALIGN,QSI1D_AL_CP,X_TIP,Y_TIP,Z_TIP,PHI1D_AL,DPHI)
                DO K=1,N_K
                    Z_AUX=Z_AUX+PHI1D_AL(K)*(PSI_TIP-PSI_TIP_COLLOCPNT(:,:,K))
                ENDDO  
!
                CALL DMRRRR(3,3,ROTATION_MATRIX,3,3,3,Z_AUX,3,3,3,Z_FUN,3)   
                !WRITE(8,*)Z_FUN(2,1),QSI1D_AL,QSI1D_CT,";"
!               ------------------------------------------------------------------------------------------------------
!               S*Z AND INTERPOLATION TO CRACK TIPS CONSIDERED
                CALL DMRRRR (3,3,S_ENR_AUX2,3,3,3,Z_FUN,3,3,3,S_ENR_0,3)
                !WRITE(10,*)S_ENR_0(2,1),QSI1D_AL,QSI1D_CT,";"
!       
                DO K=1,N_K
                    DO L=1,3
                        S_ENR(3*K-2,L)=S_ENR(3*K-2,L)+PHI1D_CT(K)*S_ENR_0(1,L)
                        S_ENR(3*K-1,L)=S_ENR(3*K-1,L)+PHI1D_CT(K)*S_ENR_0(2,L)
                        S_ENR(3*K,L)=S_ENR(3*K,L)+PHI1D_CT(K)*S_ENR_0(3,L)
                        S_ENR_0_TESTE(3*K-2,L)=S_ENR_0_TESTE(3*K-2,L)+PHI1D_CT(K)*S_ENR_0(1,L)
                        S_ENR_0_TESTE(3*K-1,L)=S_ENR_0_TESTE(3*K-1,L)+PHI1D_CT(K)*S_ENR_0(2,L)
                        S_ENR_0_TESTE(3*K,L)=S_ENR_0_TESTE(3*K,L)+PHI1D_CT(K)*S_ENR_0(3,L)
                    ENDDO
                ENDDO
!               ------------------------------------------------------------------------------------------------------
!               COMPUTING THE SST TERMS	        
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
!
                JC0=VJCS
                JC1(:)=DVJCS(:,1)*DCOS(THETA)+DVJCS(:,2)*DSIN(THETA)            

                AUX_CTE=0.D0
                AUX1=0.D0
                AUX2=0.D0
                AUX3=0.D0
                AUX4=0.D0
                AUX_CTE=(S_3/RHO)*C2*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*((THETAF(III)-THETAI(III))/2)*(RHO_BOUND/2.D0)
                DO II=1,3  
			        DO JJ=1,3 
			            DO KK=1,3                             
                            AUX1=(3*Nu/A**2)*(JC0(II)*Va(JJ)*Va(KK)+JC0(JJ)*Va(II)*Va(KK))   ! F0(teta)
!
                            AUX2=3*JC0(KK)*Va(II)*Va(JJ)/A**2
                            AUX2=AUX2+JC0(JJ)*KR(II,KK)+JC0(II)*KR(JJ,KK)
                            AUX2=(1-2*Nu)*AUX2  ! G0(teta) parte 1
!                            
                            AUX3=-(1-4*Nu)*JC0(KK)*KR(II,JJ) ! G0(teta) parte 2
!                            
                            AUX4=AUX1+AUX2+AUX3 !H0(teta)	
					        S_ENR_AUX1(KK,II,JJ)=AUX4*AUX_CTE
                        ENDDO
			        ENDDO
                ENDDO	        
		        
!		        INCLUDING ETA (II) SOURCE POINT 
                S_ENR_AUX2=0.D0
                DO KK=1,3
                    DO JJ=1,3
                        DO II=1,3
                            S_ENR_AUX2(KK,JJ)=S_ENR_AUX2(KK,JJ)+S_ENR_AUX1(KK,II,JJ)*ETAS(II)
                        ENDDO
                    ENDDO
                ENDDO
!               ------------------------------------------------------------------------------------------------------
!               Z_FUN DERIVATIVE AT FIELD POINT
                IF(QSI_REF.EQ.1)THEN
                    Z_FUN_D=Z_AUX_DAL1*DCOS(THETA)
                    Z_FUN_TEST=PSI_TIP_DERIV*DCOS(THETA)
                ELSE
                    Z_FUN_D=Z_AUX_DAL1*DSIN(THETA)
                    Z_FUN_TEST=PSI_TIP_DERIV*DSIN(THETA)
                ENDIF
!        		MULTIPLYING BY DERIVATIVE OF ENRICHED FUNCTIONS
                Z_FUN_TEST=Z_FUN_TEST*RHO+PSI_SP
                CALL DMRRRR (3,3,S_ENR_AUX2,3,3,3,Z_FUN_D,3,3,3,S_ENR_1,3)
                !WRITE(11,*)S_ENR_1(2,1),QSI1D_AL,QSI1D_CT,";"
		        !WRITE(9,*)Z_FUN_D(2,1)*RHO,QSI1D_AL,QSI1D_CT,";"
                !WRITE(7,*)Z_FUN_TEST(2,1),QSI1D_AL,QSI1D_CT,";"		          	             	
!               INTERPOLATION TO CRACK TIPS PHI1DS_CT=PHI1D(KSI_SOURCE POINT)                          
                DO K=1,N_K
                    DO L=1,3
                        S_ENR(3*K-2,L)=S_ENR(3*K-2,L)-PHI1DS_CT(K)*S_ENR_1(1,L)
                        S_ENR(3*K-1,L)=S_ENR(3*K-1,L)-PHI1DS_CT(K)*S_ENR_1(2,L)
                        S_ENR(3*K,L)=S_ENR(3*K,L)-PHI1DS_CT(K)*S_ENR_1(3,L)
                        S_ENR_1_TESTE(3*K-2,L)=S_ENR_1_TESTE(3*K-2,L)+PHI1DS_CT(K)*S_ENR_1(1,L)
                        S_ENR_1_TESTE(3*K-1,L)=S_ENR_1_TESTE(3*K-1,L)+PHI1DS_CT(K)*S_ENR_1(2,L)
                        S_ENR_1_TESTE(3*K,L)=S_ENR_1_TESTE(3*K,L)+PHI1DS_CT(K)*S_ENR_1(3,L)
                    ENDDO
                ENDDO
            ENDDO
!           ------------------------------------------------------------------------------------------------------
!           COMPUTING SEMIANALYTICALLY THE VPC
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
            !S_2=-3.D0*VaVb/(A**5.D0)    !F-1(teta)
            S_3=1/(A**3)
!
            JC0=VJCS
            JC1(:)=DVJCS(:,1)*DCOS(THETA)+DVJCS(:,2)*DSIN(THETA) 
! 
            CPV=DLOG(ABS(RHO_BOUND/Beta))    
!              
            AUX_CTE=0.D0
            AUX_CTE=(S_3*CPV)*C2*QSIW_GAUSS(I,2)*((THETAF(III)-THETAI(III))/2)
            DO II=1,3
			    DO JJ=1,3
			        DO KK=1,3 
                        AUX1=0.D0
                        AUX1=3*Nu*(JC0(II)*Va(JJ)*Va(KK)/A**2+JC0(JJ)*Va(II)*Va(KK)/A**2)   ! F0(teta)
!
                        AUX2=0.D0
                        AUX2=3*JC0(KK)*Va(II)*Va(JJ)/A**2
                        AUX2=AUX2+JC0(JJ)*KR(II,KK)+JC0(II)*KR(JJ,KK)
                        AUX2=(1-2*Nu)*AUX2  ! G0(teta) parte 1
!                            
                        AUX3=0.D0
                        AUX3=-(1-4*Nu)*JC0(KK)*KR(II,JJ) ! G0(teta) parte 2
!                            
                        AUX4=AUX1+AUX2+AUX3 !H0(teta)	
					    S_ENR_AUX1(KK,II,JJ)=AUX4*AUX_CTE
			        ENDDO
			    ENDDO
		    ENDDO    
!		        
!		    INCLUDING ETA (II) SOURCE POINT
            S_ENR_AUX2=0.D0
            DO KK=1,3
                DO JJ=1,3
                    DO II=1,3
                        S_ENR_AUX2(KK,JJ)=S_ENR_AUX2(KK,JJ)+S_ENR_AUX1(KK,II,JJ)*ETAS(II)
                    ENDDO
                ENDDO
            ENDDO
!           ------------------------------------------------------------------------------------------------------
!           WILLIAMS DISPLACEMENT DERIVATIVE AT ADIMENSIONAL DIRECTIONS EVALUATED AT FIELD POINT             
!           Z_FUN DERIVATIVE AT FIELD POINT
            IF(QSI_REF.EQ.1)THEN
                Z_FUN_D=Z_AUX_DAL1*DCOS(THETA)
            ELSE
                Z_FUN_D=Z_AUX_DAL1*DSIN(THETA)
            ENDIF
!           ------------------------------------------------------------------------------------------------------            
!        	MULTIPLYING BY DERIVATIVE ENRICHED FUNCTIONS AND INTERPOLATION TO CRACK TIPS PHIS=PHI(KSI_SOURCE POINT)  
            CALL DMRRRR (3,3,S_ENR_AUX2,3,3,3,Z_FUN_D,3,3,3,S_ENR_1,3)
!
            DO K=1,N_K
                DO L=1,3
                    S_ENR(3*K-2,L)=S_ENR(3*K-2,L)+PHI1DS_CT(K)*S_ENR_1(1,L)
                    S_ENR(3*K-1,L)=S_ENR(3*K-1,L)+PHI1DS_CT(K)*S_ENR_1(2,L)
                    S_ENR(3*K,L)=S_ENR(3*K,L)+PHI1DS_CT(K)*S_ENR_1(3,L)
                    S_ENR_VPC_TESTE(3*K-2,L)=S_ENR_VPC_TESTE(3*K-2,L)+PHI1DS_CT(K)*S_ENR_1(1,L)
                    S_ENR_VPC_TESTE(3*K-1,L)=S_ENR_VPC_TESTE(3*K-1,L)+PHI1DS_CT(K)*S_ENR_1(2,L)
                    S_ENR_VPC_TESTE(3*K,L)=S_ENR_VPC_TESTE(3*K,L)+PHI1DS_CT(K)*S_ENR_1(3,L)
                ENDDO
            ENDDO
!           ------------------------------------------------------------------------------------------------------
	    ENDDO        
    ENDDO  
!
    DEALLOCATE(VALUES)!,QSIW_GAUSS,COEFFICIENTS,PHI,DPHI,D2PHI)  
!	
!    DO I=6,11
        !WRITE(I,*)"];"
        !CLOSE(I)
!    ENDDO
!
    END SUBROUTINE
    