    SUBROUTINE XBEM_K_QUADRILATERAL_HP(NODE,ELEM,ELEMS,DOMAIN,S_ENR,N)
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
	INTEGER::NODE,NODES,ELEM,ELEMS,DOMAIN,I,J,K,L,M,NEN,N,N_CRACKTIP_PAIR,POS_CT_PAIR,QSI_REF,N_GAUSS,&
    POS_CT(N_K),II,JJ,KK,POS_ELEM,ELEM_CRACKTIP,NEN_CT_ELEM,ELEM_CT_S,NEN_CT_S,NEN_TEST
!
    REAL*8::PI,X,Y,Z,DX,DY,DZ,R,DRDN,VJC(3),DVJC(3,2),JC,Mu,Nu,C1,C2,AUX1,AUX2,AUX3,DR(3),&
	ETA(3),KR(3,3),PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),&
    VALUES[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:),&
    S_ENR(3*N_K,3),Z_FUN(3,3),Z_AUX(3,3),S_ENR_AUX1(3,3,3),S_ENR_AUX2(3,3),COORDS_TIP(3),&
    VALUES_K_ALIGN(N_K,3),COORDS_CP_PROJECTION(3),ETAS(3),SIGNAL,QSI_PNT(2),VALUES_ELEM_CRACKTIP[ALLOCATABLE](:,:)
	REAL*8::QSIW_GAUSS[ALLOCATABLE](:,:),QSI1D_CT,QSI1D_AL,QSI1D_AL_CP(N_K),COEFFICIENTS_K_ALIGN(N_K)
    REAL*8::X_TIP,Y_TIP,Z_TIP,PHI1D_CT(N_K),PHI1D_AL(N_K),R_TIP,THETA_TIP,ROTATION_MATRIX(3,3),PSI_TIP(3,3),&
    R_TIP_COLLOC,THETA_TIP_COLLOC,PSI_TIP_COLLOCPNT(3,3,N_K),&
    DR_DQSI1_CT(3),DR_DQSI2_CT(3),DS_ENR(3,3)
    REAL*8::VALUES_CT_S[ALLOCATABLE](:,:),VALUES_TEST[ALLOCATABLE](:,:),X_CT_S,Y_CT_S,Z_CT_S,&
    X_TEST,Y_TEST,Z_TEST,AUX4,DIST_AUX(N_K)
!
	PI=DACOS(-1.D0)
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
	N_GAUSS=N
    S_ENR=0.D0
    DS_ENR=0.D0
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
!   SOURCE POINT NORMAL OUTWARD MATRIX
    NEN=ELEM_TYPE(ELEMS)*ORDER_ELEM(ELEMS)+(ELEM_TYPE(ELEMS)-3)*(ORDER_ELEM(ELEMS)-1)*POL_FAMILY(ELEMS)
	ALLOCATE(VALUES(NEN,3))
	NODES=0
	DO I=1,NEN
		IF(COLLOCPOINTS_CONNECTIVITY(ELEMS,I).EQ.NODE) THEN
			NODES=I
		ENDIF
    ENDDO	
	DO I=1,NEN
		VALUES(I,1)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,I),1)
		VALUES(I,2)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,I),2)
		VALUES(I,3)=COORD_NODES(NODES_CONNECTIVITY(ELEMS,I),3)
    ENDDO
	CALL NORMAL_OUTWARD(NEN,VALUES,QSI(ELEMS,NODES,1),QSI(ELEMS,NODES,2),ETAS,VJC,DVJC,JC)
    DEALLOCATE(VALUES)	
!   ------------------------------------------------------------------------------------------------------
!   ELEMENT NODES COORDINATES
    NEN=ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM)
	ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2))
	CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM,NEN,COEFFICIENTS)
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
!   NUMERICAL INTEGRATION
    DO I=1,N_GAUSS
        DO J=1,N_GAUSS
!           FUNDAMENTAL SOLUTION S*
		    DR=0.D0
		    DRDN=0.D0
		    ETA=0.D0
!
            CALL SHAPE_FUNCTIONS(1,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!        
		    DX=X-COORD_COLLOCPOINTS(NODE,1)
		    DY=Y-COORD_COLLOCPOINTS(NODE,2)
		    DZ=Z-COORD_COLLOCPOINTS(NODE,3)
		    R=DSQRT(DX*DX+DY*DY+DZ*DZ)
		    DR(1)=DX/R
		    DR(2)=DY/R
		    DR(3)=DZ/R
!
		    CALL NORMAL_OUTWARD(NEN,VALUES,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),ETA,VJC,DVJC,JC)
!
		    DRDN=DR(1)*ETA(1)+DR(2)*ETA(2)+DR(3)*ETA(3)
!
            AUX1=0.D0
            AUX2=0.D0
            AUX3=0.D0
            AUX1=(1.D0/R**3)*C2*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*JC
		    DO KK=1,3
			    DO II=1,3
			        DO JJ=1,3
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
					    S_ENR_AUX1(KK,II,JJ)=AUX1*AUX2
				    ENDDO
			    ENDDO
		    ENDDO
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
!           ENRICHED FUNCTION Z_FUN 
!           CALCULATING CRACK TIP COORDINATES ALIGNED TO GAUSS POINT
            IF(QSI_REF.EQ.1)THEN
                QSI1D_CT=QSIW_GAUSS(J,1)
                QSI1D_AL=QSIW_GAUSS(I,1)  
            ELSE
                QSI1D_CT=QSIW_GAUSS(I,1)
                QSI1D_AL=QSIW_GAUSS(J,1)
            ENDIF
!           
            CALL SHAPE_FUNCTIONS_1D(2,QSI1D_CT,N_K-1,VALUES_K_ALIGN,COEFFICIENTS_K_ALIGN,X_TIP,Y_TIP,Z_TIP,PHI1D_CT,DPHI)           
!           CRACK TIP COORDINATES ALIGNED TO FIELD POINT
            COORDS_TIP=0.D0
            DO K=1,N_K
                COORDS_TIP(1)=COORDS_TIP(1)+PHI1D_CT(K)*COORD_CRACKTIP(POS_CT(K),1)
                COORDS_TIP(2)=COORDS_TIP(2)+PHI1D_CT(K)*COORD_CRACKTIP(POS_CT(K),2)
                COORDS_TIP(3)=COORDS_TIP(3)+PHI1D_CT(K)*COORD_CRACKTIP(POS_CT(K),3)
            ENDDO
!           ------------------------------------------------------------------------------------------------------
!           CRACK TIP ROTATION MATRIX
!            IF(QSI_REF.EQ.1)THEN
!	            CALL SHAPE_FUNCTIONS(3,CT_PAIR_QSI(POS_CT_PAIR),QSI1D_CT,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!	            CALL NORMAL_OUTWARD(NEN,VALUES,CT_PAIR_QSI(POS_CT_PAIR),QSI1D_CT,ETA,VJC,DVJC,JC)
!            ELSE
!	            CALL SHAPE_FUNCTIONS(3,QSI1D_CT,CT_PAIR_QSI(POS_CT_PAIR),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!	            CALL NORMAL_OUTWARD(NEN,VALUES,QSI1D_CT,CT_PAIR_QSI(POS_CT_PAIR),ETA,VJC,DVJC,JC)
!            ENDIF
!            DR_DQSI1_CT=0.D0
!            DR_DQSI2_CT=0.D0
!            DO K=1,3
!                DO L=1,NEN
!                    DR_DQSI1_CT(K)=DR_DQSI1_CT(K)+DPHI(L,1)*VALUES(L,K)
!                    DR_DQSI2_CT(K)=DR_DQSI2_CT(K)+DPHI(L,2)*VALUES(L,K)
!                ENDDO
!            ENDDO
!            DR_DQSI1_CT=DR_DQSI1_CT/(DSQRT(DR_DQSI1_CT(1)**2+DR_DQSI1_CT(2)**2+DR_DQSI1_CT(3)**2))
!            DR_DQSI2_CT=DR_DQSI2_CT/(DSQRT(DR_DQSI2_CT(1)**2+DR_DQSI2_CT(2)**2+DR_DQSI2_CT(3)**2))
!!           CHANGE DIRECTION OF NORMAL OUTWARD BECAUSE DISPLACEMENT EQUATION
!            IF(DUAL_BEM(COLLOCPOINTS_CONNECTIVITY(ELEM,1)).EQ."S")THEN
!                ETA=-ETA
!            !    DR_DQSI1_CT=-DR_DQSI1_CT
!            !    DR_DQSI2_CT=-DR_DQSI2_CT
!            ENDIF
!            IF(QSI_REF.EQ.1)THEN
!                ROTATION_MATRIX(:,1)=DR_DQSI1_CT
!                ROTATION_MATRIX(:,2)=ETA
!                ROTATION_MATRIX(:,3)=DR_DQSI2_CT
!            ELSE
!                ROTATION_MATRIX(:,1)=DR_DQSI2_CT
!                ROTATION_MATRIX(:,2)=ETA
!                ROTATION_MATRIX(:,3)=DR_DQSI1_CT
!            ENDIF
!           ------------------------------------------------------------------------------------------------------
!           WILLIAMS DISPL FUNCTIONS OF EACH COLLOCATION POINT ALIGNED TO FIELD POINT
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
!           ------------------------------------------------------------------------------------------------------
!           WILLIAMS FUNCTIONS OF FIELD POINT AND CRACK TIP
            DX=X-COORDS_TIP(1)
		    DY=Y-COORDS_TIP(2)
		    DZ=Z-COORDS_TIP(3)
            R_TIP=DSQRT(DX*DX+DY*DY+DZ*DZ)
            THETA_TIP=SIGNAL*ACOS((ROTATION_MATRIX(1,1)*DX+ROTATION_MATRIX(2,1)*DY+ROTATION_MATRIX(3,1)*DZ)/R_TIP)
            CALL WILLIAMS_DISPL_FUNCTIONS(R_TIP,THETA_TIP,Mu,Nu,PSI_TIP)
!           ------------------------------------------------------------------------------------------------------
!           ENRICHED Z FUNCTION
            Z_AUX=0.D0
!           VALUES_K_ALIGN IS USELESS IN FUNCTION BELOW
            CALL SHAPE_FUNCTIONS_1D(2,QSI1D_AL,N_K-1,VALUES_K_ALIGN,QSI1D_AL_CP,X_TIP,Y_TIP,Z_TIP,PHI1D_AL,DPHI)
            DO K=1,N_K
                Z_AUX=Z_AUX+PHI1D_AL(K)*(PSI_TIP-PSI_TIP_COLLOCPNT(:,:,K))
            ENDDO   
!
            CALL DMRRRR (3,3,ROTATION_MATRIX,3,3,3,Z_AUX,3,3,3,Z_FUN,3)
!           ------------------------------------------------------------------------------------------------------
!           S*Z AND INTERPOLATION TO CRACK TIPS CONSIDERED
		    CALL DMRRRR (3,3,S_ENR_AUX2,3,3,3,Z_FUN,3,3,3,DS_ENR,3)
!
            DO K=1,N_K
                DO L=1,3
                    S_ENR(3*K-2,L)=S_ENR(3*K-2,L)+PHI1D_CT(K)*DS_ENR(1,L)
                    S_ENR(3*K-1,L)=S_ENR(3*K-1,L)+PHI1D_CT(K)*DS_ENR(2,L)
                    S_ENR(3*K,L)=S_ENR(3*K,L)+PHI1D_CT(K)*DS_ENR(3,L)
                ENDDO
            ENDDO
!           ------------------------------------------------------------------------------------------------------
	    ENDDO
    ENDDO  
!
    DEALLOCATE(VALUES)!,QSIW_GAUSS,COEFFICIENTS,PHI,DPHI,D2PHI)
!	
    END SUBROUTINE