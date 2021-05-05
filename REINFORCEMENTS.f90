SUBROUTINE REINFORCEMENTS_COUPLING
    
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE REINFORCEMENTS
    USE OMP_LIB
    USE VISCO
        
    IMPLICIT NONE
    
    INTERFACE 
        SUBROUTINE INTERNAL_POINTS_MATRIX(H_IC,G_IC,COORDS,PROP,DOMAIN)
            INTEGER::DOMAIN    
            REAL*8,DIMENSION(:,:)::H_IC,G_IC
            REAL*8::COORDS(3),PROP(5)
        END SUBROUTINE
    END INTERFACE
    
    INTEGER::I,J,K,II,JJ,KK,M,P_INT_MEF,NDI_MEF,NOF,NUM_LEI,cont,NIAUX[ALLOCATABLE](:,:),NEN,P_INTEGRACAO,N_MATERIALS,&
        IOS
    
    INTEGER,DIMENSION(:),ALLOCATABLE::MATERIAL,MODEL_MATERIALS
    
    REAL*8::AUX1,AUX2,AUX3,AUX4,X,Y,Z,VET_PROP(5),VET_COORD(3),ETA(3),VJC(3),DVJC(3,2),JC,MAT_ETAF(3,9),&
        Nu,Mu,JF,DR0(3),GAMA,K1,K2,K3,K4,EVE,EE,DIST1(5)
    
    REAL*8,PARAMETER::PI=DACOS(-1.D0)
            
    REAL*8,DIMENSION(:),ALLOCATABLE::D_PHI_MEF,PHI
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::GLOCAL,GLOCALAUX,VALORES,KLOCAL,DG_CF,QSIW,DG_AUX_FF,&
        PROP_MATERIALS,PROP_EP_MATERIALS
    
    LOGICAL::INPUTFILE,LINE_HAS_ZEROS,INT_2D
    
    !CALL OMP_SET_NUM_THREADS(24)
    
    !WRITE(*,*)FILE_REINFORCEMENTS
    
    ! OPENING REINFORCEMENTS FILE
    INQUIRE(FILE=TRIM('Input_data\'//FILE_REINFORCEMENTS),EXIST=INPUTFILE) 
    IF(INPUTFILE)THEN
        OPEN(1,FILE='Input_data\'//FILE_REINFORCEMENTS,STATUS='UNKNOWN')
	ELSE
        WRITE(*,*)'THE REINFORCEMENTS FILE DOES NOT EXIST'
        WRITE(*,*)'VERIFY THE INPUT FOLDER'
        WRITE(*,*)'THE PROGRAM IS FINISHED'
        PAUSE
        STOP
    ENDIF
    
    ! DEFINF NUMBER OF INTEGRATION POINTS TO INTEGRATE SINGULAR FE
    P_INTEGRACAO=INT(0.6*N_GAUSSPOINTS)
    IF (P_INTEGRACAO.LT.30) P_INTEGRACAO=30
    ALLOCATE(QSIW(P_INTEGRACAO,2))
    CALL GAUSS_POINTS(P_INTEGRACAO,QSIW)
    ! ------------------------------------------------------------------------------------
    
    !
    !READING REINFORCEMENTS FILE ______________________________________________________________________________________________________
    !
    READ(1,*)!**************** BEM_FEM DADA DEFINITIONS *****************
    READ(1,*)!
    READ(1,*)!WARNING: RESPECT THIS FILE FORMAT.THE SPACE BETWEEN DATA
    READ(1,*)!AND THE NUMBER FORMAT
	READ(1,*)!
	READ(1,*)!ARE THERE REINFORCEMENTS INTO THE DOMAIN (Y FOR YES AND N FOR NOT):
	READ(1,*)EXISTANCE_REINFORCEMENTS
	READ(1,*)! 
	READ(1,*)!METHOD UTILIZED TO MODEL THE REINFORCEMENTS (BEM1D OR FEM):
	READ(1,*)REINFORCEMENTS_METHOD
	READ(1,*)! 
    READ(1,*)!ARE THE REINFORCEMENTS ELASTOPLASTIC? (Y OR N)
	READ(1,*)REINFORCEMENTS_ELASTOPLASTIC
	READ(1,*)!
    READ(1,*)!ARE THERE BOND-SLIP EFFECTS? (Y OR N)
    READ(1,*)REINFORCEMENTS_SLIP
    READ(1,*)
    
    IF(EXISTANCE_REINFORCEMENTS.EQ.'N') THEN
        GO TO 1
    ENDIF
    
    READ(1,*)!NUMBER OF FE NODES:
	READ(1,*)N_NOS_MEF
	READ(1,*)!
	READ(1,*)!NUMBER OF FE ELEMENTS:
	READ(1,*)N_ELEMENTOS_MEF
	READ(1,*)!
    READ(1,*)!NUMBER OF DIFFERENT MATERIALS:
    READ(1,*)N_MATERIALS
    
    !STARTING VARIABLES ------------
    ALLOCATE(COORDNOS_MEF(N_NOS_MEF,3),ORDEM_MEF(N_ELEMENTOS_MEF),ND_MEF(N_ELEMENTOS_MEF),KODENODUPLO_MEF(N_NOS_MEF),&
        MP_MEF(N_ELEMENTOS_MEF,2),N_INT_P_MEF(N_ELEMENTOS_MEF),U_MEF(3*N_NOS_MEF),P_MEF(3*N_NOS_MEF),MATERIAL(N_ELEMENTOS_MEF),&
        PROP_MATERIALS(N_MATERIALS,4),MODEL_MATERIALS(N_MATERIALS)) 
    IF (VISCO_ANALYSIS) ALLOCATE (U_ANTERIOR_MEF(3*N_NOS_MEF),P_ANTERIOR_MEF(3*N_NOS_MEF))
    COORDNOS_MEF=0.D0   
	ORDEM_MEF=0
	ND_MEF=0
	KODENODUPLO_MEF=0
	MP_MEF=0.D0
    !-------------------------------
    
    READ(1,*)!
    READ(1,*)! NODAL COORDINATES:
	READ(1,*)! NODE		X		Y       z
    DO I=1,N_NOS_MEF
        READ(1,*)J,COORDNOS_MEF(J,1),COORDNOS_MEF(J,2),COORDNOS_MEF(J,3)
    ENDDO
    READ(1,*)!
	READ(1,*)! INFORM THE POLYNOMINAL APPROXIMATION FOR THE ALL FE
	READ(1,*)! AND THE DOMAIN THAT THEY BELONG
	READ(1,*)! ELEMENT		DEGREE		DOMAIN  MATERIAL
    DO I=1,N_ELEMENTOS_MEF
        READ(1,*)J,ORDEM_MEF(J),ND_MEF(J),MATERIAL(J)
    ENDDO
    
    !FIDING MAX ELEMENT ORDER ------
    MAXORDER_MEF=0
	DO I=1,N_ELEMENTOS_MEF
		IF(ORDEM_MEF(I).GT.MAXORDER_MEF) THEN
			MAXORDER_MEF=ORDEM_MEF(I)
		ENDIF
    ENDDO
    !-------------------------------
    
    !STARTING VARIABLES ------------
    I=MAXORDER_MEF+1
    ALLOCATE(CONECTI_MEF(N_ELEMENTOS_MEF,I),QSI_MEF(N_ELEMENTOS_MEF,I),KKGLOBAL(N_ELEMENTOS_MEF,(3*I),(3*I)),&
    GGGLOBAL(N_ELEMENTOS_MEF,(3*I),(3*I)))
    CONECTI_MEF=0 
	QSI_MEF=0.D0 
	KKGLOBAL=0.D0
	GGGLOBAL=0.D0
    
    JJ=3*N_COLLOCPOINTS+6*N_NOS_MEF
    ALLOCATE(A_GLOBAL_ELASTIC(JJ,JJ))
    A_GLOBAL_ELASTIC=0.0D0
    II=3*N_COLLOCPOINTS+3*N_NOS_MEF
    JJ=3*N_COLLOCPOINTS
    KGLOBAL=>A_GLOBAL_ELASTIC(II+1:II+3*N_NOS_MEF,JJ+1:JJ+3*N_NOS_MEF)
    
    II=3*N_COLLOCPOINTS+3*N_NOS_MEF
    JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
    GGLOBAL=>A_GLOBAL_ELASTIC(II+1:II+3*N_NOS_MEF,JJ+1:JJ+3*N_NOS_MEF)
        
    READ(1,*)!
	READ(1,*)! INFORM THE CONECTIVITY FOR ALL FE:
    DO I=1,N_ELEMENTOS_MEF
		READ(1,*)! ELEMENT		NODES
		DO J=1,(ORDEM_MEF(I)+1)
			READ(1,*)II,CONECTI_MEF(II,J)
		ENDDO
    ENDDO
    
    !DEFINING THE DOUBLE NODES INTO THE FE MESH -------------
    N_INTERFACE_MEF=0
	DO I=1,N_ELEMENTOS_MEF
        DO J=1,N_ELEMENTOS_MEF
			IF(I.GT.J) THEN
      			DO K=1,(ORDEM_MEF(I)+1)
      				DO II=1,(ORDEM_MEF(J)+1)
      					IF(CONECTI_MEF(I,K).NE.CONECTI_MEF(J,II)) THEN
	  					    AUX1=COORDNOS_MEF(CONECTI_MEF(I,K),1)-COORDNOS_MEF(CONECTI_MEF(J,II),1)
                            AUX2=COORDNOS_MEF(CONECTI_MEF(I,K),2)-COORDNOS_MEF(CONECTI_MEF(J,II),2)
                            AUX3=COORDNOS_MEF(CONECTI_MEF(I,K),3)-COORDNOS_MEF(CONECTI_MEF(J,II),3)
                            AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
                            IF(AUX4.LE.TOLER) THEN
		    						    				    
                                IF((K.EQ.1).OR.(K.EQ.(ORDEM_MEF(I)+1))) THEN
                                    KODENODUPLO_MEF(CONECTI_MEF(I,K))=1
                                ENDIF
                                IF((II.EQ.1).OR.(II.EQ.(ORDEM_MEF(J)+1))) THEN
                                    KODENODUPLO_MEF(CONECTI_MEF(J,II))=1
                                ENDIF
      							      							
                                IF(ND_MEF(I) .NE. ND_MEF(J)) THEN
                                    IF(N_INTERFACE_MEF .EQ. 0) THEN
                                        N_INTERFACE_MEF = N_INTERFACE_MEF + 1
                                        ALLOCATE(NI_MEF(1,2))
                                        NI_MEF(1,1)=CONECTI_MEF(J,II)
                                        NI_MEF(1,2)=CONECTI_MEF(I,K)
                                    ELSE
                                        ALLOCATE(NIAUX(N_INTERFACE_MEF,2))
							            NIAUX = NI_MEF
                                        DEALLOCATE(NI_MEF)
                                        N_INTERFACE_MEF = N_INTERFACE_MEF + 1
                                        ALLOCATE(NI_MEF(N_INTERFACE_MEF,2))
                                        DO KK=1,N_INTERFACE_MEF-1
                                            NI_MEF(KK,1)=NIAUX(KK,1)
                                            NI_MEF(KK,2)=NIAUX(KK,2)
                                        ENDDO
                                        NI_MEF(N_INTERFACE_MEF,1)=CONECTI_MEF(J,II)
                                        NI_MEF(N_INTERFACE_MEF,2)=CONECTI_MEF(I,K)
                                        DEALLOCATE(NIAUX)
                                    ENDIF     
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDDO
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    ALLOCATE(K_CONECT(N_INTERFACE_MEF,6,6),G_CONECT(N_INTERFACE_MEF,6,6))
    ! -------------------------------------------------------
    
    !COMPUTING THE QSI COORDINATES FOR EACH FE COLOCATION POINT ---------
    DO I=1,N_ELEMENTOS_MEF
		AUX1=2.D0/ORDEM_MEF(I)
		DO J=1,(ORDEM_MEF(I)+1)
			IF((J.NE.1).AND.(J.NE.(ORDEM_MEF(I)+1))) THEN
				QSI_MEF(I,J)=-1.D0+AUX1*(J-1)
			ELSE
				IF(KODENODUPLO_MEF(CONECTI_MEF(I,J)).EQ.0) THEN
					QSI_MEF(I,J)=-1.D0+AUX1*(J-1)
				ELSE
					IF(J.EQ.1) THEN
						!QSI_MEF(I,J)=-1.D0+((AUX1)/(4.D0))
                        IF (ORDEM_MEF(I).GE.2) THEN
                            QSI_MEF(I,J)=-1.D0+((AUX1)/(2.D0))
                        ELSE
                            QSI_MEF(I,J)=-1.D0+((AUX1)/(3.D0))
                        ENDIF
                    ELSE
						!QSI_MEF(I,J)=1.D0-((AUX1)/(4.D0))
                        IF (ORDEM_MEF(I).GE.2) THEN
                            QSI_MEF(I,J)=1.D0-((AUX1)/(2.D0))
                        ELSE
                            QSI_MEF(I,J)=1.D0-((AUX1)/(3.D0))
                        ENDIF
					ENDIF
				ENDIF
			ENDIF
        ENDDO
        
        ! COLOCANDO DESCONTINUOS NOS DE EXTREMIDADE DE FIBRA
        II=0;
        JJ=0;
        DO J=1,N_ELEMENTOS_MEF
            IF (J.NE.I) THEN
                DO K=1,ORDEM_MEF(J)+1
                    AUX1=COORDNOS_MEF(CONECTI_MEF(J,K),1)-COORDNOS_MEF(CONECTI_MEF(I,1),1)
                    AUX2=COORDNOS_MEF(CONECTI_MEF(J,K),2)-COORDNOS_MEF(CONECTI_MEF(I,1),2)
                    AUX3=COORDNOS_MEF(CONECTI_MEF(J,K),3)-COORDNOS_MEF(CONECTI_MEF(I,1),3)
                    AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
                    IF (DABS(AUX4).LE.0.00000010D0) II=II+1
                
                    AUX1=COORDNOS_MEF(CONECTI_MEF(J,K),1)-COORDNOS_MEF(CONECTI_MEF(I,ORDEM_MEF(I)+1),1)
                    AUX2=COORDNOS_MEF(CONECTI_MEF(J,K),2)-COORDNOS_MEF(CONECTI_MEF(I,ORDEM_MEF(I)+1),2)
                    AUX3=COORDNOS_MEF(CONECTI_MEF(J,K),3)-COORDNOS_MEF(CONECTI_MEF(I,ORDEM_MEF(I)+1),3)
                    AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
                    IF (DABS(AUX4).LE.0.00000010D0) JJ=JJ+1
                ENDDO
            ENDIF
        ENDDO
        
        IF (KODENODUPLO_MEF(CONECTI_MEF(I,1)).EQ.0) THEN
            IF (ORDEM_MEF(I).GE.2) THEN
                IF (II.EQ.0) QSI_MEF(I,1)=QSI_MEF(I,1)+(QSI_MEF(I,2)-QSI_MEF(I,1))/2
            ELSE
                IF (II.EQ.0) QSI_MEF(I,1)=QSI_MEF(I,1)+(QSI_MEF(I,2)-QSI_MEF(I,1))/3
            ENDIF 
        ENDIF
        IF (KODENODUPLO_MEF(CONECTI_MEF(I,ORDEM_MEF(I)+1)).EQ.0) THEN
            IF (ORDEM_MEF(I).GE.2) THEN
                IF (JJ.EQ.0) QSI_MEF(I,ORDEM_MEF(I)+1)=QSI_MEF(I,ORDEM_MEF(I)+1)-(QSI_MEF(I,2)-QSI_MEF(I,1))/2
            ELSE
                IF (JJ.EQ.0) QSI_MEF(I,ORDEM_MEF(I)+1)=QSI_MEF(I,ORDEM_MEF(I)+1)-(QSI_MEF(I,2)-QSI_MEF(I,1))/3
            ENDIF
        ENDIF
                
    ENDDO
    !--------------------------------------------------------------------
    
    VISCO_ANALYSIS_REINF=.FALSE.
    EP_ANALYSIS_REINF=.FALSE.
    SLIP_REINF=.FALSE.
    
    READ(1,*)!
	READ(1,*)!MATERIAL FE PROPERTIES
    READ(1,*)!MATERIAL     YOUNG MODULUS         CROSS SECTION AREA      VISCOELASTICITY      MOD.BOLTZMANN     LEI CONST(1-KELVIN,2-BOLTZMANN,3-MAXWELL,4-HOOKE):
    DO I=1,N_MATERIALS
	    !READ(1,*)II,MP_MEF(II,1),MP_MEF(II,2)
        READ(1,*)II,PROP_MATERIALS(II,1),PROP_MATERIALS(II,2),PROP_MATERIALS(II,3),PROP_MATERIALS(II,4),MODEL_MATERIALS(II)
        !DEFINING VISCO REINFORCEMENTS
        IF (MODEL_MATERIALS(II).NE.4) VISCO_ANALYSIS_REINF=.TRUE.
    ENDDO
    
    IF (REINFORCEMENTS_ELASTOPLASTIC .EQ. 'Y') THEN
        EP_ANALYSIS_REINF=.TRUE.
	    ALLOCATE(PROP_EP_MATERIALS(N_MATERIALS,2))
	    READ(1,*)!
	    READ(1,*)!IF ELASTOPLASTIC REINFORCEMENTS, GIVE PLASTIC PROPERTIES
	    READ(1,*)!MATERIAL    	SIGMA Y				K
	    DO I=1,N_MATERIALS
	        READ(1,*)II,PROP_EP_MATERIALS(II,1),PROP_EP_MATERIALS(II,2)
	    ENDDO
	    READ(1,*)!
	    READ(1,*)!NUMBER OF LOAD STEPS
	    READ(1,*)N_PASSOS
	    READ(1,*)!
	    READ(1,*)!TOLERANCE FOR ELASTOPLASTIC ANALYSIS
	    READ(1,*)EP_TOL
    ENDIF
        
    IF (REINFORCEMENTS_SLIP .EQ. 'Y') THEN
        
        SLIP_REINF=.TRUE.
        CALL READ_SLIP_PARAMETERS
        
        ! Initiating variables
        ALLOCATE(NODAL_FORCES_PLACE(N_NODAL_FORCES,2),NODAL_FORCES_VALUE(N_NODAL_FORCES),&
            F_APP_DISP(N_NODAL_FORCES))
        F_APP_DISP=0.0D0
        
        READ(1,*)!
        READ(1,*)!N_NO      DIR     DISPLACEMENT   
        DO I=1,N_NODAL_FORCES
            READ(1,*)NODAL_FORCES_PLACE(I,1),NODAL_FORCES_PLACE(I,2),NODAL_FORCES_VALUE(I)
        ENDDO
        
    ENDIF

    ! STARTING REINFORCEMENTS VARIABLES --------    
    IF (VISCO_ANALYSIS) THEN
        II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
        ALLOCATE(B_GLOBAL(JJ,II+6*N_NOS_MEF+6*N_COLLOCPOINTS))
        B_GLOBAL=0.0D0
        IF (VISCO_ANALYSIS_REINF) THEN  
            !ALLOCATE(TIME_MODEL_REINF(N_ELEMENTOS_MEF),MP_MEF_VISCO(N_ELEMENTOS_MEF,2),&
            !K_KGLOBAL((3*N_NOS_MEF),(3*N_NOS_MEF)),G_GGLOBAL((3*N_NOS_MEF),(3*N_NOS_MEF)))
            !K_KGLOBAL=0.D0
            !G_GGLOBAL=0.D0
            ALLOCATE(TIME_MODEL_REINF(N_ELEMENTOS_MEF),MP_MEF_VISCO(N_ELEMENTOS_MEF,2))
            
            II=3*N_COLLOCPOINTS+3*N_NOS_MEF
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+3*N_COLLOCPOINTS
            K_KGLOBAL=>B_GLOBAL(II+1:II+3*N_NOS_MEF,JJ+1:JJ+3*N_NOS_MEF)
            
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+6*N_COLLOCPOINTS+3*N_NOS_MEF
            G_GGLOBAL=>B_GLOBAL(II+1:II+3*N_NOS_MEF,JJ+1:JJ+3*N_NOS_MEF)
        ENDIF
    ENDIF
    
    IF (EP_ANALYSIS_REINF) THEN
        ALLOCATE (MP_EP_MEF(N_ELEMENTOS_MEF,3))
    ENDIF
    IF (SLIP_REINF) THEN
        !I=MAXORDER_MEF+1
        !ALLOCATE(KGLOBAL_SLIP(N_NOS_MEF,N_NOS_MEF))
        !KGLOBAL_SLIP=0.0D0      
    ENDIF
    ! ----------------------------------------
    
    ! GIVING FE PROPERTIES:
    DO I=1,N_ELEMENTOS_MEF
        MP_MEF(I,1)=PROP_MATERIALS(MATERIAL(I),1)           !YOUNG
        MP_MEF(I,2)=PROP_MATERIALS(MATERIAL(I),2)           !AREA
        IF (VISCO_ANALYSIS_REINF) THEN
            MP_MEF_VISCO(I,1)=PROP_MATERIALS(MATERIAL(I),3)     !VISCO COEFICIENT (GAMA)
            MP_MEF_VISCO(I,2)=PROP_MATERIALS(MATERIAL(I),4)     !MOD. ELAST. BOLTZMANN
            TIME_MODEL_REINF(I)=MODEL_MATERIALS(MATERIAL(I))    !LEI CONSTITUTIVA
        ENDIF
        IF (EP_ANALYSIS_REINF) THEN
            MP_EP_MEF(I,1)=PROP_EP_MATERIALS(MATERIAL(I),1)     ! Sigma Y
            MP_EP_MEF(I,2)=PROP_EP_MATERIALS(MATERIAL(I),2)     ! Kp
            MP_EP_MEF(I,3)=MP_MEF(I,1)*MP_EP_MEF(I,2)/(MP_MEF(I,1)+MP_EP_MEF(I,2))   ! E_EP
        ENDIF
    ENDDO
   
    ! READING INFORMATION TO GIVE RESULTS: FORCE X DISPLACEMENT
    IF (EP_ANALYSIS_REINF) THEN
        READ(1,*,IOSTAT=IOS)!
        READ(1,*,IOSTAT=IOS)!PLOT NODE HISTORY
        READ(1,*,IOSTAT=IOS)!COORD_X    COORD_Y     COORD_Z
        READ(1,*,IOSTAT=IOS)AUX1,AUX2,AUX3
        IF (IOS.EQ.0) THEN
            NODE_HIST='B'
        ELSE
            NODE_HIST='N'
        ENDIF
        IF (NODE_HIST.EQ.'B') THEN
            DIST1=10.0D0
            DO J=1,N_COLLOCPOINTS
                DIST1(1)=COORD_COLLOCPOINTS(J,1)-AUX1
			    DIST1(2)=COORD_COLLOCPOINTS(J,2)-AUX2
	            DIST1(3)=COORD_COLLOCPOINTS(J,3)-AUX3
			    DIST1(4)=DSQRT(DIST1(1)**2.0D0+DIST1(2)**2.0D0+DIST1(3)**2.0D0)
                IF (DIST1(4).LT.DIST1(5)) THEN
                    DIST1(5)=DIST1(4)
                    N_NODE_HIST=J
                ENDIF
             ENDDO
        ENDIF
    ENDIF
    IF (SLIP_REINF) THEN
        READ(1,*,IOSTAT=IOS)!
        READ(1,*,IOSTAT=IOS)!PLOT NODE HISTORY
        READ(1,*,IOSTAT=IOS)!COORD_X    COORD_Y     COORD_Z
        READ(1,*,IOSTAT=IOS)AUX1,AUX2,AUX3
        IF (IOS.EQ.0) THEN
            NODE_HIST='Y'
        ELSE
            NODE_HIST='N'
        ENDIF
        IF (NODE_HIST.EQ.'Y') THEN
            DIST1=10.0D0
            DO J=1,N_NOS_MEF
                DIST1(1)=COORDNOS_MEF(J,1)-AUX1
			    DIST1(2)=COORDNOS_MEF(J,2)-AUX2
	            DIST1(3)=COORDNOS_MEF(J,3)-AUX3
			    DIST1(4)=DSQRT(DIST1(1)**2.0D0+DIST1(2)**2.0D0+DIST1(3)**2.0D0)
                IF (DIST1(4).LT.DIST1(5)) THEN
                    DIST1(5)=DIST1(4)
                    N_NODE_HIST=J
                ENDIF
             ENDDO
        ENDIF
    ENDIF
    
    CLOSE (1) !END OF READING REINFORCEMENTS FILE ________________________________________________________________________________    
    
    !
    ! CREATING THE FE STIFFNESS MATRIX AND THE LUMPING MATRIX ____________________________________________________________________
    !
    DO I=1,N_ELEMENTOS_MEF
        
        ! STARTING VARIABLES --------------
        ALLOCATE(KLOCAL(ORDEM_MEF(I)+1,ORDEM_MEF(I)+1),GLOCAL(ORDEM_MEF(I)+1,ORDEM_MEF(I)+1))
        
        IF (REINFORCEMENTS_METHOD.EQ.'FEM') THEN ! LOCAL MATRIXES CALCULATION _____________________________________________________
            
            CALL LOCAL_MATRIXES_FEM(I,KLOCAL,GLOCAL)
            !CALL LOCAL_MATRIXES_1DBEM(I,KLOCAL,GLOCAL)
            
        ELSE IF (REINFORCEMENTS_METHOD.EQ.'1DBEM') THEN ! LOCAL MATRIXES CALCULATION WITH 1DBEM __________________________________
            
            CALL LOCAL_MATRIXES_1DBEM(I,KLOCAL,GLOCAL)
            
        ELSE ! END OF CASE 1DBEM _________________________________________________________________________________________________ 
            
            WRITE(*,*)"THE REINFORCEMENT METHOD ISN'T VALID"
            WRITE(*,*)'VERIFY THE VALUE'
            WRITE(*,*)'THE PROGRAM IS FINISHED'
            PAUSE
            STOP
        
        ENDIF ! END OF LOCAL MATRIXES CALCULATION ________________________________________________________________________________
                
        ! MATRIX ROTATION: KLOCAL -> KKGLOBAL AND GLOCAL -> GGGLOBAL
        CALL ROTATE_MATRIXES(I,KLOCAL,GLOCAL)
                
        ! FINDING VISCOELASTIC COEFICIENTS FOR THIS FE
        IF (VISCO_ANALYSIS_REINF) THEN
            NUM_LEI = TIME_MODEL_REINF(I)
            GAMA = MP_MEF_VISCO(I,1)
            Eve = MP_MEF_VISCO(I,2)
            Ee   = MP_MEF(I,1)
            SELECT CASE(NUM_LEI)
            CASE(1)!KELVIN'S MODEL
	          K1= 1.0D0 + GAMA/DELTA_T_VISCO
	          K2= 1.0D0
	          K3= GAMA/DELTA_T_VISCO
	          K4= 0.0D0  
            CASE(2)!BOLTZMANN'S MODEL
              K1= 1.0D0 + GAMA/DELTA_T_VISCO
	          K2= GAMA/DELTA_T_VISCO + (Ee + Eve)/Eve
	          K3= GAMA/DELTA_T_VISCO
	          K4= -GAMA/DELTA_T_VISCO    
            CASE(3)!MAXWELL'S MODEL
              K1= 1.0D0
	          K2= 1.0D0 + DELTA_T_VISCO/GAMA
	          K3= 1.0D0
	          K4=-1.0D0 
            CASE(4)!HOOKE'S MODEL
              K1= 1.0D0
	          K2= 1.0D0
	          K3= 0.0D0
	          K4= 0.0D0  
            ENDSELECT
        ENDIF
        
        !TRANSFER THE LOCAL STIFNESS AND LUMPING FOR THE TOTAL ONES. KKGLOBAL -> KGLOBAL AND GGGLOBAL -> GGLOBAL
        DO J=1,ORDEM_MEF(I)+1
            DO JJ=1,3
                DO K=1,ORDEM_MEF(I)+1
                    
                    DO KK=1,3
                        IF (VISCO_ANALYSIS_REINF) THEN
                            KGLOBAL(3*(CONECTI_MEF(I,J)-1)+JJ,3*(CONECTI_MEF(I,K)-1)+KK)=&
                            KGLOBAL(3*(CONECTI_MEF(I,J)-1)+JJ,3*(CONECTI_MEF(I,K)-1)+KK)+&
                            KKGLOBAL(I,3*(J-1)+JJ,3*(K-1)+KK)*K1
                        
                            GGLOBAL(3*(CONECTI_MEF(I,J)-1)+JJ,3*(CONECTI_MEF(I,K)-1)+KK)=&
                            GGLOBAL(3*(CONECTI_MEF(I,J)-1)+JJ,3*(CONECTI_MEF(I,K)-1)+KK)+&
                            GGGLOBAL(I,3*(J-1)+JJ,3*(K-1)+KK)*K2
                            !------------------------------------------------------------------
                            K_KGLOBAL(3*(CONECTI_MEF(I,J)-1)+JJ,3*(CONECTI_MEF(I,K)-1)+KK)=&
                            K_KGLOBAL(3*(CONECTI_MEF(I,J)-1)+JJ,3*(CONECTI_MEF(I,K)-1)+KK)+&
                            KKGLOBAL(I,3*(J-1)+JJ,3*(K-1)+KK)*K3
                        
                            G_GGLOBAL(3*(CONECTI_MEF(I,J)-1)+JJ,3*(CONECTI_MEF(I,K)-1)+KK)=&
                            G_GGLOBAL(3*(CONECTI_MEF(I,J)-1)+JJ,3*(CONECTI_MEF(I,K)-1)+KK)+&
                            GGGLOBAL(I,3*(J-1)+JJ,3*(K-1)+KK)*K4
                        ELSE
                            KGLOBAL(3*(CONECTI_MEF(I,J)-1)+JJ,3*(CONECTI_MEF(I,K)-1)+KK)=&
                            KGLOBAL(3*(CONECTI_MEF(I,J)-1)+JJ,3*(CONECTI_MEF(I,K)-1)+KK)+&
                            KKGLOBAL(I,3*(J-1)+JJ,3*(K-1)+KK)
                        
                            GGLOBAL(3*(CONECTI_MEF(I,J)-1)+JJ,3*(CONECTI_MEF(I,K)-1)+KK)=&
                            GGLOBAL(3*(CONECTI_MEF(I,J)-1)+JJ,3*(CONECTI_MEF(I,K)-1)+KK)+&
                            GGGLOBAL(I,3*(J-1)+JJ,3*(K-1)+KK)
                        ENDIF
                    ENDDO
                    
                ENDDO
            ENDDO
        ENDDO
        ! ___________________________________________________________________________________________________
        
        DEALLOCATE(GLOCAL,KLOCAL)
        
    ENDDO ! END OF STIFFNES MATRIX AND LUMPING MATRIX ____________________________________________________________________________
    
    !
    !COMPUTING THE COORDINATES FOR ALL THE FE COLLOCATION POINTS _________________________________________________________________
    !
    ALLOCATE (COORDPCOLOC_MEF(N_NOS_MEF,3))
	COORDPCOLOC_MEF=0.D0
    DO I=1,N_NOS_MEF
        !IF(KODENODUPLO_MEF(I).EQ.1) THEN

		    II=0
			DO J=1,N_ELEMENTOS_MEF
				DO K=1,(ORDEM_MEF(J)+1)
					IF(CONECTI_MEF(J,K).EQ.I) THEN
						II=J
						AUX2=QSI_MEF(J,K)
						ALLOCATE(VALORES((ORDEM_MEF(J)+1),3),PHI(ORDEM_MEF(J)+1))
						VALORES=0.D0
						PHI=0.D0
						DO KK=1,(ORDEM_MEF(J)+1)
						    VALORES(KK,1)=COORDNOS_MEF(CONECTI_MEF(J,KK),1)
						    VALORES(KK,2)=COORDNOS_MEF(CONECTI_MEF(J,KK),2)
                            VALORES(KK,3)=COORDNOS_MEF(CONECTI_MEF(J,KK),3)
						ENDDO
						CALL FUNCOES_DE_FORMA_MEC_MEF(5,AUX2,(ORDEM_MEF(J)+1),J,VALORES,X,Y,Z,PHI)
                        COORDPCOLOC_MEF(I,1)=X
			            COORDPCOLOC_MEF(I,2)=Y
                        COORDPCOLOC_MEF(I,3)=Z
        
						GO TO 2
					ENDIF
				ENDDO
            ENDDO
2			DEALLOCATE(VALORES,PHI)
		!ELSE
		!	COORDPCOLOC_MEF(I,1)=COORDNOS_MEF(I,1)
		!	COORDPCOLOC_MEF(I,2)=COORDNOS_MEF(I,2)
  !          COORDPCOLOC_MEF(I,3)=COORDNOS_MEF(I,3)
		!ENDIF
    ENDDO
    !_____________________________________________________________________________________________________________________________
    
    !
    ! CREATING COONECTION ELEMENT MATRIXES _______________________________________________________________________________________
    !
    CALL CONNECTION_ELEMENT
    
    ! STARTING VARIABLES FOR MATRIXES
    ALLOCATE(H_FC(3*N_NOS_MEF,3*N_COLLOCPOINTS),G_FC(3*N_NOS_MEF,3*N_COLLOCPOINTS))
    MAT_ETAF=0.0D0
        
    JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
    G_CF=>A_GLOBAL_ELASTIC(1:3*N_COLLOCPOINTS,JJ+1:JJ+3*N_NOS_MEF)
    
    II=3*N_COLLOCPOINTS
    JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
    G_FF=>A_GLOBAL_ELASTIC(II+1:II+3*N_NOS_MEF,JJ+1:JJ+3*N_NOS_MEF)
    
    IF (VISCO_ANALYSIS) THEN
        !ALLOCATE(HH_FC(3*N_NOS_MEF,3*N_COLLOCPOINTS),GG_FC(3*N_NOS_MEF,3*N_COLLOCPOINTS),&
        !GG_FF(3*N_NOS_MEF,N_NOS_MEF),GG_CF(3*N_COLLOCPOINTS,N_NOS_MEF))
        !GG_FF=0.0D0
        !GG_CF=0.0D0
        
        II=3*N_COLLOCPOINTS
        JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
        HH_FC=>B_GLOBAL(II+1:II+3*N_NOS_MEF,JJ+1:JJ+3*N_COLLOCPOINTS)
        
        II=3*N_COLLOCPOINTS
        JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+3*N_COLLOCPOINTS+3*N_NOS_MEF
        GG_FC=>B_GLOBAL(II+1:II+3*N_NOS_MEF,JJ+1:JJ+3*N_COLLOCPOINTS)
        
        II=0
        JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+6*N_COLLOCPOINTS+3*N_NOS_MEF
        GG_CF=>B_GLOBAL(II+1:II+3*N_COLLOCPOINTS,JJ+1:JJ+3*N_NOS_MEF)
        
        II=3*N_COLLOCPOINTS
        JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+6*N_COLLOCPOINTS+3*N_NOS_MEF
        GG_FF=>B_GLOBAL(II+1:II+3*N_NOS_MEF,JJ+1:JJ+3*N_NOS_MEF)
    ENDIF
        
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(NDI_MEF,NUM_LEI,GAMA,EVE,EE,K1,K2,K3,K4,VET_PROP,VET_COORD,NEN,VALORES,NOF,MAT_ETAF,NU,MU,&
    DG_CF,KK,II,AUX1,AUX2,AUX3,JF,DR0,D_PHI_MEF,DG_AUX_FF,INT_2D)
    
    ! -----------------------------------------------------------------------------------------------------------------------------
    ! CALCULATING THE MATRIXES RELATED TO THE INTERNAL POINTS H_FC AND G_FC -------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------------------------------
    
    !$omp do schedule(dynamic)
    DO I=1,N_NOS_MEF
        
        ! FINDING ELEMENT OF NODE I AND ITS 2D DOMAIN
        NDI_MEF=0
        DO J=1,N_ELEMENTOS_MEF
            DO K=1,ORDEM_MEF(J)+1
                IF(CONECTI_MEF(J,K).EQ.I) THEN
                    NDI_MEF=ND_MEF(J)            
                ENDIF
            ENDDO
        ENDDO
        
        IF (VISCO_ANALYSIS) THEN
            NUM_LEI = TIME_MODEL(NDI_MEF)
            GAMA = EMP_VISCO(NDI_MEF,1)
            Eve = EMP_VISCO(NDI_MEF,2)
            Ee   = EMP(NDI_MEF,1)
            SELECT CASE(NUM_LEI)
            CASE(1)!KELVIN'S MODEL
	          K1= 1.0D0 + GAMA/DELTA_T_VISCO
	          K2= 1.0D0
	          K3= GAMA/DELTA_T_VISCO
	          K4= 0.0D0  
            CASE(2)!BOLTZMANN'S MODEL
              K1= 1.0D0 + GAMA/DELTA_T_VISCO
	          K2= GAMA/DELTA_T_VISCO + (Ee + Eve)/Eve
	          K3= GAMA/DELTA_T_VISCO
	          K4= -GAMA/DELTA_T_VISCO    
            CASE(3)!MAXWELL'S MODEL
              K1= 1.0D0
	          K2= 1.0D0 + DELTA_T_VISCO/GAMA
	          K3= 1.0D0
	          K4=-1.0D0 
            CASE(4)!HOOKE'S MODEL
              K1= 1.0D0
	          K2= 1.0D0
	          K3= 0.0D0
	          K4= 0.0D0  
            ENDSELECT
        ENDIF
        
        VET_PROP(1)=EMP(NDI_MEF,1)/(2*(1+EMP(NDI_MEF,2)))                       !Mu
	    VET_PROP(2)=EMP(NDI_MEF,2)                                              !Nu
	    VET_PROP(3)=((1.D0)/(16.D0*PI*VET_PROP(1)*(1.D0-VET_PROP(2))))          !C1
	    VET_PROP(4)=((-1.D0)/(8.D0*PI*(1.D0-VET_PROP(2))))                      !C2
	    VET_PROP(5)=VET_PROP(1)/(4*PI*(1.D0-VET_PROP(2)))                       !C3
        
        VET_COORD(1)=COORDPCOLOC_MEF(I,1)
        VET_COORD(2)=COORDPCOLOC_MEF(I,2)
        VET_COORD(3)=COORDPCOLOC_MEF(I,3)
        
        !CALL INTERNAL_POINTS_MATRIX(HH_AUX,GG_AUX,VET_COORD,VET_PROP,NDI_MEF)
        CALL INTERNAL_POINTS_MATRIX(H_FC(3*I-2:3*I,:),G_FC(3*I-2:3*I,:),VET_COORD,VET_PROP,NDI_MEF)
        
        IF (VISCO_ANALYSIS) THEN
            DO K=1,3         
                HH_FC(3*(I-1)+K,:) = H_FC(3*(I-1)+K,:)*K3
                GG_FC(3*(I-1)+K,:) = G_FC(3*(I-1)+K,:)*K4
                    
                H_FC(3*(I-1)+K,:) = H_FC(3*(I-1)+K,:)*K1
                G_FC(3*(I-1)+K,:) = G_FC(3*(I-1)+K,:)*K2
            ENDDO
        ENDIF
        
    ENDDO ! _______________________ END OF MATRIXES RELATED TO INTERNAL POINTS ___________________________________________________
    !$omp end do
    
    ! ---------------------------------------------------------------------------------------------------------------------------
    !COMPUTING THE MATRIX G_CF. THIS MATRIX IS CALCULATED CONSIDERING THE INTEGRAL EQUATION CHOICE ------------------------------
    !(EITHER SINGULAR OR HYPERSINGULAR) OVER EACH COLOCATION POINT --------------------------------------------------------------
         
    !$OMP DO SCHEDULE(DYNAMIC)	
    DO I=1,N_COLLOCPOINTS

        DO J=1,N_ELEM  ! FINDING NORMAL VECTOR FOR THIS COLLOCPOINT I __________________________________________
            NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
            ALLOCATE(VALORES(NEN,3))          
            DO K=1,NEN
                IF (COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
                    NOF=K
                    NDI_MEF=ND(J)                
                    DO KK=1,NEN
                        VALORES(K,1)=COORD_NODES(NODES_CONNECTIVITY(J,K),1)
		                VALORES(K,2)=COORD_NODES(NODES_CONNECTIVITY(J,K),2)
		                VALORES(K,3)=COORD_NODES(NODES_CONNECTIVITY(J,K),3)
                    ENDDO
                    CALL NORMAL_OUTWARD(NEN,VALORES,QSI(J,K,1),QSI(J,K,2),ETA,VJC,DVJC,JC)              
                    MAT_ETAF(1,1)=ETA(1)
	                MAT_ETAF(1,4)=ETA(2)
	                MAT_ETAF(1,7)=ETA(3)
                    
                    MAT_ETAF(2,2)=ETA(1)
	                MAT_ETAF(2,5)=ETA(2)
	                MAT_ETAF(2,8)=ETA(3)
                    
                    MAT_ETAF(3,3)=ETA(1)
	                MAT_ETAF(3,6)=ETA(2)
	                MAT_ETAF(3,9)=ETA(3)
                ENDIF
            ENDDO
            DEALLOCATE(VALORES)
        ENDDO  ! END OF FINDING NORMAL VECTOR __________________________________________________________________

        ! ________________________ INTEGRATING EACH REINFORCEMENT ELEMENT ______________________________________       

        IF (VISCO_ANALYSIS) THEN    ! VISCOELASTIC PARAMETERS
            NUM_LEI = TIME_MODEL(NDI_MEF)
            GAMA = EMP_VISCO(NDI_MEF,1)
            Eve = EMP_VISCO(NDI_MEF,2)
            Ee   = EMP(NDI_MEF,1)
            SELECT CASE(NUM_LEI)
            CASE(1)!KELVIN'S MODEL
	          K1= 1.0D0 + GAMA/DELTA_T_VISCO
	          K2= 1.0D0
	          K3= GAMA/DELTA_T_VISCO
	          K4= 0.0D0  
            CASE(2)!BOLTZMANN'S MODEL
              K1= 1.0D0 + GAMA/DELTA_T_VISCO
	          K2= GAMA/DELTA_T_VISCO + (Ee + Eve)/Eve
	          K3= GAMA/DELTA_T_VISCO
	          K4= -GAMA/DELTA_T_VISCO    
            CASE(3)!MAXWELL'S MODEL
              K1= 1.0D0
	          K2= 1.0D0 + DELTA_T_VISCO/GAMA
	          K3= 1.0D0
	          K4=-1.0D0 
            CASE(4)!HOOKE'S MODEL
              K1= 1.0D0
	          K2= 1.0D0
	          K3= 0.0D0
	          K4= 0.0D0  
            ENDSELECT
        ENDIF   
        
        VET_COORD(1)=COORD_COLLOCPOINTS(I,1)
        VET_COORD(2)=COORD_COLLOCPOINTS(I,2)
        VET_COORD(3)=COORD_COLLOCPOINTS(I,3)
        
        Nu=EMP(NDI_MEF,2)                          ! POISSON (Nu)
        Mu=EMP(NDI_MEF,1)/(2*(1+EMP(NDI_MEF,2)))   ! G (Mu)
        
        DO J=1,N_ELEMENTOS_MEF
            IF (ND_MEF(J).EQ.NDI_MEF) THEN
                            
                ! STRATING VARIABLES
                ALLOCATE(DG_CF(3,3*(ORDEM_MEF(J)+1)))
                
                CALL FIBER_INT_DECISION(INT_2D,J,VET_COORD)
                !INT_2D=.TRUE.
                                
                IF (EQ_TYPE(I).EQ.'S') THEN                                 ! SINGULAR EQUATION IS APPLIED
                    
                    IF (INT_2D) THEN                                        ! INTEGRING CYLINDER
                        CALL FIBER_SING_S(J,Nu,Mu,P_INTEGRACAO,QSIW,VET_COORD,DG_CF)
                    ELSE                                                    ! INTEGRING LINE
                        CALL FIBER_NONSING_S(J,Nu,Mu,VET_COORD,DG_CF)
                    ENDIF
    
                ELSE                      ! HIPERSINGULAR EQUATION IS APPLIED
                            ! ---------------------------------------------------
                            ! ---------------- UNDER CONSTRUCTION ---------------
                            ! ---------------------------------------------------
                ENDIF
                
                ! CONTRUIBUITION IN THE GLOBAL MATRIX FROM THE ELEMENT J
                DO K=1,ORDEM_MEF(J)+1
                    KK=CONECTI_MEF(J,K)
                    
                    IF (VISCO_ANALYSIS) THEN
                        DO M=1,3
                            IF (DABS(DG_CF(M,3*K-2)).LT.TOLER) DG_CF(M,3*K-2)=0.0D0
                            IF (DABS(DG_CF(M,3*K-1)).LT.TOLER) DG_CF(M,3*K-1)=0.0D0
                            IF (DABS(DG_CF(M,3*K)).LT.TOLER) DG_CF(M,3*K)=0.0D0
                            
                            G_CF(3*(I-1)+M,3*KK-2)=G_CF(3*(I-1)+M,3*KK-2)-DG_CF(M,3*K-2)*K2
                            G_CF(3*(I-1)+M,3*KK-1)=G_CF(3*(I-1)+M,3*KK-1)-DG_CF(M,3*K-1)*K2
                            G_CF(3*(I-1)+M,3*KK)=G_CF(3*(I-1)+M,3*KK)-DG_CF(M,3*K)*K2
                            
                            GG_CF(3*(I-1)+M,3*KK-2)=GG_CF(3*(I-1)+M,3*KK-2)-DG_CF(M,3*K-2)*K4
                            GG_CF(3*(I-1)+M,3*KK-1)=GG_CF(3*(I-1)+M,3*KK-1)-DG_CF(M,3*K-1)*K4
                            GG_CF(3*(I-1)+M,3*KK)=GG_CF(3*(I-1)+M,3*KK)-DG_CF(M,3*K)*K4
                        ENDDO
                    ELSE
                        DO M=1,3
                            IF (DABS(DG_CF(M,3*K-2)).LT.TOLER) DG_CF(M,3*K-2)=0.0D0
                            IF (DABS(DG_CF(M,3*K-1)).LT.TOLER) DG_CF(M,3*K-1)=0.0D0
                            IF (DABS(DG_CF(M,3*K)).LT.TOLER) DG_CF(M,3*K)=0.0D0
                            
                            G_CF(3*(I-1)+M,3*KK-2)=G_CF(3*(I-1)+M,3*KK-2)-DG_CF(M,3*K-2)
                            G_CF(3*(I-1)+M,3*KK-1)=G_CF(3*(I-1)+M,3*KK-1)-DG_CF(M,3*K-1)
                            G_CF(3*(I-1)+M,3*KK)=G_CF(3*(I-1)+M,3*KK)-DG_CF(M,3*K)
                        ENDDO                        
                    ENDIF
                ENDDO
                
                DEALLOCATE(DG_CF)
                
            ENDIF
        ENDDO ! END OF INTEGRATION EACH REINFORCEMENT ELEMENT __________________________________________________
        
    ENDDO
    !$OMP END DO  
    
!   ----------------------------------------------------------------------------------------------------------------
!   -- DETERMINATION OF G_FF MATRIX -> INTEGRATING REINF POINTS OVER REINF ELEMENTS ---------------------------------
!   ---------------------------------------------------------------------------------------------------------------
    
    !$omp do schedule(dynamic)
    DO I=1,N_NOS_MEF
    
        ! FINDING 2D DOMAIN THAT NODE I BELONGS
        NDI_MEF=0
        DO J=1,N_ELEMENTOS_MEF
	      DO K=1,(ORDEM_MEF(J)+1)
	          IF(CONECTI_MEF(J,K).EQ.I) THEN
	              NDI_MEF=ND_MEF(J)
                ENDIF
            ENDDO
        ENDDO
        
        IF (VISCO_ANALYSIS) THEN    ! VISCOELASTIC PARAMETERS
            NUM_LEI = TIME_MODEL(NDI_MEF)
            GAMA = EMP_VISCO(NDI_MEF,1)
            Eve = EMP_VISCO(NDI_MEF,2)
            Ee   = EMP(NDI_MEF,1)
            SELECT CASE(NUM_LEI)
            CASE(1)!KELVIN'S MODEL
	          K1= 1.0D0 + GAMA/DELTA_T_VISCO
	          K2= 1.0D0
	          K3= GAMA/DELTA_T_VISCO
	          K4= 0.0D0  
            CASE(2)!BOLTZMANN'S MODEL
              K1= 1.0D0 + GAMA/DELTA_T_VISCO
	          K2= GAMA/DELTA_T_VISCO + (Ee + Eve)/Eve
	          K3= GAMA/DELTA_T_VISCO
	          K4= -GAMA/DELTA_T_VISCO    
            CASE(3)!MAXWELL'S MODEL
              K1= 1.0D0
	          K2= 1.0D0 + DELTA_T_VISCO/GAMA
	          K3= 1.0D0
	          K4=-1.0D0 
            CASE(4)!HOOKE'S MODEL
              K1= 1.0D0
	          K2= 1.0D0
	          K3= 0.0D0
	          K4= 0.0D0  
            ENDSELECT
        ENDIF        
        
        Nu=EMP(NDI_MEF,2)                          ! POISSON (Nu)
        Mu=EMP(NDI_MEF,1)/(2*(1+EMP(NDI_MEF,2)))   ! G (Mu)
        VET_COORD(1)=COORDPCOLOC_MEF(I,1)   ! coords of source point
        VET_COORD(2)=COORDPCOLOC_MEF(I,2)
        VET_COORD(3)=COORDPCOLOC_MEF(I,3)   
        
        DO J=1,N_ELEMENTOS_MEF
            IF(ND_MEF(J).EQ.NDI_MEF) THEN
                                
                ALLOCATE(DG_AUX_FF(3,3*(ORDEM_MEF(J)+1)))
                DG_AUX_FF=0.0D0           
                
                CALL FIBER_INT_DECISION(INT_2D,J,VET_COORD)
                !INT_2D=.TRUE.
                
                IF (INT_2D) THEN !______ PERFORMING INTEGRATION OVER CILINDER ______________________________
                
                    CALL FIBER_SING_S(J,Nu,Mu,P_INTEGRACAO,QSIW,VET_COORD,DG_AUX_FF)
                    
                ELSE ! ______ PERFORMING INTEGRATION OVER LINE ______________________________________________
                    
                    CALL FIBER_NONSING_S(J,Nu,Mu,VET_COORD,DG_AUX_FF)  
                    
                ENDIF
                                
                ! SAVING INTO GLOBAL MATRIX
                DO K=1,(ORDEM_MEF(J)+1)
                    
                    KK=CONECTI_MEF(J,K)
                    
                    IF (VISCO_ANALYSIS) THEN
                        DO M=1,3
                            IF (DABS(DG_AUX_FF(M,3*K-2)).LT.TOLER) DG_AUX_FF(M,3*K-2)=0.0D0
                            IF (DABS(DG_AUX_FF(M,3*K-1)).LT.TOLER) DG_AUX_FF(M,3*K-1)=0.0D0
                            IF (DABS(DG_AUX_FF(M,3*K)).LT.TOLER) DG_AUX_FF(M,3*K)=0.0D0
                            
                            G_FF(3*(I-1)+M,3*KK-2)=G_FF(3*(I-1)+M,3*KK-2)-DG_AUX_FF(M,3*K-2)*K2
                            G_FF(3*(I-1)+M,3*KK-1)=G_FF(3*(I-1)+M,3*KK-1)-DG_AUX_FF(M,3*K-1)*K2
                            G_FF(3*(I-1)+M,3*KK)=G_FF(3*(I-1)+M,3*KK)-DG_AUX_FF(M,3*K)*K2
                            
                            GG_FF(3*(I-1)+M,3*KK-2)=GG_FF(3*(I-1)+M,3*KK-2)-DG_AUX_FF(M,3*K-2)*K4
                            GG_FF(3*(I-1)+M,3*KK-1)=GG_FF(3*(I-1)+M,3*KK-1)-DG_AUX_FF(M,3*K-1)*K4
                            GG_FF(3*(I-1)+M,3*KK)=GG_FF(3*(I-1)+M,3*KK)-DG_AUX_FF(M,3*K)*K4
                        ENDDO                        
                    ELSE
                        DO M=1,3
                            IF (DABS(DG_AUX_FF(M,3*K-2)).LT.TOLER) DG_AUX_FF(M,3*K-2)=0.0D0
                            IF (DABS(DG_AUX_FF(M,3*K-1)).LT.TOLER) DG_AUX_FF(M,3*K-1)=0.0D0
                            IF (DABS(DG_AUX_FF(M,3*K)).LT.TOLER) DG_AUX_FF(M,3*K)=0.0D0
                            
                            G_FF(3*(I-1)+M,3*KK-2)=G_FF(3*(I-1)+M,3*KK-2)-DG_AUX_FF(M,3*K-2)
                            G_FF(3*(I-1)+M,3*KK-1)=G_FF(3*(I-1)+M,3*KK-1)-DG_AUX_FF(M,3*K-1)
                            G_FF(3*(I-1)+M,3*KK)=G_FF(3*(I-1)+M,3*KK)-DG_AUX_FF(M,3*K)
                        ENDDO                     
                    ENDIF
 
                ENDDO
            
                DEALLOCATE(DG_AUX_FF)
                                
            ENDIF
        ENDDO  ! END OF REINFORC ELEMENTS LOOP
            
    ENDDO ! END OF REINFORC MESH NODES LOOP
    !$OMP END DO
    !$OMP END PARALLEL     
        
    CLOSE(1)
    
1   IF(EXISTANCE_REINFORCEMENTS.EQ.'N') THEN
	    CLOSE(1)
    ENDIF
    
END SUBROUTINE