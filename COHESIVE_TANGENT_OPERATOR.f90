	SUBROUTINE COHESIVE_TANGENT_OPERATOR
!
	USE ISOPARAMETRIC_MESH
	USE CRACKS_DUAL_BEM
    USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE PROPAGATION
    USE COHESIVE
    USE XBEM_FORCE_VARIABLES
    USE XBEM_SUPPORT_VARIABLES
    USE XBEM_CRACKFRONT_VARIABLES
!
	IMPLICIT NONE 
!
	INTEGER::I,J,K,II,JJ,KK,LL,IK,NUM_CHS,ELEM,NEN,IPIV[ALLOCATABLE](:),INFO,DIM_XBEM,&
    N_CHS_POINTS,NUM_CHS_CHS,ND_NODE,LOCAL_NODE,AUX_NUM_CHS,CONT,INTERFACE_NODE,NI_INTERFACE
!
	REAL*8::A[ALLOCATABLE](:,:),B[ALLOCATABLE](:,:),F[ALLOCATABLE](:),SOLUTION[ALLOCATABLE](:),&
    A_AUX[ALLOCATABLE](:,:),F_AUX[ALLOCATABLE](:),VALUES_CON_SUPPS_BACKUP[ALLOCATABLE](:),DP_DCOD,P_N
    REAL*8::VALUES[ALLOCATABLE](:,:),ETA(3),VJC(3),DVJC(3,2),JC
    REAL*8::T0,T1
    REAL*8::R_MAT_S(3,3),R_MAT_H(3,3),R_INV_S(3,3),U_GLOBAL(3),T_GLOBAL(3),U_LOCAL(3),T_LOCAL(3),R_OPENINGS(3)
    REAL*8::CHS_T_LOCAL(3),ROTATION_MATRIX_S(3,3),ROTATION_MATRIX_H(3,3) 
    REAL*8::CHS_T_BACKUP[ALLOCATABLE](:),CHS_U_BACKUP[ALLOCATABLE](:)
!
    T0=SECNDS(0.0)
	II=3*N_COLLOCPOINTS-2*(3*N_CHS_INTERFACE_POINTS+N_CRACK_POINTS)
    N_CHS_POINTS=SIZE(CHS_INT,1)
    DIM_XBEM=3*N_COLLOCPOINTS+NUM_CON_SUPPS
!   alocando B e F de tamanho cheio, sabendo que o produto deles vai dar só a influência da carga desequilibrada
!   no vetor de forças equivalentes (meio gambiarrado)    
	ALLOCATE(A(DIM_XBEM,DIM_XBEM),B(3*N_COLLOCPOINTS,II),F(II),SOLUTION(DIM_XBEM))
    ALLOCATE(CHS_T_BACKUP(SIZE(CHS_T,1)),CHS_U_BACKUP(SIZE(CHS_U,1)))
!
    CHS_T_BACKUP=CHS_T
    CHS_U_BACKUP=CHS_U
	A=0.D0
	B=0.D0
	SOLUTION=0.D0
	F=0.D0
	NUM_CHS=0
    NUM_CHS_CHS=0
    CONT=0
!   
    IF(N_CHS_POINTS.GT.0)THEN
        CALL HG_ROTATION  
    ENDIF
!
	DO I=1,N_COLLOCPOINTS
!
		DO K=1,N_ELEM
		    NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			DO KK=1,NEN
				IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.I) THEN
					JJ=K
                    EXIT
				ENDIF
			ENDDO
	    ENDDO
	    ND_NODE=ND(JJ)
!
        II=0
		JJ=0
        KK=0
        LL=0
!
		DO J=1,N_INTERFACE_POINTS
			IF(I.EQ.NI(J,1)) THEN
				II=-1
				JJ=NI(J,2)
                NI_INTERFACE=NI_INTERFACE_NUMBER(J)
            ELSE
                IF(I.EQ.NI(J,2)) THEN
					II=NI(J,1)
					JJ=-1
                    NI_INTERFACE=NI_INTERFACE_NUMBER(J)
				ENDIF
            ENDIF
        ENDDO 
        DO J=1,N_CHS_POINTS
            IF(I.EQ.CHS_INT(J,1)) THEN
				KK=-1
                LL=CHS_INT(J,2)      
                NUM_CHS=J
			ELSE
				IF(I.EQ.CHS_INT(J,2)) THEN
					LL=-1
                    KK=CHS_INT(J,1)
                    NUM_CHS=J
				ENDIF
			ENDIF
        ENDDO
		IF(II.EQ.-1) THEN
			A(1:3*N_COLLOCPOINTS,(3*I-2))=A(1:3*N_COLLOCPOINTS,(3*I-2))+H(:,(3*I-2))
			A(1:3*N_COLLOCPOINTS,3*I-1)=A(1:3*N_COLLOCPOINTS,3*I-1)+H(:,3*I-1)
			A(1:3*N_COLLOCPOINTS,3*I)=A(1:3*N_COLLOCPOINTS,3*I)+H(:,3*I)
!
			A(1:3*N_COLLOCPOINTS,(3*I-2))=A(1:3*N_COLLOCPOINTS,(3*I-2))+H(:,(3*JJ-2))
			A(1:3*N_COLLOCPOINTS,3*I-1)=A(1:3*N_COLLOCPOINTS,3*I-1)+H(:,3*JJ-1)
			A(1:3*N_COLLOCPOINTS,3*I)=A(1:3*N_COLLOCPOINTS,3*I)+H(:,3*JJ)
		ELSE
			IF(JJ.EQ.-1) THEN
				A(1:3*N_COLLOCPOINTS,(3*I-2))=A(1:3*N_COLLOCPOINTS,(3*I-2))+G(:,(3*I-2))
				A(1:3*N_COLLOCPOINTS,3*I-1)=A(1:3*N_COLLOCPOINTS,3*I-1)+G(:,3*I-1)
				A(1:3*N_COLLOCPOINTS,3*I)=A(1:3*N_COLLOCPOINTS,3*I)+G(:,3*I)
!
				A(1:3*N_COLLOCPOINTS,(3*I-2))=A(1:3*N_COLLOCPOINTS,(3*I-2))-G(:,(3*II-2))
				A(1:3*N_COLLOCPOINTS,3*I-1)=A(1:3*N_COLLOCPOINTS,3*I-1)-G(:,3*II-1)
				A(1:3*N_COLLOCPOINTS,3*I)=A(1:3*N_COLLOCPOINTS,3*I)-G(:,3*II)
            ELSE
                IF(KK.EQ.-1)THEN
                    !UF_N_DIR AND UF_T1/T2_DIR
                    AUX_NUM_CHS=NUM_CHS+N_CHS_POINTS
                    A(1:3*N_COLLOCPOINTS,3*I-2)=H_R(:,3*NUM_CHS-2)-H_R(:,(3*AUX_NUM_CHS-2))
			        A(1:3*N_COLLOCPOINTS,3*I-1)=H_R(:,3*NUM_CHS-1)+H_R(:,3*AUX_NUM_CHS-1)
			        A(1:3*N_COLLOCPOINTS,3*I)=H_R(:,3*NUM_CHS)-H_R(:,3*AUX_NUM_CHS)                 
!
                    CONT=CONT+1
                    B(:,CONT)=G_R(:,(3*NUM_CHS-2))!+G_R(:,(3*AUX_NUM_CHS-2))
                    F(CONT)=CHS_T(3*I-2)
                    !! APPLYING TANGENTIAL FORCES
                    CONT=CONT+1
                    B(:,CONT)=G_R(:,(3*NUM_CHS-1))
                    F(CONT)=CHS_T(3*I-1)
                    CONT=CONT+1
                    B(:,CONT)=G_R(:,(3*NUM_CHS))
                    F(CONT)=CHS_T(3*I)
                ELSE
                    IF(LL.EQ.-1)THEN
                        AUX_NUM_CHS=NUM_CHS+N_CHS_POINTS
                        ! COD
                        CALL COHESIVE_LAW_DERIVATIVE(NI_INTERFACE,CHS_COD(I),DP_DCOD,TANGENT_BILINEAR(NUM_CHS)) 
                        A(1:3*N_COLLOCPOINTS,3*I-2)=-H_R(:,(3*AUX_NUM_CHS-2))-(G_R(:,(3*NUM_CHS-2))+G_R(:,(3*AUX_NUM_CHS-2)))*DP_DCOD
!
                        ! CSD
 			            A(1:3*N_COLLOCPOINTS,3*I-1)=-H_R(:,3*AUX_NUM_CHS-1)
!
                        ! CTD
			            A(1:3*N_COLLOCPOINTS,3*I)=H_R(:,3*AUX_NUM_CHS)
!                        
                        ! DELTA_P BEING APPLIED
                        CONT=CONT+1
                        B(:,CONT)=G_R(:,(3*AUX_NUM_CHS-2))
                        F(CONT)=CHS_T(3*I-2)
                        ! APPLYING TANGENTIAL FORCES
                        CONT=CONT+1
                        B(:,CONT)=G_R(:,(3*AUX_NUM_CHS-1))
                        F(CONT)=CHS_T(3*I-1)
                        CONT=CONT+1
                        B(:,CONT)=G_R(:,(3*AUX_NUM_CHS))
                        F(CONT)=CHS_T(3*I)
                    ELSE
                        IF(CHS_B_CONDITIONS(3*I-2).EQ.0) THEN
					        A(1:3*N_COLLOCPOINTS,(3*I-2))=-G(:,(3*I-2))
!
					        CONT=CONT+1
					        B(:,CONT)=-H(:,(3*I-2))
					        F(CONT)=CHS_U((3*I-2))
				        ELSE
					        A(1:3*N_COLLOCPOINTS,(3*I-2))=H(:,(3*I-2))	
!
					        CONT=CONT+1
					        B(:,CONT)=G(:,(3*I-2))
					        F(CONT)=CHS_T((3*I-2))
				        ENDIF
!
				        IF(CHS_B_CONDITIONS(3*I-1).EQ.0) THEN
					        A(1:3*N_COLLOCPOINTS,3*I-1)=-G(:,3*I-1)
!
					        CONT=CONT+1
					        B(:,CONT)=-H(:,3*I-1)
					        F(CONT)=CHS_U(3*I-1)
				        ELSE
					        A(1:3*N_COLLOCPOINTS,3*I-1)=H(:,3*I-1)	
!
					        CONT=CONT+1
					        B(:,CONT)=G(:,3*I-1)
					        F(CONT)=CHS_T(3*I-1)
				        ENDIF
!
				        IF(CHS_B_CONDITIONS(3*I).EQ.0) THEN
					        A(1:3*N_COLLOCPOINTS,3*I)=-G(:,3*I)
!
					        CONT=CONT+1
				    	    B(:,CONT)=-H(:,3*I)
					        F(CONT)=CHS_U(3*I)
				        ELSE
					        A(1:3*N_COLLOCPOINTS,3*I)=H(:,3*I)	
!
					        CONT=CONT+1
					        B(:,CONT)=G(:,3*I)
					        F(CONT)=CHS_T(3*I)
				        ENDIF
                    ENDIF
                ENDIF
			ENDIF
        ENDIF
    ENDDO
!
    IF(N_CHS_POINTS.GT.0)THEN
        DEALLOCATE(H_R,G_R)
    ENDIF
!
	II=3*N_COLLOCPOINTS-2*(3*N_CHS_INTERFACE_POINTS+N_CRACK_POINTS)
	CALL DMURRV(3*N_COLLOCPOINTS,II,B,3*N_COLLOCPOINTS,II,F,1,3*N_COLLOCPOINTS,SOLUTION)
	DEALLOCATE(F)
	ALLOCATE(F(DIM_XBEM))
    F=0.D0
	F(1:3*N_COLLOCPOINTS)=SOLUTION(1:3*N_COLLOCPOINTS)
	SOLUTION=0.D0
!
   !IF (NUM_CON_SUPPS .GT. 0) THEN
   !     DO I=1,NUM_CON_SUPPS
   !         A(1:3*N_COLLOCPOINTS,3*N_COLLOCPOINTS+I)=AS_COL(:,I)
   !         A(3*N_COLLOCPOINTS+I,1:3*N_COLLOCPOINTS)=AS_LIN(I,:)
   !         F(3*N_COLLOCPOINTS+I)=VALUES_CON_SUPPS(I)
   !     ENDDO
   ! ENDIF 
    IF ((NUM_CON_SUPPS+NODES_DIST_SUPPS).GT.0) THEN
        A(1:3*N_COLLOCPOINTS,3*N_COLLOCPOINTS+1:3*N_COLLOCPOINTS+NUM_CON_SUPPS+NODES_DIST_SUPPS)=AS_COL
        A(3*N_COLLOCPOINTS+1:3*N_COLLOCPOINTS+NUM_CON_SUPPS+NODES_DIST_SUPPS,1:3*N_COLLOCPOINTS)=AS_LIN
        F(3*N_COLLOCPOINTS+1:3*N_COLLOCPOINTS+NUM_CON_SUPPS)=VALUES_CON_SUPPS
        F(3*N_COLLOCPOINTS+NUM_CON_SUPPS+1:3*N_COLLOCPOINTS+NUM_CON_SUPPS+NODES_DIST_SUPPS)=VALUES_DIST_SUPPS
    ENDIF
!	
    !OPEN(6,file='MATRIZ A.txt',status='unknown')
    !!WRITE(5,*)(CRACKTIP_COLLOCPOINTS(K),K=1,N_CRACKTIPS)
    !DO I=1,DIM_XBEM
    !    WRITE(6,'(10000F)')(A(I,K),K=1,DIM_XBEM)
    !ENDDO
!
	ALLOCATE(IPIV(DIM_XBEM))
    !OPEN(6,file='MATRIZ A.txt',status='unknown')
    !DO I=1,DIM_XBEM
    !    WRITE(6,'(10000F)')(A(I,K),K=1,DIM_XBEM)
    !ENDDO
    !CLOSE(6)
    !OPEN(6,file='TESTE_F.txt',status='unknown')
    !DO I=1,DIM_XBEM
    !    WRITE(6,'(1000F)')F(I)
    !ENDDO
    !CLOSE(6)
!    
    CALL DGESV(DIM_XBEM,1,A,DIM_XBEM,IPIV,F,DIM_XBEM,INFO)
!
	SOLUTION=F
!   SAVING XBEM SUPP INFORMATION
    IF(NUM_CON_SUPPS.GT.0)THEN
        XBEM_SUPP_REACTIONS=F(3*N_COLLOCPOINTS+1:DIM_XBEM)
    ENDIF
!
    DO I=1,N_COLLOCPOINTS
		II=0
		JJ=0
        KK=0
        LL=0
!
        DO K=1,N_ELEM
		    NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			DO KK=1,NEN
				IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.I) THEN
					JJ=K
                    LOCAL_NODE=KK
				ENDIF
			ENDDO
	    ENDDO
	    ND_NODE=ND(JJ)
        NEN=ELEM_TYPE(JJ)*ORDER_ELEM(JJ)+(ELEM_TYPE(JJ)-3)*(ORDER_ELEM(JJ)-1)*POL_FAMILY(JJ)
        ALLOCATE(VALUES(NEN,3))
        DO J=1,NEN
		    VALUES(J,1)=COORD_NODES(NODES_CONNECTIVITY(JJ,J),1)
		    VALUES(J,2)=COORD_NODES(NODES_CONNECTIVITY(JJ,J),2)
		    VALUES(J,3)=COORD_NODES(NODES_CONNECTIVITY(JJ,J),3)
        ENDDO
	    CALL NORMAL_OUTWARD(NEN,VALUES,QSI(JJ,LOCAL_NODE,1),QSI(JJ,LOCAL_NODE,2),ETA,VJC,DVJC,JC)
        DEALLOCATE(VALUES)
!
		IF(DUAL_BEM(I).EQ.'B')THEN
		    DO J=1,N_INTERFACE_POINTS
			    IF(I.EQ.NI(J,1)) THEN
				    II=-1
				    JJ=NI(J,2)
			    ELSE
				    IF(I.EQ.NI(J,2)) THEN
					    II=NI(J,1)
					    JJ=-1
				    ENDIF
			    ENDIF
            ENDDO	
!
            DO J=1,N_CHS_POINTS
                IF(I.EQ.CHS_INT(J,1)) THEN
				    KK=-1
                    LL=CHS_INT(J,2)
			    ELSE
				    IF(I.EQ.CHS_INT(J,2)) THEN
					    LL=-1
                        KK=CHS_INT(J,1)
				    ENDIF
                ENDIF
            ENDDO
!
		    IF(II.EQ.-1) THEN
			    CHS_U(3*I-2)=SOLUTION(3*I-2)
			    CHS_U(3*I-1)=SOLUTION(3*I-1)
			    CHS_U(3*I)=SOLUTION(3*I)
!
			    CHS_U(3*JJ-2)=SOLUTION(3*I-2)
			    CHS_U(3*JJ-1)=SOLUTION(3*I-1)
			    CHS_U(3*JJ)=SOLUTION(3*I)
		    ELSE
			    IF(JJ.EQ.-1) THEN
				    CHS_T(3*I-2)=-SOLUTION(3*I-2)
				    CHS_T(3*I-1)=-SOLUTION(3*I-1)
			        CHS_T(3*I)=-SOLUTION(3*I)
!
				    CHS_T(3*II-2)=SOLUTION(3*I-2)
				    CHS_T(3*II-1)=SOLUTION(3*I-1)
				    CHS_T(3*II)=SOLUTION(3*I)
                ELSE
                    IF(KK.EQ.-1)THEN
                        ! DIR ANSWERS
                        ! ROTATING ANSWERS
                        DO J=1,N_ELEM
		                    NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
			                DO K=1,NEN
				                IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
					                ELEM=J
					                LOCAL_NODE=K
				                ENDIF
			                ENDDO
                        ENDDO
!
                        CALL LOCAL_COODINATE_SYSTEM(LOCAL_NODE,ELEM,R_MAT_S,R_MAT_H)
!                       INVERTING ROTATION MATRIX
                        CALL DLINRG(3,R_MAT_S,3,R_INV_S,3) 
!                       LOCAL DISPLACEMENTS
                        U_LOCAL(1)=SOLUTION(3*I-2)
                        U_LOCAL(2)=SOLUTION(3*I-1)
                        U_LOCAL(3)=SOLUTION(3*I)
!                       FINDING GLOBAL DISPLACEMENTS
                        CALL DMURRV(3,3,R_INV_S,3,3,U_LOCAL,1,3,U_GLOBAL)
!
                        CHS_U(3*I-2)=U_GLOBAL(1)
			            CHS_U(3*I-1)=U_GLOBAL(2)
			            CHS_U(3*I)=U_GLOBAL(3)                  
!         
                        T_LOCAL(1)=CHS_T(3*I-2)
                        T_LOCAL(2)=CHS_T(3*I-1)
                        T_LOCAL(3)=CHS_T(3*I)
               !         
               !!        FINDING GLOBAL TRACTIONS
                        CALL DMURRV(3,3,R_INV_S,3,3,T_LOCAL,1,3,T_GLOBAL)
               !
                        CHS_T(3*I-2)=T_GLOBAL(1)
                        CHS_T(3*I-1)=T_GLOBAL(2)
                        CHS_T(3*I)=T_GLOBAL(3)
                    ELSE
                        IF(LL.EQ.-1)THEN
                            ! LEFT ANSWERS
                            ! ROTATING ANSWERS
                            DO J=1,N_ELEM
		                        NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
			                    DO K=1,NEN
				                    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
					                    ELEM=J
					                    LOCAL_NODE=K
				                    ENDIF
			                    ENDDO
                            ENDDO
!
                            CALL LOCAL_COODINATE_SYSTEM(LOCAL_NODE,ELEM,R_MAT_S,R_MAT_H)
!                           INVERTING ROTATION MATRIX
                            CALL DLINRG(3,R_MAT_S,3,R_INV_S,3) 
!                           LOCAL DISPLACEMENTS
                            U_LOCAL(1)=-SOLUTION(3*I-2)-SOLUTION(3*KK-2)
                            U_LOCAL(2)=SOLUTION(3*KK-1)-SOLUTION(3*I-1)
                            U_LOCAL(3)=SOLUTION(3*I)-SOLUTION(3*KK)
!                           FINDING GLOBAL DISPLACEMENTS
                            CALL DMURRV(3,3,R_INV_S,3,3,U_LOCAL,1,3,U_GLOBAL)
!
                            CHS_U(3*I-2)=U_GLOBAL(1)
			                CHS_U(3*I-1)=U_GLOBAL(2)
			                CHS_U(3*I)=U_GLOBAL(3)
!                   
                            T_LOCAL(1)=CHS_T(3*I-2)
                            T_LOCAL(2)=CHS_T(3*I-1)
                            T_LOCAL(3)=CHS_T(3*I)
                            
!                           FINDING GLOBAL TRACTIONS
                            CALL DMURRV(3,3,R_INV_S,3,3,T_LOCAL,1,3,T_GLOBAL)
!
                            CHS_T(3*I-2)=T_GLOBAL(1)
                            CHS_T(3*I-1)=T_GLOBAL(2)
                            CHS_T(3*I)=T_GLOBAL(3) 
!                           COD NODES
                            CHS_COD(I)=CHS_COD(I)+SOLUTION(3*I-2)
                            CHS_COD(KK)=CHS_COD(KK)+SOLUTION(3*I-2)
                        ELSE
				            IF (CHS_B_CONDITIONS(3*I-2).EQ.0) THEN
					            CHS_T(3*I-2)=SOLUTION(3*I-2)
				            ELSE
					            CHS_U(3*I-2)=SOLUTION(3*I-2)
				            ENDIF
!
				            IF (CHS_B_CONDITIONS(3*I-1).EQ.0) THEN
					            CHS_T(3*I-1)=SOLUTION(3*I-1)
				            ELSE
					            CHS_U(3*I-1)=SOLUTION(3*I-1)
				            ENDIF
!
				            IF (CHS_B_CONDITIONS(3*I).EQ.0) THEN
					            CHS_T(3*I)=SOLUTION(3*I)
				            ELSE
					            CHS_U(3*I)=SOLUTION(3*I)
                            ENDIF
                        ENDIF
                    ENDIF
			    ENDIF
            ENDIF
		ENDIF
    ENDDO  	    
!    ! UPDATING TRACTION TAKING INTO ACCOUNT STEP COD's
!    N_CHS_POINTS=SIZE(CHS_INT,1)
!    CHS_T=CHS_T_BACKUP
!    DO II=1,2
!        DO I=1,N_CHS_POINTS
!            INTERFACE_NODE=CHS_INT(I,II)
!            !POS=I+N_INTERFACE_POINTS*(II-1)
!            DO J=1,N_ELEM
!                NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
!                DO K=1,NEN
!                    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.INTERFACE_NODE)THEN
!                        ELEM=J
!                        LOCAL_NODE=K
!                    ENDIF
!                ENDDO
!            ENDDO
!!                    
!            CALL COHESIVE_LAW(1,CHS_COD(INTERFACE_NODE),1.D100,CHS_T_LOCAL(1))
!!                    
!            ! LOCAL ROTATION MATRIX
!            CALL LOCAL_COODINATE_SYSTEM(LOCAL_NODE,ELEM,ROTATION_MATRIX_S,ROTATION_MATRIX_H)
!!                    
!            CHS_T_LOCAL(2)=0.D0
!            CHS_T_LOCAL(3)=0.D0
!            CALL DGESV(3,1,ROTATION_MATRIX_S,3,IPIV,CHS_T_LOCAL,3,INFO)
!            T_GLOBAL=CHS_T_LOCAL
!!
!            CHS_T(3*INTERFACE_NODE-2)=T_GLOBAL(1)
!            CHS_T(3*INTERFACE_NODE-1)=T_GLOBAL(2)
!            CHS_T(3*INTERFACE_NODE)=T_GLOBAL(3)
!        ENDDO
!    ENDDO
!    
    T1=SECNDS(T0)
    WRITE(*,*)"         TIME AT COHESIVE_TANGENT_OPERATOR:",T1
!
	END SUBROUTINE