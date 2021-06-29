    SUBROUTINE COHESIVE_SYSTEM
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
	INTEGER::I,J,K,II,JJ,KK,LL,CONT,NEN,IPIV[ALLOCATABLE](:),INFO,DIM_XBEM,DIM_XBEM_CRACKFRONT
    INTEGER::AUX_NUM_CHS,ND_NODE,N_CHS_POINTS,NUM_CHS,ELEM,LOCAL_NODE
!
	REAL*8::A[ALLOCATABLE](:,:),B[ALLOCATABLE](:,:),F[ALLOCATABLE](:),SOLUTION[ALLOCATABLE](:),&
    A_AUX[ALLOCATABLE](:,:),F_AUX[ALLOCATABLE](:),VALUES_CON_SUPPS_BACKUP[ALLOCATABLE](:),&
    T0,T1
    REAL*8::R_MAT_S(3,3),R_MAT_H(3,3),R_INV_S(3,3),U_GLOBAL(3),T_GLOBAL(3),U_LOCAL(3),T_LOCAL(3)
!
    N_CHS_POINTS=SIZE(CHS_INT,1)
    T0=SECNDS(0.0)
    ! SOLTANDO DESLOCAMENTO RELATIVO TANGENTE
	II=3*N_COLLOCPOINTS-2*(3*N_CHS_INTERFACE_POINTS+N_CRACK_POINTS)
    !PRENDENDO DESLOCAMENTO RELATIVO TANGENTE
	!II=3*N_COLLOCPOINTS-2*(N_CHS_INTERFACE_POINTS+N_CRACK_POINTS)
!
	ALLOCATE(SOLUTION(3*N_COLLOCPOINTS+NUM_CON_SUPPS),F(II))
    ALLOCATE(VALUES_CON_SUPPS_BACKUP(SIZE(VALUES_CON_SUPPS,1)))
!
    DIM_XBEM=3*N_COLLOCPOINTS+NUM_CON_SUPPS
    B=0.D0
	SOLUTION=0.D0
	F=0.D0
	CONT=0
!
	DO I=1,N_COLLOCPOINTS
!
		II=0
		JJ=0
        KK=0
        LL=0
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
		    IF(II.EQ.-1) THEN
		    ELSE
			    IF(JJ.EQ.-1) THEN
                ELSE                          
				    IF(CHS_B_CONDITIONS(3*I-2).EQ.0) THEN
					    CONT=CONT+1
					    F(CONT)=CHS_U((3*I-2))
				    ELSE
					    CONT=CONT+1
					    F(CONT)=CHS_T((3*I-2))
				    ENDIF
!
				    IF(CHS_B_CONDITIONS(3*I-1).EQ.0) THEN
					    CONT=CONT+1
					    F(CONT)=CHS_U(3*I-1)
				    ELSE
					    CONT=CONT+1
					    F(CONT)=CHS_T(3*I-1)
				    ENDIF
!
				    IF(CHS_B_CONDITIONS(3*I).EQ.0) THEN
					    CONT=CONT+1
					    F(CONT)=CHS_U(3*I)
				    ELSE
					    CONT=CONT+1
					    F(CONT)=CHS_T(3*I)
                    ENDIF
			    ENDIF
		    ENDIF
		ELSE
		    DO J=1,N_CRACK_POINTS
			    IF(I.EQ.NC(J,1)) THEN
				    II=-1
				    JJ=NC(J,2)
			    ELSE
				    IF(I.EQ.NC(J,2)) THEN
					    II=NC(J,1)
					    JJ=-1
				    ENDIF
			    ENDIF
			ENDDO        
		    IF(II.EQ.-1) THEN
		        IF(CHS_B_CONDITIONS(3*I-1).EQ.0) THEN
					CONT=CONT+1
					F(CONT)=CHS_U((3*I-1))
				ELSE
					CONT=CONT+1
					F(CONT)=CHS_T((3*I-1))
				ENDIF
!
				IF(CHS_B_CONDITIONS(3*I).EQ.0) THEN
					CONT=CONT+1
					F(CONT)=CHS_U(3*I)
				ELSE
					CONT=CONT+1
					F(CONT)=CHS_T(3*I)
				ENDIF
!		    
		    ELSE
			    IF(JJ.EQ.-1) THEN
!				    			    		    
		            IF(CHS_B_CONDITIONS(3*I-1).EQ.0) THEN
					    CONT=CONT+1
				 	    F(CONT)=CHS_U((3*I-1))
				    ELSE
					    CONT=CONT+1
					    F(CONT)=CHS_T((3*I-1))
				    ENDIF
!
				    IF(CHS_B_CONDITIONS(3*I).EQ.0) THEN
					    CONT=CONT+1
					    F(CONT)=CHS_U(3*I)
				    ELSE
					    CONT=CONT+1
					    F(CONT)=CHS_T(3*I)
				    ENDIF			    
!			    
			    ELSE
				    IF(CHS_B_CONDITIONS(3*I-2).EQ.0) THEN
					    CONT=CONT+1
					    F(CONT)=CHS_U((3*I-2))
				    ELSE
					    CONT=CONT+1
					    F(CONT)=CHS_T((3*I-2))
				    ENDIF
!
				    IF(CHS_B_CONDITIONS(3*I-1).EQ.0) THEN
					    CONT=CONT+1
					    F(CONT)=CHS_U(3*I-1)
				    ELSE
					    CONT=CONT+1
					    F(CONT)=CHS_T(3*I-1)
				    ENDIF
!
				    IF(CHS_B_CONDITIONS(3*I).EQ.0) THEN
					    CONT=CONT+1
					    F(CONT)=CHS_U(3*I)
				    ELSE
					    CONT=CONT+1
					    F(CONT)=CHS_T(3*I)
				    ENDIF
			    ENDIF
		    ENDIF
		ENDIF
    ENDDO
!
	II=3*N_COLLOCPOINTS-2*(3*N_CHS_INTERFACE_POINTS+N_CRACK_POINTS)
	CALL DMURRV(3*N_COLLOCPOINTS,II,CHS_B,3*N_COLLOCPOINTS,II,F,1,3*N_COLLOCPOINTS,SOLUTION)
	DEALLOCATE(F)
	ALLOCATE(F(DIM_XBEM))
    F=0.D0
	F(1:3*N_COLLOCPOINTS)=SOLUTION(1:3*N_COLLOCPOINTS)
	SOLUTION=0.D0
!
    IF (NUM_CON_SUPPS .GT. 0) THEN
        DO I=1,NUM_CON_SUPPS
            F(3*N_COLLOCPOINTS+I)=VALUES_CON_SUPPS(I)
        ENDDO
    ENDIF 
!	
    IF (NUM_CON_LOADS .GT. 0) THEN
        !WRITE(*,*)'XBEM_FORCE'
        ALLOCATE(F_ENRICHED(3*N_COLLOCPOINTS))
        CALL XBEM_FORCE
        F(1:3*N_COLLOCPOINTS)=F(1:3*N_COLLOCPOINTS)+F_ENRICHED
        DEALLOCATE(F_ENRICHED)
    ENDIF
!
    DIM_XBEM=SIZE(F)
    ALLOCATE(IPIV(DIM_XBEM))
!
!   CALL DMURRV(DIM_XBEM,DIM_XBEM,CHS_INV_A,DIM_XBEM,DIM_XBEM,F,1,DIM_XBEM,SOLUTION)
!   CALL DLFSRG(DIM_XBEM,CHS_INV_A,DIM_XBEM,IPVT_INV_A,F,1,SOLUTION)
    CALL DGETRS('N',DIM_XBEM,1,CHS_INV_A,DIM_XBEM,IPVT_INV_A,F,DIM_XBEM,INFO)
    SOLUTION=F
!
    IF(NUM_CON_SUPPS.GT.0)THEN
        XBEM_SUPP_REACTIONS=F(3*N_COLLOCPOINTS+1:DIM_XBEM)
    ENDIF
!
200 FORMAT(4x,i5,4x,E13.6,4x,E13.6,4x,E13.6) 
!
	DO I=1,N_COLLOCPOINTS
		II=0
		JJ=0
        KK=0
        LL=0
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
!           COESIVO ROTACIONADO
       !     DO J=1,N_CHS_POINTS
       !         IF(I.EQ.CHS_INT(J,1)) THEN
				   ! KK=-1
       !             LL=CHS_INT(J,2)
			    !ELSE
				   ! IF(I.EQ.CHS_INT(J,2)) THEN
					  !  LL=-1
       !                 KK=CHS_INT(J,1)
				   ! ENDIF
       !         ENDIF
       !     ENDDO
!           ACABA COESIVO ROTACIONADO
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
                    ! COESIVO ROTACIONADO
!                    IF(KK.EQ.-1)THEN
!                        ! DIR ANSWERS
!                        ! ROTATING ANSWERS
!                        DO J=1,N_ELEM
!		                    NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
!			                DO K=1,NEN
!				                IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
!					                ELEM=J
!					                LOCAL_NODE=K
!                                    EXIT
!				                ENDIF
!			                ENDDO
!                        ENDDO
!!
!                        CALL LOCAL_COODINATE_SYSTEM(LOCAL_NODE,ELEM,R_MAT_S,R_MAT_H)
!!                       INVERTING ROTATION MATRIX
!                        CALL DLINRG(3,R_MAT_S,3,R_INV_S,3) 
!!                       LOCAL DISPLACEMENTS
!                        U_LOCAL(1)=SOLUTION(3*I-2)
!                        U_LOCAL(2)=SOLUTION(3*I-1)
!                        U_LOCAL(3)=SOLUTION(3*I)
!!                       FINDING GLOBAL DISPLACEMENTS
!                        CALL DMURRV(3,3,R_INV_S,3,3,U_LOCAL,1,3,U_GLOBAL)
!!
!                        CHS_U(3*I-2)=U_GLOBAL(1)
!			            CHS_U(3*I-1)=U_GLOBAL(2)
!			            CHS_U(3*I)=U_GLOBAL(3)                  
!!  
!                        T_LOCAL(1)=CHS_T(3*I-2)
!                        T_LOCAL(2)=CHS_T(3*I-1)
!                        T_LOCAL(3)=CHS_T(3*I)
!                            
!!                       FINDING GLOBAL TRACTIONS
!                        CALL DMURRV(3,3,R_INV_S,3,3,T_LOCAL,1,3,T_GLOBAL)
!!
!                        CHS_T(3*I-2)=T_GLOBAL(1)
!                        CHS_T(3*I-1)=T_GLOBAL(2)
!                        CHS_T(3*I)=T_GLOBAL(3) 
!                    ELSE
!                        IF(LL.EQ.-1)THEN
!                            ! LEFT ANSWERS
!                            ! ROTATING ANSWERS
!                            DO J=1,N_ELEM
!		                        NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
!			                    DO K=1,NEN
!				                    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
!					                    ELEM=J
!					                    LOCAL_NODE=K
!                                        EXIT
!				                    ENDIF
!			                    ENDDO
!                            ENDDO
!!
!                            CALL LOCAL_COODINATE_SYSTEM(LOCAL_NODE,ELEM,R_MAT_S,R_MAT_H)
!!                           INVERTING ROTATION MATRIX
!                            CALL DLINRG(3,R_MAT_S,3,R_INV_S,3) 
!!                           LOCAL DISPLACEMENTS
!                            U_LOCAL(1)=-SOLUTION(3*I-2)-SOLUTION(3*KK-2)
!                            !IF(CHS_B_CONDITIONS(3*I-1).EQ.1)THEN
!                                !CTD
!                                U_LOCAL(2)=SOLUTION(3*KK-1)-SOLUTION(3*I-1)
!                                !CSD
!                                U_LOCAL(3)=SOLUTION(3*I)-SOLUTION(3*KK)
!                                T_LOCAL(2)=CHS_T(3*I-1)
!                                T_LOCAL(3)=CHS_T(3*I)
!                            !ELSE
!                                !CTD
!                                !U_LOCAL(2)=SOLUTION(3*KK-1)-CHS_U(3*I-1)
!                                !!CSD
!                                !U_LOCAL(3)=CHS_U(3*I)-SOLUTION(3*KK)   
!                                !T_LOCAL(2)=SOLUTION(3*I-1)
!                                !T_LOCAL(3)=SOLUTION(3*I)
!                            !ENDIF
!!                           FINDING GLOBAL DISPLACEMENTS
!                            CALL DMURRV(3,3,R_INV_S,3,3,U_LOCAL,1,3,U_GLOBAL)
!!
!                            CHS_U(3*I-2)=U_GLOBAL(1)
!			                CHS_U(3*I-1)=U_GLOBAL(2)
!			                CHS_U(3*I)=U_GLOBAL(3)
!!                   
!                            T_LOCAL(1)=CHS_T(3*I-2)
!                            T_LOCAL(2)=CHS_T(3*I-1)
!                            T_LOCAL(3)=CHS_T(3*I)
!                            
!!                           FINDING GLOBAL TRACTIONS
!                            CALL DMURRV(3,3,R_INV_S,3,3,T_LOCAL,1,3,T_GLOBAL)
!!
!                            CHS_T(3*I-2)=T_GLOBAL(1)
!                            CHS_T(3*I-1)=T_GLOBAL(2)
!                            CHS_T(3*I)=T_GLOBAL(3) 
!!
!                        ELSE
                            ! ACABA COESIVO ROTACIONADO
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
            !COESIVO ROTACIONADO
            !    ENDIF
            !ENDIF
            !ACABA COESIVO ROTACIONADO
        ELSE
		    DO J=1,N_CRACK_POINTS
			    IF(I.EQ.NC(J,1)) THEN
				    II=-1
				    JJ=NC(J,2)
			    ELSE
				    IF(I.EQ.NC(J,2)) THEN
					    II=NC(J,1)
					    JJ=-1
				    ENDIF
			    ENDIF
		    ENDDO	
!
		    IF(II.EQ.-1) THEN
			    U(3*I-2)=SOLUTION(3*I-2)
			    U(3*JJ-2)=-SOLUTION(3*I-2)
!			    		    
			    IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					T(3*I-1)=SOLUTION(3*I-1)
				ELSE
					U(3*I-1)=SOLUTION(3*I-1)
				ENDIF
!
				IF (B_CONDITIONS(3*I).EQ.0) THEN
					T(3*I)=SOLUTION(3*I)
				ELSE
					U(3*I)=SOLUTION(3*I)
				ENDIF
!				    		    
		    ELSE
			    IF(JJ.EQ.-1) THEN
			        T(3*I-2)=SOLUTION(3*I-2)
				    T(3*II-2)=SOLUTION(3*I-2)
!				    			    
			        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					    T(3*I-1)=SOLUTION(3*I-1)
				    ELSE
					    U(3*I-1)=SOLUTION(3*I-1)
				    ENDIF
!
				    IF (B_CONDITIONS(3*I).EQ.0) THEN
					    T(3*I)=SOLUTION(3*I)
				    ELSE
					    U(3*I)=SOLUTION(3*I)
				    ENDIF			    
!			    
			    ELSE
				    IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					    T(3*I-2)=SOLUTION(3*I-2)
				    ELSE
					    U(3*I-2)=SOLUTION(3*I-2)
				    ENDIF
!
				    IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					    T(3*I-1)=SOLUTION(3*I-1)
				    ELSE
					    U(3*I-1)=SOLUTION(3*I-1)
				    ENDIF
!
				    IF (B_CONDITIONS(3*I).EQ.0) THEN
					    T(3*I)=SOLUTION(3*I)
				    ELSE
					    U(3*I)=SOLUTION(3*I)
				    ENDIF
			    ENDIF
		    ENDIF			
		ENDIF
	ENDDO  	     
!
    T1=SECNDS(T0)
    WRITE(*,*)"         TIME AT COHESIVE_SYSTEM:",T1
!    
	END SUBROUTINE