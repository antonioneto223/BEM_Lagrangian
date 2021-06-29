SUBROUTINE HG_MATRICES_VISCO
    
    USE VISCO
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    
    IMPLICIT NONE
    
    INTEGER::I,J,K,II,JJ,KK,JJJ,NEN,N,KKJ,NUM_LEI
    
    REAL*8::C(3,3),DX,DY,DZ,DR,R_S(3,3),R_H(3,3),H_LOCAL(3),G_LOCAL(3),VALUES(9,3),&
        LENGTH(4),AVERAGE_LENGTH,DISTANCE,DAUX,GAMA,EVE,EE,K1,K2,K3,K4,MP(3)
    
    REAL*8,DIMENSION(:),ALLOCATABLE::SOMA_H
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::DG,DH
    
    ! INICIATING VISCO MATRICES
    ALLOCATE(HB(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS),GB(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS)&
        ,HBB(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS),GBB(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS))
    HB=0.0D0
    GB=0.0D0
    HBB=0.0D0
    GBB=0.0D0
    
    !CALL OMP_SET_NUM_THREADS(24)
!$OMP PARALLEL PRIVATE(I,J,K,II,JJ,KK,JJJ,NEN,N,DG,DH,DX,DY,DZ,DR,VALUES,LENGTH,AVERAGE_LENGTH,DISTANCE,DAUX,K1,K2,K3,K4,MP,EE,NUM_LEI,GAMA,EVE)
!$OMP DO SCHEDULE(DYNAMIC)	
    DO I=1,N_COLLOCPOINTS
        
        ! FINDING DOMAIN OF COLLOC POINT I
        JJ=0
		DO K=1,N_ELEM
		    NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			DO KK=1,NEN
				IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.I) THEN
					JJ=K !! JJ guarda elemento do ponto fonte i
				ENDIF
			ENDDO
	    ENDDO
	    JJJ=ND(JJ) !! JJJ guarda qual dominio está o ponto fonte i
        
        ! GETTING MATERIAL PROPERTIES
        MP(1)=EMP(ND(JJ),1) ! E
        MP(2)=EMP(ND(JJ),2) ! Poisson
        Ee = MP(1)
            
        NUM_LEI = TIME_MODEL(ND(JJ))
        GAMA = EMP_VISCO(ND(JJ),1)
        Eve = EMP_VISCO(ND(JJ),2)
            
        SELECT CASE(NUM_LEI)
        CASE(1) ! KELVIN MODEL
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
        
        ! START INTEGRATION OF ELEMENTS (J)
		DO J=1,N_ELEM
			IF(ND(JJ).EQ.ND(J)) THEN !  integrates of elements of the same domain
				
                II=0
				! NEN guarda o numero de nós do elemento integrado J
                NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
				DO K=1,NEN
					IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
						II=1 
					ENDIF
				ENDDO
!				
				IF(II.EQ.0) THEN
		            DO K=1,NEN
			            DX=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),1)-COORD_COLLOCPOINTS(I,1)
			            DY=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),2)-COORD_COLLOCPOINTS(I,2)
	                    DZ=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),3)-COORD_COLLOCPOINTS(I,3)
			            DR=DSQRT(DX*DX+DY*DY+DZ*DZ)
			            IF(DR.LE.1.E-12) THEN
				            II=1   ! II=1 se integracao singular, senao 0
				        ENDIF
		            ENDDO
	            ENDIF
!------------------------------------------------------------------------------------------------------------------------------------------------------
!               ADAPTATIVE GAUSS POINTS
!------------------------------------------------------------------------------------------------------------------------------------------------------    
			    VALUES=0.D0
!
	            DO K=1,NEN
		            VALUES(K,1)=COORD_NODES(NODES_CONNECTIVITY(J,K),1)
		            VALUES(K,2)=COORD_NODES(NODES_CONNECTIVITY(J,K),2)
		            VALUES(K,3)=COORD_NODES(NODES_CONNECTIVITY(J,K),3)
	            ENDDO
	            SELECT CASE(ELEM_TYPE(J))
	            CASE(3)
	                LENGTH(1)=DSQRT((VALUES(1,1)-VALUES(2,1))**2+(VALUES(1,2)-VALUES(2,2))**2+(VALUES(1,3)-VALUES(2,3))**2)
                    LENGTH(2)=DSQRT((VALUES(1,1)-VALUES(3,1))**2+(VALUES(1,2)-VALUES(3,2))**2+(VALUES(1,3)-VALUES(3,3))**2)
                    LENGTH(3)=DSQRT((VALUES(2,1)-VALUES(3,1))**2+(VALUES(2,2)-VALUES(3,2))**2+(VALUES(2,3)-VALUES(3,3))**2)
                    AVERAGE_LENGTH=(LENGTH(1)+LENGTH(2)+LENGTH(3))/3
                CASE(4)
                	LENGTH(1)=DSQRT((VALUES(1,1)-VALUES(2,1))**2+(VALUES(1,2)-VALUES(2,2))**2+(VALUES(1,3)-VALUES(2,3))**2)
                    LENGTH(2)=DSQRT((VALUES(2,1)-VALUES(3,1))**2+(VALUES(2,2)-VALUES(3,2))**2+(VALUES(2,3)-VALUES(3,3))**2)
                    LENGTH(3)=DSQRT((VALUES(3,1)-VALUES(4,1))**2+(VALUES(3,2)-VALUES(4,2))**2+(VALUES(3,3)-VALUES(4,3))**2)
                    LENGTH(4)=DSQRT((VALUES(4,1)-VALUES(1,1))**2+(VALUES(4,2)-VALUES(1,2))**2+(VALUES(4,3)-VALUES(1,3))**2)
                    AVERAGE_LENGTH=(LENGTH(1)+LENGTH(2)+LENGTH(3)+LENGTH(4))/4
                ENDSELECT    
                DISTANCE=0.D0
                DO K=1,NEN
                    DAUX=DSQRT((COORD_COLLOCPOINTS(I,1)-VALUES(K,1))**2+(COORD_COLLOCPOINTS(I,2)-VALUES(K,2))**2+(COORD_COLLOCPOINTS(I,3)-VALUES(K,3))**2)
                    IF(K.EQ.1)THEN
                        DISTANCE=DAUX
                    ELSE
                        IF(DAUX.LT.DISTANCE)THEN
                            DISTANCE=DAUX
                        ENDIF
                    ENDIF
                ENDDO
                N=DNINT(DSQRT(10/(DISTANCE/AVERAGE_LENGTH)))	
                IF(N.EQ.0)THEN
                    N=1
                ENDIF		    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!               END OF ADAPTATIVE GAUSS POINTS
!------------------------------------------------------------------------------------------------------------------------------------------------------ 			    
                   
                ! --------------- FOI PARA HG_FREETERM_VISCO (adicao free terms) -----------------------------------
       !         IF(ORDEM(JJ).EQ.ORDEM(J)) THEN
				   ! IF(JJ.NE.J) THEN
					  !  AUX3=0.D0
					  !  DO K=1,(ORDEM(JJ)+1)		
						 !   AUX1=COORDPCOLOC(CONECTI(JJ,K),1) - COORDPCOLOC(CONECTI(J,(ORDEM(J)+2-K)),1)
						 !   AUX2=COORDPCOLOC(CONECTI(JJ,K),2) - COORDPCOLOC(CONECTI(J,(ORDEM(J)+2-K)),2)
						 !   AUX3=AUX3+DSQRT(AUX1*AUX1+AUX2*AUX2)
					  !  ENDDO
					  !  IF(AUX3.LE.LTOL) THEN
						 !   III=1
       !
						 !   IF(TIPOEQ(I).EQ.'HP') THEN
							!    DO KK=1,(ORDEM(J)+1)
       !
							!	    AUX1=COORDPCOLOC(CONECTI(J,KK),1) - COORDPCOLOC(I,1)
							!	    AUX2=COORDPCOLOC(CONECTI(J,KK),2) - COORDPCOLOC(I,2)
							!	    AUX3=DSQRT(AUX1*AUX1+AUX2*AUX2)
							!	    IF(AUX3.LE.1.E-12) THEN
							!	       GB((2*I-1),(2*CONECTI(J,KK)-1)) = GB((2*I-1),(2*CONECTI(J,KK)-1)) +0.5D0*K2
							!	       GB((2*I),(2*CONECTI(J,KK)))     = GB((2*I),(2*CONECTI(J,KK)))     +0.5D0*K2
							!	   
						 !              GBB((2*I-1),(2*CONECTI(J,KK)-1))= GBB((2*I-1),(2*CONECTI(J,KK)-1))+0.5D0*K4
							!	       GBB((2*I),(2*CONECTI(J,KK)))    = GBB((2*I),(2*CONECTI(J,KK)))    +0.5D0*K4
							!	    ENDIF
							!    ENDDO
						 !   ENDIF
       !
					  !  ELSE
						 !   III=0
					  !  ENDIF
				   ! ELSE
					  !  III=0
				   ! ENDIF
			    !ELSE
				   ! III=0
			    !ENDIF
                ! --------------- FIM DA PARTE DE ADICAO DE FREE TERMS ---------------------------------------------
                
                ALLOCATE(DG(3,3*NEN),DH(3,3*NEN))
                DG=0.0D0
                DH=0.0D0
                
                IF (DIST.EQ.0.0D0) THEN                     ! pontos de colocacao sobre os elementos
                    
                    IF (II.EQ.1) THEN                       ! fazer integracao singular
                        IF (EQ_TYPE(I).EQ.'S') THEN         ! equacao singular
                            CALL SINGULAR_ELEMENT_S(I,J,JJJ,DG,DH)
                        ELSE                                ! equacao hipersingular
                            CALL SINGULAR_ELEMENT_HP(I,J,JJJ,DG,DH)
                        ENDIF
                    ELSE                                    ! fazer integracao nao singular
                        IF(ELEM_TYPE(J).EQ.3) THEN          ! elemento triangular
						    IF(EQ_TYPE(I).EQ.'S')THEN
						        CALL TRIANGULAR_ELEMENT_S(I,J,JJJ,DG,DH,N)
						    ELSE
						        CALL TRIANGULAR_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
						    ENDIF   
					    ELSE
                            IF(EQ_TYPE(I).EQ.'S') THEN      ! elemento quadrilateral
						        CALL QUADRILATERAL_ELEMENT_S(I,J,JJJ,DG,DH,N)
						    ELSE
						        CALL QUADRILATERAL_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
						    ENDIF   
                        ENDIF
                    ENDIF
                    
                ELSE                                        ! pontos de colocacao fora dos elementos
                    
                    IF(ELEM_TYPE(J).EQ.3) THEN 
						IF(EQ_TYPE(I).EQ.'S')THEN
						    CALL TRIANGULAR_ELEMENT_S(I,J,JJJ,DG,DH,N)
					    ELSE
						    CALL TRIANGULAR_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
						ENDIF
				    ELSE
                        IF(EQ_TYPE(I).EQ.'S')THEN
						    CALL QUADRILATERAL_ELEMENT_S(I,J,JJJ,DG,DH,N)
						ELSE
						    CALL QUADRILATERAL_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
						ENDIF   
                    ENDIF
                    
                ENDIF ! FIM DA INTEGRACAO DO ELEMENTO
                
                DO K=1,NEN                                  ! espalhando valores integrados DG e DH
                    KK=COLLOCPOINTS_CONNECTIVITY(J,K)

                    !-------------------------------------------------------
                    GB((3*I-2),(3*KK-2)) =GB((3*I-2),(3*KK-2)) +DG(1,(3*K-2))*K2
				    GB((3*I-2),(3*KK-1)) =GB((3*I-2),(3*KK-1)) +DG(1,(3*K-1))*K2
                    GB((3*I-2),(3*KK)) =GB((3*I-2),(3*KK)) +DG(1,(3*K))*K2
                    
                    GB((3*I-1),(3*KK-2)) =GB((3*I-1),(3*KK-2)) +DG(2,(3*K-2))*K2
				    GB((3*I-1),(3*KK-1)) =GB((3*I-1),(3*KK-1)) +DG(2,(3*K-1))*K2
                    GB((3*I-1),(3*KK)) =GB((3*I-1),(3*KK)) +DG(2,(3*K))*K2
                    
                    GB((3*I),(3*KK-2)) =GB((3*I),(3*KK-2)) +DG(3,(3*K-2))*K2
				    GB((3*I),(3*KK-1)) =GB((3*I),(3*KK-1)) +DG(3,(3*K-1))*K2
                    GB((3*I),(3*KK)) =GB((3*I),(3*KK)) +DG(3,(3*K))*K2
                    !-------------------------------------------------------
                    HB((3*I-2),(3*KK-2)) =HB((3*I-2),(3*KK-2)) +DH(1,(3*K-2))*K1
				    HB((3*I-2),(3*KK-1)) =HB((3*I-2),(3*KK-1)) +DH(1,(3*K-1))*K1
                    HB((3*I-2),(3*KK)) =HB((3*I-2),(3*KK)) +DH(1,(3*K))*K1
                    
                    HB((3*I-1),(3*KK-2)) =HB((3*I-1),(3*KK-2)) +DH(2,(3*K-2))*K1
				    HB((3*I-1),(3*KK-1)) =HB((3*I-1),(3*KK-1)) +DH(2,(3*K-1))*K1
                    HB((3*I-1),(3*KK)) =HB((3*I-1),(3*KK)) +DH(2,(3*K))*K1
                    
                    HB((3*I),(3*KK-2)) =HB((3*I),(3*KK-2)) +DH(3,(3*K-2))*K1
				    HB((3*I),(3*KK-1)) =HB((3*I),(3*KK-1)) +DH(3,(3*K-1))*K1
                    HB((3*I),(3*KK)) =HB((3*I),(3*KK)) +DH(3,(3*K))*K1
                    !-------------------------------------------------------
				    GBB((3*I-2),(3*KK-2)) =GBB((3*I-2),(3*KK-2)) +DG(1,(3*K-2))*K4
				    GBB((3*I-2),(3*KK-1)) =GBB((3*I-2),(3*KK-1)) +DG(1,(3*K-1))*K4
                    GBB((3*I-2),(3*KK)) =GBB((3*I-2),(3*KK)) +DG(1,(3*K))*K4
                    
                    GBB((3*I-1),(3*KK-2)) =GBB((3*I-1),(3*KK-2)) +DG(2,(3*K-2))*K4
				    GBB((3*I-1),(3*KK-1)) =GBB((3*I-1),(3*KK-1)) +DG(2,(3*K-1))*K4
                    GBB((3*I-1),(3*KK)) =GBB((3*I-1),(3*KK)) +DG(2,(3*K))*K4
                    
                    GBB((3*I),(3*KK-2)) =GBB((3*I),(3*KK-2)) +DG(3,(3*K-2))*K4
				    GBB((3*I),(3*KK-1)) =GBB((3*I),(3*KK-1)) +DG(3,(3*K-1))*K4
                    GBB((3*I),(3*KK)) =GBB((3*I),(3*KK)) +DG(3,(3*K))*K4
                    !-------------------------------------------------------
                    HBB((3*I-2),(3*KK-2)) =HBB((3*I-2),(3*KK-2)) +DH(1,(3*K-2))*K3
				    HBB((3*I-2),(3*KK-1)) =HBB((3*I-2),(3*KK-1)) +DH(1,(3*K-1))*K3
                    HBB((3*I-2),(3*KK)) =HBB((3*I-2),(3*KK)) +DH(1,(3*K))*K3
                    
                    HBB((3*I-1),(3*KK-2)) =HBB((3*I-1),(3*KK-2)) +DH(2,(3*K-2))*K3
				    HBB((3*I-1),(3*KK-1)) =HBB((3*I-1),(3*KK-1)) +DH(2,(3*K-1))*K3
                    HBB((3*I-1),(3*KK)) =HBB((3*I-1),(3*KK)) +DH(2,(3*K))*K3
                    
                    HBB((3*I),(3*KK-2)) =HBB((3*I),(3*KK-2)) +DH(3,(3*K-2))*K3
				    HBB((3*I),(3*KK-1)) =HBB((3*I),(3*KK-1)) +DH(3,(3*K-1))*K3
                    HBB((3*I),(3*KK)) =HBB((3*I),(3*KK)) +DH(3,(3*K))*K3

                ENDDO
                
                DEALLOCATE(DG,DH)
            ENDIF ! end of if clause for domain
        ENDDO ! end of integration element j
    ENDDO ! end of colloc point i   
    !$OMP END DO
    !$OMP END PARALLEL 
    
END SUBROUTINE
    
    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    

   
SUBROUTINE HG_FREETERM_VISCO(FIRST_STEP)

    !   SUBROUTINE THAT INCLUDES THE FREE TERM AT THE DIAGONAL OF H AND G VISCOELASTIC MATRIXES
    !   USED IN THE EDGE_CRACKS_HG_UPDATE SUBROUTINE 
    !
    
    USE ISOPARAMETRIC_MESH
    USE CRACKS_DUAL_BEM
    USE SUB_REGIONS_INTERFACES
    USE ANALYSIS
    USE VISCO
    
    IMPLICIT NONE
    
    INTEGER::I,J,JJ,K,KK,NEN
    
    REAL*8::K1,K2,K3,K4,Ee,Eve,GAMA,C_BB(3,3),C(3,3),DX,DY,DZ,DR
    
    LOGICAL::FIRST_STEP
    
    IF(FIRST_STEP)THEN
        DO I=1,N_COLLOCPOINTS
                
            JJ=0  ! encontrando elemento jj do ponto de colocacao i
		    DO K=1,N_ELEM
		        NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			    DO KK=1,NEN
				    IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.I) THEN
					    JJ=K 
				    ENDIF
			    ENDDO
            ENDDO
                
            GAMA = EMP_VISCO(ND(JJ),1)
            Eve  = EMP_VISCO(ND(JJ),2)
            SELECT CASE(TIME_MODEL(ND(JJ)))
            CASE(1) ! KELVIN MODEL
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
                
            IF(EQ_TYPE(I).EQ.'S') THEN
		        C=0.D0
	            DO J=1,N_COLLOCPOINTS
			        C(1,1)=C(1,1)-HB(3*I-2,3*J-2)
			        C(1,2)=C(1,2)-HB(3*I-2,3*J-1)
				    C(1,3)=C(1,3)-HB(3*I-2,3*J)
				    C(2,1)=C(2,1)-HB(3*I-1,3*J-2)
				    C(2,2)=C(2,2)-HB(3*I-1,3*J-1)
				    C(2,3)=C(2,3)-HB(3*I-1,3*J)
				    C(3,1)=C(3,1)-HB(3*I,3*J-2)
				    C(3,2)=C(3,2)-HB(3*I,3*J-1)
				    C(3,3)=C(3,3)-HB(3*I,3*J)
                ENDDO
    !
                C_BB=0.D0
	            DO J=1,N_COLLOCPOINTS
			        C_BB(1,1)=C_BB(1,1)-HBB(3*I-2,3*J-2)
			        C_BB(1,2)=C_BB(1,2)-HBB(3*I-2,3*J-1)
				    C_BB(1,3)=C_BB(1,3)-HBB(3*I-2,3*J)
				    C_BB(2,1)=C_BB(2,1)-HBB(3*I-1,3*J-2)
				    C_BB(2,2)=C_BB(2,2)-HBB(3*I-1,3*J-1)
				    C_BB(2,3)=C_BB(2,3)-HBB(3*I-1,3*J)
				    C_BB(3,1)=C_BB(3,1)-HBB(3*I,3*J-2)
				    C_BB(3,2)=C_BB(3,2)-HBB(3*I,3*J-1)
				    C_BB(3,3)=C_BB(3,3)-HBB(3*I,3*J)
			    ENDDO
    !
                IF(DUAL_BEM(I).EQ."B")THEN
				    HB(3*I-2,3*I-2)=HB(3*I-2,3*I-2)+C(1,1)
				    HB(3*I-2,3*I-1)=HB(3*I-2,3*I-1)+C(1,2)
				    HB(3*I-2,3*I)=HB(3*I-2,3*I)+C(1,3)
				    HB(3*I-1,3*I-2)=HB(3*I-1,3*I-2)+C(2,1)
				    HB(3*I-1,3*I-1)=HB(3*I-1,3*I-1)+C(2,2)
				    HB(3*I-1,3*I)=HB(3*I-1,3*I)+C(2,3)
				    HB(3*I,3*I-2)=HB(3*I,3*I-2)+C(3,1)
				    HB(3*I,3*I-1)=HB(3*I,3*I-1)+C(3,2)
				    HB(3*I,3*I)=HB(3*I,3*I)+C(3,3)
                        
                    HBB(3*I-2,3*I-2)=HBB(3*I-2,3*I-2)+C_BB(1,1)
				    HBB(3*I-2,3*I-1)=HBB(3*I-2,3*I-1)+C_BB(1,2)
				    HBB(3*I-2,3*I)=HBB(3*I-2,3*I)+C_BB(1,3)
				    HBB(3*I-1,3*I-2)=HBB(3*I-1,3*I-2)+C_BB(2,1)
				    HBB(3*I-1,3*I-1)=HBB(3*I-1,3*I-1)+C_BB(2,2)
				    HBB(3*I-1,3*I)=HBB(3*I-1,3*I)+C_BB(2,3)
				    HBB(3*I,3*I-2)=HBB(3*I,3*I-2)+C_BB(3,1)
				    HBB(3*I,3*I-1)=HBB(3*I,3*I-1)+C_BB(3,2)
				    HBB(3*I,3*I)=HBB(3*I,3*I)+C_BB(3,3)
		        ELSE 
		            DO J=1,N_COLLOCPOINTS
		                IF(DUAL_BEM(J).EQ."H")THEN
				            DX=COORD_COLLOCPOINTS(I,1)-COORD_COLLOCPOINTS(J,1)
				            DY=COORD_COLLOCPOINTS(I,2)-COORD_COLLOCPOINTS(J,2)
				            DZ=COORD_COLLOCPOINTS(I,3)-COORD_COLLOCPOINTS(J,3)
				            DR=DSQRT(DX*DX+DY*DY+DZ*DZ)
				            IF(DR.LE.1.E-12)THEN
				                K=J
				            ENDIF
				        ENDIF
                    ENDDO
                        
				    HB(3*I-2,3*I-2)=HB(3*I-2,3*I-2)+C(1,1)/2
				    HB(3*I-2,3*I-1)=HB(3*I-2,3*I-1)+C(1,2)/2
				    HB(3*I-2,3*I)=HB(3*I-2,3*I)+C(1,3)/2
				    HB(3*I-1,3*I-2)=HB(3*I-1,3*I-2)+C(2,1)/2
				    HB(3*I-1,3*I-1)=HB(3*I-1,3*I-1)+C(2,2)/2
				    HB(3*I-1,3*I)=HB(3*I-1,3*I)+C(2,3)/2
				    HB(3*I,3*I-2)=HB(3*I,3*I-2)+C(3,1)/2
				    HB(3*I,3*I-1)=HB(3*I,3*I-1)+C(3,2)/2
				    HB(3*I,3*I)=HB(3*I,3*I)+C(3,3)/2
    				         
				    HB(3*I-2,3*K-2)=HB(3*I-2,3*K-2)+C(1,1)/2
				    HB(3*I-2,3*K-1)=HB(3*I-2,3*K-1)+C(1,2)/2
				    HB(3*I-2,3*K)=HB(3*I-2,3*K)+C(1,3)/2
				    HB(3*I-1,3*K-2)=HB(3*I-1,3*K-2)+C(2,1)/2
				    HB(3*I-1,3*K-1)=HB(3*I-1,3*K-1)+C(2,2)/2
				    HB(3*I-1,3*K)=HB(3*I-1,3*K)+C(2,3)/2
				    HB(3*I,3*K-2)=HB(3*I,3*K-2)+C(3,1)/2
				    HB(3*I,3*K-1)=HB(3*I,3*K-1)+C(3,2)/2
				    HB(3*I,3*K)=HB(3*I,3*K)+C(3,3)/2
                        
                        				        
                    HBB(3*I-2,3*I-2)=HBB(3*I-2,3*I-2)+C_BB(1,1)/2
				    HBB(3*I-2,3*I-1)=HBB(3*I-2,3*I-1)+C_BB(1,2)/2
				    HBB(3*I-2,3*I)=HBB(3*I-2,3*I)+C_BB(1,3)/2
				    HBB(3*I-1,3*I-2)=HBB(3*I-1,3*I-2)+C_BB(2,1)/2
				    HBB(3*I-1,3*I-1)=HBB(3*I-1,3*I-1)+C_BB(2,2)/2
				    HBB(3*I-1,3*I)=HBB(3*I-1,3*I)+C_BB(2,3)/2
				    HBB(3*I,3*I-2)=HBB(3*I,3*I-2)+C_BB(3,1)/2
				    HBB(3*I,3*I-1)=HBB(3*I,3*I-1)+C_BB(3,2)/2
				    HBB(3*I,3*I)=HBB(3*I,3*I)+C_BB(3,3)/2
    				         
				    HBB(3*I-2,3*K-2)=HBB(3*I-2,3*K-2)+C_BB(1,1)/2
				    HBB(3*I-2,3*K-1)=HBB(3*I-2,3*K-1)+C_BB(1,2)/2
				    HBB(3*I-2,3*K)=HBB(3*I-2,3*K)+C_BB(1,3)/2
				    HBB(3*I-1,3*K-2)=HBB(3*I-1,3*K-2)+C_BB(2,1)/2
				    HBB(3*I-1,3*K-1)=HBB(3*I-1,3*K-1)+C_BB(2,2)/2
				    HBB(3*I-1,3*K)=HBB(3*I-1,3*K)+C_BB(2,3)/2
				    HBB(3*I,3*K-2)=HBB(3*I,3*K-2)+C_BB(3,1)/2
				    HBB(3*I,3*K-1)=HBB(3*I,3*K-1)+C_BB(3,2)/2
				    HBB(3*I,3*K)=HBB(3*I,3*K)+C_BB(3,3)/2
                        
		        ENDIF
	        ELSE
			    IF(DUAL_BEM(I).EQ."B")THEN
				    GB(3*I-2,3*I-2)=GB(3*I-2,3*I-2)-0.5D0*K2
				    GB(3*I-1,3*I-1)=GB(3*I-1,3*I-1)-0.5D0*K2
				    GB(3*I,3*I)=GB(3*I,3*I)-0.5D0*K2
                        
                    GBB(3*I-2,3*I-2)=GBB(3*I-2,3*I-2)-0.5D0*K4
				    GBB(3*I-1,3*I-1)=GBB(3*I-1,3*I-1)-0.5D0*K4
				    GBB(3*I,3*I)=GBB(3*I,3*I)-0.5D0*K4
		        ELSE
		            DO J=1,N_COLLOCPOINTS
				        IF(DUAL_BEM(J).EQ."S")THEN
				            DX=COORD_COLLOCPOINTS(I,1)-COORD_COLLOCPOINTS(J,1)
				            DY=COORD_COLLOCPOINTS(I,2)-COORD_COLLOCPOINTS(J,2)
				            DZ=COORD_COLLOCPOINTS(I,3)-COORD_COLLOCPOINTS(J,3)
				            DR=DSQRT(DX*DX+DY*DY+DZ*DZ)
				            IF(DR.LE.1.E-12)THEN
				                K=J
				            ENDIF
				        ENDIF
				    ENDDO
 				    GB(3*I-2,3*I-2)=GB(3*I-2,3*I-2)-0.50D0*K2
 				    GB(3*I-1,3*I-1)=GB(3*I-1,3*I-1)-0.50D0*K2
 				    GB(3*I,3*I)=GB(3*I,3*I)-0.50D0*K2
                        
                    GBB(3*I-2,3*I-2)=GBB(3*I-2,3*I-2)-0.50D0*K4
 				    GBB(3*I-1,3*I-1)=GBB(3*I-1,3*I-1)-0.50D0*K4
 				    GBB(3*I,3*I)=GBB(3*I,3*I)-0.50D0*K4
    				         
 				    GB(3*I-2,3*K-2)=GB(3*I-2,3*K-2)+0.5D0*K2
 				    GB(3*I-1,3*K-1)=GB(3*I-1,3*K-1)+0.5D0*K2
 				    GB(3*I,3*K)=GB(3*I,3*K)+0.5D0*K2
                        
                    GBB(3*I-2,3*K-2)=GBB(3*I-2,3*K-2)+0.5D0*K4
 				    GBB(3*I-1,3*K-1)=GBB(3*I-1,3*K-1)+0.5D0*K4
 				    GBB(3*I,3*K)=GBB(3*I,3*K)+0.5D0*K4
	            ENDIF
    !				
		        C=0.D0
			    DO J=1,N_COLLOCPOINTS
			        C(1,1)=C(1,1)-HB(3*I-2,3*J-2)
			        C(1,2)=C(1,2)-HB(3*I-2,3*J-1)
				    C(1,3)=C(1,3)-HB(3*I-2,3*J)
				    C(2,1)=C(2,1)-HB(3*I-1,3*J-2)
				    C(2,2)=C(2,2)-HB(3*I-1,3*J-1)
				    C(2,3)=C(2,3)-HB(3*I-1,3*J)
				    C(3,1)=C(3,1)-HB(3*I,3*J-2)
				    C(3,2)=C(3,2)-HB(3*I,3*J-1)
				    C(3,3)=C(3,3)-HB(3*I,3*J)
		        ENDDO
    !
			    HB(3*I-2,3*I-2)=HB(3*I-2,3*I-2)+C(1,1)
			    HB(3*I-2,3*I-1)=HB(3*I-2,3*I-1)+C(1,2)
			    HB(3*I-2,3*I)=HB(3*I-2,3*I)+C(1,3)
			    HB(3*I-1,3*I-2)=HB(3*I-1,3*I-2)+C(2,1)
			    HB(3*I-1,3*I-1)=HB(3*I-1,3*I-1)+C(2,2)
			    HB(3*I-1,3*I)=HB(3*I-1,3*I)+C(2,3)
			    HB(3*I,3*I-2)=HB(3*I,3*I-2)+C(3,1)
			    HB(3*I,3*I-1)=HB(3*I,3*I-1)+C(3,2)
			    HB(3*I,3*I)=HB(3*I,3*I)+C(3,3)
                    
                C=0.D0
			    DO J=1,N_COLLOCPOINTS
			        C(1,1)=C(1,1)-HBB(3*I-2,3*J-2)
			        C(1,2)=C(1,2)-HBB(3*I-2,3*J-1)
				    C(1,3)=C(1,3)-HBB(3*I-2,3*J)
				    C(2,1)=C(2,1)-HBB(3*I-1,3*J-2)
				    C(2,2)=C(2,2)-HBB(3*I-1,3*J-1)
				    C(2,3)=C(2,3)-HBB(3*I-1,3*J)
				    C(3,1)=C(3,1)-HBB(3*I,3*J-2)
				    C(3,2)=C(3,2)-HBB(3*I,3*J-1)
				    C(3,3)=C(3,3)-HBB(3*I,3*J)
		        ENDDO
    !
			    HBB(3*I-2,3*I-2)=HBB(3*I-2,3*I-2)+C(1,1)
			    HBB(3*I-2,3*I-1)=HBB(3*I-2,3*I-1)+C(1,2)
			    HBB(3*I-2,3*I)=HBB(3*I-2,3*I)+C(1,3)
			    HBB(3*I-1,3*I-2)=HBB(3*I-1,3*I-2)+C(2,1)
			    HBB(3*I-1,3*I-1)=HBB(3*I-1,3*I-1)+C(2,2)
			    HBB(3*I-1,3*I)=HBB(3*I-1,3*I)+C(2,3)
			    HBB(3*I,3*I-2)=HBB(3*I,3*I-2)+C(3,1)
			    HBB(3*I,3*I-1)=HBB(3*I,3*I-1)+C(3,2)
			    HBB(3*I,3*I)=HBB(3*I,3*I)+C(3,3)
                    
		    ENDIF
        ENDDO
        
    ELSE  ! ________ NON FIRST STEP _____________________
            
        DO I=1,N_COLLOCPOINTS	    
		    C=0.D0
            DO J=1,N_COLLOCPOINTS
		        C(1,1)=C(1,1)-HB(3*I-2,3*J-2)
		        C(1,2)=C(1,2)-HB(3*I-2,3*J-1)
                C(1,3)=C(1,3)-HB(3*I-2,3*J)
	            C(2,1)=C(2,1)-HB(3*I-1,3*J-2)
		        C(2,2)=C(2,2)-HB(3*I-1,3*J-1)
		        C(2,3)=C(2,3)-HB(3*I-1,3*J)
		        C(3,1)=C(3,1)-HB(3*I,3*J-2)
		        C(3,2)=C(3,2)-HB(3*I,3*J-1)
		        C(3,3)=C(3,3)-HB(3*I,3*J)
	        ENDDO
!        	        
            HB(3*I-2,3*I-2)=HB(3*I-2,3*I-2)+C(1,1)
            HB(3*I-2,3*I-1)=HB(3*I-2,3*I-1)+C(1,2)
            HB(3*I-2,3*I)=HB(3*I-2,3*I)+C(1,3)
            HB(3*I-1,3*I-2)=HB(3*I-1,3*I-2)+C(2,1)
            HB(3*I-1,3*I-1)=HB(3*I-1,3*I-1)+C(2,2)
            HB(3*I-1,3*I)=HB(3*I-1,3*I)+C(2,3)
            HB(3*I,3*I-2)=HB(3*I,3*I-2)+C(3,1)
            HB(3*I,3*I-1)=HB(3*I,3*I-1)+C(3,2)
            HB(3*I,3*I)=HB(3*I,3*I)+C(3,3) 
                
            C=0.D0
            DO J=1,N_COLLOCPOINTS
		        C(1,1)=C(1,1)-HBB(3*I-2,3*J-2)
		        C(1,2)=C(1,2)-HBB(3*I-2,3*J-1)
                C(1,3)=C(1,3)-HBB(3*I-2,3*J)
	            C(2,1)=C(2,1)-HBB(3*I-1,3*J-2)
		        C(2,2)=C(2,2)-HBB(3*I-1,3*J-1)
		        C(2,3)=C(2,3)-HBB(3*I-1,3*J)
		        C(3,1)=C(3,1)-HBB(3*I,3*J-2)
		        C(3,2)=C(3,2)-HBB(3*I,3*J-1)
		        C(3,3)=C(3,3)-HBB(3*I,3*J)
	        ENDDO
!        	        
            HBB(3*I-2,3*I-2)=HBB(3*I-2,3*I-2)+C(1,1)
            HBB(3*I-2,3*I-1)=HBB(3*I-2,3*I-1)+C(1,2)
            HBB(3*I-2,3*I)=HBB(3*I-2,3*I)+C(1,3)
            HBB(3*I-1,3*I-2)=HBB(3*I-1,3*I-2)+C(2,1)
            HBB(3*I-1,3*I-1)=HBB(3*I-1,3*I-1)+C(2,2)
            HBB(3*I-1,3*I)=HBB(3*I-1,3*I)+C(2,3)
            HBB(3*I,3*I-2)=HBB(3*I,3*I-2)+C(3,1)
            HBB(3*I,3*I-1)=HBB(3*I,3*I-1)+C(3,2)
            HBB(3*I,3*I)=HBB(3*I,3*I)+C(3,3) 
        ENDDO          
    ENDIF 
        
END SUBROUTINE
    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    
