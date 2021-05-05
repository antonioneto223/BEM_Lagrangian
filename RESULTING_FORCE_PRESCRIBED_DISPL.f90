SUBROUTINE RESULTING_FORCE_PRESCRIBED_DISPL(P_RES)
    
!
!   THIS SUBROUTINE CALCULATES THE RESULTING FORCE ASSOCIATED ALL PATCHES WITH PRESCRIBED DISPLACEMENT
!
    USE ISOPARAMETRIC_MESH
    USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS	 
!
    IMPLICIT NONE       
!
    LOGICAL::NEW_ELEM,ELEMS_ALREADY_CALCULATED[allocatable](:),INTEGRA(3)
    INTEGER::I,J,K,L,II,KK,ELEM,NEN,CONT,DIR,N_GAUSS
    REAL*8,INTENT(OUT)::P_RES(3)
    REAL*8::TOL_0,VJC(3),DVJC(3,2),JC,ETA(3),AUX_P(3)
    REAL*8::COEFFICIENTS[ALLOCATABLE](:,:),AUX_DBL,T_ELEM[ALLOCATABLE](:,:),&
    PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),QSIW_GAUSS[ALLOCATABLE](:,:),&
    X,Y,Z,VALUES[ALLOCATABLE](:,:)
!
    TOL_0=1.D-8
    P_RES=0.D0
    ALLOCATE(ELEMS_ALREADY_CALCULATED(N_ELEM))
    ELEMS_ALREADY_CALCULATED=.FALSE.
    N_GAUSS=4
    ALLOCATE(QSIW_GAUSS(N_GAUSS,2))
    CALL GAUSS_POINTS(N_GAUSS,QSIW_GAUSS)
!
    ALLOCATE(T_ELEM(1,1),COEFFICIENTS(1,1),PHI(1),DPHI(1,2),D2PHI(1,2,2),VALUES(1,3))
    
    !OPEN(14,FILE='Output/RESULTING_F_EQUIVALENT.TXT',STATUS='UNKNOWN')
!
    DO I=1,N_COLLOCPOINTS
        
        INTEGRA=.FALSE.
        DO DIR=1,3
            IF((B_CONDITIONS((3*I-3)+DIR).EQ.0).AND.(DABS(U((3*I-3)+DIR)).GT.TOL_0)) INTEGRA(DIR)=.TRUE.
            !IF((B_CONDITIONS((3*I-3)+DIR).EQ.0)) INTEGRA(DIR)=.TRUE.
        ENDDO
        
        IF((INTEGRA(1)).OR.(INTEGRA(2)).OR.(INTEGRA(3)))THEN
!           FINDING ELEMENT THAT CONTAINS COLLOCATION POINT
!           ELEM = ELEMENT CONTAINING COLLOCATION POINT  
            DO K=1,N_ELEM
		        NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			    DO KK=1,NEN
				    IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.I) THEN
					    ELEM=K
                        IF(.NOT.ELEMS_ALREADY_CALCULATED(ELEM))THEN
                            ELEMS_ALREADY_CALCULATED(ELEM)=.TRUE.
                        
                            IF(SIZE(PHI).NE.NEN)THEN
                                DEALLOCATE(COEFFICIENTS,PHI,DPHI,D2PHI,T_ELEM,VALUES)
                                ALLOCATE(COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2),T_ELEM(NEN,3),VALUES(NEN,3))
                            ENDIF                           
                !
                            CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM,NEN,COEFFICIENTS)
                !
	                        DO J=1,NEN
		                        VALUES(J,1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,J),1)
		                        VALUES(J,2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,J),2)
		                        VALUES(J,3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,J),3)
                            ENDDO
                !
                            DO DIR=1,3
                                DO J=1,NEN
                                    T_ELEM(J,DIR)=T((3*COLLOCPOINTS_CONNECTIVITY(ELEM,J)-3)+DIR)
                                ENDDO
                            ENDDO
            
                            ! INTEGRAL
                            DO J=1,N_GAUSS
                                DO II=1,N_GAUSS
                                    CALL SHAPE_FUNCTIONS(2,QSIW_GAUSS(J,1),QSIW_GAUSS(II,1),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
                                    CALL NORMAL_OUTWARD(NEN,VALUES,QSIW_GAUSS(J,1),QSIW_GAUSS(II,1),ETA,VJC,DVJC,JC)  
                                    AUX_P=0.D0
                                    DO DIR=1,3
                                        DO L=1,NEN    
                                            AUX_P(DIR)=AUX_P(DIR)+T_ELEM(L,DIR)*PHI(L)
                                        ENDDO
                                    ENDDO
                                    P_RES(1)=P_RES(1)+AUX_P(1)*JC*QSIW_GAUSS(J,2)*QSIW_GAUSS(II,2)
                                    P_RES(2)=P_RES(2)+AUX_P(2)*JC*QSIW_GAUSS(J,2)*QSIW_GAUSS(II,2)
                                    P_RES(3)=P_RES(3)+AUX_P(3)*JC*QSIW_GAUSS(J,2)*QSIW_GAUSS(II,2)
                                ENDDO
                            ENDDO                            
                            
                        ENDIF
				    ENDIF
			    ENDDO
            ENDDO
!               
        ENDIF
                
    ENDDO
!   
    !WRITE(14,*)'RESULTING FORCE:'
    !WRITE(14,*)''
    !WRITE(14,*)'FX       FY      FZ'
    !WRITE(14,*)P_RES(1),P_RES(2),P_RES(3)
    
    !CLOSE(14)
!    
    
END SUBROUTINE
    
!    
!SUBROUTINE RESULTING_FORCE_PRESCRIBED_DISPL(DIR,P_RES)
!
!!   THIS SUBROUTINE CALCULATES THE RESULTING FORCE ASSOCIATED ALL PATCHES WITH PRESCRIBED DISPLACEMENT
!
!    USE ISOPARAMETRIC_MESH
!    USE MATERIALS
!	USE CRACKS_DUAL_BEM
!	USE SUB_REGIONS_INTERFACES
!	USE ANALYSIS	 
!!
!    IMPLICIT NONE       
!!
!    LOGICAL::NEW_ELEM
!    INTEGER::I,J,K,L,II,KK,DIR,ELEMS_ALREADY_CALCULATED[ALLOCATABLE](:),ELEM,NEN,CONT
!    REAL*8::TOL_0,P_RES,VJC(3),DVJC(3,2),JC,ETA(3),AUX_P
!    REAL*8::COEFFICIENTS[ALLOCATABLE](:,:),AUX_DBL,T_ELEM[ALLOCATABLE](:),&
!    PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),QSIW_GAUSS[ALLOCATABLE](:,:),&
!    X,Y,Z,VALUES[ALLOCATABLE](:,:)
!!
!    CONT=0
!    TOL_0=1.D-8
!    P_RES=0.D0
!    ALLOCATE(ELEMS_ALREADY_CALCULATED(N_ELEM))
!    ELEMS_ALREADY_CALCULATED=0
!    !N_GAUSSPOINTS=2
!    ALLOCATE(QSIW_GAUSS(N_GAUSSPOINTS,2))
!    CALL GAUSS_POINTS(N_GAUSSPOINTS,QSIW_GAUSS)
!!
!    ALLOCATE(T_ELEM(1),COEFFICIENTS(1,1),PHI(1),DPHI(1,2),D2PHI(1,2,2),VALUES(1,3))
!!
!    DO I=1,N_COLLOCPOINTS
!        NEW_ELEM=.TRUE.
!        !IF((B_CONDITIONS((3*I-3)+DIR).EQ.0).AND.(DABS(U((3*I-3)+DIR)).GT.TOL_0))THEN
!        IF((B_CONDITIONS((3*I-3)+DIR).EQ.0))THEN
!!           FINDING ELEMENT THAT CONTAINS COLLOCATION POINT
!!           ELEM = ELEMENT CONTAINING COLLOCATION POINT  
!            DO K=1,N_ELEM
!		        NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
!			    DO KK=1,NEN
!				    IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.I) THEN
!					    ELEM=K
!                        EXIT
!				    ENDIF
!			    ENDDO
!            ENDDO
!!
!            DO K=1,N_ELEM
!                IF(ELEM.EQ.ELEMS_ALREADY_CALCULATED(K))THEN
!                    NEW_ELEM=.FALSE.
!                    EXIT
!                ENDIF
!            ENDDO
!!
!            IF(NEW_ELEM)THEN
!                CONT=CONT+1
!                ELEMS_ALREADY_CALCULATED(CONT)=ELEM
!!                
!                IF(SIZE(PHI).NE.NEN)THEN
!                    DEALLOCATE(COEFFICIENTS,PHI,DPHI,D2PHI,T_ELEM,VALUES)
!                    ALLOCATE(COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2),T_ELEM(NEN),VALUES(NEN,3))
!                ENDIF                           
!!
!                CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM,NEN,COEFFICIENTS)
!!
!	            DO J=1,NEN
!		            VALUES(J,1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,J),1)
!		            VALUES(J,2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,J),2)
!		            VALUES(J,3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,J),3)
!	            ENDDO
!!
!                DO K=1,NEN
!                    T_ELEM(K)=T((3*COLLOCPOINTS_CONNECTIVITY(ELEM,K)-3)+DIR)
!                ENDDO
!                ! INTEGRAL
!                DO J=1,N_GAUSSPOINTS
!                    DO K=1,N_GAUSSPOINTS
!                        CALL SHAPE_FUNCTIONS(2,QSIW_GAUSS(J,1),QSIW_GAUSS(K,1),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
!                        CALL NORMAL_OUTWARD(NEN,VALUES,QSIW_GAUSS(J,1),QSIW_GAUSS(K,1),ETA,VJC,DVJC,JC)  
!                        AUX_P=0.D0
!                        DO L=1,NEN    
!                            AUX_P=AUX_P+T_ELEM(L)*PHI(L)
!                        ENDDO
!                        P_RES=P_RES+AUX_P*JC*QSIW_GAUSS(J,2)*QSIW_GAUSS(K,2)
!                    ENDDO
!                ENDDO
!            ENDIF
!        ENDIF
!    ENDDO
!!   
!    WRITE(*,*)'RESULTING FORCE:',P_RES
!!    
!END SUBROUTINE