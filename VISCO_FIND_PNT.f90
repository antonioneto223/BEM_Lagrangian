    
SUBROUTINE FIND_BOUND_POINT_VISCO

    !
    ! ENCONTRANDO PONTO NA MALHA DE ELEMENTOS DE CONTORNO
    !

    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE VISCO
    USE REINFORCEMENTS
    
    IMPLICIT NONE

    INTEGER::I,J,K,ELEM_A
    REAL*8::COORDS_A(3),QSI_A(2),TOLER_P,DIST1(4)
    
    ! INITIATING VARIABLES
    TOLER_P=0.00010D0
    ALLOCATE(COLLOCPOINTS_VISCO(N_NODES_BOUND_VISCO),BOUND_PNT_VISCO_COLLOC(N_NODES_BOUND_VISCO),&
        ELEMS_BOUND_PNT_VISCO(N_NODES_BOUND_VISCO),QSIS_BOUND_PNT_VISCO(N_NODES_BOUND_VISCO,2))
    ELEMS_BOUND_PNT_VISCO=0
    COLLOCPOINTS_VISCO=0
    QSIS_BOUND_PNT_VISCO=0.0D0
        
    DO I=1,N_NODES_BOUND_VISCO
        
        ! TENTANDO ENCONTRAR PTO DE COLOCACAO
        DO J=1,N_COLLOCPOINTS
            DIST1(1)=COORD_COLLOCPOINTS(J,1)-COORDS_BOUND_VISCO(I,1)
			DIST1(2)=COORD_COLLOCPOINTS(J,2)-COORDS_BOUND_VISCO(I,2)
	        DIST1(3)=COORD_COLLOCPOINTS(J,3)-COORDS_BOUND_VISCO(I,3)
			DIST1(4)=DSQRT(DIST1(1)**2.0D0+DIST1(2)**2.0D0+DIST1(3)**2.0D0)
            IF (DIST1(4).LT.TOLER_P) THEN
                COLLOCPOINTS_VISCO(I)=J
            ENDIF
        ENDDO
                
        IF (COLLOCPOINTS_VISCO(I).EQ.0) THEN  ! CASO NAO TENHA PTO DE COLOCACAO, ENCONTRAR EL E QSI
            
            COORDS_A(:)=COORDS_BOUND_VISCO(I,:)
            CALL FIND_ELEMENT_WITH_PNT(COORDS_A,ELEM_A,QSI_A)
            
            BOUND_PNT_VISCO_COLLOC(I)=.FALSE.
            ELEMS_BOUND_PNT_VISCO(I)=ELEM_A
            QSIS_BOUND_PNT_VISCO(I,:)=QSI_A(:)
        
        ELSE 
            
            BOUND_PNT_VISCO_COLLOC(I)=.TRUE.
            
        ENDIF
    
    ENDDO    
    

END SUBROUTINE FIND_BOUND_POINT_VISCO