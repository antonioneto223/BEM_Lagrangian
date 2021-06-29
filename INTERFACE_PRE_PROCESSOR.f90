SUBROUTINE INTERFACE_PRE_PROCESSOR
!	
    USE ISOPARAMETRIC_MESH
	USE CRACKS_DUAL_BEM	
	USE ANALYSIS
	USE REMESHING
    USE SUB_REGIONS_INTERFACES
    USE COHESIVE
!
	IMPLICIT NONE 
!
	INTEGER::STEP,I,J,K,II,JJ,KK,L,NEN,C,N_CRACKED_ELEM,CONT,NI_NODE_ELEM[ALLOCATABLE](:),NI_ELEM[ALLOCATABLE](:)
    INTEGER::INTERFACE_NODE,ELEM,LOCAL_NODE,POS
    INTEGER::ND_INTERFACE(2)
!
!------------------------------------------------------------------------------------------------------------------------------------------------------
!       FINDING ELEMENTS THAT CONTAINS INTERFACE NODES
!------------------------------------------------------------------------------------------------------------------------------------------------------
    ALLOCATE(NI_NODE_ELEM(4*N_INTERFACE_POINTS))
    NI_NODE_ELEM=0
    CONT=0
    DO KK=1,N_INTERFACE_POINTS
        DO I=1,N_COLLOCPOINTS
	        IF(I.EQ.NI(KK,1))THEN
	            DO J=1,N_ELEM
		            NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
			        DO K=1,NEN
				        IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
                            CONT=CONT+1
                            NI_NODE_ELEM(CONT)=J
				        ENDIF
			        ENDDO
                ENDDO
            ENDIF
        ENDDO
    ENDDO
!------------------------------------------------------------------------------------------------------------------------------------------------------
!       ELIMINATING REPEATED ELEMENTS
!------------------------------------------------------------------------------------------------------------------------------------------------------
    N_ELEM_INT=SIZE(NI_NODE_ELEM)
    DO I=1,N_ELEM_INT
        IF(NI_NODE_ELEM(I).GT.0)THEN
            DO J=1,N_ELEM_INT
                IF(J.GT.I)THEN
                    IF(NI_NODE_ELEM(I).EQ.NI_NODE_ELEM(J))THEN
                        NI_NODE_ELEM(J)=0
                    ENDIF
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    DO I=1,4*N_INTERFACE_POINTS
        IF(NI_NODE_ELEM(I).EQ.0)THEN
            N_ELEM_INT=N_ELEM_INT-1
        ENDIF
    ENDDO    
    ALLOCATE(NI_ELEM(N_ELEM_INT),INTERFACE_CONNECTIVITY(N_ELEM_INT,9))
!------------------------------------------------------------------------------------------------------------------------------------------------------
!       STOCKING ELEMENTS
!------------------------------------------------------------------------------------------------------------------------------------------------------
    CONT=0
    DO I=1,4*N_INTERFACE_POINTS
        IF(NI_NODE_ELEM(I).GT.0)THEN
            CONT=CONT+1
            NI_ELEM(CONT)=NI_NODE_ELEM(I)
        ENDIF
    ENDDO
!------------------------------------------------------------------------------------------------------------------------------------------------------
!       ORGANIZING ELEMENTS INTERFACE CONNECTIVITY
!------------------------------------------------------------------------------------------------------------------------------------------------------
    INTERFACE_CONNECTIVITY=0
    DO KK=1,N_ELEM_INT
        DO J=1,N_ELEM
            IF(NI_ELEM(KK).EQ.J)THEN
                NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                DO K=1,NEN
                    DO I=1,N_INTERFACE_POINTS
                        IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.NI(I,1))THEN
                            INTERFACE_CONNECTIVITY(KK,K)=I
                        ENDIF
                    ENDDO
                ENDDO
            ENDIF
        ENDDO
    ENDDO
!------------------------------------------------------------------------------------------------------------------------------------------------------
!       CALCULATING DATA FOR COHESIVE OPERATORS AND STOCKING NUMBER OF COHESIVE LAW ASSOCIATED TO INTERFACE POINTS
!------------------------------------------------------------------------------------------------------------------------------------------------------
    ALLOCATE(ROTATION_INTERFACE_S(3,3,2*N_INTERFACE_POINTS),ROTATION_INTERFACE_H(3,3,2*N_INTERFACE_POINTS),&
        NI_INTERFACE_NUMBER(N_INTERFACE_POINTS))
    NI_INTERFACE_NUMBER=0
    DO I=1,N_INTERFACE_POINTS
        DO II=1,2        
            INTERFACE_NODE=NI(I,II)
            POS=I+N_INTERFACE_POINTS*(II-1)
            ! ONE OF ELEMENTS THAT CONTAIN INTERFACE NODE
            DO J=1,N_ELEM
                NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                DO K=1,NEN
                    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.INTERFACE_NODE)THEN
                        ELEM=J
                        LOCAL_NODE=K
                    ENDIF
                ENDDO
            ENDDO
            ND_INTERFACE(II)=ND(ELEM)
            CALL LOCAL_COODINATE_SYSTEM(LOCAL_NODE,ELEM,ROTATION_INTERFACE_S(:,:,POS),ROTATION_INTERFACE_H(:,:,POS))
            CONTINUE
        ENDDO
        ! FINDING INTERFACE POINT COHESIVE INTERFACE NUMBER
        DO J=1,NUMBER_INTERFACES
            IF(ND_INTERFACE(1).EQ.ND_INT(J,1))THEN
                IF(ND_INTERFACE(2).EQ.ND_INT(J,2))THEN
                    NI_INTERFACE_NUMBER(I)=J
                ENDIF
            ELSEIF(ND_INTERFACE(1).EQ.ND_INT(J,2))THEN
                IF(ND_INTERFACE(2).EQ.ND_INT(J,1))THEN
                    NI_INTERFACE_NUMBER(I)=J
                ENDIF
            ENDIF
        ENDDO
        IF(NI_INTERFACE_NUMBER(I).EQ.0)THEN
            WRITE(*,*)'BREAK'
        ENDIF
    ENDDO  
!
END SUBROUTINE