	SUBROUTINE HG_ROTATION
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
    INTEGER::INT_N1,INT_N2,NEN,LOCAL_N1,LOCAL_N2,ELEM_N1,ELEM_N2
    INTEGER::I,J,K,I_ALT
    REAL*8::R_MAT(3,3),R_MAT_H(3,3),R_MAT_INV(3,3)
    REAL*8,ALLOCATABLE::MAT_H(:,:),MAT_G(:,:),MAT_H_R(:,:),MAT_G_R(:,:)
!
    N_CHS_POINTS=SIZE(CHS_INT,1)
!
    ALLOCATE(H_R(3*N_COLLOCPOINTS,3*2*N_CHS_POINTS),G_R(3*N_COLLOCPOINTS,3*2*N_CHS_POINTS))
    ALLOCATE(MAT_H(3*N_COLLOCPOINTS,3),MAT_G(3*N_COLLOCPOINTS,3),MAT_H_R(3*N_COLLOCPOINTS,3),MAT_G_R(3*N_COLLOCPOINTS,3))
!
    H_R=0.D0
    G_R=0.D0
!
    DO I=1,N_CHS_POINTS
        INT_N1=CHS_INT(I,1)
        INT_N2=CHS_INT(I,2)
!       FINDING ELEM THAT CONTAINS INT_N1 AND INT_N2, AND L_N1 AND L_N2
        DO J=1,N_ELEM
            NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
            DO K=1,NEN
                IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.INT_N1)THEN
                    ELEM_N1=J
                    LOCAL_N1=K
                ENDIF
                IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.INT_N2)THEN
                    ELEM_N2=J
                    LOCAL_N2=K
                ENDIF               
            ENDDO
        ENDDO
!
        ! INTERFACE NODE 1
        ! STORING H AND G TO ROTATE
        MAT_H=H(:,3*INT_N1-2:3*INT_N1)
        MAT_G=G(:,3*INT_N1-2:3*INT_N1)
        ! ROTATION MATRIX
        CALL LOCAL_COODINATE_SYSTEM(LOCAL_N1,ELEM_N1,R_MAT,R_MAT_H)
        ! ROTATION MATRIX INVERSE
        CALL DLINRG(3,R_MAT,3,R_MAT_INV,3)
        ! R^-1*MAT_H/G 
        !DMRRRR (3,9,MAT_ETAS,3,9,3,G_LOCAL,9,3,3,DGAUX2,3)
        CALL DMRRRR(3*N_COLLOCPOINTS,3,MAT_H,3*N_COLLOCPOINTS,3,3,R_MAT_INV,3,3*N_COLLOCPOINTS,3,MAT_H_R,3*N_COLLOCPOINTS)
        CALL DMRRRR(3*N_COLLOCPOINTS,3,MAT_G,3*N_COLLOCPOINTS,3,3,R_MAT_INV,3,3*N_COLLOCPOINTS,3,MAT_G_R,3*N_COLLOCPOINTS)
        ! STOCKING IN H/G_R
        H_R(:,3*I-2:3*I)=MAT_H_R(:,:)
        G_R(:,3*I-2:3*I)=MAT_G_R(:,:)
!
        I_ALT=I+N_CHS_POINTS
        ! INTERFACE NODE 2
        ! STORING H AND G TO ROTATE
        MAT_H=H(:,3*INT_N2-2:3*INT_N2)
        MAT_G=G(:,3*INT_N2-2:3*INT_N2)
        ! ROTATION MATRIX
        CALL LOCAL_COODINATE_SYSTEM(LOCAL_N2,ELEM_N2,R_MAT,R_MAT_H)
        ! ROTATION MATRIX INVERSE
        CALL DLINRG(3,R_MAT,3,R_MAT_INV,3)
        ! R^-1*MAT_H/G 
        !DMRRRR (3,9,MAT_ETAS,3,9,3,G_LOCAL,9,3,3,DGAUX2,3)
        CALL DMRRRR(3*N_COLLOCPOINTS,3,MAT_H,3*N_COLLOCPOINTS,3,3,R_MAT_INV,3,3*N_COLLOCPOINTS,3,MAT_H_R,3*N_COLLOCPOINTS)
        CALL DMRRRR(3*N_COLLOCPOINTS,3,MAT_G,3*N_COLLOCPOINTS,3,3,R_MAT_INV,3,3*N_COLLOCPOINTS,3,MAT_G_R,3*N_COLLOCPOINTS)
        ! STOCKING IN H/G_R
        H_R(:,3*I_ALT-2:3*I_ALT)=MAT_H_R(:,:)
        G_R(:,3*I_ALT-2:3*I_ALT)=MAT_G_R(:,:)
    ENDDO
!
    END SUBROUTINE