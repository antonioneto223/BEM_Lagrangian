    SUBROUTINE NEWTON_RAPHSON_ADIMENSIONAL_COORDINATES(NEN,VALUES,COORDS_PNT,QSI_PNT)
!   IN THIS SUBROUTINE, THE USER INFORMS NEN (NUMBER OF ELEMENT NODES), COORDINATES X, Y AND Z OF EACH NODE
!   IN THE VALUES MATRIX, THE POINT COORDINATES OF DESIRED ADIMENSIONAL COORDINATES, STORED IN QSI_PNT
    USE ISOPARAMETRIC_MESH
    USE XBEM_SUPPORT_VARIABLES
!    
    IMPLICIT NONE
!    
    INTEGER::NEN,IPIV(3),INFO,I,TOL1,J,K,M,II,JJ,KK,RANK!,ELEMS,NODES
!
    REAL*8::VALUES(NEN,3),QSI_PNT(2),KN(3,2),DQSI(2),X_TENT(3),&
    COORDS_PNT(3),DPHI(NEN,2),DIF_X(3),COEFFICIENTS(NEN,NEN),PHI(NEN),&
    D2PHI(NEN,2,2),AUXDBL,TOL2,NORM,RCOND,WORK(10),AUX1(3),AUX2(2)
!
    QSI_PNT=0.D0
    TOL1=100
    TOL2=1E-5
    I=1
    NORM=1E5
    CALL SHAPE_FUNCTIONS(1,QSI_PNT(1),QSI_PNT(2),NEN,VALUES,COEFFICIENTS,X_TENT(1),X_TENT(2),X_TENT(3),PHI,DPHI,D2PHI)
    DIF_X=COORDS_PNT-X_TENT
!
    DO WHILE (I .LT. TOL1 .AND. NORM .GT. TOL2)
        CALL SHAPE_FUNCTIONS(3,QSI_PNT(1),QSI_PNT(2),NEN,VALUES,COEFFICIENTS,AUXDBL,AUXDBL,AUXDBL,PHI,DPHI,D2PHI)
        KN=0.D0
        DO II=1,3
            DO JJ=1,2
                DO KK=1,NEN
                    KN(II,JJ)=KN(II,JJ)+VALUES(KK,II)*DPHI(KK,JJ)
                ENDDO
            ENDDO
        ENDDO
        AUX1=DIF_X
        CALL DGELSS(3,2,1,KN,3,AUX1,3,AUX2,RCOND,RANK,WORK,10,INFO)
        DQSI(1:2)=AUX1(1:2)
        QSI_PNT=QSI_PNT+DQSI
        I=I+1
        CALL SHAPE_FUNCTIONS(1,QSI_PNT(1),QSI_PNT(2),NEN,VALUES,COEFFICIENTS,X_TENT(1),X_TENT(2),X_TENT(3),PHI,DPHI,D2PHI)
        DIF_X=COORDS_PNT-X_TENT
        NORM = NORM2(DIF_X)
    ENDDO
!
    I=2
!
    END SUBROUTINE NEWTON_RAPHSON_ADIMENSIONAL_COORDINATES