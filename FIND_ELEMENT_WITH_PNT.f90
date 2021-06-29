    SUBROUTINE FIND_ELEMENT_WITH_PNT(COORDS_PNT,ELEM_NODE,QSI_PNT)
!   COORDS_PNT: COORDINATE POINTS
!   ELEM_PNT: ELEMENT THAT CONTAINS POINT
!   QSI_PNT(2): ADIMENSIONAL COORDINATES OF POINT IN ELEMENT
!
    USE ISOPARAMETRIC_MESH
    USE XBEM_SUPPORT_VARIABLES
    USE REMESHING
	USE PROPAGATION
!
    IMPLICIT NONE
!    
    INTEGER::END_FIND,I,J,NEN,ELEM_NODE,K!,M,II,JJ,KK,RANK!,ELEMS,NODES,,NEN,IPIV(3),INFO,I,TOL1,
!
    REAL*8::X_BAR(3),VALUES[ALLOCATABLE](:,:),QSI_PNT(2),DIST[ALLOCATABLE](:),&
    MED_DIST,DX,DY,DZ,DIST_BAR_PNT,COORDS_PNT(3),TOL1,TOL2,COEFFICIENTS[ALLOCATABLE](:,:),&
    PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),X_TEST,Y_TEST,&
    Z_TEST,DIST_TEST
    !,KN(3,2),DQSI(2),X_TENT(3),&
    !COORDS_PNT(3),DPHI(NEN,2),DIF_X(3),COEFFICIENTS(NEN,NEN),PHI(NEN),&
    !D2PHI(NEN,2,2),AUXDBL,TOL2,NORM,RCOND,WORK(10),AUX1(3),AUX2(2)
!    
    END_FIND=0  ! 0 MEANS ELEMENT WAS NOT FIND YET
    I=1
    QSI_PNT=0.D0
    TOL1=1.E-3
    TOL2=1.E-5
!
    DO WHILE (I .LT. N_ELEM+1 .AND. END_FIND .LT. 1)
        IF(CRACKED_ELEM(I).EQ."UNCRACKED")THEN
            NEN=ELEM_TYPE(I)*ORDER_ELEM(I)+(ELEM_TYPE(I)-3)*(ORDER_ELEM(I)-1)*POL_FAMILY(I)
            ALLOCATE(VALUES(NEN,3),DIST(NEN),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2))
            X_BAR=0.D0
            DIST=0.D0
            MED_DIST=0.D0
            DO J=1,NEN
                VALUES(J,:)=COORD_NODES(NODES_CONNECTIVITY(I,J),:)
                DO K=1,3
                    X_BAR(K)=X_BAR(K)+VALUES(J,K)
                ENDDO
            ENDDO
    !
            X_BAR = X_BAR/NEN
    !
            DO J=1,NEN
                DX=X_BAR(1)-VALUES(J,1)
                DY=X_BAR(2)-VALUES(J,2)
                DZ=X_BAR(3)-VALUES(J,3)
                DIST(J)=DSQRT(DX*DX+DY*DY+DZ*DZ)
                MED_DIST=MED_DIST+DIST(J)
            ENDDO
    !
            MED_DIST=2*MED_DIST/NEN
    !       EVALUATING DISTANCE BETWEEN POINT AND ELEMENT 'BARICENTER'
            DX=X_BAR(1)-COORDS_PNT(1)
            DY=X_BAR(2)-COORDS_PNT(2)
            DZ=X_BAR(3)-COORDS_PNT(3)
            DIST_BAR_PNT=DSQRT(DX*DX+DY*DY+DZ*DZ)
            IF(DIST_BAR_PNT.GT.MED_DIST)THEN
                I=I+1
            ELSE
                CALL NEWTON_RAPHSON_ADIMENSIONAL_COORDINATES(NEN,VALUES,COORDS_PNT,QSI_PNT)
                CALL SHAPE_FUNCTIONS(1,QSI_PNT(1),QSI_PNT(2),NEN,VALUES,COEFFICIENTS,X_TEST,Y_TEST,Z_TEST,PHI,DPHI,D2PHI)
                DX=X_TEST-COORDS_PNT(1)
                DY=Y_TEST-COORDS_PNT(2)
                DZ=Z_TEST-COORDS_PNT(3)      
                DIST_TEST=DSQRT(DX*DX+DY*DY+DZ*DZ)
                IF(DIST_TEST.GT.TOL2)THEN
                    I=I+1
                ELSE
                    IF (NEN.EQ.3.OR.NEN.EQ.6)THEN
                        IF (QSI_PNT(1).GT.(1.D0-TOL1).OR.QSI_PNT(1).LT.(0.D0+TOL1).OR.QSI_PNT(2).GT.(1.D0-TOL1).OR.QSI_PNT(2).LT.(0.D0+TOL1).OR.(QSI_PNT(1)+QSI_PNT(2).GT.(1.D0-TOL1)))THEN
                            I=I+1
                        ELSE
                            END_FIND=1
                            ELEM_NODE=I
                        ENDIF
                    ELSEIF(NEN.EQ.4.OR.NEN.EQ.8.OR.NEN.EQ.9)THEN
                        IF (QSI_PNT(1).GT.(1.D0-TOL1).OR.QSI_PNT(1).LT.(-1.D0+TOL1).OR.QSI_PNT(2).GT.(1.D0-TOL1).OR.QSI_PNT(2).LT.(-1.D0+TOL1)) THEN
                            I=I+1
                        ELSE
                            END_FIND=1
                            ELEM_NODE=I
                        ENDIF
                    ENDIF
                ENDIF
            ENDIF
            DEALLOCATE(VALUES,DIST,COEFFICIENTS,PHI,DPHI,D2PHI)
        ELSE
            I=I+1
        ENDIF 
    ENDDO
    IF(END_FIND.EQ.0)THEN
        WRITE(*,*)"POINT DOES NOT BELONG TO ANY ELEMENT! PLEASE CHECK!"
        WRITE(*,*)"X            Y           Z"
        WRITE(*,*)COORDS_PNT(1),COORDS_PNT(2),COORDS_PNT(3)
        READ(*,*)
    ENDIF
    END SUBROUTINE FIND_ELEMENT_WITH_PNT