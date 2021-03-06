    SUBROUTINE WILLIAMS_DISPL_FUNCTIONS_DERIV(QSI1,QSI2,DIR,COORDS_TIP,ROTATION_MATRIX,VALUES,NEN,ELEM,Mu,Nu,SIGNAL,SIGNAL_QSI_AL_CT,PSI_TIP_DERIV)
!   R_TIP: DISTANCE BETWEEN POINT AND CRACK FRONT
!   THETA_TIPO: ANGLE BETWEEN POINT AND CRACK FRONT COORDINATE SYSTEM
!   IT MUST ENTER THE VALUES_COLLOC (ASSOCIATED WITH COLLOCATION POINTS)
!   THIS SUBROUTINE CALCULATES WILLIAMS FUNCTION DERIVATIVE IN ADIMENSIONAL COORDINATE WITH DIRECTION DIR (DPSI/DQSI_DIR)
!   DIR: DIRECTION IS PERPENDICULAR TO CRACK FRONT
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE XBEM_SUPPORT_VARIABLES
!   
    IMPLICIT NONE
    CHARACTER(40)::FLAG
    INTEGER::I,J,ELEM,NEN,DIR
!
    REAL*8::PSI_TIP_DERIV(3,3),F_THETA(3,3),Mu,Nu,Kos,PI,R_TIP,THETA_TIP,DR_DQSIAL,VEC_R_TIP(3),&
    DR(3),AUX1,X,Y,Z,COORDS_TIP(3),ROTATION_MATRIX(3,3),QSI1,QSI2,SIGNAL,PSI_TIP(3,3),SIGNAL_QSI_AL_CT
!
    REAL*8::COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),&
    D2PHI(NEN,2,2),AUX_DBL,VALUES(NEN,3),AUXD1,AUXD2,AUXD3,DR_DQSIAL_TESTE
!
    PI=DACOS(-1.D0)
!   KOSOLOV CONSTANT FOR PLANE STRAIN
    Kos=3.D0-4.D0*Nu
!   COORDINATES OF FIELD POINT AND SHAPE FUNCTIONS DERIVATIVE
    CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM,NEN,COEFFICIENTS)
    CALL SHAPE_FUNCTIONS(1,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI) 
    CALL SHAPE_FUNCTIONS(2,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI) 
    CALL SHAPE_FUNCTIONS(3,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,AUX_DBL,AUX_DBL,AUX_DBL,PHI,DPHI,D2PHI)
!    
	DR(1)=X-COORDS_TIP(1)
	DR(2)=Y-COORDS_TIP(2)
	DR(3)=Z-COORDS_TIP(3)
    R_TIP=DSQRT(DR(1)*DR(1)+DR(2)*DR(2)+DR(3)*DR(3))
    THETA_TIP=SIGNAL*ACOS((ROTATION_MATRIX(1,1)*DR(1)+ROTATION_MATRIX(2,1)*DR(2)+ROTATION_MATRIX(3,1)*DR(3))/R_TIP)
    CALL WILLIAMS_DISPL_FUNCTIONS(R_TIP,THETA_TIP,Mu,Nu,PSI_TIP)
!
    DR_DQSIAL=0.D0
    DO I=1,3
        AUX1=0.D0
        DO J=1,NEN
            AUX1=AUX1+DPHI(J,DIR)*VALUES(J,I)
            AUXD1=AUXD1+DPHI(J,DIR)*VALUES(J,1)
            AUXD2=AUXD2+DPHI(J,DIR)*VALUES(J,2)
            AUXD3=AUXD3+DPHI(J,DIR)*VALUES(J,3)
        ENDDO
        AUX1=AUX1*DR(I)
        DR_DQSIAL=DR_DQSIAL+AUX1
    ENDDO
!
    DR_DQSIAL=DR_DQSIAL/R_TIP!*SIGNAL_QSI_AL_CT
    PSI_TIP_DERIV=(1.D0/(2*R_TIP))*PSI_TIP*DR_DQSIAL
!
    END SUBROUTINE WILLIAMS_DISPL_FUNCTIONS_DERIV