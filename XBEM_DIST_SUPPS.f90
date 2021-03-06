    SUBROUTINE XBEM_DIST_SUPPS
!
!
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE XBEM_FORCE_VARIABLES
!    
    IMPLICIT NONE
!    
    LOGICAL::SINGULAR_INT   
!    
    INTEGER::I,J,K,L,M,II,JJ,KK,LL,NEN,ELEMS,NODES,ELEM_NODE,&
    POS_COLLOCPOINT,N_GAUSS,N_DIV
!
    REAL*8::DX,DY,DZ,DR(3),U_AST(3,3), D_AST(9,3),R,&
    Mu,Nu,C1,C2,AUX1,KR(3,3),PI,&
    ETA(3),ETAS(3),MAT_ETAS(3,9),VJC(3),DVJC(3,2),JC,AUX_D(3,3),&
    QSI_PNT(2),TOLDIST,QSIW_GAUSS[ALLOCATABLE](:,:),N_P(3),TOL1,AUX_CTE2
    REAL*8::LADO
    REAL*8::THETA,THETAI[ALLOCATABLE](:),THETAF[ALLOCATABLE](:),&
    RHO_BOUND,RHO_BOUND_AUX1,RHO_BOUND_AUX2(4),QSI1,QSI2,RHO,AUX_CTE1 
    REAL*8::X,Y,Z,X_TEST,Y_TEST,Z_TEST,COEFFICIENTS_AREA[ALLOCATABLE](:,:),VALUES[ALLOCATABLE](:,:),&
    PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),QSI_COLLOC(2)
! 
    N_GAUSS=N_GAUSSPOINTS
    !N_GAUSS=2
    ALLOCATE(QSIW_GAUSS(N_GAUSS,2),COEFFICIENTS_AREA(4,4))
    CALL GAUSS_POINTS(N_GAUSS,QSIW_GAUSS)
    PI=DACOS(-1.D0)
    !F_ENRICHED=0.0D0
    TOLDIST=1.0E-5
    TOL1=1.D-3
!
    KR(1,1)=1.D0
	KR(1,2)=0.D0
    KR(1,3)=0.D0
	KR(2,1)=0.D0
	KR(2,2)=1.D0
	KR(2,3)=0.D0
	KR(3,1)=0.D0
	KR(3,2)=0.D0
	KR(3,3)=1.D0
!

!
END SUBROUTINE