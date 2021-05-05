SUBROUTINE CONCENTRATED_F_1DBEM(VET,SIZ)
    
    !
    !   CALCULATE 1DBEM ELEMENT APPLIED FORCE VECTOR 
    !   VET: STORE THE RESULT VECTOR, SIZ: SIZE OF VET
    !
    
    USE REINFORCEMENTS
    
    IMPLICIT NONE

    INTEGER,INTENT(IN)::SIZ
    REAL*8,INTENT(OUT)::VET(SIZ)
    LOGICAL::FINDD
    INTEGER::I,J,K,II,JJ,NO_G,NO_L,ELEM,QSI_P,N_NO
    REAL*8,DIMENSION(:,:),ALLOCATABLE::B_MATRIX,B_MATRIX_INV,VALUES,QSIW_MEF,MR,MAT,MAT1,MR_AUX,VET_LOCAL
    REAL*8,DIMENSION(:),ALLOCATABLE::D_PHI_MEF,VET_AUX
    REAL*8::X,Y,Z,TANG(4),L_MEF,AUX1,AUX2,AUX3,AUX4,F(3)
    
    VET=0.0D0
    
    DO I=1,N_NODAL_FORCES
        FINDD=.FALSE.
        NO_G=NODAL_FORCES_PLACE(I,1)
        DO J=1,N_ELEMENTOS_MEF
            DO K=1,ORDEM_MEF(J)+1
                IF (CONECTI_MEF(J,K).EQ.NO_G) THEN
                    ELEM=J
                    NO_L=K
                    FINDD=.TRUE.
                    EXIT
                ENDIF
            ENDDO
            IF (FINDD) EXIT
        ENDDO        
        !WRITE(*,*)'ENCONTROU'
        !WRITE(*,*)'ELEM:',ELEM,'NO LOCAL',NO_L,'NO GLOBAL:',NO_G
        !READ(*,*)
        
        ! INICIALIZANDO VARIAVEIS
        N_NO=ORDEM_MEF(ELEM)+1
        ALLOCATE(B_MATRIX(N_NO,N_NO),B_MATRIX_INV(N_NO,N_NO),VALUES(N_NO,3),D_PHI_MEF(N_NO),MR(N_NO,3*N_NO),&
            MR_AUX(1,3),VET_LOCAL(N_NO,1),MAT(3*N_NO,3),MAT1(N_NO,3),VET_AUX(3*N_NO))
        B_MATRIX=0.0D0
        MR=0.0D0
        MR_AUX=0.0D0
        B_MATRIX_INV=0.0D0
        VET_LOCAL=0.0D0
        
        ! PERFORMING NUMERICAL INTEGRATION FOR JACOBIAN
        II=CEILING((ORDEM_MEF(ELEM)+1+4+1)/2.0D0)  
        ALLOCATE(QSIW_MEF(II,2))
        CALL GAUSS_POINTS(II,QSIW_MEF)
        L_MEF=0.0D0
        DO J=1,II 
            CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSIW_MEF(J,1),N_NO,ELEM,VALUES,X,Y,Z,D_PHI_MEF)    
            ! COMPUTING THE JACOBIAN (TAN_NORM AND 2/L_MEF) --------
            TANG(1)=0.D0
	        TANG(2)=0.D0
            TANG(3)=0.0D0
	        DO K=1,N_NO
		        TANG(1)=TANG(1)+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(ELEM,K),1) !TANGENT VECTORS
		        TANG(2)=TANG(2)+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(ELEM,K),2) !TANGENT VECTORS
                TANG(3)=TANG(3)+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(ELEM,K),3) !TANGENT VECTORS
            ENDDO
	        TANG(4)=DSQRT(TANG(1)*TANG(1)+TANG(2)*TANG(2)+TANG(3)*TANG(3)) !TANGENT VECTOR NORM, JACOBIAN OF THE POINT
	        L_MEF=L_MEF+TANG(4)*QSIW_MEF(J,2)  !THE LENGHT OF THE ELEMENT IS THE SUM OF JACOBIANS, WEIGHTED BY THE GAUSS POINT WEIGHTS
        ENDDO
                 
        ! CALCULANDO MATRIZ G E G-1
        DO K=1,N_NO
            AUX4=DABS(QSI_MEF(ELEM,1)-QSI_MEF(ELEM,K))*L_MEF/2.0D0               
            B_MATRIX(K,1)=-AUX4/(2.0d0*MP_MEF(ELEM,1)*MP_MEF(ELEM,2))
                
            AUX4=DABS(QSI_MEF(ELEM,N_NO)-QSI_MEF(ELEM,K))*L_MEF/2.0D0
            B_MATRIX(K,N_NO)=-AUX4/(2.0d0*MP_MEF(ELEM,1)*MP_MEF(ELEM,2))
        ENDDO
        DO K=2,N_NO-1
            B_MATRIX(K,K)=L_MEF/(MP_MEF(ELEM,1)*MP_MEF(ELEM,2))
        ENDDO
        CALL DLINRG(SIZE(B_MATRIX,2),B_MATRIX,SIZE(B_MATRIX,1),B_MATRIX_INV,SIZE(B_MATRIX_INV,1))
        
        ! Calculando vetor U:
        QSI_P=-1.0D0+(2.0D0/ORDEM_MEF(ELEM))*DFLOAT(NO_L-1)
        DO K=1,N_NO
            AUX4=DABS(QSI_P-QSI_MEF(ELEM,K))*L_MEF/2.0D0 
            VET_LOCAL(K,1)=-AUX4/(2.0d0*MP_MEF(ELEM,1)*MP_MEF(ELEM,2))
        ENDDO        
        
        ! CALCULA MR PARA FAZER A ROTACAO
        DO J=1,N_NO
            AUX4=-1.0D0+(2.0D0/ORDEM_MEF(ELEM))*DFLOAT(J-1)
            CALL FUNCOES_DE_FORMA_MEC_MEF(6,AUX4,N_NO,ELEM,VALUES,X,Y,Z,D_PHI_MEF)  ! MUDOU: MAIS COERENCIA (VALIDAR)
            AUX1=0.0D0
            AUX2=0.0D0
            AUX3=0.0D0
            DO K=1,N_NO
		        AUX1=AUX1+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(ELEM,K),1)
		        AUX2=AUX2+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(ELEM,K),2)
                AUX3=AUX3+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(ELEM,K),3)
            ENDDO
            AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
            AUX1=AUX1/AUX4
            AUX2=AUX2/AUX4
            AUX3=AUX3/AUX4
            !TO AVOID NUMERICAL ERRORS:            
            IF (ABS(AUX1).LT.TOLER) AUX1=0.D0
            IF (ABS(AUX2).LT.TOLER) AUX2=0.D0
            IF (ABS(AUX3).LT.TOLER) AUX3=0.D0
            MR(J,3*J-2) = AUX1
            MR(J,3*J-1) = AUX2
            MR(J,3*J) = AUX3
            IF (J.EQ.NO_L) THEN
                MR_AUX(1,1)=AUX1
                MR_AUX(1,2)=AUX2
                MR_AUX(1,3)=AUX3
            ENDIF
        ENDDO
                
        ! PEGANDO VALOR DA FORÇA APLICADA E COLOCANDO EM F
        F=0.0D0
        F(NODAL_FORCES_PLACE(I,2))=NODAL_FORCES_VALUE(I)
        
        ! FAZENDO MULTIPLICACOES:
        VET_LOCAL=MATMUL(B_MATRIX_INV,VET_LOCAL)
        MAT1=MATMUL(VET_LOCAL,MR_AUX)
        MAT=MATMUL(TRANSPOSE(MR),MAT1)
        VET_AUX=MATMUL(MAT,F)
        
        ! JOGANDO PARA O VETOR GLOBAL
        II=SIZ-3*N_NOS_MEF
        DO J=1,N_NO
            JJ=CONECTI_MEF(ELEM,J)
            VET(II+3*JJ-2)=VET(II+3*JJ-2)+VET_AUX(3*J-2)
            VET(II+3*JJ-1)=VET(II+3*JJ-1)+VET_AUX(3*J-1)
            VET(II+3*JJ)=VET(II+3*JJ)+VET_AUX(3*J)
        ENDDO
        
        DEALLOCATE(B_MATRIX,B_MATRIX_INV,VALUES,D_PHI_MEF,MR,MR_AUX,VET_LOCAL,MAT,MAT1,VET_AUX,QSIW_MEF)
    ENDDO    
    
END SUBROUTINE
    
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

    
    
SUBROUTINE CONCENTRATED_U_CONNECTION(SOL,SIZ)
    
    !
    !   CALCULATES SYSTEM TO APPLY DISPLACEMENT AT THE SECOND NODE OF THE CONECTION ELEMENT 
    !   VET: STORE THE RESULT VECTOR, SIZ: SIZE OF SYSTEM WITHOUT CONNECTION ELEMENT
    !
    
    USE ISOPARAMETRIC_MESH
    USE SUB_REGIONS_INTERFACES
    USE CRACKS_DUAL_BEM
    USE REINFORCEMENTS
    
    IMPLICIT NONE

    INTEGER,INTENT(IN)::SIZ
    REAL*8,INTENT(IN)::SOL(SIZ+3*N_NODAL_FORCES)
    LOGICAL::FINDD
    INTEGER::I,J,K,KK,II,JJ,ELEM,NO_L,NO_G,N_NO,CONT_F,INFO,NO_CONECT(2)
    REAL*8,DIMENSION(:,:,:),ALLOCATABLE::KLOCAL,GLOCAL
    REAL*8,DIMENSION(:,:),ALLOCATABLE::A_AUX,B_F
    REAL*8,DIMENSION(:),ALLOCATABLE::D_PHI_MEF,VET_AUX,F
    REAL*8::X,Y,Z,TANG(4),L_MEF,AUX1,AUX2,AUX3,AUX4,COORDENADAS(2,3)
    
    ! INICIALIZANDO VARIAVEIS
    N_NO=2
    II=N_NODAL_FORCES
    ALLOCATE(KLOCAL(II,3*N_NO,3*N_NO),GLOCAL(II,3*N_NO,3*N_NO),KGLOBAL_NEW(3*N_NOS_MEF+3*II,3*N_NOS_MEF+3*II),&
    GGLOBAL_NEW(3*N_NOS_MEF+3*II,3*N_NOS_MEF+3*II),B_F(SIZ+3*II,3*II),VET_AUX(3),F(3*II))
    CONT_F=0
    KLOCAL=0.0D0
    GLOCAL=0.0D0
    KGLOBAL_NEW=0.0D0
    GGLOBAL_NEW=0.0D0
    B_F=0.0D0
    F=0.0D0
    
    DO I=1,N_NODAL_FORCES
        FINDD=.FALSE.
        NO_G=NODAL_FORCES_PLACE(I,1)
        DO J=1,N_ELEMENTOS_MEF
            DO K=1,ORDEM_MEF(J)+1
                IF (CONECTI_MEF(J,K).EQ.NO_G) THEN
                    ELEM=J
                    NO_L=K
                    FINDD=.TRUE.
                    EXIT
                ENDIF
            ENDDO
            IF (FINDD) EXIT
        ENDDO        
        !WRITE(*,*)'ENCONTROU'
        !WRITE(*,*)'ELEM:',ELEM,'NO LOCAL',NO_L,'NO GLOBAL:',NO_G
        !READ(*,*)
                
        ! CALÇULANDO RIGIDEZ DO CONNECTION ELEMENT
        DO K=1,3
            COORDENADAS(1,K)=COORDPCOLOC_MEF(NODAL_FORCES_PLACE(I,1),K)
        ENDDO
        AUX1=COORDNOS_MEF(NODAL_FORCES_PLACE(I,1),1)-COORDPCOLOC_MEF(NODAL_FORCES_PLACE(I,1),1)
        AUX2=COORDNOS_MEF(NODAL_FORCES_PLACE(I,1),2)-COORDPCOLOC_MEF(NODAL_FORCES_PLACE(I,1),2)
        AUX3=COORDNOS_MEF(NODAL_FORCES_PLACE(I,1),3)-COORDPCOLOC_MEF(NODAL_FORCES_PLACE(I,1),3)
        AUX1=2.0D0*AUX1
        AUX2=2.0D0*AUX2
        AUX3=2.0D0*AUX3
        COORDENADAS(2,1)=COORDENADAS(1,1)+AUX1
        COORDENADAS(2,2)=COORDENADAS(1,2)+AUX2
        COORDENADAS(2,3)=COORDENADAS(1,3)+AUX3
        CALL CONNECTION_ELEMENT_FOR_F(KLOCAL(I,:,:),GLOCAL(I,:,:),COORDENADAS)
        
    ENDDO
    
    ! ADICIONADO RIGIDEZ DO CONECTION ELEMENT NAS KGLOBAL E GGLOBAL
    DO J=1,3*N_NOS_MEF
        DO I=1,3*N_NOS_MEF
            KGLOBAL_NEW(I,J)=KGLOBAL(I,J)
            GGLOBAL_NEW(I,J)=GGLOBAL(I,J)
        ENDDO
    ENDDO
    DO I=1,N_NODAL_FORCES
        NO_CONECT(1)=NODAL_FORCES_PLACE(I,1)
        NO_CONECT(2)=N_NOS_MEF+I
        DO J=1,2
            DO JJ=1,3
                DO K=1,2
                    DO KK=1,3
                        KGLOBAL_NEW(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)=&
                        KGLOBAL_NEW(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)+&
                        KLOCAL(I,3*(J-1)+JJ,3*(K-1)+KK)
                        
                        GGLOBAL_NEW(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)=&
                        GGLOBAL_NEW(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)+&
                        GLOCAL(I,3*(J-1)+JJ,3*(K-1)+KK)
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
                
    ! RECRIANDO MATRIZ AGLOBAL COM OS NÓS DO CONNECTION ELEMENT
    JJ=SIZ+3*N_NODAL_FORCES
    ALLOCATE(A_AUX(JJ,JJ))
    A_AUX=0.0D0
    DO J=1,3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
        DO I=1,3*N_COLLOCPOINTS+3*N_NOS_MEF
            A_AUX(I,J)=A_GLOBAL_SLIP_F(I,J)
        ENDDO
    ENDDO
        
    II=3*N_COLLOCPOINTS+3*N_NOS_MEF
    JJ=3*N_COLLOCPOINTS
    DO J=1,3*N_NOS_MEF
        DO I=1,3*N_NOS_MEF
            A_AUX((II+I),(JJ+J))=KGLOBAL_NEW(I,J)  ! PARTE DE KGLOBAL ORIGINAL
        ENDDO
        DO I=1,3*N_NODAL_FORCES
            A_AUX((3*N_NOS_MEF+II+I),(JJ+J))=KGLOBAL_NEW(3*N_NOS_MEF+I,J)  ! PARTE DE KGLOBAL CRUZADO INFERIOR
        ENDDO
    ENDDO
        
    II=3*N_COLLOCPOINTS+3*N_NOS_MEF
    JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
    DO J=1,3*N_NOS_MEF
        DO I=1,3*N_NOS_MEF
            A_AUX((II+I),(JJ+J))=GGLOBAL_NEW(I,J) ! PARTE DO GGLOBAL ORIGINAL
        ENDDO
        DO I=1,3*N_NODAL_FORCES
            A_AUX((3*N_NOS_MEF+II+I),(JJ+J))=GGLOBAL_NEW(3*N_NOS_MEF+I,J) ! PARTE DO GGLOBAL CRUZADO INFERIOR
        ENDDO
    ENDDO
    
    II=3*N_COLLOCPOINTS+3*N_NOS_MEF
    JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
    DO J=1,3*N_NODAL_FORCES
        DO I=1,3*N_NOS_MEF
            A_AUX((II+I),(JJ+J))=KGLOBAL_NEW(I,3*N_NOS_MEF+J) ! PARTE DO KGLOBAL CRUZADO SUPERIOR
        ENDDO
        DO I=1,3*N_NODAL_FORCES
            A_AUX((3*N_NOS_MEF+II+I),(JJ+J))=KGLOBAL_NEW(3*N_NOS_MEF+I,3*N_NOS_MEF+J) ! PARTE DO KGLOBAL DIAGONAL 
        ENDDO
    ENDDO
            
    ! CRIANDO MATRIZ B PARA MULTIPLICAR PELO VETOR F A DIREITA
    II=SIZ
    DO I=1,N_NODAL_FORCES
        DO K=1,3
            B_F(II+3*(I-1)+K,3*(I-1)+K)=1.0D0
            !B_F(II+3*(I-1)+K,3*(I-1)+K)=1000.0D0   
        ENDDO
    ENDDO
        
    ! TROCANDO COLUNAS PARA APLICAR DESLOCAMENTO PRESCRITO NO NÓ DO CONNECTION
    II=SIZ
    DO J=1,N_NODAL_FORCES
        DO I=1,SIZ+3*N_NODAL_FORCES
            VET_AUX(1)=A_AUX(I,II+3*J-2)
            VET_AUX(2)=A_AUX(I,II+3*J-1)
            VET_AUX(3)=A_AUX(I,II+3*J)
                                
            A_AUX(I,II+3*J-2)=-B_F(I,3*J-2)
            A_AUX(I,II+3*J-1)=-B_F(I,3*J-1)
            A_AUX(I,II+3*J)=-B_F(I,3*J)
                    
            B_F(I,3*J-2)=-VET_AUX(1)
            B_F(I,3*J-1)=-VET_AUX(2)
            B_F(I,3*J)=-VET_AUX(3)   
            !if (vet_aux(1).ne.0.0 ) then
            !    write(*,*)
            !endif
        ENDDO
        F(3*(J-1)+NODAL_FORCES_PLACE(J,2))=NODAL_FORCES_VALUE(J)
    ENDDO
        
    ! SALVANDO MATRIZ A_GLOBAL E MATRIZ QUE MULTIPLICA O NODAL_FORCES_VALUES
    JJ=SIZ+3*N_NODAL_FORCES
    DEALLOCATE(A_GLOBAL_SLIP_F,IPIV_A_GLOBAL_SLIP_F,A_GLOBAL_SLIP_F_LU)
    ALLOCATE(A_GLOBAL_SLIP_F(JJ,JJ),IPIV_A_GLOBAL_SLIP_F(JJ),A_GLOBAL_SLIP_F_LU(JJ,JJ),&
        B_CONC_F(SIZ+3*N_NODAL_FORCES,3*N_NODAL_FORCES))
    A_GLOBAL_SLIP_F=A_AUX
    B_CONC_F=B_F

    ! MULTIPLICANDO B_F E F
    II=3*N_NODAL_FORCES
    JJ=SIZ+II
    CALL DGEMV('N',JJ,II,1.0D0,B_F,JJ,F,1,0.0D0,SOL,1)
            
    ! INVERTENDO A MATRIZ
    CALL DGETRF(JJ,JJ,A_AUX,JJ,IPIV_A_GLOBAL_SLIP_F,INFO)
    A_GLOBAL_SLIP_F_LU=A_AUX   
    
    IF (INFO.NE.0) THEN
        WRITE(*,'(A,I7,A)')'*** PROBLEM IN THE SYSTEM SOLVER, LINE:',INFO,'*****'
        WRITE(*,'(A)')'VERIFY DGETRF MATRIXES IN THE CONCENTRATED_U_CONNECTION SUBROUTINE'
        WRITE(*,'(A)')'*********************************************************'
        READ(*,*)
    ENDIF
        
END SUBROUTINE