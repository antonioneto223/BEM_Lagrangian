SUBROUTINE TXTS_VISCO

    !
    !   SUBROUTINE TO INITIATE TXT'S THAT WILL CONTAIN THE TIME RESUTLS FOR VISCOELASTIC ANALYSIS
    !
    
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE VISCO
    USE REINFORCEMENTS
    
    INTEGER::I,J,KK,NEN
    REAL*8::DIST1(4),MENOR_DIST
    CHARACTER*100::FILE
    CHARACTER*1024::FORMAT_STRING
        
    !   ABRINDO TXT'S PARA SAIDA DE DADOS DE BOUND NODES E INT NODES
    DO I=1,N_NODES_BOUND_VISCO
        IF (I .LT. 10) THEN
            FORMAT_STRING="(A13,I1)"
        ELSE IF (I.LT.100) THEN
            FORMAT_STRING="(A13,I2)"
        ELSE IF (I.LT.1000) THEN
            FORMAT_STRING="(A13,I3)"
        ENDIF
        WRITE(FILE,FORMAT_STRING)"PTO_CONTORNO_",I
        
        OPEN(I+COUNT_TXT,FILE='Output_data/Visco/'//TRIM(FILE)//'.TXT',STATUS='UNKNOWN')
    ENDDO
    
    DO I=1,N_NODES_INT_VISCO
        IF (I .LT. 10) THEN
            FORMAT_STRING="(A11,I1)"
        ELSE IF (I.LT.100) THEN
            FORMAT_STRING="(A11,I2)"
        ELSE IF (I.LT.1000) THEN
            FORMAT_STRING="(A11,I3)"
        ENDIF
        WRITE(FILE,FORMAT_STRING)"PTO_INTERN_",I
        
        OPEN(I+2*COUNT_TXT,FILE='Output_data/Visco/'//TRIM(FILE)//'.TXT',STATUS='UNKNOWN')
    ENDDO
       
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        IF (VISCO_ANALYSIS_REINF) THEN
            OPEN(3*COUNT_TXT+1,FILE='Output_data/Visco/Reinforcements/AXIAL_FORCE_TOTAL_FIBERS.TXT',STATUS='UNKNOWN')
            OPEN(3*COUNT_TXT+2,FILE='Output_data/Visco/Reinforcements/UX_FIBERS_VISCO.TXT',STATUS='UNKNOWN')
            OPEN(3*COUNT_TXT+3,FILE='Output_data/Visco/Reinforcements/UY_FIBERS_VISCO.TXT',STATUS='UNKNOWN')
            OPEN(3*COUNT_TXT+4,FILE='Output_data/Visco/Reinforcements/UZ_FIBERS_VISCO.TXT',STATUS='UNKNOWN')
            OPEN(3*COUNT_TXT+5,FILE='Output_data/Visco/Reinforcements/AXIAL_FORCE_EL_FIBERS.TXT',STATUS='UNKNOWN')
            OPEN(3*COUNT_TXT+6,FILE='Output_data/Visco/Reinforcements/AXIAL_FORCE_VISCOSA_FIBERS.TXT',STATUS='UNKNOWN')
        ELSE
            OPEN(3*COUNT_TXT+1,FILE='Output_data/Visco/Reinforcements/AXIAL_FORCE_FIBERS.TXT',STATUS='UNKNOWN')
            OPEN(3*COUNT_TXT+2,FILE='Output_data/Visco/Reinforcements/UX_FIBERS_VISCO.TXT',STATUS='UNKNOWN')
            OPEN(3*COUNT_TXT+3,FILE='Output_data/Visco/Reinforcements/UY_FIBERS_VISCO.TXT',STATUS='UNKNOWN')
            OPEN(3*COUNT_TXT+4,FILE='Output_data/Visco/Reinforcements/UZ_FIBERS_VISCO.TXT',STATUS='UNKNOWN')
            !IF (REINFORCEMENTS_ELASTOPLASTIC .EQ. 'Y') THEN
            !    OPEN(18,FILE='SAIDA_DE_DADOS_VISCO/DEF_PLAST_FIBRA_VISCO.TXT')
            !ENDIF
        ENDIF
    ENDIF
    
    !   ENCONTRANDO PONTOS DE CONTORNO MAIS PROXIMOS DOS PONTOS SELECIONADOS
    CALL FIND_BOUND_POINT_VISCO
    
    !   CABEÇALHO DOS ARQUIVOS DE SAIDA DE DADOS VISCOELASTICO-----------------------------------
    DO I=1,N_NODES_BOUND_VISCO
        WRITE(I+COUNT_TXT,*)'TIME RESULTS FOR SELECTED BOUNDARY NODE'
        WRITE(I+COUNT_TXT,*)''
        IF (BOUND_PNT_VISCO_COLLOC(I)) THEN
            WRITE(I+COUNT_TXT,*)'NODE        X               Y             Z'
            WRITE(I+COUNT_TXT,'(i5,3f15.3)')COLLOCPOINTS_VISCO(I),COORD_COLLOCPOINTS(COLLOCPOINTS_VISCO(I),1),COORD_COLLOCPOINTS(COLLOCPOINTS_VISCO(I),2),COORD_COLLOCPOINTS(COLLOCPOINTS_VISCO(I),3)
        ELSE
            WRITE(I+COUNT_TXT,'(A)')'COORDINATES        X       Y       Z'
            WRITE(I+COUNT_TXT,'(3F15.3)')COORDS_BOUND_VISCO(I,1),COORDS_BOUND_VISCO(I,2),COORDS_BOUND_VISCO(I,3)
            WRITE(I+COUNT_TXT,*)''
            WRITE(I+COUNT_TXT,'(A,I7)')'ELEMENT MESH CONTAINING THIS NODE:',ELEMS_BOUND_PNT_VISCO(I)
            WRITE(I+COUNT_TXT,*)''
            WRITE(I+COUNT_TXT,*)'NODES POSITION IN THIS ELEMENT:'
            WRITE(I+COUNT_TXT,*)'NODE       X       Y       Z'
            NEN=ELEM_TYPE(ELEMS_BOUND_PNT_VISCO(I))*ORDER_ELEM(ELEMS_BOUND_PNT_VISCO(I))+&
                (ELEM_TYPE(ELEMS_BOUND_PNT_VISCO(I))-3)*(ORDER_ELEM(ELEMS_BOUND_PNT_VISCO(I))-1)*POL_FAMILY(ELEMS_BOUND_PNT_VISCO(I))
            DO J=1,NEN
                WRITE(I+COUNT_TXT,'(I5,3F15.3)')COLLOCPOINTS_CONNECTIVITY(ELEMS_BOUND_PNT_VISCO(I),J),&
                    COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEMS_BOUND_PNT_VISCO(I),J),1),&
                    COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEMS_BOUND_PNT_VISCO(I),J),2),&
                    COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(ELEMS_BOUND_PNT_VISCO(I),J),3)
            ENDDO
        ENDIF
        WRITE(I+COUNT_TXT,*)'' 
        WRITE(I+COUNT_TXT,*)''
        WRITE(I+COUNT_TXT,'(a)')'   TEMPO            Ux             Uy             Uz               Px               Py               Pz'  
    ENDDO
    DO I=1,N_NODES_INT_VISCO
        WRITE(I+2*COUNT_TXT,*)'TIME RESULTS FOR SELECTED INTERNAL NODE'
        WRITE(I+2*COUNT_TXT,*)''
        WRITE(I+2*COUNT_TXT,*)'NODE        X                Y               Z'
        WRITE(I+2*COUNT_TXT,'(i5,3f15.3)')I,COORDS_INT_VISCO(I,1),COORDS_INT_VISCO(I,2),COORDS_INT_VISCO(I,3)
        WRITE(I+2*COUNT_TXT,*)'' 
        WRITE(I+2*COUNT_TXT,'(a)')'   TEMPO             Ux             Uy               Uz             Sx-tot            Txy-tot            Txz-tot           Sy-tot            Tyz-tot          Sz-tot              &
        Sx-el            Txy-el            Txz-el             Sy-el            Tyz-el            Sz-el             Sx-vis            Txy-vis            Txz-vis           Sy-vis            Tyz-vis         Sz-vis' 
    ENDDO
    
    !   CABEÇALHO DOS ARQUIVOS DE SAIDA DE DADOS DOS PONTOS DE FIBRA ---------------------------
    IF (VISCO_ANALYSIS_REINF) THEN
        WRITE(3*COUNT_TXT+1,*)'TIME RESULTS FOR REINFORCEMENT TOTAL STRESS'
        WRITE(3*COUNT_TXT+1,*)''
        WRITE(3*COUNT_TXT+1,*)'NORMAL EFFORT TOTAL / TOTAL AXIAL FORCE FOR ALL NODES'
        WRITE(3*COUNT_TXT+1,*)''
        WRITE(3*COUNT_TXT+1,'(a)',advance='no')' TIME'
        DO i=1,N_NOS_MEF
            IF (I .EQ. 1) THEN
                WRITE(3*COUNT_TXT+1,'(i8)',advance='no')i
            ELSE
                WRITE(3*COUNT_TXT+1,'(i14)',advance='no')i
            ENDIF
        ENDDO
        WRITE(3*COUNT_TXT+1,*)
        
        WRITE(3*COUNT_TXT+5,*)'TIME RESULTS FOR REINFORCEMENT ELASTIC STRESS PORTION'
        WRITE(3*COUNT_TXT+5,*)''
        WRITE(3*COUNT_TXT+5,*)'NORMAL EFFORT ELASTIC / ELASTIC AXIAL FORCE FOR ALL NODES'
        WRITE(3*COUNT_TXT+5,*)''
        WRITE(3*COUNT_TXT+5,'(a)',advance='no')' TIME'
        DO i=1,N_NOS_MEF
            IF (I .EQ. 1) THEN
                WRITE(3*COUNT_TXT+5,'(i8)',advance='no')i
            ELSE
                WRITE(3*COUNT_TXT+5,'(i14)',advance='no')i
            ENDIF
        ENDDO
        WRITE(3*COUNT_TXT+5,*)
        
        WRITE(3*COUNT_TXT+6,*)'TIME RESULTS FOR REINFORCEMENT VISCOUS STRESS PORTION'
        WRITE(3*COUNT_TXT+6,*)''
        WRITE(3*COUNT_TXT+6,*)'NORMAL EFFORT VISCOUS / VISCOUS AXIAL FORCE FOR ALL NODES'
        WRITE(3*COUNT_TXT+6,*)''
        WRITE(3*COUNT_TXT+6,'(a)',advance='no')' TIME'
        DO i=1,N_NOS_MEF
            IF (I .EQ. 1) THEN
                WRITE(3*COUNT_TXT+6,'(i8)',advance='no')i
            ELSE
                WRITE(3*COUNT_TXT+6,'(i14)',advance='no')i
            ENDIF
        ENDDO
        WRITE(3*COUNT_TXT+6,*)
    ELSE
        WRITE(3*COUNT_TXT+1,*)'TIME RESULTS FOR REINFORCEMENT STRESS'
        WRITE(3*COUNT_TXT+1,*)''
        WRITE(3*COUNT_TXT+1,*)'NORMAL EFFORT / AXIAL FORCE FOR ALL NODES'
        WRITE(3*COUNT_TXT+1,*)''
        WRITE(3*COUNT_TXT+1,'(a)',advance='no')' TIME'
        DO i=1,N_NOS_MEF
            IF (I .EQ. 1) THEN
                WRITE(3*COUNT_TXT+1,'(i8)',advance='no')i
            ELSE
                WRITE(3*COUNT_TXT+1,'(i14)',advance='no')i
            ENDIF
        ENDDO
        WRITE(3*COUNT_TXT+1,*)
    ENDIF
    
    !WRITE(18,*)'TIME RESULTS FOR REINFORCEMENT PLASTIC STRAIN'
    !WRITE(18,*)''
    !WRITE(18,*)'PLASTIC STRAIN FOR ALL NODES'
    !WRITE(18,*)''
    !WRITE(18,'(a)',advance='no')' TIME'
    !KK = 0
    !DO i=1,N_ELEMENTOS_MEF
    !    DO j=1,ORDEM_MEF(I)+1
    !        IF (I .NE. 1) THEN
    !            IF((j .EQ. 1) .AND. (CONECTI_MEF(i,j) .EQ. CONECTI_MEF(i-1,ORDEM_MEF(i)+1)) ) CYCLE
    !        ENDIF
    !        KK = KK + 1
    !        IF(KK .EQ. 1) THEN
    !            WRITE(18,'(i8)',advance='no')KK
    !        ELSE
    !            WRITE(18,'(i14)',advance='no')KK
    !        ENDIF
    !    ENDDO
    !ENDDO
    !WRITE(18,*)
    
    WRITE(3*COUNT_TXT+2,*)'TIME RESULTS FOR REINFORCEMENT NODAL DISPLACEMENT X'
    WRITE(3*COUNT_TXT+2,*)''
    WRITE(3*COUNT_TXT+2,*)'Ux FOR ALL NODES'
    WRITE(3*COUNT_TXT+2,*)''
    WRITE(3*COUNT_TXT+2,'(a)',advance='no')' TIME'
    DO i=1,N_NOS_MEF
        IF (I .EQ. 1) THEN
            WRITE(3*COUNT_TXT+2,'(i8)',advance='no')i
        ELSE
            WRITE(3*COUNT_TXT+2,'(i14)',advance='no')i
        ENDIF
    ENDDO
    WRITE(3*COUNT_TXT+2,*)
    
    WRITE(3*COUNT_TXT+3,*)'TIME RESULTS FOR REINFORCEMENT NODAL DISPLACEMENT y'
    WRITE(3*COUNT_TXT+3,*)''
    WRITE(3*COUNT_TXT+3,*)'Uy FOR ALL NODES'
    WRITE(3*COUNT_TXT+3,*)''
    WRITE(3*COUNT_TXT+3,'(a)',advance='no')' TIME'
    DO i=1,N_NOS_MEF
        IF (i .EQ. 1) THEN 
            WRITE(3*COUNT_TXT+3,'(i8)',advance='no')i
        ELSE
            WRITE(3*COUNT_TXT+3,'(i14)',advance='no')i
        ENDIF
    ENDDO
    WRITE(3*COUNT_TXT+3,*)
    
    WRITE(3*COUNT_TXT+4,*)'TIME RESULTS FOR REINFORCEMENT NODAL DISPLACEMENT Z'
    WRITE(3*COUNT_TXT+4,*)''
    WRITE(3*COUNT_TXT+4,*)'Uz FOR ALL NODES'
    WRITE(3*COUNT_TXT+4,*)''
    WRITE(3*COUNT_TXT+4,'(a)',advance='no')' TIME'
    DO i=1,N_NOS_MEF
        IF (i .EQ. 1) THEN 
            WRITE(3*COUNT_TXT+4,'(i8)',advance='no')i
        ELSE
            WRITE(3*COUNT_TXT+4,'(i14)',advance='no')i
        ENDIF
    ENDDO
    WRITE(3*COUNT_TXT+4,*)
    ! ----------------------------------------------------------------------------------------------
    
END SUBROUTINE TXTS_VISCO
    
    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    
 
    
    
SUBROUTINE BOUND_POINT_VISCO_OUTPUT(TIME_STEP)
    
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE VISCO
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER(4),INTENT(IN)::TIME_STEP
    INTEGER(4)::GL1,GL2,GL3,i,j,k,kk,CONT,NEN,ELEM_A
    REAL*8::AUX,AUX1,PHI[ALLOCATABLE](:),D2PHI[ALLOCATABLE](:,:,:),U_A(3),QSI_A(2),T_A(3)
    REAL*8,DIMENSION(:,:),ALLOCATABLE::VALUES,COEFFICIENTS,DPHI
    
    DO I=1,N_NODES_BOUND_VISCO
        
        IF (BOUND_PNT_VISCO_COLLOC(I)) THEN ! COINCIDENTE COM PONTO DE COLOCACAO
        
            GL1= 3*COLLOCPOINTS_VISCO(I)-2
            GL2= 3*COLLOCPOINTS_VISCO(I)-1
            GL3= 3*COLLOCPOINTS_VISCO(I)
            WRITE(I+COUNT_TXT,'(F7.1,6E18.6)')DFLOAT(TIME_STEP)*DELTA_T_VISCO,U(GL1),U(GL2),U(GL3),T(GL1),T(GL2),T(GL3)
        
        ELSE ! DENTRO DE UM ELEMENTO
            
            ELEM_A=ELEMS_BOUND_PNT_VISCO(I)
            QSI_A(:)=QSIS_BOUND_PNT_VISCO(I,:)
            
            NEN=ELEM_TYPE(ELEM_A)*ORDER_ELEM(ELEM_A)+(ELEM_TYPE(ELEM_A)-3)*(ORDER_ELEM(ELEM_A)-1)*POL_FAMILY(ELEM_A)
            ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2))
    
            DO J=1,NEN
                VALUES(J,:)=COORD_NODES(NODES_CONNECTIVITY(ELEM_A,J),:)
            ENDDO
        
            CALL SHAPE_FUNCTIONS_COEFFICIENTS(ELEM_A,NEN,COEFFICIENTS)
            CALL SHAPE_FUNCTIONS(2,QSI_A(1),QSI_A(2),NEN,VALUES,COEFFICIENTS,AUX,AUX,AUX,PHI,DPHI,D2PHI)
            
            U_A=0.D0
            T_A=0.0D0
            DO K=1,3
                DO J=1,NEN
                    U_A(K)=U_A(K)+U(3*COLLOCPOINTS_CONNECTIVITY(ELEM_A,J)+K-3)*PHI(J)
                    T_A(K)=T_A(K)+T(3*COLLOCPOINTS_CONNECTIVITY(ELEM_A,J)+K-3)*PHI(J)
                ENDDO
            ENDDO
            
            WRITE(I+COUNT_TXT,'(F7.1,6E18.6)')DFLOAT(TIME_STEP)*DELTA_T_VISCO,U_A(1),U_A(2),U_A(3),T_A(1),T_A(2),T_A(3)
            
            DEALLOCATE(VALUES,COEFFICIENTS,PHI,DPHI,D2PHI)
            
        ENDIF
        
    ENDDO
                 
    IF(EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN    
        WRITE(3*COUNT_TXT+1,'(F7.1)',advance='no')DFLOAT(TIME_STEP)*DELTA_T_VISCO
        WRITE(3*COUNT_TXT+2,'(F7.1)',advance='no')DFLOAT(TIME_STEP)*DELTA_T_VISCO
        WRITE(3*COUNT_TXT+3,'(F7.1)',advance='no')DFLOAT(TIME_STEP)*DELTA_T_VISCO
        WRITE(3*COUNT_TXT+4,'(F7.1)',advance='no')DFLOAT(TIME_STEP)*DELTA_T_VISCO
        
        IF (VISCO_ANALYSIS_REINF) THEN
            WRITE(3*COUNT_TXT+5,'(F7.1)',advance='no')DFLOAT(TIME_STEP)*DELTA_T_VISCO
            WRITE(3*COUNT_TXT+6,'(F7.1)',advance='no')DFLOAT(TIME_STEP)*DELTA_T_VISCO
        ENDIF
            
        DO I=1,N_NOS_MEF
            WRITE(3*COUNT_TXT+1,'(E14.6)',advance='no')E_NORMAL_NOD(I)
            WRITE(3*COUNT_TXT+2,'(E14.6)',advance='no')U_MEF(3*I-2)
            WRITE(3*COUNT_TXT+3,'(E14.6)',advance='no')U_MEF(3*I-1)
            WRITE(3*COUNT_TXT+4,'(E14.6)',advance='no')U_MEF(3*I)
            IF (VISCO_ANALYSIS_REINF) THEN
                !WRITE(3*COUNT_TXT+5,'(E14.6)',advance='no')E_NORMAL_NOD_EL(I)
                !WRITE(3*COUNT_TXT+6,'(E14.6)',advance='no')E_NORMAL_NOD_VISCOSA(I)
            ENDIF
        ENDDO
        WRITE(3*COUNT_TXT+1,*)
        WRITE(3*COUNT_TXT+2,*)
        WRITE(3*COUNT_TXT+3,*)
        WRITE(3*COUNT_TXT+4,*)
        IF (VISCO_ANALYSIS_REINF) THEN
            WRITE(3*COUNT_TXT+5,*)
            WRITE(3*COUNT_TXT+6,*)
        ENDIF
        
        !IF (REINFORCEMENTS_ELASTOPLASTIC .EQ. 'Y') THEN  ! SAIDA DE DADOS PARA PLASTICIDADE
        !    
        !    WRITE(18,'(i5)',advance='no')PASSO_TEMPO
        !          
        !    CONT=0
        !    KK = 0
        !    DO i=1,N_ELEMENTOS_MEF
        !        DO j=1,ORDEM_MEF(I)+1
        !            KK = KK + 1
        !            AUX1 = GRAND_ACUM_NOD(KK,2)
        !            IF (I .NE. 1) THEN
        !                IF((j .EQ. 1) .AND. (CONECTI_MEF(i,j) .EQ. CONECTI_MEF(i-1,ORDEM_MEF(i)+1)) ) THEN
        !                    CYCLE
        !                    CONT=1
        !                ELSE
        !                    CONT=0
        !                ENDIF
        !            ENDIF
        !            IF (CONT.EQ.1) AUX1=(AUX1+GRAND_ACUM_NOD(KK-1,2))/2.0D0                    
        !            WRITE(18,'(E15.5)',advance='no')AUX1
        !        ENDDO
        !    ENDDO
        !    WRITE(18,*)
        !    
        !ENDIF
        
    ENDIF
    

END SUBROUTINE BOUND_POINT_VISCO_OUTPUT
    
    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    

    
    
SUBROUTINE INT_POINT_VISCO_OUTPUT(TIME_STEP)
        
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE VISCO
    
    IMPLICIT NONE
    
    INTEGER(4),INTENT(IN)::TIME_STEP
    INTEGER(4)::i,j
    
    DO I=1,N_NODES_INT_VISCO
                
        WRITE(I+2*COUNT_TXT,'(F7.1,21E18.6)')DFLOAT(TIME_STEP)*DELTA_T_VISCO,U_INT_VISCO(3*I-2),U_INT_VISCO(3*I-1),U_INT_VISCO(3*I),&
            S_INT_VISCO(9*I-8),S_INT_VISCO(9*I-7),S_INT_VISCO(9*I-6),S_INT_VISCO(9*I-4),S_INT_VISCO(9*I-3),S_INT_VISCO(9*I),&
            S_EL(9*I-8),S_EL(9*I-7),S_EL(9*I-6),S_EL(9*I-4),S_EL(9*I-3),S_EL(9*I),S_VISCOSA(9*I-8),S_VISCOSA(9*I-7),S_VISCOSA(9*I-6),&
            S_VISCOSA(9*I-4),S_VISCOSA(9*I-3),S_VISCOSA(9*I)
        
    ENDDO
    
END SUBROUTINE INT_POINT_VISCO_OUTPUT
    
    
    
     
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    

    
    
    
SUBROUTINE ACADVIEW_OUTPUT_VISCO(TEMPO)
!
	USE ISOPARAMETRIC_MESH
	USE CRACKS_DUAL_BEM
    USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
	USE REMESHING
    USE REINFORCEMENTS
!
	IMPLICIT NONE 
!
    INTEGER,INTENT(IN)::TEMPO
!
	INTEGER::I,J,K,II,JJ,KK,L,NEN,C,N_CRACKED_ELEM,CONT,NNUM,NNODES
!
    CHARACTER*50::FILE
    CHARACTER*100::FORMAT_STRING
!	
	REAL*8::QSI1,QSI2,U_AUX(3),T_AUX(3),D2PHI[ALLOCATABLE](:,:,:),X,Y,Z,DX,DY,DZ,R,R_S(3,3),R_H(3,3)
!    
    REAL*8,DIMENSION(:),ALLOCATABLE::U_NODES,T_NODES,PHI
    REAL*8,DIMENSION(:,:),ALLOCATABLE::VALUES_U,VALUES_T,COEFFICIENTS,DPHI
!    
    ALLOCATE(U_NODES(3*N_COLLOCPOINTS),T_NODES(3*N_COLLOCPOINTS))
!------------------------------------------------------------------------------------------------------------------------------------------------------
!         U AND T GLOBAL COORDINATE SYSTEM
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
!    DO I=1,N_COLLOCPOINTS
!	    IF(DUAL_BEM(I).EQ."S")THEN
!	        DO J=1,N_ELEM
!		        NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
!			    DO K=1,NEN
!				    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
!					    JJ=J
!					    KK=K
!				    ENDIF
!			    ENDDO
!			ENDDO
!!			    
!	        CALL LOCAL_COODINATE_SYSTEM(KK,JJ,R_S,R_H)
!	        CALL DLINRG(3,R_S,3,R_S,3)
!!	        
!	        DO K=1,3
!	            U_AUX(K)=R_S(K,1)*U(3*I-2)+R_S(K,2)*U(3*I-1)+R_S(K,3)*U(3*I)
!	        	T_AUX(K)=R_S(K,1)*T(3*I-2)+R_S(K,2)*T(3*I-1)+R_S(K,3)*T(3*I) 
!	        ENDDO
!	        U(3*I-2)=U_AUX(1) 
!	        U(3*I-1)=U_AUX(2) 
!	        U(3*I)=U_AUX(3) 
!	        T(3*I-2)=T_AUX(1) 
!	        T(3*I-1)=T_AUX(2) 
!	        T(3*I)=T_AUX(3) 
!!
!	        DO J=1,N_COLLOCPOINTS
!		        IF(DUAL_BEM(J).EQ."H")THEN
!				    DX=COORD_COLLOCPOINTS(I,1)-COORD_COLLOCPOINTS(J,1)
!				    DY=COORD_COLLOCPOINTS(I,2)-COORD_COLLOCPOINTS(J,2)
!				    DZ=COORD_COLLOCPOINTS(I,3)-COORD_COLLOCPOINTS(J,3)
!				    R=DSQRT(DX*DX+DY*DY+DZ*DZ)
!				    IF(R.LE.1.E-12)THEN
!				        II=J
!				    ENDIF
!				ENDIF
!		    ENDDO
!!	        
!	        CALL DLINRG(3,R_H,3,R_H,3)
!!	        
!	        DO K=1,3
!	            U_AUX(K)=R_H(K,1)*U(3*II-2)+R_H(K,2)*U(3*II-1)+R_H(K,3)*U(3*II)
!	        	T_AUX(K)=R_H(K,1)*T(3*II-2)+R_H(K,2)*T(3*II-1)+R_H(K,3)*T(3*II) 
!	        ENDDO
!	        U(3*II-2)=U_AUX(1) 
!	        U(3*II-1)=U_AUX(2) 
!	        U(3*II)=U_AUX(3) 
!	        T(3*II-2)=T_AUX(1) 
!	        T(3*II-1)=T_AUX(2) 
!	        T(3*II)=T_AUX(3) 		    	        	        	           	   	    
!	    ENDIF
!	ENDDO	
!    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!         U AND T GLOBAL COORDINATE SYSTEM
!------------------------------------------------------------------------------------------------------------------------------------------------------		
!    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   EXTRAPOLATING DISPLACEMENTS AND TRACTIONS
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
    U_NODES=0.D0
    T_NODES=0.D0
    DO I=1,N_COLLOCPOINTS
    	C=0
	    DO J=1,N_ELEM
	        IF(CRACKED_ELEM(J).EQ."UNCRACKED")THEN
	            NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
	            DO K=1,NEN
	                IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
	                    C=C+1
	                    IF(C.LE.1)THEN
	                        ALLOCATE(VALUES_U(NEN,3),VALUES_T(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN))
	                        DO L=1,NEN
		                        VALUES_U(L,1)=U(3*COLLOCPOINTS_CONNECTIVITY(J,L)-2)
		                        VALUES_U(L,2)=U(3*COLLOCPOINTS_CONNECTIVITY(J,L)-1)
		                        VALUES_U(L,3)=U(3*COLLOCPOINTS_CONNECTIVITY(J,L))
!
		                        VALUES_T(L,1)=T(3*COLLOCPOINTS_CONNECTIVITY(J,L)-2)
		                        VALUES_T(L,2)=T(3*COLLOCPOINTS_CONNECTIVITY(J,L)-1)
		                        VALUES_T(L,3)=T(3*COLLOCPOINTS_CONNECTIVITY(J,L))		                    
	                        ENDDO 
	                        IF(ELEM_TYPE(J).EQ.3) THEN
	                            SELECTCASE(K)
	                            CASE(1)
	                                QSI1=1.D0
	                                QSI2=0.D0
	                            CASE(2)
	                                QSI1=0.D0
	                                QSI2=1.D0	                
	                            CASE(3)
	                                QSI1=0.D0
	                                QSI2=0.D0	                
	                            CASE(4)
	                                QSI1=0.5D0
	                                QSI2=0.5D0	                
	                            CASE(5)
	                                QSI1=0.D0
	                                QSI2=0.5D0	                
	                            CASE(6)
	                                QSI1=0.5D0
	                                QSI2=0.D0	                
	                            ENDSELECT
	                        ELSE
	                            SELECTCASE(K)
	                            CASE(1)
	                                QSI1=-1.D0
	                                QSI2=-1.D0
	                            CASE(2)
	                                QSI1=1.D0
	                                QSI2=-1.D0	                
	                            CASE(3)
	                                QSI1=1.D0
	                                QSI2=1.D0	                
	                            CASE(4)
	                                QSI1=-1.D0
	                                QSI2=1.D0	                
	                            CASE(5)
	                                QSI1=0.D0
	                                QSI2=-1.D0	                
	                            CASE(6)
	                                QSI1=1.D0
	                                QSI2=0.D0	                
	                            CASE(7)
	                                QSI1=0.D0
	                                QSI2=1.D0	                
	                            CASE(8)
	                                QSI1=-1.D0
	                                QSI2=0.D0	                
	                            CASE(9)
	                                QSI1=0.D0
	                                QSI2=0.D0	                
	                            ENDSELECT	            
	                        ENDIF
	                        CALL SHAPE_FUNCTIONS_COEFFICIENTS(J,NEN,COEFFICIENTS)
	                        CALL SHAPE_FUNCTIONS(2,QSI1,QSI2,NEN,VALUES_U,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
	                        DO L=1,NEN
	                            U_NODES(3*I-2)=U_NODES(3*I-2)+PHI(L)*VALUES_U(L,1)    
	                            U_NODES(3*I-1)=U_NODES(3*I-1)+PHI(L)*VALUES_U(L,2)  
	                            U_NODES(3*I)=U_NODES(3*I)+PHI(L)*VALUES_U(L,3) 
!
	                            T_NODES(3*I-2)=T_NODES(3*I-2)+PHI(L)*VALUES_T(L,1)    
	                            T_NODES(3*I-1)=T_NODES(3*I-1)+PHI(L)*VALUES_T(L,2)  
	                            T_NODES(3*I)=T_NODES(3*I)+PHI(L)*VALUES_T(L,3)	                   
	                        ENDDO
	                        DEALLOCATE(VALUES_U,VALUES_T,COEFFICIENTS,PHI)
	                    ENDIF
	                ENDIF
	            ENDDO
	        ENDIF
	    ENDDO    
    ENDDO
!    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   EXTRAPOLATING DISPLACEMENTS AND TRACTIONS
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
	N_CRACKED_ELEM=0
	DO I=1,N_ELEM
	    IF(CRACKED_ELEM(I).NE."UNCRACKED")THEN
	        N_CRACKED_ELEM=N_CRACKED_ELEM+1
	    ENDIF
    ENDDO
    
    !IF (TEMPO .LT. 10) THEN
    !    FORMAT_STRING="(A16,I1)"
    !ELSE IF (TEMPO.LT.100) THEN
    !    FORMAT_STRING="(A16,I2)"
    !ELSE IF (TEMPO.LT.1000) THEN
    !    FORMAT_STRING="(A16,I3)"
    !ENDIF
    FORMAT_STRING="(A16,I0)"
    WRITE(FILE,FORMAT_STRING)"BOUNDARY_OUTPUT_",TEMPO
        
     ! ------------------------------------------------------------------------------------------------------
    OPEN(80,FILE='Output_data\Paraview\'//TRIM(FILE)//'.vtk',STATUS='UNKNOWN')
    
    WRITE(80,'(a)')'# vtk DataFile Version 2.0'
    WRITE(80,'(a)')'Viscoelastic boundary data'
    WRITE(80,'(a)')'ASCII'
    WRITE(80,'(a)')
    WRITE(80,'(a)')'DATASET UNSTRUCTURED_GRID'
    nnodes=N_COLLOCPOINTS
    WRITE(80,'(a,1x,i0,1x,a)')'POINTS ',nnodes,' float'
	CONT=0
	DO I=1,N_COLLOCPOINTS
	    C=0
	    DO J=1,N_ELEM
	        IF(CRACKED_ELEM(J).EQ."UNCRACKED")THEN
	            NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
	            DO K=1,NEN
	                IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
	                    C=C+1
	                    IF(C.LE.1)THEN
	                        CONT=CONT+1
	                        !WRITE(80,'(3(f0.10,1x))')COORD_NODES(NODES_CONNECTIVITY(J,K),1)+U_NODES(3*NODES_CONNECTIVITY(J,K)-2),&
                         !       COORD_NODES(NODES_CONNECTIVITY(J,K),2)+U_NODES(3*NODES_CONNECTIVITY(J,K)-1),&
                         !       COORD_NODES(NODES_CONNECTIVITY(J,K),3)+U_NODES(3*NODES_CONNECTIVITY(J,K))
                            WRITE(80,'(3(f0.10,1x))')COORD_NODES(NODES_CONNECTIVITY(J,K),1),COORD_NODES(NODES_CONNECTIVITY(J,K),2),&
                                COORD_NODES(NODES_CONNECTIVITY(J,K),3)
	                    ENDIF
	                ENDIF
	            ENDDO
	        ENDIF    
	    ENDDO
    ENDDO
    nnum=0
    DO J=1,N_ELEM
        NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
        NNUM = NNUM + NEN +1
    ENDDO  
    WRITE(80,'(a,1x,i0,1x,i0)')'CELLS ',N_ELEM,nnum
    
	DO I=1,N_ELEM
	    IF(CRACKED_ELEM(I).EQ."UNCRACKED")THEN	
	        NEN=ELEM_TYPE(I)*ORDER_ELEM(I)+(ELEM_TYPE(I)-3)*(ORDER_ELEM(I)-1)*POL_FAMILY(I)
            WRITE(80,'(100(I0,1X))')NEN,(COLLOCPOINTS_CONNECTIVITY(I,J)-1,J=1,NEN)
        ENDIF
    ENDDO
    
    WRITE(80,'(a,1x,i0)')'CELL_TYPES ',N_ELEM
    DO I=1,N_ELEM
        IF(CRACKED_ELEM(I).EQ."UNCRACKED")THEN	
            NEN=ELEM_TYPE(I)*ORDER_ELEM(I)+(ELEM_TYPE(I)-3)*(ORDER_ELEM(I)-1)*POL_FAMILY(I)
            SELECT CASE(NEN)
            CASE(3,6)
                WRITE(80,'(I0)') 69
            CASE (4,9)
                WRITE(80,'(I0)') 70
            CASE (8)
                WRITE(80,'(I0)') 23  ! NAO VALIDADO
            END SELECT
        ENDIF
    ENDDO
    
    WRITE(80,'(a,1x,i0)')'POINT_DATA ',nnodes
    WRITE(80,'(a,1x,i0)')'SCALARS DISPLACEMENT FLOAT ',3
    WRITE(80,'(a)')'LOOKUP_TABLE default'
	DO I=1,N_COLLOCPOINTS
	    WRITE(80,'(4(f0.10,1x))')U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I)   
    ENDDO
               
    !WRITE(80,'(a,1x,i0)')'POINT_DATA ',N_COLLOCPOINTS
    WRITE(80,'(a,1x,i0)')'SCALARS TRACTIONS FLOAT ',3
    WRITE(80,'(a)')'LOOKUP_TABLE default'
	DO I=1,N_COLLOCPOINTS
	    WRITE(80,'(3(f0.10,1x))')T_NODES(3*I-2),T_NODES(3*I-1),T_NODES(3*I)   
    ENDDO
!	
	CLOSE(80)
    
    
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        OPEN(80,file='Output_data\Paraview\FIBER_OUTPUT.1.vtk',status='unknown')
    !
        WRITE(80,'(a)')'# vtk DataFile Version 2.0'
        WRITE(80,'(a)')'Animation data'
        WRITE(80,'(a)')'ASCII'
        WRITE(80,'(a)')
        WRITE(80,'(a)')'DATASET UNSTRUCTURED_GRID'
        nnodes=N_NOS_MEF
        WRITE(80,'(a,1x,i0,1x,a)')'POINTS ',nnodes,' float'
	    CONT=0
        DO I=1,N_NOS_MEF
            IF (SLIP_REINF) THEN
                WRITE(80,'(3(f0.10,1x))')COORDNOS_MEF(I,1)+U_MEF(3*I-2)-S_ACUM(3*I-2),&
                    COORDNOS_MEF(I,2)+U_MEF(3*I-1)-S_ACUM(3*I-1),&
                    COORDNOS_MEF(I,3)+U_MEF(3*I)-S_ACUM(3*I)
            ELSE
                WRITE(80,'(3(f0.10,1x))')COORDNOS_MEF(I,1)+U_MEF(3*I-2),&
                    COORDNOS_MEF(I,2)+U_MEF(3*I-1),&
                    COORDNOS_MEF(I,3)+U_MEF(3*I)
            ENDIF
        ENDDO
        nnum=0
        IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
            DO I=1,N_ELEMENTOS_MEF
                NNUM = NNUM + ORDEM_MEF(I)+2
            ENDDO
        ENDIF
        WRITE(80,'(a,1x,i0,1x,i0)')'CELLS ',N_ELEMENTOS_MEF,nnum
    
        DO I=1,N_ELEMENTOS_MEF
            NEN=ORDEM_MEF(I)+1
            WRITE(80,'(3(I0,1X))',ADVANCE='NO')NEN,CONECTI_MEF(I,1)-1,CONECTI_MEF(I,NEN)-1
            WRITE(80,'(100(I0,1X))')(CONECTI_MEF(I,J)-1,J=2,NEN-1)
        ENDDO
        
        WRITE(80,'(a,1x,i0)')'CELL_TYPES ',N_ELEMENTOS_MEF
        DO I=1,N_ELEMENTOS_MEF
            WRITE(80,'(I0)') 68
        ENDDO
    
        IF (SLIP_REINF) THEN
            WRITE(80,'(a,1x,i0)')'POINT_DATA ',nnodes
            WRITE(80,'(a,1x,i0)')'SCALARS DISPLACEMENT FLOAT ',3
            WRITE(80,'(a)')'LOOKUP_TABLE default'
            DO I=1,N_NOS_MEF
                WRITE(80,'(4(f0.10,1x))')U_MEF(3*I-2)-S_ACUM(3*I-2),U_MEF(3*I-1)-S_ACUM(3*I-1),U_MEF(3*I)-S_ACUM(3*I)
            ENDDO
        ELSE
            WRITE(80,'(a,1x,i0)')'POINT_DATA ',nnodes
            WRITE(80,'(a,1x,i0)')'SCALARS DISPLACEMENT FLOAT ',3
            WRITE(80,'(a)')'LOOKUP_TABLE default'
            DO I=1,N_NOS_MEF
                WRITE(80,'(4(f0.10,1x))')U_MEF(3*I-2),U_MEF(3*I-1),U_MEF(3*I)
            ENDDO
        ENDIF
                   
        !WRITE(80,'(a,1x,i0)')'POINT_DATA ',NNODES
        WRITE(80,'(a,1x,i0)')'SCALARS ADERHENCE_FORCE FLOAT ',3
        WRITE(80,'(a)')'LOOKUP_TABLE default'
        DO I=1,N_NOS_MEF
            WRITE(80,'((f0.10,1x))')P_MEF(3*I-2),P_MEF(3*I-1),P_MEF(3*I)
        ENDDO
        
        !WRITE(80,'(a,1x,i0)')'POINT_DATA ',NNODES
        WRITE(80,'(a,1x,i0)')'SCALARS NORMAL FLOAT ',1
        WRITE(80,'(a)')'LOOKUP_TABLE default'
        DO I=1,N_NOS_MEF
            WRITE(80,'((f0.10,1x))')E_NORMAL_NOD(I) 
        ENDDO
        
        IF (EP_ANALYSIS_REINF) THEN
            !WRITE(80,'(a,1x,i0)')'POINT_DATA ',NNODES
            WRITE(80,'(a,1x,i0)')'SCALARS PLASTIC_STRAIN FLOAT ',1
            WRITE(80,'(a)')'LOOKUP_TABLE default'
            DO I=1,N_NOS_MEF
                WRITE(80,'((f0.10,1x))')PL_STRAIN_NOD(I) 
            ENDDO
        ENDIF
    !	
	    CLOSE(80)
    ENDIF
    
10  FORMAT(1h#)
100	FORMAT(2x,i5,12X,i5,12X,i5)
200	FORMAT(8x,F14.9,16x,F14.9,16x,F14.9,16x,F14.9,16X,F14.9,16x,F14.9)
300	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7X,i5)
301	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5)
302	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5)
303	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5)
304 FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5)

305	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7X,i5,3x,i1)
306	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,3x,i1)
307	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,3x,i1)
308 FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,3x,i1)
309 FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,3x,i1)
    
400	FORMAT(8x,E13.6,4x,E13.6,4x,E13.6,15x,E13.6)
500 FORMAT(a,i4,a)
!    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!         U AND T LOCAL COORDINATE SYSTEM
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
!    DO I=1,N_COLLOCPOINTS
!	    IF(DUAL_BEM(I).EQ."S")THEN
!	        DO J=1,N_ELEM
!		        NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
!			    DO K=1,NEN
!				    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
!					    JJ=J
!					    KK=K
!				    ENDIF
!			    ENDDO
!			ENDDO
!!			    
!	        CALL LOCAL_COODINATE_SYSTEM(KK,JJ,R_S,R_H)
!!	        
!	        DO K=1,3
!	            U_AUX(K)=R_S(K,1)*U(3*I-2)+R_S(K,2)*U(3*I-1)+R_S(K,3)*U(3*I)
!	        	T_AUX(K)=R_S(K,1)*T(3*I-2)+R_S(K,2)*T(3*I-1)+R_S(K,3)*T(3*I) 
!	        ENDDO
!	        U(3*I-2)=U_AUX(1) 
!	        U(3*I-1)=U_AUX(2) 
!	        U(3*I)=U_AUX(3) 
!	        T(3*I-2)=T_AUX(1) 
!	        T(3*I-1)=T_AUX(2) 
!	        T(3*I)=T_AUX(3) 
!!
!	        DO J=1,N_COLLOCPOINTS
!		        IF(DUAL_BEM(J).EQ."H")THEN
!				    DX=COORD_COLLOCPOINTS(I,1)-COORD_COLLOCPOINTS(J,1)
!				    DY=COORD_COLLOCPOINTS(I,2)-COORD_COLLOCPOINTS(J,2)
!				    DZ=COORD_COLLOCPOINTS(I,3)-COORD_COLLOCPOINTS(J,3)
!				    R=DSQRT(DX*DX+DY*DY+DZ*DZ)
!				    IF(R.LE.1.E-12)THEN
!				        II=J
!				    ENDIF
!				ENDIF
!		    ENDDO
!!	        
!	        DO K=1,3
!	            U_AUX(K)=R_H(K,1)*U(3*II-2)+R_H(K,2)*U(3*II-1)+R_H(K,3)*U(3*II)
!	        	T_AUX(K)=R_H(K,1)*T(3*II-2)+R_H(K,2)*T(3*II-1)+R_H(K,3)*T(3*II) 
!	        ENDDO
!	        U(3*II-2)=U_AUX(1) 
!	        U(3*II-1)=U_AUX(2) 
!	        U(3*II)=U_AUX(3) 
!	        T(3*II-2)=T_AUX(1) 
!	        T(3*II-1)=T_AUX(2) 
!	        T(3*II)=T_AUX(3) 		    	        	        	           	   	    
!	    ENDIF
!	ENDDO		
!    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!         U AND T LOCAL COORDINATE SYSTEM
!------------------------------------------------------------------------------------------------------------------------------------------------------
!	
END SUBROUTINE ACADVIEW_OUTPUT_VISCO
