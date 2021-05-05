SUBROUTINE REINF_INPUT(FILE_REINF)
        
    USE PRE_REINFORCEMENTS 
    
    IMPLICIT NONE
    
    CHARACTER::CASE_F
    
    CHARACTER*30,INTENT(IN)::FILE_REINF
    
1   CONTINUE
    WRITE(*,*)
    WRITE(*,'(A)')'FIBERS INPUT BY MESH (M) OR BY PRE REINFORCEMENTS FILE (F) OR (N) FOR NO FIBERS?'
    READ(*,*)CASE_F
    
    IF (CASE_F.EQ."M") THEN
        GO TO 9
    ELSE IF (CASE_F.EQ."N") THEN
        GO TO 9
    ELSE IF (CASE_F.EQ."F") THEN
        CALL PRE_REINF(FILE_REINF)
    ELSE 
        WRITE(*,*)'WRITE SOME VALID CHACTER!!!'
        GO TO 1
    ENDIF
        
9   CONTINUE
    
    ! DEALLOCATING VARIABLES
    IF (TIPO_F.EQ."A") THEN
        DEALLOCATE(COORD_NODES,CONECTI)
    ELSE IF (TIPO_F.EQ."D") THEN
        DEALLOCATE(INITIAL_POINTS,FINAL_POINTS,NUM_EL,FIBER_DOMAIN,MAT_PROP,MODEL,MATERIALS,COORD_NODES,CONECTI,FE_MATERIAL,FE_DOMAIN)
        IF (ELASTOPLASTIC.EQ.1) DEALLOCATE(MAT_PROP_EP)
    ENDIF
    
END SUBROUTINE

    
! *******************************************************************************************************************
! *******************************************************************************************************************
! *******************************************************************************************************************
! *******************************************************************************************************************
! *******************************************************************************************************************
    
    
SUBROUTINE PRE_REINF(FILE_REINF)

    USE PRE_REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,J,NODE,ELEM,IOS
    CHARACTER*30,INTENT(IN)::FILE_REINF
    
    OPEN(1,FILE='Input_data\FIBERS.TXT',STATUS='OLD')
    READ(1,*,IOSTAT=IOS)!'RANDOM FIBER (A) OR DISCRETE FIBER (D)?'
    READ(1,*,IOSTAT=IOS)TIPO_F
    CLOSE(1)
    
    IF (IOS.NE.0) THEN
        WRITE(*,*)'YOU FORGET THE (A) OR (D) IN THE BEGNING OF FIBERS FILE'
        WRITE(*,*)'TRY AGAIN, FOR MORE DETAILS SEE PRE_REINF SUBROUTINE'
        READ(*,*)
    ELSE
        IF (TIPO_F.EQ.'D') THEN
            CALL INPUT_DATA_F
            NODE=0
            ELEM=0
            DO I=1,N_FIBER
                CALL CREATE_FIBER(I,NODE,ELEM)
            ENDDO
            N_NODES=NODE
            N_ELEMS=ELEM
            CALL CREATE_FILE(FILE_REINF)
            WRITE(*,*)'WROTE THE FILE:',FILE_REINF
            WRITE(*,*)''
        ELSE IF (TIPO_F.EQ.'A') THEN
            CALL INPUT_DATA_ALEAT
            NODE=0
            ELEM=0
            CALL CREATE_RANDOM(NODE,ELEM)
            CALL CREATE_FILE_RANDOM(FILE_REINF)
            WRITE(*,*)'WROTE THE FILE:',FILE_REINF
            WRITE(*,*)''
        ELSE
            WRITE(*,*)'ERROR IN THE FIBERS FILE (NOT A OR D IN THE BEGING)'
            WRITE(*,*)'TRY AGAIN, FOR MORE DETAILS SEE PRE_REINF SUBROUTINE'
            READ(*,*)
        ENDIF
    ENDIF
    
END SUBROUTINE
        
        
! *******************************************************************************************************************
! *******************************************************************************************************************
! *******************************************************************************************************************
! *******************************************************************************************************************
! *******************************************************************************************************************
  
    
SUBROUTINE INPUT_DATA_F
    
    USE PRE_REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,J,MAX_NUM_EL,IOS
    
    OPEN(1,FILE='Input_data\FIBERS.TXT',STATUS='UNKNOWN')
    
    READ(1,*)!RANDOM FIBER (A) OR DISCRETE FIBER (D)
    READ(1,*)!A
    READ(1,*)!
    READ(1,*)!NUMBER OF FIBERS
    READ(1,*)N_FIBER
    READ(1,*)!
    ALLOCATE(INITIAL_POINTS(N_FIBER,3),FINAL_POINTS(N_FIBER,3),NUM_EL(N_FIBER),FIBER_DOMAIN(N_FIBER))
    READ(1,*)!GIVE FIBERS POINTS
    DO I=1,N_FIBER
        READ(1,*)!FIBER NODE  X   Y    Z
        READ(1,*)J,INITIAL_POINTS(J,1),INITIAL_POINTS(J,2),INITIAL_POINTS(J,3)
        READ(1,*)J,FINAL_POINTS(J,1),FINAL_POINTS(J,2),FINAL_POINTS(J,3)
    ENDDO
    READ(1,*)!
    READ(1,*)! NUMBER OF DIFFERENT MATERIALS
    READ(1,*)N_MATERIALS
    ALLOCATE(MAT_PROP(N_MATERIALS,4),MODEL(N_MATERIALS),MATERIALS(N_FIBER))
    READ(1,*)!
    READ(1,*)! MAT PROPERTIES
    READ(1,*)! MAT    AREA    YOUNG     VISCOELASTICITY     MOD.BOLTZMANN       LEI CONST(1-KELVIN,2-BOLTZMANN,3-MAXWELL,4-HOOKE)
    DO I=1,N_MATERIALS
        READ(1,*)J,MAT_PROP(J,1),MAT_PROP(J,2),MAT_PROP(J,3),MAT_PROP(J,4),MODEL(J)
    ENDDO
    READ(1,*)!
    READ(1,*)! FIBERS MATERIALS
    READ(1,*)!FIBER    MATERIAL     DOMAIN
    DO I=1,N_FIBER
        READ(1,*)J,MATERIALS(J),FIBER_DOMAIN(J)
    ENDDO
    READ(1,*)!
    READ(1,*)! ELEMENTS ORDER
    READ(1,*)ORD
    READ(1,*)!
    READ(1,*)! FIBER NUMBER OF ELEMENTS
    DO I=1,N_FIBER
         READ(1,*)J,NUM_EL(J)
    ENDDO 
    
    READ(1,*)!
    READ(1,*)! ARE THE ANALYSIS ELASTOPLASTIC? (NO - 0 AND YES - 1)
    READ(1,*)ELASTOPLASTIC
    
    ! VERIFYING BOND-SLIP
    READ(1,*,IOSTAT=IOS)!
    READ(1,*,IOSTAT=IOS)! ARE THE ANALYSIS BOND-SLIP?
    READ(1,*,IOSTAT=IOS)BOND_SLIP
    
    IF (IOS.NE.0) BOND_SLIP=0
    
    IF (ELASTOPLASTIC.EQ.1) THEN
        ALLOCATE(MAT_PROP_EP(N_MATERIALS,2))
        READ(1,*)!
        READ(1,*)!************* IF ELASTOPLASTIC, GIVE PLASTIC PROPERTIES ***************
        READ(1,*)!MAT		SIGMA Y		KP
        DO I=1,N_MATERIALS
            READ(1,*)J,MAT_PROP_EP(J,1),MAT_PROP_EP(J,2)
        ENDDO
        READ(1,*)!
        READ(1,*)!NUMBER OF LOAD STEPS
        READ(1,*)LOAD_STEPS
        READ(1,*)!
        READ(1,*)!TOLERANCE
        READ(1,*)TOL
    ENDIF
    ! FINDING MAX NUM_EL
    MAX_NUM_EL=0
    DO I=1,N_FIBER
        IF (NUM_EL(I).GT.MAX_NUM_EL) MAX_NUM_EL=NUM_EL(I)
    ENDDO
    
    ! FINDING NUMBER OF NODES AND ELEMENTS
    N_NODES_MAX=0
    N_ELEMS=0
    DO I=1,N_FIBER
        DO J=1,NUM_EL(I)
            IF (J.EQ.1) N_NODES_MAX=N_NODES_MAX+(ORD+1)
            IF (J.NE.1) N_NODES_MAX=N_NODES_MAX+(ORD)
            N_ELEMS=N_ELEMS+1
        ENDDO
    ENDDO
    
    ALLOCATE(COORD_NODES(N_NODES_MAX,3),CONECTI(N_ELEMS,MAX_NUM_EL),FE_MATERIAL(N_ELEMS),FE_DOMAIN(N_ELEMS))
    
    CLOSE(1)
    
END SUBROUTINE
    
    
    
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! ***************************************************************************************************************************************************** 
   
    
SUBROUTINE CREATE_FIBER(FIBER,NODE,ELEM)
    
    USE PRE_REINFORCEMENTS
    
    IMPLICIT NONE
    
    LOGICAL::FIRST_NODE_FIBER
    INTEGER::FIBER,I,J,K,II,JJ,ELEM,NODE,GRAVOU
    REAL*8::DELTA(3),XI(3),EL_DELTA(3),COORD_AUX(3),AUX(4)
    
    GRAVOU=0
    
    DO J=1,3
        DELTA(J)=FINAL_POINTS(FIBER,J)-INITIAL_POINTS(FIBER,J)
        EL_DELTA(J)=DELTA(J)/(NUM_EL(FIBER))
    ENDDO
    
    DO I=1,NUM_EL(FIBER)
        
        ELEM=ELEM+1
        DO K=1,3
            XI(K)=INITIAL_POINTS(FIBER,K)+(I-1)*EL_DELTA(K)
        ENDDO
        
        IF (ELEM.EQ.308) THEN
            CONTINUE
        ENDIF
        
        DO J=1,ORD+1
        
            NODE=NODE+1
            IF (((J.EQ.1).AND.(I.NE.1)).OR.(GRAVOU.EQ.1)) NODE=NODE-1
            DO K=1,3
                COORD_AUX(K)=XI(K)+(J-1)*EL_DELTA(K)/(ORD)
            ENDDO
            
            ! VERIFICANDO SE JÁ TEM UM NÓ NESSA POSÇÃO
            AUX=0.0D0
            GRAVOU=0
            FIRST_NODE_FIBER=.FALSE.
            DO II=1,NODE-1
                DO K=1,3
                    AUX(K)=COORD_NODES(II,K)-COORD_AUX(K)
                ENDDO
                AUX(4)=DSQRT(AUX(1)**2.0D0+AUX(2)**2.0D0+AUX(3)**2.0D0)
                
                IF ((I.EQ.1).AND.(J.EQ.1)) FIRST_NODE_FIBER=.TRUE.
                !IF ((AUX(4).LT.TOLER).AND.(.NOT.FIRST_NODE_FIBER)) THEN
                IF (AUX(4).LT.TOLER) THEN
                    CONECTI(ELEM,J)=II
                    GRAVOU=1
                    EXIT
                ENDIF
            ENDDO
            
            IF (GRAVOU .EQ. 0) THEN
                DO K=1,3
                    COORD_NODES(NODE,K)=COORD_AUX(K)
                ENDDO
                CONECTI(ELEM,J)=NODE
            ENDIF
        ENDDO 
        
        FE_MATERIAL(ELEM)=MATERIALS(FIBER)
        FE_DOMAIN(ELEM)=FIBER_DOMAIN(FIBER)
        
    ENDDO
    
    IF (GRAVOU.EQ.1) NODE=NODE-1
    
END SUBROUTINE
    
    
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! ***************************************************************************************************************************************************** 
 

SUBROUTINE CREATE_FILE(FILE_REINF)

    USE PRE_REINFORCEMENTS
    
    IMPLICIT NONE
    
    CHARACTER*30,INTENT(IN)::FILE_REINF
    
    INTEGER::I,J
    
    OPEN(2,FILE='Input_data\'//TRIM(FILE_REINF),STATUS='UNKNOWN')

    WRITE(2,'(A)')'*************** BEM_FEM DATA DEFINITIONS *******************************'
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'WARNING: RESPECT THIS FILE FORMAT.THE SPACE BETWEEN DATA'
    WRITE(2,'(A)')'AND THE NUMBER FORMAT'
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'ARE THERE REINFORCEMENTS INTO THE DOMAIN (Y FOR YES AND N FOR NOT)'
    WRITE(2,'(A)')'Y'
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'METHOD UTILIZED TO MODEL THE REINFORCEMENTS (BEM1D OR FEM)'
    WRITE(2,'(A)')'1DBEM'
    WRITE(2,'(A)')''
    IF (ELASTOPLASTIC.EQ.1) THEN
        WRITE(2,'(A)')'ARE THERE REINFORCEMENTS ELASTOPLASTIC ? (Y OR N)'
        WRITE(2,'(A)')'Y'
    ELSE    
        WRITE(2,'(A)')'ARE THERE REINFORCEMENTS ELASTOPLASTIC ? (Y OR N)'
        WRITE(2,'(A)')'N'
    ENDIF
    WRITE(2,'(A)')''
    IF (BOND_SLIP.EQ.1) THEN
        WRITE(2,'(A)')'ARE THERE BOND-SLIP EFFECTS ? (Y OR N)'
        WRITE(2,'(A)')'Y'
    ELSE
        WRITE(2,'(A)')'ARE THERE BOND-SLIP EFFECTS? (Y OR N)'
        WRITE(2,'(A)')'N'
    ENDIF
    
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'NUMBER OF FE NODES'
    WRITE(2,'(I0)')N_NODES
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'NUMBER OF FE ELEMENTS'
    WRITE(2,'(I0)')N_ELEMS
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'NUMBER OF DIFFERENT MATERIALS'
    WRITE(2,'(I0)')N_MATERIALS
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'NODAL COORDINATES'
    WRITE(2,'(A)')'NODE X   Y   Z'
    DO I=1,N_NODES
        WRITE(2,'(I0,X,F11.6,A,F11.6,A,F11.6,A)')I,COORD_NODES(I,1),'0D0',COORD_NODES(I,2),'0D0',COORD_NODES(I,3),'0D0'
    ENDDO
    WRITE(2,'(A)')''  
    WRITE(2,'(A)')'INFORM THE POLYNOMIAL APPROXIMATION FOR ALL FE AND'
    WRITE(2,'(A)')'THE DOMAIN THAT THEY BELONG'
    WRITE(2,'(A)')'ELEMENT		DEGREE	        DOMAIN      MATERIAL'
    DO I=1,N_ELEMS
        WRITE(2,'(I0,X,I3,I3,I3)')I,ORD,FE_DOMAIN(I),FE_MATERIAL(I)    
    ENDDO
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'INFORM THE CONECTIVITY FOR ALL FE'
    DO I=1,N_ELEMS
        WRITE(2,'(A)')'ELEMENT  NODE'
        DO J=1,ORD+1
            WRITE(2,'(I0,X,I0,X)')I,CONECTI(I,J)
        ENDDO
    ENDDO
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'MATERIAL FE PROPERTIES'
    WRITE(2,'(A)')'MATERIAL         YOUNG        CROSS SECTION AREA      VISCOELASTICITY      MOD.BOLTZMANN     LEI CONST(1-KELVIN,2-BOLTZMANN,3-MAXWELL,4-HOOKE)'
    DO J=1,N_MATERIALS
        WRITE(2,'(I0,X,F18.6,A,F18.6,A,F18.6,A,F18.6,A,I8)')J,MAT_PROP(J,2),'0D0',MAT_PROP(J,1),'0D0',MAT_PROP(J,3),'0D0',MAT_PROP(J,4),'0D0',MODEL(J)
    ENDDO
    WRITE(2,'(A)')''
    
    IF (ELASTOPLASTIC.EQ.1) THEN
        WRITE(2,'(A)')'IF ELASTOPLASITC REINFORCEMENTS, GIVE PLASTIC PROPERTIES'
        WRITE(2,'(A)')'MAT          SIGMA Y              KP'
        DO J=1,N_MATERIALS
            WRITE(2,'(I0,F18.6,A,F18.6,A)')J,MAT_PROP_EP(J,1),'0D0',MAT_PROP_EP(J,2),'0D0'
        ENDDO
        WRITE(2,'(A)')''
        WRITE(2,'(A)')'NUMBER OF LOAD STEPS'
        WRITE(2,'(I4)')LOAD_STEPS
        WRITE(2,'(A)')''
        WRITE(2,'(A)')'TOLERANCE FOR ELASTOPLASTIC ANALYSIS'
        WRITE(2,'(F9.6,A)')TOL,'0D0'
    ENDIF
    
    IF (BOND_SLIP.EQ.1) THEN
        WRITE(2,'(A)')'IF CONSIDER BOND-SLIP EFFECTS, GIVE PROPETIES'
        WRITE(2,'(A)')'INFORM SLIP LAW?'
        WRITE(2,'(A)')'CONST'
        WRITE(2,'(A)')''
        WRITE(2,'(A)')'INFOR SLIP LAW PARAMETERS'
        WRITE(2,'(A)')'PARAM    VALUE'
        WRITE(2,'(I4,F25.2,A)')1,10000000000.0,'0D0'
        WRITE(2,'(A)')''
        WRITE(2,'(A)')'NUMBER OF LOAD STEPS'
        WRITE(2,'(I4)')1
        WRITE(2,'(A)')''
        WRITE(2,'(A)')'TOLERANCE FOR BOND-SLIP ANALYSIS'
        WRITE(2,'(A)')' 0.0010D0'
    ENDIF
    CLOSE(2)
    
    !SAIDA ALTERNATIVA PARA ACADVIEW COM PONTOS
    OPEN(3,FILE='Output_data/FIBERS_VIEW.OGL',STATUS='UNKNOWN')
    WRITE(3,'(A)')' Arquivo de pós-processamento'
    WRITE(3,'(A)')''
    WRITE(3,'(A)')'N. de nos    N. de elementos     N. de listas'
    WRITE(3,'(A)')'#'
    WRITE(3,'(3(I0,X))')N_NODES,N_ELEMS,0
    WRITE(3,'(A)')''
    WRITE(3,'(A)')'COORDX   COORDY   COORDZ  DX   DY   DZ'
    WRITE(3,'(A)')'#'
    DO I=1,N_NODES
        WRITE(3,'(6E18.6)')COORD_NODES(I,1),COORD_NODES(I,2),COORD_NODES(I,3),0.0D0,0.0D0,0.0D0
    ENDDO 
    WRITE(3,'(A)')' tpelem (1 - barra / 2 - triang / 3 - quad) grauaprox nó1 nó2...non'
    WRITE(3,'(A)')'#'
    DO I=1,N_ELEMS
        WRITE(3,'(I9,X,I0,X)',ADVANCE='NO')1,ORD
        DO J=1,ORD+1
            WRITE(3,'(I0,X)',ADVANCE='NO')CONECTI(I,J)
        ENDDO
        WRITE(3,'(I2)')1
        ENDDO
    CLOSE(3)
    
END SUBROUTINE