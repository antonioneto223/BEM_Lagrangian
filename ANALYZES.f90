!Esta subrotina é responsável por direcionar o fluxo do programa.
!Esta versão da subrotina possibilita escolher entre análises 
!elástica, fratura frágii ou fadiga.

SUBROUTINE ANALYZES
!
USE ISOPARAMETRIC_MESH
USE MATERIALS
USE CRACKS_DUAL_BEM
USE SUB_REGIONS_INTERFACES
USE ANALYSIS
USE REMESHING
USE PROPAGATION
USE FATIGUE
USE XBEM_FORCE_VARIABLES
USE XBEM_SUPPORT_VARIABLES
USE COHESIVE
USE REINFORCEMENTS
USE VISCO
!
IMPLICIT NONE
!
INTEGER::I,J
!
REAL*8::U_INITIAL[ALLOCATABLE](:),T_INITIAL[ALLOCATABLE](:),AUX[ALLOCATABLE](:),P_RES(3),CMOD_AB  
!
IF (VISCO_ANALYSIS) THEN
    
    CALL REINF_INPUT(FILE_REINFORCEMENTS) 
    
    WRITE(*,*)!
    WRITE(*,*)'COLLOCATION OUTPUT'
    CALL COLLOCATION_OUTPUT(-1)
        
    WRITE(*,*)!
    WRITE(*,*)'NORMAL OUTWARD'
    CALL NORMAL_OUTWARD_OUTPUT
    
    WRITE(*,*)!
    WRITE(*,*)'PROCESSING VISCOELASTIC ANALYSIS ...'
    
    WRITE(*,*)'HG_MATRICES'
    CALL HG_MATRICES_VISCO
    
    WRITE(*,*)'REINFORCEMENTS_COUPLING'
    CALL REINFORCEMENTS_COUPLING
    
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        WRITE(*,*)!
        WRITE(*,*)'     PROCESSING VISCOELASTIC ANALYSIS WITH REINFORCEMENTS ...'
        !IF (REINFORCEMENTS_ELASTOPLASTIC .EQ. 'Y')  THEN
        !    WRITE(*,*)!
        !    WRITE(*,*)'              PROCESSING ELASTOPLASTIC ANALYSIS FOR REINFORCEMENTS ...'
        !ENDIF
    ENDIF
    
    WRITE(*,*)'EDGE_CRACKS_REMESHING'
    CALL EDGE_CRACKS_REMESHING ! Atualiza HG com remalhamento e insere termos livres

    J=1 ! Marca a fase de carregamento
    
    CALL TXTS_VISCO
    WRITE(*,*)!
    WRITE(*,*)'SOLVING_ELASTIC_VISCO_PROBLEM'
    
    DO I=1,NP_TIME ! INICIATE TIME STEPS ___________________________________________________________
        
        IF ((MOD(I,5).EQ.0).OR.(I.EQ.1))  WRITE(*,*)'TIMPE STEP:',I
        
        ! Verifica se é necessario mudar fase de carregamento
        IF ((I .EQ. FINAL_STEP_NUMBER(J)) .AND. (I .NE. NP_TIME)) THEN
            CALL UPDATE_LOAD_FASE
            J= J + 1    
        ENDIF
                
        ! Atualizando variaveis
        CALL UPDATE_VISCO_VARIABLES(I)
        
        ! Resolvendo sistema no passo de tempo atual
        CALL SOLVE_ELASTIC_CRACK_PROBLEM   
        IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') CALL REINFORCEMENTS_STRESS
        
        !IF (REINFORCEMENTS_ELASTOPLASTIC .EQ. 'Y') THEN
        !    CALL EP_PROCESS_VISC(I)
        !ELSE
        !    CALL TROCA_COLUNAS_VISCO
        !ENDIF    
        
        ! Internal points calculation
        CALL INTERNAL_POINTS_VISCO(I)
        
        ! Rotinas de saida de dados do visco
        CALL BOUND_POINT_VISCO_OUTPUT(I)
        CALL INT_POINT_VISCO_OUTPUT(I) 
        CALL ACADVIEW_OUTPUT_VISCO(I)
        
    ENDDO ! END OF ALL TIME STEPS _________________________________________________________________
    
    ! USUAL OUTPUT DATAS
    WRITE(*,*)'COLLOCATION_OUTPUT'
    CALL COLLOCATION_OUTPUT(0)
    
    WRITE(*,*)'INTERNAL POINTS'
    !CALL INTERNAL_POINTS
    
    WRITE(*,*)'ACADVIEW_OUTPUT'
    CALL ACADVIEW_OUTPUT
    
    CALL NORMAL_OUTWARD_OUTPUT
    
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        CALL REINFORCEMENTS_OUTPUT
        DO I=1,4   ! CLOSING TXT VISCO-REINFORCEMENTS
            CLOSE(3*COUNT_TXT+I)
        ENDDO
    ENDIF
    
    ! CLOSING TXT'S VISCO
    DO I=1,N_NODES_BOUND_VISCO
        CLOSE(I+COUNT_TXT)
    ENDDO
    DO I=1,N_NODES_INT_VISCO
        CLOSE(I+2*COUNT_TXT)
    ENDDO
        
ELSE
    SELECT CASE(TRIM(ANALYSIS_TYPE))
    CASE("ELASTIC")
    !------ELASTIC ANALYSIS-------------
        WRITE(*,*)!
        WRITE(*,*)'PROCESSING ELASTIC ANALYSIS'
    !
        CALL COLLOCATION_OUTPUT(-1)
    !   
        WRITE(*,*)'NORMAL_OUTWARD_OUTPUT'
        CALL NORMAL_OUTWARD_OUTPUT
    !    
        IF((NUM_DIST_LOADS+NUM_DIST_SUPPS).NE.0) CALL DIST_LOADS_SUPPORTS_OUTPUT
        !IF(NUM_CON_SUPPS.NE.0) CALL FIND_ENRICHED_SUPPS
    !
        !
        !IF((NUM_CON_SUPPS.NE.0).OR.(NUM_DIST_SUPPS.NE.0))THEN
        !    WRITE(*,*)'XBEM_SUPPORT'
        !    CALL FIND_ENRICHED_SUPPS
        !    CALL XBEM_SUPPORT
        !ENDIF
    !!    
        WRITE(*,*)'HG_MATRICES'
        CALL HG_MATRICES 
    !
        IF (EXISTANCE_REINFORCEMENTS.NE.'N') THEN
            WRITE(*,*)'REINFORCEMENTS_COUPLING'
            CALL REINFORCEMENTS_COUPLING
        ENDIF
    !
        WRITE(*,*)'EDGE_CRACKS_REMESHING'
        CALL EDGE_CRACKS_REMESHING ! avalia distâncias dos nós das faces da fissura com os nós do contorno externo e faz remalhamento - remalha e atualiza
    !   
        WRITE(*,*)'NORMAL_OUTWARD_OUTPUT'
        CALL NORMAL_OUTWARD_OUTPUT
    !
        IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
            IF (EP_ANALYSIS_REINF) THEN
                WRITE(*,*)'CALCULING ELASTOPLASTIC NON LINEAR PROCESS'
                CALL REINFORCEMENTS_EP_PROCESS
                !CALL REINFORCEMENTS_EP_PROCESS_2 !TESTE
                WRITE(*,*)'WRITING REINFORCEMENTS OUTPUT'
                CALL REINFORCEMENTS_OUTPUT
            ELSE IF (SLIP_REINF) THEN
                WRITE(*,*)'CALCULING BOND-SLIP NON LINEAR PROCESS'
                CALL REINFORCEMENTS_SLIP_PROCESS
                WRITE(*,*)'WRITING REINFORCEMENTS OUTPUT'
                CALL REINFORCEMENTS_OUTPUT
            ELSE
                WRITE(*,*)'SOLVE_ELASTIC_REINFORCED_PROBLEM'
                CALL REINFORCEMENTS_SOLVER
                CALL REINFORCEMENTS_STRESS
                WRITE(*,*)'WRITING REINFORCEMENTS OUTPUT'
                CALL REINFORCEMENTS_OUTPUT
            ENDIF
        ELSE
        !
            IF((NUM_CON_SUPPS.NE.0).OR.(NUM_DIST_SUPPS.NE.0))THEN
                WRITE(*,*)'XBEM_SUPPORT'
                CALL FIND_ENRICHED_SUPPS
                CALL XBEM_SUPPORT
            ENDIF
        !  
            WRITE(*,*)'SOLVE_ELASTIC_CRACK_PROBLEM'
        !
            CALL SOLVE_ELASTIC_CRACK_PROBLEM    
        ENDIF
    !
        WRITE(*,*)'OUTPUT'
        CALL COLLOCATION_OUTPUT(0)   
    !    
        !WRITE(*,*)'INTERNAL POINTS'
        !CALL INTERNAL_POINTS
        !WRITE(*,*)'PASSOU PELO INTERNAL POINTS - TECLE ENTER'
        !READ(*,*)
    !
        CALL ACADVIEW_OUTPUT
    !    
        !WRITE(*,*) 'RESULTING P IN PATH WITH DISPLACEMENT PRESCRIBED'
        !CALL RESULTING_FORCE_PRESCRIBED_DISPL(F_RES)
    !
        WRITE(*,*)'SIT_EXTRACTION'
        SELECT CASE(SITEXTRACTION_METHOD)
        CASE("DET")
            CALL DISPLACEMENT_EXTRAPOLATION_TECHNIQUE(0)
        CASE("DFT")    
            CALL DISPLACEMENT_FITTING_TECHNIQUE(0) 
        ENDSELECT 
    !    CALL CONVERGENCE      
    !        
    CASE("BRITTLE_FRACTURE")
    !------BRITTLE FRACTURE ANALYSIS-------------
        WRITE(*,*)!
        WRITE(*,*)'PROCESSING BRITTLE FRACTURE ANALYSIS'
    !   INITIAL MESH
        CALL COLLOCATION_OUTPUT(-10)
    !
        WRITE(*,*)'HG_MATRICES'
        CALL HG_MATRICES
    !
        WRITE(*,*)'EDGE_CRACKS_REMESHING'
        CALL EDGE_CRACKS_REMESHING     
    !   MESH AFTER EDGE_CRACKS_REMESHING     
        CALL COLLOCATION_OUTPUT(-1)
    !
    !   XBEM SUPP
        IF((NUM_CON_SUPPS.NE.0).OR.(NUM_DIST_SUPPS.NE.0))THEN
            WRITE(*,*)'XBEM_SUPPORT'
            CALL FIND_ENRICHED_SUPPS
            CALL XBEM_SUPPORT
        ENDIF
    !
        ALLOCATE(U_INITIAL(3*(N_COLLOCPOINTS)),T_INITIAL(3*(N_COLLOCPOINTS)))      
        U_INITIAL=U
        T_INITIAL=T
    
        N_CRACK_POINTS=0
        WRITE(*,*)'SOLVE_BVP'
        CALL SOLVE_BVP 
       
        SELECT CASE(SITEXTRACTION_METHOD)
        CASE("DET")
            CALL DISPLACEMENT_EXTRAPOLATION_TECHNIQUE(0)
        CASE("DFT")    
            CALL DISPLACEMENT_FITTING_TECHNIQUE(0) 
        ENDSELECT
        ! CALIBRAÇÃO TORÇÃO
        !OPEN(9,file='Output_data\EX4_CALIB.txt',status='unknown')
        !CALL CMOD_EX4(CMOD_AB)
        !WRITE(9,*)CMOD_AB,XBEM_SUPP_REACTIONS(1)  
        !CLOSE(9)
    !
        CALL PROPAGATION_OUTPUT(0)
        CALL COLLOCATION_OUTPUT(0)
        CALL CRACK_OUTPUT(0)       
    
        DO I=1,N_STEPS     
            WRITE(*,*)'CRACK_PROPAGATION',I
            CALL CRACK_PROPAGATION   
        
            CALL EDGE_CRACKS_REMESHING 
                
            U=0.D0
            T=0.D0
            DO J=1,3*(N_COLLOCPOINTS-N_CRACKTIPCOLLOCPOINTS-N_ADDCOLLOCPOINTS)
                U(J)=U_INITIAL(J)
                T(J)=T_INITIAL(J)
            ENDDO
            DEALLOCATE(U_INITIAL,T_INITIAL)
            ALLOCATE(U_INITIAL(3*(N_COLLOCPOINTS)),T_INITIAL(3*(N_COLLOCPOINTS)))      
            U_INITIAL=U
            T_INITIAL=T

            N_CRACK_POINTS=0
            IF((NUM_CON_SUPPS.NE.0).OR.(NUM_DIST_SUPPS.NE.0))THEN
                WRITE(*,*)'XBEM_SUPPORT'
                CALL XBEM_SUPPORT
            ENDIF
            CALL SOLVE_BVP       
        
            SELECT CASE(SITEXTRACTION_METHOD)
            CASE("DET")
                CALL DISPLACEMENT_EXTRAPOLATION_TECHNIQUE(I)
            CASE("DFT")    
                CALL DISPLACEMENT_FITTING_TECHNIQUE(I) 
            ENDSELECT 
        
            CALL PROPAGATION_OUTPUT(I)
            CALL COLLOCATION_OUTPUT(I) 
            CALL CRACK_OUTPUT(I)                                        
        ENDDO    
    !          
    CASE("FATIGUE")
    !------FATIGUE ANALYSIS-------------
        WRITE(*,*)!
        WRITE(*,*)'PROCESSING FATIGUE ANALYSIS'
    !    
        CALL HG_MATRICES
    
        CALL EDGE_CRACKS_REMESHING     
    !    
        ALLOCATE(U_INITIAL(3*(N_COLLOCPOINTS)),T_INITIAL(3*(N_COLLOCPOINTS)))      
        U_INITIAL=U
        T_INITIAL=T
    ! 
        MINIMUM_MAXIMUM_LOAD="MINIMUM"   
        U=LOAD_RATIO*U_INITIAL
        T=LOAD_RATIO*T_INITIAL 
    
        N_CRACK_POINTS=0
        CALL SOLVE_BVP 
    !    CALL SOLVE_ELASTIC_CRACK_PROBLEM
    
        SELECT CASE(SITEXTRACTION_METHOD)
        CASE("DET")
            CALL DISPLACEMENT_EXTRAPOLATION_TECHNIQUE(0)
        CASE("DFT")    
            CALL DISPLACEMENT_FITTING_TECHNIQUE(0) 
        ENDSELECT
    
        CALL PROPAGATION_OUTPUT(0)    
    
        CALL MAXIMUM_KEQ
    
        DKi=Keqmax 

        MINIMUM_MAXIMUM_LOAD="MAXIMUM"      
        U=U_INITIAL
        T=T_INITIAL 

        N_CRACK_POINTS=0
        CALL SOLVE_BVP       
    !    CALL SOLVE_ELASTIC_CRACK_PROBLEM 
    
        SELECT CASE(SITEXTRACTION_METHOD)
        CASE("DET")
            CALL DISPLACEMENT_EXTRAPOLATION_TECHNIQUE(0)
        CASE("DFT")    
            CALL DISPLACEMENT_FITTING_TECHNIQUE(0) 
        ENDSELECT
    
        CALL PROPAGATION_OUTPUT(0)
    
        CALL MAXIMUM_KEQ
    
        DKi=Keqmax-DKi         
    !
        IF(DKi.GE.FFMP(1,4))THEN
            NCICLESi=0.D0
            INITIAL_CRACKSIZE=15.D0
            NCICLES=NCICLESi
            CRACKSIZE=INITIAL_CRACKSIZE
        !  
            OPEN(23,FILE='FATIGUE.TXT',STATUS='UNKNOWN')
	        WRITE(23,*)'         N                a                DKi                DKi+1'	     
	        WRITE(23,100)NCICLES,CRACKSIZE,DKi,DKf     
            DO I=1,N_STEPS    
                CALL FATIGUE_CRACK_PROPAGATION   
            
                CALL EDGE_CRACKS_REMESHING 
                    
                U=0.D0
                T=0.D0
                DO J=1,N_COLLOCPOINTS-N_CRACKTIPCOLLOCPOINTS-N_ADDCOLLOCPOINTS
                    U(J)=U_INITIAL(J)
                    T(J)=T_INITIAL(J)
                ENDDO
                DEALLOCATE(U_INITIAL,T_INITIAL)
                ALLOCATE(U_INITIAL(3*(N_COLLOCPOINTS)),T_INITIAL(3*(N_COLLOCPOINTS)))      
                U_INITIAL=U
                T_INITIAL=T

                MINIMUM_MAXIMUM_LOAD="MINIMUM"               
                U=LOAD_RATIO*U_INITIAL
                T=LOAD_RATIO*T_INITIAL

                N_CRACK_POINTS=0
                CALL SOLVE_BVP
    !            CALL SOLVE_ELASTIC_CRACK_PROBLEM            
                   
                SELECT CASE(SITEXTRACTION_METHOD)
                CASE("DET")
                    CALL DISPLACEMENT_EXTRAPOLATION_TECHNIQUE(I)
                CASE("DFT")    
                    CALL DISPLACEMENT_FITTING_TECHNIQUE(I) 
                ENDSELECT
            
                CALL PROPAGATION_OUTPUT(I) 
              
                CALL MAXIMUM_KEQ
            
                DKf=Keqmax 
            
                MINIMUM_MAXIMUM_LOAD="MAXIMUM"  
                U=U_INITIAL
                T=T_INITIAL

                N_CRACK_POINTS=0
                CALL SOLVE_BVP
    !            CALL SOLVE_ELASTIC_CRACK_PROBLEM
                              
                SELECT CASE(SITEXTRACTION_METHOD)
                CASE("DET")
                    CALL DISPLACEMENT_EXTRAPOLATION_TECHNIQUE(I)
                CASE("DFT")    
                    CALL DISPLACEMENT_FITTING_TECHNIQUE(I) 
                ENDSELECT 
            
                CALL PROPAGATION_OUTPUT(I)
              
                CALL MAXIMUM_KEQ
            
                DKf=Keqmax-DKf
      
                DNCICLES=((MAX_CRACKINCREMENT*Keqmax/Keqmax)*((DKf**(1-FFMP(1,3)))-(DKi**(1-FFMP(1,3)))))/(FFMP(1,2)*(1-FFMP(1,3))*(DKf-DKi))
                CRACKSIZE=CRACKSIZE+MAX_CRACKINCREMENT*Keqmax/Keqmax
                NCICLES=NCICLES+DNCICLES 
             
                WRITE(23,100)NCICLES,CRACKSIZE,DKi,DKf  
                     
                DKi=DKf	                
            ENDDO
            CLOSE(23)   
        ENDIF    
    !  
    CASE("COHESIVE")
    !------COHESIVE ANALYSIS WITH SUB-REGIONS-------------
        WRITE(*,*)!
        WRITE(*,*)'PROCESSING COHESIVE ANALYSIS'
    !
        CALL COLLOCATION_OUTPUT(-1)
    !   
        WRITE(*,*)'NORMAL_OUTWARD_OUTPUT'
        CALL NORMAL_OUTWARD_OUTPUT
    !
        CALL INTERFACE_PRE_PROCESSOR
    !
        ALLOCATE(AUX(2*N_INTERFACE_POINTS))
        AUX=0.D0   
        CALL INTERFACE_OUTPUT(-1,AUX)
    !
        WRITE(*,*)'HG_MATRICES'
        CALL HG_MATRICES 
    !
        WRITE(*,*)'EDGE_CRACKS_REMESHING'
        CALL EDGE_CRACKS_REMESHING ! avalia distâncias dos nós das faces da fissura com os nós do contorno externo e faz remalhamento - remalha e atualiza
    !   
        IF(NUM_CON_SUPPS.GT.0)THEN  
            WRITE(*,*)'XBEM_SUPPORT'
            CALL FIND_ENRICHED_SUPPS    
            CALL XBEM_SUPPORT
        ENDIF
    ! 
        SELECT CASE(TRIM(COHESIVE_OPERATOR))
        CASE('CONSTANT')
            WRITE(*,*)'SOLVE_COHESIVE_CRACK_CONSTANT_OPERATOR'
            CALL SOLVE_COHESIVE_CRACK_CO
    !
        CASE('TANGENT')
            WRITE(*,*)'SOLVE_COHESIVE_CRACK_TANGENT_OPERATOR'
            CALL SOLVE_COHESIVE_CRACK_TO
        END SELECT
    !   
        CALL NORMAL_OUTWARD_OUTPUT
    !
        WRITE(*,*)'OUTPUT'
        CALL COLLOCATION_OUTPUT(0)
    !    
        CALL ACADVIEW_OUTPUT  
    !        
        END SELECT
    ENDIF
!
100	FORMAT(3x,E15.7,3x,E15.7,3x,E15.7,3x,E15.7)       
END SUBROUTINE ANALYZES
    





