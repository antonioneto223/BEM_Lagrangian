FUNCTION F_MODELO(S,P_ATUANTE)
    
    !
    !   FUNCTION THAT RETURNS F_MODELO FOR A GIVEN NODE THAT CONTAINS A CURRENT SLIP VALUE EQUALS TO S
    !  
    !
    
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    REAL*8,INTENT(IN)::S,P_ATUANTE
    REAL*8::F_MODELO
    
    REAL*8::F_MAX,F_RES,S_1,S_2,S_3,S_4,ANG1,F1
    
    SELECT CASE(SLIP_LAW)
    CASE('CONST')
        
        !F_MODELO=SLIP_PARAMETERS(1)
        F_MODELO=DABS(P_ATUANTE)
        IF (F_MODELO.GT.SLIP_PARAMETERS(1)) F_MODELO=SLIP_PARAMETERS(1)
    
    CASE('LINEAR')
        
        F_MAX=SLIP_PARAMETERS(1)
        F_RES=SLIP_PARAMETERS(2)
        S_1=SLIP_PARAMETERS(3)
        S_2=SLIP_PARAMETERS(4)
        
        IF (DABS(S).LE.TOLER) THEN          ! TREHCO VERTICAL
        
            !F_MODELO=F_MAX
            F_MODELO=DABS(P_ATUANTE)
            IF (F_MODELO.GT.F_MAX) F_MODELO=F_MAX
            
        ELSE IF (DABS(S).LE.S_1) THEN            ! PRIMEIRO TRECHO
            
            F_MODELO=F_MAX
            
        ELSE IF (DABS(S).LT.S_2) THEN       ! SEGUNDO TRECHO
            
            F_MODELO=(F_MAX-F_RES)/(S_1-S_2)*DABS(S)
            F_MODELO=F_MODELO+(F_RES*S_1-F_MAX*S_2)/(S_1-S_2)
            
        ELSE                                ! TERCEIRO TRECHO
            
            F_MODELO=F_RES
            
        ENDIF
        
    CASE('HUANG')
        
        F_MAX=0.450D0*SLIP_PARAMETERS(1)
        F_RES=0.40D0*F_MAX
        S_1=SLIP_PARAMETERS(3)
        S_2=SLIP_PARAMETERS(4)
        S_3=SLIP_PARAMETERS(5)
        S_4=3.0D0*S_3
        ANG1=SLIP_PARAMETERS(2)
                
        IF (DABS(S).LT.S_1) THEN            ! PRIMEIRO TRECHO
            
            F_MODELO=F_MAX*(DABS(S)/S_1)**ANG1
            
        ELSE IF (DABS(S).LE.S_2) THEN       ! SEGUNDO TRECHO
            
            F_MODELO=F_MAX
            
        ELSE IF (DABS(S).LE.S_3) THEN       ! TERCEIRO TRECHO
            
            F_MODELO=F_MAX-(F_MAX-F_RES)*(DABS(S)-S_2)/(S_3-S_2)
            
        ELSE IF (DABS(S).LT.S_4) THEN       ! QUARTO TRECHO
            
            F_MODELO=F_RES-F_RES*(DABS(S)-S_3)/(S_4-S_3)
            
        ELSE                                ! QUINTO TRECHO
            
            F_MODELO=0.0D0
            
        ENDIF
            
    CASE('HUANG2')
        
        F_MAX=0.4502D0*SLIP_PARAMETERS(1)
        F_RES=0.3997D0*F_MAX
        S_1=SLIP_PARAMETERS(3)
        S_2=SLIP_PARAMETERS(4)
        S_3=SLIP_PARAMETERS(5)
        S_4=3.0D0*S_3
        ANG1=SLIP_PARAMETERS(2)
        
        F1=(F_MAX-ANG1*F_MAX)/(1.0D0+ANG1)
        
        IF (DABS(S).LE.TOLER) THEN          ! TRECHO ZERO
        
            F_MODELO=DABS(P_ATUANTE)
            IF (F_MODELO.GT.F1) F_MODELO=F1
            
        ELSE IF (DABS(S).LT.S_1) THEN       ! PRIMEIRO TRECHO
            
            F_MODELO=F1+((F_MAX-F1)/S_1)*DABS(S)
            
        ELSE IF (DABS(S).LE.S_2) THEN       ! SEGUNDO TRECHO
            
            F_MODELO=F_MAX
            
        ELSE IF (DABS(S).LE.S_3) THEN       ! TERCEIRO TRECHO
            
            F_MODELO=F_MAX-(F_MAX-F_RES)*(DABS(S)-S_2)/(S_3-S_2)
            
        ELSE IF (DABS(S).LT.S_4) THEN       ! QUARTO TRECHO
            
            F_MODELO=F_RES-F_RES*(DABS(S)-S_3)/(S_4-S_3)
            !F_MODELO=F_RES
            
        ELSE                                ! QUINTO TRECHO             
            
            F_MODELO=0.0d0
            
        ENDIF
        
    ENDSELECT
    
    
END FUNCTION
   
    
    
! **********************************************************************************************************
! **********************************************************************************************************
! **********************************************************************************************************
! **********************************************************************************************************
    
    
    
SUBROUTINE READ_SLIP_PARAMETERS

    ! 
    ! SUBOROUTINE TO READ MODEL'S PARAMETERS
    !
    
    USE REINFORCEMENTS
    
    IMPLICIT NONE
        
    INTEGER::K,IOS
    
    READ(1,*)!
    READ(1,*)! IF CONSIDERING BOND-SLIP EFFECTS, GIVE PRPERTIES
    READ(1,*)!INFORM COESIVE LAW FOR BON-SLIP? (CONST OR ...)
    READ(1,*)SLIP_LAW
    READ(1,*)
    READ(1,*)! INFORM COESVI LAW PARAMETERS
    READ(1,*)!PARAM     VALUE
    
    SELECT CASE (SLIP_LAW)
    CASE ('CONST')
        ALLOCATE(SLIP_PARAMETERS(1))
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)             ! F_MAX DE ADERENCIA
    CASE ('LINEAR')
        ALLOCATE(SLIP_PARAMETERS(4))
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)           ! F_MAX DE ADERENCIA 
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)           ! F_RES DE ADERENCIA
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)           ! S_1 DE DESLOCAMENTO RELATIVO
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)           ! S_2 DE DESLOCAMENTO RELATIVO
    CASE ('HUANG')
        ALLOCATE(SLIP_PARAMETERS(5))
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)        ! F_C DO CONCRETO
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)        ! ALFA DA CURVA
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)        ! S_1 DE DESLOCAMENTO RELATIVO
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)        ! S_2 DE DESLOCAMENTO RELATIVO
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)        ! DISTANCIA ENTRE NERVURAS DA FIBRA
    CASE ('HUANG2')
        ALLOCATE(SLIP_PARAMETERS(5))
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)        ! F_C DO CONCRETO
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)        ! ALFA DA CURVA
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)        ! S_1 DE DESLOCAMENTO RELATIVO
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)        ! S_2 DE DESLOCAMENTO RELATIVO
        READ(1,*,IOSTAT=IOS)K,SLIP_PARAMETERS(K)        ! DISTANCIA ENTRE NERVURAS DA FIBRA
    END SELECT
    
    IF (IOS.NE.0) THEN
        WRITE(*,*)'**************** ERROR IN REINFORCEMENTS FILE READING *****************'
        WRITE(*,*)'       THE NUMBER OF PARAMETERS GIVEN FOR SLIP_LAW IS INCORRECT'
        WRITE(*,*)'            PLEASE VERIFY TXT FILE AND RESTART PROGRAM'
        WRITE(*,*)'***********************************************************************'
        READ(*,*)
    ENDIF

    READ(1,*)!
	READ(1,*)!NUMBER OF LOAD STEPS
	READ(1,*)N_PASSOS
    READ(1,*)!
    READ(1,*)!TOLERANCE FOR BOND-SLIP ANALYSIS
    READ(1,*)SLIP_TOL
    READ(1,*)!
    READ(1,*)!NUMBER OF NODAL DISPLACEMENTS APPLIED
    READ(1,*)N_NODAL_FORCES
    
END SUBROUTINE