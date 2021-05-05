SUBROUTINE CREATE_RANDOM(NODE,ELEM)
    
    USE PRE_REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER,INTENT(INOUT)::NODE,ELEM
    
    LOGICAL::DENTRO,EFETIVO
    INTEGER::I,J,K,N_TOTAL_FIBRAS,CONT
    REAL*8::X(3),Y(3),XI(2),FRACAO_ATUAL,PI,AUX,VOLUME_TOTAL,COMPRIMENTOS,ANGS(3),AUX1(4),ETA(3),L_EFETIVO
    
    PI=DACOS(-1.0D0)
    EFETIVO=.FALSE.
    AREA_F=PI*F_RAIO**2.0D0
    CALL RANDOM_SEED()
    
    ! CONSIDERAR COMPRIMENTO EFETIVO PARA RESISTENCIA DAS PONTAS?
    EFETIVO=.TRUE.
    IF (EFETIVO) THEN
        L_EFETIVO=F_COMP+F_RAIO
    ELSE
        L_EFETIVO=F_COMP
    ENDIF
    ! ----------------------------------------------------
    
    FRACAO_ATUAL=0.0D0
    IF (ESPACO.EQ.'CILINDRO-X') THEN
        VOLUME_TOTAL=PI*RAIO**2.0D0*(LX-LX0)
        ! TIRANDO DAS EXTREMIDADES
        LX0=LX0+F_RAIO*1.10D0
        LX=LX-F_RAIO*1.10D0
        RAIO=RAIO-F_RAIO*1.10D0
    ELSE IF (ESPACO.EQ.'CUBO') THEN
        VOLUME_TOTAL=(LX-LX0)*(LY-LY0)*(LZ-LZ0)
        ! TIRANDO DAS EXTREMIDADES
        LX0=LX0+F_RAIO*1.10D0
        LX=LX-F_RAIO*1.10D0
        LY0=LY0+F_RAIO*1.10D0
        LY=LY-F_RAIO*1.10D0
        LZ0=LZ0+F_RAIO*1.10D0
        LZ=LZ-F_RAIO*1.10D0
    ENDIF
    
    N_TOTAL_FIBRAS=0
    COMPRIMENTOS=0.0D0
    
    DO WHILE((FRACAO-FRACAO_ATUAL).GT.0.001) 
        
        IF (ESPACO.EQ.'CILINDRO-X') THEN
            ! SORTEANDO PRIMEIRO PONTO ---------------------------

            ! SORTENADO COORDENADA X
1           CONTINUE
            CALL RANDOM_NUMBER(X(1))
            X(1)=LX0+X(1)*(LX-LX0)
            ! SORTEANDO COORDENADAS Y, DENTRO DO CILINDRO NÉ
2           CONTINUE
            CALL RANDOM_NUMBER(XI)
            X(2)=-RAIO+2.0D0*RAIO*XI(1)
            X(3)=-RAIO+2.0D0*RAIO*XI(2)
            DENTRO=.TRUE.
            AUX=X(2)**2.0D0+X(3)**2.0D0
            IF (AUX.GT.RAIO**2.0D0) DENTRO=.FALSE.
            IF (.NOT.DENTRO) GO TO 2
        
            ! SORTEANDO SEGUNDO PONTO ----------------------------
            CONT=0
4           CONTINUE
            CONT=CONT+1
            CALL RANDOM_NUMBER(XI)
            !XI(1)=XI(1)*2.0d0*PI
            !XI(2)=XI(2)*PI  
            
            !XI(1)=0.0d0            ! Para gerar fibras alinhadas com X
            !XI(2)=PI/2.0D0
            !XI(1)=PI/2.0d0            ! Para gerar fibras transversais a X
            !XI(2)=PI/2.0D0
            
            IF (DIST_F) THEN
                CALL PSEUDORANDOM_NUMBER(XI,MEDIA_INC_F,DESV_PAD_INC_F,TIPO_DIST_F)
                DO I=1,2
                    IF ((XI(I).GT.(2.0D0*PI)).OR.(XI(I).LT.(-2.0D0*PI))) GO TO 4
                ENDDO
            ELSE
                XI(1)=XI(1)*2.0d0*PI
                XI(2)=XI(2)*PI
            ENDIF
                
            Y(1)=X(1)+L_EFETIVO*DSIN(XI(2))*DCOS(XI(1))
            Y(2)=X(2)+L_EFETIVO*DSIN(XI(2))*DSIN(XI(1))
            Y(3)=X(3)+L_EFETIVO*DCOS(XI(2))
            
            ! VERIFICANDO SE O PONTO CAIU DENTRO DO CILINDRO:
            DENTRO=.TRUE.
            IF ((Y(1).LT.LX0).OR.(Y(1).GT.LX)) DENTRO=.FALSE.
            AUX=Y(2)**2.0D0+Y(3)**2.0D0
            IF (AUX.GT.RAIO**2.0D0) DENTRO=.FALSE.
            IF (CONT.GT.50) GO TO 1
            IF (.NOT.DENTRO) GO TO 4
            
        ELSE IF (ESPACO.EQ.'CUBO') THEN
            
            ! SORTENADO PONTO INICIAL
5           CONTINUE
            CALL RANDOM_NUMBER(ETA)
            X(1)=LX0+ETA(1)*(LX-LX0)
            X(2)=LY0+ETA(2)*(LY-LY0)
            X(3)=LZ0+ETA(3)*(LZ-LZ0)
            
            ! SORTEANDO SEGUNDO PONTO ----------------------------
            CONT=0
6           CONTINUE
            CONT=CONT+1
            CALL RANDOM_NUMBER(XI)
            XI(1)=XI(1)*2.0d0*PI
            XI(2)=XI(2)*PI
            !XI(1)=0.0D0            ! Para gerar fibras alinhadas com X
            !XI(2)=PI/2
            
            Y(1)=X(1)+L_EFETIVO*DSIN(XI(2))*DCOS(XI(1))
            Y(2)=X(2)+L_EFETIVO*DSIN(XI(2))*DSIN(XI(1))
            Y(3)=X(3)+L_EFETIVO*DCOS(XI(2))
            ! VERIFICANDO SE O PONTO CAIU DENTRO DO CUBO:
            DENTRO=.TRUE.
            IF ((Y(1).LT.LX0).OR.(Y(1).GT.LX)) DENTRO=.FALSE.
            IF ((Y(2).LT.LY0).OR.(Y(2).GT.LY)) DENTRO=.FALSE.
            IF ((Y(3).LT.LZ0).OR.(Y(3).GT.LZ)) DENTRO=.FALSE.
            
            IF (CONT.GT.50) GO TO 5
            IF (.NOT.DENTRO) GO TO 6
            
        ENDIF
        
        ! CRIANDO FIBRA ENTRE ESTES PONTOS
        CALL CREATE_FIBER_RANDOM(X,Y)
        
        ! CALCULANDO FRACAO DE VOLUME FEITA ATE AGORA
        FRACAO_ATUAL=FRACAO_ATUAL+AREA_F*F_COMP/VOLUME_TOTAL*100.0D0
                
        ! TESTES: SALVANDO COMPRIMENTO E ANGULOS:
        N_TOTAL_FIBRAS=N_TOTAL_FIBRAS+1        
        AUX1(1)=Y(1)-X(1)
        AUX1(2)=Y(2)-X(2)
        AUX1(3)=Y(3)-X(3)
        AUX1(4)=DSQRT(AUX1(1)**2.0D0+AUX1(2)**2.0D0+AUX1(3)**2.0D0)
        COMPRIMENTOS=COMPRIMENTOS+AUX1(4)
        AUX1(1)=AUX1(1)/AUX1(4)
        AUX1(2)=AUX1(2)/AUX1(4)
        AUX1(3)=AUX1(3)/AUX1(4)
        ANGS(1)=ANGS(1)+AUX1(1)
        ANGS(2)=ANGS(2)+AUX1(2)
        ANGS(3)=ANGS(3)+AUX1(3)
               
    ENDDO   
         
    COMPRIMENTOS=COMPRIMENTOS/DFLOAT(N_TOTAL_FIBRAS)
    ANGS(:)=ANGS(:)/DFLOAT(N_TOTAL_FIBRAS)
    !WRITE(*,*)''
    !WRITE(*,*)'MEDIA DO COMPRIMENTO:',COMPRIMENTOS
    !WRITE(*,*)''
    !WRITE(*,*)'MEDIA DE ALFA 1:',ANGS(1)
    !WRITE(*,*)'MEDIA DE ALFA 2:',ANGS(2)
    !WRITE(*,*)'MEDIA DE ALFA 3:',ANGS(3)
    !WRITE(*,*)''
    !WRITE(*,*)'NUMERO DE FIBRAS:',N_TOTAL_FIBRAS
    !WRITE(*,*)''
    !WRITE(*,*)''
    
    WRITE(*,*)'Total number of fibers:',N_TOTAL_FIBRAS 
    WRITE(*,*)'FIBERS FRACTION VOLUME OBTAINED (%):',FRACAO_ATUAL
    WRITE(*,*)''
        
END SUBROUTINE
    
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! ***************************************************************************************************************************************************** 
 
SUBROUTINE CREATE_FIBER_RANDOM(INITIAL_POINT,FINAL_POINT)

    USE PRE_REINFORCEMENTS    
    
    IMPLICIT NONE
    
    REAL*8,DIMENSION(3),INTENT(IN)::INITIAL_POINT,FINAL_POINT
    
    INTEGER::I,J,K,II,KK,JJ,GRAVOU,CONECTI_AUX(NUM_EL_F,ORD_F+1),NODE_L,ELEM_L,MAT_AUX_INT[ALLOCATABLE](:,:),SIZES(4)
    REAL*8::DELTA(3),XI(3),EL_DELTA(3),COORD_AUX(3),AUX(4),COORD_NODES_AUX(NUM_EL_F*ORD_F+1,3),MAT_AUX_REAL[ALLOCATABLE](:,:)
    
    GRAVOU=0
    
    DO J=1,3
        DELTA(J)=FINAL_POINT(J)-INITIAL_POINT(J)
        EL_DELTA(J)=DELTA(J)/(NUM_EL_F)
    ENDDO
    
    NODE_L=0
    ELEM_L=0
    
    DO I=1,NUM_EL_F
        
        ELEM_L=ELEM_L+1
        DO K=1,3
            XI(K)=INITIAL_POINT(K)+(I-1)*EL_DELTA(K)
        ENDDO
        
        DO J=1,ORD_F+1
        
            NODE_L=NODE_L+1
            IF (((J.EQ.1).AND.(I.NE.1)).OR.(GRAVOU.EQ.1)) NODE_L=NODE_L-1
            DO K=1,3
                COORD_AUX(K)=XI(K)+(J-1)*EL_DELTA(K)/(ORD_F)
            ENDDO
            
             ! VERIFICANDO SE JÁ TEM UM NÓ NESSA POSÇÃO
            AUX=0.0D0
            GRAVOU=0
            DO II=1,NODE_L-1
                DO K=1,3
                    AUX(K)=COORD_NODES_AUX(II,K)-COORD_AUX(K)
                ENDDO
                AUX(4)=DSQRT(AUX(1)**2.0D0+AUX(2)**2.0D0+AUX(3)**2.0D0)
                
                IF ((AUX(4).LT.TOLER).AND.(I.NE.1)) THEN
                    CONECTI_AUX(ELEM_L,J)=II
                    GRAVOU=1
                    EXIT
                ENDIF
            ENDDO
            
            IF (GRAVOU .EQ. 0) THEN
                DO K=1,3
                    COORD_NODES_AUX(NODE_L,K)=COORD_AUX(K)
                ENDDO
                CONECTI_AUX(ELEM_L,J)=NODE_L
            ENDIF
        ENDDO 
        
        !FE_MATERIAL(ELEM)=MATERIALS(FIBER)
        !FE_DOMAIN(ELEM)=FIBER_DOMAIN(FIBER)
        
    ENDDO
    
    ! AJUSTANDO DADOS NAS MATRIZES GLOBAIS
    IF (.NOT.ALLOCATED(COORD_NODES)) THEN
        ALLOCATE(COORD_NODES(NODE_L,3),CONECTI(ELEM_L,ORD_F+1))
        COORD_NODES=COORD_NODES_AUX
        CONECTI=CONECTI_AUX
    ELSE
        SIZES(1)=SIZE(CONECTI,1)
        SIZES(2)=SIZE(CONECTI,2)
        SIZES(3)=SIZE(COORD_NODES,1)
        SIZES(4)=SIZE(COORD_NODES,2)
        ALLOCATE(MAT_AUX_INT(SIZES(1),SIZES(2)),MAT_AUX_REAL(SIZES(3),SIZES(4)))
        MAT_AUX_INT=CONECTI
        MAT_AUX_REAL=COORD_NODES
        DEALLOCATE(CONECTI,COORD_NODES)
        ALLOCATE(CONECTI(SIZES(1)+ELEM_L,ORD_F+1),COORD_NODES(SIZES(3)+NODE_L,3))
        DO I=1,SIZES(1)
            CONECTI(I,:)=MAT_AUX_INT(I,:)
        ENDDO
        DO I=1,ELEM_L
            CONECTI(SIZES(1)+I,:)=CONECTI_AUX(I,:)+MAT_AUX_INT(SIZES(1),SIZES(2))
        ENDDO
        DO I=1,SIZES(3)
            COORD_NODES(I,:)=MAT_AUX_REAL(I,:)
        ENDDO
        DO I=1,NODE_L
            COORD_NODES(SIZES(3)+I,:)=COORD_NODES_AUX(I,:)
        ENDDO
    ENDIF
    

END SUBROUTINE
        
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! ***************************************************************************************************************************************************** 
 
SUBROUTINE CREATE_FILE_RANDOM(FILE_REINF)

    USE PRE_REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,J
    CHARACTER*30,INTENT(IN)::FILE_REINF
    
    ELASTOPLASTIC=0
    BOND_SLIP=0
    
    
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
    WRITE(2,'(I0)')SIZE(COORD_NODES,1)
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'NUMBER OF FE ELEMENTS'
    WRITE(2,'(I0)')SIZE(CONECTI,1)
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'NUMBER OF DIFFERENT MATERIALS'
    WRITE(2,'(I4)')1
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'NODAL COORDINATES'
    WRITE(2,'(A)')'NODE X   Y   Z'
    DO I=1,SIZE(COORD_NODES,1)
        WRITE(2,'(I0,D20.10,D20.10,D20.10)')I,COORD_NODES(I,1),COORD_NODES(I,2),COORD_NODES(I,3)
    ENDDO
    WRITE(2,'(A)')''  
    WRITE(2,'(A)')'INFORM THE POLYNOMIAL APPROXIMATION FOR ALL FE AND'
    WRITE(2,'(A)')'THE DOMAIN THAT THEY BELONG'
    WRITE(2,'(A)')'ELEMENT		DEGREE	        DOMAIN      MATERIAL'
    DO I=1,SIZE(CONECTI,1)
        WRITE(2,'(I0,X,I3,I3,I3)')I,ORD_F,1,1   
    ENDDO
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'INFORM THE CONECTIVITY FOR ALL FE'
    DO I=1,SIZE(CONECTI,1)
        WRITE(2,'(A)')'ELEMENT  NODE'
        DO J=1,ORD_F+1
            WRITE(2,'(I0,X,I0,X)')I,CONECTI(I,J)
        ENDDO
    ENDDO
    WRITE(2,'(A)')''
    WRITE(2,'(A)')'MATERIAL FE PROPERTIES'
    WRITE(2,'(A)')'MATERIAL         YOUNG        CROSS SECTION AREA      VISCOELASTICITY      MOD.BOLTZMANN     LEI CONST(1-KELVIN,2-BOLTZMANN,3-MAXWELL,4-HOOKE)'
    WRITE(2,'(I0,X,F18.6,A,F18.6,A,F18.6,A,F18.6,A,I8)')1,YOUNG_F,'0D0',AREA_F,'0D0',45.45454540D0,'0D0',10.0000D0,'0D0',4
    !DO J=1,N_MATERIALS
    !    WRITE(2,'(I4,F18.6,A,F18.6,A,F18.6,A,F18.6,A,I8)')J,MAT_PROP(J,2),'0D0',MAT_PROP(J,1),'0D0',MAT_PROP(J,3),'0D0',MAT_PROP(J,4),'0D0',MODEL(J)
    !ENDDO
    WRITE(2,'(A)')''
    
    IF (ELASTOPLASTIC.EQ.1) THEN
        WRITE(2,'(A)')'IF ELASTOPLASITC REINFORCEMENTS, GIVE PLASTIC PROPERTIES'
        WRITE(2,'(A)')'MAT          SIGMA Y              KP'
        DO J=1,N_MATERIALS
            WRITE(2,'(I0,X,F18.6,A,F18.6,A)')J,MAT_PROP_EP(J,1),'0D0',MAT_PROP_EP(J,2),'0D0'
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
    WRITE(3,'(A)')'Output_data/FIBERS_VIEW.TXT'
    WRITE(3,'(A)')''
    WRITE(3,'(A)')'N. de nos    N. de elementos     N. de listas'
    WRITE(3,'(A)')'#'
    WRITE(3,'(3(I0,X))')SIZE(COORD_NODES,1),SIZE(CONECTI,1),0
    WRITE(3,'(A)')''
    WRITE(3,'(A)')'COORDX   COORDY   COORDZ  DX   DY   DZ'
    WRITE(3,'(A)')'#'
    DO I=1,SIZE(COORD_NODES,1)
        WRITE(3,'(6E18.6)')COORD_NODES(I,1),COORD_NODES(I,2),COORD_NODES(I,3),0.0D0,0.0D0,0.0D0
    ENDDO
    WRITE(3,'(A)')' tpelem (1 - barra / 2 - triang / 3 - quad) grauaprox nó1 nó2...non'
    WRITE(3,'(A)')'#'
    DO I=1,SIZE(CONECTI,1)
        WRITE(3,'(I0,X,I0,X)',ADVANCE='NO')1,ORD_F
        DO J=1,ORD_F+1
            WRITE(3,'(I0,X)',ADVANCE='NO')CONECTI(I,J)
        ENDDO
        WRITE(3,'(I2)')1
        ENDDO
    CLOSE(3)

END SUBROUTINE
    
    
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! ***************************************************************************************************************************************************** 
    
    
SUBROUTINE INPUT_DATA_ALEAT
    
    USE PRE_REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,J,MAX_NUM_EL,IOS
    
    OPEN(1,FILE='Input_data\FIBERS.TXT',STATUS='OLD')
    
    READ(1,*)!RANDOM FIBER (A) OR DISCRETE FIBER (D)
    READ(1,*)!A
    READ(1,*)!
    READ(1,*)!ESPACO QUE VAI SORTEAR
    READ(1,*)ESPACO
    IF (ESPACO.EQ.'CILINDRO-X') THEN
        READ(1,*)!
        READ(1,*)! INICIO EM X
        READ(1,*)LX0
        READ(1,*)!
        READ(1,*)! FINAL EM X
        READ(1,*)LX
        READ(1,*)!
        READ(1,*)!RAIO EM YZ:
        READ(1,*)RAIO
    ELSE IF (ESPACO.EQ.'CUBO') THEN
        READ(1,*)!
        READ(1,*)! INICIO EM X
        READ(1,*)LX0
        READ(1,*)!
        READ(1,*)! FINAL EM X
        READ(1,*)LX
        READ(1,*)!
        READ(1,*)! INICIO EM Y
        READ(1,*)LY0
        READ(1,*)!
        READ(1,*)! FINAL EM Y
        READ(1,*)LY
        READ(1,*)!
        READ(1,*)! INICIO EM Z
        READ(1,*)LZ0
        READ(1,*)!
        READ(1,*)! FINAL EM Z
        READ(1,*)LZ
    ELSE
        WRITE(*,*)'ESCREVEU UM ESPACO QUE NAO EXISTE NO ARQUIVO FIBERS!'
        READ(*,*)
    ENDIF
    READ(1,*)!
    READ(1,*)! COMPRIMENTO DAS FIBRAS
    READ(1,*)F_COMP
    READ(1,*)!
    READ(1,*)! RAIO DAS FIBRAS
    READ(1,*)F_RAIO
    READ(1,*)!
    READ(1,*)! MODULO DE ELASTICIDAE DAS FIBRAS
    READ(1,*)YOUNG_F
    READ(1,*)!
    READ(1,*)! FRACAO DE VOLUME ( EM PORCENTAGEM )
    READ(1,*)FRACAO
    READ(1,*)!
    READ(1,*)! NUMERO DE ELEMENTOS POR FIBRA
    READ(1,*)NUM_EL_F
    READ(1,*)!
    READ(1,*)! ORDEM DOS ELEMENTOS
    READ(1,*)ORD_F
    
    IOS=10
    DIST_F=.FALSE.
    ! LENDO PARTE DA INCLINACAO DAS FIBRAS
    READ(1,*,IOSTAT=IOS)!
    READ(1,*,IOSTAT=IOS)!DISTRIBUICAO
    READ(1,*,IOSTAT=IOS)TIPO_DIST_F
    READ(1,*,IOSTAT=IOS)!
    READ(1,*,IOSTAT=IOS)!MEDIA  (COM X) E (COM Z)
    READ(1,*,IOSTAT=IOS)MEDIA_INC_F(1), MEDIA_INC_F(2)
    READ(1,*,IOSTAT=IOS)!
    READ(1,*,IOSTAT=IOS)!DESVIO PADRAO(COM X) E (COM Z)
    READ(1,*,IOSTAT=IOS)DESV_PAD_INC_F(1),DESV_PAD_INC_F(2)
    
    IF (IOS.EQ.0) DIST_F=.TRUE.
    
    
    CLOSE(1)
    
END SUBROUTINE
    
    
    
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! *****************************************************************************************************************************************************
! ***************************************************************************************************************************************************** 
  
    
    
    
SUBROUTINE PSEUDORANDOM_NUMBER(UI,MEDIA,DESVPAD,TIPO)
    !
    !   SUBROUTINE TO GENERATE THE PSEUDO-RANDOM NUMBER WITH DISTRIBUITON FOR FIBER INCLINATION
    !   UI: RECEIVES UNIFORM NUMBER 0 TO 1 AND RETURNS NUMBER WITH NORMAL DISTRIBUITION
    !

    IMPLICIT NONE
    
    CHARACTER*30,INTENT(IN)::TIPO
    REAL*8,INTENT(IN)::MEDIA(2),DESVPAD(2)
    REAL*8,INTENT(INOUT),DIMENSION(2)::UI
    REAL(8)::PI,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,Z1,Y,AUX1,AUX2,VAR,VARIAVEL(2)
    INTEGER::I,J

    ! ENCONTRANDO PONTOS AMOSTRADOS AO REDOR DO PONTO DE PROJETO
    P1=-0.3222324310880D0
    P2=-1.0000000000000D0
    P3=-0.34224220885470D0
    P4=-0.02042312102450D0
    P5=-0.00004536422101480D0
    Q1=0.09934846260600D0
    Q2=0.5885815704950D0
    Q3=0.5311034623660D0
    Q4=0.1035377528500D0
    Q5=0.003856070063400D0
                    
    DO J=1,2
        !CALL random_number(R1)
        !UI = R1
        SELECT CASE(TRIM(TIPO))
            CASE ('NORMAL')
                IF (UI(J).GE.0.0 .AND. UI(J).LE.0.5) THEN
                    Z1= DSQRT(DLOG(1/(UI(J)**2)))          
                    Y=-Z1-((P1+Z1*(P2+Z1*(P3+Z1*(P4+Z1*(P5)))))/(Q1+Z1*(Q2+Z1*(Q3+Z1*(Q4+Z1*(Q5))))))
                    VARIAVEL(j)=Y*DESVPAD(J)+MEDIA(J)
                END IF
                IF (UI(J).GE.0.5 .AND. UI(J).LE.1.0) THEN
                    Z1= DSQRT(DLOG(1/((1-UI(J))**2)))          
                    Y=Z1+((P1+Z1*(P2+Z1*(P3+Z1*(P4+Z1*(P5)))))/(Q1+Z1*(Q2+Z1*(Q3+Z1*(Q4+Z1*(Q5))))))
                    VARIAVEL(j)=Y*DESVPAD(J)+MEDIA(J)
                END IF 
            CASE ('LNORMAL')
                IF (UI(J).GE.0.0 .AND. UI(J).LE.0.5) THEN
                    Z1= DSQRT(DLOG(1/(UI(J)**2)))          
                    Y=-Z1-((P1+Z1*(P2+Z1*(P3+Z1*(P4+Z1*(P5)))))/(Q1+Z1*(Q2+Z1*(Q3+Z1*(Q4+Z1*(Q5))))))
                    AUX1=DSQRT(DLOG(1.0D0+(DESVPAD(J)/MEDIA(J))**2))
                    AUX2=DLOG(MEDIA(J))-0.50D0*(AUX1**2)
                    VARIAVEL(j)=EXP(Y*AUX1+AUX2)
                END IF
                IF (UI(J).GE.0.5 .AND. UI(J).LE.1.0) THEN
                    Z1= DSQRT(DLOG(1/((1-UI(J))**2)))          
                    Y=Z1+((P1+Z1*(P2+Z1*(P3+Z1*(P4+Z1*(P5)))))/(Q1+Z1*(Q2+Z1*(Q3+Z1*(Q4+Z1*(Q5))))))
                    AUX1=DSQRT(DLOG(1+(DESVPAD(J)/MEDIA(J))**2))
                    AUX2=DLOG(MEDIA(J))-0.5*(AUX1**2)
                    VARIAVEL(j)=EXP(Y*AUX1+AUX2)
                END IF
        END SELECT
    ENDDO
    
    UI(:)=VARIAVEL(:)

END SUBROUTINE