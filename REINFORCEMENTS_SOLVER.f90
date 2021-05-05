SUBROUTINE REINFORCEMENTS_SOLVER
    
    !
    !   SUBROUTINE FOR COLUMS EXCHANGE AND SYSTEM SOLVER IN THE CASE
    !   OF REINFORCED DOMAINS (SO FOR WITHOUT CRACKS)
    !
    
    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,J,II,JJ,CONT,KK,INFO,IPIV[ALLOCATABLE](:),IWORK[allocatable](:),ITER
    
    REAL*8,DIMENSION(:,:),POINTER::A,A_LINHA,B,B_LINHA
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE,TARGET::BGLOBAL
    
    REAL*8::SOLUCAO[ALLOCATABLE](:,:),F[ALLOCATABLE](:),SOLUTION[ALLOCATABLE](:)
    
    REAL*8::ANORM,AUX,MAT[ALLOCATABLE](:,:),RCOND,WORK[ALLOCATABLE](:)
    
    REAL::SWORK[ALLOCATABLE](:)
    
    character*1::norm
    
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        
        !STARTING VARIABLES
        II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
	    JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
    
        !ALLOCATE(A(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS),A_LINHA(3*N_NOS_MEF,3*N_COLLOCPOINTS),&
        !F(II),B(3*N_COLLOCPOINTS,II),B_LINHA(3*N_NOS_MEF,II),BGLOBAL(JJ,II),SOLUTION(JJ))
        
        ALLOCATE(F(II),BGLOBAL(JJ,II),SOLUTION(JJ))
        BGLOBAL=0.0D0
        
        A=>A_GLOBAL_ELASTIC(1:3*N_COLLOCPOINTS,1:3*N_COLLOCPOINTS)
        A_LINHA=>A_GLOBAL_ELASTIC(3*N_COLLOCPOINTS+1:3*N_COLLOCPOINTS+3*N_NOS_MEF,1:3*N_COLLOCPOINTS)
        
        B=>BGLOBAL(1:3*N_COLLOCPOINTS,1:II)
        B_LINHA=>BGLOBAL(3*N_COLLOCPOINTS+1:3*N_COLLOCPOINTS+3*N_NOS_MEF,1:II)
        
        !A=0.D0
	    !B=0.D0
	    !A_LINHA=0.D0
	    !B_LINHA=0.D0
	    SOLUCAO=0.D0
	    F=0.D0
	    CONT=0
	    U_MEF=0.D0
	    P_MEF=0.D0        
        
        DO I=1,N_COLLOCPOINTS
	        II=0
	        JJ=0
            IF (DUAL_BEM(I).EQ.'B') THEN
            
	            DO J=1,N_INTERFACE_POINTS
		            IF(I.EQ.NI(J,1)) THEN
			            II=-1
			            JJ=NI(J,2)
		            ELSE
			            IF(I.EQ.NI(J,2)) THEN
				            II=NI(J,1)
				            JJ=-1
			            ENDIF
		            ENDIF
	            ENDDO

	            IF(II.EQ.-1) THEN
			        A(:,(3*I-2))=A(:,(3*I-2))+H(:,(3*I-2))
			        A(:,3*I-1)=A(:,3*I-1)+H(:,3*I-1)
			        A(:,3*I)=A(:,3*I)+H(:,3*I)

			        A(:,(3*I-2))=A(:,(3*I-2))+H(:,(3*JJ-2))
			        A(:,3*I-1)=A(:,3*I-1)+H(:,3*JJ-1)
			        A(:,3*I)=A(:,3*I)+H(:,3*JJ)

                    A_LINHA(:,(3*I-2))=A_LINHA(:,(3*I-2))+H_FC(:,(3*I-2))
			        A_LINHA(:,3*I-1)=A_LINHA(:,3*I-1)+H_FC(:,3*I-1)
			        A_LINHA(:,3*I)=A_LINHA(:,3*I)+H_FC(:,3*I)

			        A_LINHA(:,(3*I-2))=A_LINHA(:,(3*I-2))+H_FC(:,(3*JJ-2))
			        A_LINHA(:,3*I-1)=A_LINHA(:,3*I-1)+H_FC(:,3*JJ-1)
			        A_LINHA(:,3*I)=A_LINHA(:,3*I)+H_FC(:,3*JJ)	    
	            ELSE
		            IF(JJ.EQ.-1) THEN
                        A(:,(3*I-2))=A(:,(3*I-2))+G(:,(3*I-2))
				        A(:,3*I-1)=A(:,3*I-1)+G(:,3*I-1)
				        A(:,3*I)=A(:,3*I)+G(:,3*I)

				        A(:,(3*I-2))=A(:,(3*I-2))-G(:,(3*II-2))
				        A(:,3*I-1)=A(:,3*I-1)-G(:,3*II-1)
				        A(:,3*I)=A(:,3*I)-G(:,3*II)

                        A_LINHA(:,(3*I-2))=A_LINHA(:,(3*I-2))+G_FC(:,(3*I-2))
				        A_LINHA(:,3*I-1)=A_LINHA(:,3*I-1)+G_FC(:,3*I-1)
				        A_LINHA(:,3*I)=A_LINHA(:,3*I)+G_FC(:,3*I)

				        A_LINHA(:,(3*I-2))=A_LINHA(:,(3*I-2))-G_FC(:,(3*II-2))
				        A_LINHA(:,3*I-1)=A_LINHA(:,3*I-1)-G_FC(:,3*II-1)
				        A_LINHA(:,3*I)=A_LINHA(:,3*I)-G_FC(:,3*II)
		            ELSE
			            IF(B_CONDITIONS(3*I-2).EQ.0) THEN
				            A(:,(3*I-2))=-G(:,(3*I-2))
				            A_LINHA(:,(3*I-2))=-G_FC(:,(3*I-2))

				            CONT=CONT+1
				            B(:,CONT)=-H(:,(3*I-2))
				            B_LINHA(:,CONT)=-H_FC(:,(3*I-2))
				            F(CONT)=U((3*I-2))
			            ELSE
				            A(:,(3*I-2))=H(:,(3*I-2))
				            A_LINHA(:,(3*I-2))=H_FC(:,(3*I-2))

				            CONT=CONT+1
				            B(:,CONT)=G(:,(3*I-2))
				            B_LINHA(:,CONT)=G_FC(:,(3*I-2))
				            F(CONT)=T((3*I-2))
			            ENDIF

			            IF(B_CONDITIONS(3*I-1).EQ.0) THEN
				            A(:,(3*I-1))=-G(:,(3*I-1))
				            A_LINHA(:,(3*I-1))=-G_FC(:,(3*I-1))

				            CONT=CONT+1
				            B(:,CONT)=-H(:,(3*I-1))
				            B_LINHA(:,CONT)=-H_FC(:,(3*I-1))
				            F(CONT)=U((3*I-1))
			            ELSE
				            A(:,(3*I-1))=H(:,(3*I-1))
				            A_LINHA(:,(3*I-1))=H_FC(:,(3*I-1))

				            CONT=CONT+1
				            B(:,CONT)=G(:,(3*I-1))
				            B_LINHA(:,CONT)=G_FC(:,(3*I-1))
				            F(CONT)=T((3*I-1))
                        ENDIF
                    
                        IF(B_CONDITIONS(3*I).EQ.0) THEN
				            A(:,(3*I))=-G(:,(3*I))
				            A_LINHA(:,(3*I))=-G_FC(:,(3*I))

				            CONT=CONT+1
				            B(:,CONT)=-H(:,(3*I))
				            B_LINHA(:,CONT)=-H_FC(:,(3*I))
				            F(CONT)=U((3*I))
			            ELSE
				            A(:,(3*I))=H(:,(3*I))
				            A_LINHA(:,(3*I))=H_FC(:,(3*I))

				            CONT=CONT+1
				            B(:,CONT)=G(:,(3*I))
				            B_LINHA(:,CONT)=G_FC(:,(3*I))
				            F(CONT)=T((3*I))
			            ENDIF
		            ENDIF
                ENDIF
            
            ELSE  ! PONTO SOBRE FISSURA
                !
                ! UNDER CONSTRUCTION
                !
            ENDIF
        ENDDO
    
        ! IDENTITY PART INSIDE AGLOBAL
        II=3*N_COLLOCPOINTS
        JJ=3*N_COLLOCPOINTS
        DO I=1,(3*N_NOS_MEF)
            A_GLOBAL_ELASTIC((II+I),(JJ+I))=1.0D0
        ENDDO

        ! DETERMINANDO VETOR A DIREITA ___________________________________________
        II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
 	    JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
	    CALL DMURRV(JJ,II,BGLOBAL,JJ,II,F,1,JJ,SOLUTION)
	    DEALLOCATE(F)
	  
	    ALLOCATE(F(JJ))
	    F(:)=SOLUTION(:)

        ! PARTES DESATIVADAS DO CODIGO ______________________
            ! VERIFICANDO NUMERO DE CONDIÇÃO
        !ANORM=0.0D0
        !DO  J=1,JJ
        !    AUX=0.0D0
        !    DO I=1,JJ
        !        AUX=AUX+DABS(AGLOBAL(I,J))
        !    ENDDO
        !    IF (AUX.GT.ANORM) ANORM=AUX
        !ENDDO
        !norm='1'
        !ALLOCATE(MAT(JJ,JJ),IPIV(JJ),WORK(4*JJ),IWORK(JJ))
        !MAT=AGLOBAL
        !call dgetrf(JJ,JJ,MAT,JJ,IPIV,INFO)
        !
        !call dgecon(norm,JJ,MAT,JJ,ANORM,RCOND,WORK,IWORK,INFO)
        !WRITE(*,*)'NUMBER CONDITION:',RCOND
        !READ(*,*)
        ! __________________________________________________
                
        IF (EP_ANALYSIS_REINF) THEN
            
            ALLOCATE(IPIV_A_GLOBAL_EP(JJ),A_GLOBAL_EP_LU(JJ,JJ),B_GLOBAL_EP(JJ,II))
            
            A_GLOBAL_EP_LU = A_GLOBAL_ELASTIC
                        
            ! OBTENDO A_GLOBAL_EP_LU EMCIMA DA MATRIZ AGLOBAL
            CALL DGETRF(JJ,JJ,A_GLOBAL_EP_LU,JJ,IPIV_A_GLOBAL_EP,INFO)

            B_GLOBAL_EP=BGLOBAL
            F_ALVO(:)=F(:)
            
            DEALLOCATE(BGLOBAL)
            
            SOLUTION(:)=F(:)
            CALL DGETRS('N',JJ,1,A_GLOBAL_EP_LU,JJ,IPIV_A_GLOBAL_EP,SOLUTION,JJ,INFO)
                        
        ELSE IF (SLIP_REINF) THEN
            
            ALLOCATE(IPIV_A_GLOBAL_SLIP_F(JJ),A_GLOBAL_SLIP_F_LU(JJ,JJ),B_GLOBAL_SLIP(JJ,II),A_GLOBAL_SLIP_F(JJ,JJ))
            
            A_GLOBAL_SLIP_F = A_GLOBAL_ELASTIC
            
            ! OBTENDO MATRIZ COM FATORACAO LU
            A_GLOBAL_SLIP_F_LU=A_GLOBAL_ELASTIC
            CALL DGETRF(JJ,JJ,A_GLOBAL_SLIP_F_LU,JJ,IPIV_A_GLOBAL_SLIP_F,INFO)
            B_GLOBAL_SLIP=BGLOBAL
            
            IF (N_NODAL_FORCES.NE.0) THEN  ! APLICANDO FORÇAS NODAIS CONCENTRADAS
                
                ! PRIMEIRA FORMA DE APLICAR: SÓ COLOCA NO LADO DIREITO DO SISTEMA
                !SOLUTION(:)=0.0D0
                !II=3*N_COLLOCPOINTS+3*N_NOS_MEF
                !DO I=1,N_NODAL_FORCES
                !    KK=3*(NODAL_FORCES_PLACE(I,1)-1)+NODAL_FORCES_PLACE(I,2)
                !    SOLUTION(II+KK)=NODAL_FORCES_VALUE(I)
                !ENDDO
                !SOLUTION(:)=SOLUTION(:)+F(:)
                !CALL DGETRS('N',JJ,1,AGLOBAL,JJ,IPIV_A_GLOBAL_SLIP_F,SOLUTION,JJ,INFO)
                
                ! SEGUNDA FORMA DE APLICAR: MULTIPLICA PELA SOLUCAO FUNDAMENTAL (COMO F CONCENTRADO DE ENRIQUECIMENTO)
                !CALL CONCENTRATED_F_1DBEM(SOLUTION,JJ)
                !SOLUTION(:)=SOLUTION(:)+F(:)
                !CALL DGETRS('N',JJ,1,AGLOBAL,JJ,IPIV_A_GLOBAL_SLIP_F,SOLUTION,JJ,INFO)
                
                ! TERCEIRA FORMA DE APLICAR: COLOCA NO CONNECTION ELEMENT (DESLOCAMENTO APLICADO)
                DEALLOCATE(SOLUTION)
                ALLOCATE(SOLUTION(JJ+3*N_NODAL_FORCES))
                CALL CONCENTRATED_U_CONNECTION(SOLUTION,JJ)
                DO I=1,JJ
                    SOLUTION(I)=SOLUTION(I)+F(I)
                ENDDO
                JJ=JJ+3*N_NODAL_FORCES
                CALL DGETRS('N',JJ,1,A_GLOBAL_SLIP_F_LU,JJ,IPIV_A_GLOBAL_SLIP_F,SOLUTION,JJ,INFO)  
                
            ELSE
                SOLUTION(:)=F(:)
                CALL DGETRS('N',JJ,1,A_GLOBAL_SLIP_F_LU,JJ,IPIV_A_GLOBAL_SLIP_F,SOLUTION,JJ,INFO)
            ENDIF
            
            DEALLOCATE(BGLOBAL)
            
        ELSE
            ! USING IMSL SOLVER
            !CALL DLSARG(JJ,AGLOBAL,JJ,F,1,SOLUTION)

            !USING DGESV 
            WRITE(*,*)"     SOLVING THE PROBLEM"
            ALLOCATE(IPIV(JJ))
            SOLUTION(:)=F(:)
            CALL DGESV(JJ,1,A_GLOBAL_ELASTIC,JJ,IPIV,SOLUTION,JJ,INFO)
            
            DEALLOCATE(BGLOBAL)
        ENDIF
                
        IF (INFO.NE.0) THEN
            WRITE(*,'(A,I7,A)')'*** PROBLEM IN THE SYSTEM SOLVER, LINE:',INFO,'*****'
            WRITE(*,'(A)')'VERIFY SYSTEM SOLVER (DGESV OR DGETRS) IN REINF_SOLVER SUBROUTINE'
            WRITE(*,'(A)')'*********************************************************'
            READ(*,*)
        ENDIF
        
        ! _______________________________________________________________
    
        ! SAVING U AND T FROM SOLUTION __________________________________
        DO I=1,N_COLLOCPOINTS
	        II=0
	        JJ=0
	        IF(DUAL_BEM(I).EQ.'B')THEN  ! POINT AT BOUNDARY
		        DO J=1,N_INTERFACE_POINTS
			        IF(I.EQ.NI(J,1)) THEN
				        II=-1
				        JJ=NI(J,2)
			        ELSE
				        IF(I.EQ.NI(J,2)) THEN
					        II=NI(J,1)
					        JJ=-1
				        ENDIF
			        ENDIF
		        ENDDO	
    
		        IF(II.EQ.-1) THEN
			        U(3*I-2)=SOLUTION(3*I-2)
			        U(3*I-1)=SOLUTION(3*I-1)
			        U(3*I)=SOLUTION(3*I)
    
			        U(3*JJ-2)=SOLUTION(3*I-2)
			        U(3*JJ-1)=SOLUTION(3*I-1)
			        U(3*JJ)=SOLUTION(3*I)
		        ELSE
			        IF(JJ.EQ.-1) THEN
				        T(3*I-2)=-SOLUTION(3*I-2)
				        T(3*I-1)=-SOLUTION(3*I-1)
			            T(3*I)=-SOLUTION(3*I)
    
				        T(3*II-2)=SOLUTION(3*I-2)
				        T(3*II-1)=SOLUTION(3*I-1)
				        T(3*II)=SOLUTION(3*I)
			        ELSE
				        IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					        T(3*I-2)=SOLUTION(3*I-2)
				        ELSE
					        U(3*I-2)=SOLUTION(3*I-2)
				        ENDIF
    
				        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					        T(3*I-1)=SOLUTION(3*I-1)
				        ELSE
					        U(3*I-1)=SOLUTION(3*I-1)
				        ENDIF
    
				        IF (B_CONDITIONS(3*I).EQ.0) THEN
					        T(3*I)=SOLUTION(3*I)
				        ELSE
					        U(3*I)=SOLUTION(3*I)
				        ENDIF
			        ENDIF
		        ENDIF
            ELSE   ! POINT AT CRACK CONTOURN
                !
                ! UNDER CONSTRUCTION
                !
            ENDIF
        ENDDO
    
        II=3*N_COLLOCPOINTS
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
        DO I=1,(3*N_NOS_MEF)
            U_MEF(I)=SOLUTION(II+I)
            P_MEF(I)=SOLUTION(JJ+I)
        ENDDO
        
        IF (SLIP_REINF.AND.(N_NODAL_FORCES.NE.0)) THEN         ! FORÇA DE REAÇÃO NOS NÓS PUCHADOS
            JJ=3*N_COLLOCPOINTS+6*N_NOS_MEF
            DO I=1,3*N_NODAL_FORCES
                F_APP_DISP(I)=SOLUTION(JJ+I)
                !WRITE(*,*)SOLUTION(JJ+I)
            ENDDO
        ENDIF
        
        DEALLOCATE(F,SOLUTION)
        ! _______________________________________________________________
    
    ELSE    ! CASO SEM REINFORCEMENTS
            
        WRITE(*,*)''
        WRITE(*,*)'****************************** ERRO *********************************'
        WRITE(*,*)' SUBROTINA REINFORCEMENTS_SOLVER USADA PARA CASO SEM REINFORCEMENTS'
        WRITE(*,*)' PARTE NAO PROGRAMADA - VERIFICAR PORQUE ESTÁ PASSANDO POR ESSA SUB'
        WRITE(*,*)''
        READ(*,*)
        
    ENDIF

    
END SUBROUTINE


! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

    
SUBROUTINE REINF_SOLVER_PLAST(N_PASSO,N_ITER)
    
    !
    ! SUBROUTINE TO APPLY THE LOAD STEP, SOLVE THE SYSTEM AND FIND DELTA U AND DELTA P FOR ELASTOPLASTIC ANALYSIS
    ! 
    ! N_PASSO: CURRENT NUMBER OF LOAD STEP
    ! N_ITER: CURRENT NUMBER OF ITERATION (INSIDE THE LOAD STEP N_PASSO)
    !

    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER,INTENT(IN)::N_PASSO,N_ITER
    
    INTEGER::I,J,K,II,JJ,KK,CONT,INFO
    
    REAL*8::AUX1,AUX2,AUX3
    
    REAL*8,DIMENSION(:),ALLOCATABLE::F,SOLUTION
    
    
    IF ((N_PASSO.EQ.1).AND.(N_ITER.EQ.1)) THEN  ! FIRST OF ALL, FIND MATRIXES ______________
        
        CALL REINFORCEMENTS_SOLVER
        
    ELSE IF (N_ITER.EQ.1) THEN ! FIRST ITERATION WITH MATRIXES READY, APPLY U AND T ________
        
        ! Starting variables
        II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
        ALLOCATE(SOLUTION(JJ),F(II))
	    SOLUTION=0.D0
	    F=0.D0
	    CONT=0
        
        DO I=1,N_COLLOCPOINTS
		    II=0
		    JJ=0
		    IF(DUAL_BEM(I).EQ.'B')THEN
		        DO J=1,N_INTERFACE_POINTS
			        IF(I.EQ.NI(J,1)) THEN
				        II=-1
				        JJ=NI(J,2)
			        ELSE
				        IF(I.EQ.NI(J,2)) THEN
					        II=NI(J,1)
					        JJ=-1
				        ENDIF
			        ENDIF
                ENDDO 
                
                IF(II.EQ.-1) THEN
                    !
                ELSE
                    IF(JJ.EQ.-1) THEN
                        !
                    ELSE  
                        IF(B_CONDITIONS(3*I-2).EQ.0) THEN
				            CONT=CONT+1
				            F(CONT)=U((3*I-2))
			            ELSE
				            CONT=CONT+1
				            F(CONT)=T((3*I-2))
			            ENDIF

			            IF(B_CONDITIONS(3*I-1).EQ.0) THEN
				            CONT=CONT+1
				            F(CONT)=U((3*I-1))
			            ELSE
				            CONT=CONT+1
				            F(CONT)=T((3*I-1))
                        ENDIF
                    
                        IF(B_CONDITIONS(3*I).EQ.0) THEN
				            CONT=CONT+1
				            F(CONT)=U((3*I))
			            ELSE
				            CONT=CONT+1
				            F(CONT)=T((3*I))
			            ENDIF
			        ENDIF
		        ENDIF
            ELSE  ! PONTO SOBRE FISSURA
                !
                ! UNDER CONSTRUCTION
                !
            ENDIF
        ENDDO  
        
        ! MULTIPLICANDO B*F
 	    II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
        !CALL DMURRV(JJ,II,BGLOBAL,JJ,II,FGLOBAL,1,JJ,SOLUTION)
        CALL DGEMV('N',JJ,II,1.0D0,B_GLOBAL_EP,JJ,F,1,0.0D0,SOLUTION,1)
        
        F_ALVO(:)=SOLUTION(:)
        
        ! RESOLVENDO O SISTEMA COM LU JA FATORADO:  
        CALL DGETRS('N',JJ,1,A_GLOBAL_EP_LU,JJ,IPIV_A_GLOBAL_EP,SOLUTION,JJ,INFO)
        
        ! SAVING U AND T FROM SOLUTION __________________________________
        DO I=1,N_COLLOCPOINTS
	        II=0
	        JJ=0
	        IF(DUAL_BEM(I).EQ.'B')THEN  ! POINT AT BOUNDARY
		        DO J=1,N_INTERFACE_POINTS
			        IF(I.EQ.NI(J,1)) THEN
				        II=-1
				        JJ=NI(J,2)
			        ELSE
				        IF(I.EQ.NI(J,2)) THEN
					        II=NI(J,1)
					        JJ=-1
				        ENDIF
			        ENDIF
		        ENDDO	
    
		        IF(II.EQ.-1) THEN
			        U(3*I-2)=SOLUTION(3*I-2)
			        U(3*I-1)=SOLUTION(3*I-1)
			        U(3*I)=SOLUTION(3*I)
    
			        U(3*JJ-2)=SOLUTION(3*I-2)
			        U(3*JJ-1)=SOLUTION(3*I-1)
			        U(3*JJ)=SOLUTION(3*I)
		        ELSE
			        IF(JJ.EQ.-1) THEN
				        T(3*I-2)=-SOLUTION(3*I-2)
				        T(3*I-1)=-SOLUTION(3*I-1)
			            T(3*I)=-SOLUTION(3*I)
    
				        T(3*II-2)=SOLUTION(3*I-2)
				        T(3*II-1)=SOLUTION(3*I-1)
				        T(3*II)=SOLUTION(3*I)
			        ELSE
				        IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					        T(3*I-2)=SOLUTION(3*I-2)
				        ELSE
					        U(3*I-2)=SOLUTION(3*I-2)
				        ENDIF
    
				        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					        T(3*I-1)=SOLUTION(3*I-1)
				        ELSE
					        U(3*I-1)=SOLUTION(3*I-1)
				        ENDIF
    
				        IF (B_CONDITIONS(3*I).EQ.0) THEN
					        T(3*I)=SOLUTION(3*I)
				        ELSE
					        U(3*I)=SOLUTION(3*I)
				        ENDIF
			        ENDIF
		        ENDIF
            ELSE   ! POINT AT CRACK CONTOURN
                !
                ! UNDER CONSTRUCTION
                !
            ENDIF
        ENDDO
        II=3*N_COLLOCPOINTS
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
        DO I=1,(3*N_NOS_MEF)
            U_MEF(I)=SOLUTION(II+I)
            P_MEF(I)=SOLUTION(JJ+I)
        ENDDO 
        
        DEALLOCATE(F,SOLUTION)
        
    ELSE ! NON FIRST ITERATION, APPLY F_DES ________________________________________________
        
        ! Starting variables
        II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
        ALLOCATE(SOLUTION(JJ))

        SOLUTION(:)=F_DES(:)
        
        ! RESOLVENDO O SISTEMA COM LU JA FATORADO:
        CALL DGETRS('N',JJ,1,A_GLOBAL_EP_LU,JJ,IPIV_A_GLOBAL_EP,SOLUTION,JJ,INFO)
        
        ! MUDOU
        U=0.0D0
        T=0.0D0
        
        ! SAVING U AND T FROM SOLUTION __________________________________
        DO I=1,N_COLLOCPOINTS
	        II=0
	        JJ=0
	        IF(DUAL_BEM(I).EQ.'B')THEN  ! POINT AT BOUNDARY
		        DO J=1,N_INTERFACE_POINTS
			        IF(I.EQ.NI(J,1)) THEN
				        II=-1
				        JJ=NI(J,2)
			        ELSE
				        IF(I.EQ.NI(J,2)) THEN
					        II=NI(J,1)
					        JJ=-1
				        ENDIF
			        ENDIF
		        ENDDO	
    
		        IF(II.EQ.-1) THEN
			        U(3*I-2)=SOLUTION(3*I-2)
			        U(3*I-1)=SOLUTION(3*I-1)
			        U(3*I)=SOLUTION(3*I)
    
			        U(3*JJ-2)=SOLUTION(3*I-2)
			        U(3*JJ-1)=SOLUTION(3*I-1)
			        U(3*JJ)=SOLUTION(3*I)
		        ELSE
			        IF(JJ.EQ.-1) THEN
				        T(3*I-2)=-SOLUTION(3*I-2)
				        T(3*I-1)=-SOLUTION(3*I-1)
			            T(3*I)=-SOLUTION(3*I)
    
				        T(3*II-2)=SOLUTION(3*I-2)
				        T(3*II-1)=SOLUTION(3*I-1)
				        T(3*II)=SOLUTION(3*I)
			        ELSE
				        IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					        T(3*I-2)=SOLUTION(3*I-2)
				        ELSE
					        U(3*I-2)=SOLUTION(3*I-2)
				        ENDIF
    
				        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					        T(3*I-1)=SOLUTION(3*I-1)
				        ELSE
					        U(3*I-1)=SOLUTION(3*I-1)
				        ENDIF
    
				        IF (B_CONDITIONS(3*I).EQ.0) THEN
					        T(3*I)=SOLUTION(3*I)
				        ELSE
					        U(3*I)=SOLUTION(3*I)
				        ENDIF
			        ENDIF
		        ENDIF
            ELSE   ! POINT AT CRACK CONTOURN
                !
                ! UNDER CONSTRUCTION
                !
            ENDIF
        ENDDO
        II=3*N_COLLOCPOINTS
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
        DO I=1,(3*N_NOS_MEF)
            U_MEF(I)=SOLUTION(II+I)
            P_MEF(I)=SOLUTION(JJ+I)
        ENDDO 
        
        DEALLOCATE(SOLUTION)
        
    ENDIF ! ________________________________________________________________________________
    
    ! Aculumando grandezas
    U_ACUM_PL(:)=U_ACUM_PL(:)+U(:)
    P_ACUM_PL(:)=P_ACUM_PL(:)+T(:)
    U_MEF_ACUM(:) = U_MEF_ACUM(:) + U_MEF(:)
    
    ! Calcular esforços normais nodais
    CALL REINFORCEMENTS_STRESS
    

END SUBROUTINE REINF_SOLVER_PLAST
    
    
    

! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

    
SUBROUTINE REINF_SOLVER_SLIP(N_PASSO,N_ITER,APPLIED_C,CHANGE_MATRIX)
    
    !
    ! SUBROUTINE TO APPLY THE LOAD STEP, SOLVE THE SYSTEM AND FIND DELTA U AND DELTA P FOR BOND-SLIP ANALYSIS
    ! 
    ! N_PASSO: CURRENT NUMBER OF LOAD STEP
    ! N_ITER: CURRENT NUMBER OF ITERATION (INSIDE THE LOAD STEP N_PASSO)
    ! APLLIED_C: VETOR APPLIED IN THE ITERATIONS OF BOND-SLIP ANALYSIS (F_C OR SLIP VALUES) 
    ! CHANGE_MATRIX: .TRUE. IF IT IS NECESSARY TO CHANGE MATRIX AND .FALSE. IF THE PREVIOUS MATRIX CAN BE RE-UTILIZED    
    !

    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    LOGICAL,INTENT(IN)::CHANGE_MATRIX
    INTEGER,INTENT(IN)::N_PASSO,N_ITER
    REAL*8,INTENT(IN)::APPLIED_C(3*N_NOS_MEF)
    
    INTEGER::I,J,K,II,JJ,KK,CONT,INFO,NI_SLIP_G[ALLOCATABLE](:)
    
    LOGICAL::COL_HAS_ZEROS(3)
    
    REAL*8::AUX1,AUX2,AUX3
        
    REAL*8,DIMENSION(:),ALLOCATABLE::F,SOLUTION,IPIV,FGLOBAL,F_C,DELTA_S_GLOBAL,DELTA_S_LOCAL,SOL
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::AGLOBAL,BGLOBAL
    
    
    IF (N_ITER.EQ.1) THEN                   ! FIRST ITERATION, APPLY BOUNDARY CONDIDITIONS AND FIN P_MEF ______________
        
        IF (N_PASSO.EQ.1) THEN              ! CREATE MATRIXES
            
            CALL REINFORCEMENTS_SOLVER
            
        ELSE                                ! USE MATRIXES ALREADY CREATED
            
            ! Starting variables
            II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
            ALLOCATE(SOLUTION(JJ),F(II),F_C(3*N_NODAL_FORCES),SOL(JJ+3*N_NODAL_FORCES))
	        SOLUTION=0.D0
            SOL=0.0D0
	        F=0.D0
	        CONT=0
        
            DO I=1,N_COLLOCPOINTS
		        II=0
		        JJ=0
		        IF(DUAL_BEM(I).EQ.'B')THEN
		            DO J=1,N_INTERFACE_POINTS
			            IF(I.EQ.NI(J,1)) THEN
				            II=-1
				            JJ=NI(J,2)
			            ELSE
				            IF(I.EQ.NI(J,2)) THEN
					            II=NI(J,1)
					            JJ=-1
				            ENDIF
			            ENDIF
                    ENDDO 
                
                    IF(II.EQ.-1) THEN
                        !
                    ELSE
                        IF(JJ.EQ.-1) THEN
                            !
                        ELSE  
                            IF(B_CONDITIONS(3*I-2).EQ.0) THEN
				                CONT=CONT+1
				                F(CONT)=U((3*I-2))
			                ELSE
				                CONT=CONT+1
				                F(CONT)=T((3*I-2))
			                ENDIF

			                IF(B_CONDITIONS(3*I-1).EQ.0) THEN
				                CONT=CONT+1
				                F(CONT)=U((3*I-1))
			                ELSE
				                CONT=CONT+1
				                F(CONT)=T((3*I-1))
                            ENDIF
                    
                            IF(B_CONDITIONS(3*I).EQ.0) THEN
				                CONT=CONT+1
				                F(CONT)=U((3*I))
			                ELSE
				                CONT=CONT+1
				                F(CONT)=T((3*I))
			                ENDIF
			            ENDIF
		            ENDIF
                ELSE  ! PONTO SOBRE FISSURA
                    !
                    ! UNDER CONSTRUCTION
                    !
                ENDIF
            ENDDO  
        
            ! MULTIPLICANDO B*F
 	        II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
            !CALL DMURRV(JJ,II,BGLOBAL,JJ,II,FGLOBAL,1,JJ,SOLUTION)
            CALL DGEMV('N',JJ,II,1.0D0,B_GLOBAL_SLIP,JJ,F,1,0.0D0,SOLUTION,1)
            
            ! SOMANDO PARCELA DO DESLOCAMENTO APLICADO
            IF (N_NODAL_FORCES.NE.0) THEN
                F_C=0.0D0
                DO J=1,N_NODAL_FORCES
                    F_C(3*(J-1)+NODAL_FORCES_PLACE(J,2))=NODAL_FORCES_VALUE(J)
                ENDDO
                II=3*N_NODAL_FORCES
                JJ=JJ+II
                CALL DGEMV('N',JJ,II,1.0D0,B_CONC_F,JJ,F_C,1,0.0D0,SOL,1)
                DO J=1,3*N_COLLOCPOINTS+6*N_NOS_MEF
                    SOL(J)=SOL(J)+SOLUTION(J)
                ENDDO         
                DEALLOCATE(SOLUTION)
                ALLOCATE(SOLUTION(JJ))
                SOLUTION(:)=SOL(:)
            ENDIF
        
            ! RESOLVENDO O SISTEMA COM LU JA FATORADO:  
            CALL DGETRS('N',JJ,1,A_GLOBAL_SLIP_F_LU,JJ,IPIV_A_GLOBAL_SLIP_F,SOLUTION,JJ,INFO)
            
            IF (INFO.NE.0) THEN
                WRITE(*,'(A,I7,A)')'*** PROBLEM IN THE SYSTEM SOLVER, LINE:',INFO,'*****'
                WRITE(*,'(A)')'VERIFY DGETRS MATRIXES IN THE REINF_SOLVER_SLIP SUBROUTINE'
                WRITE(*,'(A)')'*********************************************************'
                READ(*,*)
            ENDIF
        
            ! SAVING U AND T FROM SOLUTION __________________________________
            DO I=1,N_COLLOCPOINTS
	            II=0
	            JJ=0
	            IF(DUAL_BEM(I).EQ.'B')THEN  ! POINT AT BOUNDARY
		            DO J=1,N_INTERFACE_POINTS
			            IF(I.EQ.NI(J,1)) THEN
				            II=-1
				            JJ=NI(J,2)
			            ELSE
				            IF(I.EQ.NI(J,2)) THEN
					            II=NI(J,1)
					            JJ=-1
				            ENDIF
			            ENDIF
		            ENDDO	
    
		            IF(II.EQ.-1) THEN
			            U(3*I-2)=SOLUTION(3*I-2)
			            U(3*I-1)=SOLUTION(3*I-1)
			            U(3*I)=SOLUTION(3*I)
    
			            U(3*JJ-2)=SOLUTION(3*I-2)
			            U(3*JJ-1)=SOLUTION(3*I-1)
			            U(3*JJ)=SOLUTION(3*I)
		            ELSE
			            IF(JJ.EQ.-1) THEN
				            T(3*I-2)=-SOLUTION(3*I-2)
				            T(3*I-1)=-SOLUTION(3*I-1)
			                T(3*I)=-SOLUTION(3*I)
    
				            T(3*II-2)=SOLUTION(3*I-2)
				            T(3*II-1)=SOLUTION(3*I-1)
				            T(3*II)=SOLUTION(3*I)
			            ELSE
				            IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					            T(3*I-2)=SOLUTION(3*I-2)
				            ELSE
					            U(3*I-2)=SOLUTION(3*I-2)
				            ENDIF
    
				            IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					            T(3*I-1)=SOLUTION(3*I-1)
				            ELSE
					            U(3*I-1)=SOLUTION(3*I-1)
				            ENDIF
    
				            IF (B_CONDITIONS(3*I).EQ.0) THEN
					            T(3*I)=SOLUTION(3*I)
				            ELSE
					            U(3*I)=SOLUTION(3*I)
				            ENDIF
			            ENDIF
		            ENDIF
                ELSE   ! POINT AT CRACK CONTOURN
                    !
                    ! UNDER CONSTRUCTION
                    !
                ENDIF
            ENDDO
            II=3*N_COLLOCPOINTS
            JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
            DO I=1,(3*N_NOS_MEF)
                U_MEF(I)=SOLUTION(II+I)
                P_MEF(I)=SOLUTION(JJ+I)
            ENDDO 
            
            IF (N_NODAL_FORCES.NE.0) THEN         ! FORÇA DE REAÇÃO NOS NÓS PUCHADOS
                JJ=3*N_COLLOCPOINTS+6*N_NOS_MEF
                DO I=1,3*N_NODAL_FORCES
                    F_APP_DISP(I)=SOLUTION(JJ+I)
                    !WRITE(*,*)SOLUTION(JJ+I)
                ENDDO
            ENDIF
        
            DEALLOCATE(F,SOLUTION,F_C,SOL)
            
        ENDIF   
       
    ELSE ! ______________________________ NON FIRST ITERATION, APPLY F_C AND FIND DELTA_S AND P_MEF _________________________________________
        
        IF (.NOT.ALLOCATED(B_GLOBAL_SLIP_F)) THEN              ! CREATE B MATRIX
            JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
            ALLOCATE(B_GLOBAL_SLIP_F(JJ+3*N_NODAL_FORCES,3*N_NOS_MEF))
            ALLOCATE(B_GLOBAL_SLIP_S(JJ+3*N_NODAL_FORCES,3*N_NOS_MEF))
            B_GLOBAL_SLIP_F=0.0D0 
            ! MONTANDO MATRIZ A DIREITA
            DO J=1,3*N_NOS_MEF
                IF (N_NODAL_FORCES.EQ.0) THEN
                    II=3*N_COLLOCPOINTS
                    KK=3*N_NOS_MEF
                    DO I=1,3*N_NOS_MEF
                        B_GLOBAL_SLIP_F(II+KK+I,J)=KGLOBAL(I,J)
                    ENDDO
                ELSE
                    II=3*N_COLLOCPOINTS
                    KK=3*N_NOS_MEF
                    DO I=1,3*N_NOS_MEF+3*N_NODAL_FORCES
                        B_GLOBAL_SLIP_F(II+KK+I,J)=KGLOBAL_NEW(I,J)
                    ENDDO
                ENDIF
            ENDDO 
            
            ! ALOCANDO MATRIZES A
            JJ=3*N_COLLOCPOINTS+6*N_NOS_MEF+3*N_NODAL_FORCES
            ALLOCATE(A_GLOBAL_SLIP_S(JJ,JJ),A_GLOBAL_SLIP_S_LU(JJ,JJ),IPIV_A_GLOBAL_SLIP_S(JJ))
    
        ENDIF   ! ________ PROCEDIMENTO ABAIXO AGORA FAZ PARA TDS AS ITERACOES ____________
          
        ! PASSA NI_SLIP_MEF PARA GLOBAL         
        ALLOCATE(NI_SLIP_G(3*N_NOS_MEF))
        CALL ROTATE_NI_TO_GLOBAL(NI_SLIP_MEF,NI_SLIP_G)
        IF (CHANGE_MATRIX) THEN ! CALCULA MATRIZES INTEIRAS
            
            JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
            II=3*N_NODAL_FORCES
            ALLOCATE(AGLOBAL(JJ+II,JJ+II),BGLOBAL(JJ+II,3*N_NOS_MEF),FGLOBAL(JJ+II),DELTA_S_GLOBAL(3*N_NOS_MEF),&
                DELTA_S_LOCAL(TOTAL_ELEM_NODES_MEF),SOLUTION(JJ+II))
            BGLOBAL=0.0D0
            FGLOBAL=0.0D0
            SOLUTION=0.0D0
            
            ! COPIANDO TERMOS EM COMUM COM A_GLOBAL_SLIP_F:
            DO J=1,3*N_COLLOCPOINTS+3*N_NOS_MEF
                DO I=1,3*N_COLLOCPOINTS+6*N_NOS_MEF+3*N_NODAL_FORCES
                    AGLOBAL(I,J)=A_GLOBAL_SLIP_F(I,J)
                ENDDO
            ENDDO
            JJ=3*N_COLLOCPOINTS+6*N_NOS_MEF
            DO J=1,3*N_NODAL_FORCES
                DO I=1,3*N_COLLOCPOINTS+6*N_NOS_MEF+3*N_NODAL_FORCES
                    AGLOBAL(I,JJ+J)=A_GLOBAL_SLIP_F(I,JJ+J)
                ENDDO
            ENDDO
        
            ! MONTANDO MATRIZES COM TROCA DE COLUNAS (BASEADO EM NI_SLIP_MEF)
            II=3*N_COLLOCPOINTS+3*N_NOS_MEF
            JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF+3*N_NODAL_FORCES
            DO J=1,3*N_NOS_MEF        
                IF (NI_SLIP_G(J).EQ.0) THEN ! NÓ COLADO 
                    
                    AGLOBAL(:,II+J)=A_GLOBAL_SLIP_F(:,II+J)
            
                    BGLOBAL(:,J)=B_GLOBAL_SLIP_F(:,J)

                ELSE IF (NI_SLIP_G(J).EQ.1) THEN ! NÓ DESLOCADO
                   
                    AGLOBAL(:,II+J)=-B_GLOBAL_SLIP_F(:,J)

                    BGLOBAL(:,J)=-A_GLOBAL_SLIP_F(:,II+J)

                ELSE 
                    WRITE(*,*)'**** ERROR IN BOND-SLIP ITERATIONS, VERIFY NI_SLIP_G *****'
                    WRITE(*,*)'************************************************************'
                ENDIF
                !THE STIFFENER MAY NOT HAVE RESISTENCE IN ONE OF THE DIRECTIONS, IN THIS CASE BOTH MATRIX K AND G WILL HAVE A COLUM OF ZEROS
                !TO SOLVE THE FINAL SYSTEM, IT WILL BE ADOPTED A VALUE OF ZERO IN THE DIAGONAL (UNDERSTAND LINE = COLUM)
                !COL_HAS_ZEROS=.TRUE.
                !DO I=1,JJ
                !    IF (AGLOBAL(I,II+J).NE.0.0D0) THEN
                !        COL_HAS_ZEROS(1)=.FALSE.
                !    ENDIF
                !    IF (.NOT.COL_HAS_ZEROS(1)) EXIT
                !    !IF ((.NOT.COL_HAS_ZEROS(1)).AND.(.NOT.COL_HAS_ZEROS(2)).AND.(.NOT.COL_HAS_ZEROS(3))) EXIT
                !ENDDO
                !IF (COL_HAS_ZEROS(1)) THEN
                !    !WRITE(*,*)'COL:',II+3*J-2,'HAS ZEROS'
                !    !AGLOBAL(II+3*J-2,II+3*J-2)=1000.0D0
                !    AGLOBAL(II+J,II+J)=1000.0D0
                !ENDIF
               
            ENDDO
            
            ! SALVANDO BGLOBAL
            B_GLOBAL_SLIP_S=BGLOBAL
            
            ! MULTIPLICANDO POR APPLIED_C:
            II=3*N_NOS_MEF
            JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF+3*N_NODAL_FORCES
            CALL DGEMV('N',JJ,II,1.0D0,BGLOBAL,JJ,APPLIED_C,1,0.0D0,FGLOBAL,1)
        
            ! OBTENDO MATRIZ COM FATORACAO LU
            JJ=3*N_COLLOCPOINTS+6*N_NOS_MEF+3*N_NODAL_FORCES
            A_GLOBAL_SLIP_S=AGLOBAL
            CALL DGETRF(JJ,JJ,AGLOBAL,JJ,IPIV_A_GLOBAL_SLIP_S,INFO)
            A_GLOBAL_SLIP_S_LU=AGLOBAL
                    
            IF (INFO.NE.0) THEN
                WRITE(*,'(A,I7,A)')'*** PROBLEM IN THE LU FACTORIZATION, LINE:',INFO,'*****'
                WRITE(*,'(A)')'VERIFY DGETRF MATRIXES IN THE REINF_SOLVER_SLIP SUBROUTINE'
                WRITE(*,'(A)')'*********************************************************'
                READ(*,*)
            ENDIF
        
            ! RESOLVENDO O SISTEMA COM LU JA FATORADO:
            SOLUTION(:)=FGLOBAL(:)
            CALL DGETRS('N',JJ,1,A_GLOBAL_SLIP_S_LU,JJ,IPIV_A_GLOBAL_SLIP_S,SOLUTION,JJ,INFO)

            ! RESOLVENDO DIRETO
            !ALLOCATE(IPIV(JJ))
            !SOLUTION(:)=FGLOBAL(:)
            !CALL DGESV(JJ,1,AGLOBAL,JJ,IPIV,SOLUTION,JJ,INFO) 
        
            IF (INFO.NE.0) THEN
                WRITE(*,'(A,I7,A)')'*** PROBLEM IN THE SYSTEM SOLVER, LINE:',INFO,'*****'
                WRITE(*,'(A)')'VERIFY DGETRS MATRIXES IN THE REINF_SOLVER_SLIP SUBROUTINE'
                WRITE(*,'(A)')'*********************************************************'
                READ(*,*)
            ENDIF
        
        ELSE ! NAO MUDA MATRIZES, UTILIZA JÁ INVERTIDAS  
            
            JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
            II=3*N_NODAL_FORCES
            ALLOCATE(FGLOBAL(JJ+II),DELTA_S_GLOBAL(3*N_NOS_MEF),DELTA_S_LOCAL(TOTAL_ELEM_NODES_MEF),&
                SOLUTION(JJ+II))
            FGLOBAL=0.0D0
            SOLUTION=0.0D0
                        
            ! MULTIPLICANDO POR APPLIED_C:
            II=3*N_NOS_MEF
            JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF+3*N_NODAL_FORCES
            CALL DGEMV('N',JJ,II,1.0D0,B_GLOBAL_SLIP_S,JJ,APPLIED_C,1,0.0D0,FGLOBAL,1)
                
            ! RESOLVENDO O SISTEMA COM LU JA FATORADO:
            SOLUTION(:)=FGLOBAL(:)
            CALL DGETRS('N',JJ,1,A_GLOBAL_SLIP_S_LU,JJ,IPIV_A_GLOBAL_SLIP_S,SOLUTION,JJ,INFO)
        
            IF (INFO.NE.0) THEN
                WRITE(*,'(A,I7,A)')'*** PROBLEM IN THE SYSTEM SOLVER, LINE:',INFO,'*****'
                WRITE(*,'(A)')'VERIFY DGETRS MATRIXES IN THE REINF_SOLVER_SLIP SUBROUTINE'
                WRITE(*,'(A)')'*********************************************************'
                READ(*,*)
            ENDIF
        
        ENDIF
                
        ! MUDOU
        U=0.0D0
        T=0.0D0
        
        ! SAVING U AND T FROM SOLUTION __________________________________
        DO I=1,N_COLLOCPOINTS
	        II=0
	        JJ=0
	        IF(DUAL_BEM(I).EQ.'B')THEN  ! POINT AT BOUNDARY
		        DO J=1,N_INTERFACE_POINTS
			        IF(I.EQ.NI(J,1)) THEN
				        II=-1
				        JJ=NI(J,2)
			        ELSE
				        IF(I.EQ.NI(J,2)) THEN
					        II=NI(J,1)
					        JJ=-1
				        ENDIF
			        ENDIF
		        ENDDO	
    
		        IF(II.EQ.-1) THEN
			        U(3*I-2)=SOLUTION(3*I-2)
			        U(3*I-1)=SOLUTION(3*I-1)
			        U(3*I)=SOLUTION(3*I)
    
			        U(3*JJ-2)=SOLUTION(3*I-2)
			        U(3*JJ-1)=SOLUTION(3*I-1)
			        U(3*JJ)=SOLUTION(3*I)
		        ELSE
			        IF(JJ.EQ.-1) THEN
				        T(3*I-2)=-SOLUTION(3*I-2)
				        T(3*I-1)=-SOLUTION(3*I-1)
			            T(3*I)=-SOLUTION(3*I)
    
				        T(3*II-2)=SOLUTION(3*I-2)
				        T(3*II-1)=SOLUTION(3*I-1)
				        T(3*II)=SOLUTION(3*I)
			        ELSE
				        IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					        T(3*I-2)=SOLUTION(3*I-2)
				        ELSE
					        U(3*I-2)=SOLUTION(3*I-2)
				        ENDIF
    
				        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					        T(3*I-1)=SOLUTION(3*I-1)
				        ELSE
					        U(3*I-1)=SOLUTION(3*I-1)
				        ENDIF
    
				        IF (B_CONDITIONS(3*I).EQ.0) THEN
					        T(3*I)=SOLUTION(3*I)
				        ELSE
					        U(3*I)=SOLUTION(3*I)
				        ENDIF
			        ENDIF
		        ENDIF
            ELSE   ! POINT AT CRACK CONTOURN
                !
                ! UNDER CONSTRUCTION
                !
            ENDIF
        ENDDO
        
        DELTA_S_GLOBAL=0.0D0
        P_MEF=0.0D0
        II=3*N_COLLOCPOINTS
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
        !DO I=1,N_NOS_MEF
        !    U_MEF(3*I-2)=SOLUTION(II+3*I-2)
        !    U_MEF(3*I-1)=SOLUTION(II+3*I-1)
        !    U_MEF(3*I)=SOLUTION(II+3*I)
        !    IF (NI_SLIP_MEF(I).EQ.0) THEN ! NÓ COLADO, INCOGNITA = FORÇA
        !        P_MEF(3*I-2)=SOLUTION(JJ+3*I-2)
        !        P_MEF(3*I-1)=SOLUTION(JJ+3*I-1)
        !        P_MEF(3*I)=SOLUTION(JJ+3*I)
        !    ELSE IF (NI_SLIP_MEF(I).EQ.1) THEN ! NÓ DESLOCADO, INCOGNITA = ESCORREGAMENTO
        !        DELTA_S_GLOBAL(3*I-2)=SOLUTION(JJ+3*I-2)
        !        DELTA_S_GLOBAL(3*I-1)=SOLUTION(JJ+3*I-1)
        !        DELTA_S_GLOBAL(3*I)=SOLUTION(JJ+3*I)
        !    ENDIF
        !ENDDO 
         DO I=1,3*N_NOS_MEF
            U_MEF(I)=SOLUTION(II+I)
            IF (NI_SLIP_G(I).EQ.0) THEN ! NÓ COLADO, INCOGNITA = FORÇA
                P_MEF(I)=SOLUTION(JJ+I)
            ELSE IF (NI_SLIP_G(I).EQ.1) THEN ! NÓ DESLOCADO, INCOGNITA = ESCORREGAMENTO
                DELTA_S_GLOBAL(I)=SOLUTION(JJ+I)
            ENDIF
        ENDDO 
        
        ! FORÇA DE REAÇÃO NOS NÓS PUCHADOS
        JJ=3*N_COLLOCPOINTS+6*N_NOS_MEF
        DO I=1,3*N_NODAL_FORCES
            F_APP_DISP(I)=SOLUTION(JJ+I)
            !WRITE(*,*)SOLUTION(JJ+I)
        ENDDO
        
        DEALLOCATE(SOLUTION)
        
        ! Acumulando delta_s
        CALL ROTATE_TO_LOCAL(DELTA_S_GLOBAL,DELTA_S_LOCAL)
        S_ACUM(:)=S_ACUM(:)+DELTA_S_LOCAL(:)        
        
    ENDIF ! ________________________________________________________________________________
    
    ! Aculumando grandezas
    U_ACUM_PL(:)=U_ACUM_PL(:)+U(:)
    P_ACUM_PL(:)=P_ACUM_PL(:)+T(:)
    U_MEF_ACUM(:) = U_MEF_ACUM(:) + U_MEF(:)
    F_APP_DISP_ACUM(:) = F_APP_DISP_ACUM(:) + F_APP_DISP(:) 
        

END SUBROUTINE REINF_SOLVER_SLIP
    
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************      
    

SUBROUTINE REINF_SOLVER_PLAST_2(N_PASSO,N_ITER)
      
    !
    ! SUBROUTINE TO APPLY THE LOAD STEP, SOLVE THE SYSTEM AND FIND DELTA U AND DELTA P FOR ELASTOPLASTIC ANALYSIS
    ! 
    ! N_PASSO: CURRENT NUMBER OF LOAD STEP
    ! N_ITER: CURRENT NUMBER OF ITERATION (INSIDE THE LOAD STEP N_PASSO)
    !

    USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER,INTENT(IN)::N_PASSO,N_ITER
    
    INTEGER::I,J,K,II,JJ,KK,CONT,INFO
    
    REAL*8::AUX1,AUX2,AUX3
    
    REAL*8,DIMENSION(:),ALLOCATABLE::F,SOLUTION
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::AGLOBAL,BGLOBAL
    
    
    IF ((N_PASSO.EQ.1).AND.(N_ITER.EQ.1)) THEN  ! FIRST OF ALL, FIND MATRIXES ______________
        
        CALL REINFORCEMENTS_SOLVER
        
    ELSE IF (N_ITER.EQ.1) THEN ! FIRST ITERATION WITH MATRIXES READY, APPLY U AND T ________
        
        ! Starting variables
        II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
        ALLOCATE(SOLUTION(JJ),F(II))
	    SOLUTION=0.D0
	    F=0.D0
	    CONT=0
        
        DO I=1,N_COLLOCPOINTS
		    II=0
		    JJ=0
		    IF(DUAL_BEM(I).EQ.'B')THEN
		        DO J=1,N_INTERFACE_POINTS
			        IF(I.EQ.NI(J,1)) THEN
				        II=-1
				        JJ=NI(J,2)
			        ELSE
				        IF(I.EQ.NI(J,2)) THEN
					        II=NI(J,1)
					        JJ=-1
				        ENDIF
			        ENDIF
                ENDDO 
                
                IF(II.EQ.-1) THEN
                    !
                ELSE
                    IF(JJ.EQ.-1) THEN
                        !
                    ELSE  
                        IF(B_CONDITIONS(3*I-2).EQ.0) THEN
				            CONT=CONT+1
				            F(CONT)=U((3*I-2))
			            ELSE
				            CONT=CONT+1
				            F(CONT)=T((3*I-2))
			            ENDIF

			            IF(B_CONDITIONS(3*I-1).EQ.0) THEN
				            CONT=CONT+1
				            F(CONT)=U((3*I-1))
			            ELSE
				            CONT=CONT+1
				            F(CONT)=T((3*I-1))
                        ENDIF
                    
                        IF(B_CONDITIONS(3*I).EQ.0) THEN
				            CONT=CONT+1
				            F(CONT)=U((3*I))
			            ELSE
				            CONT=CONT+1
				            F(CONT)=T((3*I))
			            ENDIF
			        ENDIF
		        ENDIF
            ELSE  ! PONTO SOBRE FISSURA
                !
                ! UNDER CONSTRUCTION
                !
            ENDIF
        ENDDO  
        
        ! MULTIPLICANDO B*F
 	    II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
        !CALL DMURRV(JJ,II,BGLOBAL,JJ,II,FGLOBAL,1,JJ,SOLUTION)
        CALL DGEMV('N',JJ,II,1.0D0,B_GLOBAL_EP,JJ,F,1,0.0D0,SOLUTION,1)
        
        F_ALVO(:)=SOLUTION(:)
        
        ! RESOLVENDO O SISTEMA COM LU JA FATORADO:  
        CALL DGETRS('N',JJ,1,A_GLOBAL_EP_LU,JJ,IPIV_A_GLOBAL_EP,SOLUTION,JJ,INFO)
        
        ! SAVING U AND T FROM SOLUTION __________________________________
        DO I=1,N_COLLOCPOINTS
	        II=0
	        JJ=0
	        IF(DUAL_BEM(I).EQ.'B')THEN  ! POINT AT BOUNDARY
		        DO J=1,N_INTERFACE_POINTS
			        IF(I.EQ.NI(J,1)) THEN
				        II=-1
				        JJ=NI(J,2)
			        ELSE
				        IF(I.EQ.NI(J,2)) THEN
					        II=NI(J,1)
					        JJ=-1
				        ENDIF
			        ENDIF
		        ENDDO	
    
		        IF(II.EQ.-1) THEN
			        U(3*I-2)=SOLUTION(3*I-2)
			        U(3*I-1)=SOLUTION(3*I-1)
			        U(3*I)=SOLUTION(3*I)
    
			        U(3*JJ-2)=SOLUTION(3*I-2)
			        U(3*JJ-1)=SOLUTION(3*I-1)
			        U(3*JJ)=SOLUTION(3*I)
		        ELSE
			        IF(JJ.EQ.-1) THEN
				        T(3*I-2)=-SOLUTION(3*I-2)
				        T(3*I-1)=-SOLUTION(3*I-1)
			            T(3*I)=-SOLUTION(3*I)
    
				        T(3*II-2)=SOLUTION(3*I-2)
				        T(3*II-1)=SOLUTION(3*I-1)
				        T(3*II)=SOLUTION(3*I)
			        ELSE
				        IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					        T(3*I-2)=SOLUTION(3*I-2)
				        ELSE
					        U(3*I-2)=SOLUTION(3*I-2)
				        ENDIF
    
				        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					        T(3*I-1)=SOLUTION(3*I-1)
				        ELSE
					        U(3*I-1)=SOLUTION(3*I-1)
				        ENDIF
    
				        IF (B_CONDITIONS(3*I).EQ.0) THEN
					        T(3*I)=SOLUTION(3*I)
				        ELSE
					        U(3*I)=SOLUTION(3*I)
				        ENDIF
			        ENDIF
		        ENDIF
            ELSE   ! POINT AT CRACK CONTOURN
                !
                ! UNDER CONSTRUCTION
                !
            ENDIF
        ENDDO
        II=3*N_COLLOCPOINTS
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
        DO I=1,(3*N_NOS_MEF)
            U_MEF(I)=SOLUTION(II+I)
            P_MEF(I)=SOLUTION(JJ+I)
        ENDDO 
        
        DEALLOCATE(F,SOLUTION)
        
    ELSE ! NON FIRST ITERATION, APPLY F_DES ________________________________________________
        
        ! Starting variables
        II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
        ALLOCATE(SOLUTION(JJ),AGLOBAL(JJ,JJ),BGLOBAL(JJ,II))

        ! MONTANDO AGLOBAL PARA APLICACAO DE F_DES
        AGLOBAL=0.0D0
        BGLOBAL=0.0D0
        SOLUTION=0.0D0
        DO J=1,3*N_COLLOCPOINTS+3*N_NOS_MEF
            DO I=1,3*N_COLLOCPOINTS+3*N_NOS_MEF
            AGLOBAL(I,J)=A_GLOBAL_ELASTIC(I,J)
            ENDDO
        ENDDO
        CONT=0
        II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
        DO J=1,N_NOS_MEF
            IF (NI_PLAST_MEF(J).EQ.0) THEN  ! NÓ SEM FORÇA REAPLICADA
                CONT=CONT+1
                DO I=1,3*N_NOS_MEF
                    AGLOBAL(II,JJ+3*CONT-2)=A_GLOBAL_ELASTIC(II,JJ+3*J-2)
                    AGLOBAL(II,JJ+3*CONT-1)=A_GLOBAL_ELASTIC(II,JJ+3*J-1)
                    AGLOBAL(II,JJ+3*CONT)=A_GLOBAL_ELASTIC(II,JJ+3*J)
                ENDDO
            ELSE  ! NÓ COM FORÇA REAPLICADA
                
            ENDIF
        ENDDO
        
            
        SOLUTION(:)=F_DES(:)
        
        ! RESOLVENDO O SISTEMA COM LU JA FATORADO:
        CALL DGETRS('N',JJ,1,A_GLOBAL_EP_LU,JJ,IPIV_A_GLOBAL_EP,SOLUTION,JJ,INFO)
        
        ! MUDOU
        U=0.0D0
        T=0.0D0
        
        ! SAVING U AND T FROM SOLUTION __________________________________
        DO I=1,N_COLLOCPOINTS
	        II=0
	        JJ=0
	        IF(DUAL_BEM(I).EQ.'B')THEN  ! POINT AT BOUNDARY
		        DO J=1,N_INTERFACE_POINTS
			        IF(I.EQ.NI(J,1)) THEN
				        II=-1
				        JJ=NI(J,2)
			        ELSE
				        IF(I.EQ.NI(J,2)) THEN
					        II=NI(J,1)
					        JJ=-1
				        ENDIF
			        ENDIF
		        ENDDO	
    
		        IF(II.EQ.-1) THEN
			        U(3*I-2)=SOLUTION(3*I-2)
			        U(3*I-1)=SOLUTION(3*I-1)
			        U(3*I)=SOLUTION(3*I)
    
			        U(3*JJ-2)=SOLUTION(3*I-2)
			        U(3*JJ-1)=SOLUTION(3*I-1)
			        U(3*JJ)=SOLUTION(3*I)
		        ELSE
			        IF(JJ.EQ.-1) THEN
				        T(3*I-2)=-SOLUTION(3*I-2)
				        T(3*I-1)=-SOLUTION(3*I-1)
			            T(3*I)=-SOLUTION(3*I)
    
				        T(3*II-2)=SOLUTION(3*I-2)
				        T(3*II-1)=SOLUTION(3*I-1)
				        T(3*II)=SOLUTION(3*I)
			        ELSE
				        IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					        T(3*I-2)=SOLUTION(3*I-2)
				        ELSE
					        U(3*I-2)=SOLUTION(3*I-2)
				        ENDIF
    
				        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					        T(3*I-1)=SOLUTION(3*I-1)
				        ELSE
					        U(3*I-1)=SOLUTION(3*I-1)
				        ENDIF
    
				        IF (B_CONDITIONS(3*I).EQ.0) THEN
					        T(3*I)=SOLUTION(3*I)
				        ELSE
					        U(3*I)=SOLUTION(3*I)
				        ENDIF
			        ENDIF
		        ENDIF
            ELSE   ! POINT AT CRACK CONTOURN
                !
                ! UNDER CONSTRUCTION
                !
            ENDIF
        ENDDO
        II=3*N_COLLOCPOINTS
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
        DO I=1,(3*N_NOS_MEF)
            U_MEF(I)=SOLUTION(II+I)
            P_MEF(I)=SOLUTION(JJ+I)
        ENDDO 
        
        DEALLOCATE(SOLUTION)
        
    ENDIF ! ________________________________________________________________________________
    
    ! Aculumando grandezas
    U_ACUM_PL(:)=U_ACUM_PL(:)+U(:)
    P_ACUM_PL(:)=P_ACUM_PL(:)+T(:)
    U_MEF_ACUM(:) = U_MEF_ACUM(:) + U_MEF(:)
    
    ! Calcular esforços normais nodais
    CALL REINFORCEMENTS_STRESS
    
END SUBROUTINE REINF_SOLVER_PLAST_2
    
