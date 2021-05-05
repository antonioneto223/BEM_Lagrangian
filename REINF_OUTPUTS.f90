SUBROUTINE REINFORCEMENTS_OUTPUT

    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,J,K
    REAL*8::AUX,DIST
        
    OPEN(2,file='Output_data/RESULTS_REINFORCEMENT.txt',status='unknown')

	  WRITE(2,*)
	  WRITE(2,*)'RESULTS AT THE REINFORCEMENTS.BEM/FEM COUPLING'
	  WRITE(2,*)
      WRITE(2,*)'NUMERICAL METHOD FOR REINFORCEMENTS:'
      WRITE(2,*)REINFORCEMENTS_METHOD
      WRITE(2,*)
	  ! ESCREVENDO COORDENADAS NODAIS DO REFORÇO
	  !IF (REINFORCEMENTS_METHOD.EQ.'FEM') THEN
	  !  WRITE(2,*)'COORDINATES OF NODAL POINTS'
	  !  WRITE(2,*)'   NODE          X                Y          Z'
   !     DO I=1,N_NOS_MEF
   !         WRITE(2,'(I3,3E14.6)')i,COORDNOS_MEF(i,1),COORDNOS_MEF(i,2),COORDNOS_MEF(i,3)       
   !     ENDDO
   !   ELSE
        WRITE(2,*)'COORDINATES OF NODAL AND COLLOCATION POINTS'
	    WRITE(2,*)'                     NODAL POINT                     COLLOCATION POINT'
        WRITE(2,'(A)')'   NODE    X             Y          Z                 X           Y             Z'
        DO I=1,N_NOS_MEF
            WRITE(2,'(i5,6E14.6)')i,COORDNOS_MEF(i,1),COORDNOS_MEF(i,2),COORDNOS_MEF(i,3),COORDPCOLOC_MEF(I,1),COORDPCOLOC_MEF(i,2),COORDPCOLOC_MEF(i,3)
        ENDDO
      !ENDIF
	  
	  ! ESCREVENDO CONECTIVIDADE DOS ELEMENTOS DO REFORÇO
	  WRITE(2,*)''
	  WRITE(2,*)'CONECTIVITY OF REINFORCEMENTS ELEMENTS'
	  DO I=1,N_ELEMENTOS_MEF
		    WRITE(2,*)'   ELEMENT       DOMAIN	      NODES'
		    DO J=1,(ORDEM_MEF(I)+1)
			    WRITE(2,'(I4,I3,I6)')I,ND_MEF(I),CONECTI_MEF(I,J)
		    ENDDO
	  ENDDO
	  
	  ! ESCREVENDO PROPRIEDADES DO REFORÇO
	  WRITE(2,*)''	  
	  WRITE(2,*)'MATERIAL PROPERTIES ALONG THE REINFORCEMENTS'
	  WRITE(2,*)'ELEMENT         DOMAIN		      YOUNG S MODULUS   CROSS SECTIONAL AREA'
	  DO I=1,N_ELEMENTOS_MEF
      		WRITE(2,'(I4,I3,2F16.6)')I,ND_MEF(I),MP_MEF(I,1),MP_MEF(I,2)
	  ENDDO
        
      ! RESULTADOS DE DESLOCAMENTO E FORÇA DE SUPERFICIE
      IF (SLIP_REINF) THEN
          WRITE(2,*)''
          WRITE(2,*)'VALUES CALCULATED AND PRESCRIBED ALONG THE REINFORCEMENTS + SLIP NODAL VALUES'
          WRITE(2,'(A)')'  NODE         UX               UY              UZ                   PX               PY                PZ                 SX                SY               SZ'
          DO I=1,N_NOS_MEF
            WRITE(2,'(i5,9E18.6)')I,U_MEF(3*I-2)-S_ACUM(3*I-2),U_MEF(3*I-1)-S_ACUM(3*I-1),U_MEF(3*I)-S_ACUM(3*I),&
                P_MEF(3*I-2),P_MEF(3*I-1),P_MEF(3*I),S_ACUM(3*I-2),S_ACUM(3*I-1),S_ACUM(3*I)
          ENDDO      
      ELSE
          WRITE(2,*)''
          WRITE(2,*)'VALUES CALCULATED AND PRESCRIBED ALONG THE REINFORCEMENTS'
          WRITE(2,'(A)')'     NODE         UX               UY              UZ                P_X           P_Y         P_Z'
          DO I=1,N_NOS_MEF
            WRITE(2,'(i5,6E18.6)')I,U_MEF(3*I-2),U_MEF(3*I-1),U_MEF(3*I),P_MEF(3*I-2),P_MEF(3*I-1),P_MEF(3*I)
          ENDDO
      ENDIF
      
      ! RESULTADOS DE TENSAO NORMAL NOS ELEMENTOS
      IF (EP_ANALYSIS_REINF) THEN  
          WRITE(2,*)
          WRITE(2,*)'AXIAL FORCE AND PLASTIC STRAIN AT NODES'
          WRITE(2,'(A)')'   NODE      DIST_S        AXIAL_FORCE         PLASTIC_STRAIN'
          DO I=1,N_NOS_MEF
              DIST=0.0D0
              DO J=1,3
                  DIST=DIST+(COORDNOS_MEF(I,J)-COORDPCOLOC_MEF(I,J))**2.0D0
              ENDDO
              DIST=DSQRT(DIST)
              IF (DIST.GT.TOLER) THEN
                  AUX=DIST
              ELSE
                  AUX=AUX+DSQRT((COORDPCOLOC_MEF(I,1)-COORDPCOLOC_MEF(I-1,1))**2.0D0+(COORDPCOLOC_MEF(I,2)-COORDPCOLOC_MEF(I-1,2))**2.0D0&
                  +(COORDPCOLOC_MEF(I,3)-COORDPCOLOC_MEF(I-1,3))**2.0D0)
              ENDIF  
              
              WRITE(2,'(I8,F14.4,2E18.6)')I,AUX,E_NORMAL_NOD(I),PL_STRAIN_NOD(I)
          ENDDO
      ELSE 
          WRITE(2,*)
          WRITE(2,*)'AXIAL FORCE AT REINFORCEMENTS ELEMENTS NODES'
          WRITE(2,'(A)')'   NODE      DIST_S        AXIAL_FORCE'
          DO I=1,N_NOS_MEF
              DIST=0.0D0
              DO J=1,3
                  DIST=DIST+(COORDNOS_MEF(I,J)-COORDPCOLOC_MEF(I,J))**2.0D0
              ENDDO
              DIST=DSQRT(DIST)
              IF (DIST.GT.TOLER) THEN
                  AUX=DIST
              ELSE
                  AUX=AUX+DSQRT((COORDPCOLOC_MEF(I,1)-COORDPCOLOC_MEF(I-1,1))**2.0D0+(COORDPCOLOC_MEF(I,2)-COORDPCOLOC_MEF(I-1,2))**2.0D0&
                  +(COORDPCOLOC_MEF(I,3)-COORDPCOLOC_MEF(I-1,3))**2.0D0)
              ENDIF  
              
              WRITE(2,'(I8,F14.4,F14.6)')I,AUX,E_NORMAL_NOD(I)
          ENDDO
          !WRITE(2,*)
          !WRITE(2,*)'AXIAL FORCE AT REINFORCEMENTS ELEMENTS NODES'
          !WRITE(2,'(A)')'TOTAL_NODE   ELE  LOCAL_NODE   AXIAL FORCE'
          !K=0
          !DO I=1,N_ELEMENTOS_MEF
          !    DO J=1,ORDEM_MEF(I)+1
          !      K=K+1
          !      WRITE(2,'(I8,I6,I8,E18.6)')K,I,J,E_NORMAL(I,J)
          !    ENDDO
          !ENDDO
      ENDIF

END SUBROUTINE
    
    
        
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  
    

SUBROUTINE NODE_HISTORIC_PRINT(STEP)
    
    !
    ! SUBROUTINE TO PLOT NODE HISTORIC DISPLACEMENTS
    ! STEP: number of current load step of non-linear process
    !

    USE REINFORCEMENTS
    USE ISOPARAMETRIC_MESH
    
    IMPLICIT NONE
    
    INTEGER,INTENT(IN)::STEP
    INTEGER::I,J,K,N_NO
    
    IF (EP_ANALYSIS_REINF) THEN             ! SAIDA PARA ELASTOPLASTICIDADE
        
        IF (STEP.EQ.1) THEN
            OPEN(18,FILE='Output_data/NODE_HISTORY.TXT',STATUS='UNKNOWN')
            WRITE(18,*)'N NODE:',N_NODE_HIST
            WRITE(18,*)'COORD X:',COORD_COLLOCPOINTS(N_NODE_HIST,1)
            WRITE(18,*)'COORD Y:',COORD_COLLOCPOINTS(N_NODE_HIST,2)
            WRITE(18,*)'COORD Z:',COORD_COLLOCPOINTS(N_NODE_HIST,3)
            WRITE(18,*)''
            WRITE(18,'(A)')' LOAD_STEP  UX		      UY		      UZ'
        ENDIF
    
        WRITE(18,'(I4,3E15.5)')STEP,U_ACUM_PL(3*N_NODE_HIST-2),U_ACUM_PL(3*N_NODE_HIST-1),U_ACUM_PL(3*N_NODE_HIST)
    
        IF (STEP.EQ.N_PASSOS) CLOSE (18)
        
    ENDIF
    
    IF (SLIP_REINF) THEN                    ! SAIDA PARA ESCORREGAMENTO
        
        IF (STEP.EQ.1) THEN
            IF (NODE_HIST.EQ.'Y') THEN
                OPEN(19,FILE='Output_data/NODE_HISTORY_SLIP.TXT',STATUS='UNKNOWN')
                WRITE(19,*)'N NODE:',N_NODE_HIST
                WRITE(19,*)'COORD X:',COORDPCOLOC_MEF(N_NODE_HIST,1)
                WRITE(19,*)'COORD Y:',COORDPCOLOC_MEF(N_NODE_HIST,2)
                WRITE(19,*)'COORD Z:',COORDPCOLOC_MEF(N_NODE_HIST,3)
                WRITE(19,*)''
                WRITE(19,'(A)')'LOAD_STEP  UX				UY			UZ      		P_LOCAL     S_LOCAL'
            ENDIF
            IF (N_NODAL_FORCES.NE.0) THEN
                OPEN(20,FILE='Output_data/NODE_HISTORY_REACTIONS.TXT',STATUS='UNKNOWN')
                WRITE(20,*)''
                WRITE(20,'(A)',ADVANCE='NO')'STEP'
                DO I=1,N_NODAL_FORCES
                    WRITE(20,'(A,I1,A,I1,A,I1)',ADVANCE='NO')'          FX_',I,'    		FY_',I,'    		FZ_',I
                ENDDO
                WRITE(20,*)''
            ENDIF
        ENDIF
            
        IF (NODE_HIST.EQ.'Y') THEN
            K=0
            DO I=1,N_ELEMENTOS_MEF
                DO J=1,ORDEM_MEF(I)+1
                    K=K+1
                    IF (CONECTI_MEF(I,J).EQ.N_NODE_HIST) N_NO=K                    
                ENDDO
            ENDDO
            WRITE(19,'(I4,5E15.5)')STEP,U_MEF_ACUM(3*N_NODE_HIST-2),U_MEF_ACUM(3*N_NODE_HIST-1),U_MEF_ACUM(3*N_NODE_HIST),&
                P_MEF_ACUM_LOCAL(N_NO),S_ACUM(N_NO)  
        ENDIF
        
        IF (N_NODAL_FORCES.NE.0) THEN
            WRITE(20,'(I4)',ADVANCE='NO')STEP
            DO I=1,N_NODAL_FORCES
                WRITE(20,'(3E15.5)',ADVANCE='NO')F_APP_DISP_ACUM(3*I-2),F_APP_DISP_ACUM(3*I-1),F_APP_DISP_ACUM(3*I)
            ENDDO
            WRITE(20,'(A)')''
        ENDIF
    
        IF (STEP.EQ.N_PASSOS) THEN
            IF (NODE_HIST.EQ.'Y')   CLOSE (19)
            IF (N_NODAL_FORCES.NE.0)    CLOSE(20)
        ENDIF
        
    ENDIF
    
END SUBROUTINE
    
        
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  
    

SUBROUTINE TXT_OUTPUT_BOND_SLIP(STEP,ITRS,ERR)

    !
    !   STEP: NUMBER OF CURRENT LOAD STEP, ITRS: NUMBER OF ITERATIONS FOR THIS STEP, ERR: ERROR NUMBER IN THE FINAL ITERATION
    !
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER,INTENT(IN)::STEP,ITRS
    REAL*8,INTENT(IN)::ERR
    INTEGER::I,J,K
    
    ! Escreve resultados no txt
    WRITE(17,*)''
    WRITE(17,'(a,i3,a)')'____________________|| PASSO DE CARGA:',STEP,'||____________________'
    WRITE(17,*)''
    WRITE(17,*)'ERRO AO FINAL DO PASSO:',ERR
    WRITE(17,*)'NUMERO DE ITERACOES:', ITRS
    WRITE(17,*)''
    WRITE(17,*)'GRANDEZAS ATUAIS AO FINAL DO PASSO NOS NÓS DOS ELEMENTOS:'
    WRITE(17,'(A)')'NO_TOTAL   NO_GLOBAL    ELEM   NO_LOCAL      P_MEF_ACUM        S_ACUM'
    K=0
    DO I=1,N_ELEMENTOS_MEF
        DO J=1,ORDEM_MEF(I)+1
            K=K+1
            WRITE(17,'(I8,I8,I6,I6,E18.6,E18.6)')K,CONECTI_MEF(I,J),I,J,P_MEF_ACUM_LOCAL(K),S_ACUM(K)
        ENDDO
    ENDDO
    WRITE(17,*)''
    WRITE(17,*)'REACOES NO NÓS PUCHADOS:'
    WRITE(17,*)'NO_GLOBAL   FX      FY      FZ'    
    DO I=1,N_NODAL_FORCES
        WRITE(17,'(I8,3E18.6)')NODAL_FORCES_PLACE(I,1),F_APP_DISP_ACUM(3*I-2),F_APP_DISP_ACUM(3*I-1),F_APP_DISP_ACUM(3*I)
    ENDDO        
    WRITE(17,*)''
        
END SUBROUTINE