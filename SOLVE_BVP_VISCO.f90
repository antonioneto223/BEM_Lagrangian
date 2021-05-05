SUBROUTINE SOLVE_BVP_VISCO

    !
    !   SUBROUTINE TO ASSEMBLE GLOBAL MATRIXES AND SOLVE LINEAR SYSTEM FOR VISCOELASTIC PROBLEMS (REINFORCED OR NOT)
    !
    
    USE ISOPARAMETRIC_MESH
	USE CRACKS_DUAL_BEM
    USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE VISCO
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,J,K,II,JJ,CONT,INFO,NUM_LEI
    
    REAL*8::K1,K2,K3,K4,GAMA
    REAL*8,DIMENSION(:,:),POINTER::A,A_LINHA,B,B_LINHA
    REAL*8,DIMENSION(:,:),ALLOCATABLE::AGLOBAL
    REAL*8,DIMENSION(:),ALLOCATABLE::SOLUTION,F,FGLOBAL
    
    IF(EXISTANCE_REINFORCEMENTS.EQ.'N')THEN  ! ______________________ NON REINFORCED MEDIA __________________________
        

        IF (.NOT.ALLOCATED(A_GLOBAL_LU)) THEN
         
            ! Starting variables
            II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            ALLOCATE(AGLOBAL(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS),SOLUTION(3*N_COLLOCPOINTS),F(II),&
            B_GLOBAL(3*N_COLLOCPOINTS,II))
            AGLOBAL=0.D0
            B_GLOBAL=0.0D0
	        B=>B_GLOBAL(1:3*N_COLLOCPOINTS,1:II)
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
			            AGLOBAL(:,(3*I-2))=AGLOBAL(:,(3*I-2))+HB(:,(3*I-2))
			            AGLOBAL(:,3*I-1)=AGLOBAL(:,3*I-1)+HB(:,3*I-1)
			            AGLOBAL(:,3*I)=AGLOBAL(:,3*I)+HB(:,3*I)
        !
			            AGLOBAL(:,(3*I-2))=AGLOBAL(:,(3*I-2))+HB(:,(3*JJ-2))
			            AGLOBAL(:,3*I-1)=AGLOBAL(:,3*I-1)+HB(:,3*JJ-1)
			            AGLOBAL(:,3*I)=AGLOBAL(:,3*I)+HB(:,3*JJ)
                    ELSE
                        IF(JJ.EQ.-1) THEN
				            AGLOBAL(:,(3*I-2))=AGLOBAL(:,(3*I-2))+GB(:,(3*I-2))
				            AGLOBAL(:,3*I-1)=AGLOBAL(:,3*I-1)+GB(:,3*I-1)
				            AGLOBAL(:,3*I)=AGLOBAL(:,3*I)+GB(:,3*I)
            !
				            AGLOBAL(:,(3*I-2))=AGLOBAL(:,(3*I-2))-GB(:,(3*II-2))
				            AGLOBAL(:,3*I-1)=AGLOBAL(:,3*I-1)-GB(:,3*II-1)
				            AGLOBAL(:,3*I)=AGLOBAL(:,3*I)-GB(:,3*II)
                        ELSE  
                            IF(B_CONDITIONS(3*I-2).EQ.0) THEN
					            AGLOBAL(:,(3*I-2))=-GB(:,(3*I-2))
        !
					            CONT=CONT+1
					            B(:,CONT)=-HB(:,(3*I-2))
					            F(CONT)=U((3*I-2))
				            ELSE
        					    AGLOBAL(:,(3*I-2))=HB(:,(3*I-2))	
    !
					            CONT=CONT+1
					            B(:,CONT)=GB(:,(3*I-2))
					            F(CONT)=T((3*I-2))
				            ENDIF
                            IF(B_CONDITIONS(3*I-1).EQ.0) THEN
					            AGLOBAL(:,3*I-1)=-GB(:,3*I-1)
        !
					            CONT=CONT+1
					            B(:,CONT)=-HB(:,3*I-1)
					            F(CONT)=U(3*I-1)
				            ELSE
					            AGLOBAL(:,3*I-1)=HB(:,3*I-1)	
        !
					            CONT=CONT+1
					            B(:,CONT)=GB(:,3*I-1)
					            F(CONT)=T(3*I-1)
				            ENDIF
        !
				            IF(B_CONDITIONS(3*I).EQ.0) THEN
					            AGLOBAL(:,3*I)=-GB(:,3*I)
        !
					            CONT=CONT+1
				    	        B(:,CONT)=-HB(:,3*I)
					            F(CONT)=U(3*I)
				            ELSE
					            AGLOBAL(:,3*I)=HB(:,3*I)	
        !
					            CONT=CONT+1
					            B(:,CONT)=GB(:,3*I)
					            F(CONT)=T(3*I)
				            ENDIF
			            ENDIF
		            ENDIF
		        ELSE
		            DO J=1,N_CRACK_POINTS
			            IF(I.EQ.NC(J,1)) THEN
				            II=-1
				            JJ=NC(J,2)
			            ELSE
				            IF(I.EQ.NC(J,2)) THEN
					            II=NC(J,1)
					            JJ=-1
				            ENDIF
			            ENDIF
			        ENDDO        
		            IF(II.EQ.-1) THEN
		    	        AGLOBAL(:,(3*I-2))=AGLOBAL(:,(3*I-2))+HB(:,(3*I-2))
			            AGLOBAL(:,(3*I-2))=AGLOBAL(:,(3*I-2))-HB(:,(3*JJ-2))			    		    
        !			    
		                IF(B_CONDITIONS(3*I-1).EQ.0) THEN
			                AGLOBAL(:,(3*I-1))=-GB(:,(3*I-1))
        !
					        CONT=CONT+1
					        B(:,CONT)=-HB(:,(3*I-1))
					        F(CONT)=U((3*I-1))
				        ELSE
					        AGLOBAL(:,(3*I-1))=HB(:,(3*I-1))	
        !
					        CONT=CONT+1
					        B(:,CONT)=GB(:,(3*I-1))
					        F(CONT)=T((3*I-1))
				        ENDIF
        !
				        IF(B_CONDITIONS(3*I).EQ.0) THEN
					        AGLOBAL(:,3*I)=-GB(:,3*I)
        !
					        CONT=CONT+1
					        B(:,CONT)=-HB(:,3*I)
					        F(CONT)=U(3*I)
				        ELSE
					        AGLOBAL(:,3*I)=HB(:,3*I)	
        !
					        CONT=CONT+1
					        B(:,CONT)=GB(:,3*I)
					        F(CONT)=T(3*I)
				        ENDIF
        !		    
		            ELSE
			            IF(JJ.EQ.-1) THEN
				            AGLOBAL(:,3*I-2)=AGLOBAL(:,3*I-2)-GB(:,3*I-2)
				            AGLOBAL(:,3*I-2)=AGLOBAL(:,3*I-2)-GB(:,3*II-2)
        !				    			    		    
		                    IF(B_CONDITIONS(3*I-1).EQ.0) THEN
			                    AGLOBAL(:,(3*I-1))=-GB(:,(3*I-1))
        !
					            CONT=CONT+1
					            B(:,CONT)=-HB(:,(3*I-1))
				 	            F(CONT)=U((3*I-1))
				            ELSE
					            AGLOBAL(:,(3*I-1))=HB(:,(3*I-1))	
        !
					            CONT=CONT+1
					            B(:,CONT)=GB(:,(3*I-1))
					            F(CONT)=T((3*I-1))
				            ENDIF
        !
				            IF(B_CONDITIONS(3*I).EQ.0) THEN
					            AGLOBAL(:,3*I)=-GB(:,3*I)
        !
					            CONT=CONT+1
					            B(:,CONT)=-HB(:,3*I)
					            F(CONT)=U(3*I)
				            ELSE
					            AGLOBAL(:,3*I)=HB(:,3*I)	
        !
					            CONT=CONT+1
					            B(:,CONT)=GB(:,3*I)
					            F(CONT)=T(3*I)
				            ENDIF			    
        !			    
			            ELSE
				            IF(B_CONDITIONS(3*I-2).EQ.0) THEN
					            AGLOBAL(:,(3*I-2))=-GB(:,(3*I-2))
        !
					            CONT=CONT+1
					            B(:,CONT)=-HB(:,(3*I-2))
					            F(CONT)=U((3*I-2))
				            ELSE
					            AGLOBAL(:,(3*I-2))=HB(:,(3*I-2))	
        !
					            CONT=CONT+1
					            B(:,CONT)=GB(:,(3*I-2))
					            F(CONT)=T((3*I-2))
				            ENDIF
        !
				            IF(B_CONDITIONS(3*I-1).EQ.0) THEN
					            AGLOBAL(:,3*I-1)=-GB(:,3*I-1)
        !
					            CONT=CONT+1
					            B(:,CONT)=-HB(:,3*I-1)
					            F(CONT)=U(3*I-1)
				            ELSE
					            AGLOBAL(:,3*I-1)=HB(:,3*I-1)	
        !
					            CONT=CONT+1
					            B(:,CONT)=GB(:,3*I-1)
					            F(CONT)=T(3*I-1)
				            ENDIF
        !
				            IF(B_CONDITIONS(3*I).EQ.0) THEN
					            AGLOBAL(:,3*I)=-GB(:,3*I)
        !
					            CONT=CONT+1
				    	        B(:,CONT)=-HB(:,3*I)
					            F(CONT)=U(3*I)
				            ELSE
					            AGLOBAL(:,3*I)=HB(:,3*I)	
        !
					            CONT=CONT+1
					            B(:,CONT)=GB(:,3*I)
					            F(CONT)=T(3*I)
				            ENDIF
			            ENDIF
		            ENDIF
		        ENDIF
            ENDDO  ! End of collocation points loop
                    
            ! Somando vetores livres (a direita da eq.)
            II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            !CALL DMURRV(3*N_COLLOCPOINTS,II,B_GLOBAL,3*N_COLLOCPOINTS,II,F,1,3*N_COLLOCPOINTS,SOLUTION)
            CALL DGEMV('N',3*N_COLLOCPOINTS,II,1.0D0,B_GLOBAL,3*N_COLLOCPOINTS,F,1,0.0D0,SOLUTION,1)
	        DEALLOCATE(F)
            ALLOCATE(F(3*N_COLLOCPOINTS))
	        F=SOLUTION
        
           !CALL DMURRV(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS,HBB,3*N_COLLOCPOINTS,3*N_COLLOCPOINTS,U_ANTERIOR,1,3*N_COLLOCPOINTS,SOLUTION)
            CALL DGEMV('N',3*N_COLLOCPOINTS,3*N_COLLOCPOINTS,1.0D0,HBB,3*N_COLLOCPOINTS,U_ANTERIOR,1,0.0D0,SOLUTION,1)
            F(:)=F(:)+SOLUTION(:)
            !CALL DMURRV(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS,GBB,3*N_COLLOCPOINTS,3*N_COLLOCPOINTS,T_ANTERIOR,1,3*N_COLLOCPOINTS,SOLUTION)
            CALL DGEMV('N',3*N_COLLOCPOINTS,3*N_COLLOCPOINTS,1.0D0,GBB,3*N_COLLOCPOINTS,T_ANTERIOR,1,0.0D0,SOLUTION,1)
            F(:)=F(:)+SOLUTION(:)
                
            ! INVERTENDO A MATRIZ AGLOBAL E GUARDANDO EM A_GLOBAL_LU (FATORACAO LU):
            JJ=3*N_COLLOCPOINTS
            ALLOCATE(A_GLOBAL_LU(JJ,JJ),IPVT_AG(JJ))
            !CALL DLINRG (JJ,AGLOBAL,JJ,AGLOBAL_INV,JJ) 
            !CALL DLFTRG(A,A_GLOBAL_LU,IPVT_AG) 
            A_GLOBAL_LU=AGLOBAL
            CALL DGETRF(JJ,JJ,A_GLOBAL_LU,JJ,IPVT_AG,INFO)
            
            IF (INFO.NE.0) THEN
                WRITE(*,*)'*************** CAUTION **********************'
                WRITE(*,*)'******* INFO NOT EQUAL ZERO IN DGETRF ********'
                WRITE(*,*)'VERIFY SOLVE_BVP_VISCO WITHOUT REINFORCEMENTS'
                READ(*,*)
            ENDIF
            
            DEALLOCATE(AGLOBAL)
            
        ELSE ! MATRIZ A_GLOBAL_LU JÁ CALCLULADA, SOMENTE MONTA VETOR DO LADO DIREITO
            
            ! Starting variables
            II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            ALLOCATE(SOLUTION(3*N_COLLOCPOINTS),F(II))
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
					            F(CONT)=U(3*I-1)
				            ELSE
					            CONT=CONT+1
					            F(CONT)=T(3*I-1)
				            ENDIF
        !
				            IF(B_CONDITIONS(3*I).EQ.0) THEN
					            CONT=CONT+1
					            F(CONT)=U(3*I)
				            ELSE
					            CONT=CONT+1
					            F(CONT)=T(3*I)
				            ENDIF
			            ENDIF
		            ENDIF
		        ELSE
		            DO J=1,N_CRACK_POINTS
			            IF(I.EQ.NC(J,1)) THEN
				            II=-1
				            JJ=NC(J,2)
			            ELSE
				            IF(I.EQ.NC(J,2)) THEN
					            II=NC(J,1)
					            JJ=-1
				            ENDIF
			            ENDIF
			        ENDDO        
		            IF(II.EQ.-1) THEN	    		    
                        !			    
		                IF(B_CONDITIONS(3*I-1).EQ.0) THEN
					        CONT=CONT+1
					        F(CONT)=U((3*I-1))
				        ELSE
					        CONT=CONT+1
					        F(CONT)=T((3*I-1))
				        ENDIF
        !
				        IF(B_CONDITIONS(3*I).EQ.0) THEN
					        CONT=CONT+1
					        F(CONT)=U(3*I)
				        ELSE
					        CONT=CONT+1
					        F(CONT)=T(3*I)
				        ENDIF
        !		    
		            ELSE
			            IF(JJ.EQ.-1) THEN
                            !   			    		    
		                    IF(B_CONDITIONS(3*I-1).EQ.0) THEN
					            CONT=CONT+1
				 	            F(CONT)=U((3*I-1))
				            ELSE
					            CONT=CONT+1
					            F(CONT)=T((3*I-1))
				            ENDIF
        !
				            IF(B_CONDITIONS(3*I).EQ.0) THEN
					            CONT=CONT+1
					            F(CONT)=U(3*I)
				            ELSE
					            CONT=CONT+1
					            F(CONT)=T(3*I)
				            ENDIF			    
        !			    
			            ELSE
				            IF(B_CONDITIONS(3*I-2).EQ.0) THEN
					            CONT=CONT+1
					            F(CONT)=U((3*I-2))
				            ELSE
					            CONT=CONT+1
					            F(CONT)=T((3*I-2))
				            ENDIF
        !
				            IF(B_CONDITIONS(3*I-1).EQ.0) THEN
					            CONT=CONT+1
					            F(CONT)=U(3*I-1)
				            ELSE
					            CONT=CONT+1
					            F(CONT)=T(3*I-1)
				            ENDIF
        !
				            IF(B_CONDITIONS(3*I).EQ.0) THEN
					            CONT=CONT+1
					            F(CONT)=U(3*I)
				            ELSE
					            CONT=CONT+1
					            F(CONT)=T(3*I)
				            ENDIF
			            ENDIF
		            ENDIF
		        ENDIF
            ENDDO  ! End of collocation points loop
                    
            ! Somando vetores livres (a direita da eq.)
            II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            !CALL DMURRV(3*N_COLLOCPOINTS,II,B_GLOBAL,3*N_COLLOCPOINTS,II,F,1,3*N_COLLOCPOINTS,SOLUTION)
            CALL DGEMV('N',3*N_COLLOCPOINTS,II,1.0D0,B_GLOBAL,3*N_COLLOCPOINTS,F,1,0.0D0,SOLUTION,1)
	        DEALLOCATE(F)
            ALLOCATE(F(3*N_COLLOCPOINTS))
	        F=SOLUTION
        
            !CALL DMURRV(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS,HBB,3*N_COLLOCPOINTS,3*N_COLLOCPOINTS,U_ANTERIOR,1,3*N_COLLOCPOINTS,SOLUTION)
            CALL DGEMV('N',3*N_COLLOCPOINTS,3*N_COLLOCPOINTS,1.0D0,HBB,3*N_COLLOCPOINTS,U_ANTERIOR,1,0.0D0,SOLUTION,1)
            F(:)=F(:)+SOLUTION(:)
            !CALL DMURRV(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS,GBB,3*N_COLLOCPOINTS,3*N_COLLOCPOINTS,T_ANTERIOR,1,3*N_COLLOCPOINTS,SOLUTION)
            CALL DGEMV('N',3*N_COLLOCPOINTS,3*N_COLLOCPOINTS,1.0D0,GBB,3*N_COLLOCPOINTS,T_ANTERIOR,1,0.0D0,SOLUTION,1)
            F(:)=F(:)+SOLUTION(:)
                            
        ENDIF
                
        ! RESOLVENDO O SISTEMA COM LU FATORADO:  
        JJ=3*N_COLLOCPOINTS
        !CALL DMURRV(JJ,JJ,AGLOBAL_INV,JJ,JJ,F,1,JJ,SOLUCAO)  
        !CALL DLFSRG(JJ,A_GLOBAL_LU,JJ,IPVT_AG,F,1,SOLUTION
        SOLUTION=F
        CALL DGETRS('N',JJ,1,A_GLOBAL_LU,JJ,IPVT_AG,SOLUTION,JJ,INFO)
        
        ! Salvando valores de U e T encontrados
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
    !
		        IF(II.EQ.-1) THEN
			        U(3*I-2)=SOLUTION(3*I-2)
			        U(3*I-1)=SOLUTION(3*I-1)
			        U(3*I)=SOLUTION(3*I)
    !
			        U(3*JJ-2)=SOLUTION(3*I-2)
			        U(3*JJ-1)=SOLUTION(3*I-1)
			        U(3*JJ)=SOLUTION(3*I)
		        ELSE
			        IF(JJ.EQ.-1) THEN
				        T(3*I-2)=-SOLUTION(3*I-2)
				        T(3*I-1)=-SOLUTION(3*I-1)
			            T(3*I)=-SOLUTION(3*I)
    !
				        T(3*II-2)=SOLUTION(3*I-2)
				        T(3*II-1)=SOLUTION(3*I-1)
				        T(3*II)=SOLUTION(3*I)
			        ELSE
				        IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					        T(3*I-2)=SOLUTION(3*I-2)
				        ELSE
					        U(3*I-2)=SOLUTION(3*I-2)
				        ENDIF
    !
				        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					        T(3*I-1)=SOLUTION(3*I-1)
				        ELSE
					        U(3*I-1)=SOLUTION(3*I-1)
				        ENDIF
    !
				        IF (B_CONDITIONS(3*I).EQ.0) THEN
					        T(3*I)=SOLUTION(3*I)
				        ELSE
					        U(3*I)=SOLUTION(3*I)
				        ENDIF
			        ENDIF
		        ENDIF
		    ELSE
		        DO J=1,N_CRACK_POINTS
			        IF(I.EQ.NC(J,1)) THEN
				        II=-1
				        JJ=NC(J,2)
			        ELSE
				        IF(I.EQ.NC(J,2)) THEN
					        II=NC(J,1)
					        JJ=-1
				        ENDIF
			        ENDIF
		        ENDDO	
    !
		        IF(II.EQ.-1) THEN
			        U(3*I-2)=SOLUTION(3*I-2)
			        U(3*JJ-2)=-SOLUTION(3*I-2)
    !			    		    
			        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					    T(3*I-1)=SOLUTION(3*I-1)
				    ELSE
					    U(3*I-1)=SOLUTION(3*I-1)
				    ENDIF
    !
				    IF (B_CONDITIONS(3*I).EQ.0) THEN
					    T(3*I)=SOLUTION(3*I)
				    ELSE
					    U(3*I)=SOLUTION(3*I)
				    ENDIF
    !				    		    
		        ELSE
			        IF(JJ.EQ.-1) THEN
			            T(3*I-2)=SOLUTION(3*I-2)
				        T(3*II-2)=SOLUTION(3*I-2)
    !				    			    
			            IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					        T(3*I-1)=SOLUTION(3*I-1)
				        ELSE
					        U(3*I-1)=SOLUTION(3*I-1)
				        ENDIF
    !
				        IF (B_CONDITIONS(3*I).EQ.0) THEN
					        T(3*I)=SOLUTION(3*I)
				        ELSE
					        U(3*I)=SOLUTION(3*I)
				        ENDIF			    
    !			    
			        ELSE
				        IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					        T(3*I-2)=SOLUTION(3*I-2)
				        ELSE
					        U(3*I-2)=SOLUTION(3*I-2)
				        ENDIF
    !
				        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					        T(3*I-1)=SOLUTION(3*I-1)
				        ELSE
					        U(3*I-1)=SOLUTION(3*I-1)
				        ENDIF
    !
				        IF (B_CONDITIONS(3*I).EQ.0) THEN
					        T(3*I)=SOLUTION(3*I)
				        ELSE
					        U(3*I)=SOLUTION(3*I)
				        ENDIF
			        ENDIF
		        ENDIF			
		    ENDIF
	    ENDDO  
        
        DEALLOCATE(SOLUTION,F)
        
             !--------------------------------------------------------------------------------------------------------------------
    ELSE     !------------- REINFORCED VISCOELASTIC MEDIA ------------------------------------------------------------------------
             !--------------------------------------------------------------------------------------------------------------------
        
        IF (.NOT.ALLOCATED(A_GLOBAL_LU)) THEN
            
            ! Starting variables
            II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
            !ALLOCATE(A(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS),A_LINHA(3*N_NOS_MEF,3*N_COLLOCPOINTS),&
            !    SOLUTION(JJ),B(3*N_COLLOCPOINTS,II),B_LINHA(3*N_NOS_MEF,II),F(II),&
            !    AGLOBAL(JJ,JJ),B_GLOBAL(JJ,II+6*N_NOS_MEF+6*N_COLLOCPOINTS),&
            !    FGLOBAL(II+6*N_NOS_MEF+6*N_COLLOCPOINTS))
            ALLOCATE(SOLUTION(JJ),F(II),FGLOBAL(II+6*N_NOS_MEF+6*N_COLLOCPOINTS))       
            FGLOBAL=0.0D0
	        SOLUTION=0.D0
	        F=0.D0
	        CONT=0
            
            A=>A_GLOBAL_ELASTIC(1:3*N_COLLOCPOINTS,1:3*N_COLLOCPOINTS)
            A_LINHA=>A_GLOBAL_ELASTIC(3*N_COLLOCPOINTS+1:3*N_COLLOCPOINTS+3*N_NOS_MEF,1:3*N_COLLOCPOINTS)
            
            B=>B_GLOBAL(1:3*N_COLLOCPOINTS,1:II)
            B_LINHA=>B_GLOBAL(3*N_COLLOCPOINTS+1:3*N_COLLOCPOINTS+3*N_NOS_MEF,1:II)

        
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
			            A(:,(3*I-2))=A(:,(3*I-2))+HB(:,(3*I-2))
			            A(:,3*I-1)=A(:,3*I-1)+HB(:,3*I-1)
			            A(:,3*I)=A(:,3*I)+HB(:,3*I)
        !
			            A(:,(3*I-2))=A(:,(3*I-2))+HB(:,(3*JJ-2))
			            A(:,3*I-1)=A(:,3*I-1)+HB(:,3*JJ-1)
			            A(:,3*I)=A(:,3*I)+HB(:,3*JJ)
                    
                        A_LINHA(:,(3*I-2))=A_LINHA(:,(3*I-2))+H_FC(:,(3*I-2))
			            A_LINHA(:,3*I-1)=A_LINHA(:,3*I-1)+H_FC(:,3*I-1)
			            A_LINHA(:,3*I)=A_LINHA(:,3*I)+H_FC(:,3*I)

			            A_LINHA(:,(3*I-2))=A_LINHA(:,(3*I-2))+H_FC(:,(3*JJ-2))
			            A_LINHA(:,3*I-1)=A_LINHA(:,3*I-1)+H_FC(:,3*JJ-1)
			            A_LINHA(:,3*I)=A_LINHA(:,3*I)+H_FC(:,3*JJ)
                    ELSE
                        IF(JJ.EQ.-1) THEN
				            A(:,(3*I-2))=A(:,(3*I-2))+GB(:,(3*I-2))
				            A(:,3*I-1)=A(:,3*I-1)+GB(:,3*I-1)
				            A(:,3*I)=A(:,3*I)+GB(:,3*I)
            !
				            A(:,(3*I-2))=A(:,(3*I-2))-GB(:,(3*II-2))
				            A(:,3*I-1)=A(:,3*I-1)-GB(:,3*II-1)
				            A(:,3*I)=A(:,3*I)-GB(:,3*II)
                        
                            A_LINHA(:,(3*I-2))=A_LINHA(:,(3*I-2))+G_FC(:,(3*I-2))
				            A_LINHA(:,3*I-1)=A_LINHA(:,3*I-1)+G_FC(:,3*I-1)
				            A_LINHA(:,3*I)=A_LINHA(:,3*I)+G_FC(:,3*I)

				            A_LINHA(:,(3*I-2))=A_LINHA(:,(3*I-2))-G_FC(:,(3*II-2))
				            A_LINHA(:,3*I-1)=A_LINHA(:,3*I-1)-G_FC(:,3*II-1)
				            A_LINHA(:,3*I)=A_LINHA(:,3*I)-G_FC(:,3*II)
                        ELSE  
                            IF(B_CONDITIONS(3*I-2).EQ.0) THEN
				                A(:,(3*I-2))=-GB(:,(3*I-2))
				                A_LINHA(:,(3*I-2))=-G_FC(:,(3*I-2))

				                CONT=CONT+1
				                B(:,CONT)=-HB(:,(3*I-2))
				                B_LINHA(:,CONT)=-H_FC(:,(3*I-2))
				                F(CONT)=U((3*I-2))
			                ELSE
				                A(:,(3*I-2))=HB(:,(3*I-2))
				                A_LINHA(:,(3*I-2))=H_FC(:,(3*I-2))

				                CONT=CONT+1
				                B(:,CONT)=GB(:,(3*I-2))
				                B_LINHA(:,CONT)=G_FC(:,(3*I-2))
				                F(CONT)=T((3*I-2))
			                ENDIF

			                IF(B_CONDITIONS(3*I-1).EQ.0) THEN
				                A(:,(3*I-1))=-GB(:,(3*I-1))
				                A_LINHA(:,(3*I-1))=-G_FC(:,(3*I-1))

				                CONT=CONT+1
				                B(:,CONT)=-HB(:,(3*I-1))
				                B_LINHA(:,CONT)=-H_FC(:,(3*I-1))
				                F(CONT)=U((3*I-1))
			                ELSE
				                A(:,(3*I-1))=HB(:,(3*I-1))
				                A_LINHA(:,(3*I-1))=H_FC(:,(3*I-1))

				                CONT=CONT+1
				                B(:,CONT)=GB(:,(3*I-1))
				                B_LINHA(:,CONT)=G_FC(:,(3*I-1))
				                F(CONT)=T((3*I-1))
                            ENDIF
                    
                            IF(B_CONDITIONS(3*I).EQ.0) THEN
				                A(:,(3*I))=-GB(:,(3*I))
				                A_LINHA(:,(3*I))=-G_FC(:,(3*I))

				                CONT=CONT+1
				                B(:,CONT)=-HB(:,(3*I))
				                B_LINHA(:,CONT)=-H_FC(:,(3*I))
				                F(CONT)=U((3*I))
			                ELSE
				                A(:,(3*I))=HB(:,(3*I))
				                A_LINHA(:,(3*I))=H_FC(:,(3*I))

				                CONT=CONT+1
				                B(:,CONT)=GB(:,(3*I))
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
            ENDDO  ! End of collocation points loop
        
            !
            ! ------ CONSTRUCTING THE AGLOBAL MATRIX -----------
            !
            !DO J=1,(3*N_COLLOCPOINTS)
            !    DO I=1,(3*N_COLLOCPOINTS)
            !        AGLOBAL(I,J)=A(I,J)
            !    ENDDO
            !ENDDO

            !JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
            !DO J=1,(3*N_NOS_MEF)
            !    DO I=1,(3*N_COLLOCPOINTS)
            !        AGLOBAL(I,(JJ+J))=-G_CF(I,J)
            !    ENDDO
            !ENDDO

            !JJ=3*N_COLLOCPOINTS
            !DO J=1,(3*N_COLLOCPOINTS)
            !    DO I=1,(3*N_NOS_MEF)
            !        AGLOBAL((JJ+I),J)=A_LINHA(I,J)
            !    ENDDO
            !ENDDO

            II=3*N_COLLOCPOINTS
            JJ=3*N_COLLOCPOINTS
            DO I=1,N_NOS_MEF
                CONT=0
                DO J=1,N_ELEMENTOS_MEF
                    DO K=1,ORDEM_MEF(J)+1
                        IF (CONECTI_MEF(J,K).EQ.I) CONT=ND_MEF(J)
                    ENDDO
                ENDDO
                ! VISCOELASTIC COEFFICIENT ----------------
	            NUM_LEI = TIME_MODEL(CONT)
                GAMA = EMP_VISCO(CONT,1)
                SELECT CASE(NUM_LEI)
                CASE(1)!KELVIN'S MODEL
	                  K1= 1.0D0 + GAMA/DELTA_T_VISCO
                CASE(2)!BOLTZMANN'S MODEL
    	              K1= 1.0D0 + GAMA/DELTA_T_VISCO
                CASE(3)!MAXWELL'S MODEL
                      K1= 1.0D0
                CASE(4)!HOOKE'S MODEL
       	              K1= 1.0D0
                ENDSELECT
                A_GLOBAL_ELASTIC((II+3*I-2),(JJ+3*I-2))=1.0D0*K1
                A_GLOBAL_ELASTIC((II+3*I-1),(JJ+3*I-1))=1.0D0*K1
                A_GLOBAL_ELASTIC((II+3*I),(JJ+3*I))=1.0D0*K1
            ENDDO

            !II=3*N_COLLOCPOINTS
            !JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
            !DO J=1,(3*N_NOS_MEF)
            !    DO I=1,(3*N_NOS_MEF)
            !        AGLOBAL((II+I),(JJ+J))=-G_FF(I,J)
            !    ENDDO
            !ENDDO       

            !II=3*N_COLLOCPOINTS+3*N_NOS_MEF
            !JJ=3*N_COLLOCPOINTS
            !DO J=1,(3*N_NOS_MEF)
            !    DO I=1,(3*N_NOS_MEF)
            !        AGLOBAL((II+I),(JJ+J))=KGLOBAL(I,J)
            !    ENDDO
            !ENDDO
            !II=3*N_COLLOCPOINTS+3*N_NOS_MEF
            !JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
            !DO J=1,(3*N_NOS_MEF)
            !    DO I=1,(3*N_NOS_MEF)
            !        AGLOBAL((II+I),(JJ+J))=GGLOBAL(I,J)
            !    ENDDO
            !ENDDO
               
            !
            !       CONSTRUCTING THE BGLOBAL MATRIX -------------------------------     
            ! 
            
            !II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            !DO J=1,(II)
            !    DO I=1,(3*N_COLLOCPOINTS)
            !        B_GLOBAL(I,J)=B(I,J)
            !    ENDDO
            !ENDDO
            !II=3*N_COLLOCPOINTS
            !JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            !DO J=1,(JJ)
            !    DO I=1,(3*N_NOS_MEF)
            !        B_GLOBAL((I+II),J)=B_LINHA(I,J)
            !    ENDDO
            !ENDDO           
        
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            DO J=1,3*N_COLLOCPOINTS
                DO I=1,3*N_COLLOCPOINTS
                    B_GLOBAL(I,J+JJ) = HBB(I,J)
                ENDDO
            ENDDO
            DEALLOCATE(HBB)
            
            II=0
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+3*N_COLLOCPOINTS+3*N_NOS_MEF
            DO J=1,3*N_COLLOCPOINTS
                DO I=1,3*N_COLLOCPOINTS
                    B_GLOBAL(II+I,JJ+J) = GBB(I,J)
                ENDDO
            ENDDO
            DEALLOCATE(GBB)
        
            !II=3*N_COLLOCPOINTS
            !JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            !DO j=1,3*N_COLLOCPOINTS    
            !    DO i=1,3*N_NOS_MEF
            !        B_GLOBAL(I+II,J+JJ) = HH_FC(I,J)
            !    ENDDO
            !ENDDO
        
            II=3*N_COLLOCPOINTS
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+3*N_COLLOCPOINTS
            DO I=1,N_NOS_MEF
                CONT=0
	            DO J=1,N_ELEMENTOS_MEF
	                  DO K=1,(ORDEM_MEF(J)+1)
	                      IF(CONECTI_MEF(J,K).EQ.I) CONT=ND_MEF(J) !CONT guarda qual dominio o nó atual do mef pertence
                      ENDDO
                ENDDO
                ! VISCOELASTIC COEFFICIENTS --------------
	            NUM_LEI = TIME_MODEL(CONT)
                GAMA = EMP_VISCO(CONT,1)
                SELECT CASE(NUM_LEI)
                CASE(1)!KELVIN'S MODEL
	                K3= GAMA/DELTA_T_VISCO
                CASE(2)!BOLTZMANN'S MODEL
	                K3= GAMA/DELTA_T_VISCO
                CASE(3)!MAXWELL'S MODEL
	                K3= 1.0D0
                CASE(4)!HOOKE'S MODEL
	                K3= 0.0D0
                ENDSELECT
	            ! ----------------------------------------
                B_GLOBAL(II+3*I-2,JJ+3*I-2) = 1.0D0*K3
                B_GLOBAL(II+3*I-1,JJ+3*I-1) = 1.0D0*K3
                B_GLOBAL(II+3*I,JJ+3*I) = 1.0D0*K3
            ENDDO
                
            !II=3*N_COLLOCPOINTS
            !JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+3*N_COLLOCPOINTS+3*N_NOS_MEF
            !DO J=1,3*N_COLLOCPOINTS
            !    DO I=1,3*N_NOS_MEF
            !        B_GLOBAL(II+I,JJ+J) = GG_FC(I,J)
            !    ENDDO
            !ENDDO
        
            !II=0
            !JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+6*N_COLLOCPOINTS+3*N_NOS_MEF
            !DO J=1,3*N_NOS_MEF
            !    DO I=1,3*N_COLLOCPOINTS
            !        B_GLOBAL(II+I,JJ+J) = GG_CF(I,J)
            !    ENDDO
            !ENDDO
        
            !II=3*N_COLLOCPOINTS
            !JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+6*N_COLLOCPOINTS+3*N_NOS_MEF
            !DO J=1,3*N_NOS_MEF
            !    DO I=1,3*N_NOS_MEF
            !        B_GLOBAL(II+I,JJ+J) = GG_FF(I,J)
            !    ENDDO
            !ENDDO
            
            !IF (VISCO_ANALYSIS_REINF) THEN
            !    II=3*N_COLLOCPOINTS+3*N_NOS_MEF
            !    JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+3*N_COLLOCPOINTS
            !    DO J=1,3*N_NOS_MEF
            !        DO I=1,3*N_NOS_MEF
            !            B_GLOBAL(II+I,JJ+J) = K_KGLOBAL(I,J)
            !        ENDDO
            !    ENDDO
            !    II=3*N_COLLOCPOINTS+3*N_NOS_MEF
            !    JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+6*N_COLLOCPOINTS+3*N_NOS_MEF
            !    DO J=1,3*N_NOS_MEF
            !        DO I=1,3*N_NOS_MEF
            !            B_GLOBAL(II+I,JJ+J) = G_GGLOBAL(I,J)
            !        ENDDO
            !    ENDDO
            !ENDIF
            
            !
            !     CREATING FGLOBAL WITH FREE TERMS ------------------------------
            !
            JJ = 0
            DO i=1,3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
                FGLOBAL(i+jj) = F(i)
            ENDDO
         
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            DO i=1,3*N_COLLOCPOINTS
                FGLOBAL(i+jj) = U_ANTERIOR(i) 
            ENDDO
      
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+3*N_COLLOCPOINTS
            DO i=1,3*N_NOS_MEF
                FGLOBAL(i+jj) = U_ANTERIOR_MEF(i)
            ENDDO
      
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+3*N_COLLOCPOINTS+3*N_NOS_MEF
            DO i=1,3*N_COLLOCPOINTS
                FGLOBAL(i+jj) = T_ANTERIOR(i)
            ENDDO
      
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+6*N_COLLOCPOINTS+3*N_NOS_MEF
            DO i=1,3*N_NOS_MEF
                FGLOBAL(i+jj) = P_ANTERIOR_MEF(i)
            ENDDO 
      
            !          
            !     MULTIPLYING BGLOBAL AND FGLOBAL (FREE TERMS) ---------------------------
            !
       	    JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
 	        II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+6*N_NOS_MEF+6*N_COLLOCPOINTS
            !CALL DMURRV(JJ,II,BGLOBAL,JJ,II,FGLOBAL,1,JJ,SOLUTION)
            CALL DGEMV('N',JJ,II,1.0D0,B_GLOBAL,JJ,FGLOBAL,1,0.0D0,SOLUTION,1)
                
            ! Resolvendo sistema
            !CALL DLSARG(3*N_COLLOCPOINTS,A,3*N_COLLOCPOINTS,F,1,SOLUTION)
        
            ! INVERTENDO A MATRIZ AGLOBAL E GUARDANDO EM A_GLOBAL_LU (FATORACAO LU):
            JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
            ALLOCATE(A_GLOBAL_LU(JJ,JJ),IPVT_AG(JJ))
            
            A_GLOBAL_LU=A_GLOBAL_ELASTIC
            CALL DGETRF(JJ,JJ,A_GLOBAL_LU,JJ,IPVT_AG,INFO)
            
            IF (INFO.NE.0) THEN
                WRITE(*,*)'*************** CAUTION **********************'
                WRITE(*,*)'******* INFO NOT EQUAL ZERO IN DGETRF ********'
                WRITE(*,*)'* VERIFY SOLVE_BVP_VISCO WITH REINFORCEMENTS *'
                READ(*,*)
            ENDIF
                                   
        ELSE            ! A_GLOBAL_LU E BGLOBAL JÁ CALCULADOS
        
            ! Starting variables
            II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
            ALLOCATE(SOLUTION(JJ),F(II),FGLOBAL(II+6*N_NOS_MEF+6*N_COLLOCPOINTS))
            FGLOBAL=0.0D0
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
            ENDDO  ! End of collocation points loop
            
            !
            !     CREATING FGLOBAL WITH FREE TERMS ------------------------------
            !
            JJ = 0
            DO i=1,3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
                FGLOBAL(i+jj) = F(i)
            ENDDO
         
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
            DO i=1,3*N_COLLOCPOINTS
                FGLOBAL(i+jj) = U_ANTERIOR(i) 
            ENDDO
      
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+3*N_COLLOCPOINTS
            DO i=1,3*N_NOS_MEF
                FGLOBAL(i+jj) = U_ANTERIOR_MEF(i)
            ENDDO
      
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+3*N_COLLOCPOINTS+3*N_NOS_MEF
            DO i=1,3*N_COLLOCPOINTS
                FGLOBAL(i+jj) = T_ANTERIOR(i)
            ENDDO
      
            JJ=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+6*N_COLLOCPOINTS+3*N_NOS_MEF
            DO i=1,3*N_NOS_MEF
                FGLOBAL(i+jj) = P_ANTERIOR_MEF(i)
            ENDDO 
      
            !          
            !     MULTIPLYING BGLOBAL AND FGLOBAL (FREE TERMS) ---------------------------
            !
       	    JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
 	        II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)+6*N_NOS_MEF+6*N_COLLOCPOINTS
            !CALL DMURRV(JJ,II,BGLOBAL,JJ,II,FGLOBAL,1,JJ,SOLUTION)
            CALL DGEMV('N',JJ,II,1.0D0,B_GLOBAL,JJ,FGLOBAL,1,0.0D0,SOLUTION,1)
            
        ENDIF
        
        ! RESOLVENDO O SISTEMA COM LU FATORADO:  
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF+3*N_NOS_MEF
        CALL DGETRS('N',JJ,1,A_GLOBAL_LU,JJ,IPVT_AG,SOLUTION,JJ,INFO)
        
        DEALLOCATE(F,FGLOBAL)
        
        ! Salvando valores de U e T encontrados
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
    !
		        IF(II.EQ.-1) THEN
			        U(3*I-2)=SOLUTION(3*I-2)
			        U(3*I-1)=SOLUTION(3*I-1)
			        U(3*I)=SOLUTION(3*I)
    !
			        U(3*JJ-2)=SOLUTION(3*I-2)
			        U(3*JJ-1)=SOLUTION(3*I-1)
			        U(3*JJ)=SOLUTION(3*I)
		        ELSE
			        IF(JJ.EQ.-1) THEN
				        T(3*I-2)=-SOLUTION(3*I-2)
				        T(3*I-1)=-SOLUTION(3*I-1)
			            T(3*I)=-SOLUTION(3*I)
    !
				        T(3*II-2)=SOLUTION(3*I-2)
				        T(3*II-1)=SOLUTION(3*I-1)
				        T(3*II)=SOLUTION(3*I)
			        ELSE
				        IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					        T(3*I-2)=SOLUTION(3*I-2)
				        ELSE
					        U(3*I-2)=SOLUTION(3*I-2)
				        ENDIF
    !
				        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					        T(3*I-1)=SOLUTION(3*I-1)
				        ELSE
					        U(3*I-1)=SOLUTION(3*I-1)
				        ENDIF
    !
				        IF (B_CONDITIONS(3*I).EQ.0) THEN
					        T(3*I)=SOLUTION(3*I)
				        ELSE
					        U(3*I)=SOLUTION(3*I)
				        ENDIF
			        ENDIF
		        ENDIF
		    ELSE
		        DO J=1,N_CRACK_POINTS
			        IF(I.EQ.NC(J,1)) THEN
				        II=-1
				        JJ=NC(J,2)
			        ELSE
				        IF(I.EQ.NC(J,2)) THEN
					        II=NC(J,1)
					        JJ=-1
				        ENDIF
			        ENDIF
		        ENDDO	
    !
		        IF(II.EQ.-1) THEN
			        U(3*I-2)=SOLUTION(3*I-2)
			        U(3*JJ-2)=-SOLUTION(3*I-2)
    !			    		    
			        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					    T(3*I-1)=SOLUTION(3*I-1)
				    ELSE
					    U(3*I-1)=SOLUTION(3*I-1)
				    ENDIF
    !
				    IF (B_CONDITIONS(3*I).EQ.0) THEN
					    T(3*I)=SOLUTION(3*I)
				    ELSE
					    U(3*I)=SOLUTION(3*I)
				    ENDIF
    !				    		    
		        ELSE
			        IF(JJ.EQ.-1) THEN
			            T(3*I-2)=SOLUTION(3*I-2)
				        T(3*II-2)=SOLUTION(3*I-2)
    !				    			    
			            IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					        T(3*I-1)=SOLUTION(3*I-1)
				        ELSE
					        U(3*I-1)=SOLUTION(3*I-1)
				        ENDIF
    !
				        IF (B_CONDITIONS(3*I).EQ.0) THEN
					        T(3*I)=SOLUTION(3*I)
				        ELSE
					        U(3*I)=SOLUTION(3*I)
				        ENDIF			    
    !			    
			        ELSE
				        IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					        T(3*I-2)=SOLUTION(3*I-2)
				        ELSE
					        U(3*I-2)=SOLUTION(3*I-2)
				        ENDIF
    !
				        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					        T(3*I-1)=SOLUTION(3*I-1)
				        ELSE
					        U(3*I-1)=SOLUTION(3*I-1)
				        ENDIF
    !
				        IF (B_CONDITIONS(3*I).EQ.0) THEN
					        T(3*I)=SOLUTION(3*I)
				        ELSE
					        U(3*I)=SOLUTION(3*I)
				        ENDIF
			        ENDIF
		        ENDIF			
		    ENDIF
        ENDDO  
        
        II=3*N_COLLOCPOINTS
        JJ=3*N_COLLOCPOINTS+3*N_NOS_MEF
        DO I=1,(3*N_NOS_MEF)
            U_MEF(I)=SOLUTION(II+I)
            P_MEF(I)=SOLUTION(JJ+I)
        ENDDO 
        
        DEALLOCATE(SOLUTION)
        
    ENDIF
    
END SUBROUTINE
    