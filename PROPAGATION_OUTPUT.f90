	SUBROUTINE PROPAGATION_OUTPUT(STEP)
!
	USE ISOPARAMETRIC_MESH
	USE CRACKS_DUAL_BEM
    USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
	USE REMESHING
	USE FATIGUE
!
	IMPLICIT NONE 
!
	INTEGER::I,J,K,II,JJ,KK,L,NEN,C,N_CRACKED_ELEM,STEP
!	
	REAL*8::QSI1,QSI2,U_AUX(3),T_AUX(3),U_NODES(3*N_COLLOCPOINTS),T_NODES(3*N_COLLOCPOINTS),VALUES_U[ALLOCATABLE](:,:),VALUES_T[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:),PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),X,Y,Z,DX,DY,DZ,R,R_S(3,3),R_H(3,3)
!
    CHARACTER(40)::NOME	
!
!------------------------------------------------------------------------------------------------------------------------------------------------------
!         U AND T GLOBAL COORDINATE SYSTEM
!------------------------------------------------------------------------------------------------------------------------------------------------------
    DO I=1,N_COLLOCPOINTS
	    IF(DUAL_BEM(I).EQ."S")THEN
	        DO J=1,N_ELEM
		        NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
			    DO K=1,NEN
				    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
					    JJ=J
					    KK=K
				    ENDIF
			    ENDDO
			ENDDO
!			    
	        CALL LOCAL_COODINATE_SYSTEM(KK,JJ,R_S,R_H)
!	        
	        CALL DLINRG(3,R_S,3,R_S,3)
!	        
	        DO K=1,3
	            U_AUX(K)=R_S(K,1)*U(3*I-2)+R_S(K,2)*U(3*I-1)+R_S(K,3)*U(3*I)
	        	T_AUX(K)=R_S(K,1)*T(3*I-2)+R_S(K,2)*T(3*I-1)+R_S(K,3)*T(3*I) 
	        ENDDO
	        U(3*I-2)=U_AUX(1) 
	        U(3*I-1)=U_AUX(2) 
	        U(3*I)=U_AUX(3) 
	        T(3*I-2)=T_AUX(1) 
	        T(3*I-1)=T_AUX(2) 
	        T(3*I)=T_AUX(3) 
!
	        DO J=1,N_COLLOCPOINTS
		        IF(DUAL_BEM(J).EQ."H")THEN
				    DX=COORD_COLLOCPOINTS(I,1)-COORD_COLLOCPOINTS(J,1)
				    DY=COORD_COLLOCPOINTS(I,2)-COORD_COLLOCPOINTS(J,2)
				    DZ=COORD_COLLOCPOINTS(I,3)-COORD_COLLOCPOINTS(J,3)
				    R=DSQRT(DX*DX+DY*DY+DZ*DZ)
				    IF(R.LE.1.E-12)THEN
				        II=J
				    ENDIF
				ENDIF
		    ENDDO
!	        
	        CALL DLINRG(3,R_H,3,R_H,3)
!	        
	        DO K=1,3
	            U_AUX(K)=R_H(K,1)*U(3*II-2)+R_H(K,2)*U(3*II-1)+R_H(K,3)*U(3*II)
	        	T_AUX(K)=R_H(K,1)*T(3*II-2)+R_H(K,2)*T(3*II-1)+R_H(K,3)*T(3*II) 
	        ENDDO
	        U(3*II-2)=U_AUX(1) 
	        U(3*II-1)=U_AUX(2) 
	        U(3*II)=U_AUX(3) 
	        T(3*II-2)=T_AUX(1) 
	        T(3*II-1)=T_AUX(2) 
	        T(3*II)=T_AUX(3) 		    	        	        	           	   	    
	    ENDIF
	ENDDO	   	
!    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   EXTRAPOLATING DISPLACEMENTS AND TRACTIONS
!------------------------------------------------------------------------------------------------------------------------------------------------------
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
!   WRITING PROPAGATION OUTPUT
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
	N_CRACKED_ELEM=0
	DO I=1,N_ELEM
	    IF(CRACKED_ELEM(I).NE."UNCRACKED")THEN
	        N_CRACKED_ELEM=N_CRACKED_ELEM+1
	    ENDIF
	ENDDO
!
	OPEN(90,file='LIXO',status='unknown')
		WRITE(90,500)'PROPAGATION_OUTPUT',STEP,' ',TRIM(MINIMUM_MAXIMUM_LOAD),'.ogl'
	CLOSE(90)
!	
	OPEN(90,file='LIXO',status='OLD')
		READ(90,500)NOME
	CLOSE(90, STATUS = 'DELETE')
	
	NOME=ADJUSTL(NOME)	
		
	OPEN(3,file='Output_data\'//NOME,status='unknown')
!
	WRITE(3,*)'Arquivo de p�s-processamento'
	WRITE(3,*)
	WRITE(3,*)'N� de nos     N� de elementos     N� de listas'
	WRITE(3,10)
	WRITE(3,100)N_COLLOCPOINTS,N_ELEM-N_CRACKED_ELEM,6
    WRITE(3,*)
	WRITE(3,*)'coordx coordy coordz deslx delsy deslz'
	WRITE(3,10)
	DO I=1,N_COLLOCPOINTS
	    C=0
	    DO J=1,N_ELEM
	        IF(CRACKED_ELEM(J).EQ."UNCRACKED")THEN
	            NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
	            DO K=1,NEN
	                IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
	                    C=C+1
	                    IF(C.LE.1)THEN
	                        WRITE(3,200)COORD_NODES(NODES_CONNECTIVITY(J,K),1),COORD_NODES(NODES_CONNECTIVITY(J,K),2),COORD_NODES(NODES_CONNECTIVITY(J,K),3),U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I)
	                    ENDIF
	                ENDIF
	            ENDDO
	        ENDIF    
	    ENDDO
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'tpelem (1 - barra / 2 - triang / 3 - quad) grauaprox n�1 n�2...non'
	WRITE(3,10)
	DO I=1,N_ELEM
	    IF(CRACKED_ELEM(I).EQ."UNCRACKED")THEN	
	        NEN=ELEM_TYPE(I)*ORDER_ELEM(I)+(ELEM_TYPE(I)-3)*(ORDER_ELEM(I)-1)*POL_FAMILY(I)
            SELECT CASE(NEN)
            CASE(3)
                WRITE(3,300)2,1,COLLOCPOINTS_CONNECTIVITY(I,3),COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,2)
            CASE(6)
                WRITE(3,301)2,2,COLLOCPOINTS_CONNECTIVITY(I,3),COLLOCPOINTS_CONNECTIVITY(I,6),COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,5),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,2)
            CASE(4)
                WRITE(3,302)3,1,COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,2),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,3)
            CASE(8)
                WRITE(3,303)3,2,COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,5),COLLOCPOINTS_CONNECTIVITY(I,2),COLLOCPOINTS_CONNECTIVITY(I,8),COLLOCPOINTS_CONNECTIVITY(I,6),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,7),COLLOCPOINTS_CONNECTIVITY(I,3)
            CASE(9)
                WRITE(3,304)3,2,COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,5),COLLOCPOINTS_CONNECTIVITY(I,2),COLLOCPOINTS_CONNECTIVITY(I,8),COLLOCPOINTS_CONNECTIVITY(I,9),COLLOCPOINTS_CONNECTIVITY(I,6),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,7),COLLOCPOINTS_CONNECTIVITY(I,3)              
            ENDSELECT
        ENDIF
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'Deslocamento Ux no sistema de refer�ncia global'
	WRITE(3,10)
	WRITE(3,*)'Ux'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),U_NODES(3*I-2)    
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'Deslocamento Uy no sistema de refer�ncia global'
	WRITE(3,10)
	WRITE(3,*)'Uy'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),U_NODES(3*I-1)    
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'Deslocamento Uz no sistema de refer�ncia global'
	WRITE(3,10)
	WRITE(3,*)'Uz'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),U_NODES(3*I)    
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'For�a Tx no sistema de refer�ncia global'
	WRITE(3,10)
	WRITE(3,*)'Tx'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),T_NODES(3*I-2)    
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'For�a Ty no sistema de refer�ncia global'
	WRITE(3,10)
	WRITE(3,*)'Ty'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),T_NODES(3*I-1)    
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'For�a Tz no sistema de refer�ncia global'
	WRITE(3,10)
	WRITE(3,*)'Tz'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),T_NODES(3*I)    
	ENDDO
!	
	CLOSE(3)
!
10  FORMAT(1h#)
100	FORMAT(2x,i5,12X,i5,12X,i5)
200	FORMAT(8x,F14.9,16x,F14.9,16x,F14.9,16x,F14.9,16X,F14.9,16x,F14.9)
300	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7X,i5)
301	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5)
302	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5)
303	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5)
304	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5)
400	FORMAT(8x,E13.6,4x,E13.6,4x,E13.6,15x,E13.6)
500 FORMAT(a,i4,a,a,a)
!    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!         U AND T LOCAL COORDINATE SYSTEM
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
    DO I=1,N_COLLOCPOINTS
	    IF(DUAL_BEM(I).EQ."S")THEN
	        DO J=1,N_ELEM
		        NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
			    DO K=1,NEN
				    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
					    JJ=J
					    KK=K
				    ENDIF
			    ENDDO
			ENDDO
!			    
	        CALL LOCAL_COODINATE_SYSTEM(KK,JJ,R_S,R_H)
!	        
	        DO K=1,3
	            U_AUX(K)=R_S(K,1)*U(3*I-2)+R_S(K,2)*U(3*I-1)+R_S(K,3)*U(3*I)
	        	T_AUX(K)=R_S(K,1)*T(3*I-2)+R_S(K,2)*T(3*I-1)+R_S(K,3)*T(3*I) 
	        ENDDO
	        U(3*I-2)=U_AUX(1) 
	        U(3*I-1)=U_AUX(2) 
	        U(3*I)=U_AUX(3) 
	        T(3*I-2)=T_AUX(1) 
	        T(3*I-1)=T_AUX(2) 
	        T(3*I)=T_AUX(3) 
!
	        DO J=1,N_COLLOCPOINTS
		        IF(DUAL_BEM(J).EQ."H")THEN
				    DX=COORD_COLLOCPOINTS(I,1)-COORD_COLLOCPOINTS(J,1)
				    DY=COORD_COLLOCPOINTS(I,2)-COORD_COLLOCPOINTS(J,2)
				    DZ=COORD_COLLOCPOINTS(I,3)-COORD_COLLOCPOINTS(J,3)
				    R=DSQRT(DX*DX+DY*DY+DZ*DZ)
				    IF(R.LE.1.E-12)THEN
				        II=J
				    ENDIF
				ENDIF
		    ENDDO
!	        
	        DO K=1,3
	            U_AUX(K)=R_H(K,1)*U(3*II-2)+R_H(K,2)*U(3*II-1)+R_H(K,3)*U(3*II)
	        	T_AUX(K)=R_H(K,1)*T(3*II-2)+R_H(K,2)*T(3*II-1)+R_H(K,3)*T(3*II) 
	        ENDDO
	        U(3*II-2)=U_AUX(1) 
	        U(3*II-1)=U_AUX(2) 
	        U(3*II)=U_AUX(3) 
	        T(3*II-2)=T_AUX(1) 
	        T(3*II-1)=T_AUX(2) 
	        T(3*II)=T_AUX(3) 		    	        	        	           	   	    
	    ENDIF
	ENDDO		
!    	
	END SUBROUTINE