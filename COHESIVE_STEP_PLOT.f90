	SUBROUTINE COHESIVE_STEP_PLOT(STEP)
!
	USE ISOPARAMETRIC_MESH
	USE CRACKS_DUAL_BEM	
	USE ANALYSIS
	USE REMESHING
    USE COHESIVE
!
	IMPLICIT NONE 
!
	INTEGER::STEP,I,J,K,II,JJ,KK,L,NEN,C,N_CRACKED_ELEM
!
	REAL*8::U_AUX(3),T_AUX(3),R_S(3,3),R_H(3,3),DX,DY,DZ,R
!
    CHARACTER(40)::NOME			
!    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!         U AND T GLOBAL COORDINATE SYSTEM
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
	        CALL DLINRG(3,R_S,3,R_S,3)
!	        
	        DO K=1,3
	            U_AUX(K)=R_S(K,1)*CHS_U(3*I-2)+R_S(K,2)*CHS_U(3*I-1)+R_S(K,3)*CHS_U(3*I)
	        	T_AUX(K)=R_S(K,1)*CHS_T(3*I-2)+R_S(K,2)*CHS_T(3*I-1)+R_S(K,3)*CHS_T(3*I) 
	        ENDDO
	        CHS_U(3*I-2)=U_AUX(1) 
	        CHS_U(3*I-1)=U_AUX(2) 
	        CHS_U(3*I)=U_AUX(3) 
	        CHS_T(3*I-2)=T_AUX(1) 
	        CHS_T(3*I-1)=T_AUX(2) 
	        CHS_T(3*I)=T_AUX(3) 
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
	            U_AUX(K)=R_H(K,1)*CHS_U(3*II-2)+R_H(K,2)*CHS_U(3*II-1)+R_H(K,3)*CHS_U(3*II)
	        	T_AUX(K)=R_H(K,1)*CHS_T(3*II-2)+R_H(K,2)*CHS_T(3*II-1)+R_H(K,3)*CHS_T(3*II) 
	        ENDDO
	        CHS_U(3*II-2)=U_AUX(1) 
	        CHS_U(3*II-1)=U_AUX(2) 
	        CHS_U(3*II)=U_AUX(3) 
	        CHS_T(3*II-2)=T_AUX(1) 
	        CHS_T(3*II-1)=T_AUX(2) 
	        CHS_T(3*II)=T_AUX(3) 		    	        	        	           	   	    
	    ENDIF
	ENDDO	
!    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!         U AND T GLOBAL COORDINATE SYSTEM
!------------------------------------------------------------------------------------------------------------------------------------------------------		 	
!	   
	N_CRACKED_ELEM=0
	DO I=1,N_ELEM
	    IF(CRACKED_ELEM(I).NE."UNCRACKED")THEN
	        N_CRACKED_ELEM=N_CRACKED_ELEM+1
	    ENDIF
    ENDDO
!	
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   WRITTING THE ACADVIEW FILE
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
	OPEN(90,file='LIXO',status='unknown')
		WRITE(90,500)'COHESIVE_OUTPUT_IT',STEP,'.ogl'
	CLOSE(90)
!	
	OPEN(90,file='LIXO',status='OLD')
		READ(90,500)NOME
	CLOSE(90, STATUS = 'DELETE')
	
	NOME=ADJUSTL(NOME)	
		
	OPEN(3,file='Output_data\'//NOME,status='unknown')
!
	WRITE(3,*)'Arquivo de pós-processamento'
	WRITE(3,*)
	WRITE(3,*)'Nº de nos     Nº de elementos     Nº de listas'
	WRITE(3,10)
	WRITE(3,100)N_COLLOCPOINTS,N_ELEM-N_CRACKED_ELEM,6
    WRITE(3,*)
	WRITE(3,*)'coordx coordy coordz deslx delsy deslz'
	WRITE(3,10)
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,200)COORD_COLLOCPOINTS(I,1),COORD_COLLOCPOINTS(I,2),COORD_COLLOCPOINTS(I,3),CHS_U(3*I-2),CHS_U(3*I-1),CHS_U(3*I)
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'tpelem (1 - barra / 2 - triang / 3 - quad) grauaprox nó1 nó2...non'
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
	WRITE(3,*)'Deslocamento Ux no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'Ux'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)CHS_U(3*I-2),CHS_U(3*I-1),CHS_U(3*I),CHS_U(3*I-2)    
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'Deslocamento Uy no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'Uy'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)CHS_U(3*I-2),CHS_U(3*I-1),CHS_U(3*I),CHS_U(3*I-1)    
	ENDDO	
	WRITE(3,*)
	WRITE(3,*)'Deslocamento Uz no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'Uz'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)CHS_U(3*I-2),CHS_U(3*I-1),CHS_U(3*I),CHS_U(3*I)    
	ENDDO	
	WRITE(3,*)
	WRITE(3,*)'Força Tx no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'Tx'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)CHS_U(3*I-2),CHS_U(3*I-1),CHS_U(3*I),CHS_T(3*I-2)    
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'Força Ty no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'Ty'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)CHS_U(3*I-2),CHS_U(3*I-1),CHS_U(3*I),CHS_T(3*I-1)    
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'Força Tz no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'Tz'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)CHS_U(3*I-2),CHS_U(3*I-1),CHS_U(3*I),CHS_T(3*I)    
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
500 FORMAT(a,i4,a)
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
	            U_AUX(K)=R_S(K,1)*CHS_U(3*I-2)+R_S(K,2)*CHS_U(3*I-1)+R_S(K,3)*CHS_U(3*I)
	        	T_AUX(K)=R_S(K,1)*CHS_T(3*I-2)+R_S(K,2)*CHS_T(3*I-1)+R_S(K,3)*CHS_T(3*I) 
	        ENDDO
	        CHS_U(3*I-2)=U_AUX(1) 
	        CHS_U(3*I-1)=U_AUX(2) 
	        CHS_U(3*I)=U_AUX(3) 
	        CHS_T(3*I-2)=T_AUX(1) 
	        CHS_T(3*I-1)=T_AUX(2) 
	        CHS_T(3*I)=T_AUX(3) 
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
	            U_AUX(K)=R_H(K,1)*CHS_U(3*II-2)+R_H(K,2)*CHS_U(3*II-1)+R_H(K,3)*CHS_U(3*II)
	        	T_AUX(K)=R_H(K,1)*CHS_T(3*II-2)+R_H(K,2)*CHS_T(3*II-1)+R_H(K,3)*CHS_T(3*II) 
	        ENDDO
	        CHS_U(3*II-2)=U_AUX(1) 
	        CHS_U(3*II-1)=U_AUX(2) 
	        CHS_U(3*II)=U_AUX(3) 
	        CHS_T(3*II-2)=T_AUX(1) 
	        CHS_T(3*II-1)=T_AUX(2) 
	        CHS_T(3*II)=T_AUX(3) 		    	        	        	           	   	    
	    ENDIF
	ENDDO		
!    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!         U AND T LOCAL COORDINATE SYSTEM
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
	END SUBROUTINE