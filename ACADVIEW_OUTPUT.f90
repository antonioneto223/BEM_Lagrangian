	SUBROUTINE ACADVIEW_OUTPUT
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
	INTEGER::I,J,K,II,JJ,KK,L,NEN,C,N_CRACKED_ELEM,CONT,NNUM,NNODES
!	
	REAL*8::QSI1,QSI2,U_AUX(3),T_AUX(3),U_NODES(3*N_COLLOCPOINTS),T_NODES(3*N_COLLOCPOINTS),VALUES_U[ALLOCATABLE](:,:),VALUES_T[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:),PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),X,Y,Z,DX,DY,DZ,R,R_S(3,3),R_H(3,3)
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
!			
	OPEN(3,file='Output_data\ACADVIEW_OUTPUT.ogl',status='unknown')
!
	WRITE(3,*)'Arquivo de pós-processamento'
	WRITE(3,*)
	WRITE(3,*)'Nº de nos     Nº de elementos     Nº de listas'
	WRITE(3,10)
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        IF (EP_ANALYSIS_REINF) THEN
            WRITE(3,100)N_COLLOCPOINTS+N_NOS_MEF,N_ELEM-N_CRACKED_ELEM+N_ELEMENTOS_MEF,8
        ELSE
            WRITE(3,100)N_COLLOCPOINTS+N_NOS_MEF,N_ELEM-N_CRACKED_ELEM+N_ELEMENTOS_MEF,7
        ENDIF
    ELSE
        WRITE(3,100)N_COLLOCPOINTS,N_ELEM-N_CRACKED_ELEM,6
    ENDIF
    WRITE(3,*)
	WRITE(3,*)'coordx coordy coordz deslx delsy deslz'
	WRITE(3,10)
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
	                        WRITE(3,200)COORD_NODES(NODES_CONNECTIVITY(J,K),1),COORD_NODES(NODES_CONNECTIVITY(J,K),2),COORD_NODES(NODES_CONNECTIVITY(J,K),3),U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I)
	                    ENDIF
	                ENDIF
	            ENDDO
	        ENDIF    
	    ENDDO
    ENDDO
    IF(EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        DO I=1,N_NOS_MEF
                WRITE(3,200)COORDNOS_MEF(I,1),COORDNOS_MEF(I,2),COORDNOS_MEF(I,3),U_MEF(3*I-2),U_MEF(3*I-1),U_MEF(3*I)
        ENDDO
    ENDIF

	WRITE(3,*)
	WRITE(3,*)'tpelem (1 - barra / 2 - triang / 3 - quad) grauaprox nó1 nó2...non'
	WRITE(3,10)
	DO I=1,N_ELEM
	    IF(CRACKED_ELEM(I).EQ."UNCRACKED")THEN	
	        NEN=ELEM_TYPE(I)*ORDER_ELEM(I)+(ELEM_TYPE(I)-3)*(ORDER_ELEM(I)-1)*POL_FAMILY(I)
            IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
                SELECT CASE(NEN)
                CASE(3)
                    WRITE(3,305)2,1,COLLOCPOINTS_CONNECTIVITY(I,3),COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,2),0
                CASE(6)
                    WRITE(3,306)2,2,COLLOCPOINTS_CONNECTIVITY(I,3),COLLOCPOINTS_CONNECTIVITY(I,6),COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,5),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,2),0
                CASE(4)
                    WRITE(3,307)3,1,COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,2),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,3),0
                CASE(8)
                    WRITE(3,308)3,2,COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,5),COLLOCPOINTS_CONNECTIVITY(I,2),COLLOCPOINTS_CONNECTIVITY(I,8),COLLOCPOINTS_CONNECTIVITY(I,6),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,7),COLLOCPOINTS_CONNECTIVITY(I,3),0
                CASE(9)
                    WRITE(3,309)3,2,COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,5),COLLOCPOINTS_CONNECTIVITY(I,2),COLLOCPOINTS_CONNECTIVITY(I,8),COLLOCPOINTS_CONNECTIVITY(I,9),COLLOCPOINTS_CONNECTIVITY(I,6),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,7),COLLOCPOINTS_CONNECTIVITY(I,3),0
                ENDSELECT  
            ELSE
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
        ENDIF
    ENDDO
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        DO I=1,N_ELEMENTOS_MEF
            WRITE(3,'(I9,I8)',ADVANCE='NO')1,ORDEM_MEF(I)
            DO J=1,ORDEM_MEF(I)+1
                WRITE(3,'(I12)',ADVANCE='NO')CONECTI_MEF(I,J)+N_COLLOCPOINTS
            ENDDO
            WRITE(3,'(I2)')1
        ENDDO
    ENDIF
    
	WRITE(3,*)
	WRITE(3,*)'Deslocamento Ux no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'Ux'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),U_NODES(3*I-2)    
    ENDDO
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        DO I=1,N_NOS_MEF
            WRITE(3,400)U_MEF(3*I-2),U_MEF(3*I-1),U_MEF(3*I),U_MEF(3*I-2) 
        ENDDO
    ENDIF
    
	WRITE(3,*)
	WRITE(3,*)'Deslocamento Uy no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'Uy'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),U_NODES(3*I-1)    
    ENDDO
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        DO I=1,N_NOS_MEF
            WRITE(3,400)U_MEF(3*I-2),U_MEF(3*I-1),U_MEF(3*I),U_MEF(3*I-1) 
        ENDDO
    ENDIF
    
	WRITE(3,*)
	WRITE(3,*)'Deslocamento Uz no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'Uz'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),U_NODES(3*I)    
    ENDDO
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        DO I=1,N_NOS_MEF
            WRITE(3,400)U_MEF(3*I-2),U_MEF(3*I-1),U_MEF(3*I),U_MEF(3*I) 
        ENDDO
    ENDIF
        
	WRITE(3,*)
	WRITE(3,*)'Força Tx no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'Tx'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),T_NODES(3*I-2)    
    ENDDO
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        DO I=1,N_NOS_MEF
            WRITE(3,400)U_MEF(3*I-2),U_MEF(3*I-1),U_MEF(3*I),P_MEF(I) 
        ENDDO
    ENDIF
    
	WRITE(3,*)
	WRITE(3,*)'Força Ty no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'Ty'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),T_NODES(3*I-1)    
    ENDDO
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        DO I=1,N_NOS_MEF
            WRITE(3,400)U_MEF(3*I-2),U_MEF(3*I-1),U_MEF(3*I),0.0
        ENDDO
    ENDIF
        
	WRITE(3,*)
	WRITE(3,*)'Força Tz no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'Tz'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),T_NODES(3*I)    
    ENDDO
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        DO I=1,N_NOS_MEF
            WRITE(3,400)U_MEF(3*I-2),U_MEF(3*I-1),U_MEF(3*I),0.0
        ENDDO
    ENDIF
    
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
    	WRITE(3,*)
	    WRITE(3,*)'Esforço normal em nós de fibra'
	    WRITE(3,10)
	    WRITE(3,*)'N'
	    DO I=1,N_COLLOCPOINTS
	        WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),0.0d0   
        ENDDO
        DO I=1,N_NOS_MEF
            WRITE(3,400)U_MEF(3*I-2),U_MEF(3*I-1),U_MEF(3*I),E_NORMAL_NOD(I) 
        ENDDO
        IF (EP_ANALYSIS_REINF) THEN
            	WRITE(3,*)
	            WRITE(3,*)'Deformacao plastica em nós de fibra'
	            WRITE(3,10)
	            WRITE(3,*)'Dp'
	            DO I=1,N_COLLOCPOINTS
	                WRITE(3,400)U_NODES(3*I-2),U_NODES(3*I-1),U_NODES(3*I),0.0d0   
                ENDDO
                DO I=1,N_NOS_MEF
                    WRITE(3,400)U_MEF(3*I-2),U_MEF(3*I-1),U_MEF(3*I),PL_STRAIN_NOD(I) 
                ENDDO
        ENDIF
    ENDIF   
!	
	CLOSE(3)
!
    !write(filename2,'(a,a1,i0)') trim(filename),'.',step
    
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        OPEN(80,file='Output_data\Paraview\BOUNDARY_OUTPUT_INITIAL_MESH.vtk',status='unknown')
    ELSE
        OPEN(80,file='Output_data\Paraview\BOUNDARY_OUTPUT_INITIAL_MESH.vtk',status='unknown')
    ENDIF
!
    WRITE(80,'(a)')'# vtk DataFile Version 2.0'
    WRITE(80,'(a)')'Animation data'
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
	                        WRITE(80,'(3(f0.10,1x))')COORD_NODES(NODES_CONNECTIVITY(J,K),1),COORD_NODES(NODES_CONNECTIVITY(J,K),2),COORD_NODES(NODES_CONNECTIVITY(J,K),3)
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
    
    ! ------------------------------------------------------------------------------------------------------
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        OPEN(80,file='Output_data\Paraview\BOUNDARY_OUTPUT.1.vtk',status='unknown')
    ELSE
        OPEN(80,file='Output_data\Paraview\BOUNDARY_OUTPUT.1.vtk',status='unknown')
    ENDIF
!
    WRITE(80,'(a)')'# vtk DataFile Version 2.0'
    WRITE(80,'(a)')'Animation data'
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
	                        WRITE(80,'(3(f0.10,1x))')COORD_NODES(NODES_CONNECTIVITY(J,K),1)+U_NODES(3*NODES_CONNECTIVITY(J,K)-2),&
                                COORD_NODES(NODES_CONNECTIVITY(J,K),2)+U_NODES(3*NODES_CONNECTIVITY(J,K)-1),&
                                COORD_NODES(NODES_CONNECTIVITY(J,K),3)+U_NODES(3*NODES_CONNECTIVITY(J,K))
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
!------------------------------------------------------------------------------------------------------------------------------------------------------
!         U AND T LOCAL COORDINATE SYSTEM
!------------------------------------------------------------------------------------------------------------------------------------------------------
!	
	END SUBROUTINE