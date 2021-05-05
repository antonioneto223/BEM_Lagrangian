	SUBROUTINE MESH_GENERATION
!
	USE NURBS_SURFACES
	USE TRIMMING_CURVES
	USE ISOPARAMETRIC_MESH
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
	USE DFPORT
!
	IMPLICIT NONE 

    INTEGER::I,J,K,L,M,II,N_MAX,CONT1,CONT2,NEN,ANSYS_VER
    INTEGER::N_CURVEPOINTS,N_PARAMETRICNODES,N_PARAMETRICCOLLOCPOINTS,N_PARAMETRICELEMENTS,PARAMETRIC_CONNECTIVITY[ALLOCATABLE](:,:),COMMON_ELEMENTS[ALLOCATABLE](:),REPEATED_CP[ALLOCATABLE](:)
    INTEGER::ELEM_TYPE_AUX[ALLOCATABLE](:),ORDER_ELEM_AUX[ALLOCATABLE](:),POL_FAMILY_AUX[ALLOCATABLE](:),NODES_CONNECTIVITY_AUX[ALLOCATABLE](:,:),COLLOCPOINTS_CONNECTIVITY_AUX[ALLOCATABLE](:,:),ND_AUX[ALLOCATABLE](:),CRACK_TIP_AUX[ALLOCATABLE](:),EDGE_DISCONTINUITY_AUX[ALLOCATABLE](:,:)
!
    REAL*8::DELTAQSI(2),VALUES[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:),PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),X1,X2,X3
    REAL*8::PARAMETRIC_VERTICES(4,2),TRIMMING_CURVEPOINTS[ALLOCATABLE](:,:),PARAMETRICCOORD_NODES[ALLOCATABLE](:,:),A,DELTAA,QSI1,QSI2,SIGNED_DIST_PARAMETRICNODES[ALLOCATABLE](:)
    REAL*8::COORD_NODES_AUX[ALLOCATABLE](:,:)
    REAL*8::QSI_AUX[ALLOCATABLE](:,:,:)
!
    CHARACTER(:),ALLOCATABLE ::ANSYS_IN,ANSYS_OUT
    CHARACTER(6)::KEYWORD
    CHARACTER(3)::CMD1
    CHARACTER(30)::EQ_TYPE_AUX[ALLOCATABLE](:),DUAL_BEM_AUX[ALLOCATABLE](:),DISCONTINUOUS_COLLOCPOINTS_AUX[ALLOCATABLE](:)
!
    LOGICAL::KEEP_READING         
!	
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   DEFINING THE NODES PARAMETRIC COORDINATES   
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
    N_MAX=0    
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN 
            IF(SUBDIVISION_PD1(I)+1.GE.N_MAX)THEN
                N_MAX=SUBDIVISION_PD1(I)+1
            ENDIF
            IF(SUBDIVISION_PD2(I)+1.GE.N_MAX)THEN
                N_MAX=SUBDIVISION_PD2(I)+1
            ENDIF            
        ENDIF     
    ENDDO  
    ALLOCATE(NODES(N_SURFACES,2,N_MAX))
!    
	DO I=1,N_SURFACES 
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN 		     
            NODES(I,1,:)=KNOT_VECTORS(I,1,1)
            NODES(I,2,:)=KNOT_VECTORS(I,2,1)
        ENDIF    
    ENDDO    
!    
	DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN 	 
	        DELTAQSI(1)=(KNOT_VECTORS(I,1,(N_BF(I,1)+POLYNOMIAL_ORDERS(I,1)+1))-KNOT_VECTORS(I,1,1))/SUBDIVISION_PD1(I)
	        DELTAQSI(2)=(KNOT_VECTORS(I,2,(N_BF(I,2)+POLYNOMIAL_ORDERS(I,2)+1))-KNOT_VECTORS(I,2,1))/SUBDIVISION_PD2(I)
            DO J=1,SUBDIVISION_PD1(I)
                NODES(I,1,J)=NODES(I,1,J)+(J-1)*DELTAQSI(1)
            ENDDO
            DO J=1,SUBDIVISION_PD2(I)
                NODES(I,2,J)=NODES(I,2,J)+(J-1)*DELTAQSI(2) 
            ENDDO	         
	    ENDIF    	    
	ENDDO
!	
	DO I=1,N_SURFACES  
	    IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN     
            NODES(I,1,SUBDIVISION_PD1(I)+1)=KNOT_VECTORS(I,1,(N_BF(I,1)+POLYNOMIAL_ORDERS(I,1)+1))
            NODES(I,2,SUBDIVISION_PD2(I)+1)=KNOT_VECTORS(I,2,(N_BF(I,2)+POLYNOMIAL_ORDERS(I,2)+1))
        ENDIF    
    ENDDO 
!     	   	
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   COMPUTING THE NUMBER OF NODES AND THEIR PHYSICAL COORDINATES
!------------------------------------------------------------------------------------------------------------------------------------------------------
    N_NODES=0
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN 
            N_NODES=N_NODES+(SUBDIVISION_PD1(I)+1)*(SUBDIVISION_PD2(I)+1)
        ENDIF    
    ENDDO	
!    
    ALLOCATE(COORD_NODES(N_NODES,3))
    COORD_NODES=0.D0
    CONT1=0
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN    
            ALLOCATE(VALUES(N_GBF(I),3),PHI(N_GBF(I)))            
            DO J=1,N_GBF(I)                
                VALUES(J,1)=COORD_CONTROLPOINTS(CONTROLPOINTS_CONNECTIVITY(I,J),1)
                VALUES(J,2)=COORD_CONTROLPOINTS(CONTROLPOINTS_CONNECTIVITY(I,J),2)
                VALUES(J,3)=COORD_CONTROLPOINTS(CONTROLPOINTS_CONNECTIVITY(I,J),3)                                            
            ENDDO                 
            DO J=1,SUBDIVISION_PD2(I)+1
                DO K=1,SUBDIVISION_PD1(I)+1
                    CONT1=CONT1+1 
                    CALL PHYSICAL_MAPPING(I,NODES(I,1,K),NODES(I,2,J),VALUES,PHI,X1,X2,X3) 
                    COORD_NODES(CONT1,1)=X1
                    COORD_NODES(CONT1,2)=X2
                    COORD_NODES(CONT1,3)=X3                                                                                     
                ENDDO
            ENDDO
            DEALLOCATE(VALUES,PHI)
        ENDIF                                               
    ENDDO 		
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   DEFINING THE NUMBER OF LINEAR QUADRILATERAL ELEMENTS AND THE NODE CONNECTIVITY
!------------------------------------------------------------------------------------------------------------------------------------------------------	
    N_ELEM=0
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN     
            N_ELEM=N_ELEM+SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)
            N_ELEM_SURFACES(I)=SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)
        ENDIF             
    ENDDO
!    
    ALLOCATE(ELEM_TYPE(N_ELEM),ORDER_ELEM(N_ELEM),POL_FAMILY(N_ELEM),NODES_CONNECTIVITY(N_ELEM,4),COLLOCPOINTS_CONNECTIVITY(N_ELEM,4),QSI(N_ELEM,9,2),ND(N_ELEM),EDGE_DISCONTINUITY(N_ELEM,4))
    ELEM_TYPE=4
    ORDER_ELEM=1
    POL_FAMILY=1
    EDGE_DISCONTINUITY=0
    CONT1=0
    CONT2=0
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN      
            DO K=1,SUBDIVISION_PD2(I)
                DO L=1,SUBDIVISION_PD1(I)
                    CONT2=CONT2+1
                    NODES_CONNECTIVITY(CONT2,1)=CONT1+K*(SUBDIVISION_PD1(I)+1)+L-(SUBDIVISION_PD1(I)+1)
                    NODES_CONNECTIVITY(CONT2,2)=CONT1+K*(SUBDIVISION_PD1(I)+1)+L-(SUBDIVISION_PD1(I)+1)+1
                    NODES_CONNECTIVITY(CONT2,3)=CONT1+K*(SUBDIVISION_PD1(I)+1)+L+1   
                    NODES_CONNECTIVITY(CONT2,4)=CONT1+K*(SUBDIVISION_PD1(I)+1)+L
                ENDDO
            ENDDO
            CONT1=CONT1+(SUBDIVISION_PD1(I)+1)*(SUBDIVISION_PD2(I)+1)
        ENDIF    
    ENDDO
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   COLLOCATION POINTS GENERATION: COMPUTING THE NUMBER OF COLLOCATION POINTS
!------------------------------------------------------------------------------------------------------------------------------------------------------
    CONT1=0
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN     
            SELECT CASE (SURFACE_TYPE(I))
            CASE("B")
                IF (DESCONT_TYPE(I).EQ."DB") THEN
                    CONT1=CONT1+4*SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)
                    N_COLLOCPOINTS_SURFACES(I)=4*SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)    
                ELSE
                    CONT1=CONT1+(SUBDIVISION_PD1(I)+1)*(SUBDIVISION_PD2(I)+1) 
                    N_COLLOCPOINTS_SURFACES(I)=(SUBDIVISION_PD1(I)+1)*(SUBDIVISION_PD2(I)+1)                      
                ENDIF        
            CASE("CD")
                CONT1=CONT1+4*SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)
                N_COLLOCPOINTS_SURFACES(I)=4*SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)
            CASE("CT")
                CONT1=CONT1+4*SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)
                N_COLLOCPOINTS_SURFACES(I)=4*SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)   
            CASE("CI")
                CONT1=CONT1+4*SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)
                N_COLLOCPOINTS_SURFACES(I)=4*SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)                   
            ENDSELECT
        ENDIF    
    ENDDO
    N_COLLOCPOINTS=CONT1
!    
    ALLOCATE(COORD_COLLOCPOINTS(N_COLLOCPOINTS,3),B_CONDITIONS(3*N_COLLOCPOINTS),U(3*N_COLLOCPOINTS),T(3*N_COLLOCPOINTS),EQ_TYPE(N_COLLOCPOINTS),DUAL_BEM(N_COLLOCPOINTS),H(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS),G(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS),DISCONTINUOUS_COLLOCPOINTS(N_COLLOCPOINTS),CRACK_TIP(N_COLLOCPOINTS))  
    COORD_COLLOCPOINTS=0.D0
    B_CONDITIONS=1
    DISCONTINUOUS_COLLOCPOINTS="CONTINUOUS"
    DUAL_BEM="B"
    EQ_TYPE="S"
    CRACK_TIP=0
    U=0.D0
    T=0.D0
    H=0.D0
    G=0.D0
    CONT1=0
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN     
            SELECT CASE (SURFACE_TYPE(I))
            CASE("B")
                IF (DESCONT_TYPE(I).EQ."DB") THEN
                    DO J=1,4*SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)
                        CONT1=CONT1+1
                    ENDDO  
                ELSE
                    CONT1=CONT1+(SUBDIVISION_PD1(I)+1)*(SUBDIVISION_PD2(I)+1)
                ENDIF
            CASE("CD")
                DO J=1,4*SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)
                    CONT1=CONT1+1
                    DUAL_BEM(CONT1)="S"
                    EQ_TYPE(CONT1)="S"
                ENDDO    
            CASE("CT")
                DO J=1,4*SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)
                    CONT1=CONT1+1
                    DUAL_BEM(CONT1)="H"
                    EQ_TYPE(CONT1)="H"
                ENDDO  
            CASE("CI")
                DO J=1,4*SUBDIVISION_PD1(I)*SUBDIVISION_PD2(I)
                    CONT1=CONT1+1
                ENDDO  
            ENDSELECT
        ENDIF        
    ENDDO   
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   COLLOCATION POINTS GENERATION: CONTINUOUS, DISCONTINUOUS AND EDGE DISCONTINUOUS ELEMENTS   
!------------------------------------------------------------------------------------------------------------------------------------------------------
    CONT1=0
    CONT2=0
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN     
            SELECT CASE(SURFACE_TYPE(I))
            CASE("B")
                IF (DESCONT_TYPE(I).EQ."DB") THEN
                    DO J=1,SUBDIVISION_PD2(I)
                        DO K=1,SUBDIVISION_PD1(I)
                            CONT1=CONT1+1
                            DO L=1,4
                                CONT2=CONT2+1
                                COLLOCPOINTS_CONNECTIVITY(CONT1,L)=CONT2
                            ENDDO                                                          
                        ENDDO
                    ENDDO 
                ELSE
                    DO J=1,SUBDIVISION_PD2(I)
                        DO K=1,SUBDIVISION_PD1(I)
                            CONT1=CONT1+1
                            COLLOCPOINTS_CONNECTIVITY(CONT1,1)=CONT2+J*(SUBDIVISION_PD1(I)+1)+K-(SUBDIVISION_PD1(I)+1)
                            COLLOCPOINTS_CONNECTIVITY(CONT1,2)=CONT2+J*(SUBDIVISION_PD1(I)+1)+K-(SUBDIVISION_PD1(I)+1)+1
                            COLLOCPOINTS_CONNECTIVITY(CONT1,3)=CONT2+J*(SUBDIVISION_PD1(I)+1)+K+1   
                            COLLOCPOINTS_CONNECTIVITY(CONT1,4)=CONT2+J*(SUBDIVISION_PD1(I)+1)+K                                                           
                        ENDDO
                    ENDDO
                    CONT2=CONT2+(SUBDIVISION_PD1(I)+1)*(SUBDIVISION_PD2(I)+1)
                ENDIF
            CASE("CD")
                DO J=1,SUBDIVISION_PD2(I)
                    DO K=1,SUBDIVISION_PD1(I)
                        CONT1=CONT1+1
                        DO L=1,4
                            CONT2=CONT2+1
                            COLLOCPOINTS_CONNECTIVITY(CONT1,L)=CONT2
                        ENDDO                                                          
                    ENDDO
                ENDDO        
            CASE("CT")        
                DO J=1,SUBDIVISION_PD2(I)
                    DO K=1,SUBDIVISION_PD1(I)
                        CONT1=CONT1+1
                        DO L=1,4
                            CONT2=CONT2+1
                            COLLOCPOINTS_CONNECTIVITY(CONT1,L)=CONT2
                        ENDDO                                                          
                    ENDDO
                ENDDO   
            CASE("CI")        
                DO J=1,SUBDIVISION_PD2(I)
                    DO K=1,SUBDIVISION_PD1(I)
                        CONT1=CONT1+1
                        DO L=1,4
                            CONT2=CONT2+1
                            COLLOCPOINTS_CONNECTIVITY(CONT1,L)=CONT2
                        ENDDO                                                          
                    ENDDO
                ENDDO 
            ENDSELECT
        ENDIF      
    ENDDO        
    CONT1=0
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN     
            SELECT CASE(SURFACE_TYPE(I))
            CASE("B")
                IF (DESCONT_TYPE(I).EQ."DB") THEN
                    DO J=1,SUBDIVISION_PD2(I)
                        DO K=1,SUBDIVISION_PD1(I)
                            CONT1=CONT1+1 
                            ND(CONT1)=DOMAINS(I)       
                            DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,1))="DISCONTINUOUS"
                            DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,2))="DISCONTINUOUS"
                            DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,3))="DISCONTINUOUS"
                            DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,4))="DISCONTINUOUS"
                        ENDDO
                    ENDDO 
                ELSE                  
                    DO J=1,SUBDIVISION_PD2(I)
                        DO K=1,SUBDIVISION_PD1(I)
                            CONT1=CONT1+1
                            ND(CONT1)=DOMAINS(I)  
                            IF (J.EQ.1)THEN
                                EDGE_DISCONTINUITY(CONT1,1)=4
                                EDGE_DISCONTINUITY(CONT1,2)=2                            
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,1))="EDGEDISCONTINUOUS"
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,2))="EDGEDISCONTINUOUS"
                                IF(K.EQ.1)THEN
                                    EDGE_DISCONTINUITY(CONT1,2)=2
                                    EDGE_DISCONTINUITY(CONT1,4)=3                                
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,1))="DISCONTINUOUS"
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,2))="EDGEDISCONTINUOUS"
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,4))="EDGEDISCONTINUOUS"
                                ELSE IF(K.EQ.SUBDIVISION_PD1(I))THEN
                                    EDGE_DISCONTINUITY(CONT1,1)=4
                                    EDGE_DISCONTINUITY(CONT1,3)=3                                
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,1))="EDGEDISCONTINUOUS"
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,2))="DISCONTINUOUS"
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,3))="EDGEDISCONTINUOUS"
                                ENDIF
                            ELSE IF(J.EQ.SUBDIVISION_PD2(I))THEN
                                EDGE_DISCONTINUITY(CONT1,3)=2
                                EDGE_DISCONTINUITY(CONT1,4)=4                            
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,3))="EDGEDISCONTINUOUS"
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,4))="EDGEDISCONTINUOUS"
                                IF(K.EQ.1)THEN
                                    EDGE_DISCONTINUITY(CONT1,1)=1
                                    EDGE_DISCONTINUITY(CONT1,3)=2                                
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,1))="EDGEDISCONTINUOUS"
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,3))="EDGEDISCONTINUOUS"
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,4))="DISCONTINUOUS"
                                ELSE IF(K.EQ.SUBDIVISION_PD1(I))THEN
                                    EDGE_DISCONTINUITY(CONT1,2)=1
                                    EDGE_DISCONTINUITY(CONT1,4)=4                                
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,2))="EDGEDISCONTINUOUS"
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,3))="DISCONTINUOUS"
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,4))="EDGEDISCONTINUOUS"
                                ENDIF
                            ELSE
                                IF(K.EQ.1)THEN
                                    EDGE_DISCONTINUITY(CONT1,1)=1
                                    EDGE_DISCONTINUITY(CONT1,4)=3                                
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,1))="EDGEDISCONTINUOUS"
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,4))="EDGEDISCONTINUOUS"
                                ELSE IF(K.EQ.SUBDIVISION_PD1(I))THEN
                                    EDGE_DISCONTINUITY(CONT1,2)=1
                                    EDGE_DISCONTINUITY(CONT1,3)=3                                
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,2))="EDGEDISCONTINUOUS"
                                    DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,3))="EDGEDISCONTINUOUS"                        
                                ENDIF                                
                            ENDIF
                        ENDDO
                    ENDDO 
                ENDIF 
            CASE("CD")
                DO J=1,SUBDIVISION_PD2(I)
                    DO K=1,SUBDIVISION_PD1(I)
                        CONT1=CONT1+1 
                        ND(CONT1)=DOMAINS(I)       
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,1))="DISCONTINUOUS"
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,2))="DISCONTINUOUS"
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,3))="DISCONTINUOUS"
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,4))="DISCONTINUOUS"
                    ENDDO
                ENDDO        
            CASE("CT")
                 DO J=1,SUBDIVISION_PD2(I)
                    DO K=1,SUBDIVISION_PD1(I)
                        CONT1=CONT1+1
                        ND(CONT1)=DOMAINS(I)                             
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,1))="DISCONTINUOUS"
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,2))="DISCONTINUOUS"
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,3))="DISCONTINUOUS"
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,4))="DISCONTINUOUS"
                    ENDDO
                 ENDDO 
            CASE("CI")
                DO J=1,SUBDIVISION_PD2(I)
                    DO K=1,SUBDIVISION_PD1(I)
                        CONT1=CONT1+1 
                        ND(CONT1)=DOMAINS(I)       
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,1))="DISCONTINUOUS"
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,2))="DISCONTINUOUS"
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,3))="DISCONTINUOUS"
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(CONT1,4))="DISCONTINUOUS"
                    ENDDO
                ENDDO 
            ENDSELECT  
        ENDIF      
    ENDDO   
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   DEFINING CRACKTIP COLLOCATION POINTS   
!------------------------------------------------------------------------------------------------------------------------------------------------------  
    CONT1=0
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN     
            SELECT CASE(SURFACE_TYPE(I))
            CASE("B")
                CONT1=CONT1+(SUBDIVISION_PD1(I))*(SUBDIVISION_PD2(I))
            ! CAUTION: PODE DAR PROBLEMA DEPOIS
            CASE("CI")
                CONT1=CONT1+(SUBDIVISION_PD1(I))*(SUBDIVISION_PD2(I))
            CASE("CD")
                DO J=1,SUBDIVISION_PD2(I)
                    DO K=1,SUBDIVISION_PD1(I)
                        CONT1=CONT1+1
                        IF(J.EQ.1.AND.K.EQ.1)THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,1))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,2))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,4))=1                                                
                        ELSE IF(J.EQ.1.AND.K.EQ.SUBDIVISION_PD1(I))THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,1))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,2))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,3))=1     
                        ELSE IF(J.EQ.SUBDIVISION_PD2(I).AND.K.EQ.1)THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,1))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,3))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,4))=1     
                        ELSE IF(J.EQ.SUBDIVISION_PD2(I).AND.K.EQ.SUBDIVISION_PD1(I))THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,2))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,3))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,4))=1  
                        ELSE IF(J.EQ.1.AND.K.NE.1.AND.K.NE.SUBDIVISION_PD1(I))THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,1))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,2))=1     
                        ELSE IF(J.EQ.SUBDIVISION_PD2(I).AND.K.NE.1.AND.K.NE.SUBDIVISION_PD1(I))THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,3))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,4))=1          
                        ELSE IF(K.EQ.1.AND.J.NE.1.AND.J.NE.SUBDIVISION_PD2(I))THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,1))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,4))=1     
                        ELSE IF(K.EQ.SUBDIVISION_PD1(I).AND.J.NE.1.AND.J.NE.SUBDIVISION_PD2(I))THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,2))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,3))=1                                                             
                        ENDIF                                                               
                    ENDDO
                ENDDO
            CASE("CT")
                DO J=1,SUBDIVISION_PD2(I)
                    DO K=1,SUBDIVISION_PD1(I)
                        CONT1=CONT1+1
                        IF(J.EQ.1.AND.K.EQ.1)THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,1))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,2))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,4))=1                                                
                        ELSE IF(J.EQ.1.AND.K.EQ.SUBDIVISION_PD1(I))THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,1))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,2))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,3))=1     
                        ELSE IF(J.EQ.SUBDIVISION_PD2(I).AND.K.EQ.1)THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,1))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,3))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,4))=1     
                        ELSE IF(J.EQ.SUBDIVISION_PD2(I).AND.K.EQ.SUBDIVISION_PD1(I))THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,2))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,3))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,4))=1  
                        ELSE IF(J.EQ.1.AND.K.NE.1.AND.K.NE.SUBDIVISION_PD1(I))THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,1))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,2))=1     
                        ELSE IF(J.EQ.SUBDIVISION_PD2(I).AND.K.NE.1.AND.K.NE.SUBDIVISION_PD1(I))THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,3))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,4))=1          
                        ELSE IF(K.EQ.1.AND.J.NE.1.AND.J.NE.SUBDIVISION_PD2(I))THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,1))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,4))=1     
                        ELSE IF(K.EQ.SUBDIVISION_PD1(I).AND.J.NE.1.AND.J.NE.SUBDIVISION_PD2(I))THEN
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,2))=1
                            CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(CONT1,3))=1                                                             
                        ENDIF   
                    ENDDO
                ENDDO        
            ENDSELECT
        ENDIF    
    ENDDO
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   MESH GENERATION: TRIMMED NURBS SURFACES   
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   COHESIVE INTERFACE CAN ONLY BE UNTRIMMED
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).NE."UNTRIMMED")THEN 
!            
            PARAMETRIC_VERTICES(1,1)=KNOT_VECTORS(I,1,1)
            PARAMETRIC_VERTICES(1,2)=KNOT_VECTORS(I,2,1)
            PARAMETRIC_VERTICES(2,1)=KNOT_VECTORS(I,1,N_BF(I,1)+POLYNOMIAL_ORDERS(I,1)+1)
            PARAMETRIC_VERTICES(2,2)=KNOT_VECTORS(I,2,1)            
            PARAMETRIC_VERTICES(3,1)=KNOT_VECTORS(I,1,N_BF(I,1)+POLYNOMIAL_ORDERS(I,1)+1)
            PARAMETRIC_VERTICES(3,2)=KNOT_VECTORS(I,2,N_BF(I,2)+POLYNOMIAL_ORDERS(I,2)+1)
            PARAMETRIC_VERTICES(4,1)=KNOT_VECTORS(I,1,1)
            PARAMETRIC_VERTICES(4,2)=KNOT_VECTORS(I,2,N_BF(I,2)+POLYNOMIAL_ORDERS(I,2)+1)
!
            N_CURVEPOINTS=0
            DO J=1,N_TRIMMINGCURVES
                IF(SURFACES(J).EQ.I)THEN
                    N_CURVEPOINTS=N_CURVEPOINTS+SUBDIVISION(J)+1	
                ENDIF
            ENDDO
            ALLOCATE(TRIMMING_CURVEPOINTS(N_CURVEPOINTS,2))
            TRIMMING_CURVEPOINTS=0.D0
            N_CURVEPOINTS=0    
	        DO J=1,N_TRIMMINGCURVES
                IF(SURFACES(J).EQ.I)THEN
                    ALLOCATE(VALUES(N_BFCURVES(J),2),PHI(N_BFCURVES(J)))            
                    DO K=1,N_BFCURVES(J)                
                        VALUES(K,1)=COORD_PARAMETRICCONTROLPOINTS(PARAMETRICCONTROLPOINTS_CONNECTIVITY(J,K),1)
                        VALUES(K,2)=COORD_PARAMETRICCONTROLPOINTS(PARAMETRICCONTROLPOINTS_CONNECTIVITY(J,K),2)                                            
                    ENDDO                 	 
	                DELTAA=(KNOT_VECTOR(J,(N_BFCURVES(J)+POLYNOMIAL_ORDER(J)+1))-KNOT_VECTOR(J,1))/SUBDIVISION(J)
	                A=KNOT_VECTOR(J,1)
	                DO K=1,SUBDIVISION(J)
	                    N_CURVEPOINTS=N_CURVEPOINTS+1
	                    CALL PARAMETRICAL_MAPPING(J,A,VALUES,PHI,QSI1,QSI2)
	                    TRIMMING_CURVEPOINTS(N_CURVEPOINTS,1)=QSI1
	                    TRIMMING_CURVEPOINTS(N_CURVEPOINTS,2)=QSI2
	                    A=A+DELTAA
	                ENDDO
	                A=KNOT_VECTOR(J,(N_BFCURVES(J)+POLYNOMIAL_ORDER(J)+1))
	                N_CURVEPOINTS=N_CURVEPOINTS+1
	                CALL PARAMETRICAL_MAPPING(J,A,VALUES,PHI,QSI1,QSI2)
	                TRIMMING_CURVEPOINTS(N_CURVEPOINTS,1)=QSI1
	                TRIMMING_CURVEPOINTS(N_CURVEPOINTS,2)=QSI2	           
	                DEALLOCATE(VALUES,PHI)	         
	            ENDIF    	    
	        ENDDO             
!                       
            IF(TRIMMED_SURFACE(I).EQ."TRIMMEDINNER")THEN
                CALL WRITE_TRIMMEDINNERFILE(I,PARAMETRIC_VERTICES,N_CURVEPOINTS,TRIMMING_CURVEPOINTS,SUBDIVISION_PD1(I))  
            ENDIF        
            IF(TRIMMED_SURFACE(I).EQ."TRIMMEDOUTER")THEN
                CALL WRITE_TRIMMEDOUTERFILE(I,N_CURVEPOINTS,TRIMMING_CURVEPOINTS,SUBDIVISION_PD1(I))                
            ENDIF     
!
            ANSYS_IN  = 'P:\filename.log'
            ANSYS_OUT  = 'P:\results.out'             
            OPEN(5,FILE='Input_data\ANSYS_VER.txt')
            READ(5,*)ANSYS_VER
            CLOSE(5)
            WRITE(CMD1,"(I0)")ANSYS_VER
            II = SYSTEM('"'//'"%ANSYS'//CMD1//'_DIR%\bin\winx64\ANSYS'//CMD1//'.exe" -B -I '//ANSYS_IN//' -O '//ANSYS_OUT//'"')
            !PC ROCHA
            !II = SYSTEM('"C:\Program Files\ANSYS Inc\v191\ansys\bin\winx64\ANSYS191.exe" -B -I '//ANSYS_IN//' -O '//ANSYS_OUT)
            ! 
	        OPEN(14,file='parametric_mesh.txt',status='unknown')
	        KEEP_READING=.TRUE.
	        DO WHILE(KEEP_READING)
	            READ(14,10)KEYWORD
	            IF(TRIM(KEYWORD).EQ.'*ELSE')THEN
	                KEEP_READING=.FALSE.
	            ENDIF
	        ENDDO
	        READ(14,20)N_PARAMETRICNODES
	        READ(14,20)N_PARAMETRICELEMENTS
	        N_ELEM_SURFACES(I)=N_PARAMETRICELEMENTS
	        ALLOCATE(PARAMETRICCOORD_NODES(N_PARAMETRICNODES,2),PARAMETRIC_CONNECTIVITY(N_PARAMETRICELEMENTS,4),SIGNED_DIST_PARAMETRICNODES(N_PARAMETRICNODES),COMMON_ELEMENTS(N_PARAMETRICNODES),REPEATED_CP(N_PARAMETRICNODES))
	        PARAMETRICCOORD_NODES=0.D0
	        PARAMETRIC_CONNECTIVITY=0
	        COMMON_ELEMENTS=0	        
	        KEEP_READING=.TRUE.
	        DO WHILE(KEEP_READING)
	            READ(14,10)KEYWORD
	            IF(TRIM(KEYWORD).EQ.'NBLOCK')THEN
	                KEEP_READING=.FALSE.
	                READ(14,*)
	            ENDIF
	        ENDDO
!	        
	        DO J=1,N_PARAMETRICNODES
	            READ(14,30)PARAMETRICCOORD_NODES(J,1),PARAMETRICCOORD_NODES(J,2)   
	        ENDDO
	        READ(14,*)
	        READ(14,*)
	        READ(14,*)	        	        
	        DO J=1,N_PARAMETRICELEMENTS
	            READ(14,*)II,II,II,II,II,II,II,II,II,II,II,PARAMETRIC_CONNECTIVITY(J,1),PARAMETRIC_CONNECTIVITY(J,2),PARAMETRIC_CONNECTIVITY(J,3),PARAMETRIC_CONNECTIVITY(J,4)
	        ENDDO	        
!
	        CLOSE(14)
!
            CALL BOUNDARY_PARAMETRICNODES(I,N_PARAMETRICNODES,N_PARAMETRICELEMENTS,PARAMETRIC_CONNECTIVITY,PARAMETRICCOORD_NODES,SIGNED_DIST_PARAMETRICNODES,N_CURVEPOINTS,TRIMMING_CURVEPOINTS,COMMON_ELEMENTS)	        
!
            IF(SURFACE_TYPE(I).EQ."B")THEN
                N_PARAMETRICCOLLOCPOINTS=N_PARAMETRICNODES
	            DO J=1,N_PARAMETRICNODES
	                IF(COMMON_ELEMENTS(J).GE.3)THEN
	                    N_PARAMETRICCOLLOCPOINTS=N_PARAMETRICCOLLOCPOINTS+COMMON_ELEMENTS(J)-1
	                ENDIF
	            ENDDO
                N_COLLOCPOINTS_SURFACES(I)=N_PARAMETRICCOLLOCPOINTS
            ELSE
                CONT1=0       
                DO J=1,N_PARAMETRICELEMENTS
                    IF(PARAMETRIC_CONNECTIVITY(J,4).EQ.0.OR.PARAMETRIC_CONNECTIVITY(J,3).EQ.PARAMETRIC_CONNECTIVITY(J,4))THEN
                        CONT1=CONT1+3
                    ELSE
                        CONT1=CONT1+4                
                    ENDIF            
                ENDDO
                N_PARAMETRICCOLLOCPOINTS=CONT1
                N_COLLOCPOINTS_SURFACES(I)=N_PARAMETRICCOLLOCPOINTS
            ENDIF       
!	        	                  
            DEALLOCATE(TRIMMING_CURVEPOINTS)
!------------------------------------------------------------------------------------------------------------------------------------------------------
!           COMPUTING THE PHYSICAL COORDINATES OF THE NODES
!------------------------------------------------------------------------------------------------------------------------------------------------------
            ALLOCATE(COORD_NODES_AUX(N_NODES,3))
            COORD_NODES_AUX=COORD_NODES
            DEALLOCATE(COORD_NODES)
            ALLOCATE(COORD_NODES(N_NODES+N_PARAMETRICNODES,3))
            DO J=1,N_NODES
                COORD_NODES(J,:)=COORD_NODES_AUX(J,:)
            ENDDO
            ALLOCATE(VALUES(N_GBF(I),3),PHI(N_GBF(I)))            
            DO J=1,N_GBF(I)                
                VALUES(J,1)=COORD_CONTROLPOINTS(CONTROLPOINTS_CONNECTIVITY(I,J),1)
                VALUES(J,2)=COORD_CONTROLPOINTS(CONTROLPOINTS_CONNECTIVITY(I,J),2)
                VALUES(J,3)=COORD_CONTROLPOINTS(CONTROLPOINTS_CONNECTIVITY(I,J),3)                                            
            ENDDO                        
            DO J=1,N_PARAMETRICNODES
                CALL PHYSICAL_MAPPING(I,PARAMETRICCOORD_NODES(J,1),PARAMETRICCOORD_NODES(J,2),VALUES,PHI,X1,X2,X3)
                COORD_NODES(N_NODES+J,1)=X1
                COORD_NODES(N_NODES+J,2)=X2
                COORD_NODES(N_NODES+J,3)=X3
            ENDDO
            DEALLOCATE(VALUES,PHI) 
            DEALLOCATE(COORD_NODES_AUX)                    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!           DEFINING THE NEW ELEMENTS
!------------------------------------------------------------------------------------------------------------------------------------------------------    
            ALLOCATE(ELEM_TYPE_AUX(N_ELEM),ORDER_ELEM_AUX(N_ELEM),POL_FAMILY_AUX(N_ELEM),NODES_CONNECTIVITY_AUX(N_ELEM,4),COLLOCPOINTS_CONNECTIVITY_AUX(N_ELEM,4),QSI_AUX(N_ELEM,9,2),ND_AUX(N_ELEM),EDGE_DISCONTINUITY_AUX(N_ELEM,4))
            ELEM_TYPE_AUX=ELEM_TYPE
            ORDER_ELEM_AUX=ORDER_ELEM
            POL_FAMILY_AUX=POL_FAMILY
            NODES_CONNECTIVITY_AUX=NODES_CONNECTIVITY 
            COLLOCPOINTS_CONNECTIVITY_AUX=COLLOCPOINTS_CONNECTIVITY 
            QSI_AUX=QSI 
            ND_AUX=ND
            EDGE_DISCONTINUITY_AUX=EDGE_DISCONTINUITY
            DEALLOCATE(ELEM_TYPE,ORDER_ELEM,POL_FAMILY,NODES_CONNECTIVITY,COLLOCPOINTS_CONNECTIVITY,QSI,ND,EDGE_DISCONTINUITY)
            ALLOCATE(ELEM_TYPE(N_ELEM+N_PARAMETRICELEMENTS),ORDER_ELEM(N_ELEM+N_PARAMETRICELEMENTS),POL_FAMILY(N_ELEM+N_PARAMETRICELEMENTS),NODES_CONNECTIVITY(N_ELEM+N_PARAMETRICELEMENTS,4),COLLOCPOINTS_CONNECTIVITY(N_ELEM+N_PARAMETRICELEMENTS,4),QSI(N_ELEM+N_PARAMETRICELEMENTS,9,2),ND(N_ELEM+N_PARAMETRICELEMENTS),EDGE_DISCONTINUITY(N_ELEM+N_PARAMETRICELEMENTS,4))                                                                    
            ORDER_ELEM=1
            POL_FAMILY=1
            EDGE_DISCONTINUITY=0
            ND=DOMAINS(I)
            DO J=1,N_ELEM
                ELEM_TYPE(J)=ELEM_TYPE_AUX(J)
                ORDER_ELEM(J)=ORDER_ELEM_AUX(J)
                POL_FAMILY(J)=POL_FAMILY_AUX(J)
                NODES_CONNECTIVITY(J,1)=NODES_CONNECTIVITY_AUX(J,1)
                NODES_CONNECTIVITY(J,2)=NODES_CONNECTIVITY_AUX(J,2)  
                NODES_CONNECTIVITY(J,3)=NODES_CONNECTIVITY_AUX(J,3)  
                NODES_CONNECTIVITY(J,4)=NODES_CONNECTIVITY_AUX(J,4)                                               
                COLLOCPOINTS_CONNECTIVITY(J,1)=COLLOCPOINTS_CONNECTIVITY_AUX(J,1)
                COLLOCPOINTS_CONNECTIVITY(J,2)=COLLOCPOINTS_CONNECTIVITY_AUX(J,2)   
                COLLOCPOINTS_CONNECTIVITY(J,3)=COLLOCPOINTS_CONNECTIVITY_AUX(J,3)
                COLLOCPOINTS_CONNECTIVITY(J,4)=COLLOCPOINTS_CONNECTIVITY_AUX(J,4)                              
                ND(J)=ND_AUX(J) 
                EDGE_DISCONTINUITY(J,:)=EDGE_DISCONTINUITY_AUX(J,:)               
            ENDDO
            DEALLOCATE(ELEM_TYPE_AUX,ORDER_ELEM_AUX,POL_FAMILY_AUX,NODES_CONNECTIVITY_AUX,COLLOCPOINTS_CONNECTIVITY_AUX,QSI_AUX,ND_AUX,EDGE_DISCONTINUITY_AUX)            
!
            ALLOCATE(EQ_TYPE_AUX(N_COLLOCPOINTS),DUAL_BEM_AUX(N_COLLOCPOINTS),DISCONTINUOUS_COLLOCPOINTS_AUX(N_COLLOCPOINTS),CRACK_TIP_AUX(N_COLLOCPOINTS)) 
            EQ_TYPE_AUX=EQ_TYPE
            DUAL_BEM_AUX=DUAL_BEM
            DISCONTINUOUS_COLLOCPOINTS_AUX=DISCONTINUOUS_COLLOCPOINTS
            CRACK_TIP_AUX=CRACK_TIP
            DEALLOCATE(COORD_COLLOCPOINTS,B_CONDITIONS,U,T,EQ_TYPE,DUAL_BEM,H,G,DISCONTINUOUS_COLLOCPOINTS,CRACK_TIP)
            ALLOCATE(COORD_COLLOCPOINTS((N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS),3),B_CONDITIONS(3*(N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS)),U(3*(N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS)),T(3*(N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS)),EQ_TYPE(N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS),DUAL_BEM(N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS),H(3*(N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS),3*(N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS)),G(3*(N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS),3*(N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS)),DISCONTINUOUS_COLLOCPOINTS(N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS),CRACK_TIP(N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS)) 
            DISCONTINUOUS_COLLOCPOINTS="CONTINUOUS"
            DO J=1,N_COLLOCPOINTS
                EQ_TYPE(J)=EQ_TYPE_AUX(J)
                DUAL_BEM(J)=DUAL_BEM_AUX(J)
                DISCONTINUOUS_COLLOCPOINTS(J)=DISCONTINUOUS_COLLOCPOINTS_AUX(J)
                CRACK_TIP(J)=CRACK_TIP_AUX(J)
            ENDDO
            DEALLOCATE(EQ_TYPE_AUX,DUAL_BEM_AUX,DISCONTINUOUS_COLLOCPOINTS_AUX,CRACK_TIP_AUX)
!                        
            CONT1=0            
            REPEATED_CP=0
!
            SELECT CASE(SURFACE_TYPE(I))
            CASE("B")
                DO J=1,N_PARAMETRICCOLLOCPOINTS
                    DUAL_BEM(N_COLLOCPOINTS+J)="B"                
                    EQ_TYPE(N_COLLOCPOINTS+J)="S"
                    CRACK_TIP(N_COLLOCPOINTS+J)=0
                    DISCONTINUOUS_COLLOCPOINTS(N_COLLOCPOINTS+J)="CONTINUOUS"  
                ENDDO
                DO J=1,N_PARAMETRICELEMENTS
                    EDGE_DISCONTINUITY(N_ELEM+J,:)=0
                    IF(PARAMETRIC_CONNECTIVITY(J,4).EQ.0.OR.PARAMETRIC_CONNECTIVITY(J,3).EQ.PARAMETRIC_CONNECTIVITY(J,4))THEN
                    !--------TRIANGULAR ELEMENTS-------------                
                        ELEM_TYPE(N_ELEM+J)=3
                        NODES_CONNECTIVITY(N_ELEM+J,1)=N_NODES+PARAMETRIC_CONNECTIVITY(J,1)  
                        NODES_CONNECTIVITY(N_ELEM+J,2)=N_NODES+PARAMETRIC_CONNECTIVITY(J,2)   
                        NODES_CONNECTIVITY(N_ELEM+J,3)=N_NODES+PARAMETRIC_CONNECTIVITY(J,3)
!
                        IF (COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,1)).LT.3)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,1)                                    
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,1)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,1)).EQ.0)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,1)  
                            REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,1))=1
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,1)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,1)).EQ.1)THEN
                            CONT1=CONT1+1
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)=N_COLLOCPOINTS+N_PARAMETRICNODES+CONT1                             
                        ENDIF
!
                        IF (COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,2)).LT.3)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,2)                                    
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,2)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,2)).EQ.0)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,2)  
                            REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,2))=1
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,2)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,2)).EQ.1)THEN
                            CONT1=CONT1+1
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)=N_COLLOCPOINTS+N_PARAMETRICNODES+CONT1                             
                        ENDIF 
!
                        IF (COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,3)).LT.3)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,3)                                    
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,3)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,3)).EQ.0)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,3)  
                            REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,3))=1
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,3)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,3)).EQ.1)THEN
                            CONT1=CONT1+1
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)=N_COLLOCPOINTS+N_PARAMETRICNODES+CONT1                             
                        ENDIF                                                                              
!                                 
                    ELSE
                    !--------QUADRILATERAL ELEMENTS-------------                     
                        ELEM_TYPE(N_ELEM+J)=4
                        NODES_CONNECTIVITY(N_ELEM+J,1)=N_NODES+PARAMETRIC_CONNECTIVITY(J,1)  
                        NODES_CONNECTIVITY(N_ELEM+J,2)=N_NODES+PARAMETRIC_CONNECTIVITY(J,2)   
                        NODES_CONNECTIVITY(N_ELEM+J,3)=N_NODES+PARAMETRIC_CONNECTIVITY(J,3)    
                        NODES_CONNECTIVITY(N_ELEM+J,4)=N_NODES+PARAMETRIC_CONNECTIVITY(J,4)
!
                        IF (COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,1)).LT.3)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,1)                                    
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,1)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,1)).EQ.0)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,1)  
                            REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,1))=1
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,1)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,1)).EQ.1)THEN
                            CONT1=CONT1+1
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)=N_COLLOCPOINTS+N_PARAMETRICNODES+CONT1                             
                        ENDIF
!
                        IF (COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,2)).LT.3)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,2)                                    
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,2)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,2)).EQ.0)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,2)  
                            REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,2))=1
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,2)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,2)).EQ.1)THEN
                            CONT1=CONT1+1
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)=N_COLLOCPOINTS+N_PARAMETRICNODES+CONT1                             
                        ENDIF 
!
                        IF (COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,3)).LT.3)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,3)                                    
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,3)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,3)).EQ.0)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,3)  
                            REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,3))=1
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,3)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,3)).EQ.1)THEN
                            CONT1=CONT1+1
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)=N_COLLOCPOINTS+N_PARAMETRICNODES+CONT1                             
                        ENDIF
!
                        IF (COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,4)).LT.3)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,4)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,4)                                    
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,4)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,4)).EQ.0)THEN
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,4)=N_COLLOCPOINTS+PARAMETRIC_CONNECTIVITY(J,4)  
                            REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,4))=1
                        ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,4)).GE.3.AND.REPEATED_CP(PARAMETRIC_CONNECTIVITY(J,4)).EQ.1)THEN
                            CONT1=CONT1+1
                            COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,4)=N_COLLOCPOINTS+N_PARAMETRICNODES+CONT1                             
                        ENDIF                                                                          
!
                    ENDIF
                ENDDO                  
            CASE("CD")
                DO J=1,N_PARAMETRICCOLLOCPOINTS
                    DUAL_BEM(N_COLLOCPOINTS+J)="S"                
                    EQ_TYPE(N_COLLOCPOINTS+J)="S"
                    CRACK_TIP(N_COLLOCPOINTS+J)=0
                    DISCONTINUOUS_COLLOCPOINTS(N_COLLOCPOINTS+J)="DISCONTINUOUS" 
                ENDDO            
                DO J=1,N_PARAMETRICELEMENTS
                    EDGE_DISCONTINUITY(N_ELEM+J,:)=0                  
                    IF(PARAMETRIC_CONNECTIVITY(J,4).EQ.0.OR.PARAMETRIC_CONNECTIVITY(J,3).EQ.PARAMETRIC_CONNECTIVITY(J,4))THEN
                    !--------TRIANGULAR ELEMENTS------------- 
                        ELEM_TYPE(N_ELEM+J)=3
                        NODES_CONNECTIVITY(N_ELEM+J,1)=N_NODES+PARAMETRIC_CONNECTIVITY(J,1)  
                        NODES_CONNECTIVITY(N_ELEM+J,2)=N_NODES+PARAMETRIC_CONNECTIVITY(J,2)   
                        NODES_CONNECTIVITY(N_ELEM+J,3)=N_NODES+PARAMETRIC_CONNECTIVITY(J,3) 
!                     
                        CONT1=CONT1+3
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)=N_COLLOCPOINTS+CONT1-2  
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)=N_COLLOCPOINTS+CONT1-1   
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)=N_COLLOCPOINTS+CONT1 
                    ELSE
                    !--------QUADRILATERAL ELEMENTS-------------                     
                        ELEM_TYPE(N_ELEM+J)=4
                        NODES_CONNECTIVITY(N_ELEM+J,1)=N_NODES+PARAMETRIC_CONNECTIVITY(J,1)  
                        NODES_CONNECTIVITY(N_ELEM+J,2)=N_NODES+PARAMETRIC_CONNECTIVITY(J,2)   
                        NODES_CONNECTIVITY(N_ELEM+J,3)=N_NODES+PARAMETRIC_CONNECTIVITY(J,3)    
                        NODES_CONNECTIVITY(N_ELEM+J,4)=N_NODES+PARAMETRIC_CONNECTIVITY(J,4)   
!
                        CONT1=CONT1+4
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)=N_COLLOCPOINTS+CONT1-3  
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)=N_COLLOCPOINTS+CONT1-2   
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)=N_COLLOCPOINTS+CONT1-1    
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,4)=N_COLLOCPOINTS+CONT1                                               
                    ENDIF
                ENDDO                                                      
            CASE("CT")
                DO J=1,N_PARAMETRICCOLLOCPOINTS
                    DUAL_BEM(N_COLLOCPOINTS+J)="H"                
                    EQ_TYPE(N_COLLOCPOINTS+J)="H"
                    CRACK_TIP(N_COLLOCPOINTS+J)=0
                    DISCONTINUOUS_COLLOCPOINTS(N_COLLOCPOINTS+J)="DISCONTINUOUS" 
                ENDDO             
                DO J=1,N_PARAMETRICELEMENTS
                    EDGE_DISCONTINUITY(N_ELEM+J,:)=0                  
                    IF(PARAMETRIC_CONNECTIVITY(J,4).EQ.0.OR.PARAMETRIC_CONNECTIVITY(J,3).EQ.PARAMETRIC_CONNECTIVITY(J,4))THEN
                    !--------TRIANGULAR ELEMENTS------------- 
                        ELEM_TYPE(N_ELEM+J)=3
                        NODES_CONNECTIVITY(N_ELEM+J,1)=N_NODES+PARAMETRIC_CONNECTIVITY(J,1)  
                        NODES_CONNECTIVITY(N_ELEM+J,2)=N_NODES+PARAMETRIC_CONNECTIVITY(J,2)   
                        NODES_CONNECTIVITY(N_ELEM+J,3)=N_NODES+PARAMETRIC_CONNECTIVITY(J,3) 
!                                        
                        CONT1=CONT1+3
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)=N_COLLOCPOINTS+CONT1-2  
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)=N_COLLOCPOINTS+CONT1-1   
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)=N_COLLOCPOINTS+CONT1                         
                    ELSE
                    !--------QUADRILATERAL ELEMENTS------------- 
                        ELEM_TYPE(N_ELEM+J)=4
                        NODES_CONNECTIVITY(N_ELEM+J,1)=N_NODES+PARAMETRIC_CONNECTIVITY(J,1)  
                        NODES_CONNECTIVITY(N_ELEM+J,2)=N_NODES+PARAMETRIC_CONNECTIVITY(J,2)   
                        NODES_CONNECTIVITY(N_ELEM+J,3)=N_NODES+PARAMETRIC_CONNECTIVITY(J,3)    
                        NODES_CONNECTIVITY(N_ELEM+J,4)=N_NODES+PARAMETRIC_CONNECTIVITY(J,4)   
!
                        CONT1=CONT1+4
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)=N_COLLOCPOINTS+CONT1-3  
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)=N_COLLOCPOINTS+CONT1-2   
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)=N_COLLOCPOINTS+CONT1-1    
                        COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,4)=N_COLLOCPOINTS+CONT1                                                                      
                    ENDIF
                ENDDO             
            ENDSELECT
!
            SELECT CASE(SURFACE_TYPE(I))
            CASE("B")
                DO J=1,N_PARAMETRICELEMENTS
                    IF(PARAMETRIC_CONNECTIVITY(J,4).EQ.0.OR.PARAMETRIC_CONNECTIVITY(J,3).EQ.PARAMETRIC_CONNECTIVITY(J,4))THEN 
                    !--------TRIANGULAR ELEMENTS------------- 
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,1)).EQ.0.D0)THEN
                            IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,1)).GE.3)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1))="DISCONTINUOUS"
                            ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,1)).EQ.1)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1))="DISCONTINUOUS"                            
                            ELSE
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1))="EDGEDISCONTINUOUS"
                            ENDIF
                        ENDIF    
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,2)).EQ.0.D0)THEN
                            IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,2)).GE.3)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2))="DISCONTINUOUS"
                            ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,2)).EQ.1)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2))="DISCONTINUOUS"                                 
                            ELSE
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2))="EDGEDISCONTINUOUS"
                            ENDIF
                        ENDIF                    
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,3)).EQ.0.D0)THEN
                            IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,3)).GE.3)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3))="DISCONTINUOUS"
                            ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,3)).EQ.1)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3))="DISCONTINUOUS"                                 
                            ELSE
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3))="EDGEDISCONTINUOUS"
                            ENDIF
                        ENDIF                    
                    
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,1)).EQ.0.D0.AND.SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,2)).EQ.0.D0)THEN
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)).EQ."EDGEDISCONTINUOUS")THEN
                                EDGE_DISCONTINUITY(N_ELEM+J,1)=3
                            ENDIF
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)).EQ."EDGEDISCONTINUOUS")THEN    
                                EDGE_DISCONTINUITY(N_ELEM+J,2)=2
                            ENDIF                   
                        ENDIF
                        
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,2)).EQ.0.D0.AND.SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,3)).EQ.0.D0)THEN
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)).EQ."EDGEDISCONTINUOUS")THEN                        
                                EDGE_DISCONTINUITY(N_ELEM+J,2)=1
                            ENDIF
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)).EQ."EDGEDISCONTINUOUS")THEN                                
                                EDGE_DISCONTINUITY(N_ELEM+J,3)=3
                            ENDIF                       
                        ENDIF                        
                    
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,3)).EQ.0.D0.AND.SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,1)).EQ.0.D0)THEN
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)).EQ."EDGEDISCONTINUOUS")THEN                        
                                EDGE_DISCONTINUITY(N_ELEM+J,3)=2
                            ENDIF
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)).EQ."EDGEDISCONTINUOUS")THEN                                
                                EDGE_DISCONTINUITY(N_ELEM+J,1)=1
                            ENDIF                       
                        ENDIF                                               
!               !                                                            
                    ELSE
                    !--------QUADRILATERAL ELEMENTS-------------
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,1)).EQ.0.D0)THEN
                            IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,1)).GE.3)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1))="DISCONTINUOUS"
                            ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,1)).EQ.1)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1))="DISCONTINUOUS"                                 
                            ELSE
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1))="EDGEDISCONTINUOUS"
                            ENDIF
                        ENDIF    
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,2)).EQ.0.D0)THEN
                            IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,2)).GE.3)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2))="DISCONTINUOUS"
                            ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,2)).EQ.1)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2))="DISCONTINUOUS"                                 
                            ELSE
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2))="EDGEDISCONTINUOUS"
                            ENDIF
                        ENDIF                    
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,3)).EQ.0.D0)THEN
                            IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,3)).GE.3)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3))="DISCONTINUOUS"
                            ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,3)).EQ.1)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3))="DISCONTINUOUS"                                 
                            ELSE
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3))="EDGEDISCONTINUOUS"
                            ENDIF
                        ENDIF 
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,4)).EQ.0.D0)THEN
                            IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,4)).GE.3)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,4))="DISCONTINUOUS"
                            ELSE IF(COMMON_ELEMENTS(PARAMETRIC_CONNECTIVITY(J,4)).EQ.1)THEN
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,4))="DISCONTINUOUS"                                 
                            ELSE
                                DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,4))="EDGEDISCONTINUOUS"
                            ENDIF
                        ENDIF                                            
                    
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,1)).EQ.0.D0.AND.SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,2)).EQ.0.D0)THEN
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)).EQ."EDGEDISCONTINUOUS")THEN
                                EDGE_DISCONTINUITY(N_ELEM+J,1)=4
                            ENDIF
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)).EQ."EDGEDISCONTINUOUS")THEN    
                                EDGE_DISCONTINUITY(N_ELEM+J,2)=2
                            ENDIF                   
                        ENDIF
                    
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,2)).EQ.0.D0.AND.SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,3)).EQ.0.D0)THEN
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2)).EQ."EDGEDISCONTINUOUS")THEN                        
                                EDGE_DISCONTINUITY(N_ELEM+J,2)=1
                            ENDIF
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)).EQ."EDGEDISCONTINUOUS")THEN                                
                                EDGE_DISCONTINUITY(N_ELEM+J,3)=3
                            ENDIF                       
                        ENDIF    
                    
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,3)).EQ.0.D0.AND.SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,4)).EQ.0.D0)THEN
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3)).EQ."EDGEDISCONTINUOUS")THEN                        
                                EDGE_DISCONTINUITY(N_ELEM+J,3)=2
                            ENDIF
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,4)).EQ."EDGEDISCONTINUOUS")THEN                                
                                EDGE_DISCONTINUITY(N_ELEM+J,4)=4
                            ENDIF                       
                        ENDIF           
                        
                        IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,4)).EQ.0.D0.AND.SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,1)).EQ.0.D0)THEN
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,4)).EQ."EDGEDISCONTINUOUS")THEN                        
                                EDGE_DISCONTINUITY(N_ELEM+J,4)=3
                            ENDIF
                            IF(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1)).EQ."EDGEDISCONTINUOUS")THEN                                
                                EDGE_DISCONTINUITY(N_ELEM+J,1)=1
                            ENDIF                       
                        ENDIF                                                              
                    ENDIF
                ENDDO    
                                                       
            CASE("CD")
                 DO J=1,N_PARAMETRICELEMENTS
                    IF(PARAMETRIC_CONNECTIVITY(J,4).EQ.0.OR.PARAMETRIC_CONNECTIVITY(J,3).EQ.PARAMETRIC_CONNECTIVITY(J,4))THEN
                    !--------TRIANGULAR ELEMENTS-------------  
                         IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,1)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1))=1
                         ENDIF
                         IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,2)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2))=1
                         ENDIF
                          IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,3)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3))=1
                         ENDIF                                                 
                    ELSE
                    !--------QUADRILATERAL ELEMENTS-------------                    
                         IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,1)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1))=1
                         ENDIF
                         IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,2)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2))=1
                         ENDIF
                         IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,3)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3))=1
                         ENDIF 
                         IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,4)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,4))=1
                         ENDIF                                              
                    ENDIF
                ENDDO           
            CASE("CT")
                DO J=1,N_PARAMETRICELEMENTS
                    IF(PARAMETRIC_CONNECTIVITY(J,4).EQ.0.OR.PARAMETRIC_CONNECTIVITY(J,3).EQ.PARAMETRIC_CONNECTIVITY(J,4))THEN
                    !--------TRIANGULAR ELEMENTS-------------                      
                         IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,1)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1))=1
                         ENDIF
                         IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,2)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2))=1
                         ENDIF
                          IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,3)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3))=1
                         ENDIF                                                 
                    ELSE
                    !--------QUADRILATERAL ELEMENTS-------------                    
                         IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,1)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,1))=1
                         ENDIF
                         IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,2)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,2))=1
                         ENDIF
                         IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,3)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,3))=1
                         ENDIF 
                         IF(SIGNED_DIST_PARAMETRICNODES(PARAMETRIC_CONNECTIVITY(J,4)).EQ.0.D0)THEN
                             CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+J,4))=1
                         ENDIF                                              
                    ENDIF
                ENDDO                      
            ENDSELECT            
!------------------------------------------------------------------------------------------------------------------------------------------------------
!           UPDATING THE MESH SIZE
!------------------------------------------------------------------------------------------------------------------------------------------------------                            
            N_NODES=N_NODES+N_PARAMETRICNODES
            N_COLLOCPOINTS=N_COLLOCPOINTS+N_PARAMETRICCOLLOCPOINTS
            N_ELEM=N_ELEM+N_PARAMETRICELEMENTS 
            DEALLOCATE(PARAMETRICCOORD_NODES,PARAMETRIC_CONNECTIVITY,SIGNED_DIST_PARAMETRICNODES,COMMON_ELEMENTS,REPEATED_CP)                      
        ENDIF
    ENDDO
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   INVERTING SURFACES
!------------------------------------------------------------------------------------------------------------------------------------------------------ 
    ALLOCATE(EDGE_DISCONTINUITY_AUX(1,4),NODES_CONNECTIVITY_AUX(1,4),COLLOCPOINTS_CONNECTIVITY_AUX(1,4),CRACK_TIP_AUX(4),DISCONTINUOUS_COLLOCPOINTS_AUX(4))
    CONT1=0
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).EQ."UNTRIMMED")THEN    
            IF(SURFACE_INVERTION(I).EQ.-1)THEN
                DO J=CONT1+1,CONT1+N_ELEM_SURFACES(I)
                    NODES_CONNECTIVITY_AUX(1,:)=NODES_CONNECTIVITY(J,:)  
                    COLLOCPOINTS_CONNECTIVITY_AUX(1,:)=COLLOCPOINTS_CONNECTIVITY(J,:)
                    EDGE_DISCONTINUITY_AUX(1,:)=EDGE_DISCONTINUITY(J,:) 
                    SELECT CASE(ELEM_TYPE(J))
                    CASE(3)
                        CRACK_TIP_AUX(1)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,1))
                        CRACK_TIP_AUX(2)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,2))
                        CRACK_TIP_AUX(3)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,3))
    !
                        DISCONTINUOUS_COLLOCPOINTS_AUX(1)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,1)) 
                        DISCONTINUOUS_COLLOCPOINTS_AUX(2)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,2))
                        DISCONTINUOUS_COLLOCPOINTS_AUX(3)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,3))
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,1).EQ.1)THEN
                            EDGE_DISCONTINUITY(J,1)=3
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,1).EQ.3)THEN
                             EDGE_DISCONTINUITY(J,1)=1                   
                        ENDIF
    !                    
                        NODES_CONNECTIVITY(J,2)=NODES_CONNECTIVITY_AUX(1,3) 
                        COLLOCPOINTS_CONNECTIVITY(J,2)=COLLOCPOINTS_CONNECTIVITY_AUX(1,3)
                        CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,2))=CRACK_TIP_AUX(3)    
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,2))=DISCONTINUOUS_COLLOCPOINTS_AUX(3)
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,3).EQ.2)THEN
                            EDGE_DISCONTINUITY(J,2)=2
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,3).EQ.3)THEN
                             EDGE_DISCONTINUITY(J,2)=1                   
                        ENDIF                    
    !       
                        NODES_CONNECTIVITY(J,3)=NODES_CONNECTIVITY_AUX(1,2) 
                        COLLOCPOINTS_CONNECTIVITY(J,3)=COLLOCPOINTS_CONNECTIVITY_AUX(1,2)
                        CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,3))=CRACK_TIP_AUX(2)    
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,3))=DISCONTINUOUS_COLLOCPOINTS_AUX(2) 
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,2).EQ.1)THEN
                            EDGE_DISCONTINUITY(J,3)=3
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,2).EQ.2)THEN
                             EDGE_DISCONTINUITY(J,3)=2                   
                        ENDIF                    
    !                                                                          
                    CASE(4)
                        CRACK_TIP_AUX(1)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,1))
                        CRACK_TIP_AUX(2)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,2))
                        CRACK_TIP_AUX(3)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,3))
                        CRACK_TIP_AUX(4)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,4))                    
    !
                        DISCONTINUOUS_COLLOCPOINTS_AUX(1)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,1)) 
                        DISCONTINUOUS_COLLOCPOINTS_AUX(2)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,2))
                        DISCONTINUOUS_COLLOCPOINTS_AUX(3)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,3))
                        DISCONTINUOUS_COLLOCPOINTS_AUX(4)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,4))                    
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,1).EQ.1)THEN
                            EDGE_DISCONTINUITY(J,1)=4
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,1).EQ.4)THEN
                             EDGE_DISCONTINUITY(J,1)=1                   
                        ENDIF
    !                    
                        NODES_CONNECTIVITY(J,2)=NODES_CONNECTIVITY_AUX(1,4) 
                        COLLOCPOINTS_CONNECTIVITY(J,2)=COLLOCPOINTS_CONNECTIVITY_AUX(1,4)
                        CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,2))=CRACK_TIP_AUX(4)    
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,2))=DISCONTINUOUS_COLLOCPOINTS_AUX(4)
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,4).EQ.3)THEN
                            EDGE_DISCONTINUITY(J,2)=2
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,4).EQ.4)THEN
                             EDGE_DISCONTINUITY(J,2)=1                   
                        ENDIF
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,3).EQ.2)THEN
                            EDGE_DISCONTINUITY(J,3)=3
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,3).EQ.3)THEN
                             EDGE_DISCONTINUITY(J,3)=2                   
                        ENDIF
    !                    
                        NODES_CONNECTIVITY(J,4)=NODES_CONNECTIVITY_AUX(1,2) 
                        COLLOCPOINTS_CONNECTIVITY(J,4)=COLLOCPOINTS_CONNECTIVITY_AUX(1,2)
                        CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,4))=CRACK_TIP_AUX(2)    
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,4))=DISCONTINUOUS_COLLOCPOINTS_AUX(2)
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,2).EQ.1)THEN
                            EDGE_DISCONTINUITY(J,4)=4
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,2).EQ.2)THEN
                             EDGE_DISCONTINUITY(J,4)=3                   
                        ENDIF
    !                                                                            
                    ENDSELECT                             
                ENDDO
            ENDIF
            CONT1=CONT1+N_ELEM_SURFACES(I)
        ENDIF
    ENDDO
    DO I=1,N_SURFACES
        IF(TRIMMED_SURFACE(I).NE."UNTRIMMED")THEN    
            IF(SURFACE_INVERTION(I).EQ.-1)THEN
                DO J=CONT1+1,CONT1+N_ELEM_SURFACES(I)
                    NODES_CONNECTIVITY_AUX(1,:)=NODES_CONNECTIVITY(J,:)  
                    COLLOCPOINTS_CONNECTIVITY_AUX(1,:)=COLLOCPOINTS_CONNECTIVITY(J,:)
                    EDGE_DISCONTINUITY_AUX(1,:)=EDGE_DISCONTINUITY(J,:) 
                    SELECT CASE(ELEM_TYPE(J))
                    CASE(3)
                        CRACK_TIP_AUX(1)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,1))
                        CRACK_TIP_AUX(2)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,2))
                        CRACK_TIP_AUX(3)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,3))
    !
                        DISCONTINUOUS_COLLOCPOINTS_AUX(1)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,1)) 
                        DISCONTINUOUS_COLLOCPOINTS_AUX(2)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,2))
                        DISCONTINUOUS_COLLOCPOINTS_AUX(3)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,3))
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,1).EQ.1)THEN
                            EDGE_DISCONTINUITY(J,1)=3
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,1).EQ.3)THEN
                             EDGE_DISCONTINUITY(J,1)=1                   
                        ENDIF
    !                    
                        NODES_CONNECTIVITY(J,2)=NODES_CONNECTIVITY_AUX(1,3) 
                        COLLOCPOINTS_CONNECTIVITY(J,2)=COLLOCPOINTS_CONNECTIVITY_AUX(1,3)
                        CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,2))=CRACK_TIP_AUX(3)    
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,2))=DISCONTINUOUS_COLLOCPOINTS_AUX(3)
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,3).EQ.2)THEN
                            EDGE_DISCONTINUITY(J,2)=2
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,3).EQ.3)THEN
                             EDGE_DISCONTINUITY(J,2)=1                   
                        ENDIF                    
    !       
                        NODES_CONNECTIVITY(J,3)=NODES_CONNECTIVITY_AUX(1,2) 
                        COLLOCPOINTS_CONNECTIVITY(J,3)=COLLOCPOINTS_CONNECTIVITY_AUX(1,2)
                        CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,3))=CRACK_TIP_AUX(2)    
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,3))=DISCONTINUOUS_COLLOCPOINTS_AUX(2) 
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,2).EQ.1)THEN
                            EDGE_DISCONTINUITY(J,3)=3
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,2).EQ.2)THEN
                             EDGE_DISCONTINUITY(J,3)=2                   
                        ENDIF                    
    !                                                                          
                    CASE(4)
                        CRACK_TIP_AUX(1)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,1))
                        CRACK_TIP_AUX(2)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,2))
                        CRACK_TIP_AUX(3)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,3))
                        CRACK_TIP_AUX(4)=CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,4))                    
    !
                        DISCONTINUOUS_COLLOCPOINTS_AUX(1)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,1)) 
                        DISCONTINUOUS_COLLOCPOINTS_AUX(2)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,2))
                        DISCONTINUOUS_COLLOCPOINTS_AUX(3)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,3))
                        DISCONTINUOUS_COLLOCPOINTS_AUX(4)=DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,4))                    
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,1).EQ.1)THEN
                            EDGE_DISCONTINUITY(J,1)=4
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,1).EQ.4)THEN
                             EDGE_DISCONTINUITY(J,1)=1                   
                        ENDIF
    !                    
                        NODES_CONNECTIVITY(J,2)=NODES_CONNECTIVITY_AUX(1,4) 
                        COLLOCPOINTS_CONNECTIVITY(J,2)=COLLOCPOINTS_CONNECTIVITY_AUX(1,4)
                        CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,2))=CRACK_TIP_AUX(4)    
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,2))=DISCONTINUOUS_COLLOCPOINTS_AUX(4)
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,4).EQ.3)THEN
                            EDGE_DISCONTINUITY(J,2)=2
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,4).EQ.4)THEN
                             EDGE_DISCONTINUITY(J,2)=1                   
                        ENDIF
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,3).EQ.2)THEN
                            EDGE_DISCONTINUITY(J,3)=3
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,3).EQ.3)THEN
                             EDGE_DISCONTINUITY(J,3)=2                   
                        ENDIF
    !                    
                        NODES_CONNECTIVITY(J,4)=NODES_CONNECTIVITY_AUX(1,2) 
                        COLLOCPOINTS_CONNECTIVITY(J,4)=COLLOCPOINTS_CONNECTIVITY_AUX(1,2)
                        CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,4))=CRACK_TIP_AUX(2)    
                        DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,4))=DISCONTINUOUS_COLLOCPOINTS_AUX(2)
    !                                                                                                                       
                        IF(EDGE_DISCONTINUITY_AUX(1,2).EQ.1)THEN
                            EDGE_DISCONTINUITY(J,4)=4
                        ELSE IF(EDGE_DISCONTINUITY_AUX(1,2).EQ.2)THEN
                             EDGE_DISCONTINUITY(J,4)=3                   
                        ENDIF
    !                                                                            
                    ENDSELECT                             
                ENDDO
            ENDIF
            CONT1=CONT1+N_ELEM_SURFACES(I)            
        ENDIF
    ENDDO    
    DEALLOCATE(EDGE_DISCONTINUITY_AUX,NODES_CONNECTIVITY_AUX,COLLOCPOINTS_CONNECTIVITY_AUX,CRACK_TIP_AUX,DISCONTINUOUS_COLLOCPOINTS_AUX)
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   COLLOCATION POINTS GENERATION: DIMENSIONLESS COORDENATES   
!------------------------------------------------------------------------------------------------------------------------------------------------------
    DO I=1,N_ELEM   
        SELECT CASE(ELEM_TYPE(I))
        CASE(3)
            SELECT CASE(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(I,1)))
            CASE("CONTINUOUS")
                QSI(I,1,1)=1.D0 
                QSI(I,1,2)=0.D0
            CASE("DISCONTINUOUS")                
                QSI(I,1,1)=0.8D0 
                QSI(I,1,2)=0.1D0
            CASE("EDGEDISCONTINUOUS")
                IF(EDGE_DISCONTINUITY(I,1).EQ.1)THEN                   
                    QSI(I,1,1)=0.8D0 
                    QSI(I,1,2)=0.2D0
                ELSE IF(EDGE_DISCONTINUITY(I,1).EQ.3)THEN
                    QSI(I,1,1)=0.8D0 
                    QSI(I,1,2)=0.D0                
                ENDIF                             
            END SELECT     
!
            SELECT CASE(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(I,2)))
            CASE("CONTINUOUS")
                QSI(I,2,1)=0.D0
                QSI(I,2,2)=1.D0
            CASE("DISCONTINUOUS")                
                QSI(I,2,1)=0.1D0
                QSI(I,2,2)=0.8D0
            CASE("EDGEDISCONTINUOUS")
                IF(EDGE_DISCONTINUITY(I,2).EQ.1)THEN                   
                    QSI(I,2,1)=0.2D0 
                    QSI(I,2,2)=0.8D0
                ELSE IF(EDGE_DISCONTINUITY(I,2).EQ.2)THEN
                    QSI(I,2,1)=0.D0 
                    QSI(I,2,2)=0.8D0                
                ENDIF                             
            END SELECT         
!
            SELECT CASE(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(I,3)))
            CASE("CONTINUOUS")
                QSI(I,3,1)=0.D0
                QSI(I,3,2)=0.D0 
            CASE("DISCONTINUOUS")                
                QSI(I,3,1)=0.1D0
                QSI(I,3,2)=0.1D0 
            CASE("EDGEDISCONTINUOUS")
                IF(EDGE_DISCONTINUITY(I,3).EQ.2)THEN                   
                    QSI(I,3,1)=0.D0 
                    QSI(I,3,2)=0.2D0
                ELSE IF(EDGE_DISCONTINUITY(I,3).EQ.3)THEN
                    QSI(I,3,1)=0.2D0
                    QSI(I,3,2)=0.D0                
                ENDIF                             
            END SELECT         
!
        CASE(4)       
            SELECT CASE(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(I,1)))
            CASE("CONTINUOUS")
                QSI(I,1,1)=-1.D0
                QSI(I,1,2)=-1.D0
            CASE("DISCONTINUOUS")                
                QSI(I,1,1)=-0.8D0
                QSI(I,1,2)=-0.8D0  
            CASE("EDGEDISCONTINUOUS")
                IF(EDGE_DISCONTINUITY(I,1).EQ.1)THEN                   
                    QSI(I,1,1)=-0.8D0
                    QSI(I,1,2)=-1.D0
                ELSE IF(EDGE_DISCONTINUITY(I,1).EQ.4)THEN
                    QSI(I,1,1)=-1.D0
                    QSI(I,1,2)=-0.8D0                
                ENDIF                             
            END SELECT
!
            SELECT CASE(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(I,2)))
            CASE("CONTINUOUS")
                QSI(I,2,1)=1.D0
                QSI(I,2,2)=-1.D0
            CASE("DISCONTINUOUS")                
                QSI(I,2,1)=0.8D0
                QSI(I,2,2)=-0.8D0  
            CASE("EDGEDISCONTINUOUS")
                IF(EDGE_DISCONTINUITY(I,2).EQ.1)THEN                   
                    QSI(I,2,1)=0.8D0
                    QSI(I,2,2)=-1.D0
                ELSE IF(EDGE_DISCONTINUITY(I,2).EQ.2)THEN
                    QSI(I,2,1)=1.D0
                    QSI(I,2,2)=-0.8D0                
                ENDIF                             
            END SELECT
!
            SELECT CASE(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(I,3)))
            CASE("CONTINUOUS")
                QSI(I,3,1)=1.D0
                QSI(I,3,2)=1.D0
            CASE("DISCONTINUOUS")                
                QSI(I,3,1)=0.8D0
                QSI(I,3,2)=0.8D0  
            CASE("EDGEDISCONTINUOUS")
                IF(EDGE_DISCONTINUITY(I,3).EQ.2)THEN                   
                    QSI(I,3,1)=1.D0
                    QSI(I,3,2)=0.8D0
                ELSE IF(EDGE_DISCONTINUITY(I,3).EQ.3)THEN
                    QSI(I,3,1)=0.8D0
                    QSI(I,3,2)=1.D0               
                ENDIF                             
            END SELECT
!
            SELECT CASE(DISCONTINUOUS_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(I,4)))
            CASE("CONTINUOUS")
                QSI(I,4,1)=-1.D0
                QSI(I,4,2)=1.D0
            CASE("DISCONTINUOUS")                
                QSI(I,4,1)=-0.8D0
                QSI(I,4,2)=0.8D0  
            CASE("EDGEDISCONTINUOUS")
                IF(EDGE_DISCONTINUITY(I,4).EQ.3)THEN                   
                    QSI(I,4,1)=-0.8D0
                    QSI(I,4,2)=1.D0 
                ELSE IF(EDGE_DISCONTINUITY(I,4).EQ.4)THEN
                    QSI(I,4,1)=-1.D0
                    QSI(I,4,2)=0.8D0            
                ENDIF                             
            END SELECT            
        ENDSELECT       
    ENDDO  
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   COLLOCATION POINTS GENERATION: PHYSICAL COORDENATES   
!------------------------------------------------------------------------------------------------------------------------------------------------------
    COORD_COLLOCPOINTS=0.D0
    DO I=1,N_COLLOCPOINTS
        DO J=1,N_ELEM
            NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
            ALLOCATE(VALUES(NEN,3))
            DO K=1,NEN
                VALUES(K,1)=COORD_NODES(NODES_CONNECTIVITY(J,K),1)
	            VALUES(K,2)=COORD_NODES(NODES_CONNECTIVITY(J,K),2)
	            VALUES(K,3)=COORD_NODES(NODES_CONNECTIVITY(J,K),3)
            ENDDO
            DO K=1,NEN
                IF (COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
                    CALL SHAPE_FUNCTIONS(1,QSI(J,K,1),QSI(J,K,2),NEN,VALUES,COEFFICIENTS,X1,X2,X3,PHI,DPHI,D2PHI)
                    COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),1)=X1
                    COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),2)=X2
                    COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),3)=X3
                ENDIF
            ENDDO
            DEALLOCATE(VALUES)
        ENDDO
    ENDDO    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   OTHER VARIABLES   
!------------------------------------------------------------------------------------------------------------------------------------------------------    
    B_CONDITIONS=1
    U=0.D0
    T=0.D0
    H=0.D0
    G=0.D0    
!  
    10  FORMAT(a)
    20  FORMAT(12X,i)
    30  FORMAT(28x,E21.18,E20.18)          
!          
	END SUBROUTINE