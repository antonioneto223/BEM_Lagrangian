	SUBROUTINE NORMAL_OUTWARD_OUTPUT
!
	USE ISOPARAMETRIC_MESH
	USE CRACKS_DUAL_BEM
    USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
	USE REMESHING
!
	IMPLICIT NONE 
!
	INTEGER::I,J,K,II,JJ,KK,L,NEN,C,N_CRACKED_ELEM,CONT
!	
	REAL*8::QSI1,QSI2,ETA(3),ETA_NODES(3*N_COLLOCPOINTS),VALUES[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:),PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),VJC(3),DVJC(3,2),JC,X,Y,Z,DX,DY,DZ,R	
!    
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   EXTRAPOLATING DISPLACEMENTS AND TRACTIONS
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
    ETA_NODES=0.D0
    DO I=1,N_COLLOCPOINTS
    	C=0
	    DO J=1,N_ELEM
	        IF(CRACKED_ELEM(J).EQ."UNCRACKED")THEN
	            NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
	            DO K=1,NEN
	                IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
	                    C=C+1
	                    IF(C.LE.1)THEN
	                        ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN))
	                        DO L=1,NEN
		                        VALUES(L,1)=COORD_NODES(NODES_CONNECTIVITY(J,L),1)
		                        VALUES(L,2)=COORD_NODES(NODES_CONNECTIVITY(J,L),2)
		                        VALUES(L,3)=COORD_NODES(NODES_CONNECTIVITY(J,L),3)                  
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
	                        CALL NORMAL_OUTWARD(NEN,VALUES,QSI1,QSI2,ETA,VJC,DVJC,JC)

                            ETA_NODES(3*I-2)=ETA(1)   
                            ETA_NODES(3*I-1)=ETA(2)   
                            ETA_NODES(3*I)=ETA(3)  	 
	                                              
	                        DEALLOCATE(VALUES,COEFFICIENTS,PHI)
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
	OPEN(3,file='Output_data\NORMAL_OUTWARD_OUTPUT.ogl',status='unknown')
!
	WRITE(3,*)'Arquivo de pós-processamento'
	WRITE(3,*)
	WRITE(3,*)'Nº de nos     Nº de elementos     Nº de listas'
	WRITE(3,10)
	WRITE(3,100)N_COLLOCPOINTS,N_ELEM-N_CRACKED_ELEM,3
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
	                        WRITE(3,200)COORD_NODES(NODES_CONNECTIVITY(J,K),1),COORD_NODES(NODES_CONNECTIVITY(J,K),2),COORD_NODES(NODES_CONNECTIVITY(J,K),3),0.D0,0.D0,0.D0
	                    ENDIF
	                ENDIF
	            ENDDO
	        ENDIF    
	    ENDDO
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'tpelem (1 - barra / 2 - triang / 3 - quad) grauaprox nó1 nó2...non'
	WRITE(3,10)
	DO I=1,N_ELEM
	    IF(CRACKED_ELEM(I).EQ."UNCRACKED")THEN	
	        NEN=ELEM_TYPE(I)*ORDER_ELEM(I)+(ELEM_TYPE(I)-3)*(ORDER_ELEM(I)-1)*POL_FAMILY(I)
            SELECT CASE(NEN)
            CASE(3)
                WRITE(3,300)2,1,COLLOCPOINTS_CONNECTIVITY(I,3),COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,2),ND(I)
            CASE(6)
                WRITE(3,301)2,2,COLLOCPOINTS_CONNECTIVITY(I,3),COLLOCPOINTS_CONNECTIVITY(I,6),COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,5),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,2),ND(I)
            CASE(4)
                WRITE(3,302)3,1,COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,2),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,3),ND(I)
            CASE(8)
                WRITE(3,303)3,2,COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,5),COLLOCPOINTS_CONNECTIVITY(I,2),COLLOCPOINTS_CONNECTIVITY(I,8),COLLOCPOINTS_CONNECTIVITY(I,6),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,7),COLLOCPOINTS_CONNECTIVITY(I,3),ND(I)
            CASE(9)
                WRITE(3,304)3,2,COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,5),COLLOCPOINTS_CONNECTIVITY(I,2),COLLOCPOINTS_CONNECTIVITY(I,8),COLLOCPOINTS_CONNECTIVITY(I,9),COLLOCPOINTS_CONNECTIVITY(I,6),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,7),COLLOCPOINTS_CONNECTIVITY(I,3),ND(I)      
            ENDSELECT
        ENDIF
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'Componente ETAx do versor normal no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'ETAx'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)0.D0,0.D0,0.D0,ETA_NODES(3*I-2)    
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'Componente ETAy do versor normal no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'ETAy'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)0.D0,0.D0,0.D0,ETA_NODES(3*I-1)    
	ENDDO
	WRITE(3,*)
	WRITE(3,*)'Componente ETAz do versor normal no sistema de referência global'
	WRITE(3,10)
	WRITE(3,*)'ETAz'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)0.D0,0.D0,0.D0,ETA_NODES(3*I)    
	ENDDO
!	
	CLOSE(3)
!
10  FORMAT(1h#)
100	FORMAT(2x,i5,12X,i5,12X,i5)
200	FORMAT(8x,F14.9,16x,F14.9,16x,F14.9,16x,F14.9,16X,F14.9,16x,F14.9)
300	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7X,i5,7x,i2)
301	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i2)
302	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i2)
303	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i2)
304	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i2)
400	FORMAT(8x,E13.6,4x,E13.6,4x,E13.6,15x,E13.6)
500 FORMAT(a,i4,a)
!	
	END SUBROUTINE











