	SUBROUTINE INTERCEPTED_EDGES(ELEM,EDGE,INTERCEPTED_EDGE)
!
	USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
	USE REMESHING
!
	IMPLICIT NONE 
!
	INTEGER::I,J,K,L,JJ,KK,ELEM,EDGE,NEN,IPIV(2),INFO
!
    REAL*8::D1,D2,PI1(3),PI2(3),PF1(3),PF2(3),V1(3),V2(3),A1,A2,X1(3),X2(3),R,BOUNDARYCRACK_TOL,A,B,C,D,E,F		
    REAL*8::VALUES[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:),PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),X,Y,Z,QSI1,QSI2    	
!    
    LOGICAL::INTERCEPTED_EDGE	
!
    BOUNDARYCRACK_TOL=1.D-8
    INTERCEPTED_EDGE=.FALSE.
!    
    SELECT CASE(ELEM_TYPE(ELEM))
    CASE(3)  
        SELECT CASE(EDGE)
        CASE(1)
            PI1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,1),1)
            PI1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,1),2)
            PI1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,1),3)
! 
            PF1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,2),1)
            PF1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,2),2)
            PF1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,2),3)                       
        CASE(2)
            PI1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,2),1)
            PI1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,2),2)
            PI1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,2),3)
! 
            PF1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,3),1)
            PF1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,3),2)
            PF1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,3),3)    
        CASE(3)
            PI1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,3),1)
            PI1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,3),2)
            PI1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,3),3)
! 
            PF1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,1),1)
            PF1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,1),2)
            PF1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,1),3)       
        ENDSELECT    
    CASE(4)  
        SELECT CASE(EDGE)
        CASE(1)
            PI1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,1),1)
            PI1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,1),2)
            PI1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,1),3)
! 
            PF1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,2),1)
            PF1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,2),2)
            PF1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,2),3)                       
        CASE(2)
            PI1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,2),1)
            PI1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,2),2)
            PI1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,2),3)
! 
            PF1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,3),1)
            PF1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,3),2)
            PF1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,3),3)    
        CASE(3)
            PI1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,3),1)
            PI1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,3),2)
            PI1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,3),3)
! 
            PF1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,4),1)
            PF1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,4),2)
            PF1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,4),3)    
        CASE(4)
            PI1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,4),1)
            PI1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,4),2)
            PI1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,4),3)
! 
            PF1(1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,1),1)
            PF1(2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,1),2)
            PF1(3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,1),3)   
        ENDSELECT
    ENDSELECT     
    D1=DSQRT((PF1(1)-PI1(1))**2+(PF1(2)-PI1(2))**2+(PF1(3)-PI1(3))**2)
    V1=(PF1-PI1)/D1
    
    DO J=1,N_ELEM
        IF(DUAL_BEM(COLLOCPOINTS_CONNECTIVITY(J,1)).EQ."S")THEN
            DO K=1,ELEM_TYPE(J)
                SELECT CASE(ELEM_TYPE(J))
                CASE(3)
                    SELECT CASE(K)
                    CASE(1)
                        PI2(1)=COORD_NODES(NODES_CONNECTIVITY(J,1),1)
                        PI2(2)=COORD_NODES(NODES_CONNECTIVITY(J,1),2)
                        PI2(3)=COORD_NODES(NODES_CONNECTIVITY(J,1),3)
! 
                        PF2(1)=COORD_NODES(NODES_CONNECTIVITY(J,2),1)
                        PF2(2)=COORD_NODES(NODES_CONNECTIVITY(J,2),2)
                        PF2(3)=COORD_NODES(NODES_CONNECTIVITY(J,2),3)                 
                    CASE(2)
                        PI2(1)=COORD_NODES(NODES_CONNECTIVITY(J,2),1)
                        PI2(2)=COORD_NODES(NODES_CONNECTIVITY(J,2),2)
                        PI2(3)=COORD_NODES(NODES_CONNECTIVITY(J,2),3)
!  
                        PF2(1)=COORD_NODES(NODES_CONNECTIVITY(J,3),1)
                        PF2(2)=COORD_NODES(NODES_CONNECTIVITY(J,3),2)
                        PF2(3)=COORD_NODES(NODES_CONNECTIVITY(J,3),3)                 
                    CASE(3)
                        PI2(1)=COORD_NODES(NODES_CONNECTIVITY(J,3),1)
                        PI2(2)=COORD_NODES(NODES_CONNECTIVITY(J,3),2)
                        PI2(3)=COORD_NODES(NODES_CONNECTIVITY(J,3),3)
! 
                        PF2(1)=COORD_NODES(NODES_CONNECTIVITY(J,1),1)
                        PF2(2)=COORD_NODES(NODES_CONNECTIVITY(J,1),2)
                        PF2(3)=COORD_NODES(NODES_CONNECTIVITY(J,1),3)                                 
                    ENDSELECT                
                CASE(4)
                    SELECT CASE(K)
                    CASE(1)
                        PI2(1)=COORD_NODES(NODES_CONNECTIVITY(J,1),1)
                        PI2(2)=COORD_NODES(NODES_CONNECTIVITY(J,1),2)
                        PI2(3)=COORD_NODES(NODES_CONNECTIVITY(J,1),3)
! 
                        PF2(1)=COORD_NODES(NODES_CONNECTIVITY(J,2),1)
                        PF2(2)=COORD_NODES(NODES_CONNECTIVITY(J,2),2)
                        PF2(3)=COORD_NODES(NODES_CONNECTIVITY(J,2),3)                 
                    CASE(2)
                       PI2(1)=COORD_NODES(NODES_CONNECTIVITY(J,2),1)
                       PI2(2)=COORD_NODES(NODES_CONNECTIVITY(J,2),2)
                       PI2(3)=COORD_NODES(NODES_CONNECTIVITY(J,2),3)
! 
                        PF2(1)=COORD_NODES(NODES_CONNECTIVITY(J,3),1)
                        PF2(2)=COORD_NODES(NODES_CONNECTIVITY(J,3),2)
                        PF2(3)=COORD_NODES(NODES_CONNECTIVITY(J,3),3)                 
                    CASE(3)
                        PI2(1)=COORD_NODES(NODES_CONNECTIVITY(J,3),1)
                        PI2(2)=COORD_NODES(NODES_CONNECTIVITY(J,3),2)
                        PI2(3)=COORD_NODES(NODES_CONNECTIVITY(J,3),3)
! 
                        PF2(1)=COORD_NODES(NODES_CONNECTIVITY(J,4),1)
                        PF2(2)=COORD_NODES(NODES_CONNECTIVITY(J,4),2)
                        PF2(3)=COORD_NODES(NODES_CONNECTIVITY(J,4),3)                 
                    CASE(4)
                        PI2(1)=COORD_NODES(NODES_CONNECTIVITY(J,4),1)
                        PI2(2)=COORD_NODES(NODES_CONNECTIVITY(J,4),2)
                        PI2(3)=COORD_NODES(NODES_CONNECTIVITY(J,4),3)
! 
                        PF2(1)=COORD_NODES(NODES_CONNECTIVITY(J,1),1)
                        PF2(2)=COORD_NODES(NODES_CONNECTIVITY(J,1),2)
                        PF2(3)=COORD_NODES(NODES_CONNECTIVITY(J,1),3)                 
                    ENDSELECT
                ENDSELECT
                D2=DSQRT((PF2(1)-PI2(1))**2+(PF2(2)-PI2(2))**2+(PF2(3)-PI2(3))**2)
                V2=(PF2-PI2)/D2
!                
                A=(V1(1)**2+V1(2)**2+V1(3)**2)
                B=-(V2(1)*V1(1)+V2(2)*V1(2)+V2(3)*V1(3))
                C=(PI2(1)*V1(1)+PI2(2)*V1(2)+PI2(3)*V1(3))-(PI1(1)*V1(1)+PI1(2)*V1(2)+PI1(3)*V1(3))
                D=(V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3))
                E=-(V2(1)**2+V2(2)**2+V2(3)**2)
                F=(PI2(1)*V2(1)+PI2(2)*V2(2)+PI2(3)*V2(3))-(PI1(1)*V2(1)+PI1(2)*V2(2)+PI1(3)*V2(3))
!
                A2=(F*A-D*C)/(E*A-D*B)
                A1=(C-B*A2)/A
!
                IF(0.D0.LE.A1.AND.A1.LE.D1)THEN
                    IF(0.D0.LE.A2.AND.A2.LE.D2)THEN
                        X2=PI2+V2*A2   
                        X1=PI1+V1*A1
                        R=DSQRT((X2(1)-X1(1))**2+(X2(2)-X1(2))**2+(X2(3)-X1(3))**2)
                        IF(R/(D1+D2).LE.BOUNDARYCRACK_TOL)THEN 
                            COORD_INTERCEPTPOINTS(ELEM,EDGE,:)=X2
                            INTERCEPTED_EDGE=.TRUE.
!                           
                            NEN=ELEM_TYPE(ELEM)*ORDER_ELEM(ELEM)+(ELEM_TYPE(ELEM)-3)*(ORDER_ELEM(ELEM)-1)*POL_FAMILY(ELEM)
                            ALLOCATE(VALUES(NEN,3))
                            DO L=1,NEN
                                VALUES(L,1)=COORD_NODES(NODES_CONNECTIVITY(ELEM,L),1)
	                            VALUES(L,2)=COORD_NODES(NODES_CONNECTIVITY(ELEM,L),2)
	                            VALUES(L,3)=COORD_NODES(NODES_CONNECTIVITY(ELEM,L),3)
	                        ENDDO               
!
                            IF(EDGE.EQ.1)THEN
                                QSI1=(2.D0/D1)*A1-1.D0
                                CALL SHAPE_FUNCTIONS(1,QSI1,1.D0,4,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI) 
                                COORD_EXTRAINTERCEPTPOINTS(ELEM,1)=X
                                COORD_EXTRAINTERCEPTPOINTS(ELEM,2)=Y 
                                COORD_EXTRAINTERCEPTPOINTS(ELEM,3)=Z                               
                            ENDIF                           
                            IF(EDGE.EQ.2)THEN
                                QSI2=(2.D0/D1)*A1-1.D0
                                CALL SHAPE_FUNCTIONS(1,-1.D0,QSI2,4,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI) 
                                COORD_EXTRAINTERCEPTPOINTS(ELEM,1)=X
                                COORD_EXTRAINTERCEPTPOINTS(ELEM,2)=Y 
                                COORD_EXTRAINTERCEPTPOINTS(ELEM,3)=Z                               
                            ENDIF
                            IF(EDGE.EQ.3)THEN
                                QSI1=(-2.D0/D1)*A1+1.D0
                                CALL SHAPE_FUNCTIONS(1,QSI1,-1.D0,4,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI) 
                                COORD_EXTRAINTERCEPTPOINTS(ELEM,1)=X
                                COORD_EXTRAINTERCEPTPOINTS(ELEM,2)=Y 
                                COORD_EXTRAINTERCEPTPOINTS(ELEM,3)=Z                               
                            ENDIF    
                            IF(EDGE.EQ.4)THEN
                                QSI2=(-2.D0/D1)*A1+1.D0
                                CALL SHAPE_FUNCTIONS(1,1.D0,QSI2,4,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI) 
                                COORD_EXTRAINTERCEPTPOINTS(ELEM,1)=X
                                COORD_EXTRAINTERCEPTPOINTS(ELEM,2)=Y 
                                COORD_EXTRAINTERCEPTPOINTS(ELEM,3)=Z                               
                            ENDIF                                                                                
!                            
                            DEALLOCATE(VALUES)                                                   
                        ENDIF
                    ENDIF    
                ENDIF                                                            
!                                              
            ENDDO    
        ENDIF
    ENDDO
! 	   
    END SUBROUTINE INTERCEPTED_EDGES