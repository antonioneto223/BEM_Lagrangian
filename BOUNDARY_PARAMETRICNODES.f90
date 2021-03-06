	SUBROUTINE BOUNDARY_PARAMETRICNODES(SURFACE,N_PARAMETRICNODES,N_PARAMETRICELEMENTS,PARAMETRIC_CONNECTIVITY,PARAMETRICCOORD_NODES,SIGNED_DIST_PARAMETRICNODES,N_CURVEPOINTS,TRIMMING_CURVEPOINTS,COMMON_ELEMENTS)  
!
	USE NURBS_SURFACES
	USE TRIMMING_CURVES
!
	IMPLICIT NONE 

    INTEGER::I,J,K,CONT1,CONT2,SURFACE,N_PARAMETRICNODES,N_PARAMETRICELEMENTS,PARAMETRIC_CONNECTIVITY(N_PARAMETRICELEMENTS,4),N_CURVEPOINTS,COMMON_ELEMENTS(N_PARAMETRICNODES),N_BOUNDARYPOINTS,N_BOUNDARYLINES,BOUNDARY_LINES[ALLOCATABLE](:,:)
    REAL*8::PARAMETRICCOORD_NODES(N_PARAMETRICNODES,2),SIGNED_DIST_PARAMETRICNODES(N_PARAMETRICNODES),TRIMMING_CURVEPOINTS(N_CURVEPOINTS,2),A,B,SIZE,D,TOL,BOUNDARY_POINTS[ALLOCATABLE](:,:)  	
!
    TOL=1.D-3
    SIGNED_DIST_PARAMETRICNODES=1.D0
!
    N_BOUNDARYLINES=0
    N_BOUNDARYPOINTS=0
    DO I=1,N_TRIMMINGCURVES
        IF(SURFACES(I).EQ.SURFACE)THEN
            DO J=1,SUBDIVISION(I)
                N_BOUNDARYLINES=N_BOUNDARYLINES+1
                N_BOUNDARYPOINTS=N_BOUNDARYPOINTS+1             
            ENDDO
        ENDIF
    ENDDO      
!    
    ALLOCATE(BOUNDARY_POINTS(N_BOUNDARYPOINTS,2),BOUNDARY_LINES(N_BOUNDARYLINES,2))  
!    
    CONT1=0
    CONT2=0
    DO I=1,N_TRIMMINGCURVES
        IF(SURFACES(I).EQ.SURFACE)THEN
            DO J=1,SUBDIVISION(I)
                CONT1=CONT1+1
                CONT2=CONT2+1
                BOUNDARY_POINTS(CONT2,1)=TRIMMING_CURVEPOINTS(CONT1,1)
                BOUNDARY_POINTS(CONT2,2)=TRIMMING_CURVEPOINTS(CONT1,2)                
            ENDDO
            CONT1=CONT1+1                   
        ENDIF
    ENDDO     
!   
    CONT1=0
    N_BOUNDARYLINES=0
    DO I=1,N_TRIMMINGCURVES
        IF(SURFACES(I).EQ.SURFACE)THEN
            DO J=1,SUBDIVISION(I)
                N_BOUNDARYLINES=N_BOUNDARYLINES+1
                CONT1=CONT1+1
                BOUNDARY_LINES(N_BOUNDARYLINES,1)=CONT1              
            ENDDO
            CONT1=CONT1+1
        ENDIF
    ENDDO
!    
    CONT1=0
    CONT2=0 
    DO I=1,N_TRIMMINGCURVES
        IF(SURFACES(I).EQ.SURFACE)THEN
            IF(OPENED_CLOSED_CURVE(I).EQ."O".AND.LAST_OPENED_CURVE(I).EQ."N")THEN           
                DO J=1,SUBDIVISION(I)
                    CONT1=CONT1+1
                    CONT2=CONT2+1
                    BOUNDARY_LINES(CONT1,1)=CONT1
                    BOUNDARY_LINES(CONT1,2)=CONT1+1      
                ENDDO                   
            ELSE IF(OPENED_CLOSED_CURVE(I).EQ."O".AND.LAST_OPENED_CURVE(I).EQ."S")THEN                 
                DO J=1,SUBDIVISION(I)
                    CONT1=CONT1+1
                    CONT2=CONT2+1
                    IF(J.LT.SUBDIVISION(I))THEN
                        BOUNDARY_LINES(CONT1,1)=CONT1
                        BOUNDARY_LINES(CONT1,2)=CONT1+1 
                    ELSE
                        BOUNDARY_LINES(CONT1,1)=CONT1
                        BOUNDARY_LINES(CONT1,2)=CONT1+1-CONT2 
                    ENDIF
                ENDDO    
                CONT2=0                                         
            ELSE IF(OPENED_CLOSED_CURVE(I).EQ."C")THEN
                CONT2=0
                DO J=1,SUBDIVISION(I)
                    CONT1=CONT1+1
                    CONT2=CONT2+1
                    IF(J.LT.SUBDIVISION(I))THEN
                        BOUNDARY_LINES(CONT1,1)=CONT1
                        BOUNDARY_LINES(CONT1,2)=CONT1+1 
                    ELSE
                        BOUNDARY_LINES(CONT1,1)=CONT1
                        BOUNDARY_LINES(CONT1,2)=CONT1+1-CONT2 
                    ENDIF    
                ENDDO      
            ENDIF
        ENDIF                 
    ENDDO        
!    
    DO I=1,N_PARAMETRICNODES

        SIZE=ABS(KNOT_VECTORS(SURFACE,1,1)-KNOT_VECTORS(SURFACE,1,(N_BF(SURFACE,1)+POLYNOMIAL_ORDERS(SURFACE,1)+1)))
        D=ABS((PARAMETRICCOORD_NODES(I,2)-(KNOT_VECTORS(SURFACE,2,1)))/SIZE)
        IF(D.LT.TOL)THEN
            SIGNED_DIST_PARAMETRICNODES(I)=0.D0
        ENDIF 
        
        SIZE=ABS(KNOT_VECTORS(SURFACE,2,1)-KNOT_VECTORS(SURFACE,2,(N_BF(SURFACE,2)+POLYNOMIAL_ORDERS(SURFACE,2)+1)))
        D=ABS((PARAMETRICCOORD_NODES(I,1)-(KNOT_VECTORS(SURFACE,1,(N_BF(SURFACE,1)+POLYNOMIAL_ORDERS(SURFACE,1)+1))))/SIZE)
        IF(D.LT.TOL)THEN
            SIGNED_DIST_PARAMETRICNODES(I)=0.D0
        ENDIF
        
        SIZE=ABS(KNOT_VECTORS(SURFACE,1,1)-KNOT_VECTORS(SURFACE,1,(N_BF(SURFACE,1)+POLYNOMIAL_ORDERS(SURFACE,1)+1)))
        D=ABS((PARAMETRICCOORD_NODES(I,2)-(KNOT_VECTORS(SURFACE,2,(N_BF(SURFACE,2)+POLYNOMIAL_ORDERS(SURFACE,2)+1))))/SIZE)
        IF(D.LT.TOL)THEN
            SIGNED_DIST_PARAMETRICNODES(I)=0.D0
        ENDIF 
        
        SIZE=ABS(KNOT_VECTORS(SURFACE,2,1)-KNOT_VECTORS(SURFACE,2,(N_BF(SURFACE,2)+POLYNOMIAL_ORDERS(SURFACE,2)+1)))
        D=ABS((PARAMETRICCOORD_NODES(I,1)-(KNOT_VECTORS(SURFACE,1,1)))/SIZE)
        IF(D.LT.TOL)THEN
            SIGNED_DIST_PARAMETRICNODES(I)=0.D0
        ENDIF                                
        
        DO J=1,N_BOUNDARYLINES
            IF((BOUNDARY_POINTS(BOUNDARY_LINES(J,1),1)-BOUNDARY_POINTS(BOUNDARY_LINES(J,2),1)).NE.0.D0)THEN
                A=(BOUNDARY_POINTS(BOUNDARY_LINES(J,1),2)-BOUNDARY_POINTS(BOUNDARY_LINES(J,2),2))/(BOUNDARY_POINTS(BOUNDARY_LINES(J,1),1)-BOUNDARY_POINTS(BOUNDARY_LINES(J,2),1))
                B=BOUNDARY_POINTS(BOUNDARY_LINES(J,1),2)-A*BOUNDARY_POINTS(BOUNDARY_LINES(J,1),1)
                SIZE=ABS(BOUNDARY_POINTS(BOUNDARY_LINES(J,1),1)-BOUNDARY_POINTS(BOUNDARY_LINES(J,2),1))+ABS(BOUNDARY_POINTS(BOUNDARY_LINES(J,1),2)-BOUNDARY_POINTS(BOUNDARY_LINES(J,2),2))
                D=ABS((PARAMETRICCOORD_NODES(I,2)-(A*PARAMETRICCOORD_NODES(I,1)+B))/SIZE)
                IF(D.LT.TOL)THEN
                    IF(BOUNDARY_POINTS(BOUNDARY_LINES(J,1),1).LT.BOUNDARY_POINTS(BOUNDARY_LINES(J,2),1))THEN
                        IF(BOUNDARY_POINTS(BOUNDARY_LINES(J,1),1)*(1-TOL).LE.PARAMETRICCOORD_NODES(I,1).AND.PARAMETRICCOORD_NODES(I,1).LE.BOUNDARY_POINTS(BOUNDARY_LINES(J,2),1)*(1+TOL))THEN
                            SIGNED_DIST_PARAMETRICNODES(I)=0.D0
                        ENDIF    
                    ELSE
                         IF(BOUNDARY_POINTS(BOUNDARY_LINES(J,2),1)*(1-TOL).LE.PARAMETRICCOORD_NODES(I,1).AND.PARAMETRICCOORD_NODES(I,1).LE.BOUNDARY_POINTS(BOUNDARY_LINES(J,1),1)*(1+TOL))THEN
                            SIGNED_DIST_PARAMETRICNODES(I)=0.D0
                        ENDIF                     
                    ENDIF
                ENDIF
            ELSE    
                SIZE=ABS(BOUNDARY_POINTS(BOUNDARY_LINES(J,1),1)-BOUNDARY_POINTS(BOUNDARY_LINES(J,2),1))+ABS(BOUNDARY_POINTS(BOUNDARY_LINES(J,1),2)-BOUNDARY_POINTS(BOUNDARY_LINES(J,2),2))
                D=ABS((PARAMETRICCOORD_NODES(I,1)-(BOUNDARY_POINTS(BOUNDARY_LINES(J,1),1)))/SIZE)
                IF(D.LT.TOL)THEN
                    IF(BOUNDARY_POINTS(BOUNDARY_LINES(J,1),2).LT.BOUNDARY_POINTS(BOUNDARY_LINES(J,2),2))THEN
                        IF(BOUNDARY_POINTS(BOUNDARY_LINES(J,1),2)*(1-TOL).LE.PARAMETRICCOORD_NODES(I,2).AND.PARAMETRICCOORD_NODES(I,2).LE.BOUNDARY_POINTS(BOUNDARY_LINES(J,2),2)*(1+TOL))THEN
                            SIGNED_DIST_PARAMETRICNODES(I)=0.D0
                        ENDIF    
                    ELSE
                         IF(BOUNDARY_POINTS(BOUNDARY_LINES(J,2),2)*(1-TOL).LE.PARAMETRICCOORD_NODES(I,2).AND.PARAMETRICCOORD_NODES(I,2).LE.BOUNDARY_POINTS(BOUNDARY_LINES(J,1),2)*(1+TOL))THEN
                            SIGNED_DIST_PARAMETRICNODES(I)=0.D0
                        ENDIF                     
                    ENDIF
                ENDIF           
            ENDIF     
        ENDDO
    ENDDO 
!
    DO I=1,N_PARAMETRICNODES    
        IF(SIGNED_DIST_PARAMETRICNODES(I).EQ.0.D0)THEN
            DO J=1,N_PARAMETRICELEMENTS
                IF(PARAMETRIC_CONNECTIVITY(J,4).EQ.0.OR.PARAMETRIC_CONNECTIVITY(J,3).EQ.PARAMETRIC_CONNECTIVITY(J,4))THEN
                    DO K=1,3
                        IF(PARAMETRIC_CONNECTIVITY(J,K).EQ.I)THEN
                            COMMON_ELEMENTS(I)=COMMON_ELEMENTS(I)+1
                        ENDIF
                    ENDDO              
                ELSE
                    DO K=1,4
                        IF(PARAMETRIC_CONNECTIVITY(J,K).EQ.I)THEN
                            COMMON_ELEMENTS(I)=COMMON_ELEMENTS(I)+1
                        ENDIF
                    ENDDO
                ENDIF    
            ENDDO  
        ENDIF
    ENDDO        
!            
	END SUBROUTINE    