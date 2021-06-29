	SUBROUTINE BOUNDARY_CRACK_PROPAGATION
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
	INTEGER::I,J,K,II,NEN
!
    REAL*8::X,Y,Z,R,P(3),T1(3),T2(3),X1(3),X3(3),V1(3),V2(3),V3(3),V4(3),a,b,c,d,e,A1,A2,CRACKTIPQSI(2),PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),VALUES[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:)
    REAL*8::L,BOUNDARYCRACK_TOL	
    BOUNDARYCRACK_TOL=1.D-2 
!
    DO I=1,N_ELEM       
        !------------------------------------------------------------------------------------------------------------------------------------------------------
        !   DEFINING NEW COORDINATES FOR THE NEW ELEMENTS AT THE BOUNDARY CRACKTIP 
        !------------------------------------------------------------------------------------------------------------------------------------------------------ 
        IF(CRACKED_ELEM(I).EQ."CRACKTIP")THEN 
            X1=COORD_NODES(NODES_CONNECTIVITY(I,1),:)
            X3=COORD_NODES(NODES_CONNECTIVITY(I,3),:)
            V1=COORD_NODES(NODES_CONNECTIVITY(I,2),:)-COORD_NODES(NODES_CONNECTIVITY(I,1),:) 
            V2=COORD_NODES(NODES_CONNECTIVITY(I,4),:)-COORD_NODES(NODES_CONNECTIVITY(I,1),:)      
            V3=COORD_NODES(NODES_CONNECTIVITY(I,2),:)-COORD_NODES(NODES_CONNECTIVITY(I,3),:) 
            V4=COORD_NODES(NODES_CONNECTIVITY(I,4),:)-COORD_NODES(NODES_CONNECTIVITY(I,3),:)
            L=1.D0/4.D0*(DSQRT(V1(1)**2+V1(2)**2+V1(3)**2)+DSQRT(V2(1)**2+V2(2)**2+V2(3)**2)+DSQRT(V3(1)**2+V3(2)**2+V3(3)**2)+DSQRT(V4(1)**2+V4(2)**2+V4(3)**2))         
            DO J=1,N_ELEM
                IF(DUAL_BEM(COLLOCPOINTS_CONNECTIVITY(J,1)).EQ."S")THEN
                    NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                    DO K=1,NEN
                        IF(CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(J,K)).EQ.1)THEN
                            P=COORD_NODES(NODES_CONNECTIVITY(J,K),:)
                            
                            a=V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)
                            b=V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
                            c=V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3)
                            d=V1(1)*(X1(1)-P(1))+V1(2)*(X1(2)-P(2))+V1(3)*(X1(3)-P(3))
                            e=V2(1)*(X1(1)-P(1))+V2(2)*(X1(2)-P(2))+V2(3)*(X1(3)-P(3))
                            A1=(b*e-c*d)/(a*c-b**2) 
                            A2=(b*d-a*e)/(a*c-b**2)
                            IF(0.d0.LE.A1.AND.A1.LE.1.D0)THEN
                                IF(0.d0.LE.A2.AND.A2.LE.1.D0)THEN
                                    IF(A1+A2.LE.1.D0)THEN
                                        T1=X1+V1*A1+V2*A2
                                        R=DSQRT((T1(1)-P(1))**2+(T1(2)-P(2))**2+(T1(3)-P(3))**2)
                                        IF(R/(L).LE.BOUNDARYCRACK_TOL)THEN
                                            COORD_BOUNDARYCRACKTIP(I,:)=P 
                                            CRACKTIPQSI(1)=-1+2*A1	         
                                            CRACKTIPQSI(2)=-1+2*A2
                                        ENDIF    
                                    ENDIF
                                ENDIF
                            ENDIF
                            a=V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3)
                            b=V3(1)*V4(1)+V3(2)*V4(2)+V3(3)*V4(3)
                            c=V4(1)*V4(1)+V4(2)*V4(2)+V4(3)*V4(3)
                            d=V3(1)*(X3(1)-P(1))+V3(2)*(X3(2)-P(2))+V3(3)*(X3(3)-P(3))
                            e=V4(1)*(X3(1)-P(1))+V4(2)*(X3(2)-P(2))+V4(3)*(X3(3)-P(3))
                            A1=(b*d-a*e)/(a*c-b**2)                                
                            A2=(b*e-c*d)/(a*c-b**2) 
                            IF(0.d0.LE.A1.AND.A1.LE.1.D0)THEN
                                IF(0.d0.LE.A2.AND.A2.LE.1.D0)THEN
                                    IF(A1+A2.LE.1.D0)THEN
                                        T2=X3+V3*A2+V4*A1
                                        R=DSQRT((T2(1)-P(1))**2+(T2(2)-P(2))**2+(T2(3)-P(3))**2)
                                        IF(R/(L).LE.BOUNDARYCRACK_TOL)THEN
                                            COORD_BOUNDARYCRACKTIP(I,:)=P 
                                            CRACKTIPQSI(1)=1-2*A1	         
                                            CRACKTIPQSI(2)=1-2*A2
                                        ENDIF    
                                    ENDIF
                                ENDIF
                            ENDIF                                                
                        ENDIF
                    ENDDO
                ENDIF
            ENDDO          
    ! 
            NEN=ELEM_TYPE(I)*ORDER_ELEM(I)+(ELEM_TYPE(I)-3)*(ORDER_ELEM(I)-1)*POL_FAMILY(I)
            ALLOCATE(VALUES(NEN,3))
            DO J=1,NEN
	            VALUES(J,1)=COORD_NODES(NODES_CONNECTIVITY(I,J),1)
	            VALUES(J,2)=COORD_NODES(NODES_CONNECTIVITY(I,J),2)
	            VALUES(J,3)=COORD_NODES(NODES_CONNECTIVITY(I,J),3)
            ENDDO
    !
            SELECT CASE(ELEM_TYPE(I))
            CASE(3)
        
            CASE(4)	  	
                SELECT CASE(CRACKED_EDGES(I,1))
                CASE(1)
                    CALL SHAPE_FUNCTIONS(1,1.D0,CRACKTIPQSI(2),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)    
                    COORD_BOUNDARYCRACKTIPSTOP(I,1,1)=X
                    COORD_BOUNDARYCRACKTIPSTOP(I,1,2)=Y
                    COORD_BOUNDARYCRACKTIPSTOP(I,1,3)=Z                        
                    CALL SHAPE_FUNCTIONS(1,-1.D0,CRACKTIPQSI(2),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)  
                    COORD_BOUNDARYCRACKTIPSTOP(I,2,1)=X
                    COORD_BOUNDARYCRACKTIPSTOP(I,2,2)=Y
                    COORD_BOUNDARYCRACKTIPSTOP(I,2,3)=Z                 
                CASE(2)
                    CALL SHAPE_FUNCTIONS(1,CRACKTIPQSI(1),1.D0,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)    
                    COORD_BOUNDARYCRACKTIPSTOP(I,1,1)=X
                    COORD_BOUNDARYCRACKTIPSTOP(I,1,2)=Y
                    COORD_BOUNDARYCRACKTIPSTOP(I,1,3)=Z                        
                    CALL SHAPE_FUNCTIONS(1,CRACKTIPQSI(1),-1.D0,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)  
                    COORD_BOUNDARYCRACKTIPSTOP(I,2,1)=X
                    COORD_BOUNDARYCRACKTIPSTOP(I,2,2)=Y
                    COORD_BOUNDARYCRACKTIPSTOP(I,2,3)=Z     
                CASE(3)
                    CALL SHAPE_FUNCTIONS(1,-1.D0,CRACKTIPQSI(2),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)    
                    COORD_BOUNDARYCRACKTIPSTOP(I,1,1)=X
                    COORD_BOUNDARYCRACKTIPSTOP(I,1,2)=Y
                    COORD_BOUNDARYCRACKTIPSTOP(I,1,3)=Z                        
                    CALL SHAPE_FUNCTIONS(1,1.D0,CRACKTIPQSI(2),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)  
                    COORD_BOUNDARYCRACKTIPSTOP(I,2,1)=X
                    COORD_BOUNDARYCRACKTIPSTOP(I,2,2)=Y
                    COORD_BOUNDARYCRACKTIPSTOP(I,2,3)=Z     
                CASE(4)
                    CALL SHAPE_FUNCTIONS(1,CRACKTIPQSI(1),-1.D0,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)    
                    COORD_BOUNDARYCRACKTIPSTOP(I,1,1)=X
                    COORD_BOUNDARYCRACKTIPSTOP(I,1,2)=Y
                    COORD_BOUNDARYCRACKTIPSTOP(I,1,3)=Z                        
                    CALL SHAPE_FUNCTIONS(1,CRACKTIPQSI(1),1.D0,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)  
                    COORD_BOUNDARYCRACKTIPSTOP(I,2,1)=X
                    COORD_BOUNDARYCRACKTIPSTOP(I,2,2)=Y
                    COORD_BOUNDARYCRACKTIPSTOP(I,2,3)=Z         
                ENDSELECT        
            ENDSELECT
!                
            DEALLOCATE(VALUES)
        ENDIF        
    ENDDO
!     	   
    END SUBROUTINE BOUNDARY_CRACK_PROPAGATION