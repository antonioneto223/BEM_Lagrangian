	SUBROUTINE MINIMUM_DISTANCE_POINT_ELEMENT(CRACKTIP,ELEMENT,Dmin,A1,A2)
!
	USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
	USE REMESHING
	USE PROPAGATION
!
	IMPLICIT NONE 
!
	INTEGER::CRACKTIP,ELEMENT,J,K,II,NEN
!
    REAL*8::Dmin,Dmin1,Dmin2,A1,A2,A1_AUX,A2_AUX,X,Y,Z,R,P(3),T1(3),T2(3),X1(3),X3(3),V1(3),V2(3),V3(3),V4(3),a,b,c,d,e,PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),VALUES[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:)	
!
    SELECT CASE(ELEM_TYPE(ELEMENT))
    CASE(3)
        X1=COORD_NODES(NODES_CONNECTIVITY(ELEMENT,3),:)
        V1=COORD_NODES(NODES_CONNECTIVITY(ELEMENT,1),:)-COORD_NODES(NODES_CONNECTIVITY(ELEMENT,3),:)
        V2=COORD_NODES(NODES_CONNECTIVITY(ELEMENT,2),:)-COORD_NODES(NODES_CONNECTIVITY(ELEMENT,3),:)       
    !
        P=COORD_CRACKTIP(CRACKTIP,:)
    !            
        a=V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)
        b=V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
        c=V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3)
        d=V1(1)*(X1(1)-P(1))+V1(2)*(X1(2)-P(2))+V1(3)*(X1(3)-P(3))
        e=V2(1)*(X1(1)-P(1))+V2(2)*(X1(2)-P(2))+V2(3)*(X1(3)-P(3))
        A1_AUX=(b*e-c*d)/(a*c-b**2) 
        A2_AUX=(b*d-a*e)/(a*c-b**2)
        IF(0.d0.LE.A1_AUX.AND.A1_AUX.LE.1.D0)THEN
            IF(0.d0.LE.A2_AUX.AND.A2_AUX.LE.1.D0)THEN
                IF(A1_AUX+A2_AUX.LE.1.D0)THEN
                    T1=X1+V1*A1_AUX+V2*A2_AUX
                    Dmin=DSQRT((T1(1)-P(1))**2+(T1(2)-P(2))**2+(T1(3)-P(3))**2)
                    A1=A1_AUX  
                    A2=A2_AUX                  
                ENDIF
            ENDIF
        ENDIF   
    CASE(4)
        X1=COORD_NODES(NODES_CONNECTIVITY(ELEMENT,1),:)
        X3=COORD_NODES(NODES_CONNECTIVITY(ELEMENT,3),:)
        V1=COORD_NODES(NODES_CONNECTIVITY(ELEMENT,2),:)-COORD_NODES(NODES_CONNECTIVITY(ELEMENT,1),:)
        V2=COORD_NODES(NODES_CONNECTIVITY(ELEMENT,4),:)-COORD_NODES(NODES_CONNECTIVITY(ELEMENT,1),:)      
        V3=COORD_NODES(NODES_CONNECTIVITY(ELEMENT,2),:)-COORD_NODES(NODES_CONNECTIVITY(ELEMENT,3),:) 
        V4=COORD_NODES(NODES_CONNECTIVITY(ELEMENT,4),:)-COORD_NODES(NODES_CONNECTIVITY(ELEMENT,3),:)  
    !
        P=COORD_CRACKTIP(CRACKTIP,:)
    !            
        a=V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)
        b=V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
        c=V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3)
        d=V1(1)*(X1(1)-P(1))+V1(2)*(X1(2)-P(2))+V1(3)*(X1(3)-P(3))
        e=V2(1)*(X1(1)-P(1))+V2(2)*(X1(2)-P(2))+V2(3)*(X1(3)-P(3))
        A1_AUX=(b*e-c*d)/(a*c-b**2) 
        A2_AUX=(b*d-a*e)/(a*c-b**2)
        IF(0.d0.LE.A1_AUX.AND.A1_AUX.LE.1.D0)THEN
            IF(0.d0.LE.A2_AUX.AND.A2_AUX.LE.1.D0)THEN
                IF(A1_AUX+A2_AUX.LE.1.D0)THEN
                    T1=X1+V1*A1_AUX+V2*A2_AUX
                    Dmin1=DSQRT((T1(1)-P(1))**2+(T1(2)-P(2))**2+(T1(3)-P(3))**2)
                    T1=X1+V1*A2_AUX+V2*A1_AUX
                    Dmin2=DSQRT((T1(1)-P(1))**2+(T1(2)-P(2))**2+(T1(3)-P(3))**2)                    
                    IF(Dmin1.LT.Dmin2)THEN
                        Dmin=Dmin1
                        A1=-1+2*A1_AUX	         
                        A2=-1+2*A2_AUX
                    ELSE
                        Dmin=Dmin2                    
                        A1=-1+2*A2_AUX	         
                        A2=-1+2*A1_AUX                   
                    ENDIF                                     
                ENDIF
            ENDIF
        ENDIF
        a=V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3)
        b=V3(1)*V4(1)+V3(2)*V4(2)+V3(3)*V4(3)
        c=V4(1)*V4(1)+V4(2)*V4(2)+V4(3)*V4(3)
        d=V3(1)*(X3(1)-P(1))+V3(2)*(X3(2)-P(2))+V3(3)*(X3(3)-P(3))
        e=V4(1)*(X3(1)-P(1))+V4(2)*(X3(2)-P(2))+V4(3)*(X3(3)-P(3))
        A1_AUX=(b*e-c*d)/(a*c-b**2)         
        A2_AUX=(b*d-a*e)/(a*c-b**2)        
        IF(0.d0.LE.A1_AUX.AND.A1_AUX.LE.1.D0)THEN
            IF(0.d0.LE.A2_AUX.AND.A2_AUX.LE.1.D0)THEN
                IF(A1_AUX+A2_AUX.LE.1.D0)THEN
                    T2=X3+V3*A1_AUX+V4*A2_AUX
                    Dmin1=DSQRT((T2(1)-P(1))**2+(T2(2)-P(2))**2+(T2(3)-P(3))**2)
                    T2=X3+V3*A2_AUX+V4*A1_AUX
                    Dmin2=DSQRT((T2(1)-P(1))**2+(T2(2)-P(2))**2+(T2(3)-P(3))**2)                    
                    IF(Dmin1.LT.Dmin2)THEN
                        Dmin=Dmin1
                        A1=1-2*A2_AUX	         
                        A2=1-2*A1_AUX
                    ELSE
                        Dmin=Dmin2                    
                        A1=1-2*A1_AUX	         
                        A2=1-2*A2_AUX                    
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
    END SELECT                                                     
!     	   
    END SUBROUTINE MINIMUM_DISTANCE_POINT_ELEMENT