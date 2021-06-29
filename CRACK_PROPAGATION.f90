	SUBROUTINE CRACK_PROPAGATION
!	
    USE ISOPARAMETRIC_MESH
    USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
	USE REMESHING
	USE PROPAGATION	   
	USE FATIGUE
!
    IMPLICIT NONE
    INTEGER::I,J,K,II,JJ,KK,ELEMENT,LOCAL_NUMBER,LOCAL_NUMBER1,LOCAL_NUMBER2,NEN,CONT
    INTEGER::ELEM_TYPE_AUX[ALLOCATABLE](:),ORDER_ELEM_AUX[ALLOCATABLE](:),POL_FAMILY_AUX[ALLOCATABLE](:),NODES_CONNECTIVITY_AUX[ALLOCATABLE](:,:),COLLOCPOINTS_CONNECTIVITY_AUX[ALLOCATABLE](:,:),ND_AUX[ALLOCATABLE](:),CRACKED_EDGES_AUX[ALLOCATABLE](:,:),B_CONDITIONS_AUX[ALLOCATABLE](:),CRACK_TIP_AUX[ALLOCATABLE](:),NEW_ELEMENTS_AUX[ALLOCATABLE](:,:),EDGE_DISCONTINUITY_AUX[ALLOCATABLE](:,:) 
!    
    REAL*8::PI,GTHETA,GTHETA_AUX,THETA,E,Nu,KIC,Dmin,Dmin_AUX,A1,A2,A1_AUX,A2_AUX,DX,DY,DZ,R,F
    REAL*8::ETAX(3),ETAY(3),ETAZ(3),DIRECTION1(3),DIRECTION2(3),R_S(3,3),R_H(3,3),H_LOCAL(3),G_LOCAL(3),HG_AUX[ALLOCATABLE](:,:),UT_AUX[ALLOCATABLE](:)  
    REAL*8::VALUES[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:),PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:) 
    REAL*8::COORD_NODES_AUX[ALLOCATABLE](:,:),COORD_COLLOCPOINTS_AUX[ALLOCATABLE](:,:),QSI_AUX[ALLOCATABLE](:,:,:),COORD_INTERCEPTPOINTS_AUX[ALLOCATABLE](:,:,:),COORD_EXTRAINTERCEPTPOINTS_AUX[ALLOCATABLE](:,:),COORD_BOUNDARYCRACKTIP_AUX[ALLOCATABLE](:,:),COORD_BOUNDARYCRACKTIPSTOP_AUX[ALLOCATABLE](:,:,:),NC_AUX[ALLOCATABLE](:,:)
    REAL*8::X,Y,Z
!
    CHARACTER*30::DISCONTINUOUS_COLLOCPOINTS_AUX[ALLOCATABLE](:),CRACKED_ELEM_AUX[ALLOCATABLE](:),CHANGEABLE_ELEM_AUX[ALLOCATABLE](:),DUAL_BEM_AUX[ALLOCATABLE](:),EQ_TYPE_AUX[ALLOCATABLE](:),CHANGEABLE_COLLOCPOINTS_AUX[ALLOCATABLE](:)  
!	
!------------------------------------------------------------------------------------------------------------------------------------------------------
!         H AND G GLOBAL COORDINATE SYSTEM
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
	        DO J=1,3*N_COLLOCPOINTS
	            H_LOCAL=0.D0
	            G_LOCAL=0.D0
	            DO K=1,3	            
	                H_LOCAL(K)=H(J,3*I-2)*R_S(1,K)+H(J,3*I-1)*R_S(2,K)+H(J,3*I)*R_S(3,K)
	                G_LOCAL(K)=G(J,3*I-2)*R_S(1,K)+G(J,3*I-1)*R_S(2,K)+G(J,3*I)*R_S(3,K)
	            ENDDO    
                H(J,3*I-2)=H_LOCAL(1)
                H(J,3*I-1)=H_LOCAL(2)
                H(J,3*I)=H_LOCAL(3)
!
                G(J,3*I-2)=G_LOCAL(1)
                G(J,3*I-1)=G_LOCAL(2)
                G(J,3*I)=G_LOCAL(3)                              
	        ENDDO
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
	        DO J=1,3*N_COLLOCPOINTS
	            H_LOCAL=0.D0
	            G_LOCAL=0.D0
	            DO K=1,3	            
	                H_LOCAL(K)=H(J,3*II-2)*R_H(1,K)+H(J,3*II-1)*R_H(2,K)+H(J,3*II)*R_H(3,K)
	                G_LOCAL(K)=G(J,3*II-2)*R_H(1,K)+G(J,3*II-1)*R_H(2,K)+G(J,3*II)*R_H(3,K)
	            ENDDO    
                H(J,3*II-2)=H_LOCAL(1)
                H(J,3*II-1)=H_LOCAL(2)
                H(J,3*II)=H_LOCAL(3)
!
                G(J,3*II-2)=G_LOCAL(1)
                G(J,3*II-1)=G_LOCAL(2)
                G(J,3*II)=G_LOCAL(3)                              
	        ENDDO		    		    
	    ENDIF    
	ENDDO	      
!
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   COMPUTING THE PROPAGATION ANGLES AND THE EQUIVALENT STRESS INTENSITY FACTORS 
!------------------------------------------------------------------------------------------------------------------------------------------------------
!  
    PI=DACOS(-1.D0)
    CONT=0    
    DO I=1,N_COLLOCPOINTS
        IF(DUAL_BEM(I).EQ."S")THEN
            IF(CRACK_TIP(I).EQ.1.OR.CRACK_TIP(I).EQ.3)THEN
                CONT=CONT+1 
                DO J=1,N_ELEM
                    NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                    DO K=1,NEN
                        IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
                            ELEMENT=J
                            LOCAL_NUMBER=K                    
                        ENDIF
                    ENDDO
                ENDDO  
                E=EMP(ND(ELEMENT),1)
                Nu=EMP(ND(ELEMENT),2)

                SELECT CASE(PROPAGATION_CRITERIA)  
                CASE("MERR") ! MAXIMUM ENERGY RELEASE RATE CRITERIA                        
                    CALL MAXIMUM_ENERGY_RELEASE_RATE(Nu,E,KI(CONT),KII(CONT),KIII(CONT),THETAP(CONT))
                    PHIP(CONT)=0.D0  
                    Gmax(CONT)=(1+Nu)/E*COS(THETAP(CONT))**2*((1-Nu**2)/(2*(1+Nu))*(KI(CONT)**2*(1+COS(THETAP(CONT)))-4*KI(CONT)*KII(CONT)*SIN(THETAP(CONT))+KII(CONT)**2*(5-3*COS(THETAP(CONT))))+KIII(CONT)**2) 
                    Keq(CONT)=DSQRT(E*Gmax(CONT)/(1-Nu**2))
                CASE("MPS") ! MAXIMUM PRINCIPAL STRESS CRITERIA
                    CALL MAXIMUM_PRINCIPAL_STRESS(KI(CONT),KII(CONT),KIII(CONT),THETAP(CONT))
                    F=(1-Nu)*KI(CONT)*(3*COS(THETAP(CONT)/2.D0)+COS(3*THETAP(CONT)/2.D0))/4.D0-(1-Nu)*KII(CONT)*(3*SIN(THETAP(CONT)/2.D0)+3*SIN(3*THETAP(CONT)/2.D0))/4.D0-Nu*KI(CONT)*(5*COS(THETAP(CONT)/2.D0)-COS(3*THETAP(CONT)/2.D0))/4.D0+Nu*KII(CONT)*(5*SIN(THETAP(CONT)/2.D0)-3*SIN(3*THETAP(CONT)/2.D0))/4.D0
                    PHIP(CONT)=0.5D0*ATAN(2*KIII(CONT)*COS(THETAP(CONT)/2.D0)/F)
                    Keq(CONT)=0.5D0*COS(THETAP(CONT)/2.D0)*(KI(CONT)*COS(THETAP(CONT))**2-1.5D0*KII(CONT)*SIN(THETAP(CONT))+DSQRT((KI(CONT)*COS(THETAP(CONT)/2.D0)**2-1.5D0*KII(CONT)*SIN(THETAP(CONT)))**2+4*KIII(CONT)**2)) 
                ENDSELECT                      
            ENDIF
        ENDIF         
    ENDDO   
!
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   CHECKING BRITTLE FRACTURE CRITERIA 
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
    CONT=0
    Keqmax=0.D0
    PROPAGATE=.FALSE.
    DO I=1,N_COLLOCPOINTS
        IF(CRACK_TIP(I).EQ.1.AND.DUAL_BEM(I).EQ."S")THEN
            CONT=CONT+1 
            DO J=1,N_ELEM
                NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                DO K=1,NEN
                    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
                        ELEMENT=J
                        LOCAL_NUMBER=K                    
                    ENDIF
                ENDDO
            ENDDO 
            KIC=FFMP(ND(ELEMENT),1) 
            IF(Keq(CONT).GE.KIC)THEN
                PROPAGATE=.TRUE.
            ENDIF
            IF(Keq(CONT).GE.Keqmax)THEN
                Keqmax=Keq(CONT)
            ENDIF
        ELSE IF(CRACK_TIP(I).EQ.3.AND.DUAL_BEM(I).EQ."S")THEN
            CONT=CONT+1 
            DO J=1,N_ELEM
                NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                DO K=1,NEN
                    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
                        ELEMENT=J
                        LOCAL_NUMBER=K                    
                    ENDIF
                ENDDO
            ENDDO 
            KIC=FFMP(ND(ELEMENT),1) 
            IF(Keq(CONT).GE.KIC)THEN
                PROPAGATE=.TRUE.
            ENDIF
            IF(Keq(CONT).GE.Keqmax)THEN
                Keqmax=Keq(CONT)
            ENDIF
        ENDIF      
    ENDDO                                      
!
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   CREATING THE NEW CRACK TIP ELEMENTS 
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
    IF(PROPAGATE)THEN
        N_CRACKTIPELEMENTS=N_CRACKTIPS
        N_CRACKTIPCOLLOCPOINTS=4*N_CRACKTIPELEMENTS
        N_CRACKTIPNODES=N_CRACKTIPS+2
        N_INTPOINTS=2*N_CRACKTIPELEMENTS
        ALLOCATE(COORD_INTPOINTS(N_INTPOINTS,3),NDI(N_INTPOINTS)) 
        N_CRACK_POINTS=0
	    DO I=1,N_COLLOCPOINTS
            IF(DUAL_BEM(I).EQ.'S')THEN
                N_CRACK_POINTS=N_CRACK_POINTS+1
            ENDIF
	    ENDDO
!  
        CONT=0
        DO I=1,N_COLLOCPOINTS
            IF(DUAL_BEM(I).EQ."S")THEN
                IF(CRACK_TIP(I).EQ.1.OR.CRACK_TIP(I).EQ.3)THEN
                    CONT=CONT+1 
 	                IF(Keq(CONT).GE.0.05D0*Keqmax)THEN
 	                    CRACK_INCREMENT(CONT)=MAX_CRACKINCREMENT*Keq(CONT)/Keqmax
 	                ELSE
 	                    CRACK_INCREMENT(CONT)=MAX_CRACKINCREMENT*0.05D0
 	                ENDIF                    
                ENDIF
            ENDIF
        ENDDO 
!        
        CALL VERTICAL_INCREMENTS               
!   
        CONT=0
        DO I=1,N_COLLOCPOINTS
            IF(DUAL_BEM(I).EQ."S")THEN
                IF(CRACK_TIP(I).EQ.1.OR.CRACK_TIP(I).EQ.3)THEN
                    CONT=CONT+1 
                    DO J=1,N_ELEM
                        NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                        DO K=1,NEN
                            IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
                                ELEMENT=J
                                LOCAL_NUMBER=K                    
                            ENDIF
                        ENDDO
                    ENDDO        
 	                CALL LOCAL_COODINATE_SYSTEM(LOCAL_NUMBER,ELEMENT,R_S,R_H)
 	                ETAX(1)=R_S(1,1)
 	                ETAX(2)=R_S(1,2)
 	                ETAX(3)=R_S(1,3) 	
 	                ETAY(1)=R_S(2,1)
 	                ETAY(2)=R_S(2,2)
 	                ETAY(3)=R_S(2,3)
 	                ETAZ(1)=R_S(3,1)
 	                ETAZ(2)=R_S(3,2)
 	                ETAZ(3)=R_S(3,3)
 	                
 	                DIRECTION1=COS(THETAP(CONT))*ETAY+SIN(THETAP(CONT))*ETAX
 	                DIRECTION1=DIRECTION1/DSQRT(DIRECTION1(1)**2+DIRECTION1(2)**2+DIRECTION1(3)**2) 
 	                	             	                
 	                DIRECTION2=-SIN(THETAP(CONT))*ETAY+COS(THETAP(CONT))*ETAX
 	                DIRECTION2=DIRECTION2/DSQRT(DIRECTION2(1)**2+DIRECTION2(2)**2+DIRECTION2(3)**2) 
 	                     	                 	      	                 	        
 	                COORD_CRACKTIP(CONT,:)=COORD_CRACKTIP(CONT,:)+CRACK_INCREMENT(CONT)*DIRECTION1+VERTICAL_INCREMENT(CONT)*DIRECTION2       
 	            ENDIF
 	        ENDIF                        
        ENDDO
!        
        ALLOCATE(COORD_NODES_AUX(N_NODES,3),COORD_COLLOCPOINTS_AUX(N_COLLOCPOINTS,3),DISCONTINUOUS_COLLOCPOINTS_AUX(N_COLLOCPOINTS),QSI_AUX(N_ELEM,9,2),ELEM_TYPE_AUX(N_ELEM),ORDER_ELEM_AUX(N_ELEM),POL_FAMILY_AUX(N_ELEM),NODES_CONNECTIVITY_AUX(N_ELEM,4),COLLOCPOINTS_CONNECTIVITY_AUX(N_ELEM,4),ND_AUX(N_ELEM),EDGE_DISCONTINUITY_AUX(N_ELEM,4),CRACKED_EDGES_AUX(N_ELEM,2),NEW_ELEMENTS_AUX(N_ELEM,3),CRACKED_ELEM_AUX(N_ELEM),CHANGEABLE_ELEM_AUX(N_ELEM),COORD_INTERCEPTPOINTS_AUX(N_ELEM,4,3),COORD_EXTRAINTERCEPTPOINTS_AUX(N_ELEM,3),COORD_BOUNDARYCRACKTIP_AUX(N_ELEM,3),COORD_BOUNDARYCRACKTIPSTOP_AUX(N_ELEM,2,3),DUAL_BEM_AUX(N_COLLOCPOINTS),EQ_TYPE_AUX(N_COLLOCPOINTS),B_CONDITIONS_AUX(3*N_COLLOCPOINTS),CRACK_TIP_AUX(N_COLLOCPOINTS),CHANGEABLE_COLLOCPOINTS_AUX(N_COLLOCPOINTS),NC_AUX(N_CRACK_POINTS,2))  
        COORD_NODES_AUX=COORD_NODES
        COORD_COLLOCPOINTS_AUX=COORD_COLLOCPOINTS
        QSI_AUX=QSI
        DISCONTINUOUS_COLLOCPOINTS_AUX=DISCONTINUOUS_COLLOCPOINTS
        ELEM_TYPE_AUX=ELEM_TYPE
        ORDER_ELEM_AUX=ORDER_ELEM
        POL_FAMILY_AUX=POL_FAMILY
        NODES_CONNECTIVITY_AUX=NODES_CONNECTIVITY
        COLLOCPOINTS_CONNECTIVITY_AUX=COLLOCPOINTS_CONNECTIVITY
        ND_AUX=ND
        EDGE_DISCONTINUITY_AUX=EDGE_DISCONTINUITY        
        CRACKED_EDGES_AUX=CRACKED_EDGES 
        NEW_ELEMENTS_AUX=NEW_ELEMENTS         
        CRACKED_ELEM_AUX=CRACKED_ELEM
        CHANGEABLE_ELEM_AUX=CHANGEABLE_ELEM
        COORD_INTERCEPTPOINTS_AUX=COORD_INTERCEPTPOINTS
        COORD_EXTRAINTERCEPTPOINTS_AUX=COORD_EXTRAINTERCEPTPOINTS
        COORD_BOUNDARYCRACKTIP_AUX=COORD_BOUNDARYCRACKTIP
        COORD_BOUNDARYCRACKTIPSTOP_AUX=COORD_BOUNDARYCRACKTIPSTOP
        DUAL_BEM_AUX=DUAL_BEM 
        EQ_TYPE_AUX=EQ_TYPE   
        B_CONDITIONS_AUX=B_CONDITIONS
        CRACK_TIP_AUX=CRACK_TIP
        CHANGEABLE_COLLOCPOINTS_AUX=CHANGEABLE_COLLOCPOINTS
        NC_AUX=NC
        DEALLOCATE(COORD_NODES,COORD_COLLOCPOINTS,DISCONTINUOUS_COLLOCPOINTS,QSI,ELEM_TYPE,ORDER_ELEM,POL_FAMILY,NODES_CONNECTIVITY,COLLOCPOINTS_CONNECTIVITY,ND,EDGE_DISCONTINUITY,CRACKED_EDGES,NEW_ELEMENTS,CRACKED_ELEM,CHANGEABLE_ELEM,COORD_INTERCEPTPOINTS,COORD_EXTRAINTERCEPTPOINTS,COORD_BOUNDARYCRACKTIP,COORD_BOUNDARYCRACKTIPSTOP,DUAL_BEM,EQ_TYPE,B_CONDITIONS,CRACK_TIP,CHANGEABLE_COLLOCPOINTS,NC) 
        ALLOCATE(COORD_NODES(N_NODES+N_CRACKTIPNODES,3),COORD_COLLOCPOINTS(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS,3),DISCONTINUOUS_COLLOCPOINTS(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS),QSI(N_ELEM+N_CRACKTIPELEMENTS,9,2),ELEM_TYPE(N_ELEM+N_CRACKTIPELEMENTS),ORDER_ELEM(N_ELEM+N_CRACKTIPELEMENTS),POL_FAMILY(N_ELEM+N_CRACKTIPELEMENTS),NODES_CONNECTIVITY(N_ELEM+N_CRACKTIPELEMENTS,4),COLLOCPOINTS_CONNECTIVITY(N_ELEM+N_CRACKTIPELEMENTS,4),ND(N_ELEM+N_CRACKTIPELEMENTS),EDGE_DISCONTINUITY(N_ELEM+N_CRACKTIPELEMENTS,4),CRACKED_EDGES(N_ELEM+N_CRACKTIPELEMENTS,2),NEW_ELEMENTS(N_ELEM+N_CRACKTIPELEMENTS,3),CRACKED_ELEM(N_ELEM+N_CRACKTIPELEMENTS),CHANGEABLE_ELEM(N_ELEM+N_CRACKTIPELEMENTS),COORD_INTERCEPTPOINTS(N_ELEM+N_CRACKTIPELEMENTS,4,3),COORD_EXTRAINTERCEPTPOINTS(N_ELEM+N_CRACKTIPELEMENTS,3),COORD_BOUNDARYCRACKTIP(N_ELEM+N_CRACKTIPELEMENTS,3),COORD_BOUNDARYCRACKTIPSTOP(N_ELEM+N_CRACKTIPELEMENTS,2,3),DUAL_BEM(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS),EQ_TYPE(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS),B_CONDITIONS(3*(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS)),CRACK_TIP(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS),CHANGEABLE_COLLOCPOINTS(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS),NC(N_CRACK_POINTS+2*N_CRACKTIPELEMENTS,2))       
        COORD_NODES=0.D0  
        COORD_COLLOCPOINTS=0.D0
        QSI=0.D0 
        DISCONTINUOUS_COLLOCPOINTS="DISCONTINUOUS"
        ELEM_TYPE=4 
        ORDER_ELEM=1
        POL_FAMILY=1
        NODES_CONNECTIVITY=0
        COLLOCPOINTS_CONNECTIVITY=0
        ND=0
        EDGE_DISCONTINUITY=0        
        CRACKED_EDGES=0     
        NEW_ELEMENTS=0           
        CRACKED_ELEM="UNCRACKED"
        CHANGEABLE_ELEM="UNCHANGEABLE"
        COORD_INTERCEPTPOINTS=0.D0 
        COORD_EXTRAINTERCEPTPOINTS=0.D0  
        COORD_BOUNDARYCRACKTIP=0.D0
        COORD_BOUNDARYCRACKTIPSTOP=0.D0     
        DUAL_BEM="S"
        EQ_TYPE="S"
        B_CONDITIONS=1
        CRACK_TIP=0  
        CHANGEABLE_COLLOCPOINTS="UNCHANGEABLE"
        DO I=1,N_NODES
            DO J=1,3
                COORD_NODES(I,J)=COORD_NODES_AUX(I,J)
            ENDDO
        ENDDO
        DO I=1,N_COLLOCPOINTS
            DISCONTINUOUS_COLLOCPOINTS(I)=DISCONTINUOUS_COLLOCPOINTS_AUX(I)
            DUAL_BEM(I)=DUAL_BEM_AUX(I)
            EQ_TYPE(I)=EQ_TYPE_AUX(I)
            B_CONDITIONS(3*I-2)=B_CONDITIONS_AUX(3*I-2)
            B_CONDITIONS(3*I-1)=B_CONDITIONS_AUX(3*I-1)  
            B_CONDITIONS(3*I)=B_CONDITIONS_AUX(3*I)
            CRACK_TIP(I)=CRACK_TIP_AUX(I) 
            CHANGEABLE_COLLOCPOINTS(I)=CHANGEABLE_COLLOCPOINTS_AUX(I)             
            DO J=1,3
                COORD_COLLOCPOINTS(I,J)=COORD_COLLOCPOINTS_AUX(I,J) 
            ENDDO     
        ENDDO
        DO I=1,N_ELEM
            ELEM_TYPE(I)=ELEM_TYPE_AUX(I)
            ORDER_ELEM(I)=ORDER_ELEM_AUX(I)
            POL_FAMILY(I)=POL_FAMILY_AUX(I)
            ND(I)=ND_AUX(I)
            EDGE_DISCONTINUITY(I,:)=EDGE_DISCONTINUITY_AUX(I,:)            
            CRACKED_ELEM(I)=CRACKED_ELEM_AUX(I)
            CHANGEABLE_ELEM(I)=CHANGEABLE_ELEM_AUX(I)
            CRACKED_EDGES(I,1)=CRACKED_EDGES_AUX(I,1)
            CRACKED_EDGES(I,2)=CRACKED_EDGES_AUX(I,2) 
            NEW_ELEMENTS(I,1)=NEW_ELEMENTS_AUX(I,1)        
            NEW_ELEMENTS(I,2)=NEW_ELEMENTS_AUX(I,2)        
            NEW_ELEMENTS(I,3)=NEW_ELEMENTS_AUX(I,3)  
            COORD_EXTRAINTERCEPTPOINTS(I,1)=COORD_EXTRAINTERCEPTPOINTS_AUX(I,1)    
            COORD_EXTRAINTERCEPTPOINTS(I,2)=COORD_EXTRAINTERCEPTPOINTS_AUX(I,2) 
            COORD_EXTRAINTERCEPTPOINTS(I,3)=COORD_EXTRAINTERCEPTPOINTS_AUX(I,3)                                               
            COORD_BOUNDARYCRACKTIP(I,1)=COORD_BOUNDARYCRACKTIP_AUX(I,1)
            COORD_BOUNDARYCRACKTIP(I,2)=COORD_BOUNDARYCRACKTIP_AUX(I,2)
            COORD_BOUNDARYCRACKTIP(I,3)=COORD_BOUNDARYCRACKTIP_AUX(I,3)    
            COORD_BOUNDARYCRACKTIPSTOP(I,1,1)=COORD_BOUNDARYCRACKTIPSTOP_AUX(I,1,1)
            COORD_BOUNDARYCRACKTIPSTOP(I,1,2)=COORD_BOUNDARYCRACKTIPSTOP_AUX(I,1,2)
            COORD_BOUNDARYCRACKTIPSTOP(I,1,3)=COORD_BOUNDARYCRACKTIPSTOP_AUX(I,1,3)   
            COORD_BOUNDARYCRACKTIPSTOP(I,2,1)=COORD_BOUNDARYCRACKTIPSTOP_AUX(I,2,1)
            COORD_BOUNDARYCRACKTIPSTOP(I,2,2)=COORD_BOUNDARYCRACKTIPSTOP_AUX(I,2,2)
            COORD_BOUNDARYCRACKTIPSTOP(I,2,3)=COORD_BOUNDARYCRACKTIPSTOP_AUX(I,2,3) 
            DO J=1,9
                QSI(I,J,1)=QSI_AUX(I,J,1)
                QSI(I,J,2)=QSI_AUX(I,J,2)
            ENDDO                    
            DO J=1,4
                NODES_CONNECTIVITY(I,J)=NODES_CONNECTIVITY_AUX(I,J)
                COLLOCPOINTS_CONNECTIVITY(I,J)=COLLOCPOINTS_CONNECTIVITY_AUX(I,J)
                COORD_INTERCEPTPOINTS(I,J,1)=COORD_INTERCEPTPOINTS_AUX(I,J,1)
                COORD_INTERCEPTPOINTS(I,J,2)=COORD_INTERCEPTPOINTS_AUX(I,J,2)
                COORD_INTERCEPTPOINTS(I,J,3)=COORD_INTERCEPTPOINTS_AUX(I,J,3)            
            ENDDO
        ENDDO
        DO I=1,N_CRACK_POINTS
            NC(I,:)=NC_AUX(I,:)
        ENDDO
        DEALLOCATE(COORD_NODES_AUX,COORD_COLLOCPOINTS_AUX,QSI_AUX,DISCONTINUOUS_COLLOCPOINTS_AUX,ELEM_TYPE_AUX,ORDER_ELEM_AUX,POL_FAMILY_AUX,NODES_CONNECTIVITY_AUX,COLLOCPOINTS_CONNECTIVITY_AUX,ND_AUX,EDGE_DISCONTINUITY_AUX,CRACKED_EDGES_AUX,NEW_ELEMENTS_AUX,CRACKED_ELEM_AUX,CHANGEABLE_ELEM_AUX,COORD_INTERCEPTPOINTS_AUX,COORD_EXTRAINTERCEPTPOINTS_AUX,COORD_BOUNDARYCRACKTIP_AUX,COORD_BOUNDARYCRACKTIPSTOP_AUX,DUAL_BEM_AUX,EQ_TYPE_AUX,B_CONDITIONS_AUX,CRACK_TIP_AUX,CHANGEABLE_COLLOCPOINTS_AUX,NC_AUX) 
!        
!------------------------------------------------------------------------------------------------------------------------------------------------------
!       DEFINING THE COORDINATES OF THE NEW NODES 
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
        CONT=0
        DO I=1,N_CRACKTIPS
            IF(CRACK_TIP(CRACKTIP_COLLOCPOINTS(I)).EQ.3)THEN
                CONT=CONT+1
                Dmin=1.D14
                Dmin_AUX=1.D16
                DO J=1,N_ELEM
                    IF(DUAL_BEM(COLLOCPOINTS_CONNECTIVITY(J,1)).EQ."B")THEN
                        CALL MINIMUM_DISTANCE_POINT_ELEMENT(I,J,Dmin_AUX,A1_AUX,A2_AUX)
                        IF(Dmin_aux.LT.Dmin)THEN
                            Dmin=Dmin_aux    
                            A1=A1_AUX
                            A2=A2_AUX
                            ELEMENT=J  
                        ENDIF
                    ENDIF
                ENDDO
!                
                NEN=ELEM_TYPE(ELEMENT)*ORDER_ELEM(ELEMENT)+(ELEM_TYPE(ELEMENT)-3)*(ORDER_ELEM(ELEMENT)-1)*POL_FAMILY(ELEMENT)
                ALLOCATE(VALUES(NEN,3))
                DO J=1,NEN
		            VALUES(J,1)=COORD_NODES(NODES_CONNECTIVITY(ELEMENT,J),1)
		            VALUES(J,2)=COORD_NODES(NODES_CONNECTIVITY(ELEMENT,J),2)
		            VALUES(J,3)=COORD_NODES(NODES_CONNECTIVITY(ELEMENT,J),3)
	            ENDDO
	            CALL SHAPE_FUNCTIONS(1,A1,A2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
	            COORD_NODES(N_NODES+CONT,1)=X
	            COORD_NODES(N_NODES+CONT,2)=Y
	            COORD_NODES(N_NODES+CONT,3)=Z 
	            COORD_NODES(N_NODES+DINT((N_CRACKTIPS+2)/2.D0)+CONT,:)=COORD_NODES(N_NODES+CONT,:)
	            DEALLOCATE(VALUES)	            
            ELSE IF(CRACK_TIP(CRACKTIP_COLLOCPOINTS(I)).EQ.1)THEN
               IF(MOD(I,2).EQ.0.D0)THEN
                   CONT=CONT+1
                   IF(I.LT.N_CRACKTIPS)THEN
                       COORD_NODES(N_NODES+CONT,:)=(COORD_CRACKTIP(I,:)+COORD_CRACKTIP(I+1,:))/2.D0
                   ELSE
                       COORD_NODES(N_NODES+CONT,:)=(COORD_CRACKTIP(I,:)+COORD_CRACKTIP(I-1,:))/2.D0           
                   ENDIF  
	               COORD_NODES(N_NODES+DINT((N_CRACKTIPS+2)/2.D0)+CONT,:)=COORD_NODES(N_NODES+CONT,:)                  
               ENDIF
            ENDIF        
        ENDDO
!        
!------------------------------------------------------------------------------------------------------------------------------------------------------
!       DEFINING THE ELEMENTS OF THE NEW CRACKTIP 
!------------------------------------------------------------------------------------------------------------------------------------------------------
!
        CONT=0
        DO I=1,N_CRACKTIPS
            IF(MOD(I,2).NE.0.D0)THEN
                CONT=CONT+1
                DO J=1,N_ELEM
                    NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                    DO K=1,NEN
                        IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.CRACKTIP_COLLOCPOINTS(I))THEN
                            ELEMENT=J
                            LOCAL_NUMBER1=K
                        ENDIF
                        IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.CRACKTIP_COLLOCPOINTS(I+1))THEN
                            ELEMENT=J
                            LOCAL_NUMBER2=K
                        ENDIF                    
                    ENDDO
                ENDDO                  
!
                DO J=1,N_ELEM
                    NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                    DO K=1,NEN
                        DX=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),1)-COORD_COLLOCPOINTS(CRACKTIP_COLLOCPOINTS(I),1)
                        DY=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),2)-COORD_COLLOCPOINTS(CRACKTIP_COLLOCPOINTS(I),2)
                        DZ=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),3)-COORD_COLLOCPOINTS(CRACKTIP_COLLOCPOINTS(I),3)
                        R=DSQRT(DX**2+DY**2+DZ**2)
                        IF(R.LE.TOL.AND.DUAL_BEM(COLLOCPOINTS_CONNECTIVITY(J,K)).EQ."H")THEN
                            ELEMENT=J
                            LOCAL_NUMBER1=K                                                                        
                        ENDIF
                        DX=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),1)-COORD_COLLOCPOINTS(CRACKTIP_COLLOCPOINTS(I+1),1)
                        DY=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),2)-COORD_COLLOCPOINTS(CRACKTIP_COLLOCPOINTS(I+1),2)
                        DZ=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),3)-COORD_COLLOCPOINTS(CRACKTIP_COLLOCPOINTS(I+1),3)
                        R=DSQRT(DX**2+DY**2+DZ**2)
                        IF(R.LE.TOL.AND.DUAL_BEM(COLLOCPOINTS_CONNECTIVITY(J,K)).EQ."H")THEN
                            ELEMENT=J
                            LOCAL_NUMBER2=K                                                                        
                        ENDIF                        
                    ENDDO
                ENDDO
!
                NODES_CONNECTIVITY(N_ELEM+CONT,1)=N_NODES+CONT
                NODES_CONNECTIVITY(N_ELEM+CONT,2)=N_NODES+CONT+1  
                NODES_CONNECTIVITY(N_ELEM+CONT,3)=NODES_CONNECTIVITY(ELEMENT,LOCAL_NUMBER2) 
                NODES_CONNECTIVITY(N_ELEM+CONT,4)=NODES_CONNECTIVITY(ELEMENT,LOCAL_NUMBER1)                    
!
                COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,1)=N_COLLOCPOINTS+4*CONT-3
                COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,2)=N_COLLOCPOINTS+4*CONT-2 
                COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,3)=N_COLLOCPOINTS+4*CONT-1
                COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,4)=N_COLLOCPOINTS+4*CONT
!                
                NC(N_CRACK_POINTS+4*CONT-3,1)=N_COLLOCPOINTS+4*CONT-3
                NC(N_CRACK_POINTS+4*CONT-2,1)=N_COLLOCPOINTS+4*CONT-2
                NC(N_CRACK_POINTS+4*CONT-1,1)=N_COLLOCPOINTS+4*CONT-1
                NC(N_CRACK_POINTS+4*CONT,1)=N_COLLOCPOINTS+4*CONT
!
                IF(CONT.EQ.1)THEN
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,1))=1   
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,2))=1  
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,3))=0   
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,4))=0
                ELSE IF(CONT.EQ.DINT(N_CRACKTIPS/2.D0))THEN
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,1))=1   
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,2))=1  
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,3))=0   
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,4))=0               
                ELSE
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,1))=1   
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,2))=1  
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,3))=0   
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+CONT,4))=0                
                ENDIF 
!                
                NODES_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,1)=N_NODES+DINT((N_CRACKTIPS+2)/2.D0)+CONT  
                NODES_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,2)=NODES_CONNECTIVITY(ELEMENT,LOCAL_NUMBER1)     
                NODES_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,3)=NODES_CONNECTIVITY(ELEMENT,LOCAL_NUMBER2)
                NODES_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,4)=N_NODES+DINT((N_CRACKTIPS+2)/2.D0)+CONT+1                                                                                       
!
                COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,1)=N_COLLOCPOINTS+2*N_CRACKTIPELEMENTS+4*CONT-3
                COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,2)=N_COLLOCPOINTS+2*N_CRACKTIPELEMENTS+4*CONT-2
                COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,3)=N_COLLOCPOINTS+2*N_CRACKTIPELEMENTS+4*CONT-1
                COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,4)=N_COLLOCPOINTS+2*N_CRACKTIPELEMENTS+4*CONT
!                
                NC(N_CRACK_POINTS+4*CONT-3,2)=N_COLLOCPOINTS+2*N_CRACKTIPELEMENTS+4*CONT-3
                NC(N_CRACK_POINTS+4*CONT-2,2)=N_COLLOCPOINTS+2*N_CRACKTIPELEMENTS+4*CONT
                NC(N_CRACK_POINTS+4*CONT-1,2)=N_COLLOCPOINTS+2*N_CRACKTIPELEMENTS+4*CONT-1
                NC(N_CRACK_POINTS+4*CONT,2)=N_COLLOCPOINTS+2*N_CRACKTIPELEMENTS+4*CONT-2                
!
                IF(CONT.EQ.1)THEN
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,1))=1 
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,2))=0 
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,3))=0                                   
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,4))=1
                ELSE IF(CONT.EQ.DINT(N_CRACKTIPS/2.D0))THEN     
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,1))=1 
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,2))=0 
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,3))=0                                   
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,4))=1          
                ELSE
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,1))=1 
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,2))=0 
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,3))=0                                   
                    CRACK_TIP(COLLOCPOINTS_CONNECTIVITY(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,4))=1                
                ENDIF                                                               
!                                     
                ND(N_ELEM+CONT)=ND(ELEMENT)
                ND(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT)=ND(ELEMENT)
!
                NDI(2*I-1)=ND(ELEMENT) 
                NDI(2*I)=ND(ELEMENT)
                NDI(2*(I+1)-1)=ND(ELEMENT) 
                NDI(2*(I+1))=ND(ELEMENT)                                                                                                                 
!
                QSI(N_ELEM+CONT,1,1)=-0.8D0    
                QSI(N_ELEM+CONT,1,2)=-0.8D0    
                QSI(N_ELEM+CONT,2,1)=0.8D0 
                QSI(N_ELEM+CONT,2,2)=-0.8D0    
                QSI(N_ELEM+CONT,3,1)=0.8D0
                QSI(N_ELEM+CONT,3,2)=0.8D0    
                QSI(N_ELEM+CONT,4,1)=-0.8D0
                QSI(N_ELEM+CONT,4,2)=0.8D0    
!
                QSI(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,1,1)=-0.8D0    
                QSI(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,1,2)=-0.8D0    
                QSI(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,2,1)=0.8D0 
                QSI(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,2,2)=-0.8D0    
                QSI(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,3,1)=0.8D0
                QSI(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,3,2)=0.8D0    
                QSI(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,4,1)=-0.8D0
                QSI(N_ELEM+DINT(N_CRACKTIPS/2.D0)+CONT,4,2)=0.8D0                                   
!               
            ENDIF
        ENDDO 
!       
        DO I=N_COLLOCPOINTS+1,N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS
            DO J=N_ELEM+1,N_ELEM+N_CRACKTIPELEMENTS
                NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                ALLOCATE(VALUES(NEN,3))
                DO K=1,NEN
                    VALUES(K,1)=COORD_NODES(NODES_CONNECTIVITY(J,K),1)
	                VALUES(K,2)=COORD_NODES(NODES_CONNECTIVITY(J,K),2)
	                VALUES(K,3)=COORD_NODES(NODES_CONNECTIVITY(J,K),3)
	            ENDDO               
!                        
                DO K=1,NEN
                    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
                        CALL SHAPE_FUNCTIONS(1,QSI(J,K,1),QSI(J,K,2),4,VALUES,COEFFICIENTS,COORD_COLLOCPOINTS(I,1),COORD_COLLOCPOINTS(I,2),COORD_COLLOCPOINTS(I,3),PHI,DPHI,D2PHI)
                    ENDIF
                ENDDO            
                DEALLOCATE(VALUES)
            ENDDO
            IF(I.LE.(N_COLLOCPOINTS+DINT(N_CRACKTIPCOLLOCPOINTS/2.D0)))THEN
                DUAL_BEM(I)="S"
                EQ_TYPE(I)="S"
            ELSE
                DUAL_BEM(I)="H"
                EQ_TYPE(I)="H"            
            ENDIF
        ENDDO
!        
        DO I=1,N_INTPOINTS
            COORD_INTPOINTS(I,:)=COORD_COLLOCPOINTS(N_COLLOCPOINTS+I,:)
        ENDDO         
!------------------------------------------------------------------------------------------------------------------------------------------------------
!               ALLOCATING THE ADDITIONAL "L" AT THE MATRIX H AND G
!------------------------------------------------------------------------------------------------------------------------------------------------------ 
        ALLOCATE(HG_AUX(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS),UT_AUX(3*N_COLLOCPOINTS))
        HG_AUX=H
        UT_AUX=U
        DEALLOCATE(H,U)
        ALLOCATE(H(3*(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS),3*(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS)),U(3*(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS))) 
        H=0.D0
        U=0.D0
        DO I=1,3*N_COLLOCPOINTS
            U(I)=UT_AUX(I)
            DO J=1,3*N_COLLOCPOINTS
                H(I,J)=HG_AUX(I,J)
            ENDDO
        ENDDO
        HG_AUX=G
        UT_AUX=T
        DEALLOCATE(G,T)
        ALLOCATE(G(3*(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS),3*(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS)),T(3*(N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS))) 
        G=0.D0
        T=0.D0
        DO I=1,3*N_COLLOCPOINTS
            T(I)=UT_AUX(I)
            DO J=1,3*N_COLLOCPOINTS
                G(I,J)=HG_AUX(I,J)
            ENDDO
        ENDDO
!
        DO I=1,N_COLLOCPOINTS
            IF(CRACK_TIP(I).EQ.0)THEN
                CRACK_TIP(I)=0        
            ELSE IF(CRACK_TIP(I).EQ.1)THEN
                CRACK_TIP(I)=0
            ELSE IF(CRACK_TIP(I).EQ.2)THEN
                 CRACK_TIP(I)=1                  
            ELSE IF(CRACK_TIP(I).EQ.3)THEN 
                 CRACK_TIP(I)=1               
            ENDIF 
        ENDDO  
!                
        N_NODES=N_NODES+N_CRACKTIPNODES
        N_COLLOCPOINTS=N_COLLOCPOINTS+N_CRACKTIPCOLLOCPOINTS
	    N_ELEM=N_ELEM+N_CRACKTIPELEMENTS
	    N_CRACK_POINTS=N_CRACK_POINTS+2*N_CRACKTIPELEMENTS
!
        CALL CRACK_PROPAGATION_HG_UPDATE
!      	    
	ELSE
        N_CRACKTIPNODES=0
        N_CRACKTIPCOLLOCPOINTS=0
        N_CRACKTIPELEMENTS=0                                                                                 
    ENDIF                  
! 
    DEALLOCATE(COORD_CRACKTIP)
    DEALLOCATE(KI)
    DEALLOCATE(KII)
    DEALLOCATE(KIII)
    DEALLOCATE(COD)
    DEALLOCATE(CSD)
    DEALLOCATE(CTD)
    DEALLOCATE(THETAP)
    DEALLOCATE(PHIP)
    DEALLOCATE(Keq)
    DEALLOCATE(Gmax)
    DEALLOCATE(CRACK_INCREMENT)
    DEALLOCATE(VERTICAL_INCREMENT)
    DEALLOCATE(CRACKTIP_COLLOCPOINTS)
    IF(PROPAGATE) DEALLOCATE(COORD_INTPOINTS)
    IF(PROPAGATE) DEALLOCATE(NDI)
!           	    
    END SUBROUTINE CRACK_PROPAGATION

