	SUBROUTINE INTERNAL_DISPLACEMENTS
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
    INTEGER::I,J,K,L,II,JJ,KK,NEN,N,N_GAUSS
!
    REAL*8::AUX1,AUX2,Mu,Nu,PI,VJC(3),DVJC(3,2),JC,C1,C2,C3,X,Y,Z,DX,DY,DZ,R,DR(3),DRDN,KR(3,3),LENGTH(4),AVERAGE_LENGTH,DISTANCE,DAUX,ETA(3),VALUES[ALLOCATABLE](:,:),&
    PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),D2PHI[ALLOCATABLE](:,:,:),MATPHI[ALLOCATABLE](:,:),DGGAUX[ALLOCATABLE](:,:),DHHAUX[ALLOCATABLE](:,:),QSIW_GAUSS[ALLOCATABLE](:,:),&
    GG_LOCAL(3,3),HH_LOCAL(3,3),COEFFICIENTS[ALLOCATABLE](:,:),THETA,RHO,THETAI,THETAF,RHO_BOUND,RHO_BOUND_AUX1,RHO_BOUND_AUX2(3),QSI1,QSI2,QSI0(2),SOLUTION[ALLOCATABLE](:)
!    
    REAL*8::R_S(3,3),R_H(3,3),U_AUX(3),T_AUX(3)
! 
    REAL*8::SOMA_HH(3*N_INTPOINTS,3)
!
	ALLOCATE(HH(3*N_INTPOINTS,3*N_COLLOCPOINTS),GG(3*N_INTPOINTS,3*N_COLLOCPOINTS))
!
	ALLOCATE(U_INT(3*N_INTPOINTS))
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
	U_INT=0.D0
!
	PI=DACOS(-1.D0)
	KR(1,1)=1.D0
	KR(1,2)=0.D0
	KR(1,3)=0.D0
	KR(2,1)=0.D0
	KR(2,2)=1.D0
	KR(2,3)=0.D0
	KR(3,1)=0.D0
	KR(3,2)=0.D0
	KR(3,3)=1.D0		
!
	HH=0.D0
	GG=0.D0
!
	DO I=1,N_INTPOINTS	
	    Mu=EMP(NDI(I),1)/(2*(1+EMP(NDI(I),2)))
	    Nu=EMP(NDI(I),2)
	    C1=((1.D0)/(16.D0*PI*Mu*(1.D0-Nu)))
	    C2=((-1.D0)/(8.D0*PI*(1.D0-Nu)))
	    C3=Mu/(4*PI*(1.D0-Nu))
!
		DO J=1,N_ELEM
	        IF(CRACKED_ELEM(J).EQ."UNCRACKED")THEN			    
		        IF(NDI(I).EQ.ND(J)) THEN
		            NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
			        ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),MATPHI(3,3*NEN),DGGAUX(3,3*NEN),DHHAUX(3,3*NEN)) 
		            CALL SHAPE_FUNCTIONS_COEFFICIENTS(J,NEN,COEFFICIENTS)
!    				 
			        DO K=1,NEN
				        VALUES(K,1)=COORD_NODES(NODES_CONNECTIVITY(J,K),1)
				        VALUES(K,2)=COORD_NODES(NODES_CONNECTIVITY(J,K),2)
		                VALUES(K,3)=COORD_NODES(NODES_CONNECTIVITY(J,K),3)
			        ENDDO
    !
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    !               ADAPTATIVE GAUSS POINTS
    !------------------------------------------------------------------------------------------------------------------------------------------------------    
                    SELECT CASE(ELEM_TYPE(J))
                    CASE(3)
                        LENGTH(1)=DSQRT((VALUES(1,1)-VALUES(2,1))**2+(VALUES(1,2)-VALUES(2,2))**2+(VALUES(1,3)-VALUES(2,3))**2)
                        LENGTH(2)=DSQRT((VALUES(1,1)-VALUES(3,1))**2+(VALUES(1,2)-VALUES(3,2))**2+(VALUES(1,3)-VALUES(3,3))**2)
                        LENGTH(3)=DSQRT((VALUES(2,1)-VALUES(3,1))**2+(VALUES(2,2)-VALUES(3,2))**2+(VALUES(2,3)-VALUES(3,3))**2)
                        AVERAGE_LENGTH=(LENGTH(1)+LENGTH(2)+LENGTH(3))/3
                    CASE(4)
            	        LENGTH(1)=DSQRT((VALUES(1,1)-VALUES(2,1))**2+(VALUES(1,2)-VALUES(2,2))**2+(VALUES(1,3)-VALUES(2,3))**2)
                        LENGTH(2)=DSQRT((VALUES(2,1)-VALUES(3,1))**2+(VALUES(2,2)-VALUES(3,2))**2+(VALUES(2,3)-VALUES(3,3))**2)
                        LENGTH(3)=DSQRT((VALUES(3,1)-VALUES(4,1))**2+(VALUES(3,2)-VALUES(4,2))**2+(VALUES(3,3)-VALUES(4,3))**2)
                        LENGTH(4)=DSQRT((VALUES(4,1)-VALUES(1,1))**2+(VALUES(4,2)-VALUES(1,2))**2+(VALUES(4,3)-VALUES(1,3))**2)
                        AVERAGE_LENGTH=(LENGTH(1)+LENGTH(2)+LENGTH(3)+LENGTH(4))/4
                    ENDSELECT  
                    DISTANCE=0.D0
                    DO K=1,NEN
                        DAUX=DSQRT((COORD_INTPOINTS(I,1)-VALUES(K,1))**2+(COORD_INTPOINTS(I,2)-VALUES(K,2))**2+(COORD_INTPOINTS(I,3)-VALUES(K,3))**2)
                        IF(K.EQ.1)THEN
                            DISTANCE=DAUX
                        ELSE
                            IF(DAUX.LT.DISTANCE)THEN
                                DISTANCE=DAUX
                            ENDIF
                        ENDIF
                    ENDDO
                    N=DNINT(DSQRT(10/(DISTANCE/AVERAGE_LENGTH)))
                    IF(N.EQ.0)THEN
                        N=1
                    ENDIF	
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    !               ADAPTATIVE GAUSS POINTS
    !------------------------------------------------------------------------------------------------------------------------------------------------------      
                    SELECTCASE(ELEM_TYPE(J))
                    CASE(3)
                        N_GAUSS=15*N
                        ALLOCATE(QSIW_GAUSS(N_GAUSS,2))
                        CALL GAUSS_POINTS(N_GAUSS,QSIW_GAUSS)
                        THETAI=3*PI/4
                        THETAF=PI
                        QSI0(1)=1.D0
                        QSI0(2)=0.D0
    !                        
                        DO K=1,N_GAUSS
                            THETA=(THETAF-THETAI)/2*QSIW_GAUSS(K,1)+(THETAF+THETAI)/2
                        !***************************************************************************************************************************  
                        !       COMPUTING RHO_BOUND(THETA)     
                        !***************************************************************************************************************************             
                            RHO_BOUND_AUX2(1)=(-QSI0(1))/DCOS(THETA)
                            RHO_BOUND_AUX2(2)=(-QSI0(2))/DSIN(THETA)
                            RHO_BOUND_AUX2(3)=(1-QSI0(1)-QSI0(2))/(DCOS(THETA)+DSIN(THETA))                        
                            RHO_BOUND=DSQRT(8.D0)
                            DO L=1,3
                                IF(RHO_BOUND_AUX2(L).GT.0.D0)THEN
                                    RHO_BOUND_AUX1=RHO_BOUND_AUX2(L)
                                    IF(RHO_BOUND_AUX1.LE.RHO_BOUND)THEN
                                        RHO_BOUND=RHO_BOUND_AUX1
                                    ENDIF
                                ENDIF
                            ENDDO   
                        !***************************************************************************************************************************  
                            DO L=1,N_GAUSS
                                RHO=(RHO_BOUND/2.D0)*QSIW_GAUSS(L,1)+(RHO_BOUND/2.D0) 
                                DR=0.d0
                                DRDN=0.D0
                                ETA=0.D0
                                DPHI=0.D0
                                PHI=0.D0
                                MATPHI=0.D0
                                GG_LOCAL=0.D0
                                HH_LOCAL=0.D0
                                QSI1=RHO*DCOS(THETA)+QSI0(1)
                                QSI2=RHO*DSIN(THETA)+QSI0(2)          
    !
                                CALL SHAPE_FUNCTIONS(1,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
    !        
                                DX=X-COORD_INTPOINTS(I,1)
                                DY=Y-COORD_INTPOINTS(I,2)
                                DZ=Z-COORD_INTPOINTS(I,3)
                                R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                                DR(1)=DX/R
                                DR(2)=DY/R
                                DR(3)=DZ/R
    !
                                CALL NORMAL_OUTWARD(NEN,VALUES,QSI1,QSI2,ETA,VJC,DVJC,JC)
    !
                                DRDN=DR(1)*ETA(1)+DR(2)*ETA(2)+DR(3)*ETA(3)
    !                        
                                DO II=1,3
	                                DO JJ=1,3
		                                AUX1=0.D0
		                                AUX1=(3.D0-4.D0*Nu)*KR(II,JJ)
		                                AUX1=AUX1+DR(II)*DR(JJ)
		                                GG_LOCAL(II,JJ)=(1.D0/R)*C1*AUX1*RHO*QSIW_GAUSS(K,2)*QSIW_GAUSS(L,2)*JC*((THETAF-THETAI)/2)*(RHO_BOUND/2.D0)
    			                    
		                                AUX2=0.D0
		                                AUX2=(1.D0-2.D0*Nu)*KR(II,JJ)+3.D0*DR(II)*DR(JJ)
		                                AUX2=AUX2*DRDN
	                                    AUX2=AUX2-(1.D0-2.D0*Nu)*(DR(II)*ETA(JJ)-DR(JJ)*ETA(II))
		                                HH_LOCAL(II,JJ)=(1.D0/R**2)*C2*AUX2*RHO*QSIW_GAUSS(K,2)*QSIW_GAUSS(L,2)*JC*((THETAF-THETAI)/2)*(RHO_BOUND/2.D0)
	                                ENDDO
                                ENDDO                        
    !                
                                CALL SHAPE_FUNCTIONS(2,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
    !                                                    
                                DO II=1,NEN
	                                MATPHI(1,(3*II-2))=PHI(II)
	                                MATPHI(2,(3*II-1))=PHI(II)
	                                MATPHI(3,(3*II))=PHI(II)
                                ENDDO
    !		                    
			                    DGGAUX=0.D0
			                    DHHAUX=0.D0
    !
                                CALL DMRRRR (3,3,GG_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DGGAUX,3)
    !
                                CALL DMRRRR (3,3,HH_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHHAUX,3)
    !
			                    DO JJ=1,NEN
				                    KK=COLLOCPOINTS_CONNECTIVITY(J,JJ)
    !
                                    GG((3*I-2),(3*KK-2))=GG((3*I-2),(3*KK-2))+DGGAUX(1,(3*JJ-2))
			                        GG((3*I-2),3*KK-1)=GG((3*I-2),3*KK-1)+DGGAUX(1,3*JJ-1)
			                        GG((3*I-2),3*KK)=GG((3*I-2),3*KK)+DGGAUX(1,3*JJ)
			                        GG((3*I-1),(3*KK-2))=GG((3*I-1),(3*KK-2))+DGGAUX(2,(3*JJ-2))
			                        GG((3*I-1),3*KK-1)=GG((3*I-1),3*KK-1)+DGGAUX(2,3*JJ-1)
			                        GG((3*I-1),3*KK)=GG((3*I-1),3*KK)+DGGAUX(2,3*JJ)
			                        GG((3*I),(3*KK-2))=GG((3*I),(3*KK-2))+DGGAUX(3,(3*JJ-2))
			                        GG((3*I),3*KK-1)=GG((3*I),3*KK-1)+DGGAUX(3,3*JJ-1)
			                        GG((3*I),3*KK)=GG((3*I),3*KK)+DGGAUX(3,3*JJ)
    !
                                    HH((3*I-2),(3*KK-2))=HH((3*I-2),(3*KK-2))+DHHAUX(1,(3*JJ-2))
			                        HH((3*I-2),3*KK-1)=HH((3*I-2),3*KK-1)+DHHAUX(1,3*JJ-1)
			                        HH((3*I-2),3*KK)=HH((3*I-2),3*KK)+DHHAUX(1,3*JJ)
			                        HH((3*I-1),(3*KK-2))=HH((3*I-1),(3*KK-2))+DHHAUX(2,(3*JJ-2))
			                        HH((3*I-1),3*KK-1)=HH((3*I-1),3*KK-1)+DHHAUX(2,3*JJ-1)
			                        HH((3*I-1),3*KK)=HH((3*I-1),3*KK)+DHHAUX(2,3*JJ)
			                        HH((3*I),(3*KK-2))=HH((3*I),(3*KK-2))+DHHAUX(3,(3*JJ-2))
			                        HH((3*I),3*KK-1)=HH((3*I),3*KK-1)+DHHAUX(3,3*JJ-1)
			                        HH((3*I),3*KK)=HH((3*I),3*KK)+DHHAUX(3,3*JJ)						            				            					            				            				            					            				            				            
		   	                    ENDDO                        
                            ENDDO
                        ENDDO
                        DEALLOCATE(QSIW_GAUSS)            	                         
                    CASE(4)
                        N_GAUSS=15*N
                        ALLOCATE(QSIW_GAUSS(N_GAUSS,2))
                        CALL GAUSS_POINTS(N_GAUSS,QSIW_GAUSS)              
                        DO K=1,N_GAUSS
                            DO L=1,N_GAUSS
	                            DR=0.d0
	                            DRDN=0.D0
	                            ETA=0.D0
	                            DPHI=0.D0
	                            PHI=0.D0
	                            MATPHI=0.D0
	                            GG_LOCAL=0.D0
	                            HH_LOCAL=0.D0
    !
                                CALL SHAPE_FUNCTIONS(1,QSIW_GAUSS(K,1),QSIW_GAUSS(L,1),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
    !        
	                            DX=X-COORD_INTPOINTS(I,1)
	                            DY=Y-COORD_INTPOINTS(I,2)
	                            DZ=Z-COORD_INTPOINTS(I,3)
	                            R=DSQRT(DX*DX+DY*DY+DZ*DZ)
	                            DR(1)=DX/R
	                            DR(2)=DY/R
	                            DR(3)=DZ/R
    !
	                            CALL NORMAL_OUTWARD(NEN,VALUES,QSIW_GAUSS(K,1),QSIW_GAUSS(L,1),ETA,VJC,DVJC,JC)
    !
	                            DRDN=DR(1)*ETA(1)+DR(2)*ETA(2)+DR(3)*ETA(3)
    !
	                            DO II=1,3
		                            DO JJ=1,3
			                            AUX1=0.D0
			                            AUX1=(3.D0-4.D0*Nu)*KR(II,JJ)
			                            AUX1=AUX1+DR(II)*DR(JJ)
			                            GG_LOCAL(II,JJ)=(1.D0/R)*C1*AUX1*QSIW_GAUSS(K,2)*QSIW_GAUSS(L,2)*JC
    				                    
			                            AUX2=0.D0
			                            AUX2=(1.D0-2.D0*Nu)*KR(II,JJ)+3.D0*DR(II)*DR(JJ)
			                            AUX2=AUX2*DRDN
		                                AUX2=AUX2-(1.D0-2.D0*Nu)*(DR(II)*ETA(JJ)-DR(JJ)*ETA(II))
			                            HH_LOCAL(II,JJ)=(1.D0/R**2)*C2*AUX2*QSIW_GAUSS(K,2)*QSIW_GAUSS(L,2)*JC    				                    
		                            ENDDO
	                            ENDDO
    !
                                CALL SHAPE_FUNCTIONS(2,QSIW_GAUSS(K,1),QSIW_GAUSS(L,1),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
    !	
	                            DO II=1,NEN
		                            MATPHI(1,(3*II-2))=PHI(II)
		                            MATPHI(2,(3*II-1))=PHI(II)
		                            MATPHI(3,(3*II))=PHI(II)
	                            ENDDO
    !    		                    
				                DGGAUX=0.D0
				                DHHAUX=0.D0
    !
	                            CALL DMRRRR (3,3,GG_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DGGAUX,3)
    !
	                            CALL DMRRRR (3,3,HH_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHHAUX,3)
    		                    
				                DO JJ=1,NEN
					                KK=COLLOCPOINTS_CONNECTIVITY(J,JJ)					                
    !
	                                GG((3*I-2),(3*KK-2))=GG((3*I-2),(3*KK-2))+DGGAUX(1,(3*JJ-2))
				                    GG((3*I-2),3*KK-1)=GG((3*I-2),3*KK-1)+DGGAUX(1,3*JJ-1)
				                    GG((3*I-2),3*KK)=GG((3*I-2),3*KK)+DGGAUX(1,3*JJ)
				                    GG((3*I-1),(3*KK-2))=GG((3*I-1),(3*KK-2))+DGGAUX(2,(3*JJ-2))
				                    GG((3*I-1),3*KK-1)=GG((3*I-1),3*KK-1)+DGGAUX(2,3*JJ-1)
				                    GG((3*I-1),3*KK)=GG((3*I-1),3*KK)+DGGAUX(2,3*JJ)
				                    GG((3*I),(3*KK-2))=GG((3*I),(3*KK-2))+DGGAUX(3,(3*JJ-2))
				                    GG((3*I),3*KK-1)=GG((3*I),3*KK-1)+DGGAUX(3,3*JJ-1)
				                    GG((3*I),3*KK)=GG((3*I),3*KK)+DGGAUX(3,3*JJ)
    !
	                                HH((3*I-2),(3*KK-2))=HH((3*I-2),(3*KK-2))+DHHAUX(1,(3*JJ-2))
				                    HH((3*I-2),3*KK-1)=HH((3*I-2),3*KK-1)+DHHAUX(1,3*JJ-1)
				                    HH((3*I-2),3*KK)=HH((3*I-2),3*KK)+DHHAUX(1,3*JJ)
				                    HH((3*I-1),(3*KK-2))=HH((3*I-1),(3*KK-2))+DHHAUX(2,(3*JJ-2))
				                    HH((3*I-1),3*KK-1)=HH((3*I-1),3*KK-1)+DHHAUX(2,3*JJ-1)
				                    HH((3*I-1),3*KK)=HH((3*I-1),3*KK)+DHHAUX(2,3*JJ)
				                    HH((3*I),(3*KK-2))=HH((3*I),(3*KK-2))+DHHAUX(3,(3*JJ-2))
				                    HH((3*I),3*KK-1)=HH((3*I),3*KK-1)+DHHAUX(3,3*JJ-1)
				                    HH((3*I),3*KK)=HH((3*I),3*KK)+DHHAUX(3,3*JJ)				            				            					            				            				            					            				            				            
			   	                ENDDO		   	        
                            ENDDO
                        ENDDO
                        DEALLOCATE(QSIW_GAUSS)        
			        ENDSELECT
			        DEALLOCATE(VALUES,COEFFICIENTS,PHI,DPHI,MATPHI,DGGAUX,DHHAUX)
		        ENDIF
			ENDIF
		ENDDO
    ENDDO
!-----------------------------------------------------------------------------------------------------------
	SOMA_HH=0.D0
    DO I=1,3*N_INTPOINTS
        DO J=1,N_COLLOCPOINTS
            SOMA_HH(I,1)=SOMA_HH(I,1)+HH(I,3*J-2)
            SOMA_HH(I,2)=SOMA_HH(I,2)+HH(I,3*J-1)
            SOMA_HH(I,3)=SOMA_HH(I,3)+HH(I,3*J)
        ENDDO
        WRITE(*,*)SOMA_HH(I,1),SOMA_HH(I,2),SOMA_HH(I,3) 
    ENDDO
!    
!-----------------------------------------------------------------------------------------------------------
	ALLOCATE(SOLUTION(3*N_INTPOINTS)) 
	SOLUTION=0.D0
	CALL DMURRV(3*N_INTPOINTS,3*N_COLLOCPOINTS,GG,3*N_INTPOINTS,3*N_COLLOCPOINTS,T,1,3*N_INTPOINTS,SOLUTION)
	U_INT=SOLUTION
!
	SOLUTION=0.D0
	CALL DMURRV(3*N_INTPOINTS,3*N_COLLOCPOINTS,HH,3*N_INTPOINTS,3*N_COLLOCPOINTS,U,1,3*N_INTPOINTS,SOLUTION)
	U_INT=U_INT-SOLUTION
!	
	DEALLOCATE(SOLUTION,HH,GG)
!    
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