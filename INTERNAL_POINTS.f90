SUBROUTINE INTERNAL_POINTS
    
    ! SUBROUTINE TO CALCULATE DISPLACEMENTS AT INTERNAL POINTS GIVEN BY USER (INPUT)
    
	USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
	USE REMESHING
	USE PROPAGATION
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTERFACE 
        SUBROUTINE INTERNAL_POINTS_MATRIX(H_IC,G_IC,COORDS,PROP,DOMAIN)
            INTEGER::DOMAIN    
            REAL*8,DIMENSION(:,:)::H_IC,G_IC
            REAL*8::COORDS(3),PROP(5)
        END SUBROUTINE
    END INTERFACE
    INTERFACE 
        SUBROUTINE INTERNAL_POINTS_MATRIX_HP(H_IC,G_IC,COORDS,PROP,DOMAIN)
            INTEGER::DOMAIN    
            REAL*8,DIMENSION(:,:)::H_IC,G_IC
            REAL*8::COORDS(3),PROP(5)
        END SUBROUTINE
    END INTERFACE
    
    LOGICAL::INT_2D
    INTEGER::I,J,K,II,JJ,KK,NEN,N,N_GAUSS,P_INTEGRACAO,M
    
    REAL*8::VET_PROP(5),VET_COORD(3),SOLUTION[ALLOCATABLE](:),C1,MU,NU,JC,PI
    REAL*8,DIMENSION(:,:),ALLOCATABLE::G_IC,H_IC,HH_IC,GG_IC,HH_AUX,GG_AUX,HHH_AUX,GGG_AUX,SOMA_HH,DG_CF,DG_AUX_IF,G_IF,GG_IF,QSIW
        
    ! STARTING VARIABLES
    ALLOCATE(U_INT_INPUT(3*INPUT_INT_NODES),H_IC(3*INPUT_INT_NODES,3*N_COLLOCPOINTS),G_IC(3*INPUT_INT_NODES,3*N_COLLOCPOINTS),&
    SOMA_HH(3*INPUT_INT_NODES,3),HH_IC(9*INPUT_INT_NODES,3*N_COLLOCPOINTS),&
    GG_IC(9*INPUT_INT_NODES,3*N_COLLOCPOINTS),S_INT(9*INPUT_INT_NODES))
    PI=DACOS(-1.D0)
        
    IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
        ALLOCATE(G_IF(3*INPUT_INT_NODES,3*N_NOS_MEF),GG_IF(9*INPUT_INT_NODES,3*N_NOS_MEF))
        G_IF=0.0D0
        GG_IF=0.0D0
    ENDIF
        
            
    DO I=1,INPUT_INT_NODES !________________________________________________________________________________________________________
        
        VET_PROP(1)=EMP(NDI_INPUT(I),1)/(2*(1+EMP(NDI_INPUT(I),2)))             !Mu
	    VET_PROP(2)=EMP(NDI_INPUT(I),2)                                         !Nu
	    VET_PROP(3)=((1.D0)/(16.D0*PI*VET_PROP(1)*(1.D0-VET_PROP(2))))          !C1
	    VET_PROP(4)=((-1.D0)/(8.D0*PI*(1.D0-VET_PROP(2))))                      !C2
	    VET_PROP(5)=VET_PROP(1)/(4*PI*(1.D0-VET_PROP(2)))                       !C3
        
        VET_COORD(1)=COORDINT_NOD_INPUT(I,1)
        VET_COORD(2)=COORDINT_NOD_INPUT(I,2)
        VET_COORD(3)=COORDINT_NOD_INPUT(I,3)
        
        CALL INTERNAL_POINTS_MATRIX(H_IC(3*I-2:3*I,:),G_IC(3*I-2:3*I,:),VET_COORD,VET_PROP,NDI_INPUT(I))
        CALL INTERNAL_POINTS_MATRIX_HP(HH_IC(9*I-8:9*I,:),GG_IC(9*I-8:9*I,:),VET_COORD,VET_PROP,NDI_INPUT(I))
        
        !DO J=1,3*N_COLLOCPOINTS
        !    DO K=1,3
        !        H_IC(3*(I-1)+K,J) = HH_AUX(K,J)
        !        G_IC(3*(I-1)+K,J) = GG_AUX(K,J)
        !    ENDDO
        !    DO K=1,9
        !        HH_IC(3*(I-1)+K,J) = HHH_AUX(K,J)
        !        GG_IC(3*(I-1)+K,J) = GGG_AUX(K,J)
        !    ENDDO
        !ENDDO
        
        IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN ! ______ REINFORCED CASE ____________________
                         
            P_INTEGRACAO=INT(0.6*N_GAUSSPOINTS)
            IF (P_INTEGRACAO.LT.30) P_INTEGRACAO=30
            ALLOCATE(QSIW(P_INTEGRACAO,2))
            CALL GAUSS_POINTS(P_INTEGRACAO,QSIW)
            
            DO J=1,N_ELEMENTOS_MEF
                IF(NDI_INPUT(I).EQ.ND_MEF(J)) THEN
                    
                    ! Starting variables
                    ALLOCATE(DG_CF(3,3*(ORDEM_MEF(J)+1)))
                    Nu=EMP(NDI_INPUT(I),2)                               ! POISSON (Nu)
                    Mu=EMP(NDI_INPUT(I),1)/(2*(1+EMP(NDI_INPUT(I),2)))   ! G (Mu)

                    ! Integrating fiber (singular equation)
                    CALL FIBER_INT_DECISION(INT_2D,J,VET_COORD)
                    IF (INT_2D) THEN                                        
                        CALL FIBER_SING_S(J,Nu,Mu,P_INTEGRACAO,QSIW,VET_COORD,DG_CF)
                    ELSE                                                    
                        CALL FIBER_NONSING_S(J,Nu,Mu,VET_COORD,DG_CF)
                    ENDIF
                    
                    ! Contribuiting
                    DO K=1,ORDEM_MEF(J)+1
                        KK=CONECTI_MEF(J,K)
                        DO M=1,3
                            IF (DABS(DG_CF(M,3*K-2)).LT.TOLER) DG_CF(M,3*K-2)=0.0D0
                            IF (DABS(DG_CF(M,3*K-1)).LT.TOLER) DG_CF(M,3*K-1)=0.0D0
                            IF (DABS(DG_CF(M,3*K)).LT.TOLER) DG_CF(M,3*K)=0.0D0
                            
                            G_IF(3*(I-1)+M,3*KK-2)=G_IF(3*(I-1)+M,3*KK-2)+DG_CF(M,3*K-2)
                            G_IF(3*(I-1)+M,3*KK-1)=G_IF(3*(I-1)+M,3*KK-1)+DG_CF(M,3*K-1)
                            G_IF(3*(I-1)+M,3*KK)=G_IF(3*(I-1)+M,3*KK)+DG_CF(M,3*K)
                        ENDDO
                    ENDDO
                    DEALLOCATE(DG_CF)
                    
                    ! Integrating fiber (hipersingular equation)
                    ALLOCATE(DG_CF(9,3*(ORDEM_MEF(J)+1)))
                    IF (INT_2D) THEN                                        
                        CALL FIBER_SING_HP(J,Nu,Mu,P_INTEGRACAO,QSIW,VET_COORD,DG_CF)
                    ELSE                                                    
                        CALL FIBER_NONSING_HP(J,Nu,Mu,VET_COORD,DG_CF)
                    ENDIF
                    ! Contribuing in the global matrix for J element
                    DO K=1,ORDEM_MEF(J)+1
                        KK=CONECTI_MEF(J,K)
                        DO M=1,9
                            IF (DABS(DG_CF(M,3*K-2)).LT.TOLER) DG_CF(M,3*K-2)=0.0D0
                            IF (DABS(DG_CF(M,3*K-1)).LT.TOLER) DG_CF(M,3*K-1)=0.0D0
                            IF (DABS(DG_CF(M,3*K)).LT.TOLER) DG_CF(M,3*K)=0.0D0
                            
                            GG_IF(9*(I-1)+M,3*KK-2)=GG_IF(9*(I-1)+M,3*KK-2)+DG_CF(M,3*K-2)
                            GG_IF(9*(I-1)+M,3*KK-1)=GG_IF(9*(I-1)+M,3*KK-1)+DG_CF(M,3*K-1)
                            GG_IF(9*(I-1)+M,3*KK)=GG_IF(9*(I-1)+M,3*KK)+DG_CF(M,3*K)
                        ENDDO
                    ENDDO
                    DEALLOCATE(DG_CF)
                        
                ENDIF
            ENDDO
            
        ENDIF ! _____________ END OF REINFORECEMENT MATRIX _____________________________________
        
        DEALLOCATE(QSIW)   
    ENDDO ! END OF INTERNAL POINTS LOOP ____________________________________________________________________________________________
    
    IF (INPUT_INT_NODES.NE.0) THEN
    
        !-----------------------------------------------------------------------------------------------------------
	    SOMA_HH=0.D0
        WRITE(*,*)'MINHA SUBROTINA:'
        DO I=1,3*INPUT_INT_NODES
            DO J=1,N_COLLOCPOINTS
                SOMA_HH(I,1)=SOMA_HH(I,1)+H_IC(I,3*J-2)
                SOMA_HH(I,2)=SOMA_HH(I,2)+H_IC(I,3*J-1)
                SOMA_HH(I,3)=SOMA_HH(I,3)+H_IC(I,3*J)
            ENDDO
            WRITE(*,*)SOMA_HH(I,1),SOMA_HH(I,2),SOMA_HH(I,3) 
        ENDDO
        !    
        !-----------------------------------------------------------------------------------------------------------
    
        ! _____________________________________ APPLYING SOMIGLIANA EQUATION ___________________________________________________________
        ALLOCATE(SOLUTION(3*INPUT_INT_NODES)) 
	    SOLUTION=0.D0
	    CALL DMURRV(3*INPUT_INT_NODES,3*N_COLLOCPOINTS,G_IC,3*INPUT_INT_NODES,3*N_COLLOCPOINTS,T,1,3*INPUT_INT_NODES,SOLUTION)
	    U_INT_INPUT=SOLUTION    
    
        IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
            SOLUTION=0.0D0
            CALL DMURRV(3*INPUT_INT_NODES,3*N_NOS_MEF,G_IF,3*INPUT_INT_NODES,3*N_NOS_MEF,P_MEF,1,3*INPUT_INT_NODES,SOLUTION)
            U_INT_INPUT=U_INT_INPUT+SOLUTION
        ENDIF 
    !
	    SOLUTION=0.D0
	    CALL DMURRV(3*INPUT_INT_NODES,3*N_COLLOCPOINTS,H_IC,3*INPUT_INT_NODES,3*N_COLLOCPOINTS,U,1,3*INPUT_INT_NODES,SOLUTION)
	    U_INT_INPUT=U_INT_INPUT-SOLUTION

        
        !_____________________________________ APPLYING SOMIGLIANA HP EQUATION ___________________________________________________________
    
        DEALLOCATE(SOLUTION)
        ALLOCATE(SOLUTION(9*INPUT_INT_NODES)) 
	    SOLUTION=0.D0
	    CALL DMURRV(9*INPUT_INT_NODES,3*N_COLLOCPOINTS,GG_IC,9*INPUT_INT_NODES,3*N_COLLOCPOINTS,T,1,9*INPUT_INT_NODES,SOLUTION)
	    S_INT=SOLUTION
    
        IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') THEN
            SOLUTION=0.0D0
            CALL DMURRV(9*INPUT_INT_NODES,3*N_NOS_MEF,GG_IF,9*INPUT_INT_NODES,3*N_NOS_MEF,P_MEF,1,9*INPUT_INT_NODES,SOLUTION)
            S_INT=S_INT+SOLUTION
        ENDIF 
    
        SOLUTION=0.D0
	    CALL DMURRV(9*INPUT_INT_NODES,3*N_COLLOCPOINTS,HH_IC,9*INPUT_INT_NODES,3*N_COLLOCPOINTS,U,1,9*INPUT_INT_NODES,SOLUTION)
	    S_INT=S_INT-SOLUTION
    
    
	    DEALLOCATE(SOLUTION,H_IC,G_IC,HH_IC,GG_IC)
        IF (EXISTANCE_REINFORCEMENTS.EQ.'Y') DEALLOCATE(G_IF,GG_IF)
        ! ______________________________________________________________________________________________________________________________
    
        ! WRITING THE SOLUTION
        OPEN(12,file='Output_data/INTERNAL_POINTS.TXT',status='unknown')
        WRITE(12,*)''
        WRITE(12,*)'INTERNAL POINTS GIVEN BY USER RESULTS'
        WRITE(12,*)''
        WRITE(12,*)'POINT      X        Y         Z        DOMAIN'
        DO I=1,INPUT_INT_NODES
            WRITE(12,'(I4,3F11.3,I3)')I,COORDINT_NOD_INPUT(I,1),COORDINT_NOD_INPUT(I,2),COORDINT_NOD_INPUT(I,3),NDI_INPUT(I)
        ENDDO
        WRITE(12,*)''
        WRITE(12,*)'DISPLACEMENTS RESULTS'
        WRITE(12,*)'POINT       UX         UY          UZ'
        DO I=1,INPUT_INT_NODES
            WRITE(12,'(I4,3E15.6)')I,U_INT_INPUT(3*I-2),U_INT_INPUT(3*I-1),U_INT_INPUT(3*I)
        ENDDO
        WRITE(12,*)''
        WRITE(12,*)'STRESS RESULTS'
        WRITE(12,'(A)')' POINT    Sx-tot        Txy-tot       Txz-tot        Sy-tot           Tyz-tot          Sz-tot'
        DO I=1,INPUT_INT_NODES
            WRITE(12,'(I4,6E15.6)')I,S_INT(9*I-8),S_INT(9*I-7),S_INT(9*I-6),S_INT(9*I-4),S_INT(9*I-3),S_INT(9*I)
        ENDDO
        CLOSE(12)
    
    ENDIF
    
ENDSUBROUTINE
    
    
    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    

    
    
SUBROUTINE INTERNAL_POINTS_MATRIX(H_IC,G_IC,COORDS,PROP,DOMAIN)

    ! THIS SUBROUTINE CALCULATES THE H_IC AND G_IC MATRIXES LINES REFERENT TO A GIVEN INTERNAL POINT, WHICH WILL BE USED TO OBTAIN
    ! THE INTERNAL POINTS DISPLACEMENTS
    !
    ! H_IC: MATRIX (3,NELEM) FOR OUTPUT (INTEGRATION OF P*)
    ! G_IC: MATRIX (3,NELEM) FOR OUTPUT (INTEGRATION OF U*)
    ! COODS: X,Y,Z COORDINATES OF INTERNAL POINT IN QUESTION
    ! PROP: PROPERTIES OF INTERNAL NODE DOMAIN: (1)=MU, (2)=NU, (3)=C1, (4)=C2, (5)=C3
    ! DOMAIN: THE NUMBER OF THE 3D DOMAIN IN WHICH THE INTERNAL POINT IS POSITIONATED
    
    USE ISOPARAMETRIC_MESH
    USE SUB_REGIONS_INTERFACES
    
    IMPLICIT NONE
    
    INTEGER::DOMAIN
    INTEGER::I,II,J,JJ,K,KK,L,NEN,N,N_GAUSS
    
    !REAL*8,DIMENSION(3,3*N_COLLOCPOINTS)::H_IC,G_IC
    REAL*8,DIMENSION(:,:)::H_IC,G_IC
    
    REAL*8::COORDS(3),PROP(5),KR(3,3)  
    REAL*8::AUX,PI,VALUES[ALLOCATABLE](:,:),COEFFICIENTS[ALLOCATABLE](:,:),PHI[ALLOCATABLE](:),DPHI[ALLOCATABLE](:,:),&
    MATPHI[ALLOCATABLE](:,:),DGGAUX[ALLOCATABLE](:,:),DHHAUX[ALLOCATABLE](:,:),AVERAGE_LENGTH,DISTANCE,DAUX,LENGTH(4),&
    QSIW_GAUSS[ALLOCATABLE](:,:),THETAI,THETAF,QSI0(2),THETA,RHO_BOUND,RHO_BOUND_AUX1,RHO_BOUND_AUX2(3),RHO,DRDN,ETA(3),&
    QSI1,QSI2,D2PHI[ALLOCATABLE](:,:,:),DX,DY,DZ,R,DR(3),VJC(3),DVJC(3,2),JC,AUX1,AUX2,GG_LOCAL(3,3),HH_LOCAL(3,3),X,Y,Z
    
    PI=DACOS(-1.D0)
    H_IC=0.0D0
    G_IC=0.0D0
    KR=0.0D0
    KR(1,1)=1.0D0
    KR(2,2)=1.0D0
    KR(3,3)=1.0D0
    
    DO J=1,N_ELEM
        
        IF (DOMAIN .EQ. ND(J)) THEN
            
            NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
            ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),MATPHI(3,3*NEN),DGGAUX(3,3*NEN),DHHAUX(3,3*NEN))
            CALL SHAPE_FUNCTIONS_COEFFICIENTS(J,NEN,COEFFICIENTS)
            
            DO K=1,NEN
				VALUES(K,1)=COORD_NODES(NODES_CONNECTIVITY(J,K),1)
				VALUES(K,2)=COORD_NODES(NODES_CONNECTIVITY(J,K),2)
		        VALUES(K,3)=COORD_NODES(NODES_CONNECTIVITY(J,K),3)
            ENDDO
            
            ! -------------------- ADAPTATIVE GAUSS POINTS ------------------------------------------------------------------------
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
                DAUX=DSQRT((COORDS(1)-VALUES(K,1))**2+(COORDS(2)-VALUES(K,2))**2+(COORDS(3)-VALUES(K,3))**2)
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
            ! -------------- END OF ADAPTATIVE GAUSS POINTS -----------------------------------------------------------------------
            
            SELECTCASE(ELEM_TYPE(J))
                
            CASE(3) ! ______________ INTEGRATING TRIANGULAR ELEMENTS ______________________________________________________________
                
                N_GAUSS=15*N
                ALLOCATE(QSIW_GAUSS(N_GAUSS,2))
                CALL GAUSS_POINTS(N_GAUSS,QSIW_GAUSS)
                THETAI=3*PI/4
                THETAF=PI
                QSI0(1)=1.D0
                QSI0(2)=0.D0
                        
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

                        CALL SHAPE_FUNCTIONS(1,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
        
                        DX=X-COORDS(1)
                        DY=Y-COORDS(2)
                        DZ=Z-COORDS(3)
                        R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                        DR(1)=DX/R
                        DR(2)=DY/R
                        DR(3)=DZ/R

                        CALL NORMAL_OUTWARD(NEN,VALUES,QSI1,QSI2,ETA,VJC,DVJC,JC)

                        DRDN=DR(1)*ETA(1)+DR(2)*ETA(2)+DR(3)*ETA(3)
                        
                        DO II=1,3
	                        DO JJ=1,3
		                        AUX1=0.D0
		                        AUX1=(3.D0-4.D0*PROP(2))*KR(II,JJ)
		                        AUX1=AUX1+DR(II)*DR(JJ)
		                        GG_LOCAL(II,JJ)=(1.D0/R)*PROP(3)*AUX1*RHO*QSIW_GAUSS(K,2)*QSIW_GAUSS(L,2)*JC*((THETAF-THETAI)/2)*(RHO_BOUND/2.D0)
    			                    
		                        AUX2=0.D0
		                        AUX2=(1.D0-2.D0*PROP(2))*KR(II,JJ)+3.D0*DR(II)*DR(JJ)
		                        AUX2=AUX2*DRDN
	                            AUX2=AUX2-(1.D0-2.D0*PROP(2))*(DR(II)*ETA(JJ)-DR(JJ)*ETA(II))
		                        HH_LOCAL(II,JJ)=(1.D0/R**2)*PROP(4)*AUX2*RHO*QSIW_GAUSS(K,2)*QSIW_GAUSS(L,2)*JC*((THETAF-THETAI)/2)*(RHO_BOUND/2.D0)
	                        ENDDO
                        ENDDO                        
                
                        CALL SHAPE_FUNCTIONS(2,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
                                                   
                        DO II=1,NEN
	                        MATPHI(1,(3*II-2))=PHI(II)
	                        MATPHI(2,(3*II-1))=PHI(II)
	                        MATPHI(3,(3*II))=PHI(II)
                        ENDDO
	                    
			            DGGAUX=0.D0
			            DHHAUX=0.D0

                        CALL DMRRRR (3,3,GG_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DGGAUX,3)

                        CALL DMRRRR (3,3,HH_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHHAUX,3)
                       
			            DO JJ=1,NEN
				            KK=COLLOCPOINTS_CONNECTIVITY(J,JJ)
                            DO I=1,3
                                G_IC(I,(3*KK-2))=G_IC(I,(3*KK-2))+DGGAUX(I,(3*JJ-2))
				                G_IC(I,3*KK-1)=G_IC(I,3*KK-1)+DGGAUX(I,3*JJ-1)
				                G_IC(I,3*KK)=G_IC(I,3*KK)+DGGAUX(I,3*JJ)
                                
                                H_IC(I,(3*KK-2))=H_IC(I,(3*KK-2))+DHHAUX(I,(3*JJ-2))
				                H_IC(I,3*KK-1)=H_IC(I,3*KK-1)+DHHAUX(I,3*JJ-1)
				                H_IC(I,3*KK)=H_IC(I,3*KK)+DHHAUX(I,3*JJ)	
                            ENDDO						            				            					            				            				            					            				            				            
		   	            ENDDO                        
                    ENDDO
                ENDDO
                DEALLOCATE(QSIW_GAUSS)            	                   
                
            CASE(4) ! ______________ INTEGRATING QUADRILATERAL ELEMENTS ___________________________________________________________
                
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

                        CALL SHAPE_FUNCTIONS(1,QSIW_GAUSS(K,1),QSIW_GAUSS(L,1),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
       
	                    DX=X-COORDS(1)
	                    DY=Y-COORDS(2)
	                    DZ=Z-COORDS(3)
	                    R=DSQRT(DX*DX+DY*DY+DZ*DZ)
	                    DR(1)=DX/R
	                    DR(2)=DY/R
	                    DR(3)=DZ/R

	                    CALL NORMAL_OUTWARD(NEN,VALUES,QSIW_GAUSS(K,1),QSIW_GAUSS(L,1),ETA,VJC,DVJC,JC)

	                    DRDN=DR(1)*ETA(1)+DR(2)*ETA(2)+DR(3)*ETA(3)

	                    DO II=1,3
		                    DO JJ=1,3
			                    AUX1=0.D0
			                    AUX1=(3.D0-4.D0*PROP(2))*KR(II,JJ)
			                    AUX1=AUX1+DR(II)*DR(JJ)
			                    GG_LOCAL(II,JJ)=(1.D0/R)*PROP(3)*AUX1*QSIW_GAUSS(K,2)*QSIW_GAUSS(L,2)*JC
    				                    
			                    AUX2=0.D0
			                    AUX2=(1.D0-2.D0*PROP(2))*KR(II,JJ)+3.D0*DR(II)*DR(JJ)
			                    AUX2=AUX2*DRDN
		                        AUX2=AUX2-(1.D0-2.D0*PROP(2))*(DR(II)*ETA(JJ)-DR(JJ)*ETA(II))
			                    HH_LOCAL(II,JJ)=(1.D0/R**2)*PROP(4)*AUX2*QSIW_GAUSS(K,2)*QSIW_GAUSS(L,2)*JC    				                    
		                    ENDDO
	                    ENDDO

                        CALL SHAPE_FUNCTIONS(2,QSIW_GAUSS(K,1),QSIW_GAUSS(L,1),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
	
	                    DO II=1,NEN
		                    MATPHI(1,(3*II-2))=PHI(II)
		                    MATPHI(2,(3*II-1))=PHI(II)
		                    MATPHI(3,(3*II))=PHI(II)
	                    ENDDO
   		                    
				        DGGAUX=0.D0
				        DHHAUX=0.D0

	                    CALL DMRRRR (3,3,GG_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DGGAUX,3)

	                    CALL DMRRRR (3,3,HH_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHHAUX,3)
    		                    
				        DO JJ=1,NEN
					        KK=COLLOCPOINTS_CONNECTIVITY(J,JJ)					                
                            DO I=1,3
                                G_IC(I,(3*KK-2))=G_IC(I,(3*KK-2))+DGGAUX(I,(3*JJ-2))
				                G_IC(I,3*KK-1)=G_IC(I,3*KK-1)+DGGAUX(I,3*JJ-1)
				                G_IC(I,3*KK)=G_IC(I,3*KK)+DGGAUX(I,3*JJ)
                                
                                H_IC(I,(3*KK-2))=H_IC(I,(3*KK-2))+DHHAUX(I,(3*JJ-2))
				                H_IC(I,3*KK-1)=H_IC(I,3*KK-1)+DHHAUX(I,3*JJ-1)
				                H_IC(I,3*KK)=H_IC(I,3*KK)+DHHAUX(I,3*JJ)	
                            ENDDO			            				            					            				            				            					            				            				            
			   	        ENDDO		   	        
                    ENDDO
                ENDDO
                DEALLOCATE(QSIW_GAUSS)        
                
            END SELECT ! ______________ END OF BOUNDARY ELEMENTS INTEGRATION ______________________________________________________
            DEALLOCATE(VALUES,COEFFICIENTS,PHI,DPHI,MATPHI,DGGAUX,DHHAUX)
            
        ENDIF
        
    ENDDO
    
END SUBROUTINE
    
    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    
! ****************************************************************************************************************************************************    

    
SUBROUTINE INTERNAL_POINTS_MATRIX_HP(H_IC,G_IC,COORDS,PROP,DOMAIN)

    ! THIS SUBROUTINE CALCULATES THE H_IC AND G_IC MATRIXES LINES REFERENT TO A GIVEN INTERNAL POINT, WHICH WILL BE USED TO OBTAIN
    ! THE INTERNAL POINTS STRESSES
    !
    ! H_IC: MATRIX (9,3*N_COLLOCPOINTS) FOR OUTPUT (INTEGRATION OF S*)
    ! G_IC: MATRIX (9,3*N_COLLOCPOINTS) FOR OUTPUT (INTEGRATION OF D*)
    ! COODS: X,Y,Z COORDINATES OF INTERNAL POINT IN QUESTION
    ! PROP: PROPERTIES OF INTERNAL NODE DOMAIN: (1)=MU, (2)=NU, (3)=C1, (4)=C2, (5)=C3
    ! DOMAIN: THE NUMBER OF THE 3D DOMAIN IN WHICH THE INTERNAL POINT IS POSITIONATED
    
    USE ISOPARAMETRIC_MESH
    USE SUB_REGIONS_INTERFACES
    
    IMPLICIT NONE
    
    INTEGER::I,II,J,JJ,K,KK,L,NEN,N,M,N_GAUSS,J_EL,DOMAIN
    
    REAL*8,DIMENSION(:,:)::H_IC,G_IC
    
    REAL*8::COORDS(3),PROP(5),KR(3,3)  
    
    REAL*8::AUX,PI,AVERAGE_LENGTH,DISTANCE,DAUX,LENGTH(4),THETAI,THETAF,QSI0(2),THETA,RHO_BOUND,&
        RHO_BOUND_AUX1,RHO_BOUND_AUX2(3),RHO,DRDN,ETA(3),QSI1,QSI2,DX,DY,DZ,R,DR(3),VJC(3),DVJC(3,2),&
        JC,AUX1,AUX2,AUX3,GG_LOCAL(9,3),HH_LOCAL(9,3),X,Y,Z,C1,C2,NU
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::VALUES,COEFFICIENTS,DPHI,MATPHI,DGGAUX,DHHAUX,QSIW_GAUSS
    
    REAL*8::D2PHI[ALLOCATABLE](:,:,:),PHI[ALLOCATABLE](:)
    
    PI=DACOS(-1.D0)
    H_IC=0.0D0
    G_IC=0.0D0
    KR=0.0D0
    KR(1,1)=1.0D0
    KR(2,2)=1.0D0
    KR(3,3)=1.0D0
    C1=-PROP(4)
    C2=PROP(5)
    Nu=PROP(2)
    
    DO J_EL=1,N_ELEM
        
        IF (DOMAIN .EQ. ND(J_EL)) THEN
            
            NEN=ELEM_TYPE(J_EL)*ORDER_ELEM(J_EL)+(ELEM_TYPE(J_EL)-3)*(ORDER_ELEM(J_EL)-1)*POL_FAMILY(J_EL)
            ALLOCATE(VALUES(NEN,3),COEFFICIENTS(NEN,NEN),PHI(NEN),DPHI(NEN,2),MATPHI(3,3*NEN),DGGAUX(9,3*NEN),DHHAUX(9,3*NEN))
            CALL SHAPE_FUNCTIONS_COEFFICIENTS(J_EL,NEN,COEFFICIENTS)
            
            DO K=1,NEN
				VALUES(K,1)=COORD_NODES(NODES_CONNECTIVITY(J_EL,K),1)
				VALUES(K,2)=COORD_NODES(NODES_CONNECTIVITY(J_EL,K),2)
		        VALUES(K,3)=COORD_NODES(NODES_CONNECTIVITY(J_EL,K),3)
            ENDDO
            
            ! -------------------- ADAPTATIVE GAUSS POINTS ------------------------------------------------------------------------
            SELECT CASE(ELEM_TYPE(J_EL))
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
                DAUX=DSQRT((COORDS(1)-VALUES(K,1))**2+(COORDS(2)-VALUES(K,2))**2+(COORDS(3)-VALUES(K,3))**2)
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
            ! -------------- END OF ADAPTATIVE GAUSS POINTS -----------------------------------------------------------------------
            
            SELECTCASE(ELEM_TYPE(J_EL))
                
            CASE(3) ! ______________ INTEGRATING TRIANGULAR ELEMENTS ______________________________________________________________
                
                N_GAUSS=30*N
                ALLOCATE(QSIW_GAUSS(N_GAUSS,2))
                CALL GAUSS_POINTS(N_GAUSS,QSIW_GAUSS)
                THETAI=3*PI/4
                THETAF=PI
                QSI0(1)=1.D0
                QSI0(2)=0.D0
                        
                DO I=1,N_GAUSS
                    THETA=(THETAF-THETAI)/2*QSIW_GAUSS(K,1)+(THETAF+THETAI)/2
                !***************************************************************************************************************************  
                !       COMPUTING RHO_BOUND(THETA)     
                !***************************************************************************************************************************             
                    RHO_BOUND_AUX2(1)=(-QSI0(1))/DCOS(THETA)
                    RHO_BOUND_AUX2(2)=(-QSI0(2))/DSIN(THETA)
                    RHO_BOUND_AUX2(3)=(1-QSI0(1)-QSI0(2))/(DCOS(THETA)+DSIN(THETA))                        
                    RHO_BOUND=DSQRT(8.D0)
					DO K=1,3
						IF(RHO_BOUND_AUX2(K).GT.0.D0)THEN
							RHO_BOUND_AUX1=RHO_BOUND_AUX2(K)
							IF(RHO_BOUND_AUX1.LE.RHO_BOUND)THEN
								RHO_BOUND=RHO_BOUND_AUX1
							ENDIF
						ENDIF
					ENDDO  
                !***************************************************************************************************************************  
					DO J=1,N_GAUSS
						RHO=(RHO_BOUND/2.D0)*QSIW_GAUSS(J,1)+(RHO_BOUND/2.D0)            
						DR=0.d0
						DRDN=0.D0
						ETA=0.D0
						PHI=0.D0
						MATPHI=0.D0
						GG_LOCAL=0.D0
						HH_LOCAL=0.D0
						QSI1=RHO*DCOS(THETA)+QSI0(1)
						QSI2=RHO*DSIN(THETA)+QSI0(2)
			!
						CALL SHAPE_FUNCTIONS(1,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
        
                        DX=X-COORDS(1)
                        DY=Y-COORDS(2)
                        DZ=Z-COORDS(3)
                        R=DSQRT(DX*DX+DY*DY+DZ*DZ)
                        DR(1)=DX/R
                        DR(2)=DY/R
                        DR(3)=DZ/R

                        CALL NORMAL_OUTWARD(NEN,VALUES,QSI1,QSI2,ETA,VJC,DVJC,JC)

                        DRDN=DR(1)*ETA(1)+DR(2)*ETA(2)+DR(3)*ETA(3)
                        
						DO K=1,3
							DO L=1,3
								DO M=1,3
									AUX1=0.D0
									AUX1=KR(K,M)*DR(L)+KR(L,M)*DR(K)    
									AUX1=AUX1-KR(K,L)*DR(M)
									AUX1=AUX1*(1.D0-2.D0*Nu)
									AUX1=AUX1+3.D0*DR(K)*DR(L)*DR(M)
									GG_LOCAL((L+3*K-3),M)=(1.D0/R**2)*C1*AUX1*RHO*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*JC*((THETAF-THETAI)/2)*(RHO_BOUND/2.D0)
			!
									AUX2=0.D0
									AUX2=(1.D0-2.D0*Nu)*KR(K,L)*DR(M)
									AUX3=KR(K,M)*DR(L)+KR(L,M)*DR(K)
									AUX3=AUX3*Nu
									AUX2=AUX2+AUX3-5.D0*DR(K)*DR(L)*DR(M)
									AUX2=AUX2*3.D0*DRDN
			!
									AUX3=ETA(K)*DR(L)*DR(M)
									AUX3=AUX3+ETA(L)*DR(K)*DR(M)
									AUX3=AUX3*3.D0*Nu
									AUX2=AUX2+AUX3
			!
									AUX3=3.D0*ETA(M)*DR(K)*DR(L)
									AUX3=AUX3+ETA(L)*KR(K,M)
									AUX3=AUX3+ETA(K)*KR(L,M)
									AUX3=AUX3*(1.D0-2.D0*Nu)
									AUX2=AUX2+AUX3
			!
									AUX3=-(1.D0-4.D0*Nu)*ETA(M)*KR(K,L)
									AUX2=AUX2+AUX3
									HH_LOCAL((L+3*K-3),M)=(1.D0/R**3)*C2*AUX2*RHO*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*JC*((THETAF-THETAI)/2)*(RHO_BOUND/2.D0)
								ENDDO
							ENDDO
						ENDDO                        
                
                        CALL SHAPE_FUNCTIONS(2,QSI1,QSI2,NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
                                                   
                        DO II=1,NEN
	                        MATPHI(1,(3*II-2))=PHI(II)
	                        MATPHI(2,(3*II-1))=PHI(II)
	                        MATPHI(3,(3*II))=PHI(II)
                        ENDDO
	                    
                        ! CALL DMRRRR (3,3,GG_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DGGAUX,3)
						DGGAUX=MATMUL(GG_LOCAL,MATPHI)
						
                        ! CALL DMRRRR (3,3,HH_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHHAUX,3)
                        DHHAUX=MATMUL(HH_LOCAL,MATPHI)
					   
				        DO JJ=1,NEN
					        KK=COLLOCPOINTS_CONNECTIVITY(J_EL,JJ)
							DO M=1,9
								G_IC(M,(3*KK-2))=G_IC(M,(3*KK-2))+DGGAUX(M,(3*JJ-2))
								G_IC(M,3*KK-1)=G_IC(M,3*KK-1)+DGGAUX(M,3*JJ-1)
								G_IC(M,3*KK)=G_IC(M,3*KK)+DGGAUX(M,3*JJ)
								
								H_IC(M,(3*KK-2))=H_IC(M,(3*KK-2))+DHHAUX(M,(3*JJ-2))
								H_IC(M,3*KK-1)=H_IC(M,3*KK-1)+DHHAUX(M,3*JJ-1)
								H_IC(M,3*KK)=H_IC(M,3*KK)+DHHAUX(M,3*JJ)
							ENDDO				            				            					            				            				            					            				            				            
			   	        ENDDO                      
                    ENDDO
                ENDDO
                DEALLOCATE(QSIW_GAUSS)            	                   
                
            CASE(4) ! ______________ INTEGRATING QUADRILATERAL ELEMENTS ___________________________________________________________
                
                N_GAUSS=30*N
                ALLOCATE(QSIW_GAUSS(N_GAUSS,2))
                CALL GAUSS_POINTS(N_GAUSS,QSIW_GAUSS)    
                
                DO I=1,N_GAUSS
                    DO J=1,N_GAUSS
		                DR=0.d0
		                DRDN=0.D0
		                ETA=0.D0
		                PHI=0.D0
		                MATPHI=0.D0
		                GG_LOCAL=0.D0
		                HH_LOCAL=0.D0

                        CALL SHAPE_FUNCTIONS(1,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
       
	                    DX=X-COORDS(1)
	                    DY=Y-COORDS(2)
	                    DZ=Z-COORDS(3)
	                    R=DSQRT(DX*DX+DY*DY+DZ*DZ)
	                    DR(1)=DX/R
	                    DR(2)=DY/R
	                    DR(3)=DZ/R

	                    CALL NORMAL_OUTWARD(NEN,VALUES,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),ETA,VJC,DVJC,JC)

	                    DRDN=DR(1)*ETA(1)+DR(2)*ETA(2)+DR(3)*ETA(3)

		                DO K=1,3
			                DO L=1,3
			                    DO M=1,3
					                AUX1=0.D0
					                AUX1=KR(K,M)*DR(L)+KR(L,M)*DR(K)    
					                AUX1=AUX1-KR(K,L)*DR(M)
					                AUX1=AUX1*(1.D0-2.D0*Nu)
					                AUX1=AUX1+3.D0*DR(K)*DR(L)*DR(M)
					                GG_LOCAL((L+3*K-3),M)=(1.D0/R**2)*C1*AUX1*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*JC				        
            !
					                AUX2=0.D0
					                AUX2=(1.D0-2.D0*Nu)*KR(K,L)*DR(M)
					                AUX3=KR(K,M)*DR(L)+KR(L,M)*DR(K)
					                AUX3=AUX3*Nu
					                AUX2=AUX2+AUX3-5.D0*DR(K)*DR(L)*DR(M)
					                AUX2=AUX2*3.D0*DRDN
            !
					                AUX3=ETA(K)*DR(L)*DR(M)
					                AUX3=AUX3+ETA(L)*DR(K)*DR(M)
					                AUX3=AUX3*3.D0*Nu
					                AUX2=AUX2+AUX3
            !
					                AUX3=3.D0*ETA(M)*DR(K)*DR(L)
					                AUX3=AUX3+ETA(L)*KR(K,M)
					                AUX3=AUX3+ETA(K)*KR(L,M)
					                AUX3=AUX3*(1.D0-2.D0*Nu)
					                AUX2=AUX2+AUX3
            !
					                AUX3=-(1.D0-4.D0*Nu)*ETA(M)*KR(K,L)
					                AUX2=AUX2+AUX3
					                HH_LOCAL((L+3*K-3),M)=(1.D0/R**3)*C2*AUX2*QSIW_GAUSS(I,2)*QSIW_GAUSS(J,2)*JC
				                ENDDO
			                ENDDO
                        ENDDO
                        
            !           MULTIPLYING BY SHAPE FUNCTIONS ---------------------------------------------------------------------------------------------------	
                        CALL SHAPE_FUNCTIONS(2,QSIW_GAUSS(I,1),QSIW_GAUSS(J,1),NEN,VALUES,COEFFICIENTS,X,Y,Z,PHI,DPHI,D2PHI)
	
		                DO K=1,NEN
			                MATPHI(1,(3*K-2))=PHI(K)
			                MATPHI(2,(3*K-1))=PHI(K)
			                MATPHI(3,(3*K))=PHI(K)
		                ENDDO
   		                    
	                    !CALL DMRRRR (9,3,GG_LOCAL,9,3,(3*NEN),MATPHI,3,3,(3*NEN),DGGAUX,3)
						DGGAUX=MATMUL(GG_LOCAL,MATPHI)

	                    ! CALL DMRRRR (3,3,HH_LOCAL,3,3,(3*NEN),MATPHI,3,3,(3*NEN),DHHAUX,3)
						DHHAUX=MATMUL(HH_LOCAL,MATPHI)
    		                    
				        DO JJ=1,NEN
					        KK=COLLOCPOINTS_CONNECTIVITY(J_EL,JJ)
							DO M=1,9
								G_IC(M,(3*KK-2))=G_IC(M,(3*KK-2))+DGGAUX(M,(3*JJ-2))
								G_IC(M,3*KK-1)=G_IC(M,3*KK-1)+DGGAUX(M,3*JJ-1)
								G_IC(M,3*KK)=G_IC(M,3*KK)+DGGAUX(M,3*JJ)
								
								H_IC(M,(3*KK-2))=H_IC(M,(3*KK-2))+DHHAUX(M,(3*JJ-2))
								H_IC(M,3*KK-1)=H_IC(M,3*KK-1)+DHHAUX(M,3*JJ-1)
								H_IC(M,3*KK)=H_IC(M,3*KK)+DHHAUX(M,3*JJ)
							ENDDO				            				            					            				            				            					            				            				            
			   	        ENDDO		   	        
                    ENDDO
                ENDDO
                DEALLOCATE(QSIW_GAUSS)        
                
            END SELECT ! ______________ END OF BOUNDARY ELEMENTS INTEGRATION ______________________________________________________
            DEALLOCATE(VALUES,COEFFICIENTS,PHI,DPHI,MATPHI,DGGAUX,DHHAUX)
            
        ENDIF
        
    ENDDO
    
END SUBROUTINE