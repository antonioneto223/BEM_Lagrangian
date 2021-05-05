!SUBROUTINE LOCAL_MATRIXES_1DBEM(ELEM,KLOCAL,GLOCAL)
!    
!    !    
!    ! SUBROUTINE THAT CALCULATES KLOCAL AND GLOCAL FOR A GIVEN FIBER ELEMENT (ELEM) USING 1DBEM 
!    !
!    ! KLOCAL AND GLOCAL: MATRIXES CONTAINING THE RESULTS OF THIS SUBROUTINE
!    !
!    
!    USE REINFORCEMENTS
!    
!    IMPLICIT NONE
!    
!    INTEGER::ELEM,cont,K,J,II,JJ,P_INT_MEF,N_NO,NOF,SIDE
!    
!    REAL*8::KLOCAL(ORDEM_MEF(ELEM)+1,ORDEM_MEF(ELEM)+1),GLOCAL(ORDEM_MEF(ELEM)+1,ORDEM_MEF(ELEM)+1)
!    
!    REAL*8::GI_MEF[ALLOCATABLE](:),OME_MEF[ALLOCATABLE](:),QSIW_MEF[ALLOCATABLE](:,:),D_PHI_MEF[ALLOCATABLE](:),PHI_MEF[ALLOCATABLE](:),&
!        VALORES[ALLOCATABLE](:,:),TAN1,TAN2,TAN3,TAN_NORM,L_MEF,X,Y,Z,A_MATRIX[ALLOCATABLE](:,:),B_MATRIX[ALLOCATABLE](:,:),B_MATRIX_INV[ALLOCATABLE](:,:),&
!        C_MATRIX[ALLOCATABLE](:,:),AUX4,QSI_SOURCE,QSI_ELEMENT,JACOBIAN_SIDE(2),QSI_COORD(3),AUX1,AUX2,AUX3
!    
!    ! STARTING VARIABLES --------------
!        N_NO=ORDEM_MEF(ELEM)+1
!        ALLOCATE(D_PHI_MEF(N_NO),PHI_MEF(N_NO),VALORES(N_NO,3))
!	    PHI_MEF=0.D0
!	    D_PHI_MEF=0.D0
!	    L_MEF=0.D0
!	    VALORES=0.D0
!        ! ---------------------------------
!        
!        DO K=1,N_NO
!            VALORES(K,1)=COORDNOS_MEF(CONECTI_MEF(ELEM,K),1)
!            VALORES(K,2)=COORDNOS_MEF(CONECTI_MEF(ELEM,K),2)
!            VALORES(K,3)=COORDNOS_MEF(CONECTI_MEF(ELEM,K),3)
!        END DO
!        
!        ! IDENTIFYING DISCONTINUOUS NODES --------------
!        JJ=0
!        DO J=1,N_NO
!            JJ=JJ+KODENODUPLO_MEF(CONECTI_MEF(ELEM,J))
!        ENDDO
!        IF(JJ.EQ.0) THEN
!            NOF=0  !ONLY CONTINUOUS NODES ARE PRESENT
!        ELSE
!            NOF=1  !DISCONTINUOUS NODES ARE PRESENT
!        ENDIF
!    ! ----------------------------------------------
!    
!    
!    ! FINDING NUMBER OF QUADRATURE POINTS -----------------------------
!    P_INT_MEF=CEILING((ORDEM_MEF(ELEM)+1+4+1)/2.0D0)  !Soma 4 para considerar o Jacobiano na integral
!            
!    cont = 0
!    DO J=1,N_ELEMENTOS_MEF
!        DO K=1,ORDEM_MEF(ELEM)+1
!            IF ((CONECTI_MEF(J,K) .EQ. CONECTI_MEF(ELEM,1)) .OR.(CONECTI_MEF(J,K) .EQ. CONECTI_MEF(ELEM,ORDEM_MEF(ELEM)+1))) THEN
!                cont = cont + 1
!            ENDIF
!        ENDDO
!    ENDDO
!    cont = cont - 2  !subtraindo os encontros do proprio elemento
!    N_INT_P_MEF(ELEM) = P_INT_MEF
!    IF (cont .LT. 2) N_INT_P_MEF(ELEM) = 2*N_INT_P_MEF(ELEM)    
!    ! ------------------------------------------------------------------
!            
!    ! STARTING GAUSS INTEGRATION POINTS VARIABLES ----------------------
!    ALLOCATE(GI_MEF(P_INT_MEF),OME_MEF(P_INT_MEF),QSIW_MEF(P_INT_MEF,2))
!
!    GI_MEF=0.D0   !Array of length P_INT_MEF containing quadrature points.
!    OME_MEF=0.D0  !Array of length P_INT_MEF containing quadrature weights.
!    QSIW_MEF=0.D0 !Array with both quadrature points and weights
!
!    CALL DGQRUL (P_INT_MEF,1,0.D0,0.D0,0,GI_MEF,GI_MEF,OME_MEF)
!
!    QSIW_MEF(:,1)=GI_MEF(:)
!    QSIW_MEF(:,2)=OME_MEF(:)
!
!    DEALLOCATE(GI_MEF,OME_MEF)
!    ! ------------------------------------------------------------------
!    
!    DO J=1,P_INT_MEF ! PERFORMING NUMERICAL INTEGRATION FOR JACOBIAN 
!                
!        CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSIW_MEF(J,1),N_NO,ELEM,VALORES,X,Y,Z,D_PHI_MEF)
!                
!        ! COMPUTING THE JACOBIAN (TAN_NORM AND 2/L_MEF) --------
!        TAN1=0.D0
!	    TAN2=0.D0
!        TAN3=0.0D0
!	    DO K=1,N_NO
!		    TAN1=TAN1+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(ELEM,K),1) !TANGENT VECTORS
!		    TAN2=TAN2+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(ELEM,K),2) !TANGENT VECTORS
!            TAN3=TAN3+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(ELEM,K),3) !TANGENT VECTORS
!        ENDDO
!	    TAN_NORM=DSQRT(TAN1*TAN1+TAN2*TAN2+TAN3*TAN3) !TANGENT VECTOR NORM, JACOBIAN OF THE POINT
!	    L_MEF=L_MEF+TAN_NORM*QSIW_MEF(J,2)  !THE LENGHT OF THE ELEMENT IS THE SUM OF JACOBIANS, WEIGHTED BY THE GAUSS POINT WEIGHTS
!        ! ------------------------------------------------------
!                              
!    ENDDO 
!    
!    ! __________________ FINDING MATRIXES FOR 1DBEM _________________________
!    ALLOCATE(A_MATRIX(N_NO,N_NO),B_MATRIX(N_NO,N_NO),B_MATRIX_INV(N_NO,N_NO),C_MATRIX(N_NO,N_NO))
!            
!    B_MATRIX=0.0d0
!    A_MATRIX=0.0d0
!    A_MATRIX(:,1)=-1.0d0/2.0d0
!    A_MATRIX(:,N_NO)=-1.0d0/2.0d0
!    
!    DO K=1,N_NO
!        ! --------------------------------------------------------------------------  MUDOU
!        AUX1=COORDNOS_MEF(CONECTI_MEF(ELEM,1),1)-COORDNOS_MEF(CONECTI_MEF(ELEM,K),1)
!        AUX2=COORDNOS_MEF(CONECTI_MEF(ELEM,1),2)-COORDNOS_MEF(CONECTI_MEF(ELEM,K),2)
!        AUX3=COORDNOS_MEF(CONECTI_MEF(ELEM,1),3)-COORDNOS_MEF(CONECTI_MEF(ELEM,K),3)
!        AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
!        B_MATRIX(K,1)=-AUX4/(2.0d0*MP_MEF(ELEM,1)*MP_MEF(ELEM,2))
!        
!        AUX1=COORDNOS_MEF(CONECTI_MEF(ELEM,N_NO),1)-COORDNOS_MEF(CONECTI_MEF(ELEM,K),1)
!        AUX2=COORDNOS_MEF(CONECTI_MEF(ELEM,N_NO),2)-COORDNOS_MEF(CONECTI_MEF(ELEM,K),2)
!        AUX3=COORDNOS_MEF(CONECTI_MEF(ELEM,N_NO),3)-COORDNOS_MEF(CONECTI_MEF(ELEM,K),3)
!        AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
!        B_MATRIX(K,N_NO)=-AUX4/(2.0d0*MP_MEF(ELEM,1)*MP_MEF(ELEM,2))
!                
!        !AUX4=DABS(QSI_MEF(ELEM,1)-QSI_MEF(ELEM,K))*L_MEF/2.0D0               
!        !B_MATRIX(K,1)=-AUX4/(2.0d0*MP_MEF(ELEM,1)*MP_MEF(ELEM,2))
!        !        
!        !AUX4=DABS(QSI_MEF(ELEM,N_NO)-QSI_MEF(ELEM,K))*L_MEF/2.0D0
!        !B_MATRIX(K,N_NO)=-AUX4/(2.0d0*MP_MEF(ELEM,1)*MP_MEF(ELEM,2))
!        ! ---------------------------------------------------------------------------- MUDOU
!            
!        A_MATRIX(K,K)=A_MATRIX(K,K)+1.0d0                
!    ENDDO 
!    
!    DO K=2,N_NO-1
!        B_MATRIX(K,K)=L_MEF/(MP_MEF(ELEM,1)*MP_MEF(ELEM,2))
!    ENDDO
!            
!    CALL DLINRG(SIZE(B_MATRIX,2),B_MATRIX,SIZE(B_MATRIX,1),B_MATRIX_INV,SIZE(B_MATRIX_INV,1))
!            
!    KLOCAL=MATMUL(B_MATRIX_INV,A_MATRIX)   ! SITFFNESS MATRIX FROM 1DBEM (RESULT OF THIS SUBROUTINE)
!            
!    C_MATRIX=0.0D0
!    
!    DO K=1,N_NO
!        IF ((K.eq.1).or.(K.eq.N_NO)) THEN
!                    
!            DO J=1,P_INT_MEF !PERFORMIN NUMERICAL INTEGRATION FOR MATRIX G 
!                        
!                !CALL FUNCOES_DE_FORMA_MEC_MEF(2,QSIW_MEF(J,1),N_NO,ELEM,VALORES,X,Y,Z,PHI_MEF)
!                 CALL FUNCOES_DE_FORMA_MEC_MEF(4,QSIW_MEF(J,1),N_NO,ELEM,VALORES,X,Y,Z,PHI_MEF)
!    			
!                CALL FUNCOES_DE_FORMA_MEC_MEF(5,QSIW_MEF(J,1),N_NO,ELEM,VALORES,X,Y,Z,PHI_MEF)
!                QSI_COORD(1)=X
!                QSI_COORD(2)=Y   ! COORDENADAS DO PONTO DE INTEGRACAO
!                QSI_COORD(3)=Z
!                        
!                !CALL FUNCOES_DE_FORMA_MEC_MEF(3,qsi_element,N_NO,I,VALORES,X,Y,Z,D_PHI_MEF)
!                CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSIW_MEF(J,1),N_NO,ELEM,VALORES,X,Y,Z,D_PHI_MEF)  !VALIDADO EL RETO
!                        
!                DO II=1,N_NO
!    
!                    TAN1=0.0D0
!                    TAN2=0.0D0
!                    TAN3=0.0D0
!                    DO JJ=1,N_NO
!                        TAN1=TAN1+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(ELEM,JJ),1) !TANGENT VECTORS
!		                TAN2=TAN2+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(ELEM,JJ),2) !TANGENT VECTORS
!                        TAN3=TAN3+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(ELEM,JJ),3) !TANGENT VECTORS
!        	        ENDDO
!        	        TAN_NORM=DSQRT(TAN1*TAN1+TAN2*TAN2+TAN3*TAN3) !TANGENT VECTOR NORM, JACOBIAN OF THE POINT
!                            
!                    ! ---------------------------------------------------------------------------- MUDOU
!                    AUX1=QSI_COORD(1)-COORDNOS_MEF(CONECTI_MEF(ELEM,K),1)
!                    AUX2=QSI_COORD(2)-COORDNOS_MEF(CONECTI_MEF(ELEM,K),2)
!                    AUX3=QSI_COORD(3)-COORDNOS_MEF(CONECTI_MEF(ELEM,K),3)
!                    AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
!                            
!                    !AUX4=DABS(QSIW_MEF(J,1)-QSI_MEF(ELEM,K))*L_MEF/2.0D0
!                    ! ---------------------------------------------------------------------------- MUDOU
!                                
!                    C_MATRIX(K,II)=C_MATRIX(K,II)+AUX4*PHI_MEF(II)*TAN_NORM*QSIW_MEF(J,2)
!                ENDDO
!            ENDDO
!        ELSE
!!                Lenght_side=abs(coord(K)-coord(K))
!            qsi_source=QSI_MEF(ELEM,K)
!            Jacobian_side(1)=(qsi_source+1.0d0)/2.0d0
!            Jacobian_side(2)=(1.0d0-qsi_source)/2.0d0 
!            do SIDE=1,2             
!                do II=1,N_NO
!                    do J=1,P_INT_MEF
!                        if (SIDE.eq.1) then
!                            qsi_element=(QSIW_MEF(J,1)+1.0d0)*(1.0d0+qsi_source)/2.0d0-1.0d0
!                        elseif (SIDE.eq.2) then
!                            qsi_element=(QSIW_MEF(J,1)+1.0d0)*(1.0d0-qsi_source)/2.0d0+qsi_source
!                        end if
!                        !CALL FUNCOES_DE_FORMA_MEC_MEF(2,qsi_element,N_NO,ELEM,VALORES,X,Y,Z,PHI_MEF)
!                        CALL FUNCOES_DE_FORMA_MEC_MEF(4,qsi_element,N_NO,ELEM,VALORES,X,Y,Z,PHI_MEF)
!                        CALL FUNCOES_DE_FORMA_MEC_MEF(5,qsi_element,N_NO,ELEM,VALORES,X,Y,Z,PHI_MEF)
!                        QSI_COORD(1)=X
!                        QSI_COORD(2)=Y
!                        QSI_COORD(3)=Z
!!                            jacobian=get_qsi_jacobian(qsi_element,coord)
!                        CALL FUNCOES_DE_FORMA_MEC_MEF(6,qsi_element,N_NO,ELEM,VALORES,X,Y,Z,D_PHI_MEF)
!    
!        	            TAN1=0.D0
!	                    TAN2=0.D0
!                        TAN3=0.0D0
!	                    DO JJ=1,N_NO
!		                    TAN1=TAN1+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(ELEM,JJ),1) !TANGENT VECTORS
!		                    TAN2=TAN2+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(ELEM,JJ),2) !TANGENT VECTORS
!                            TAN3=TAN3+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(ELEM,JJ),3) !TANGENT VECTORS
!        	            ENDDO
!        	            TAN_NORM=DSQRT(TAN1*TAN1+TAN2*TAN2+TAN3*TAN3) !TANGENT VECTOR NORM, JACOBIAN OF THE POINT
!            	                
!                        ! ---------------------------------------------------------------------------- MUDOU
!                        AUX1=QSI_COORD(1)-COORDNOS_MEF(CONECTI_MEF(ELEM,K),1)
!                        AUX2=QSI_COORD(2)-COORDNOS_MEF(CONECTI_MEF(ELEM,K),2)
!                        AUX3=QSI_COORD(3)-COORDNOS_MEF(CONECTI_MEF(ELEM,K),3)
!                        AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
!                                
!                        !AUX4=DABS(qsi_element-QSI_MEF(ELEM,K))*L_MEF/2.0D0
!                        ! ---------------------------------------------------------------------------- MUDOU
!                
!                        C_MATRIX(K,II)=C_MATRIX(K,II)+AUX4*PHI_MEF(II)*Jacobian_side(SIDE)*QSIW_MEF(J,2)*TAN_NORM
!                    ENDDO
!                ENDDO
!            ENDDO
!        ENDIF
!    ENDDO
!            
!    
!    C_MATRIX=-1.0d0/(2.0d0*MP_MEF(ELEM,1)*MP_MEF(ELEM,2))*C_MATRIX
!            
!    GLOCAL=MATMUL(B_MATRIX_INV,C_MATRIX)   ! LUMPING MATRIX FOR 1DBEM (RESULT OF THIS SUBROUTINE)
!         
!    DEALLOCATE(A_MATRIX,B_MATRIX,B_MATRIX_INV,C_MATRIX,D_PHI_MEF,PHI_MEF,VALORES,QSIW_MEF)            
!    
!    ! ____________________________________________________ END OF 1DBEM MATRIXES _________________________________________________
!        
!
!    
!END SUBROUTINE
    

    
    
SUBROUTINE LOCAL_MATRIXES_1DBEM(ELEM,KLOCAL,GLOCAL)
    
    !    
    ! SUBROUTINE THAT CALCULATES KLOCAL AND GLOCAL FOR A GIVEN FIBER ELEMENT (ELEM) USING 1DBEM 
    !
    ! KLOCAL AND GLOCAL: MATRIXES CONTAINING THE RESULTS OF THIS SUBROUTINE
    !
    
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::ELEM,cont,K,J,II,JJ,P_INT_MEF,N_NO,NOF,SIDE
    
    REAL*8::KLOCAL(ORDEM_MEF(ELEM)+1,ORDEM_MEF(ELEM)+1),GLOCAL(ORDEM_MEF(ELEM)+1,ORDEM_MEF(ELEM)+1)
    
    REAL*8::GI_MEF[ALLOCATABLE](:),OME_MEF[ALLOCATABLE](:),QSIW_MEF[ALLOCATABLE](:,:),D_PHI_MEF[ALLOCATABLE](:),PHI_MEF[ALLOCATABLE](:),&
        VALORES[ALLOCATABLE](:,:),TAN1,TAN2,TAN3,TAN_NORM,L_MEF,X,Y,Z,A_MATRIX[ALLOCATABLE](:,:),B_MATRIX[ALLOCATABLE](:,:),B_MATRIX_INV[ALLOCATABLE](:,:),&
        C_MATRIX[ALLOCATABLE](:,:),AUX4,QSI_SOURCE,QSI_ELEMENT,JACOBIAN_SIDE(2),QSI_COORD(3)
    
    ! STARTING VARIABLES --------------
        N_NO=ORDEM_MEF(ELEM)+1
        ALLOCATE(D_PHI_MEF(N_NO),PHI_MEF(N_NO),VALORES(N_NO,3))
	    PHI_MEF=0.D0
	    D_PHI_MEF=0.D0
	    L_MEF=0.D0
	    VALORES=0.D0
        ! ---------------------------------
        
        DO K=1,N_NO
            VALORES(K,1)=COORDNOS_MEF(CONECTI_MEF(ELEM,K),1)
            VALORES(K,2)=COORDNOS_MEF(CONECTI_MEF(ELEM,K),2)
            VALORES(K,3)=COORDNOS_MEF(CONECTI_MEF(ELEM,K),3)
        END DO
        
        ! IDENTIFYING DISCONTINUOUS NODES --------------
        JJ=0
        DO J=1,N_NO
            JJ=JJ+KODENODUPLO_MEF(CONECTI_MEF(ELEM,J))
        ENDDO
        IF(JJ.EQ.0) THEN
            NOF=0  !ONLY CONTINUOUS NODES ARE PRESENT
        ELSE
            NOF=1  !DISCONTINUOUS NODES ARE PRESENT
        ENDIF
    ! ----------------------------------------------
    
    
    ! FINDING NUMBER OF QUADRATURE POINTS -----------------------------
    P_INT_MEF=CEILING((ORDEM_MEF(ELEM)+1+4+1)/2.0D0)  !Soma 4 para considerar o Jacobiano na integral
            
    cont = 0
    DO J=1,N_ELEMENTOS_MEF
        DO K=1,ORDEM_MEF(ELEM)+1
            IF ((CONECTI_MEF(J,K) .EQ. CONECTI_MEF(ELEM,1)) .OR.(CONECTI_MEF(J,K) .EQ. CONECTI_MEF(ELEM,ORDEM_MEF(ELEM)+1))) THEN
                cont = cont + 1
            ENDIF
        ENDDO
    ENDDO
    cont = cont - 2  !subtraindo os encontros do proprio elemento
    N_INT_P_MEF(ELEM) = P_INT_MEF
    IF (cont .LT. 2) N_INT_P_MEF(ELEM) = 2*N_INT_P_MEF(ELEM)    
    ! ------------------------------------------------------------------
            
    ! STARTING GAUSS INTEGRATION POINTS VARIABLES ----------------------
    ALLOCATE(GI_MEF(P_INT_MEF),OME_MEF(P_INT_MEF),QSIW_MEF(P_INT_MEF,2))

    GI_MEF=0.D0   !Array of length P_INT_MEF containing quadrature points.
    OME_MEF=0.D0  !Array of length P_INT_MEF containing quadrature weights.
    QSIW_MEF=0.D0 !Array with both quadrature points and weights

    CALL DGQRUL (P_INT_MEF,1,0.D0,0.D0,0,GI_MEF,GI_MEF,OME_MEF)

    QSIW_MEF(:,1)=GI_MEF(:)
    QSIW_MEF(:,2)=OME_MEF(:)

    DEALLOCATE(GI_MEF,OME_MEF)
    ! ------------------------------------------------------------------
    
    DO J=1,P_INT_MEF ! PERFORMING NUMERICAL INTEGRATION FOR JACOBIAN 
                
        CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSIW_MEF(J,1),N_NO,ELEM,VALORES,X,Y,Z,D_PHI_MEF)
                
        ! COMPUTING THE JACOBIAN (TAN_NORM AND 2/L_MEF) --------
        TAN1=0.D0
	    TAN2=0.D0
        TAN3=0.0D0
	    DO K=1,N_NO
		    TAN1=TAN1+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(ELEM,K),1) !TANGENT VECTORS
		    TAN2=TAN2+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(ELEM,K),2) !TANGENT VECTORS
            TAN3=TAN3+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(ELEM,K),3) !TANGENT VECTORS
        ENDDO
	    TAN_NORM=DSQRT(TAN1*TAN1+TAN2*TAN2+TAN3*TAN3) !TANGENT VECTOR NORM, JACOBIAN OF THE POINT
	    L_MEF=L_MEF+TAN_NORM*QSIW_MEF(J,2)  !THE LENGHT OF THE ELEMENT IS THE SUM OF JACOBIANS, WEIGHTED BY THE GAUSS POINT WEIGHTS
        ! ------------------------------------------------------
                              
    ENDDO 
    
    ! __________________ FINDING MATRIXES FOR 1DBEM _________________________
    ALLOCATE(A_MATRIX(N_NO,N_NO),B_MATRIX(N_NO,N_NO),B_MATRIX_INV(N_NO,N_NO),C_MATRIX(N_NO,N_NO))
            
    B_MATRIX=0.0d0
    A_MATRIX=0.0d0
    A_MATRIX(:,1)=-1.0d0/2.0d0
    A_MATRIX(:,N_NO)=-1.0d0/2.0d0
    
    DO K=1,N_NO
        ! --------------------------------------------------------------------------  MUDOU
        !AUX1=COORDNOS_MEF(CONECTI_MEF(I,1),1)-COORDNOS_MEF(CONECTI_MEF(I,K),1)
        !AUX2=COORDNOS_MEF(CONECTI_MEF(I,1),2)-COORDNOS_MEF(CONECTI_MEF(I,K),2)
        !AUX3=COORDNOS_MEF(CONECTI_MEF(I,1),3)-COORDNOS_MEF(CONECTI_MEF(I,K),3)
        !AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
        !B_MATRIX(K,1)=-AUX4/(2.0d0*MP_MEF(I,1)*MP_MEF(I,2))
        !
        !AUX1=COORDNOS_MEF(CONECTI_MEF(I,N_NO),1)-COORDNOS_MEF(CONECTI_MEF(I,K),1)
        !AUX2=COORDNOS_MEF(CONECTI_MEF(I,N_NO),2)-COORDNOS_MEF(CONECTI_MEF(I,K),2)
        !AUX3=COORDNOS_MEF(CONECTI_MEF(I,N_NO),3)-COORDNOS_MEF(CONECTI_MEF(I,K),3)
        !AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
        !B_MATRIX(K,N_NO)=-AUX4/(2.0d0*MP_MEF(I,1)*MP_MEF(I,2))
                
        AUX4=DABS(QSI_MEF(ELEM,1)-QSI_MEF(ELEM,K))*L_MEF/2.0D0               
        B_MATRIX(K,1)=-AUX4/(2.0d0*MP_MEF(ELEM,1)*MP_MEF(ELEM,2))
                
        AUX4=DABS(QSI_MEF(ELEM,N_NO)-QSI_MEF(ELEM,K))*L_MEF/2.0D0
        B_MATRIX(K,N_NO)=-AUX4/(2.0d0*MP_MEF(ELEM,1)*MP_MEF(ELEM,2))
        ! ---------------------------------------------------------------------------- MUDOU
            
        A_MATRIX(K,K)=A_MATRIX(K,K)+1.0d0                
    ENDDO 
    
    DO K=2,N_NO-1
        B_MATRIX(K,K)=L_MEF/(MP_MEF(ELEM,1)*MP_MEF(ELEM,2))
    ENDDO
            
    CALL DLINRG(SIZE(B_MATRIX,2),B_MATRIX,SIZE(B_MATRIX,1),B_MATRIX_INV,SIZE(B_MATRIX_INV,1))
            
    KLOCAL=MATMUL(B_MATRIX_INV,A_MATRIX)   ! SITFFNESS MATRIX FROM 1DBEM (RESULT OF THIS SUBROUTINE)
            
    C_MATRIX=0.0D0
    
    DO K=1,N_NO
        IF ((K.eq.1).or.(K.eq.N_NO)) THEN
                    
            DO J=1,P_INT_MEF !PERFORMIN NUMERICAL INTEGRATION FOR MATRIX G 
                        
                CALL FUNCOES_DE_FORMA_MEC_MEF(2,QSIW_MEF(J,1),N_NO,ELEM,VALORES,X,Y,Z,PHI_MEF)
    			
                CALL FUNCOES_DE_FORMA_MEC_MEF(5,QSIW_MEF(J,1),N_NO,ELEM,VALORES,X,Y,Z,PHI_MEF)
                QSI_COORD(1)=X
                QSI_COORD(2)=Y   ! COORDENADAS DO PONTO DE INTEGRACAO
                QSI_COORD(3)=Z
                        
                !CALL FUNCOES_DE_FORMA_MEC_MEF(3,qsi_element,N_NO,I,VALORES,X,Y,Z,D_PHI_MEF)
                CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSIW_MEF(J,1),N_NO,ELEM,VALORES,X,Y,Z,D_PHI_MEF)  !VALIDADO EL RETO
                        
                DO II=1,N_NO
    
                    TAN1=0.0D0
                    TAN2=0.0D0
                    TAN3=0.0D0
                    DO JJ=1,N_NO
                        TAN1=TAN1+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(ELEM,JJ),1) !TANGENT VECTORS
		                TAN2=TAN2+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(ELEM,JJ),2) !TANGENT VECTORS
                        TAN3=TAN3+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(ELEM,JJ),3) !TANGENT VECTORS
        	        ENDDO
        	        TAN_NORM=DSQRT(TAN1*TAN1+TAN2*TAN2+TAN3*TAN3) !TANGENT VECTOR NORM, JACOBIAN OF THE POINT
                            
                    ! ---------------------------------------------------------------------------- MUDOU
                    !AUX1=QSI_COORD(1)-COORDNOS_MEF(CONECTI_MEF(I,K),1)
                    !AUX2=QSI_COORD(2)-COORDNOS_MEF(CONECTI_MEF(I,K),2)
                    !AUX3=QSI_COORD(3)-COORDNOS_MEF(CONECTI_MEF(I,K),3)
                    !AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
                            
                    AUX4=DABS(QSIW_MEF(J,1)-QSI_MEF(ELEM,K))*L_MEF/2.0D0
                    ! ---------------------------------------------------------------------------- MUDOU
                                
                    C_MATRIX(K,II)=C_MATRIX(K,II)+AUX4*PHI_MEF(II)*TAN_NORM*QSIW_MEF(J,2)
                ENDDO
            ENDDO
        ELSE
!           Lenght_side=abs(coord(K)-coord(K))
            qsi_source=QSI_MEF(ELEM,K)
            Jacobian_side(1)=(qsi_source+1.0d0)/2.0d0
            Jacobian_side(2)=(1.0d0-qsi_source)/2.0d0 
            do SIDE=1,2             
                do II=1,N_NO
                    do J=1,P_INT_MEF
                        if (SIDE.eq.1) then
                            qsi_element=(QSIW_MEF(J,1)+1.0d0)*(1.0d0+qsi_source)/2.0d0-1.0d0
                        elseif (SIDE.eq.2) then
                            qsi_element=(QSIW_MEF(J,1)+1.0d0)*(1.0d0-qsi_source)/2.0d0+qsi_source
                        end if
                        CALL FUNCOES_DE_FORMA_MEC_MEF(2,qsi_element,N_NO,ELEM,VALORES,X,Y,Z,PHI_MEF)
                        CALL FUNCOES_DE_FORMA_MEC_MEF(5,qsi_element,N_NO,ELEM,VALORES,X,Y,Z,PHI_MEF)
                        QSI_COORD(1)=X
                        QSI_COORD(2)=Y
                        QSI_COORD(3)=Z
!                            jacobian=get_qsi_jacobian(qsi_element,coord)
                        CALL FUNCOES_DE_FORMA_MEC_MEF(6,qsi_element,N_NO,ELEM,VALORES,X,Y,Z,D_PHI_MEF)
    
        	            TAN1=0.D0
	                    TAN2=0.D0
                        TAN3=0.0D0
	                    DO JJ=1,N_NO
		                    TAN1=TAN1+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(ELEM,JJ),1) !TANGENT VECTORS
		                    TAN2=TAN2+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(ELEM,JJ),2) !TANGENT VECTORS
                            TAN3=TAN3+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(ELEM,JJ),3) !TANGENT VECTORS
        	            ENDDO
        	            TAN_NORM=DSQRT(TAN1*TAN1+TAN2*TAN2+TAN3*TAN3) !TANGENT VECTOR NORM, JACOBIAN OF THE POINT
            	                
                        ! ---------------------------------------------------------------------------- MUDOU
                        !AUX1=QSI_COORD(1)-COORDNOS_MEF(CONECTI_MEF(I,K),1)
                        !AUX2=QSI_COORD(2)-COORDNOS_MEF(CONECTI_MEF(I,K),2)
                        !AUX3=QSI_COORD(3)-COORDNOS_MEF(CONECTI_MEF(I,K),3)
                        !AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
                                
                        AUX4=DABS(qsi_element-QSI_MEF(ELEM,K))*L_MEF/2.0D0
                        ! ---------------------------------------------------------------------------- MUDOU
                
                        C_MATRIX(K,II)=C_MATRIX(K,II)+AUX4*PHI_MEF(II)*Jacobian_side(SIDE)*QSIW_MEF(J,2)*TAN_NORM
                    ENDDO
                ENDDO
            ENDDO
        ENDIF
    ENDDO
            
    
    C_MATRIX=-1.0d0/(2.0d0*MP_MEF(ELEM,1)*MP_MEF(ELEM,2))*C_MATRIX
            
    GLOCAL=MATMUL(B_MATRIX_INV,C_MATRIX)   ! LUMPING MATRIX FOR 1DBEM (RESULT OF THIS SUBROUTINE)
         
    DEALLOCATE(A_MATRIX,B_MATRIX,B_MATRIX_INV,C_MATRIX,D_PHI_MEF,PHI_MEF,VALORES,QSIW_MEF)            
    
    ! ____________________________________________________ END OF 1DBEM MATRIXES _________________________________________________
        

    
END SUBROUTINE
    
    
    
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  
    
    
     
    
    
SUBROUTINE LOCAL_MATRIXES_FEM(ELEM,KLOCAL,GLOCAL)
    
    !    
    ! SUBROUTINE THAT CALCULATES KLOCAL AND GLOCAL FOR A GIVEN FIBER ELEMENT (ELEM) USING FEM METHODO 
    !
    ! KLOCAL AND GLOCAL: MATRIXES CONTAINING THE RESULTS OF THIS SUBROUTINE
    !

    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::I,ELEM,J,K,JJ,NOF,N_NO,P_INT_MEF,cont
    
    REAL*8::KLOCAL(ORDEM_MEF(ELEM)+1,ORDEM_MEF(ELEM)+1),GLOCAL(ORDEM_MEF(ELEM)+1,ORDEM_MEF(ELEM)+1)
    
    REAL*8::B[ALLOCATABLE](:),BT[ALLOCATABLE](:),KLOCALAUX[ALLOCATABLE](:,:),GLOCALAUX[ALLOCATABLE](:,:),VALORES[ALLOCATABLE](:,:),&
        PHI_MEF[ALLOCATABLE](:),D_PHI_MEF[ALLOCATABLE](:),L_MEF,GI_MEF[ALLOCATABLE](:),OME_MEF[ALLOCATABLE](:),QSIW_MEF[ALLOCATABLE](:,:),&
        TAN1,TAN2,TAN3,TAN_NORM,X,Y,Z

    ! GIVEN THE ELEMENT NUMBER:
    I=ELEM
    
    
    ! STARTING VARIABLES --------------
    N_NO=ORDEM_MEF(I)+1
    ALLOCATE(B(N_NO),BT(N_NO),KLOCALAUX(N_NO,N_NO),D_PHI_MEF(N_NO),PHI_MEF(N_NO),GLOCALAUX(N_NO,N_NO),VALORES(N_NO,3))
    B=0.D0
	BT=0.D0
	PHI_MEF=0.D0
	D_PHI_MEF=0.D0
	KLOCAL=0.D0
	GLOCAL=0.D0
	L_MEF=0.D0
	VALORES=0.D0
    ! ---------------------------------
        
    DO K=1,N_NO
        VALORES(K,1)=COORDNOS_MEF(CONECTI_MEF(I,K),1)
        VALORES(K,2)=COORDNOS_MEF(CONECTI_MEF(I,K),2)
        VALORES(K,3)=COORDNOS_MEF(CONECTI_MEF(I,K),3)
    END DO
        
    ! IDENTIFYING DISCONTINUOUS NODES --------------
    JJ=0
    DO J=1,N_NO
        JJ=JJ+KODENODUPLO_MEF(CONECTI_MEF(I,J))
    ENDDO
    IF(JJ.EQ.0) THEN
        NOF=0  !ONLY CONTINUOUS NODES ARE PRESENT
    ELSE
        NOF=1  !DISCONTINUOUS NODES ARE PRESENT
    ENDIF
    ! ----------------------------------------------
    
    ! FINDING NUMBER OF QUADRATURE POINTS -------
    P_INT_MEF=ORDEM_MEF(I)-1
    P_INT_MEF=P_INT_MEF*P_INT_MEF
    IF(P_INT_MEF.LT.2) THEN
        P_INT_MEF=2
    ELSE
        J=(P_INT_MEF+1)/2
        P_INT_MEF=(J+MOD((P_INT_MEF+1),2))    
    ENDIF
            
    CONT = 0
    DO J=1,N_ELEMENTOS_MEF
        DO K=1,ORDEM_MEF(I)+1
            IF ((CONECTI_MEF(J,K) .EQ. CONECTI_MEF(I,1)) .OR. (CONECTI_MEF(J,K) .EQ. CONECTI_MEF(I,ORDEM_MEF(I)+1))) THEN
                cont = cont + 1
            ENDIF
        ENDDO
    ENDDO          
    cont = cont - 2  !subtraindo os encontros do proprio elemento
    N_INT_P_MEF(I) = (ORDEM_MEF(I)+1)
    IF (cont .LT. 2) N_INT_P_MEF(I) = 2*N_INT_P_MEF(I)    
!           WRITE(*,*)'ELEM:',I
!           WRITE(*,*)'CONT:',CONT,'NUMBER POINTS:',P_INT_MEF
!           READ(*,*)   
!           IF (P_INT_MEF .LT. ORDEM_MEF(I)+1) P_INT_MEF = ORDEM_MEF(I)+1   
            
    ! -------------------------------------------
                       
    ! STARTING GAUSS INTEGRATION VARIABLES-------
    ALLOCATE(GI_MEF(P_INT_MEF),OME_MEF(P_INT_MEF),QSIW_MEF(P_INT_MEF,2))
     
	GI_MEF=0.D0   !Array of length P_INT_MEF containing quadrature points.
	OME_MEF=0.D0  !Array of length P_INT_MEF containing quadrature weights.
	QSIW_MEF=0.D0 !Array with both quadrature points and weights
            
    CALL DGQRUL (P_INT_MEF,1,0.D0,0.D0,0,GI_MEF,GI_MEF,OME_MEF)
    
    QSIW_MEF(:,1)=GI_MEF(:)
	QSIW_MEF(:,2)=OME_MEF(:)   
          
	DEALLOCATE(GI_MEF,OME_MEF)
    ! --------------------------------------------
            
    DO J=1,P_INT_MEF ! PERFORMING NUMERICAL INTEGRATION ____________________________________
                
        ! COMPUTING THE DERIVATIVE VALUES OF THE SHAPE FUNCTIONS
        CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSIW_MEF(J,1),N_NO,I,VALORES,X,Y,Z,D_PHI_MEF)
                
        B=D_PHI_MEF
        BT=B
            
        KLOCALAUX=0.0D0
        CALL DMRRRR (N_NO,1,BT,N_NO,1,N_NO,B,1,N_NO,N_NO,KLOCALAUX,N_NO)
                               
        !AUX2=MP_MEF(I,1)*MP_MEF(I,2)*QSIW_MEF(J,2)
        KLOCAL=KLOCAL+KLOCALAUX*MP_MEF(I,1)*MP_MEF(I,2)*QSIW_MEF(J,2)	
        ! ------------------------------------------------------
                
        ! COMPUTING THE JACOBIAN (TAN_NORM AND 2/L_MEF) --------
        TAN1=0.0D0
        TAN2=0.0D0
        TAN3=0.0D0
        DO K=1,N_NO
		    TAN1=TAN1+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),1) !TANGENT VECTORS
		    TAN2=TAN2+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),2) !TANGENT VECTORS
            TAN3=TAN3+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),3) !TANGENT VECTORS
        ENDDO
        TAN_NORM=DSQRT(TAN1*TAN1+TAN2*TAN2+TAN3*TAN3) !TANGENT VECTOR NORM, JACOBIAN OF THE POINT
	    L_MEF=L_MEF+TAN_NORM*QSIW_MEF(J,2)  !THE LENGHT OF THE ELEMENT IS THE SUM OF JACOBIANS, WEIGHTED BY THE GAUSS POINT WEIGHTS
        ! ------------------------------------------------------
                
        ! COMPUTING THE APPROXIMATE FUNCTION VALUES FOR LUMPING MATRIX
        CALL FUNCOES_DE_FORMA_MEC_MEF(4,QSIW_MEF(J,1),N_NO,I,VALORES,X,Y,Z,PHI_MEF)
                
        B=PHI_MEF
        BT=B
                
        GLOCALAUX=0.D0
        CALL DMRRRR (N_NO,1,BT,N_NO,1,N_NO,B,1,N_NO,N_NO,GLOCALAUX,N_NO)
        GLOCAL=GLOCAL+GLOCALAUX*QSIW_MEF(J,2)	
        ! ------------------------------------------------------------
            
    ENDDO ! END OF NUMERICAL INTEGRATION ___________________________________________________
            
    !MULTIPLYING KLOCAL AND F_MEF BY RESPECTIVE JACOBIANS
    KLOCAL=KLOCAL*(2.D0/L_MEF)
	GLOCAL=GLOCAL*(L_MEF/2.D0)
    !----------------------------------------------------
    
    !CLOSING VARIABLES
    DEALLOCATE(B,BT,KLOCALAUX,D_PHI_MEF,PHI_MEF,GLOCALAUX,VALORES,QSIW_MEF)
                

END SUBROUTINE
        
    
    
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  
    
    
     
    
!SUBROUTINE ROTATE_MATRIXES(ELEM,KLOCAL,GLOCAL)
!    !
!    ! THIS SUBROUTINE PERFORMS THE ROTATION OF THE KLOCAL AND GLOCAL MATRIXES (INPUT)
!    !
!    ! ELEM IS THE FIBER ELEMENT NUMBER
!    ! OUTPUT OF THIS SUBROUTINE IS STORES IN THE GLOBAL VARIABLES KKGLOCAL AND GGGLOBAL
!    !
!
!    USE REINFORCEMENTS
!    
!    IMPLICIT NONE
!    
!    INTEGER::ELEM,I,J,K,N_NO
!    
!    REAL*8::KLOCAL(ORDEM_MEF(ELEM)+1,ORDEM_MEF(ELEM)+1),GLOCAL(ORDEM_MEF(ELEM)+1,ORDEM_MEF(ELEM)+1)
!
!    REAL*8,DIMENSION(:,:),ALLOCATABLE::KLOCALAUX,GLOCALAUX,MR,MRT,VALORES,KLOCALAUX_SLIP
!    
!    REAL*8::AUX1,AUX2,AUX3,AUX4,PHI_MEF[ALLOCATABLE](:),D_PHI_MEF[ALLOCATABLE](:),X,Y,Z
!    
!    LOGICAL::LINE_HAS_ZEROS
!    
!    ! GIVEN THE ELEMENT NUMBER:
!    I=ELEM
!    
!    ! EXPAND THE LOCAL STIFNESS MATRIX TO 3 FREEDOM DEGREES, KLOCAL -> KLOCALAUX
!    N_NO=ORDEM_MEF(I)+1
!    ALLOCATE(MR(N_NO,3*N_NO),MRT(3*N_NO,N_NO),PHI_MEF(N_NO),D_PHI_MEF(N_NO),VALORES(N_NO,3))
!    MR=0.D0     !ROTATION MATRIX
!    MRT=0.D0    !ROTATION MATRIX TRANSPOSED
!    ! ---------------------------------------------------------------------------
!    
!    ! ________________________ ROTATE THE LOCAL STIFNESS MATRIX KLOCAL -> KLOCALAUX ______________________
!    DO J=1,N_NO
!            
!        PHI_MEF=0.0D0
!        !CALL FUNCOES_DE_FORMA_MEC_MEF(3,QSI_MEF(I,J),(ORDEM_MEF(I)+1),I,VALORES,X,Y,Z,D_PHI_MEF)  
!        CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSI_MEF(I,J),(ORDEM_MEF(I)+1),I,VALORES,X,Y,Z,D_PHI_MEF)  ! MUDOU: MAIS COERENCIA (VALIDAR)(DUVIDA)
!            
!        AUX1=0.0D0
!        AUX2=0.0D0
!        AUX3=0.0D0
!        DO K=1,N_NO
!		    AUX1=AUX1+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),1)
!		    AUX2=AUX2+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),2)
!            AUX3=AUX3+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),3)
!		    !AUX1=AUX1+D_PHI_MEF(K)*COORDPCOLOC_MEF(CONECTI_MEF(I,K),1) 
!		    !AUX2=AUX2+D_PHI_MEF(K)*COORDPCOLOC_MEF(CONECTI_MEF(I,K),2)
!      !      AUX3=AUX3+D_PHI_MEF(K)*COORDPCOLOC_MEF(CONECTI_MEF(I,K),3)
!        ENDDO
!        AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
!        AUX1=AUX1/AUX4
!        AUX2=AUX2/AUX4
!        AUX3=AUX3/AUX4
!            
!        !TO AVOID NUMERICAL ERRORS:            
!        IF (ABS(AUX1).LT.TOLER) AUX1=0.D0
!        IF (ABS(AUX2).LT.TOLER) AUX2=0.D0
!        IF (ABS(AUX3).LT.TOLER) AUX3=0.D0
!
!        MR(J,3*J-2) = AUX1
!        MR(J,3*J-1) = AUX2
!        MR(J,3*J) = AUX3
!                               
!    ENDDO
!    MRT = TRANSPOSE(MR)
!        
!    ALLOCATE(KLOCALAUX(N_NO,3*N_NO))
!    KLOCALAUX=0.0D0
!        
!    KLOCALAUX=MATMUL(KLOCAL,MR)
!    
!    IF (SLIP_REINF) THEN
!        !ALLOCATE(KLOCALAUX_SLIP(3*N_NO,N_NO))
!        !KLOCALAUX_SLIP=MATMUL(MRT,KLOCAL)
!    ENDIF
!        
!    ! ____________________________ END OF ROTATING MATRIXES __________________________________________________
!        
!    ! SAVING ELEMENT MATRIXES INTO A UNIQUE VARIABLE
!    DO J=1,(N_NO)
!        DO K=1,(3*N_NO)
!            KKGLOBAL(I,J,K)=KLOCALAUX(J,K)
!        ENDDO
!        DO K=1,N_NO
!            GGGLOBAL(I,J,K)=GLOCAL(J,K)
!        ENDDO
!        !IF (SLIP_REINF) THEN
!        !    DO K=1,N_NO
!        !        KKGLOBAL_SLIP(I,J,K)=KLOCALAUX_SLIP(J,K)
!        !    ENDDO
!        !ENDIF
!    ENDDO
!        
!    ! CLOSING VARIABLES
!    DEALLOCATE(MR,MRT,PHI_MEF,D_PHI_MEF,VALORES,KLOCALAUX)
!    
!    !IF (SLIP_REINF) DEALLOCATE(KLOCALAUX_SLIP)
!
!END SUBROUTINE
    
        
    
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  

    
    
SUBROUTINE ROTATE_MATRIXES(ELEM,KLOCAL,GLOCAL)
    !
    ! THIS SUBROUTINE PERFORMS THE ROTATION OF THE KLOCAL AND GLOCAL MATRIXES (INPUT)
    !
    ! ELEM IS THE FIBER ELEMENT NUMBER
    ! OUTPUT OF THIS SUBROUTINE IS STORES IN THE GLOBAL VARIABLES KKGLOCAL AND GGGLOBAL
    !

    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::ELEM,I,J,K,N_NO
    
    REAL*8::KLOCAL(ORDEM_MEF(ELEM)+1,ORDEM_MEF(ELEM)+1),GLOCAL(ORDEM_MEF(ELEM)+1,ORDEM_MEF(ELEM)+1)

    REAL*8,DIMENSION(:,:),ALLOCATABLE::KLOCALAUX,GLOCALAUX,MR,MRT,VALORES,KLOCALAUX_SLIP
    
    REAL*8::AUX1,AUX2,AUX3,AUX4,PHI_MEF[ALLOCATABLE](:),D_PHI_MEF[ALLOCATABLE](:),X,Y,Z
    
    LOGICAL::LINE_HAS_ZEROS
    
    ! GIVEN THE ELEMENT NUMBER:
    I=ELEM
    
    ! EXPAND THE LOCAL STIFNESS MATRIX TO 3 FREEDOM DEGREES, KLOCAL -> KLOCALAUX
    N_NO=ORDEM_MEF(I)+1
    ALLOCATE(MR(3*N_NO,3*N_NO),MRT(3*N_NO,3*N_NO),PHI_MEF(N_NO),D_PHI_MEF(N_NO),VALORES(N_NO,3))
    MR=0.D0     !ROTATION MATRIX
    MRT=0.D0    !ROTATION MATRIX TRANSPOSED
    ! ---------------------------------------------------------------------------
    
    ! ________________________ ROTATE THE LOCAL STIFNESS MATRIX KLOCAL -> KLOCALAUX ______________________
    DO J=1,N_NO
            
        PHI_MEF=0.0D0
        !CALL FUNCOES_DE_FORMA_MEC_MEF(3,QSI_MEF(I,J),(ORDEM_MEF(I)+1),I,VALORES,X,Y,Z,D_PHI_MEF)  
        CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSI_MEF(I,J),(ORDEM_MEF(I)+1),I,VALORES,X,Y,Z,D_PHI_MEF)  ! MUDOU: MAIS COERENCIA (VALIDAR)(DUVIDA)
            
        AUX1=0.0D0
        AUX2=0.0D0
        AUX3=0.0D0
        DO K=1,N_NO
		    AUX1=AUX1+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),1)
		    AUX2=AUX2+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),2)
            AUX3=AUX3+D_PHI_MEF(K)*COORDNOS_MEF(CONECTI_MEF(I,K),3)
		    !AUX1=AUX1+D_PHI_MEF(K)*COORDPCOLOC_MEF(CONECTI_MEF(I,K),1) 
		    !AUX2=AUX2+D_PHI_MEF(K)*COORDPCOLOC_MEF(CONECTI_MEF(I,K),2)
      !      AUX3=AUX3+D_PHI_MEF(K)*COORDPCOLOC_MEF(CONECTI_MEF(I,K),3)
        ENDDO
        AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
        AUX1=AUX1/AUX4
        AUX2=AUX2/AUX4
        AUX3=AUX3/AUX4
            
        !TO AVOID NUMERICAL ERRORS:            
        IF (ABS(AUX1).LT.TOLER) AUX1=0.D0
        IF (ABS(AUX2).LT.TOLER) AUX2=0.D0
        IF (ABS(AUX3).LT.TOLER) AUX3=0.D0

        !MR(J,3*J-2) = AUX1
        !MR(J,3*J-1) = AUX2
        !MR(J,3*J) = AUX3
        
        if (AUX1.LT.TOLER .and. AUX3.LT.TOLER) then
            MR(3*J-2,3*J-1) = AUX2
            MR(3*J-1,3*J-2) = -AUX2
            MR(3*J,3*J) = 1.0d0
        else
            MR(3*J-2,3*J-2) = AUX1
            MR(3*J-2,3*J-1) = AUX2
            MR(3*J-2,3*J) = AUX3
            MR(3*J-1,3*J-2) = -AUX1*AUX2/dsqrt(AUX1*AUX1 + AUX3*AUX3)
            MR(3*J-1,3*J-1) = dsqrt(AUX1*AUX1 + AUX3*AUX3)
            MR(3*J-1,3*J) = -AUX2*AUX3/dsqrt(AUX1*AUX1 + AUX3*AUX3)
            MR(3*J,3*J-2) = -AUX3/dsqrt(AUX1*AUX1 + AUX3*AUX3)
            MR(3*J,3*J) = AUX1/dsqrt(AUX1*AUX1 + AUX3*AUX3)
        end if
        
                               
    ENDDO
    MRT = TRANSPOSE(MR)
        
    ALLOCATE(KLOCALAUX(3*N_NO,3*N_NO),GLOCALAUX(3*N_NO,3*N_NO))
    KLOCALAUX=0.0D0
    GLOCALAUX=0.0D0
    DO J=1,3*N_NO
        GLOCALAUX(J,J)=1.0D0
    ENDDO
    DO J=1,N_NO
        DO K=1,N_NO
            KLOCALAUX((3*J-2),(3*K-2))=KLOCAL(J,K)
            GLOCALAUX((3*J-2),(3*K-2))=GLOCAL(J,K)
        ENDDO
    ENDDO
    
    KLOCALAUX=MATMUL(MATMUL(MRT,KLOCALAUX),MR)
    GLOCALAUX=MATMUL(MATMUL(MRT,GLOCALAUX),MR)
    
    IF (SLIP_REINF) THEN
        !ALLOCATE(KLOCALAUX_SLIP(3*N_NO,N_NO))
        !KLOCALAUX_SLIP=MATMUL(MRT,KLOCAL)
    ENDIF
        
    ! ____________________________ END OF ROTATING MATRIXES __________________________________________________
        
    ! SAVING ELEMENT MATRIXES INTO A UNIQUE VARIABLE
    DO J=1,(3*N_NO)
        DO K=1,(3*N_NO)
            KKGLOBAL(I,J,K)=KLOCALAUX(J,K)
            GGGLOBAL(I,J,K)=GLOCALAUX(J,K)
        ENDDO
        IF (SLIP_REINF) THEN
            !DO K=1,N_NO
            !    KKGLOBAL_SLIP(I,J,K)=KLOCALAUX_SLIP(J,K)
            !ENDDO
        ENDIF
    ENDDO
        
    ! CLOSING VARIABLES
    DEALLOCATE(MR,MRT,PHI_MEF,D_PHI_MEF,VALORES,KLOCALAUX,GLOCALAUX)
    
    !IF (SLIP_REINF) DEALLOCATE(KLOCALAUX_SLIP)

END SUBROUTINE
    
    
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  
    
    
     
    
SUBROUTINE CONNECTION_ELEMENT
  
    USE REINFORCEMENTS
    USE VISCO
    
    IMPLICIT NONE
    
    INTEGER::I,J,JJ,K,KK,NO_CONECT(2),ELEM_CONECT(2),N_NO,NUM_LEI
    
    REAL*8::YOUNG_CONECT,AREA_CONECT,AUX1,AUX2,AUX3,L_CONECT,KLOCALAUX(2,2),GLOCALAUX(2,2),&
        K1,K2,K3,K4,GAMA,EVE,EE
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::KLOCAL,GLOCAL,MR,MRT,A_MATRIX,B_MATRIX,C_MATRIX,B_MATRIX_INV
    
    DO I=1,N_INTERFACE_MEF
            
        NO_CONECT(1) = NI_MEF(I,1)  !FIRST CONECTION ELEMENT NODE
        NO_CONECT(2) = NI_MEF(I,2)  !SECOND CONECTION ELEMENT NODE
        
        DO k=1,N_ELEMENTOS_MEF ! FINDING ELEMENTS OF NO_I AND NO_J
            IF((NO_CONECT(1) .EQ. CONECTI_MEF(K,1)) .OR. (NO_CONECT(1) .EQ. CONECTI_MEF(K,ORDEM_MEF(K)+1))) ELEM_CONECT(1)=K
     
            IF((NO_CONECT(2) .EQ. CONECTI_MEF(K,1)) .OR. (NO_CONECT(2) .EQ. CONECTI_MEF(K,ORDEM_MEF(K)+1))) ELEM_CONECT(2)=K
        ENDDO
                
        ! CONECTION ELEMENT PROPERTIES
        YOUNG_CONECT = (MP_MEF(ELEM_CONECT(1),1)+MP_MEF(ELEM_CONECT(2),1))/(2.0D0)*(1.0D0)  
        AREA_CONECT = (MP_MEF(ELEM_CONECT(1),2)+MP_MEF(ELEM_CONECT(2),2))/(2.0D0)*(1.0D0)
        
        ! CONECTION ELEMENT INCLINATION (CONSTANT)
        AUX1 = COORDPCOLOC_MEF(NO_CONECT(2),1)-COORDPCOLOC_MEF(NO_CONECT(1),1)
        AUX2 = COORDPCOLOC_MEF(NO_CONECT(2),2)-COORDPCOLOC_MEF(NO_CONECT(1),2)
        AUX3 = COORDPCOLOC_MEF(NO_CONECT(2),3)-COORDPCOLOC_MEF(NO_CONECT(1),3)
        L_CONECT = DSQRT(AUX1**2.0D0+AUX2**2.0D0+AUX3**2.0D0)
                
        ALLOCATE(KLOCAL(6,6),GLOCAL(6,6))
        
        IF (REINFORCEMENTS_METHOD.EQ.'FEM') THEN
        
            ! CALCULATION ROTATION MATRIX
            N_NO=2
            ALLOCATE(MR(N_NO,3*N_NO),MRT(3*N_NO,N_NO))
            
            ! ANGLES
            AUX1=AUX1/L_CONECT
            AUX2=AUX2/L_CONECT
            AUX3=AUX3/L_CONECT
                        
            !TO AVOID NUMERICAL ERRORS:            
            IF (ABS(AUX1).LT.TOLER) AUX1=0.D0
            IF (ABS(AUX2).LT.TOLER) AUX2=0.D0
            IF (ABS(AUX3).LT.TOLER) AUX3=0.D0
                        
            MR=0.0D0
            DO J=1,N_NO
                MR(J,3*J-2) = AUX1
                MR(J,3*J-1) = AUX2
                MR(J,3*J) = AUX3
            ENDDO
            MRT=TRANSPOSE(MR)
            
            ! ENTERING STIFFNESS MATRIX             
            KLOCALAUX(1,1) = 1.0D0
            KLOCALAUX(1,2) = -1.0D0
            KLOCALAUX(2,1) = -1.0D0
            KLOCALAUX(2,2) = 1.0D0
            
            DO K=1,2
                DO KK=1,2
                    KLOCALAUX(K,KK) = KLOCALAUX(K,KK)*YOUNG_CONECT*AREA_CONECT/L_CONECT
                ENDDO
            ENDDO
            
            ! ROTATING SITFFNESS MATRIX 
            KLOCAL=MATMUL(MRT,MATMUL(KLOCALAUX,MR))
            
            ! ENTERING LUMPING MARTRIX           
            GLOCALAUX(1,1) = L_CONECT/3.0D0
            GLOCALAUX(1,2) = L_CONECT/6.0D0
            GLOCALAUX(2,1) = L_CONECT/6.0D0
            GLOCALAUX(2,2) = L_CONECT/3.0D0
            
            ! ROTATIONG LUMPING MATRIX 
            GLOCAL=MATMUL(MRT,MATMUL(GLOCALAUX,MR))
            
            DEALLOCATE(MR,MRT)
            
        ELSE IF (REINFORCEMENTS_METHOD.EQ.'1DBEM') THEN
            
            ! CALCULATION ROTATION MATRIX 
            N_NO = 2
            ALLOCATE(MR(N_NO,3*N_NO),MRT(3*N_NO,N_NO))
            
            ! ANGLES
            AUX1=AUX1/L_CONECT
            AUX2=AUX2/L_CONECT
            AUX3=AUX3/L_CONECT
            
            !TO AVOID NUMERICAL ERRORS:            
            IF (ABS(AUX1).LT.TOLER) AUX1=0.D0
            IF (ABS(AUX2).LT.TOLER) AUX2=0.D0
            IF (ABS(AUX3).LT.TOLER) AUX3=0.D0
                        
            MR=0.0D0
            DO J=1,N_NO
                MR(J,3*J-2) = AUX1
                MR(J,3*J-1) = AUX2
                MR(J,3*J) = AUX3
            ENDDO
            MRT=TRANSPOSE(MR)
            
            ! FINDING STIFFNESS MATRIX K 
            ALLOCATE(A_MATRIX(N_NO,N_NO),B_MATRIX(N_NO,N_NO),B_MATRIX_INV(N_NO,N_NO),C_MATRIX(N_NO,N_NO))        
 
            ! ENTERING H AND G MATRIXES
            A_MATRIX(1,1)=0.50D0
            A_MATRIX(1,2)=-0.50D0
            A_MATRIX(2,1)=-0.50D0
            A_MATRIX(2,2)=0.50D0
            
            B_MATRIX(1,1)=0.0D0
            B_MATRIX(1,2)=-L_CONECT/(2.0D0*YOUNG_CONECT*AREA_CONECT)
            B_MATRIX(2,1)=-L_CONECT/(2.0D0*YOUNG_CONECT*AREA_CONECT)
            B_MATRIX(2,2)=0.0D0
            
            CALL DLINRG(SIZE(B_MATRIX,2),B_MATRIX,SIZE(B_MATRIX,1),B_MATRIX_INV,SIZE(B_MATRIX_INV,1))
            
            KLOCALAUX=MATMUL(B_MATRIX_INV,A_MATRIX)
            
            ! ENTRING G BAR MATRIX 
            C_MATRIX(1,1) = 1.0D0/6.0D0
            C_MATRIX(1,2) = 1.0D0/3.0D0
            C_MATRIX(2,1) = 1.0D0/3.0D0
            C_MATRIX(2,2) = 1.0D0/6.0D0
            
            C_MATRIX=C_MATRIX*(-L_CONECT**2.0D0/(2.0D0*YOUNG_CONECT*AREA_CONECT))
            
            GLOCALAUX=MATMUL(B_MATRIX_INV,C_MATRIX)
            
            ! ROTATING MATRIXES
            KLOCAL=MATMUL(MRT,MATMUL(KLOCALAUX,MR))
            GLOCAL=MATMUL(MRT,MATMUL(GLOCALAUX,MR))
            
            DEALLOCATE(A_MATRIX,B_MATRIX,B_MATRIX_INV,C_MATRIX,MR,MRT)
            
        ENDIF               
        
        ! SAVING CONECTION ELEMENT MATRIXES IN GLOBAL VARIABLES
        DO K=1,6
            DO J=1,6
                K_CONECT(I,K,J)=KLOCAL(K,J)
                G_CONECT(I,K,J)=GLOCAL(K,J)
            ENDDO
        ENDDO
        
        ! GETTING VISCOELASTIC CONSTRAINTS
        IF (VISCO_ANALYSIS_REINF) THEN
            NUM_LEI = TIME_MODEL_REINF(ELEM_CONECT(1))
            GAMA = MP_MEF_VISCO(ELEM_CONECT(1),1)
            Eve = MP_MEF_VISCO(ELEM_CONECT(1),2)
            Ee   = MP_MEF(ELEM_CONECT(1),1)
            SELECT CASE(NUM_LEI)
            CASE(1)!KELVIN'S MODEL
	          K1= 1.0D0 + GAMA/DELTA_T_VISCO
	          K2= 1.0D0
	          K3= GAMA/DELTA_T_VISCO
	          K4= 0.0D0  
            CASE(2)!BOLTZMANN'S MODEL
              K1= 1.0D0 + GAMA/DELTA_T_VISCO
	          K2= GAMA/DELTA_T_VISCO + (Ee + Eve)/Eve
	          K3= GAMA/DELTA_T_VISCO
	          K4= -GAMA/DELTA_T_VISCO    
            CASE(3)!MAXWELL'S MODEL
              K1= 1.0D0
	          K2= 1.0D0 + DELTA_T_VISCO/GAMA
	          K3= 1.0D0
	          K4=-1.0D0 
            CASE(4)!HOOKE'S MODEL
              K1= 1.0D0
	          K2= 1.0D0
	          K3= 0.0D0
	          K4= 0.0D0  
            ENDSELECT
        ENDIF
        
        ! TRANSFERING CONECTIOM ELEMENT'S KLOCAL AND GLOCAL TO GLOBAL ONES -> KGLOBAL E GGLOBAL
        DO J=1,2
            DO JJ=1,3
                DO K=1,2
                    DO KK=1,3
                        IF (VISCO_ANALYSIS_REINF) THEN
                            KGLOBAL(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)=&
                            KGLOBAL(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)+&
                            KKGLOBAL(I,3*(J-1)+JJ,3*(K-1)+KK)*K1
                        
                            GGLOBAL(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)=&
                            GGLOBAL(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)+&
                            GGGLOBAL(I,3*(J-1)+JJ,3*(K-1)+KK)*K2
                            !------------------------------------------------------------------
                            K_KGLOBAL(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)=&
                            K_KGLOBAL(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)+&
                            KKGLOBAL(I,3*(J-1)+JJ,3*(K-1)+KK)*K3
                        
                            G_GGLOBAL(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)=&
                            G_GGLOBAL(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)+&
                            GGGLOBAL(I,3*(J-1)+JJ,3*(K-1)+KK)*K4
                        ELSE
                            KGLOBAL(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)=&
                            KGLOBAL(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)+&
                            KLOCAL(3*(J-1)+JJ,3*(K-1)+KK)
                        
                            GGLOBAL(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)=&
                            GGLOBAL(3*(NO_CONECT(J)-1)+JJ,3*(NO_CONECT(K)-1)+KK)+&
                            GLOCAL(3*(J-1)+JJ,3*(K-1)+KK)
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
        ENDDO     
        
        DEALLOCATE(KLOCAL,GLOCAL)
        
    ENDDO ! END OF CONECTION ELEMENT _____________________________________________________________________________________________

        
  
END SUBROUTINE
    
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************  
    
    
     
    
SUBROUTINE CONNECTION_ELEMENT_FOR_F(KG,GG,COORDS)
  
    ! COORDS(I,J): COORDENADA J DO NÓ I PERTERCENTE AO CONNECTION ELEMENT

    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    REAL*8,INTENT(IN)::COORDS(2,3)
    REAL*8,INTENT(OUT)::KG(6,6),GG(6,6)
    
    INTEGER::I,J,JJ,K,KK,NO_CONECT(2),ELEM_CONECT(2),N_NO
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::A_MATRIX,B_MATRIX,C_MATRIX,B_MATRIX_INV
    REAL*8::MR(2,6),MRT(6,2),KLOCALAUX(2,2),GLOCALAUX(2,2),YOUNG_CONECT,AREA_CONECT,AUX(4),L_CONECT
    
    ! FINDING ELEMENTS CONNECTED TO THE CONNECTION ELEMENT   
    NO_CONECT=0
    ELEM_CONECT=0
    DO I=1,N_ELEMENTOS_MEF 
        DO J=1,ORDEM_MEF(I)+1
            AUX=0.0D0
            DO K=1,3
                AUX(1)=AUX(1)+(COORDPCOLOC_MEF(CONECTI_MEF(I,J),K)-COORDS(1,K))**2.0D0
                AUX(2)=AUX(2)+(COORDPCOLOC_MEF(CONECTI_MEF(I,J),K)-COORDS(2,K))**2.0D0
            ENDDO
            AUX(1)=DSQRT(AUX(1))
            AUX(2)=DSQRT(AUX(2))
            IF (AUX(1).LT.TOLER) THEN
                NO_CONECT(1)=J
                ELEM_CONECT(1)=I
            ENDIF
            IF (AUX(2).LT.TOLER) THEN
                NO_CONECT(2)=J
                ELEM_CONECT(2)=I
            ENDIF
        ENDDO
    ENDDO
                
    ! CONECTION ELEMENT PROPERTIES
    IF (ELEM_CONECT(1).NE.0) THEN
        YOUNG_CONECT = MP_MEF(ELEM_CONECT(1),1)
        AREA_CONECT = MP_MEF(ELEM_CONECT(1),2)
    ELSE IF (ELEM_CONECT(2).NE.0) THEN
        YOUNG_CONECT = MP_MEF(ELEM_CONECT(2),1)
        AREA_CONECT = MP_MEF(ELEM_CONECT(2),2)
    ENDIF
    
    ! CONECTION ELEMENT INCLINATION (CONSTANT)
    DO K=1,3
        AUX(K)=COORDS(2,K)-COORDS(1,K)
    ENDDO
    L_CONECT = DSQRT(AUX(1)**2.0D0+AUX(2)**2.0D0+AUX(3)**2.0D0)
    
    ! utiliza L_REAL:
    !L_CONECT=80.0D0*AREA_CONECT/201.0D0
                    
    IF (REINFORCEMENTS_METHOD.EQ.'FEM') THEN
        
        ! CALCULATION ROTATION MATRIX 
        AUX(1)=AUX(1)/L_CONECT
        AUX(2)=AUX(2)/L_CONECT
        AUX(3)=AUX(3)/L_CONECT   
        MR=0.0D0
        MR(1,1)=AUX(1)
        MR(1,2)=AUX(2)
        MR(1,3)=AUX(3)
        MR(2,4)=AUX(1)
        MR(2,5)=AUX(2)
        MR(2,6)=AUX(3)
        MRT=TRANSPOSE(MR)
            
        ! ENTERING STIFFNESS MATRIX             
        KLOCALAUX(1,1) = 1.0D0
        KLOCALAUX(1,2) = -1.0D0
        KLOCALAUX(2,1) = -1.0D0
        KLOCALAUX(2,2) = 1.0D0
            
        DO K=1,2
            DO KK=1,2
                KLOCALAUX(K,KK) = KLOCALAUX(K,KK)*YOUNG_CONECT*AREA_CONECT/L_CONECT
            ENDDO
        ENDDO
            
        ! ROTATING SITFFNESS MATRIX 
        KG=MATMUL(MRT,MATMUL(KLOCALAUX,MR))
            
        ! ENTERING LUMPING MARTRIX           
        GLOCALAUX(1,1) = L_CONECT/3.0D0
        GLOCALAUX(1,2) = L_CONECT/6.0D0
        GLOCALAUX(2,1) = L_CONECT/6.0D0
        GLOCALAUX(2,2) = L_CONECT/3.0D0
            
        ! ROTATIONG LUMPING MATRIX 
        GG=MATMUL(MRT,MATMUL(GLOCALAUX,MR))
                        
    ELSE IF (REINFORCEMENTS_METHOD.EQ.'1DBEM') THEN
            
        ! CALCULATION ROTATION MATRIX             
        AUX(1)=AUX(1)/L_CONECT
        AUX(2)=AUX(2)/L_CONECT
        AUX(3)=AUX(3)/L_CONECT    
        MR=0.0D0
        MR(1,1)=AUX(1)
        MR(1,2)=AUX(2)
        MR(1,3)=AUX(3)
        MR(2,4)=AUX(1)
        MR(2,5)=AUX(2)
        MR(2,6)=AUX(3)
        MRT=TRANSPOSE(MR)
            
        ! FINDING STIFFNESS MATRIX K 
        N_NO = 2
        ALLOCATE(A_MATRIX(N_NO,N_NO),B_MATRIX(N_NO,N_NO),B_MATRIX_INV(N_NO,N_NO),C_MATRIX(N_NO,N_NO))        
 
        ! ENTERING H AND G MATRIXES
        A_MATRIX(1,1)=0.50D0
        A_MATRIX(1,2)=-0.50D0
        A_MATRIX(2,1)=-0.50D0
        A_MATRIX(2,2)=0.50D0
            
        B_MATRIX(1,1)=0.0D0
        B_MATRIX(1,2)=-L_CONECT/(2.0D0*YOUNG_CONECT*AREA_CONECT)
        B_MATRIX(2,1)=-L_CONECT/(2.0D0*YOUNG_CONECT*AREA_CONECT)
        B_MATRIX(2,2)=0.0D0
            
        CALL DLINRG(SIZE(B_MATRIX,2),B_MATRIX,SIZE(B_MATRIX,1),B_MATRIX_INV,SIZE(B_MATRIX_INV,1))
            
        KLOCALAUX=MATMUL(B_MATRIX_INV,A_MATRIX)
            
        ! ENTRING G BAR MATRIX 
        C_MATRIX(1,1) = 1.0D0/6.0D0
        C_MATRIX(1,2) = 1.0D0/3.0D0
        C_MATRIX(2,1) = 1.0D0/3.0D0
        C_MATRIX(2,2) = 1.0D0/6.0D0
            
        C_MATRIX=C_MATRIX*(-L_CONECT**2.0D0/(2.0D0*YOUNG_CONECT*AREA_CONECT))
            
        GLOCALAUX=MATMUL(B_MATRIX_INV,C_MATRIX)
            
        ! ROTATING MATRIXES
        KG=MATMUL(MRT,MATMUL(KLOCALAUX,MR))
        GG=MATMUL(MRT,MATMUL(GLOCALAUX,MR))
            
        DEALLOCATE(A_MATRIX,B_MATRIX,B_MATRIX_INV,C_MATRIX)
            
    ENDIF               
    
END SUBROUTINE