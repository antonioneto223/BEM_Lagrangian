SUBROUTINE FIBER_NONSING_S(ELEM,Nu,Mu,COORDS,DG_CF)
    
    ! FUNCTION THAT INTEGRATES U* FUNDAMENTAL SOLUTION OVER A GIVEN FIBER ELEMENT (ELEM) 
    ! HAVING A SOURCE POINT THAT IS OUTSIDE THE ELEMENT
    
    ! Nu AND Mu: POISSION COEFICIENT AND TRANSVERSAL ELASTIC MODULUS
    ! COORDS: X,Y AND Z COORDINATES OF THE SOURCE POINT
    ! DG_CF: MATRIX TO STORE THE RESULTS OF THIS SUBROUTINE
    
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER::ELEM,K,KK,JJ,P_INTEGRACAO,N
    
    REAL*8::Nu,Mu,COORDS(3),DG_CF(3,3*(ORDEM_MEF(ELEM)+1))
    
    REAL*8::GG_LOCAL(3,3),KR1(3,3),AUX1,AUX2,AUX3,JC,PI1,C1,C2,X,Y,Z,DX,DY,DZ,R,DR(3),AVERAGE_LENGTH,DISTANCE
    
    REAL*8,DIMENSION(:),ALLOCATABLE::PHI
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::DG_AUX,MATPHI,VALORES,QSIW
    
    ! STARTING VARIABLES
    ALLOCATE(DG_AUX(3,3*(ORDEM_MEF(ELEM)+1)),VALORES(ORDEM_MEF(ELEM)+1,3),PHI(ORDEM_MEF(ELEM)+1),&
        MATPHI(3,3*(ORDEM_MEF(ELEM)+1)))
    DG_AUX=0.0D0
    DG_CF=0.0D0
    MATPHI=0.0D0
    PI1=DACOS(-1.D0)
    KR1=0.0D0
    KR1(1,1)=1.0D0
    KR1(2,2)=1.0D0
    KR1(3,3)=1.0D0
    
    C1=((1.D0)/(16.D0*PI1*Mu*(1.D0-Nu)))
	C2=((-1.D0)/(8.D0*PI1*(1.D0-Nu)))
    
    DO K=1,ORDEM_MEF(ELEM)+1
        VALORES(K,1)=COORDNOS_MEF(CONECTI_MEF(ELEM,K),1)
        VALORES(K,2)=COORDNOS_MEF(CONECTI_MEF(ELEM,K),2)
        VALORES(K,3)=COORDNOS_MEF(CONECTI_MEF(ELEM,K),3)
    ENDDO
    
    ! _______________________ ADAPTATIVE GAUSS POINTS ____________________________
    ! ____________________________________________________________________________
    
    AVERAGE_LENGTH=0.0D0     ! Calcula comprimento aproximado
    DO K=1,ORDEM_MEF(ELEM)
        AUX1=VALORES(K+1,1)-VALORES(K,1)
        AUX2=VALORES(K+1,2)-VALORES(K,2)
        AUX3=VALORES(K+1,3)-VALORES(K,3)
        AVERAGE_LENGTH=AVERAGE_LENGTH+DSQRT(AUX1**2.0D0+AUX2**2.0D0+AUX3**2.0D0)
    ENDDO   
    DISTANCE=0.D0
    DO K=1,ORDEM_MEF(ELEM)+1
        AUX1=DSQRT((COORDS(1)-VALORES(K,1))**2.0D0+(COORDS(2)-VALORES(K,2))**2.0D0+(COORDS(3)-VALORES(K,3))**2.0D0)
        IF(K.EQ.1)THEN
            DISTANCE=AUX1
        ELSE
            IF(AUX1.LT.DISTANCE)THEN
                DISTANCE=AUX1
            ENDIF
        ENDIF
    ENDDO
    N=DNINT(15*DSQRT(10/(DISTANCE/AVERAGE_LENGTH)))
    IF(N.LT.15)THEN
        N=15
    ENDIF	
    P_INTEGRACAO=N
    ALLOCATE(QSIW(P_INTEGRACAO,2))
    CALL GAUSS_POINTS(P_INTEGRACAO,QSIW)    
    ! _______________________ END OF ADAPTATIVE GAUSS POINTS _______________________
    ! ______________________________________________________________________________

    DO JJ=1,P_INTEGRACAO 
                        
        CALL FUNCOES_DE_FORMA_MEC_MEF(5,QSIW(JJ,1),ORDEM_MEF(ELEM)+1,ELEM,VALORES,X,Y,Z,PHI)
        DX=X-COORDS(1)
        DY=Y-COORDS(2)
        DZ=Z-COORDS(3)
        R=DSQRT(DX*DX+DY*DY+DZ*DZ)
        DR(1)=DX/R
		DR(2)=DY/R
		DR(3)=DZ/R
                        
        CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSIW(JJ,1),ORDEM_MEF(ELEM)+1,ELEM,VALORES,X,Y,Z,PHI)
        AUX1=0.0D0
        AUX2=0.0D0
        AUX3=0.0D0
        DO K=1,(ORDEM_MEF(ELEM)+1)
            AUX1=AUX1+PHI(K)*VALORES(K,1)
            AUX2=AUX2+PHI(K)*VALORES(K,2)
            AUX3=AUX3+PHI(K)*VALORES(K,3)
        ENDDO
        JC=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)

        DO K=1,3
            DO KK=1,3
                AUX1=0.D0
				AUX1=(3.D0-4.D0*Nu)*KR1(K,KK)
				AUX1=AUX1+DR(K)*DR(KK)
				GG_LOCAL(K,KK)=(1.D0/R)*C1*AUX1*QSIW(JJ,2)*JC   
            ENDDO
        ENDDO
                        
        CALL FUNCOES_DE_FORMA_MEC_MEF(2,QSIW(JJ,1),ORDEM_MEF(ELEM)+1,ELEM,VALORES,X,Y,Z,PHI)
        DO K=1,ORDEM_MEF(ELEM)+1
            MATPHI(1,(3*K-2))=PHI(K)
            MATPHI(2,(3*K-1))=PHI(K)
            MATPHI(3,(3*K))=PHI(K)
        ENDDO
                        
        DG_AUX=0.0D0
        CALL DMRRRR (3,3,GG_LOCAL,3,3,3*(ORDEM_MEF(ELEM)+1),MATPHI,3,3,3*(ORDEM_MEF(ELEM)+1),DG_AUX,3)
                         
        DG_CF=DG_CF+DG_AUX
                        
    ENDDO ! end of integration point
    
    ! CLOSING VARIABLES
    DEALLOCATE(DG_AUX,VALORES,PHI,MATPHI)
    
    !if (elem.eq.40) then
    !    write(10,*)dg_cf
    !endif
    !
    
END SUBROUTINE
    
    
        
    
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************
    
  
    
    
SUBROUTINE FIBER_SING_S(ELEM,Nu,Mu,P_INTEGRACAO,QSIW,COORDS,DG_FF)
    
    ! FUNCTION THAT INTEGRATES U* FUNDAMENTAL SOLUTION OVER A GIVEN FIBER ELEMENT (ELEM) 
    ! HAVING A SOURCE POINT INSIDE THIS FIBER ELEMENT (SINGULARITY)
    !
    ! Nu AND Mu: POISSION COEFICIENT AND TRANSVERSAL ELASTIC MODULUS
    ! P_INTEGRACAO: NUMBER OF INTEGRATION POINTS FOR GAUSS-LEGENGRE NUMERICAL INTEGRATION
    ! QSIW: INTEGRATION POINTS COORDINATES AND WEIGHTS
    ! COORDS: X,Y AND Z COORDINATES OF THE SOURCE POINT
    ! DG_FF: MATRIX TO STORE THE RESULTS OF THIS SUBROUTINE
    ! DR0: VECTOR OF DERIVATIVES IN X,Y,Z OF THE VECTOR R AT THE SOURCE POINT

    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    INTEGER,INTENT(IN)::ELEM,P_INTEGRACAO
    
    INTEGER::I,J,K,KK,JJ,P_I,LADO,P_I_THETA
    
    REAL*8,INTENT(IN)::COORDS(3),Nu,Mu,QSIW(P_INTEGRACAO,2)
    
    REAL*8,INTENT(OUT)::DG_FF(3,3*(ORDEM_MEF(ELEM)+1))
    
    REAL*8,DIMENSION(:),ALLOCATABLE::PHI,PHI_MEF,D_PHI_MEF,MAT_N
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE::VALORES,MATPHI,DG_AUX,RT,QSIW_THETA
        
    REAL*8::C1,C2,PI1,DX,DY,DZ,X,Y,Z,R,DR(3),JC,KR1(3,3),GG_LOCAL(3,3),JAC_T,THETA,&
        AUX1,AUX2,AUX3,AUX(3),XK,YK,ZK,RAIO,GG_ACUM(3,3)
    
    ! STARTING VARIABLES
    ALLOCATE(PHI_MEF(ORDEM_MEF(ELEM)+1),VALORES(ORDEM_MEF(ELEM)+1,3),PHI(ORDEM_MEF(ELEM)+1),&
        MATPHI(3,3*(ORDEM_MEF(ELEM)+1)),DG_AUX(3,3*(ORDEM_MEF(ELEM)+1)),D_PHI_MEF(ORDEM_MEF(ELEM)+1),&
        RT(3,3),MAT_N(3))
    DG_FF=0.0D0
    PI1=DACOS(-1.D0)
    C1=((1.D0)/(16.D0*PI1*Mu*(1.D0-Nu)))
	C2=((-1.D0)/(8.D0*PI1*(1.D0-Nu)))
    KR1=0.0D0
    KR1(1,1)=1.0D0
    KR1(2,2)=1.0D0
    KR1(3,3)=1.0D0
        
    DO K=1,(ORDEM_MEF(ELEM)+1)
        VALORES(K,1)=COORDNOS_MEF(CONECTI_MEF(ELEM,K),1)
        VALORES(K,2)=COORDNOS_MEF(CONECTI_MEF(ELEM,K),2)
        VALORES(K,3)=COORDNOS_MEF(CONECTI_MEF(ELEM,K),3)
    ENDDO
    
    P_I=P_INTEGRACAO
    P_I_THETA=20
    ALLOCATE(QSIW_THETA(P_I_THETA,2))
    CALL GAUSS_POINTS(P_I_THETA,QSIW_THETA)
    RAIO=DSQRT(MP_MEF(ELEM,2)/PI1)        
    !RAIO=DSQRT((MP_MEF(ELEM,2)*10.0d0)/PI1)        
    JAC_T=PI1*RAIO
    
    DO I=1,P_I   ! Starting numerical integration
               
        ! JACOBIAN (TANGENT VECTOR)
        CALL FUNCOES_DE_FORMA_MEC_MEF(6,QSIW(I,1),ORDEM_MEF(ELEM)+1,ELEM,VALORES,X,Y,Z,D_PHI_MEF)
        AUX(1)=0.0D0
        AUX(2)=0.0D0
        AUX(3)=0.0D0
        DO K=1,(ORDEM_MEF(ELEM)+1)
            AUX(1)=AUX(1)+D_PHI_MEF(K)*VALORES(K,1)
            AUX(2)=AUX(2)+D_PHI_MEF(K)*VALORES(K,2)
            AUX(3)=AUX(3)+D_PHI_MEF(K)*VALORES(K,3)
        ENDDO
        JC=DSQRT(AUX(1)*AUX(1)+AUX(2)*AUX(2)+AUX(3)*AUX(3))
        AUX(1)=AUX(1)/JC
        AUX(2)=AUX(2)/JC
        AUX(3)=AUX(3)/JC
        
        ! Posicao do ponto K (pto campo na linha central do cilindro):
        CALL FUNCOES_DE_FORMA_MEC_MEF(5,QSIW(I,1),ORDEM_MEF(ELEM)+1,ELEM,VALORES,XK,YK,ZK,D_PHI_MEF)
        
        GG_ACUM=0.0D0
        
        DO J=1,P_I_THETA  ! Starting integration in theta
                        
            ! Valo de theta:
            THETA=PI1*(QSIW_THETA(J,1)+1.0D0)
            
            ! Posicao do ponto campo:
            IF ((DABS(AUX1).LT.TOLER).AND.(DABS(AUX2).LT.TOLER)) THEN
                !MAT_N(1)=0.0D0
                !MAT_N(2)=DCOS(THETA)
                !MAT_N(3)=DSIN(THETA)
                MAT_N(1)=-DSIN(THETA)
                MAT_N(2)=DCOS(THETA)
                MAT_N(3)=0.0D0
            ELSE
                MAT_N(1)=-AUX(2)*DCOS(THETA)-AUX(1)*AUX(3)*DSIN(THETA)
                MAT_N(2)=AUX(1)*DCOS(THETA)-AUX(2)*AUX(3)*DSIN(THETA)
                MAT_N(3)=DSIN(THETA)*(AUX(1)**2.0D0+AUX(2)**2.0D0)
                MAT_N(:)=MAT_N(:)*(1.0D0/DSQRT(AUX(1)**2.0D0+AUX(2)**2.0D0))
            ENDIF
            X=XK+RAIO*MAT_N(1)
            Y=YK+RAIO*MAT_N(2)
            Z=ZK+RAIO*MAT_N(3)
            
            ! Nucleo de integracao
            DX=X-COORDS(1)
            DY=Y-COORDS(2)
            DZ=Z-COORDS(3)
            R=DSQRT(DX*DX+DY*DY+DZ*DZ)
            DR(1)=DX/R
		    DR(2)=DY/R
		    DR(3)=DZ/R                       
                        
            GG_LOCAL=0.0D0
            DO K=1,3
                DO KK=1,3
                    AUX1=0.D0
			        AUX1=(3.D0-4.D0*Nu)*KR1(K,KK)
			        AUX1=AUX1+DR(K)*DR(KK)
			        GG_LOCAL(K,KK)=(1.D0/R)*C1*AUX1*QSIW(I,2)*QSIW_THETA(J,2)*JC*JAC_T
                ENDDO
            ENDDO
            
            !if ((elem.eq.2).and.(i.eq.1)) write(10,*)qsiw_theta(j,1),gg_local(1,2)/QSIW_theta(j,2)
            
            GG_ACUM = GG_ACUM + GG_LOCAL
            
        ENDDO ! End of theta integration
                               
        !IF (elem.eq.2) write(10,*)qsiw(i,1),gg_acum(1,1)/QSIW(I,2)
        
        CALL FUNCOES_DE_FORMA_MEC_MEF(2,QSIW(I,1),ORDEM_MEF(ELEM)+1,ELEM,VALORES,X,Y,Z,PHI)
        MATPHI=0.0D0
        DO K=1,ORDEM_MEF(ELEM)+1
            MATPHI(1,(3*K-2))=PHI(K)
            MATPHI(2,(3*K-1))=PHI(K)
            MATPHI(3,(3*K))=PHI(K)
        ENDDO
                        
        DG_AUX=0.0D0     
        DG_AUX=MATMUL(GG_ACUM,MATPHI)
        DG_FF=DG_FF+DG_AUX
        
    ENDDO ! end of integration point loop
    
    DG_FF=DG_FF/(2.0D0*PI1*RAIO)
                    
    ! CLOSING VARIABLES
    DEALLOCATE(PHI_MEF,VALORES,PHI,MATPHI,DG_AUX,D_PHI_MEF,RT,MAT_N)
    
    !if (elem.eq.40) then        
    !    write(10,*)dg_ff
    !endif
        
END SUBROUTINE
    
    
        
    
    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************    
! ***********************************************************************************************************************************************************************************
! ***********************************************************************************************************************************************************************************
  
    
        
  
SUBROUTINE FIBER_INT_DECISION(CIL,ELEM,COORDS)

    !
    ! SUBROUTINE TO DECIDE IF THIS FIBER ELEMENT MUST BE INTEGRATED OVER CILINDRICAL SURFACE OR LINE 
    ! THE LINE INTEGRATION IS USED TO ELEMENTS FAR FROM SOURCE POINT
    ! CIL = OUTPUT VARIABLE (.TRUE. FOR CILINDRICAL INTEGRATION AND .FALSE. FOR LINE)
    ! ELEM: FIBER ELEMENT NUMBER
    ! COORDS: SOURCE POINT COORDINATES
    ! NOF: SOURCE POINT (GLOBAL NODE POSITION)
    !
    
    USE REINFORCEMENTS
    
    IMPLICIT NONE
    
    LOGICAL,INTENT(OUT)::CIL
    INTEGER,INTENT(IN)::ELEM
    REAL*8,INTENT(IN)::COORDS(3)
    
    INTEGER::I,J,K
    REAL*8::AUX(2),L_MEF,MENOR_DIST,DIST_PARA_INT_LINHA,RAIO
    
    ! PARAMETRO:
    DIST_PARA_INT_LINHA=100.00d0
    RAIO=DSQRT(MP_MEF(ELEM,2)/DACOS(-1.D0)) 
    
    ! Estimativa do vetor x (distancia do ponto fonte ate o elemento integrado)
    AUX=0.0D0
    DO I=1,3        
        K=CONECTI_MEF(ELEM,ORDEM_MEF(ELEM)+1)
        AUX(1)=AUX(1)+(COORDNOS_MEF(K,I)-COORDS(I))**2.0D0
        K=CONECTI_MEF(ELEM,1)
        AUX(2)=AUX(2)+(COORDNOS_MEF(K,I)-COORDS(I))**2.0D0
    ENDDO
    AUX(1)=DSQRT(AUX(1))
    AUX(2)=DSQRT(AUX(2))
        
    IF (DABS(AUX(1)).LT.DABS(AUX(2))) THEN
        MENOR_DIST=AUX(1)
    ELSE
        MENOR_DIST=AUX(2)
    ENDIF
    
    MENOR_DIST=MENOR_DIST/RAIO
    
    IF (MENOR_DIST.LT.DIST_PARA_INT_LINHA) THEN     ! INTEGRA CILINDRO
        CIL=.TRUE.
    ELSE
        CIL=.FALSE.                                 ! INTEGRA LINHA
    ENDIF
    
END SUBROUTINE