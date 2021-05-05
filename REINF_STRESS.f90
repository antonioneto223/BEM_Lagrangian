SUBROUTINE REINFORCEMENTS_STRESS
    
    !
    !   CALCULATE VALUE OF NORMAL EFFORTS AT EACH COLLOCATION POINT
    !   STORE IN E_NORMAL(I,J) VARIABLE (ELEMENT I AND LOCAL COLLOCATION POINT J)
    !
    
    USE REINFORCEMENTS    
    USE OMP_LIB
        
    IMPLICIT NONE
    
    INTEGER::I,J,JJ,K,II,TOTAL_ELEM_NODES,P_INT_MEF,MAX_ORDER
    
    INTEGER,ALLOCATABLE,DIMENSION(:)::contelem_pnod
    
    REAL*8::AUX1,AUX2,AUX3,AUX4,TANG(3),TAN_NORM,QSI_SOURCE,JACOBIAN_SIDE,QSI_ELEMENT,X,Y,Z
    
    REAL*8,ALLOCATABLE,DIMENSION(:)::F,FF,FFLOCAL,SOLUCAO1,PHI_MEF,SOLUCAO,SOLUCAO2,SOLUCAO3,E_NORMAL_ELEM,D_PHI_MEF
    
    REAL*8,ALLOCATABLE,DIMENSION(:,:)::A,B,VALORES,MR,QSIW_MEF
    
    !CALL OMP_SET_NUM_THREADS(24)
    
    ! Initiating variable e_normal
    TOTAL_ELEM_NODES=0
    MAX_ORDER=0
    DO I=1,N_ELEMENTOS_MEF
        TOTAL_ELEM_NODES=TOTAL_ELEM_NODES+ORDEM_MEF(I)+1
        IF (ORDEM_MEF(I).GT.MAX_ORDER) MAX_ORDER=ORDEM_MEF(I)
    END DO
    
    !ALLOCATE(E_NORMAL(TOTAL_ELEM_NODES,4))
    IF (.NOT.ALLOCATED(E_NORMAL)) ALLOCATE(E_NORMAL_NOD(N_NOS_MEF),E_NORMAL(N_ELEMENTOS_MEF,MAX_ORDER+1))
    E_NORMAL=0.0D0
    
    ! _______________________ FINDING VALUES OF STRESS ________________________
    
    !!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(COORDNOS_MEF,CONECTI_MEF,ORDEM_MEF,N_ELEMENTOS_MEF,N_NOS_MEF,QSI_MEF,&
    !E_NORMAL,KKGLOBAL,GGGLOBAL,U_MEF,P_MEF,COORDPCOLOC_MEF,F_DES_ANT)
    !! AUMENTOU F_DES_ANT AQUI !!!!!
    !
    !!$omp do schedule(dynamic)
    DO I=1,N_ELEMENTOS_MEF
        II=ORDEM_MEF(I)+1
        ALLOCATE(A(3*II,3*II),F(3*II),SOLUCAO(3*II),B(3*II,3*II),FF(3*II),FFLOCAL(II),SOLUCAO1(3*II),&
            PHI_MEF(II),VALORES(II,3),MR(II,3*II),SOLUCAO2(II),SOLUCAO3(II),E_NORMAL_ELEM(II),D_PHI_MEF(II))
        A=0.D0
        F=0.D0
        SOLUCAO=0.D0
        B=0.D0
        FF=0.D0
        SOLUCAO1=0.D0
        MR=0.D0
            
        A(:,:)=KKGLOBAL(I,:,:)
        B(:,:)=GGGLOBAL(I,:,:)
        DO J=1,(ORDEM_MEF(I)+1)
            K=CONECTI_MEF(I,J)
            
            IF (SLIP_REINF) THEN
                F(3*J-2)=U_MEF(3*K-2)-S_ACUM(3*K-2)
                F(3*J-1)=U_MEF(3*K-1)-S_ACUM(3*K-1)   
                F(3*J)=U_MEF(3*K)-S_ACUM(3*K)
            ELSE
                F(3*J-2)=U_MEF(3*K-2)
                F(3*J-1)=U_MEF(3*K-1)   
                F(3*J)=U_MEF(3*K)
            ENDIF
            
            FF(3*J-2)=-P_MEF(3*K-2)
            FF(3*J-1)=-P_MEF(3*K-1)    
            FF(3*J)=-P_MEF(3*K)
        ENDDO
        CALL DMURRV(3*II,3*II,A,3*II,3*II,F,1,3*II,SOLUCAO)
        CALL DMURRV(3*II,3*II,B,3*II,3*II,FF,1,3*II,SOLUCAO1)

        DO JJ=1,(ORDEM_MEF(I)+1)
            !VALORES(JJ,1)=COORDNOS_MEF(CONECTI_MEF(I,JJ),1)
            !VALORES(JJ,2)=COORDNOS_MEF(CONECTI_MEF(I,JJ),2)
            !VALORES(JJ,3)=COORDNOS_MEF(CONECTI_MEF(I,JJ),3)
            VALORES(JJ,1)=COORDPCOLOC_MEF(CONECTI_MEF(I,JJ),1)
            VALORES(JJ,2)=COORDPCOLOC_MEF(CONECTI_MEF(I,JJ),2)
            VALORES(JJ,3)=COORDPCOLOC_MEF(CONECTI_MEF(I,JJ),3)
        ENDDO
            
        DO J=1,(ORDEM_MEF(I)+1)
            D_PHI_MEF=0.D0
                
            CALL FUNCOES_DE_FORMA_MEC_MEF(3,QSI_MEF(I,J),(ORDEM_MEF(I)+1),I,VALORES,X,Y,Z,D_PHI_MEF)
            
        	AUX1=0.D0
	        AUX2=0.D0
            AUX3=0.0D0
	        DO K=1,ORDEM_MEF(I)+1
		        AUX1=AUX1+D_PHI_MEF(K)*VALORES(K,1)
		        AUX2=AUX2+D_PHI_MEF(K)*VALORES(K,2)
                AUX3=AUX3+D_PHI_MEF(K)*VALORES(K,3)
        	ENDDO
            AUX4=DSQRT(AUX1*AUX1+AUX2*AUX2+AUX3*AUX3)
            AUX1=AUX1/AUX4
            AUX2=AUX2/AUX4
            AUX3=AUX3/AUX4
            !______________TO AVOID NUMERICAL ERRORS:            
            IF (ABS(AUX1).LT.TOLER) AUX1=0.D0
            IF (ABS(AUX2).LT.TOLER) AUX2=0.D0
            IF (ABS(AUX3).LT.TOLER) AUX3=0.D0
            !______________
        
            MR(J,3*J-2) = AUX1
            MR(J,3*J-1) = AUX2
            MR(J,3*J) = AUX3
        ENDDO
        
        
        SOLUCAO2=0.D0
        SOLUCAO3=0.D0
!            SOLUCAO2=MATMUL(TRANSPOSE(MR),SOLUCAO)
        CALL DMURRV(II,3*II,MR,II,3*II,SOLUCAO,1,II,SOLUCAO2)
!            SOLUCAO3=MATMUL(TRANSPOSE(MR),SOLUCAO1)
        CALL DMURRV(II,3*II,MR,II,3*II,SOLUCAO1,1,II,SOLUCAO3)
   
        CALL DMURRV(II,3*II,MR,II,3*II,FF,1,II,FFLOCAL)
        E_NORMAL_ELEM=SOLUCAO2-SOLUCAO3
        
        ! MUDOU ------------
        IF (EP_ANALYSIS_REINF) THEN
            DO J=1,ORDEM_MEF(I)+1
                K=CONECTI_MEF(I,J)
            
                SOLUCAO1(3*J-2)=F_DES_ANT(3*K-2)
                SOLUCAO1(3*J-1)=F_DES_ANT(3*K-1)    
                SOLUCAO1(3*J)=F_DES_ANT(3*K)
            ENDDO
            SOLUCAO2=MATMUL(MR,SOLUCAO1)
            
            E_NORMAL_ELEM=E_NORMAL_ELEM-SOLUCAO2
        ENDIF
        ! FIM DE MUDOU -----
        
        !WRITE(*,*)E_NORMAL_ELEM
        !READ(*,*)
        
        P_INT_MEF=CEILING((ORDEM_MEF(I)+1)/2.0D0)     
        ALLOCATE(QSIW_MEF(P_INT_MEF,2))
        CALL GAUSS_POINTS(P_INT_MEF,QSIW_MEF)
            
        DO K=1,(ORDEM_MEF(I)+1)
                
            !E_NORMAL(KK,1)=I
            !E_NORMAL(KK,2)=CONECTI_MEF(I,K)
            IF ((K.eq.1).or.(K.eq.(ORDEM_MEF(I)+1))) then
                
                IF (K.EQ.1) E_NORMAL_ELEM(1)=-E_NORMAL_ELEM(1)
                    
                !E_NORMAL(KK,3)=DABS(E_NORMAL_ELEM(2*K-1))
                E_NORMAL(I,K)=E_NORMAL_ELEM(K)

            ELSE
                
                QSI_SOURCE=QSI_MEF(I,K)              
                JACOBIAN_SIDE=(qsi_source+1.0d0)/2.0d0
                AUX1=0.D0
                DO II=1,(ORDEM_MEF(I)+1)
                    DO J=1,P_INT_MEF
                        
                        QSI_ELEMENT=(QSIW_MEF(J,1)+1.0D0)*JACOBIAN_SIDE-1.0d0

                        CALL FUNCOES_DE_FORMA_MEC_MEF(2,QSI_ELEMENT,(ORDEM_MEF(I)+1),I,VALORES,X,Y,Z,PHI_MEF)
                                                        
                        CALL FUNCOES_DE_FORMA_MEC_MEF(3,QSI_ELEMENT,(ORDEM_MEF(I)+1),I,VALORES,X,Y,Z,D_PHI_MEF)
!
                        TANG=0.D0
                        DO JJ=1,(ORDEM_MEF(I)+1)
	                     !   TANG(1)=TANG(1)+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(I,JJ),1) !TANGENT VECTORS
		                    !TANG(2)=TANG(2)+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(I,JJ),2) 
                      !      TANG(3)=TANG(3)+D_PHI_MEF(JJ)*COORDNOS_MEF(CONECTI_MEF(I,JJ),3)  
                            TANG(1)=TANG(1)+D_PHI_MEF(JJ)*COORDPCOLOC_MEF(CONECTI_MEF(I,JJ),1) !TANGENT VECTORS
		                    TANG(2)=TANG(2)+D_PHI_MEF(JJ)*COORDPCOLOC_MEF(CONECTI_MEF(I,JJ),2) 
                            TANG(3)=TANG(3)+D_PHI_MEF(JJ)*COORDPCOLOC_MEF(CONECTI_MEF(I,JJ),3) 
        	            ENDDO
        	            TAN_NORM=DSQRT(TANG(1)**2.0+TANG(2)**2.0+TANG(3)**2.0) !TANGENT VECTOR NORM, JACOBIAN OF THE POINT
            	                
                        !ROTATE FF VECTOR:
                        AUX1=AUX1+PHI_MEF(II)*JACOBIAN_SIDE*QSIW_MEF(J,2)*TAN_NORM*FFLOCAL(II)
                    ENDDO
                ENDDO

                !E_NORMAL(KK,3)=E_NORMAL_ELEM(1)-AUX1
                E_NORMAL(I,K)=E_NORMAL_ELEM(1)-AUX1
     
                !E_NORMAL(KK,3)=DABS(E_NORMAL(KK,3))
            ENDIF
!                PAUSE
        ENDDO
        DEALLOCATE(QSIW_MEF)
!_____________________________________________________________________________            

        !WRITE(*,*)(E_NORMAL(I,J),J=1,ORDEM_MEF(I)+1)
        !READ(*,*)
   
        DEALLOCATE(A,F,SOLUCAO,B,FF,FFLOCAL,SOLUCAO1,PHI_MEF,VALORES,MR,SOLUCAO2,SOLUCAO3,D_PHI_MEF,E_NORMAL_ELEM)       
    ENDDO  
    !!$omp end do
    !!$OMP END PARALLEL 
    
    ! CALCULANDO TENSAO POR NÓ
    IF (.NOT.EP_ANALYSIS_REINF) THEN
        ALLOCATE(contelem_pnod(N_NOS_MEF))
        E_NORMAL_NOD=0.0D0
        contelem_pnod=0
        DO I=1,N_ELEMENTOS_MEF
            DO J=1,ORDEM_MEF(I)+1
                E_NORMAL_NOD(CONECTI_MEF(I,J))=E_NORMAL_NOD(CONECTI_MEF(I,J))+E_NORMAL(I,J)
                contelem_pnod(CONECTI_MEF(I,J))=contelem_pnod(CONECTI_MEF(I,J))+1            
            ENDDO
        ENDDO
        DO I=1,N_NOS_MEF
            E_NORMAL_NOD(I)=E_NORMAL_NOD(I)/DFLOAT(contelem_pnod(I))
        ENDDO
        DEALLOCATE(contelem_pnod)
    ENDIF
    
    
END SUBROUTINE