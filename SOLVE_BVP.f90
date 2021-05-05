	SUBROUTINE SOLVE_BVP
!
    !INCLUDE MKL.FI
	USE ISOPARAMETRIC_MESH
	USE CRACKS_DUAL_BEM
    USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
    USE PROPAGATION
    USE XBEM_FORCE_VARIABLES
    USE XBEM_SUPPORT_VARIABLES
    USE XBEM_CRACKFRONT_VARIABLES
!
	IMPLICIT NONE 
!
	INTEGER::I,J,K,II,JJ,KK,CONT,NEN,IPIV[ALLOCATABLE](:),INFO,DIM_XBEM,DIM_XBEM_CRACKFRONT,AUX,IDO,IRANK
!
	REAL*8::A[ALLOCATABLE](:,:),B[ALLOCATABLE](:,:),F[ALLOCATABLE](:),SOLUTION[ALLOCATABLE](:),&
    NCOND,RCOND,IPVT_A[ALLOCATABLE](:),ANORM,WORK[ALLOCATABLE](:),DLANGE,WORK2[ALLOCATABLE](:),&
    A_ALT[ALLOCATABLE](:,:),P_TEST[ALLOCATABLE](:),AUX_SOL[ALLOCATABLE](:),Z_TEST[ALLOCATABLE](:),A_RES[ALLOCATABLE](:,:)
!
	II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
    DIM_XBEM=3*N_COLLOCPOINTS+NUM_CON_SUPPS+NODES_DIST_SUPPS
!
	ALLOCATE(A(DIM_XBEM,DIM_XBEM),SOLUTION(DIM_XBEM),B(3*N_COLLOCPOINTS,II),F(II))
    IF(ALLOCATED(IPVT_A))DEALLOCATE(IPVT_A,WORK)
    ALLOCATE(IPVT_A(DIM_XBEM),WORK(DIM_XBEM),WORK2(4*DIM_XBEM))
!
	A=0.D0
	B=0.D0
	SOLUTION=0.D0
	F=0.D0
	CONT=0
!
	DO I=1,N_COLLOCPOINTS
		II=0
		JJ=0
		IF(DUAL_BEM(I).EQ.'B')THEN
		    DO J=1,SIZE(NI,1)
			    IF(I.EQ.NI(J,1)) THEN
				    II=-1
				    JJ=NI(J,2)
			    ELSE
				    IF(I.EQ.NI(J,2)) THEN
					    II=NI(J,1)
					    JJ=-1
				    ENDIF
			    ENDIF
			ENDDO        
		    IF(II.EQ.-1) THEN
			    A(1:3*N_COLLOCPOINTS,(3*I-2))=A(1:3*N_COLLOCPOINTS,(3*I-2))+H(:,(3*I-2))
			    A(1:3*N_COLLOCPOINTS,3*I-1)=A(1:3*N_COLLOCPOINTS,3*I-1)+H(:,3*I-1)
			    A(1:3*N_COLLOCPOINTS,3*I)=A(1:3*N_COLLOCPOINTS,3*I)+H(:,3*I)
!
			    A(1:3*N_COLLOCPOINTS,(3*I-2))=A(1:3*N_COLLOCPOINTS,(3*I-2))+H(:,(3*JJ-2))
			    A(1:3*N_COLLOCPOINTS,3*I-1)=A(1:3*N_COLLOCPOINTS,3*I-1)+H(:,3*JJ-1)
			    A(1:3*N_COLLOCPOINTS,3*I)=A(1:3*N_COLLOCPOINTS,3*I)+H(:,3*JJ)
		    ELSE
			    IF(JJ.EQ.-1) THEN
				    A(1:3*N_COLLOCPOINTS,(3*I-2))=A(1:3*N_COLLOCPOINTS,(3*I-2))+G(:,(3*I-2))
				    A(1:3*N_COLLOCPOINTS,3*I-1)=A(1:3*N_COLLOCPOINTS,3*I-1)+G(:,3*I-1)
				    A(1:3*N_COLLOCPOINTS,3*I)=A(1:3*N_COLLOCPOINTS,3*I)+G(:,3*I)
!
				    A(1:3*N_COLLOCPOINTS,(3*I-2))=A(1:3*N_COLLOCPOINTS,(3*I-2))-G(:,(3*II-2))
				    A(1:3*N_COLLOCPOINTS,3*I-1)=A(1:3*N_COLLOCPOINTS,3*I-1)-G(:,3*II-1)
				    A(1:3*N_COLLOCPOINTS,3*I)=A(1:3*N_COLLOCPOINTS,3*I)-G(:,3*II)
			    ELSE
				    IF(B_CONDITIONS(3*I-2).EQ.0) THEN
					    A(1:3*N_COLLOCPOINTS,(3*I-2))=-G(:,(3*I-2))
!
					    CONT=CONT+1
					    B(:,CONT)=-H(:,(3*I-2))
					    F(CONT)=U((3*I-2))
				    ELSE
					    A(1:3*N_COLLOCPOINTS,(3*I-2))=H(:,(3*I-2))	
!
					    CONT=CONT+1
					    B(:,CONT)=G(:,(3*I-2))
					    F(CONT)=T((3*I-2))
				    ENDIF
!
				    IF(B_CONDITIONS(3*I-1).EQ.0) THEN
					    A(1:3*N_COLLOCPOINTS,3*I-1)=-G(:,3*I-1)
!
					    CONT=CONT+1
					    B(:,CONT)=-H(:,3*I-1)
					    F(CONT)=U(3*I-1)
				    ELSE
					    A(1:3*N_COLLOCPOINTS,3*I-1)=H(:,3*I-1)	
!
					    CONT=CONT+1
					    B(:,CONT)=G(:,3*I-1)
					    F(CONT)=T(3*I-1)
				    ENDIF
!
				    IF(B_CONDITIONS(3*I).EQ.0) THEN
					    A(1:3*N_COLLOCPOINTS,3*I)=-G(:,3*I)
!
					    CONT=CONT+1
				    	B(:,CONT)=-H(:,3*I)
					    F(CONT)=U(3*I)
				    ELSE
					    A(1:3*N_COLLOCPOINTS,3*I)=H(:,3*I)	
!
					    CONT=CONT+1
					    B(:,CONT)=G(:,3*I)
					    F(CONT)=T(3*I)
				    ENDIF
			    ENDIF
		    ENDIF
		ELSE
		    DO J=1,N_CRACK_POINTS
			    IF(I.EQ.NC(J,1)) THEN
				    II=-1
				    JJ=NC(J,2)
			    ELSE
				    IF(I.EQ.NC(J,2)) THEN
					    II=NC(J,1)
					    JJ=-1
				    ENDIF
			    ENDIF
			ENDDO        
		    IF(II.EQ.-1) THEN
		    	A(1:3*N_COLLOCPOINTS,(3*I-2))=A(1:3*N_COLLOCPOINTS,(3*I-2))+H(:,(3*I-2))
			    A(1:3*N_COLLOCPOINTS,(3*I-2))=A(1:3*N_COLLOCPOINTS,(3*I-2))-H(:,(3*JJ-2))			    		    
!			    
		        IF(B_CONDITIONS(3*I-1).EQ.0) THEN
			        A(1:3*N_COLLOCPOINTS,(3*I-1))=-G(:,(3*I-1))
!
					CONT=CONT+1
					B(:,CONT)=-H(:,(3*I-1))
					F(CONT)=U((3*I-1))
				ELSE
					A(1:3*N_COLLOCPOINTS,(3*I-1))=H(:,(3*I-1))	
!
					CONT=CONT+1
					B(:,CONT)=G(:,(3*I-1))
					F(CONT)=T((3*I-1))
				ENDIF
!
				IF(B_CONDITIONS(3*I).EQ.0) THEN
					A(1:3*N_COLLOCPOINTS,3*I)=-G(:,3*I)
!
					CONT=CONT+1
					B(:,CONT)=-H(:,3*I)
					F(CONT)=U(3*I)
				ELSE
					A(1:3*N_COLLOCPOINTS,3*I)=H(:,3*I)	
!
					CONT=CONT+1
					B(:,CONT)=G(:,3*I)
					F(CONT)=T(3*I)
				ENDIF
!		    
		    ELSE
			    IF(JJ.EQ.-1) THEN
				    A(1:3*N_COLLOCPOINTS,3*I-2)=A(1:3*N_COLLOCPOINTS,3*I-2)-G(:,3*I-2)
				    A(1:3*N_COLLOCPOINTS,3*I-2)=A(1:3*N_COLLOCPOINTS,3*I-2)-G(:,3*II-2)
!				    			    		    
		            IF(B_CONDITIONS(3*I-1).EQ.0) THEN
			            A(1:3*N_COLLOCPOINTS,(3*I-1))=-G(:,(3*I-1))
!
					    CONT=CONT+1
					    B(:,CONT)=-H(:,(3*I-1))
				 	    F(CONT)=U((3*I-1))
				    ELSE
					    A(1:3*N_COLLOCPOINTS,(3*I-1))=H(:,(3*I-1))	
!
					    CONT=CONT+1
					    B(:,CONT)=G(:,(3*I-1))
					    F(CONT)=T((3*I-1))
				    ENDIF
!
				    IF(B_CONDITIONS(3*I).EQ.0) THEN
					    A(1:3*N_COLLOCPOINTS,3*I)=-G(:,3*I)
!
					    CONT=CONT+1
					    B(:,CONT)=-H(:,3*I)
					    F(CONT)=U(3*I)
				    ELSE
					    A(1:3*N_COLLOCPOINTS,3*I)=H(:,3*I)	
!
					    CONT=CONT+1
					    B(:,CONT)=G(:,3*I)
					    F(CONT)=T(3*I)
				    ENDIF			    
!			    
			    ELSE
				    IF(B_CONDITIONS(3*I-2).EQ.0) THEN
					    A(1:3*N_COLLOCPOINTS,(3*I-2))=-G(:,(3*I-2))
!
					    CONT=CONT+1
					    B(:,CONT)=-H(:,(3*I-2))
					    F(CONT)=U((3*I-2))
				    ELSE
					    A(1:3*N_COLLOCPOINTS,(3*I-2))=H(:,(3*I-2))	
!
					    CONT=CONT+1
					    B(:,CONT)=G(:,(3*I-2))
					    F(CONT)=T((3*I-2))
				    ENDIF
!
				    IF(B_CONDITIONS(3*I-1).EQ.0) THEN
					    A(1:3*N_COLLOCPOINTS,3*I-1)=-G(:,3*I-1)
!
					    CONT=CONT+1
					    B(:,CONT)=-H(:,3*I-1)
					    F(CONT)=U(3*I-1)
				    ELSE
					    A(1:3*N_COLLOCPOINTS,3*I-1)=H(:,3*I-1)	
!
					    CONT=CONT+1
					    B(:,CONT)=G(:,3*I-1)
					    F(CONT)=T(3*I-1)
				    ENDIF
!
				    IF(B_CONDITIONS(3*I).EQ.0) THEN
					    A(1:3*N_COLLOCPOINTS,3*I)=-G(:,3*I)
!
					    CONT=CONT+1
				    	B(:,CONT)=-H(:,3*I)
					    F(CONT)=U(3*I)
				    ELSE
					    A(1:3*N_COLLOCPOINTS,3*I)=H(:,3*I)	
!
					    CONT=CONT+1
					    B(:,CONT)=G(:,3*I)
					    F(CONT)=T(3*I)
				    ENDIF
			    ENDIF
		    ENDIF
		ENDIF
    ENDDO
!
    II=3*N_COLLOCPOINTS-2*(3*N_INTERFACE_POINTS+N_CRACK_POINTS)
	CALL DMURRV(3*N_COLLOCPOINTS,II,B,3*N_COLLOCPOINTS,II,F,1,3*N_COLLOCPOINTS,SOLUTION)
	DEALLOCATE(F)
	ALLOCATE(F(DIM_XBEM))
    F=0.D0
	F(1:3*N_COLLOCPOINTS)=SOLUTION(1:3*N_COLLOCPOINTS)
	SOLUTION=0.D0
!
    IF ((NUM_CON_SUPPS+NODES_DIST_SUPPS).GT.0) THEN
        !DO I=1,NUM_CON_SUPPS
            A(1:3*N_COLLOCPOINTS,3*N_COLLOCPOINTS+1:3*N_COLLOCPOINTS+NUM_CON_SUPPS+NODES_DIST_SUPPS)=AS_COL
            A(3*N_COLLOCPOINTS+1:3*N_COLLOCPOINTS+NUM_CON_SUPPS+NODES_DIST_SUPPS,1:3*N_COLLOCPOINTS)=AS_LIN
            F(3*N_COLLOCPOINTS+1:3*N_COLLOCPOINTS+NUM_CON_SUPPS)=VALUES_CON_SUPPS
            F(3*N_COLLOCPOINTS+NUM_CON_SUPPS+1:3*N_COLLOCPOINTS+NUM_CON_SUPPS+NODES_DIST_SUPPS)=VALUES_DIST_SUPPS
        !ENDDO
    ENDIF  
    !AUX=3*N_COLLOCPOINTS+NUM_CON_SUPPS
    !IF (NODES_DIST_SUPPS.GT.0) THEN
    !    A(1:3*N_COLLOCPOINTS,AUX+1:AUX+NODES_DIST_SUPPS)=AS_COL!(:,1:NODES_DIST_SUPPS)
    !    A(AUX+1:AUX+NODES_DIST_SUPPS,1:3*N_COLLOCPOINTS)=AS_LIN!(1:NODES_DIST_SUPPS,:)
    !    F(AUX+1:AUX+NODES_DIST_SUPPS)=VALUES_DIST_SUPPS(:)
    !ENDIF  
!
!	TIRAR DAQUI DEPOIS
    IF ((NUM_CON_LOADS+NUM_DIST_LOADS).GT.0) THEN
        WRITE(*,*)'XBEM_FORCE'
        ALLOCATE(F_ENRICHED(3*N_COLLOCPOINTS))
        F_ENRICHED=0.D0
        CALL XBEM_FORCE
        IF(NUM_DIST_LOADS.GT.0) CALL XBEM_DIST_FORCE
        F(1:3*N_COLLOCPOINTS)=F(1:3*N_COLLOCPOINTS)+F_ENRICHED
        DEALLOCATE(F_ENRICHED)
    ENDIF
!
    !IF (SITEXTRACTION_METHOD.EQ."XBEM") THEN
    !    N_K=2
    !    CALL XBEM_K_AUXILIAR_ELEMS
    !    DIM_XBEM=DIM_XBEM+3*N_CRACKTIPS
    !    ALLOCATE(A_AUX(3*N_COLLOCPOINTS,3*N_COLLOCPOINTS),F_AUX(3*N_COLLOCPOINTS))
    !    A_AUX=A
    !    F_AUX=F
    !    DEALLOCATE(A,F,SOLUTION)
    !    ALLOCATE(A(DIM_XBEM,DIM_XBEM),F(DIM_XBEM),SOLUTION(DIM_XBEM))
    !    A=0.D0
    !    F=0.D0
    !    A(1:3*N_COLLOCPOINTS,1:3*N_COLLOCPOINTS)=A_AUX
    !    F(1:3*N_COLLOCPOINTS)=F_AUX
    !    CALL XBEM_CRACK_TIP(A,DIM_XBEM)
    !ENDIF
!
!   WRITE(*,*)"SOLVING THE PROBLEM"
!
	ALLOCATE(IPIV(DIM_XBEM))
    !OPEN(6,file='F.txt',status='unknown')
    !DO I=1,DIM_XBEM
    !    WRITE(6,'(1000F)')F(I)
    !ENDDO
    !CLOSE(6)
    !
    !OPEN(6,file='MATRIZ A.txt',status='unknown')
    !DO I=1,DIM_XBEM
    !    WRITE(6,'(10000F)')(A(I,K),K=1,DIM_XBEM)
    !ENDDO

!   SOLUCAO TRADICIONAL
    DIM_XBEM=SIZE(F)
    SOLUTION=F
!
!    ALLOCATE(A_RES(DIM_XBEM,DIM_XBEM),A_ALT(DIM_XBEM,DIM_XBEM))
!!   PRE MULTIPLICANDO PELA TRANSPOSTA
!   ! DTRNRR TRANSPÕE
!    CALL DTRNRR (DIM_XBEM,DIM_XBEM,A,DIM_XBEM,DIM_XBEM,DIM_XBEM,A_ALT,DIM_XBEM)
!   ! DMURRV FAZ MATRIZ VETOR
!    CALL DMURRV (DIM_XBEM,DIM_XBEM,A_ALT,DIM_XBEM,DIM_XBEM,F,1,DIM_XBEM,SOLUTION)
!   ! DMXTXF FAZ A^T*A
!    CALL DMXTXF (DIM_XBEM,DIM_XBEM,A,DIM_XBEM,DIM_XBEM,A_RES,DIM_XBEM)
!!   !MEDINDO NUMERO DE CONDICAO tentando precondiconante
!    RCOND=0.d0
!    CALL DGETRF(DIM_XBEM,DIM_XBEM,A_RES,DIM_XBEM,IPVT_A,INFO)
!    ANORM = DLANGE('1',DIM_XBEM,DIM_XBEM,A_RES,DIM_XBEM, WORK)
!    CALL DGECON('1',DIM_XBEM,A_RES,DIM_XBEM,ANORM,RCOND,WORK2,DIM_XBEM,INFO)
!    CALL DGETRS('N',DIM_XBEM,1,A_RES,DIM_XBEM,IPVT_A,SOLUTION,DIM_XBEM,INFO)
!    NCOND=1.D0/RCOND
!   PRÉ-MULTIPLICANDO PELA DIAGONAL
    !A_ALT=0.D0
    !CALL DLSGRR (DIM_XBEM,DIM_XBEM,A,DIM_XBEM,1.d-8,IRANK,A_ALT,DIM_XBEM)
    !! A_ALT*F
    !CALL DMURRV (DIM_XBEM,DIM_XBEM,A_ALT,DIM_XBEM,DIM_XBEM,F,1,DIM_XBEM,SOLUTION)
    !! A_ALT*A
    !CALL DMRRRR (DIM_XBEM,DIM_XBEM,A_ALT,DIM_XBEM,DIM_XBEM,DIM_XBEM,A,DIM_XBEM,DIM_XBEM,DIM_XBEM,A_RES,DIM_XBEM)
!  
!   MEDINDO NUMERO DE CONDICAO
    !RCOND=0.d0
    !CALL DGETRF(DIM_XBEM,DIM_XBEM,A,DIM_XBEM,IPVT_A,INFO)
    !ANORM = DLANGE('1',DIM_XBEM,DIM_XBEM,A,DIM_XBEM, WORK)
    !CALL DGECON('1',DIM_XBEM,A,DIM_XBEM,ANORM,RCOND,WORK2,DIM_XBEM,INFO)
    !CALL DGETRS('N',DIM_XBEM,1,A,DIM_XBEM,IPVT_A,SOLUTION,DIM_XBEM,INFO)
    !NCOND=1.D0/RCOND
!   
!   SOLUCAO TRADICIONAL
    CALL DGESV(DIM_XBEM,1,A,DIM_XBEM,IPIV,SOLUTION,DIM_XBEM,INFO)
!
    !WRITE(6,*)
    !WRITE(*,*)'INFO=',INFO
    !WRITE(*,*)'RCOND=',RCOND
    !CLOSE(6)
!
    IF((NUM_CON_SUPPS+NUM_DIST_SUPPS).GT.0)THEN
        XBEM_SUPP_REACTIONS=SOLUTION(3*N_COLLOCPOINTS+1:DIM_XBEM)
        OPEN(2,FILE="Output_data\XBEM_SUPP_REACTIONS.TXT",STATUS='UNKNOWN')
        DO I=1,NUM_CON_SUPPS+NODES_DIST_SUPPS
            WRITE(2,*)XBEM_SUPP_REACTIONS(I)
        ENDDO
        CLOSE(2,STATUS='KEEP')
    ENDIF
!
200 FORMAT(4x,i5,4x,E13.6,4x,E13.6,4x,E13.6) 
    !IF(SITEXTRACTION_METHOD.EQ."XBEM")THEN
    !    OPEN(3,FILE="STRESS_INTENSITY_FACTOR_XBEM.TXT",STATUS='UNKNOWN')
    !    WRITE(3,*)"CRACK TIP    KI     		   KII			    KIII"
    !    DO I=1,N_CRACKTIPS
    !        WRITE(3,200)CRACKTIP_COLLOCPOINTS(I),F(3*N_COLLOCPOINTS+3*I-2),F(3*N_COLLOCPOINTS+3*I-1),F(3*N_COLLOCPOINTS+3*I)
    !    ENDDO
    !    CLOSE(3,STATUS='KEEP')
    !ENDIF
!
	DO I=1,N_COLLOCPOINTS
		II=0
		JJ=0
		IF(DUAL_BEM(I).EQ.'B')THEN
		    DO J=1,SIZE(NI,1)
			    IF(I.EQ.NI(J,1)) THEN
				    II=-1
				    JJ=NI(J,2)
			    ELSE
				    IF(I.EQ.NI(J,2)) THEN
					    II=NI(J,1)
					    JJ=-1
				    ENDIF
			    ENDIF
		    ENDDO	
!
		    IF(II.EQ.-1) THEN
			    U(3*I-2)=SOLUTION(3*I-2)
			    U(3*I-1)=SOLUTION(3*I-1)
			    U(3*I)=SOLUTION(3*I)
!
			    U(3*JJ-2)=SOLUTION(3*I-2)
			    U(3*JJ-1)=SOLUTION(3*I-1)
			    U(3*JJ)=SOLUTION(3*I)
		    ELSE
			    IF(JJ.EQ.-1) THEN
				    T(3*I-2)=-SOLUTION(3*I-2)
				    T(3*I-1)=-SOLUTION(3*I-1)
			        T(3*I)=-SOLUTION(3*I)
!
				    T(3*II-2)=SOLUTION(3*I-2)
				    T(3*II-1)=SOLUTION(3*I-1)
				    T(3*II)=SOLUTION(3*I)
			    ELSE
				    IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					    T(3*I-2)=SOLUTION(3*I-2)
				    ELSE
					    U(3*I-2)=SOLUTION(3*I-2)
				    ENDIF
!
				    IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					    T(3*I-1)=SOLUTION(3*I-1)
				    ELSE
					    U(3*I-1)=SOLUTION(3*I-1)
				    ENDIF
!
				    IF (B_CONDITIONS(3*I).EQ.0) THEN
					    T(3*I)=SOLUTION(3*I)
				    ELSE
					    U(3*I)=SOLUTION(3*I)
				    ENDIF
			    ENDIF
		    ENDIF
		ELSE
		    DO J=1,N_CRACK_POINTS
			    IF(I.EQ.NC(J,1)) THEN
				    II=-1
				    JJ=NC(J,2)
			    ELSE
				    IF(I.EQ.NC(J,2)) THEN
					    II=NC(J,1)
					    JJ=-1
				    ENDIF
			    ENDIF
		    ENDDO	
!
		    IF(II.EQ.-1) THEN
			    U(3*I-2)=SOLUTION(3*I-2)
			    U(3*JJ-2)=-SOLUTION(3*I-2)
!			    		    
			    IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					T(3*I-1)=SOLUTION(3*I-1)
				ELSE
					U(3*I-1)=SOLUTION(3*I-1)
				ENDIF
!
				IF (B_CONDITIONS(3*I).EQ.0) THEN
					T(3*I)=SOLUTION(3*I)
				ELSE
					U(3*I)=SOLUTION(3*I)
				ENDIF
!				    		    
		    ELSE
			    IF(JJ.EQ.-1) THEN
			        T(3*I-2)=SOLUTION(3*I-2)
				    T(3*II-2)=SOLUTION(3*I-2)
!				    			    
			        IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					    T(3*I-1)=SOLUTION(3*I-1)
				    ELSE
					    U(3*I-1)=SOLUTION(3*I-1)
				    ENDIF
!
				    IF (B_CONDITIONS(3*I).EQ.0) THEN
					    T(3*I)=SOLUTION(3*I)
				    ELSE
					    U(3*I)=SOLUTION(3*I)
				    ENDIF			    
!			    
			    ELSE
				    IF (B_CONDITIONS(3*I-2).EQ.0) THEN
					    T(3*I-2)=SOLUTION(3*I-2)
				    ELSE
					    U(3*I-2)=SOLUTION(3*I-2)
				    ENDIF
!
				    IF (B_CONDITIONS(3*I-1).EQ.0) THEN
					    T(3*I-1)=SOLUTION(3*I-1)
				    ELSE
					    U(3*I-1)=SOLUTION(3*I-1)
				    ENDIF
!
				    IF (B_CONDITIONS(3*I).EQ.0) THEN
					    T(3*I)=SOLUTION(3*I)
				    ELSE
					    U(3*I)=SOLUTION(3*I)
				    ENDIF
			    ENDIF
		    ENDIF			
		ENDIF
	ENDDO  	    
!    
	END SUBROUTINE