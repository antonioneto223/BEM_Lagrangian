	SUBROUTINE SOLVE_ELASTIC_CRACK_PROBLEM
!	
    USE ISOPARAMETRIC_MESH
    USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS	
	USE FATIGUE   	  
    USE VISCO
!
    IMPLICIT NONE   
    INTEGER::I,J,K,II,NEN,CONT,IT,N_CRACK_COLLOCPOINTS,KODE_CONSTRAINS[ALLOCATABLE](:),NCAUX[ALLOCATABLE](:,:)
!
    REAL*8::VALUES[ALLOCATABLE](:,:),DX,DY,DZ,DR,TOLCOD,TOLTn
    REAL*8::U_PRESC(3*N_COLLOCPOINTS),T_PRESC(3*N_COLLOCPOINTS),COD_AUX[ALLOCATABLE](:),Tn_AUX[ALLOCATABLE](:)
!
    LOGICAL::ADHESION_INTERPENETRATION[ALLOCATABLE](:)
!
    TOLCOD=-1.0D-16 
    TOLTn=1.0D-16        
    N_CRACK_COLLOCPOINTS=0
    N_CRACK_POINTS=0
    
    ALLOCATE(COD(N_COLLOCPOINTS),CSD(N_COLLOCPOINTS),CTD(N_COLLOCPOINTS))
       
    DO I=1,N_COLLOCPOINTS
        IF(DUAL_BEM(I).EQ.'S')THEN
            N_CRACK_COLLOCPOINTS=N_CRACK_COLLOCPOINTS+1
        ENDIF
    ENDDO
    IF(N_CRACK_COLLOCPOINTS.NE.0)THEN    
        ALLOCATE(KODE_CONSTRAINS(N_CRACK_COLLOCPOINTS),NCAUX(N_CRACK_COLLOCPOINTS,2),ADHESION_INTERPENETRATION(N_CRACK_COLLOCPOINTS))
        ALLOCATE(COD_AUX(N_CRACK_COLLOCPOINTS),Tn_AUX(N_CRACK_COLLOCPOINTS))
    ELSE
        IF (VISCO_ANALYSIS) THEN
            CALL SOLVE_BVP_VISCO
        ELSE
            CALL SOLVE_BVP  
        ENDIF
    ENDIF    
    NCAUX=NC 
    ADHESION_INTERPENETRATION=.TRUE.
    KODE_CONSTRAINS=0 
    U_PRESC=U
    T_PRESC=T
    IT=0
!    
    DO WHILE (ANY(ADHESION_INTERPENETRATION))
        IT=IT+1      
        WRITE(*,*)IT 
        U=U_PRESC
        T=T_PRESC
!        
        CALL SOLVE_BVP
!
        COD=0
        CSD=0
        CTD=0
        
        COD_AUX=0
        Tn_AUX=0
        CONT=0
        DO I=1,N_COLLOCPOINTS
            IF(DUAL_BEM(I).EQ.'S')THEN
                DO J=1,N_COLLOCPOINTS
                    IF(DUAL_BEM(J).EQ.'H')THEN
                        DX=COORD_COLLOCPOINTS(I,1)-COORD_COLLOCPOINTS(J,1)
				        DY=COORD_COLLOCPOINTS(I,2)-COORD_COLLOCPOINTS(J,2)
		                DZ=COORD_COLLOCPOINTS(I,3)-COORD_COLLOCPOINTS(J,3)
				        DR=DSQRT(DX*DX+DY*DY+DZ*DZ)
				        IF(DR.LE.TOL) THEN
				            II=J
                        ENDIF
                    ENDIF    
                ENDDO
                COD(I)=-U(3*I-2)-U(3*II-2)     
                CSD(I)=U(3*I-1)+U(3*II-1) 
                CTD(I)=U(3*I)-U(3*II) 
!                
                COD(II)=COD(I)   
                CSD(II)=CSD(I)     
                CTD(II)=CTD(I)
                
                CONT=CONT+1
                COD_AUX(CONT)=COD(I)   
                Tn_AUX(CONT)=T(3*I-2) 
            ENDIF     
        ENDDO 
!               
        DO I=1,N_CRACK_COLLOCPOINTS
            SELECT CASE(KODE_CONSTRAINS(I))
            CASE(0)
                IF(COD(NCAUX(I,1)).LT.TOLCOD)THEN
                    KODE_CONSTRAINS(I)=1
                    ADHESION_INTERPENETRATION(I)=.TRUE.
                ELSE
                    ADHESION_INTERPENETRATION(I)=.FALSE.                
                ENDIF
            CASE(1)
                IF(T(3*NCAUX(I,1)-2).GT.TOLTn)THEN
                    KODE_CONSTRAINS(I)=0
                    ADHESION_INTERPENETRATION(I)=.TRUE.
                ELSE
                    ADHESION_INTERPENETRATION(I)=.FALSE.                   
                ENDIF
            ENDSELECT             
        ENDDO
!            
        CONT=0        
        DO I=1,N_CRACK_COLLOCPOINTS
            IF(KODE_CONSTRAINS(I).EQ.1)THEN
                CONT=CONT+1        
            ENDIF        
        ENDDO
        IF(CONT.NE.0)THEN
            DEALLOCATE(NC)
            ALLOCATE(NC(CONT,2))
            CONT=0       
            DO I=1,N_CRACK_COLLOCPOINTS
                IF(KODE_CONSTRAINS(I).EQ.1)THEN
                    CONT=CONT+1
                    NC(CONT,:)=NCAUX(I,:)        
                ENDIF
            ENDDO
        ENDIF    
        N_CRACK_POINTS=CONT 
        ! forma bruta de desligar a interpenetração de faces de fissura
        ADHESION_INTERPENETRATION(:)=.FALSE.
!        
    ENDDO    
!
    OPEN(30,file='Output_data\DISPLACEMENT_DISCONTINUITIES.ogl',status='unknown')
    WRITE(30,*)'COLLOCPOINT        COD        CSD        CTD        Tn        Tt1        Tt2'
	WRITE(30,*)
	DO I=1,N_COLLOCPOINTS
	    IF(DUAL_BEM(I).EQ."B")THEN
	        WRITE(30,100) I,COD(I),CSD(I),CTD(I),0.D0,0.D0,0.D0
	    ELSE
	    	WRITE(30,100) I,COD(I),CSD(I),CTD(I),T(3*I-2),T(3*I-1),T(3*I)
	    ENDIF
	ENDDO
	CLOSE(30)
!
    DEALLOCATE(COD,CSD,CTD) 
!    		
    100	FORMAT(I5,8x,F14.9,8x,F14.9,8x,F14.9,8x,F14.9,8x,F14.9,8x,F14.9)    
!    
    END SUBROUTINE SOLVE_ELASTIC_CRACK_PROBLEM