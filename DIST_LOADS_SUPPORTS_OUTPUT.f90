	SUBROUTINE DIST_LOADS_SUPPORTS_OUTPUT
!
	USE ISOPARAMETRIC_MESH
	USE CRACKS_DUAL_BEM	
	USE ANALYSIS
	USE REMESHING
    USE SUB_REGIONS_INTERFACES
	USE XBEM_FORCE_VARIABLES
    USE XBEM_SUPPORT_VARIABLES
!
	INTEGER::I,J,K,II,JJ,KK,L,NEN,C,N_CRACKED_ELEM,CONT
!
    REAL*8::ETA_DIST_LOADS[ALLOCATABLE](:,:),VALUES(4,3),VJC(3),DVJC(3,2),JC,&
    ETA_DIST_SUPPS_X[ALLOCATABLE](:,:),ETA_DIST_SUPPS_Y[ALLOCATABLE](:,:),ETA_DIST_SUPPS_Z[ALLOCATABLE](:,:)
!
    CHARACTER(40)::NOME		
!
    N_CRACKED_ELEM=0
	DO I=1,N_ELEM
	    IF(CRACKED_ELEM(I).NE."UNCRACKED")THEN
	        N_CRACKED_ELEM=N_CRACKED_ELEM+1
	    ENDIF
    ENDDO
!
    ALLOCATE(ETA_DIST_LOADS(NUM_DIST_LOADS,3),ETA_DIST_SUPPS_X(NUM_DIST_SUPPS_X,3),&
    ETA_DIST_SUPPS_Y(NUM_DIST_SUPPS_Y,3),ETA_DIST_SUPPS_Z(NUM_DIST_SUPPS_Z,3))
!
    ! CONSIDERING ONLY PLANE ELEMENTS
    DO I=1,NUM_DIST_LOADS
        CALL NORMAL_OUTWARD(4,COORDS_DIST_LOADS(I,:,:),0.D0,0.D0,ETA_DIST_LOADS(I,:),VJC,DVJC,JC)
    ENDDO
    DO I=1,NUM_DIST_SUPPS_X
        DO J=1,4
            VALUES(J,1)=COORDS_DIST_SUPPS_X(CONEC_DIST_SUPPS_X(I,J),1)
            VALUES(J,2)=COORDS_DIST_SUPPS_X(CONEC_DIST_SUPPS_X(I,J),2)
            VALUES(J,3)=COORDS_DIST_SUPPS_X(CONEC_DIST_SUPPS_X(I,J),3)
        ENDDO
        CALL NORMAL_OUTWARD(4,VALUES,0.D0,0.D0,ETA_DIST_SUPPS_X(I,:),VJC,DVJC,JC)
    ENDDO
    DO I=1,NUM_DIST_SUPPS_Y
        DO J=1,4
            VALUES(J,1)=COORDS_DIST_SUPPS_Y(CONEC_DIST_SUPPS_Y(I,J),1)
            VALUES(J,2)=COORDS_DIST_SUPPS_Y(CONEC_DIST_SUPPS_Y(I,J),2)
            VALUES(J,3)=COORDS_DIST_SUPPS_Y(CONEC_DIST_SUPPS_Y(I,J),3)
        ENDDO
        CALL NORMAL_OUTWARD(4,VALUES,0.D0,0.D0,ETA_DIST_SUPPS_Y(I,:),VJC,DVJC,JC)
    ENDDO
    DO I=1,NUM_DIST_SUPPS_Z
        DO J=1,4
            VALUES(J,1)=COORDS_DIST_SUPPS_Z(CONEC_DIST_SUPPS_Z(I,J),1)
            VALUES(J,2)=COORDS_DIST_SUPPS_Z(CONEC_DIST_SUPPS_Z(I,J),2)
            VALUES(J,3)=COORDS_DIST_SUPPS_Z(CONEC_DIST_SUPPS_Z(I,J),3)
        ENDDO
        CALL NORMAL_OUTWARD(4,VALUES,0.D0,0.D0,ETA_DIST_SUPPS_Z(I,:),VJC,DVJC,JC)
    ENDDO
!    
    OPEN(90,file='LIXO',status='unknown')
    WRITE(90,500)'XBEM_DIST_OUTPUT.ogl'
	CLOSE(90)
!	
	OPEN(90,file='LIXO',status='OLD')
		READ(90,500)NOME
	CLOSE(90, STATUS = 'DELETE')
	
	NOME=ADJUSTL(NOME)	
		
	OPEN(3,file='Output_data\'//NOME,status='unknown')
!
	WRITE(3,*)'Arquivo de p�s-processamento'
	WRITE(3,*)
	WRITE(3,*)'N� de nos     N� de elementos     N� de listas'
	WRITE(3,10)
	WRITE(3,100)N_COLLOCPOINTS+4*NUM_DIST_LOADS+NODES_DIST_SUPPS,N_ELEM-N_CRACKED_ELEM+NUM_DIST_LOADS+NUM_DIST_SUPPS,3
    WRITE(3,*)
	WRITE(3,*)'coordx coordy coordz deslx delsy deslz'
	WRITE(3,10)
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,200)COORD_COLLOCPOINTS(I,1),COORD_COLLOCPOINTS(I,2),COORD_COLLOCPOINTS(I,3),0.D0,0.D0,0.D0
    ENDDO
    DO I=1,NUM_DIST_LOADS
        DO J=1,4
            WRITE(3,200)COORDS_DIST_LOADS(I,J,1),COORDS_DIST_LOADS(I,J,2),COORDS_DIST_LOADS(I,J,3),0.D0,0.D0,0.D0
        ENDDO
    ENDDO
    DO I=1,NODES_DIST_SUPPS_X
        WRITE(3,200)COORDS_DIST_SUPPS_X(I,1),COORDS_DIST_SUPPS_X(I,2),COORDS_DIST_SUPPS_X(I,3),0.D0,0.D0,0.D0
    ENDDO
    DO I=1,NODES_DIST_SUPPS_Y
        WRITE(3,200)COORDS_DIST_SUPPS_Y(I,1),COORDS_DIST_SUPPS_Y(I,2),COORDS_DIST_SUPPS_Y(I,3),0.D0,0.D0,0.D0
    ENDDO
    DO I=1,NODES_DIST_SUPPS_Z
        WRITE(3,200)COORDS_DIST_SUPPS_Z(I,1),COORDS_DIST_SUPPS_Z(I,2),COORDS_DIST_SUPPS_Z(I,3),0.D0,0.D0,0.D0
    ENDDO
	WRITE(3,*)
	WRITE(3,*)'tpelem (1 - barra / 2 - triang / 3 - quad) grauaprox n�1 n�2...non'
	WRITE(3,10)
	DO I=1,N_ELEM
	    IF(CRACKED_ELEM(I).EQ."UNCRACKED")THEN	
	        NEN=ELEM_TYPE(I)*ORDER_ELEM(I)+(ELEM_TYPE(I)-3)*(ORDER_ELEM(I)-1)*POL_FAMILY(I)
            SELECT CASE(NEN)
            CASE(3)
                WRITE(3,300)2,1,COLLOCPOINTS_CONNECTIVITY(I,3),COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,2),ND(I)
            CASE(6)
                WRITE(3,301)2,2,COLLOCPOINTS_CONNECTIVITY(I,3),COLLOCPOINTS_CONNECTIVITY(I,6),COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,5),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,2),ND(I)
            CASE(4)
                WRITE(3,302)3,1,COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,2),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,3),ND(I)
            CASE(8)
                WRITE(3,303)3,2,COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,5),COLLOCPOINTS_CONNECTIVITY(I,2),COLLOCPOINTS_CONNECTIVITY(I,8),COLLOCPOINTS_CONNECTIVITY(I,6),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,7),COLLOCPOINTS_CONNECTIVITY(I,3),ND(I)
            CASE(9)
                WRITE(3,304)3,2,COLLOCPOINTS_CONNECTIVITY(I,1),COLLOCPOINTS_CONNECTIVITY(I,5),COLLOCPOINTS_CONNECTIVITY(I,2),COLLOCPOINTS_CONNECTIVITY(I,8),COLLOCPOINTS_CONNECTIVITY(I,9),COLLOCPOINTS_CONNECTIVITY(I,6),COLLOCPOINTS_CONNECTIVITY(I,4),COLLOCPOINTS_CONNECTIVITY(I,7),COLLOCPOINTS_CONNECTIVITY(I,3),ND(I)            
            ENDSELECT
        ENDIF
    ENDDO
    DO I=1,NUM_DIST_LOADS
        CONT=N_COLLOCPOINTS+4*(I-1)
        WRITE(3,302)3,1,CONT+1,CONT+2,CONT+4,CONT+3,0
    ENDDO
    CONT=N_COLLOCPOINTS+4*NUM_DIST_LOADS
    DO I=1,NUM_DIST_SUPPS_X
        WRITE(3,302)3,1,CONT+CONEC_DIST_SUPPS_X(I,1),CONT+CONEC_DIST_SUPPS_X(I,2),CONT+CONEC_DIST_SUPPS_X(I,4),CONT+CONEC_DIST_SUPPS_X(I,3),0
    ENDDO
    CONT=CONT+NODES_DIST_SUPPS_X
    DO I=1,NUM_DIST_SUPPS_Y
        WRITE(3,302)3,1,CONT+CONEC_DIST_SUPPS_Y(I,1),CONT+CONEC_DIST_SUPPS_Y(I,2),CONT+CONEC_DIST_SUPPS_Y(I,4),CONT+CONEC_DIST_SUPPS_Y(I,3),0
    ENDDO
    CONT=CONT+NODES_DIST_SUPPS_Y
    DO I=1,NUM_DIST_SUPPS_Z
        WRITE(3,302)3,1,CONT+CONEC_DIST_SUPPS_Z(I,1),CONT+CONEC_DIST_SUPPS_Z(I,2),CONT+CONEC_DIST_SUPPS_Z(I,4),CONT+CONEC_DIST_SUPPS_Z(I,3),0
    ENDDO
	WRITE(3,*)
	WRITE(3,*)'ETAX_dist_supps'
	WRITE(3,10)
	WRITE(3,*)'ETAX_dist_supps'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)0.D0,0.D0,0.D0,0.D0
    ENDDO
    DO I=1,NUM_DIST_LOADS
        DO J=1,4
	        WRITE(3,400)0.D0,0.D0,0.D0,ETA_DIST_LOADS(I,1)
        ENDDO
    ENDDO
    DO I=1,NUM_DIST_SUPPS_X
        DO J=1,4
	        WRITE(3,400)0.D0,0.D0,0.D0,ETA_DIST_SUPPS_X(I,1)
        ENDDO
    ENDDO
    DO I=1,NUM_DIST_SUPPS_Y
        DO J=1,4
	        WRITE(3,400)0.D0,0.D0,0.D0,ETA_DIST_SUPPS_Y(I,1)
        ENDDO
    ENDDO
    DO I=1,NUM_DIST_SUPPS_Z
        DO J=1,4
	        WRITE(3,400)0.D0,0.D0,0.D0,ETA_DIST_SUPPS_Z(I,1)
        ENDDO
    ENDDO
	WRITE(3,*)
    	WRITE(3,*)'ETAY_dist_supps'
	WRITE(3,10)
	WRITE(3,*)'ETAY_dist_supps'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)0.D0,0.D0,0.D0,0.D0
    ENDDO
    DO I=1,NUM_DIST_LOADS
        DO J=1,4
	        WRITE(3,400)0.D0,0.D0,0.D0,ETA_DIST_LOADS(I,2)
        ENDDO
    ENDDO
    DO I=1,NUM_DIST_SUPPS_X
        DO J=1,4
	        WRITE(3,400)0.D0,0.D0,0.D0,ETA_DIST_SUPPS_X(I,2)
        ENDDO
    ENDDO
    DO I=1,NUM_DIST_SUPPS_Y
        DO J=1,4
	        WRITE(3,400)0.D0,0.D0,0.D0,ETA_DIST_SUPPS_Y(I,2)
        ENDDO
    ENDDO
    DO I=1,NUM_DIST_SUPPS_Z
        DO J=1,4
	        WRITE(3,400)0.D0,0.D0,0.D0,ETA_DIST_SUPPS_Z(I,2)
        ENDDO
    ENDDO
	WRITE(3,*)
    	WRITE(3,*)'ETAZ_dist_supps'
	WRITE(3,10)
	WRITE(3,*)'ETAZ_dist_supps'
	DO I=1,N_COLLOCPOINTS
	    WRITE(3,400)0.D0,0.D0,0.D0,0.D0
    ENDDO
    DO I=1,NUM_DIST_LOADS
        DO J=1,4
	        WRITE(3,400)0.D0,0.D0,0.D0,ETA_DIST_LOADS(I,3)
        ENDDO
    ENDDO
    DO I=1,NUM_DIST_SUPPS_X
        DO J=1,4
	        WRITE(3,400)0.D0,0.D0,0.D0,ETA_DIST_SUPPS_X(I,3)
        ENDDO
    ENDDO
    DO I=1,NUM_DIST_SUPPS_Y
        DO J=1,4
	        WRITE(3,400)0.D0,0.D0,0.D0,ETA_DIST_SUPPS_Y(I,3)
        ENDDO
    ENDDO
    DO I=1,NUM_DIST_SUPPS_Z
        DO J=1,4
	        WRITE(3,400)0.D0,0.D0,0.D0,ETA_DIST_SUPPS_Z(I,3)
        ENDDO
    ENDDO
	WRITE(3,*)
	CLOSE(3)
!
10  FORMAT(1h#)
100	FORMAT(2x,i5,12X,i5,12X,i5)
200	FORMAT(8x,F14.9,16x,F14.9,16x,F14.9,16x,F14.9,16X,F14.9,16x,F14.9)
300	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7X,i5,7x,i2)
301	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i2)
302	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i2)
303	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i2)
304	FORMAT(8x,i1,7x,i1,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i5,7x,i2)
400	FORMAT(8x,E13.6,4x,E13.6,4x,E13.6,15x,E13.6)
500 FORMAT(a,i4,a)
!
    END SUBROUTINE