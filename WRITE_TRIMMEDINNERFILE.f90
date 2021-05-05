	SUBROUTINE WRITE_TRIMMEDINNERFILE(PATCH,PARAMETRIC_VERTICES,N_CURVEPOINTS,TRIMMING_CURVEPOINTS,SMART_SIZE)  
!
	USE NURBS_SURFACES
	USE TRIMMING_CURVES
!
	IMPLICIT NONE 

    INTEGER::I,J,CONT1,CONT2,PATCH,N_CURVEPOINTS,SMART_SIZE
    REAL*8::PARAMETRIC_VERTICES(4,2),TRIMMING_CURVEPOINTS(N_CURVEPOINTS,2)   
!   ALTERA DE ACORDO COM O LOCAL QUE O ANSYS SALVA
!    OPEN(13,FILE='C:\Users\sergiogustavo\filename.log',status='unknown')		
    OPEN(13,FILE='P:\filename.log',status='unknown')	
!	
	WRITE(13,10)'/BATCH' 
	WRITE(13,10)'/COM,ANSYS RELEASE 12.0.1  UP20090415       08:56:18    04/10/2017'	   
	WRITE(13,10)'/input,menust,tmp,'',,,,,,,,,,,,,,,,1'
	WRITE(13,10)'/GRA,POWER'
	WRITE(13,10)'/GST,ON'
	WRITE(13,10)'/PLO,INFO,3'	
	WRITE(13,10)'/GRO,CURL,ON'				    
	WRITE(13,10)'/CPLANE,1'
	WRITE(13,10)'/REPLOT,RESIZE'
	WRITE(13,10)'WPSTYLE,,,,,,,,0'
	WRITE(13,10)'/PREP7'		
    WRITE(13,10)'!*'
    WRITE(13,10)'ET,1,PLANE42'
    WRITE(13,10)'!*'
    WRITE(13,10)'!*'
    WRITE(13,10)'MPTEMP,,,,,,,,'                		
    WRITE(13,10)'MPTEMP,1,0'  
    WRITE(13,10)'MPDATA,EX,1,,10' 
    WRITE(13,10)'MPDATA,PRXY,1,,0.2'  
    WRITE(13,20)'K,1,',PARAMETRIC_VERTICES(1,1),',',PARAMETRIC_VERTICES(1,2),',,' 
    WRITE(13,20)'K,2,',PARAMETRIC_VERTICES(2,1),',',PARAMETRIC_VERTICES(2,2),',,'             
    WRITE(13,20)'K,3,',PARAMETRIC_VERTICES(3,1),',',PARAMETRIC_VERTICES(3,2),',,'
    WRITE(13,20)'K,4,',PARAMETRIC_VERTICES(4,1),',',PARAMETRIC_VERTICES(4,2),',,'
    WRITE(13,10)'LSTR,       1,       2'   
    WRITE(13,10)'LSTR,       2,       3'
    WRITE(13,10)'LSTR,       3,       4'
    WRITE(13,10)'LSTR,       4,       1'
!    
    CONT1=0
    CONT2=0
    DO I=1,N_TRIMMINGCURVES
        IF(SURFACES(I).EQ.PATCH)THEN
            DO J=1,SUBDIVISION(I)
                CONT1=CONT1+1
                CONT2=CONT2+1
                WRITE(13,30)'K,',4+CONT2,',',TRIMMING_CURVEPOINTS(CONT1,1),',',TRIMMING_CURVEPOINTS(CONT1,2),',,'
            ENDDO
            CONT1=CONT1+1        
        ENDIF
    ENDDO
!                          
    CONT1=0
    CONT2=0 
    DO I=1,N_TRIMMINGCURVES
        IF(SURFACES(I).EQ.PATCH)THEN
            IF(OPENED_CLOSED_CURVE(I).EQ."O".AND.LAST_OPENED_CURVE(I).EQ."N")THEN           
                DO J=1,SUBDIVISION(I)
                    CONT1=CONT1+1
                    CONT2=CONT2+1
                    WRITE(13,40)'LSTR,     ',4+CONT1,',     ',4+CONT1+1       
                ENDDO                   
            ELSE IF(OPENED_CLOSED_CURVE(I).EQ."O".AND.LAST_OPENED_CURVE(I).EQ."S")THEN                 
                DO J=1,SUBDIVISION(I)
                    CONT1=CONT1+1
                    CONT2=CONT2+1
                    IF(J.LT.SUBDIVISION(I))THEN
                        WRITE(13,40)'LSTR,     ',4+CONT1,',     ',4+CONT1+1
                    ELSE
                        WRITE(13,40)'LSTR,     ',4+CONT1,',     ',4+CONT1+1-CONT2
                    ENDIF
                ENDDO    
                CONT2=0                                         
            ELSE IF(OPENED_CLOSED_CURVE(I).EQ."C")THEN
                CONT2=0
                DO J=1,SUBDIVISION(I)
                    CONT1=CONT1+1
                    CONT2=CONT2+1
                    IF(J.LT.SUBDIVISION(I))THEN
                        WRITE(13,40)'LSTR,     ',4+CONT1,',     ',4+CONT1+1
                    ELSE
                        WRITE(13,40)'LSTR,     ',4+CONT1,',     ',4+CONT1+1-CONT2
                    ENDIF    
                ENDDO      
            ENDIF
        ENDIF                 
    ENDDO          
!
	WRITE(13,50)'FLST,2,',4+CONT1,',4'
	DO I=1,4+CONT1
	    WRITE(13,60)'FITEM,2,',I
	ENDDO
    WRITE(13,10)'AL,P51X'	    			
	WRITE(13,60)'SMRT,',SUBDIVISION_PD1(PATCH)
	WRITE(13,50)'MSHAPE,',SUBDIVISION_PD2(PATCH),',2D'
	WRITE(13,10)'MSHKEY,0'
	WRITE(13,10)'!*' 
	WRITE(13,10)'CM,_Y,AREA'
	WRITE(13,10)'ASEL, , , ,       1'
	WRITE(13,10)'CM,_Y1,AREA'
	WRITE(13,10)'CMSEL,S,_Y'  
	WRITE(13,10)'!*' 
	WRITE(13,10)'AMESH,_Y1' 	    
    WRITE(13,10)'!*'
    WRITE(13,10)'CMDELE,_Y' 
    WRITE(13,10)'CMDELE,_Y1'  
    WRITE(13,10)'CMDELE,_Y2'        
	WRITE(13,10)'!*'     
    WRITE(13,10)'FINISH' 
    WRITE(13,10)'/POST1'         
    WRITE(13,10)'FINISH'  
    WRITE(13,10)'CDWRITE,DB,parametric_mesh,txt' 
    WRITE(13,10)'! /EXIT,MODEL' 
! 
	CLOSE(13)
!
10  FORMAT(a)    
20  FORMAT(a,F11.6,a,F11.6,a) 
30  FORMAT(a,i3,a,F11.6,a,F11.6,a)
40  FORMAT(a,i3,a,i3) 
50  FORMAT(a,i3,a)
60  FORMAT(a,i3)    
!            
	END SUBROUTINE    