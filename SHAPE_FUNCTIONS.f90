	SUBROUTINE SHAPE_FUNCTIONS(FLAG,A1,A2,NEN,VALUES,COEFFICIENTS,X1,X2,X3,PHI,DPHI,D2PHI)
!
!	IN THIS SUBROUTINE, THE USER MUST INFORM THE QSI, THE ORDER OF THE 
!	ELEMENT, THE VALUES TO BE INTERPOLATED AND THE RESULTS ARE STORED IN
!	X1, X2, X3, PHI AND DPHI VARIABLES
!
!	IF FLAG=1, A VALUE IS INTERPOLATED USING THE SHAPE FUNCTIONS
!	THE RESULTS ARE STORED IN X, Y AND Z VARIABLES.
!	IF FLAG=2, THE VALUES OF THE INTERPOLATION FUNCTIONS ARE CALCULATED FOR THE 
!	GIVEN ADIMENSIONAL COORDINATES A1 AND A2. THE RESULTS ARE STORED IN PHI VARIABLE
!	IF FLAG=3, THE FIRST DERIVATIVES OF THE SHAPE FUNCTIONS AT THE GIVEN ADIMENSIONAL
!	COORDINATES A1 AND A2 ARE CALCULATED.THE RESULTS ARE STORED IN DPHI VARIABLE
!	IF FLAG=4 THE FIRST DERIVATIVES OF THE INTERPOLATION FUNCTIONS AT THE GIVEN ADIMENSIONAL
!	COORDINATES A1 AND A2 ARE CALCULATED.THE RESULTS ARE STORED IN DPHI VARIABLE
!	IF FLAG=5, THE SECOND DERIVATIVES OF THE SHAPE FUNCTIONS AT THE GIVEN ADIMENSIONAL
!	COORDINATES A1 AND A2 ARE CALCULATED.THE RESULTS ARE STORED IN D2PHI VARIABLE
!
!   A1 = DIMENSIONLESS COORDENATE QSI1
!   A2 = DIMENSIONLESS COORDENATE QSI2
!   ELEM = ELEMENT NUMBER
!   VALUES = NODAL OR COLLOCATION POINTS VALUES THAT ARE USED TOGETHER WITH THE SHAPE FUNCTIONS TO INTERPOLATE 
!
	USE ISOPARAMETRIC_MESH
!
	IMPLICIT NONE 
!
	INTEGER::I,NEN,FLAG
!    		 
	REAL*8::A1,A2,VALUES(NEN,3),X1,X2,X3,PHI(NEN),DPHI(NEN,2),D2PHI(NEN,2,2),COEFFICIENTS(NEN,NEN)
!	
    SELECT CASE(FLAG)
    CASE(1)
		IF (NEN.EQ.3)THEN
		    X1=(A1)*VALUES(1,1)+(A2)*VALUES(2,1)+(1-A1-A2)*VALUES(3,1)
		    X2=(A1)*VALUES(1,2)+(A2)*VALUES(2,2)+(1-A1-A2)*VALUES(3,2)
		    X3=(A1)*VALUES(1,3)+(A2)*VALUES(2,3)+(1-A1-A2)*VALUES(3,3)
		ELSE IF(NEN.EQ.6)THEN
		    X1=(A1*(2*A1-1))*VALUES(1,1)+(A2*(2*A2-1))*VALUES(2,1)+((1-A1-A2)*(1-2*A1-2*A2))*VALUES(3,1)+&
		      (4*A1*A2)*VALUES(4,1)+(4*A2*(1-A1-A2))*VALUES(5,1)+(4*A1*(1-A1-A2))*VALUES(6,1)
		    X2=(A1*(2*A1-1))*VALUES(1,2)+(A2*(2*A2-1))*VALUES(2,2)+((1-A1-A2)*(1-2*A1-2*A2))*VALUES(3,2)+&
              (4*A1*A2)*VALUES(4,2)+(4*A2*(1-A1-A2))*VALUES(5,2)+(4*A1*(1-A1-A2))*VALUES(6,2)		           
		    X3=(A1*(2*A1-1))*VALUES(1,3)+(A2*(2*A2-1))*VALUES(2,3)+((1-A1-A2)*(1-2*A1-2*A2))*VALUES(3,3)+&
		      (4*A1*A2)*VALUES(4,3)+(4*A2*(1-A1-A2))*VALUES(5,3)+(4*A1*(1-A1-A2))*VALUES(6,3)
		ELSE IF(NEN.EQ.4)THEN
		    X1=((1-A1)*(1-A2)/4)*VALUES(1,1)+((1+A1)*(1-A2)/4)*VALUES(2,1)+((1+A1)*(1+A2)/4)*VALUES(3,1)+((1-A1)*(1+A2)/4)*VALUES(4,1)
		    X2=((1-A1)*(1-A2)/4)*VALUES(1,2)+((1+A1)*(1-A2)/4)*VALUES(2,2)+((1+A1)*(1+A2)/4)*VALUES(3,2)+((1-A1)*(1+A2)/4)*VALUES(4,2)
		    X3=((1-A1)*(1-A2)/4)*VALUES(1,3)+((1+A1)*(1-A2)/4)*VALUES(2,3)+((1+A1)*(1+A2)/4)*VALUES(3,3)+((1-A1)*(1+A2)/4)*VALUES(4,3)
		ELSE IF(NEN.EQ.8)THEN
			X1=((1-A1)*(A2-1)*(A2+A1+1)/4)*VALUES(1,1)+((A1+1)*(A2-1)*(A2-A1+1)/4)*VALUES(2,1)+((A1+1)*(A2+1)*(A2+A1-1)/4)*VALUES(3,1)+((1-A1)*(A2+1)*(A2-A1-1)/4)*VALUES(4,1)+&
		      ((A1-1)*(A1+1)*(A2-1)/2)*VALUES(5,1)+((-A1-1)*(A2-1)*(A2+1)/2)*VALUES(6,1)+((1-A1)*(A1+1)*(A2+1)/2)*VALUES(7,1)+((A1-1)*(A2-1)*(A2+1)/2)*VALUES(8,1)
		    X2=((1-A1)*(A2-1)*(A2+A1+1)/4)*VALUES(1,2)+((A1+1)*(A2-1)*(A2-A1+1)/4)*VALUES(2,2)+((A1+1)*(A2+1)*(A2+A1-1)/4)*VALUES(3,2)+((1-A1)*(A2+1)*(A2-A1-1)/4)*VALUES(4,2)+&
		      ((A1-1)*(A1+1)*(A2-1)/2)*VALUES(5,2)+((-A1-1)*(A2-1)*(A2+1)/2)*VALUES(6,2)+((1-A1)*(A1+1)*(A2+1)/2)*VALUES(7,2)+((A1-1)*(A2-1)*(A2+1)/2)*VALUES(8,2)
		    X3=((1-A1)*(A2-1)*(A2+A1+1)/4)*VALUES(1,3)+((A1+1)*(A2-1)*(A2-A1+1)/4)*VALUES(2,3)+((A1+1)*(A2+1)*(A2+A1-1)/4)*VALUES(3,3)+((1-A1)*(A2+1)*(A2-A1-1)/4)*VALUES(4,3)+&
		      ((A1-1)*(A1+1)*(A2-1)/2)*VALUES(5,3)+((-A1-1)*(A2-1)*(A2+1)/2)*VALUES(6,3)+((1-A1)*(A1+1)*(A2+1)/2)*VALUES(7,3)+((A1-1)*(A2-1)*(A2+1)/2)*VALUES(8,3)		
		ELSE IF(NEN.EQ.9)THEN
			X1=((A1-1)*(A2-1)*A1*A2/4)*VALUES(1,1)+((A1+1)*(A2-1)*A1*A2/4)*VALUES(2,1)+((A1+1)*(A2+1)*A1*A2/4)*VALUES(3,1)+((A1-1)*(A2+1)*A1*A2/4)*VALUES(4,1)+&
		      ((1-A1**2)*A2*(A2-1)/2)*VALUES(5,1)+((1-A2**2)*A1*(A1+1)/2)*VALUES(6,1)+((1-A1**2)*A2*(A2+1)/2)*VALUES(7,1)+((1-A2**2)*A1*(A1-1)/2)*VALUES(8,1)+((1-A1**2)*(1-A2**2))*VALUES(9,1)
            X2=((A1-1)*(A2-1)*A1*A2/4)*VALUES(1,2)+((A1+1)*(A2-1)*A1*A2/4)*VALUES(2,2)+((A1+1)*(A2+1)*A1*A2/4)*VALUES(3,2)+((A1-1)*(A2+1)*A1*A2/4)*VALUES(4,2)+&
		      ((1-A1**2)*A2*(A2-1)/2)*VALUES(5,2)+((1-A2**2)*A1*(A1+1)/2)*VALUES(6,2)+((1-A1**2)*A2*(A2+1)/2)*VALUES(7,2)+((1-A2**2)*A1*(A1-1)/2)*VALUES(8,2)+((1-A1**2)*(1-A2**2))*VALUES(9,2)
		    X3=((A1-1)*(A2-1)*A1*A2/4)*VALUES(1,3)+((A1+1)*(A2-1)*A1*A2/4)*VALUES(2,3)+((A1+1)*(A2+1)*A1*A2/4)*VALUES(3,3)+((A1-1)*(A2+1)*A1*A2/4)*VALUES(4,3)+&
		      ((1-A1**2)*A2*(A2-1)/2)*VALUES(5,3)+((1-A2**2)*A1*(A1+1)/2)*VALUES(6,3)+((1-A1**2)*A2*(A2+1)/2)*VALUES(7,3)+((1-A2**2)*A1*(A1-1)/2)*VALUES(8,3)+((1-A1**2)*(1-A2**2))*VALUES(9,3)
		ENDIF
    CASE(2)
		IF (NEN.EQ.3)THEN
		    DO I=1,NEN			    	
		        PHI(I)=COEFFICIENTS(I,1)*1.D0+COEFFICIENTS(I,2)*A1+COEFFICIENTS(I,3)*A2
		    ENDDO
		ELSE IF(NEN.EQ.6)THEN
		    DO I=1,NEN			    	
		        PHI(I)=COEFFICIENTS(I,1)*1.D0+COEFFICIENTS(I,2)*A1+COEFFICIENTS(I,3)*A2+COEFFICIENTS(I,4)*A1*A2+COEFFICIENTS(I,5)*A1*A1+COEFFICIENTS(I,6)*A2*A2
		    ENDDO
		ELSE IF(NEN.EQ.4)THEN
			DO I=1,NEN			    	
		        PHI(I)=COEFFICIENTS(I,1)*1.D0+COEFFICIENTS(I,2)*A1+COEFFICIENTS(I,3)*A2+COEFFICIENTS(I,4)*A1*A2
		    ENDDO  	    		    
		ELSE IF(NEN.EQ.8)THEN
			DO I=1,NEN			    	
		        PHI(I)=COEFFICIENTS(I,1)*1.D0+COEFFICIENTS(I,2)*A1+COEFFICIENTS(I,3)*A2+COEFFICIENTS(I,4)*A1*A2+COEFFICIENTS(I,5)*A1*A1+COEFFICIENTS(I,6)*A2*A2+COEFFICIENTS(I,7)*A1*A2*A2+COEFFICIENTS(I,8)*A2*A1*A1
		    ENDDO  
		ELSE IF(NEN.EQ.9)THEN
			DO I=1,NEN			    	
		        PHI(I)=COEFFICIENTS(I,1)*1.D0+COEFFICIENTS(I,2)*A1+COEFFICIENTS(I,3)*A2+COEFFICIENTS(I,4)*A1*A2+COEFFICIENTS(I,5)*A1*A1+COEFFICIENTS(I,6)*A2*A2+COEFFICIENTS(I,7)*A1*A2*A2+COEFFICIENTS(I,8)*A2*A1*A1+COEFFICIENTS(I,9)*A1*A1*A2*A2
		    ENDDO  		        		    
		ENDIF
	CASE(3)
		IF (NEN.EQ.3)THEN
		    DPHI(1,1)=1.D0
		    DPHI(1,2)=0.D0
		    DPHI(2,1)=0.D0
		    DPHI(2,2)=1.D0
		    DPHI(3,1)=-1.D0
		    DPHI(3,2)=-1.D0 	    
		ELSE IF(NEN.EQ.6)THEN
		    DPHI(1,1)=4.D0*A1-1.D0
		    DPHI(1,2)=0
		    DPHI(2,1)=0
		    DPHI(2,2)=4.D0*A2-1.D0
		    DPHI(3,1)=4.D0*A1+4.D0*A2-3.D0
		    DPHI(3,2)=4.D0*A1+4.D0*A2-3.D0
		    DPHI(4,1)=4.D0*A2
		    DPHI(4,2)=4.D0*A1
		    DPHI(5,1)=-4.D0*A2
		    DPHI(5,2)=-8.D0*A2-4.D0*A1+4.D0
		    DPHI(6,1)=-8.D0*A1-4.D0*A2+4.D0
		    DPHI(6,2)=-4.D0*A1	 		    			    				
		ELSE IF(NEN.EQ.4)THEN
		    DPHI(1,1)=(A2-1)/4.D0
		    DPHI(1,2)=(A1-1)/4.D0
		    DPHI(2,1)=(1-A2)/4.D0
		    DPHI(2,2)=(-1-A1)/4.D0
		    DPHI(3,1)=(A2+1)/4.D0
		    DPHI(3,2)=(A1+1)/4.D0
		    DPHI(4,1)=(-1-A2)/4.D0
		    DPHI(4,2)=(1-A1)/4.D0		
		ELSE IF(NEN.EQ.8)THEN
		    DPHI(1,1)=(-2.D0*A2*A1-A2**2.D0+2.D0*A1+A2)/4.D0
		    DPHI(1,2)=(-A1**2.D0-2.D0*A1*A2+2.D0*A2+A1)/4.D0
		    DPHI(2,1)=(-2.D0*A2*A1+A2**2.D0+2.D0*A1-A2)/4.D0
		    DPHI(2,2)=(-A1**2.D0+2.D0*A1*A2+2.D0*A2-A1)/4.D0
		    DPHI(3,1)=(2.D0*A2*A1+A2**2.D0+2.D0*A1+A2)/4.D0
		    DPHI(3,2)=(A1**2.D0+2.D0*A1*A2+2.D0*A2+A1)/4.D0
		    DPHI(4,1)=(2.D0*A2*A1-A2**2.D0+2.D0*A1-A2)/4.D0
		    DPHI(4,2)=(A1**2.D0-2.D0*A1*A2+2.D0*A2-A1)/4.D0
		    DPHI(5,1)=A1*A2-A1
		    DPHI(5,2)=(A1**2.D0-1.D0)/2.D0
		    DPHI(6,1)=(-A2**2+1.D0)/2.D0
		    DPHI(6,2)=-A1*A2-A2		
		    DPHI(7,1)=-A1*A2-A1	
		    DPHI(7,2)=(-A1**2.D0+1.D0)/2.D0	
		    DPHI(8,1)=(A2**2.D0-1.D0)/2.D0
		    DPHI(8,2)=A1*A2-A2
		ELSE IF(NEN.EQ.9)THEN
		    DPHI(1,1)=(A2-1.D0)*A2*(2.D0*A1-1.D0)/4.D0
		    DPHI(1,2)=(A1-1.D0)*A1*(2.D0*A2-1.D0)/4.D0
		    DPHI(2,1)=(A2-1.D0)*A2*(2.D0*A1+1.D0)/4.D0
		    DPHI(2,2)=(A1+1.D0)*A1*(2.D0*A2-1.D0)/4.D0
		    DPHI(3,1)=(A2+1.D0)*A2*(2.D0*A1+1.D0)/4.D0
		    DPHI(3,2)=(A1+1.D0)*A1*(2.D0*A2+1.D0)/4.D0
		    DPHI(4,1)=(A2+1.D0)*A2*(2.D0*A1-1.D0)/4.D0
		    DPHI(4,2)=(A1-1.D0)*A1*(2.D0*A2+1.D0)/4.D0	    
		    DPHI(5,1)=(1.D0-A2)*A2*A1
		    DPHI(5,2)=(1.D0-A1**2.D0)*(2.D0*A2-1.D0)/2.D0    
		    DPHI(6,1)=(1.D0-A2**2.D0)*(2.D0*A1+1.D0)/2.D0
		    DPHI(6,2)=-(1.D0+A1)*A1*A2   
		    DPHI(7,1)=-(1.D0+A2)*A2*A1
		    DPHI(7,2)=(1.D0-A1**2.D0)*(2.D0*A2+1.D0)/2.D0  
		    DPHI(8,1)=(1.D0-A2**2.D0)*(2.D0*A1-1.D0)/2.D0
		    DPHI(8,2)=(1.D0-A1)*A1*A2    		    
		    DPHI(9,1)=2.D0*A1*(A2**2.D0-1.D0)
		    DPHI(9,2)=2.D0*A2*(A1**2.D0-1.D0)   		    
		ENDIF
	CASE(4)
		IF (NEN.EQ.3)THEN
		    DO I=1,NEN
		        DPHI(I,1)=COEFFICIENTS(I,2)
		        DPHI(I,2)=COEFFICIENTS(I,3)
		    ENDDO		    
		ELSE IF(NEN.EQ.6)THEN
			DO I=1,NEN
		        DPHI(I,1)=COEFFICIENTS(I,2)+COEFFICIENTS(I,4)*A2+2.D0*COEFFICIENTS(I,5)*A1
		        DPHI(I,2)=COEFFICIENTS(I,3)+COEFFICIENTS(I,4)*A1+2.D0*COEFFICIENTS(I,6)*A2
		    ENDDO		    				
		ELSE IF(NEN.EQ.4)THEN
			DO I=1,NEN
		        DPHI(I,1)=COEFFICIENTS(I,2)+COEFFICIENTS(I,4)*A2
		        DPHI(I,2)=COEFFICIENTS(I,3)+COEFFICIENTS(I,4)*A1
		    ENDDO		
		ELSE IF(NEN.EQ.8)THEN
			DO I=1,NEN
		        DPHI(I,1)=COEFFICIENTS(I,2)+COEFFICIENTS(I,4)*A2+2.D0*COEFFICIENTS(I,5)*A1+COEFFICIENTS(I,7)*A2*A2+2.D0*COEFFICIENTS(I,8)*A2*A1
		        DPHI(I,2)=COEFFICIENTS(I,3)+COEFFICIENTS(I,4)*A1+2.D0*COEFFICIENTS(I,6)*A2+2.D0*COEFFICIENTS(I,7)*A1*A2+COEFFICIENTS(I,8)*A1*A1
		    ENDDO
	    ELSE IF(NEN.EQ.9)THEN
	    	DO I=1,NEN
		        DPHI(I,1)=COEFFICIENTS(I,2)+COEFFICIENTS(I,4)*A2+2.D0*COEFFICIENTS(I,5)*A1+COEFFICIENTS(I,7)*A2*A2+2.D0*COEFFICIENTS(I,8)*A2*A1+2.D0*COEFFICIENTS(I,9)*A1*A2*A2
		        DPHI(I,2)=COEFFICIENTS(I,3)+COEFFICIENTS(I,4)*A1+2.D0*COEFFICIENTS(I,6)*A2+2.D0*COEFFICIENTS(I,7)*A1*A2+COEFFICIENTS(I,8)*A1*A1+2.D0*COEFFICIENTS(I,9)*A1*A1*A2
		    ENDDO					  		        		    
		ENDIF
	CASE(5)
		IF (NEN.EQ.3)THEN
            D2PHI(1,1,1)=0.D0	
            D2PHI(1,1,2)=0.D0	
            D2PHI(1,2,1)=0.D0		 
            D2PHI(1,2,2)=0.D0	   	
            D2PHI(2,1,1)=0.D0		
            D2PHI(2,1,2)=0.D0		
            D2PHI(2,2,1)=0.D0		 
            D2PHI(2,2,2)=0.D0	  
            D2PHI(3,1,1)=0.D0		
            D2PHI(3,1,2)=0.D0		
            D2PHI(3,2,1)=0.D0		 
            D2PHI(3,2,2)=0.D0	  	    
		ELSE IF(NEN.EQ.6)THEN
            D2PHI(1,1,1)=4.D0	
            D2PHI(1,1,2)=0.D0	
            D2PHI(1,2,1)=0.D0	 
            D2PHI(1,2,2)=0.D0  
            D2PHI(2,1,1)=0.D0	
            D2PHI(2,1,2)=0.D0	
            D2PHI(2,2,1)=0.D0	 
            D2PHI(2,2,2)=4.D0 
            D2PHI(3,1,1)=4.D0	
            D2PHI(3,1,2)=4.D0	
            D2PHI(3,2,1)=4.D0	 
            D2PHI(3,2,2)=4.D0
            D2PHI(4,1,1)=0.D0	
            D2PHI(4,1,2)=4.D0	
            D2PHI(4,2,1)=4.D0	 
            D2PHI(4,2,2)=0.D0
            D2PHI(5,1,1)=0.D0	
            D2PHI(5,1,2)=-4.D0	
            D2PHI(5,2,1)=-4.D0	 
            D2PHI(5,2,2)=-8.D0
            D2PHI(6,1,1)=-8.D0	
            D2PHI(6,1,2)=-4.D0	
            D2PHI(6,2,1)=-4.D0	 
            D2PHI(6,2,2)=0.D0   		    			    				
		ELSE IF(NEN.EQ.4)THEN
		    D2PHI(1,1,1)=0.D0	
            D2PHI(1,1,2)=1.D0/4.D0	
            D2PHI(1,2,1)=1.D0/4.D0 
            D2PHI(1,2,2)=0.D0 
            D2PHI(2,1,1)=0.D0	
            D2PHI(2,1,2)=-1.D0/4.D0	
            D2PHI(2,2,1)=-1.D0/4.D0 
            D2PHI(2,2,2)=0.D0 
            D2PHI(3,1,1)=0.D0	
            D2PHI(3,1,2)=1.D0/4.D0	
            D2PHI(3,2,1)=1.D0/4.D0 
            D2PHI(3,2,2)=0.D0 
            D2PHI(4,1,1)=0.D0	
            D2PHI(4,1,2)=-1.D0/4.D0	
            D2PHI(4,2,1)=-1.D0/4.D0 
            D2PHI(4,2,2)=0.D0 
		ELSE IF(NEN.EQ.8)THEN
		    D2PHI(1,1,1)=(-2.D0*A2+2.D0)/4.D0	
            D2PHI(1,1,2)=(-2.D0*A1-2.D0*A2+1.D0)/4.D0	
            D2PHI(1,2,1)=(-2.D0*A1-2.D0*A2+1.D0)/4.D0	
            D2PHI(1,2,2)=(-2.D0*A1+2.D0)/4.D0
		    D2PHI(2,1,1)=(-2.D0*A2+2.D0)/4.D0	
            D2PHI(2,1,2)=(-2.D0*A1+2.D0*A2-1.D0)/4.D0	
            D2PHI(2,2,1)=(-2.D0*A1+2.D0*A2-1.D0)/4.D0	
            D2PHI(2,2,2)=(2.D0*A1+2.D0)/4.D0     
		    D2PHI(3,1,1)=(2.D0*A2+2.D0)/4.D0	
            D2PHI(3,1,2)=(2.D0*A1+2.D0*A2+1.D0)/4.D0	
            D2PHI(3,2,1)=(2.D0*A1+2.D0*A2+1.D0)/4.D0	
            D2PHI(3,2,2)=(2.D0*A1+2.D0)/4.D0  
 		    D2PHI(4,1,1)=(2.D0*A2+2.D0)/4.D0	
            D2PHI(4,1,2)=(2.D0*A1-2.D0*A2-1.D0)/4.D0	
            D2PHI(4,2,1)=(2.D0*A1-2.D0*A2-1.D0)/4.D0		
            D2PHI(4,2,2)=(-2.D0*A1+2.D0)/4.D0            
 		    D2PHI(5,1,1)=A2-1.D0
            D2PHI(5,1,2)=A1
            D2PHI(5,2,1)=A1		
            D2PHI(5,2,2)=0.D0   
            D2PHI(6,1,1)=0.D0
            D2PHI(6,1,2)=-A2
            D2PHI(6,2,1)=-A2		
            D2PHI(6,2,2)=-A1-1.D0      
 		    D2PHI(7,1,1)=-A2-1.D0
            D2PHI(7,1,2)=-A1
            D2PHI(7,2,1)=-A1		
            D2PHI(7,2,2)=0.D0    
            D2PHI(8,1,1)=0.D0
            D2PHI(8,1,2)=A2
            D2PHI(8,2,1)=A2		
            D2PHI(8,2,2)=A1-1.D0                                                        
		ELSE IF(NEN.EQ.9)THEN
		    D2PHI(1,1,1)=2.D0*(A2-1.D0)*A2/4.D0	
            D2PHI(1,1,2)=(2.D0*A2-1.D0)*(2.D0*A1-1.D0)/4.D0	
            D2PHI(1,2,1)=(2.D0*A2-1.D0)*(2.D0*A1-1.D0)/4.D0	
            D2PHI(1,2,2)=2.D0*(A1-1.D0)*A1/4.D0		
            D2PHI(2,1,1)=2.D0*(A2-1.D0)*A2/4.D0	
            D2PHI(2,1,2)=(2.D0*A2-1.D0)*(2.D0*A1+1.D0)/4.D0	
            D2PHI(2,2,1)=(2.D0*A2-1.D0)*(2.D0*A1+1.D0)/4.D0		
            D2PHI(2,2,2)=2.D0*(A1+1.D0)*A1/4.D0	    
            D2PHI(3,1,1)=2.D0*(A2+1.D0)*A2/4.D0	
            D2PHI(3,1,2)=(2.D0*A2+1.D0)*(2.D0*A1+1.D0)/4.D0	
            D2PHI(3,2,1)=(2.D0*A2+1.D0)*(2.D0*A1+1.D0)/4.D0		
            D2PHI(3,2,2)=2.D0*(A1+1.D0)*A1/4.D0	          
            D2PHI(4,1,1)=2.D0*(A2+1.D0)*A2/4.D0	
            D2PHI(4,1,2)=(2.D0*A2+1.D0)*(2.D0*A1-1.D0)/4.D0	
            D2PHI(4,2,1)=(2.D0*A2+1.D0)*(2.D0*A1-1.D0)/4.D0		
            D2PHI(4,2,2)=2.D0*(A1-1.D0)*A1/4.D0  
            D2PHI(5,1,1)=A2*(1.D0-A2)
            D2PHI(5,1,2)=A1*(1.D0-2.D0*A2)
            D2PHI(5,2,1)=A1*(1.D0-2.D0*A2)		
            D2PHI(5,2,2)=(1.D0-A1**2.D0)     
            D2PHI(6,1,1)=(1.D0-A2**2.D0)
            D2PHI(6,1,2)=-A2*(2.D0*A1+1.D0)
            D2PHI(6,2,1)=-A2*(2.D0*A1+1.D0)		
            D2PHI(6,2,2)=-A1*(1.D0+A1)       
            D2PHI(7,1,1)=-A2*(1.D0+A2)
            D2PHI(7,1,2)=-A1*(1.D0+2.D0*A2)
            D2PHI(7,2,1)=-A1*(1.D0+2.D0*A2)		
            D2PHI(7,2,2)=(1.D0-A1**2.D0)
            D2PHI(8,1,1)=(1.D0-A2**2.D0)
            D2PHI(8,1,2)=-A2*(2.D0*A1-1.D0)
            D2PHI(8,2,1)=-A2*(2.D0*A1-1.D0)		
            D2PHI(8,2,2)=A1*(1.D0-A1)  
            D2PHI(9,1,1)=2.D0*(A2**2.D0-1.D0)
            D2PHI(9,1,2)=4.D0*A1*A2
            D2PHI(9,2,1)=4.D0*A1*A2		
            D2PHI(9,2,2)=2.D0*(A1**2.D0-1.D0)                                                                  
		ENDIF			  
    END SELECT
!    
 	END SUBROUTINE SHAPE_FUNCTIONS
 
 
 
 
 	SUBROUTINE SHAPE_FUNCTIONS_COEFFICIENTS(ELEM,NEN,COEFFICIENTS)
!
!	IN THIS SUBROUTINE CALCULATES THE COEFFICIENTS OF THE INTERPOLATION FUNCTIONS CONSIDERING  
!	THE ADIMENSIONAL COLLOCATION POINTS COORDINATES.
!
	USE ISOPARAMETRIC_MESH
!	
 	IMPLICIT NONE 
!
	INTEGER::I,ELEM,NEN
!
	REAL*8::MAT_C(NEN,NEN),VET_C(NEN),COEFFICIENTS(NEN,NEN),COEFFICIENTS_AUX(NEN)
!	
	COEFFICIENTS=0.D0	
	IF (NEN.EQ.3)THEN			    
	    DO I=1,NEN		
		    MAT_C(I,1)=1.D0
		    MAT_C(I,2)=QSI(ELEM,I,1)
		    MAT_C(I,3)=QSI(ELEM,I,2)
		ENDDO
!		    
		VET_C(1)=1.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(1,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=1.D0
		VET_C(3)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(2,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=1.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(3,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
    ELSE IF(NEN.EQ.6)THEN
		DO I=1,NEN		
		    MAT_C(I,1)=1.D0
		    MAT_C(I,2)=QSI(ELEM,I,1)
		    MAT_C(I,3)=QSI(ELEM,I,2)
		    MAT_C(I,4)=QSI(ELEM,I,1)*QSI(ELEM,I,2)
		    MAT_C(I,5)=QSI(ELEM,I,1)*QSI(ELEM,I,1)
		    MAT_C(I,6)=QSI(ELEM,I,2)*QSI(ELEM,I,2)
		ENDDO
!		    
		VET_C(1)=1.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(1,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=1.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(2,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=1.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(3,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=1.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(4,I)=COEFFICIENTS_AUX(I)
		ENDDO
!	
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=1.D0
		VET_C(6)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(5,I)=COEFFICIENTS_AUX(I)
		ENDDO
!	
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=1.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(6,I)=COEFFICIENTS_AUX(I)
		ENDDO
!	
	ELSE IF(NEN.EQ.4)THEN
		DO I=1,NEN		
		    MAT_C(I,1)=1.D0
		    MAT_C(I,2)=QSI(ELEM,I,1)
		    MAT_C(I,3)=QSI(ELEM,I,2)
		    MAT_C(I,4)=QSI(ELEM,I,1)*QSI(ELEM,I,2)
		ENDDO
!		    
		VET_C(1)=1.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)
		DO I=1,NEN	
		    COEFFICIENTS(1,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=1.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(2,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=1.D0
		VET_C(4)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(3,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=1.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(4,I)=COEFFICIENTS_AUX(I)
		ENDDO
!		
    ELSE IF(NEN.EQ.8)THEN
		DO I=1,NEN		
		    MAT_C(I,1)=1.D0
		    MAT_C(I,2)=QSI(ELEM,I,1)
		    MAT_C(I,3)=QSI(ELEM,I,2)
		    MAT_C(I,4)=QSI(ELEM,I,1)*QSI(ELEM,I,2)
		    MAT_C(I,5)=QSI(ELEM,I,1)*QSI(ELEM,I,1)
		    MAT_C(I,6)=QSI(ELEM,I,2)*QSI(ELEM,I,2)
		    MAT_C(I,7)=QSI(ELEM,I,1)*QSI(ELEM,I,2)*QSI(ELEM,I,2)
		    MAT_C(I,8)=QSI(ELEM,I,2)*QSI(ELEM,I,1)*QSI(ELEM,I,1)
		ENDDO
!		    
		VET_C(1)=1.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(1,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=1.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(2,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=1.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(3,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=1.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(4,I)=COEFFICIENTS_AUX(I)
		ENDDO
!	
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=1.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(5,I)=COEFFICIENTS_AUX(I)
		ENDDO
!	
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=1.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(6,I)=COEFFICIENTS_AUX(I)
		ENDDO
!	
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=1.D0
		VET_C(8)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(7,I)=COEFFICIENTS_AUX(I)
		ENDDO
!	
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=1.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(8,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
    ELSE IF(NEN.EQ.9)THEN
		DO I=1,NEN		
		    MAT_C(I,1)=1.D0
		    MAT_C(I,2)=QSI(ELEM,I,1)
		    MAT_C(I,3)=QSI(ELEM,I,2)
		    MAT_C(I,4)=QSI(ELEM,I,1)*QSI(ELEM,I,2)
		    MAT_C(I,5)=QSI(ELEM,I,1)*QSI(ELEM,I,1)
		    MAT_C(I,6)=QSI(ELEM,I,2)*QSI(ELEM,I,2)
		    MAT_C(I,7)=QSI(ELEM,I,1)*QSI(ELEM,I,2)*QSI(ELEM,I,2)
		    MAT_C(I,8)=QSI(ELEM,I,2)*QSI(ELEM,I,1)*QSI(ELEM,I,1)
		    MAT_C(I,9)=QSI(ELEM,I,1)*QSI(ELEM,I,1)*QSI(ELEM,I,2)*QSI(ELEM,I,2)
		ENDDO
!		    
		VET_C(1)=1.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		VET_C(9)=0.D0
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(1,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=1.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		VET_C(9)=0.D0		
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(2,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=1.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		VET_C(9)=0.D0		
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(3,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=1.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		VET_C(9)=0.D0		
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(4,I)=COEFFICIENTS_AUX(I)
		ENDDO
!	
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=1.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		VET_C(9)=0.D0		
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(5,I)=COEFFICIENTS_AUX(I)
		ENDDO
!	
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=1.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		VET_C(9)=0.D0		
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(6,I)=COEFFICIENTS_AUX(I)
		ENDDO
!	
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=1.D0
		VET_C(8)=0.D0
		VET_C(9)=0.D0		
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(7,I)=COEFFICIENTS_AUX(I)
		ENDDO
!	
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=1.D0
		VET_C(9)=0.D0		
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(8,I)=COEFFICIENTS_AUX(I)
		ENDDO
!	
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		VET_C(5)=0.D0
		VET_C(6)=0.D0
		VET_C(7)=0.D0
		VET_C(8)=0.D0
		VET_C(9)=1.D0		
		CALL DLSARG(NEN,MAT_C,NEN,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,NEN	
		    COEFFICIENTS(9,I)=COEFFICIENTS_AUX(I)
		ENDDO
!									  	
 	ENDIF    	
    END SUBROUTINE SHAPE_FUNCTIONS_COEFFICIENTS	
    
    
SUBROUTINE SHAPE_FUNCTIONS_COEFFICIENTS_XBEM_DIST(ELEM,COEFFICIENTS)
!
!	IN THIS SUBROUTINE CALCULATES THE COEFFICIENTS OF THE INTERPOLATION FUNCTIONS CONSIDERING  
!	THE ADIMENSIONAL COLLOCATION POINTS COORDINATES.
!   ELEM = INTEGRATED AREA
!
	USE ISOPARAMETRIC_MESH
    USE XBEM_FORCE_VARIABLES
!	
 	IMPLICIT NONE 
!
	INTEGER::I,ELEM
!
	REAL*8::MAT_C(4,4),VET_C(4),COEFFICIENTS(4,4),COEFFICIENTS_AUX(4)
!	
	COEFFICIENTS=0.D0	
		DO I=1,4		
		    MAT_C(I,1)=1.D0
		    MAT_C(I,2)=QSI_DIST(ELEM,I,1)
		    MAT_C(I,3)=QSI_DIST(ELEM,I,2)
		    MAT_C(I,4)=QSI_DIST(ELEM,I,1)*QSI_DIST(ELEM,I,2)
		ENDDO
!		    
		VET_C(1)=1.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		CALL DLSARG(4,MAT_C,4,VET_C,1,COEFFICIENTS_AUX)
		DO I=1,4	
		    COEFFICIENTS(1,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=1.D0
		VET_C(3)=0.D0
		VET_C(4)=0.D0
		CALL DLSARG(4,MAT_C,4,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,4	
		    COEFFICIENTS(2,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=1.D0
		VET_C(4)=0.D0
		CALL DLSARG(4,MAT_C,4,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,4	
		    COEFFICIENTS(3,I)=COEFFICIENTS_AUX(I)
		ENDDO
!
		VET_C(1)=0.D0
		VET_C(2)=0.D0
		VET_C(3)=0.D0
		VET_C(4)=1.D0
		CALL DLSARG(4,MAT_C,4,VET_C,1,COEFFICIENTS_AUX)	
		DO I=1,4	
		    COEFFICIENTS(4,I)=COEFFICIENTS_AUX(I)
        ENDDO
        
    END SUBROUTINE