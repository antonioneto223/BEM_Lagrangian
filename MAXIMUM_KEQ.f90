	SUBROUTINE MAXIMUM_KEQ
!	
    USE ISOPARAMETRIC_MESH
    USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
	USE REMESHING
	USE PROPAGATION	 
	USE FATIGUE  
!
    IMPLICIT NONE
    INTEGER::I,J,K,ELEMENT,LOCAL_NUMBER,CONT,NEN
!    
    REAL*8::THETA,E,Nu,F 
!
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   COMPUTING THE MAXIMUM EQUIVALENT STRESS INTENSITY FACTORS 
!------------------------------------------------------------------------------------------------------------------------------------------------------
!  
    CONT=0    
    DO I=1,N_COLLOCPOINTS
        IF(DUAL_BEM(I).EQ."S")THEN
            IF(CRACK_TIP(I).EQ.1.OR.CRACK_TIP(I).EQ.3)THEN
                CONT=CONT+1 
                DO J=1,N_ELEM
                    NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                    DO K=1,NEN
                        IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
                            ELEMENT=J
                            LOCAL_NUMBER=K                    
                        ENDIF
                    ENDDO
                ENDDO  
                E=EMP(ND(ELEMENT),1)
                Nu=EMP(ND(ELEMENT),2)

                SELECT CASE(PROPAGATION_CRITERIA)  
                CASE("MERR") ! MAXIMUM ENERGY RELEASE RATE CRITERIA                        
                    CALL MAXIMUM_ENERGY_RELEASE_RATE(Nu,E,KI(CONT),KII(CONT),KIII(CONT),THETAP(CONT))
                    PHIP(CONT)=0.D0  
                    Gmax(CONT)=(1+Nu)/E*COS(THETAP(CONT))**2*((1-Nu**2)/(2*(1+Nu))*(KI(CONT)**2*(1+COS(THETAP(CONT)))-4*KI(CONT)*KII(CONT)*SIN(THETAP(CONT))+KII(CONT)**2*(5-3*COS(THETAP(CONT))))+KIII(CONT)**2) 
                    Keq(CONT)=DSQRT(E*Gmax(CONT)/(1-Nu**2))
                CASE("MPS") ! MAXIMUM PRINCIPAL STRESS CRITERIA
                    CALL MAXIMUM_PRINCIPAL_STRESS(KI(CONT),KII(CONT),KIII(CONT),THETAP(CONT))
                    F=(1-Nu)*KI(CONT)*(3*COS(THETAP(CONT)/2.D0)+COS(3*THETAP(CONT)/2.D0))/4.D0-(1-Nu)*KII(CONT)*(3*SIN(THETAP(CONT)/2.D0)+3*SIN(3*THETAP(CONT)/2.D0))/4.D0-Nu*KI(CONT)*(5*COS(THETAP(CONT)/2.D0)-COS(3*THETAP(CONT)/2.D0))/4.D0+Nu*KII(CONT)*(5*SIN(THETAP(CONT)/2.D0)-3*SIN(3*THETAP(CONT)/2.D0))/4.D0
                    PHIP(CONT)=0.5D0*ATAN(2*KIII(CONT)*COS(THETAP(CONT)/2.D0)/F)
                    Keq(CONT)=0.5D0*COS(THETAP(CONT)/2.D0)*(KI(CONT)*COS(THETAP(CONT))**2-1.5D0*KII(CONT)*SIN(THETAP(CONT))+DSQRT((KI(CONT)*COS(THETAP(CONT)/2.D0)**2-1.5D0*KII(CONT)*SIN(THETAP(CONT)))**2+4*KIII(CONT)**2)) 
                ENDSELECT                      
            ENDIF
        ENDIF         
    ENDDO   
!
    CONT=0
    Keqmax=0.D0
    DO I=1,N_COLLOCPOINTS
        IF(CRACK_TIP(I).EQ.1.AND.DUAL_BEM(I).EQ."S")THEN
            CONT=CONT+1 
            DO J=1,N_ELEM
                NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                DO K=1,NEN
                    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
                        ELEMENT=J
                        LOCAL_NUMBER=K                    
                    ENDIF
                ENDDO
            ENDDO 
            IF(Keq(CONT).GE.Keqmax)THEN
                Keqmax=Keq(CONT)
            ENDIF
        ELSE IF(CRACK_TIP(I).EQ.3.AND.DUAL_BEM(I).EQ."S")THEN
            CONT=CONT+1 
            DO J=1,N_ELEM
                NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
                DO K=1,NEN
                    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I)THEN
                        ELEMENT=J
                        LOCAL_NUMBER=K                    
                    ENDIF
                ENDDO
            ENDDO 
            IF(Keq(CONT).GE.Keqmax)THEN
                Keqmax=Keq(CONT)
            ENDIF
        ENDIF      
    ENDDO       
!    Keqmax=Keq(1)                                 
!
    IF(MINIMUM_MAXIMUM_LOAD.EQ."MINIMUM")THEN 
        DEALLOCATE(COORD_CRACKTIP,KI,KII,KIII,COD,CSD,CTD,THETAP,PHIP,Keq,Gmax,CRACK_INCREMENT,VERTICAL_INCREMENT,CRACKTIP_COLLOCPOINTS)
    ENDIF    
!               	    
    END SUBROUTINE MAXIMUM_KEQ

