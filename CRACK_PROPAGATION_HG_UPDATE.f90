	SUBROUTINE CRACK_PROPAGATION_HG_UPDATE
!
	USE ISOPARAMETRIC_MESH
	USE MATERIALS
	USE CRACKS_DUAL_BEM
	USE SUB_REGIONS_INTERFACES
	USE ANALYSIS
	USE REMESHING
	USE PROPAGATION
!
	IMPLICIT NONE 
!
	INTEGER::I,J,K,II,JJ,KK,JJJ,NEN,N
!
	REAL*8::DG[ALLOCATABLE](:,:),DH[ALLOCATABLE](:,:),DX,DY,DZ,DR,C(3,3)
	REAL*8::VALUES(9,3),LENGTH(4),AVERAGE_LENGTH,DISTANCE,DAUX,SOMA_H(3*N_COLLOCPOINTS)
!	REAL*8::H_AUX(3),G_AUX(3),R_S(3,3),R_H(3,3)               		 
!
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   ADDING THE NEW CRACK FRONT ELEMENTS CONTRIBUTIONS
!------------------------------------------------------------------------------------------------------------------------------------------------------ 
!
    CALL OMP_SET_NUM_THREADS(24)
!$OMP PARALLEL PRIVATE(I,J,K,II,JJ,KK,JJJ,NEN,N,DG,DH,DX,DY,DZ,DR,VALUES,LENGTH,AVERAGE_LENGTH,DISTANCE,DAUX)
!$OMP DO SCHEDULE(DYNAMIC)   
  	DO I=1,N_COLLOCPOINTS-N_CRACKTIPCOLLOCPOINTS
		JJ=0
		DO K=1,N_ELEM
		    NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			DO KK=1,NEN
				IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.I) THEN
					JJ=K
				ENDIF
			ENDDO
	    ENDDO
	    JJJ=ND(JJ)
!	    
		DO J=N_ELEM-N_CRACKTIPELEMENTS+1,N_ELEM	    
		    IF(ND(JJ).EQ.ND(J)) THEN
			    II=0
			    NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
			    DO K=1,NEN
				    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
					    II=1
				    ENDIF
			    ENDDO
!				
			    IF(II.EQ.0) THEN
	                DO K=1,NEN
		                DX=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),1)-COORD_COLLOCPOINTS(I,1)
		                DY=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),2)-COORD_COLLOCPOINTS(I,2)
                        DZ=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),3)-COORD_COLLOCPOINTS(I,3)
		                DR=DSQRT(DX*DX+DY*DY+DZ*DZ)
		                IF(DR.LE.1.E-12) THEN
			                II=1
			            ENDIF
	                ENDDO
                ENDIF
!
		        VALUES=0.D0
!
                DO K=1,NEN
	                VALUES(K,1)=COORD_NODES(NODES_CONNECTIVITY(J,K),1)
	                VALUES(K,2)=COORD_NODES(NODES_CONNECTIVITY(J,K),2)
	                VALUES(K,3)=COORD_NODES(NODES_CONNECTIVITY(J,K),3)
                ENDDO
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
                    DAUX=DSQRT((COORD_COLLOCPOINTS(I,1)-VALUES(K,1))**2+(COORD_COLLOCPOINTS(I,2)-VALUES(K,2))**2+(COORD_COLLOCPOINTS(I,3)-VALUES(K,3))**2)
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
!
		        ALLOCATE(DG(3,3*NEN),DH(3,3*NEN))
			    DG=0.D0
			    DH=0.D0
!
			    IF(DIST.EQ.0.D0) THEN
				    IF((II.EQ.1)) THEN
					    IF(EQ_TYPE(I).EQ.'S') THEN
					        CALL SINGULAR_ELEMENT_S(I,J,JJJ,DG,DH)
					    ELSE
						    CALL SINGULAR_ELEMENT_HP(I,J,JJJ,DG,DH)
					    ENDIF
				    ELSE
					    IF(ELEM_TYPE(J).EQ.3) THEN    
					        IF(EQ_TYPE(I).EQ.'S')THEN
					            CALL TRIANGULAR_ELEMENT_S(I,J,JJJ,DG,DH,N)
					        ELSE
					            CALL TRIANGULAR_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
					        ENDIF   
					    ELSE
                            IF(EQ_TYPE(I).EQ.'S')THEN
					            CALL QUADRILATERAL_ELEMENT_S(I,J,JJJ,DG,DH,N)
					        ELSE
					            CALL QUADRILATERAL_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
					        ENDIF   
					    ENDIF
				    ENDIF
			    ELSE
                    IF(ELEM_TYPE(J).EQ.3) THEN 
					    IF(EQ_TYPE(I).EQ.'S')THEN
					        CALL TRIANGULAR_ELEMENT_S(I,J,JJJ,DG,DH,N)
				        ELSE
					        CALL TRIANGULAR_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
					    ENDIF
			        ELSE
                        IF(EQ_TYPE(I).EQ.'S')THEN
					        CALL QUADRILATERAL_ELEMENT_S(I,J,JJJ,DG,DH,N)
					    ELSE
					        CALL QUADRILATERAL_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
					    ENDIF   
			        ENDIF
			    ENDIF
!
			    DO K=1,NEN
				    KK=COLLOCPOINTS_CONNECTIVITY(J,K)
!
			        G((3*I-2),(3*KK-2))=G((3*I-2),(3*KK-2))+DG(1,(3*K-2))
				    G((3*I-2),3*KK-1)=G((3*I-2),3*KK-1)+DG(1,3*K-1)
			        G((3*I-2),3*KK)=G((3*I-2),3*KK)+DG(1,3*K)
!				    
			        G((3*I-1),(3*KK-2))=G((3*I-1),(3*KK-2))+DG(2,(3*K-2))
				    G((3*I-1),3*KK-1)=G((3*I-1),3*KK-1)+DG(2,3*K-1)
			        G((3*I-1),3*KK)=G((3*I-1),3*KK)+DG(2,3*K)
!				    
			        G((3*I),(3*KK-2))=G((3*I),(3*KK-2))+DG(3,(3*K-2))
				    G((3*I),3*KK-1)=G((3*I),3*KK-1)+DG(3,3*K-1)
			        G((3*I),3*KK)=G((3*I),3*KK)+DG(3,3*K)
!
			        H((3*I-2),(3*KK-2))=H((3*I-2),(3*KK-2))+DH(1,(3*K-2))
				    H((3*I-2),3*KK-1)=H((3*I-2),3*KK-1)+DH(1,3*K-1)
			        H((3*I-2),3*KK)=H((3*I-2),3*KK)+DH(1,3*K)
!				    
			        H((3*I-1),(3*KK-2))=H((3*I-1),(3*KK-2))+DH(2,(3*K-2))
				    H((3*I-1),3*KK-1)=H((3*I-1),3*KK-1)+DH(2,3*K-1)
			        H((3*I-1),3*KK)=H((3*I-1),3*KK)+DH(2,3*K)
!				    
			        H((3*I),(3*KK-2))=H((3*I),(3*KK-2))+DH(3,(3*K-2))
				    H((3*I),3*KK-1)=H((3*I),3*KK-1)+DH(3,3*K-1)
			        H((3*I),3*KK)=H((3*I),3*KK)+DH(3,3*K)				    
			    ENDDO
!
		        DEALLOCATE(DG,DH)
            ENDIF
	    ENDDO 
	ENDDO
!$OMP END DO
!$OMP END PARALLEL  	
!
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   ADDING THE NEW COLLOCATION POINTS CONTRIBUTIONS
!------------------------------------------------------------------------------------------------------------------------------------------------------ 
! 
    CALL OMP_SET_NUM_THREADS(24)
!$OMP PARALLEL PRIVATE(I,J,K,II,JJ,KK,JJJ,NEN,N,DG,DH,DX,DY,DZ,DR,VALUES,LENGTH,AVERAGE_LENGTH,DISTANCE,DAUX)
!$OMP DO SCHEDULE(DYNAMIC)   	
  	DO I=N_COLLOCPOINTS-N_CRACKTIPCOLLOCPOINTS+1,N_COLLOCPOINTS
		JJ=0
		DO K=1,N_ELEM
		    NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			DO KK=1,NEN
				IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.I) THEN
					JJ=K
				ENDIF
			ENDDO
	    ENDDO
	    JJJ=ND(JJ)
!	    
		DO J=1,N_ELEM
		    IF(CRACKED_ELEM(J).EQ."UNCRACKED")THEN		    
			    IF(ND(JJ).EQ.ND(J)) THEN
				    II=0
				    NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
				    DO K=1,NEN
					    IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
						    II=1
					    ENDIF
				    ENDDO
!				
				    IF(II.EQ.0) THEN
		                DO K=1,NEN
			                DX=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),1)-COORD_COLLOCPOINTS(I,1)
			                DY=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),2)-COORD_COLLOCPOINTS(I,2)
	                        DZ=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),3)-COORD_COLLOCPOINTS(I,3)
			                DR=DSQRT(DX*DX+DY*DY+DZ*DZ)
			                IF(DR.LE.1.E-12) THEN
				                II=1
				            ENDIF
		                ENDDO
	                ENDIF
!
			        VALUES=0.D0
!
	                DO K=1,NEN
		                VALUES(K,1)=COORD_NODES(NODES_CONNECTIVITY(J,K),1)
		                VALUES(K,2)=COORD_NODES(NODES_CONNECTIVITY(J,K),2)
		                VALUES(K,3)=COORD_NODES(NODES_CONNECTIVITY(J,K),3)
	                ENDDO
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
                        DAUX=DSQRT((COORD_COLLOCPOINTS(I,1)-VALUES(K,1))**2+(COORD_COLLOCPOINTS(I,2)-VALUES(K,2))**2+(COORD_COLLOCPOINTS(I,3)-VALUES(K,3))**2)
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
!
			        ALLOCATE(DG(3,3*NEN),DH(3,3*NEN))
				    DG=0.D0
				    DH=0.D0
!
				    IF(DIST.EQ.0.D0) THEN
					    IF((II.EQ.1)) THEN
						    IF(EQ_TYPE(I).EQ.'S') THEN
						        CALL SINGULAR_ELEMENT_S(I,J,JJJ,DG,DH)
						    ELSE
							    CALL SINGULAR_ELEMENT_HP(I,J,JJJ,DG,DH)
						    ENDIF
					    ELSE
						    IF(ELEM_TYPE(J).EQ.3) THEN    
						        IF(EQ_TYPE(I).EQ.'S')THEN
						            CALL TRIANGULAR_ELEMENT_S(I,J,JJJ,DG,DH,N)
						        ELSE
						            CALL TRIANGULAR_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
						        ENDIF   
						    ELSE
                                IF(EQ_TYPE(I).EQ.'S')THEN
						            CALL QUADRILATERAL_ELEMENT_S(I,J,JJJ,DG,DH,N)
						        ELSE
						            CALL QUADRILATERAL_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
						        ENDIF   
						    ENDIF
					    ENDIF
				    ELSE
	                    IF(ELEM_TYPE(J).EQ.3) THEN 
						    IF(EQ_TYPE(I).EQ.'S')THEN
						        CALL TRIANGULAR_ELEMENT_S(I,J,JJJ,DG,DH,N)
					        ELSE
						        CALL TRIANGULAR_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
						    ENDIF
				        ELSE
                            IF(EQ_TYPE(I).EQ.'S')THEN
						        CALL QUADRILATERAL_ELEMENT_S(I,J,JJJ,DG,DH,N)
						    ELSE
						        CALL QUADRILATERAL_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
						    ENDIF   
				        ENDIF
				    ENDIF
!
				    DO K=1,NEN
					    KK=COLLOCPOINTS_CONNECTIVITY(J,K)
!
				        G((3*I-2),(3*KK-2))=G((3*I-2),(3*KK-2))+DG(1,(3*K-2))
					    G((3*I-2),3*KK-1)=G((3*I-2),3*KK-1)+DG(1,3*K-1)
				        G((3*I-2),3*KK)=G((3*I-2),3*KK)+DG(1,3*K)
!				    
				        G((3*I-1),(3*KK-2))=G((3*I-1),(3*KK-2))+DG(2,(3*K-2))
					    G((3*I-1),3*KK-1)=G((3*I-1),3*KK-1)+DG(2,3*K-1)
				        G((3*I-1),3*KK)=G((3*I-1),3*KK)+DG(2,3*K)
!				    
				        G((3*I),(3*KK-2))=G((3*I),(3*KK-2))+DG(3,(3*K-2))
					    G((3*I),3*KK-1)=G((3*I),3*KK-1)+DG(3,3*K-1)
				        G((3*I),3*KK)=G((3*I),3*KK)+DG(3,3*K)
!
				        H((3*I-2),(3*KK-2))=H((3*I-2),(3*KK-2))+DH(1,(3*K-2))
					    H((3*I-2),3*KK-1)=H((3*I-2),3*KK-1)+DH(1,3*K-1)
				        H((3*I-2),3*KK)=H((3*I-2),3*KK)+DH(1,3*K)
!				    
				        H((3*I-1),(3*KK-2))=H((3*I-1),(3*KK-2))+DH(2,(3*K-2))
					    H((3*I-1),3*KK-1)=H((3*I-1),3*KK-1)+DH(2,3*K-1)
				        H((3*I-1),3*KK)=H((3*I-1),3*KK)+DH(2,3*K)
!				    
				        H((3*I),(3*KK-2))=H((3*I),(3*KK-2))+DH(3,(3*K-2))
					    H((3*I),3*KK-1)=H((3*I),3*KK-1)+DH(3,3*K-1)
				        H((3*I),3*KK)=H((3*I),3*KK)+DH(3,3*K)				    
				    ENDDO
!
			        DEALLOCATE(DG,DH)
	            ENDIF
	        ENDIF
	    ENDDO 
	ENDDO
!$OMP END DO
!$OMP END PARALLEL 	
!
    OPEN(13,file='Output_data\SOMA_CRACKPROP.ogl',status='unknown')  
	SOMA_H=0.D0
	II=0
    DO I=1,N_COLLOCPOINTS
        IF(DUAL_BEM(I).EQ."B")THEN 
            II=II+1   
            DO J=1,3*N_COLLOCPOINTS
                SOMA_H(3*II-2)=SOMA_H(3*II-2)+H(3*I-2,J)
                SOMA_H(3*II-1)=SOMA_H(3*II-1)+H(3*I-1,J)
                SOMA_H(3*II)=SOMA_H(3*II)+H(3*I,J)                
            ENDDO
            WRITE(13,*)II,SOMA_H(3*II-2),SOMA_H(3*II-1),SOMA_H(3*II)
        ENDIF    
    ENDDO
    CLOSE(13)   
	SOMA_H=0.D0
	II=0
    DO I=1,N_COLLOCPOINTS
        IF(DUAL_BEM(I).EQ."S")THEN 
            II=II+1   
            DO J=1,3*N_COLLOCPOINTS
                SOMA_H(3*II-2)=SOMA_H(3*II-2)+H(3*I-2,J)
                SOMA_H(3*II-1)=SOMA_H(3*II-1)+H(3*I-1,J)
                SOMA_H(3*II)=SOMA_H(3*II)+H(3*I,J)                
            ENDDO
        ENDIF    
    ENDDO  
	SOMA_H=0.D0
	II=0
    DO I=1,N_COLLOCPOINTS
        IF(DUAL_BEM(I).EQ."H")THEN 
            II=II+1   
            DO J=1,3*N_COLLOCPOINTS
                SOMA_H(3*II-2)=SOMA_H(3*II-2)+H(3*I-2,J)
                SOMA_H(3*II-1)=SOMA_H(3*II-1)+H(3*I-1,J)
                SOMA_H(3*II)=SOMA_H(3*II)+H(3*I,J)                
            ENDDO                                    
        ENDIF    
    ENDDO
!------------------------------------------------------------------------------------------------------------------------------------------------------
!		INCLUDING THE ADDTIONAL FREE TERMS AT THE DIAGONAL H AND G
!------------------------------------------------------------------------------------------------------------------------------------------------------
	DO I=N_COLLOCPOINTS-N_CRACKTIPCOLLOCPOINTS+1,N_COLLOCPOINTS
	    IF(EQ_TYPE(I).EQ.'S') THEN
		    C=0.D0
	        DO J=1,N_COLLOCPOINTS
			    C(1,1)=C(1,1)-H(3*I-2,3*J-2)
			    C(1,2)=C(1,2)-H(3*I-2,3*J-1)
				C(1,3)=C(1,3)-H(3*I-2,3*J)
				C(2,1)=C(2,1)-H(3*I-1,3*J-2)
				C(2,2)=C(2,2)-H(3*I-1,3*J-1)
				C(2,3)=C(2,3)-H(3*I-1,3*J)
				C(3,1)=C(3,1)-H(3*I,3*J-2)
				C(3,2)=C(3,2)-H(3*I,3*J-1)
				C(3,3)=C(3,3)-H(3*I,3*J)
			ENDDO
!
            IF(DUAL_BEM(I).EQ."B")THEN
				H(3*I-2,3*I-2)=H(3*I-2,3*I-2)+C(1,1)
				H(3*I-2,3*I-1)=H(3*I-2,3*I-1)+C(1,2)
				H(3*I-2,3*I)=H(3*I-2,3*I)+C(1,3)
				H(3*I-1,3*I-2)=H(3*I-1,3*I-2)+C(2,1)
				H(3*I-1,3*I-1)=H(3*I-1,3*I-1)+C(2,2)
				H(3*I-1,3*I)=H(3*I-1,3*I)+C(2,3)
				H(3*I,3*I-2)=H(3*I,3*I-2)+C(3,1)
				H(3*I,3*I-1)=H(3*I,3*I-1)+C(3,2)
				H(3*I,3*I)=H(3*I,3*I)+C(3,3)
		    ELSE 
		        DO J=1,N_COLLOCPOINTS
		            IF(DUAL_BEM(J).EQ."H")THEN
				        DX=COORD_COLLOCPOINTS(I,1)-COORD_COLLOCPOINTS(J,1)
				        DY=COORD_COLLOCPOINTS(I,2)-COORD_COLLOCPOINTS(J,2)
				        DZ=COORD_COLLOCPOINTS(I,3)-COORD_COLLOCPOINTS(J,3)
				        DR=DSQRT(DX*DX+DY*DY+DZ*DZ)
				        IF(DR.LE.1.E-12)THEN
				            K=J
				        ENDIF
				    ENDIF
				ENDDO
				H(3*I-2,3*I-2)=H(3*I-2,3*I-2)+C(1,1)/2
				H(3*I-2,3*I-1)=H(3*I-2,3*I-1)+C(1,2)/2
				H(3*I-2,3*I)=H(3*I-2,3*I)+C(1,3)/2
				H(3*I-1,3*I-2)=H(3*I-1,3*I-2)+C(2,1)/2
				H(3*I-1,3*I-1)=H(3*I-1,3*I-1)+C(2,2)/2
				H(3*I-1,3*I)=H(3*I-1,3*I)+C(2,3)/2
				H(3*I,3*I-2)=H(3*I,3*I-2)+C(3,1)/2
				H(3*I,3*I-1)=H(3*I,3*I-1)+C(3,2)/2
				H(3*I,3*I)=H(3*I,3*I)+C(3,3)/2
				         
				H(3*I-2,3*K-2)=H(3*I-2,3*K-2)+C(1,1)/2
				H(3*I-2,3*K-1)=H(3*I-2,3*K-1)+C(1,2)/2
				H(3*I-2,3*K)=H(3*I-2,3*K)+C(1,3)/2
				H(3*I-1,3*K-2)=H(3*I-1,3*K-2)+C(2,1)/2
				H(3*I-1,3*K-1)=H(3*I-1,3*K-1)+C(2,2)/2
				H(3*I-1,3*K)=H(3*I-1,3*K)+C(2,3)/2
				H(3*I,3*K-2)=H(3*I,3*K-2)+C(3,1)/2
				H(3*I,3*K-1)=H(3*I,3*K-1)+C(3,2)/2
				H(3*I,3*K)=H(3*I,3*K)+C(3,3)/2
		    ENDIF
	    ELSE
			IF(DUAL_BEM(I).EQ."B")THEN
				G(3*I-2,3*I-2)=G(3*I-2,3*I-2)-0.5D0
				G(3*I-1,3*I-1)=G(3*I-1,3*I-1)-0.5D0
				G(3*I,3*I)=G(3*I,3*I)-0.5D0
		    ELSE
		        DO J=1,N_COLLOCPOINTS
				    IF(DUAL_BEM(J).EQ."S")THEN
				        DX=COORD_COLLOCPOINTS(I,1)-COORD_COLLOCPOINTS(J,1)
				        DY=COORD_COLLOCPOINTS(I,2)-COORD_COLLOCPOINTS(J,2)
				        DZ=COORD_COLLOCPOINTS(I,3)-COORD_COLLOCPOINTS(J,3)
				        DR=DSQRT(DX*DX+DY*DY+DZ*DZ)
				        IF(DR.LE.1.E-12)THEN
				            K=J
				        ENDIF
				    ENDIF
				ENDDO
 				G(3*I-2,3*I-2)=G(3*I-2,3*I-2)-0.50D0
 				G(3*I-1,3*I-1)=G(3*I-1,3*I-1)-0.50D0
 				G(3*I,3*I)=G(3*I,3*I)-0.50D0
				         
 				G(3*I-2,3*K-2)=G(3*I-2,3*K-2)+0.5D0
 				G(3*I-1,3*K-1)=G(3*I-1,3*K-1)+0.5D0
 				G(3*I,3*K)=G(3*I,3*K)+0.5D0
	        ENDIF
!				
		    C=0.D0
			DO J=1,N_COLLOCPOINTS
			    C(1,1)=C(1,1)-H(3*I-2,3*J-2)
			    C(1,2)=C(1,2)-H(3*I-2,3*J-1)
				C(1,3)=C(1,3)-H(3*I-2,3*J)
				C(2,1)=C(2,1)-H(3*I-1,3*J-2)
				C(2,2)=C(2,2)-H(3*I-1,3*J-1)
				C(2,3)=C(2,3)-H(3*I-1,3*J)
				C(3,1)=C(3,1)-H(3*I,3*J-2)
				C(3,2)=C(3,2)-H(3*I,3*J-1)
				C(3,3)=C(3,3)-H(3*I,3*J)
		    ENDDO
!
			H(3*I-2,3*I-2)=H(3*I-2,3*I-2)+C(1,1)
			H(3*I-2,3*I-1)=H(3*I-2,3*I-1)+C(1,2)
			H(3*I-2,3*I)=H(3*I-2,3*I)+C(1,3)
			H(3*I-1,3*I-2)=H(3*I-1,3*I-2)+C(2,1)
			H(3*I-1,3*I-1)=H(3*I-1,3*I-1)+C(2,2)
			H(3*I-1,3*I)=H(3*I-1,3*I)+C(2,3)
			H(3*I,3*I-2)=H(3*I,3*I-2)+C(3,1)
			H(3*I,3*I-1)=H(3*I,3*I-1)+C(3,2)
			H(3*I,3*I)=H(3*I,3*I)+C(3,3)
		ENDIF
    ENDDO        	 
!
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   REMOVING THE CHANGEABLE ELEMENTS CONTRIBUTIONS
!------------------------------------------------------------------------------------------------------------------------------------------------------ 
!
    CALL OMP_SET_NUM_THREADS(24)
!$OMP PARALLEL PRIVATE(I,J,K,II,JJ,KK,JJJ,NEN,N,DG,DH,DX,DY,DZ,DR,VALUES,LENGTH,AVERAGE_LENGTH,DISTANCE,DAUX)
!$OMP DO SCHEDULE(DYNAMIC)    
  	DO I=1,N_COLLOCPOINTS
		JJ=0
		DO K=1,N_ELEM
		    NEN=ELEM_TYPE(K)*ORDER_ELEM(K)+(ELEM_TYPE(K)-3)*(ORDER_ELEM(K)-1)*POL_FAMILY(K)
			DO KK=1,NEN
				IF(COLLOCPOINTS_CONNECTIVITY(K,KK).EQ.I) THEN
					JJ=K
				ENDIF
			ENDDO
	    ENDDO
	    JJJ=ND(JJ)
!	    
		DO J=1,N_ELEM
		    IF(CHANGEABLE_ELEM(J).EQ."CHANGEABLE")THEN	    
		        IF(ND(JJ).EQ.ND(J)) THEN
			        II=0
			        NEN=ELEM_TYPE(J)*ORDER_ELEM(J)+(ELEM_TYPE(J)-3)*(ORDER_ELEM(J)-1)*POL_FAMILY(J)
			        DO K=1,NEN
				        IF(COLLOCPOINTS_CONNECTIVITY(J,K).EQ.I) THEN
					        II=1
				        ENDIF
			        ENDDO
    !				
			        IF(II.EQ.0) THEN
	                    DO K=1,NEN
		                    DX=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),1)-COORD_COLLOCPOINTS(I,1)
		                    DY=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),2)-COORD_COLLOCPOINTS(I,2)
                            DZ=COORD_COLLOCPOINTS(COLLOCPOINTS_CONNECTIVITY(J,K),3)-COORD_COLLOCPOINTS(I,3)
		                    DR=DSQRT(DX*DX+DY*DY+DZ*DZ)
		                    IF(DR.LE.1.E-12) THEN
			                    II=1
			                ENDIF
	                    ENDDO
                    ENDIF
    !
		            VALUES=0.D0
    !
                    DO K=1,NEN
	                    VALUES(K,1)=COORD_NODES(NODES_CONNECTIVITY(J,K),1)
	                    VALUES(K,2)=COORD_NODES(NODES_CONNECTIVITY(J,K),2)
	                    VALUES(K,3)=COORD_NODES(NODES_CONNECTIVITY(J,K),3)
                    ENDDO
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
                        DAUX=DSQRT((COORD_COLLOCPOINTS(I,1)-VALUES(K,1))**2+(COORD_COLLOCPOINTS(I,2)-VALUES(K,2))**2+(COORD_COLLOCPOINTS(I,3)-VALUES(K,3))**2)
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
    !
		            ALLOCATE(DG(3,3*NEN),DH(3,3*NEN))
			        DG=0.D0
			        DH=0.D0
    !
			        IF(DIST.EQ.0.D0) THEN
				        IF((II.EQ.1)) THEN
					        IF(EQ_TYPE(I).EQ.'S') THEN
					            CALL SINGULAR_ELEMENT_S(I,J,JJJ,DG,DH)
					        ELSE
						        CALL SINGULAR_ELEMENT_HP(I,J,JJJ,DG,DH)
					        ENDIF
				        ELSE
					        IF(ELEM_TYPE(J).EQ.3) THEN    
					            IF(EQ_TYPE(I).EQ.'S')THEN
					                CALL TRIANGULAR_ELEMENT_S(I,J,JJJ,DG,DH,N)
					            ELSE
					                CALL TRIANGULAR_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
					            ENDIF   
					        ELSE
                                IF(EQ_TYPE(I).EQ.'S')THEN
					                CALL QUADRILATERAL_ELEMENT_S(I,J,JJJ,DG,DH,N)
					            ELSE
					                CALL QUADRILATERAL_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
					            ENDIF   
					        ENDIF
				        ENDIF
			        ELSE
                        IF(ELEM_TYPE(J).EQ.3) THEN 
					        IF(EQ_TYPE(I).EQ.'S')THEN
					            CALL TRIANGULAR_ELEMENT_S(I,J,JJJ,DG,DH,N)
				            ELSE
					            CALL TRIANGULAR_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
					        ENDIF
			            ELSE
                            IF(EQ_TYPE(I).EQ.'S')THEN
					            CALL QUADRILATERAL_ELEMENT_S(I,J,JJJ,DG,DH,N)
					        ELSE
					            CALL QUADRILATERAL_ELEMENT_HP(I,J,JJ,JJJ,DG,DH,N) 
					        ENDIF   
			            ENDIF
			        ENDIF
    !
			        DO K=1,NEN
				        KK=COLLOCPOINTS_CONNECTIVITY(J,K)
    !
			            G((3*I-2),(3*KK-2))=G((3*I-2),(3*KK-2))-DG(1,(3*K-2))
				        G((3*I-2),3*KK-1)=G((3*I-2),3*KK-1)-DG(1,3*K-1)
			            G((3*I-2),3*KK)=G((3*I-2),3*KK)-DG(1,3*K)
    !				    
			            G((3*I-1),(3*KK-2))=G((3*I-1),(3*KK-2))-DG(2,(3*K-2))
				        G((3*I-1),3*KK-1)=G((3*I-1),3*KK-1)-DG(2,3*K-1)
			            G((3*I-1),3*KK)=G((3*I-1),3*KK)-DG(2,3*K)
    !				    
			            G((3*I),(3*KK-2))=G((3*I),(3*KK-2))-DG(3,(3*K-2))
				        G((3*I),3*KK-1)=G((3*I),3*KK-1)-DG(3,3*K-1)
			            G((3*I),3*KK)=G((3*I),3*KK)-DG(3,3*K)
    !
			            H((3*I-2),(3*KK-2))=H((3*I-2),(3*KK-2))-DH(1,(3*K-2))
				        H((3*I-2),3*KK-1)=H((3*I-2),3*KK-1)-DH(1,3*K-1)
			            H((3*I-2),3*KK)=H((3*I-2),3*KK)-DH(1,3*K)
    !				    
			            H((3*I-1),(3*KK-2))=H((3*I-1),(3*KK-2))-DH(2,(3*K-2))
				        H((3*I-1),3*KK-1)=H((3*I-1),3*KK-1)-DH(2,3*K-1)
			            H((3*I-1),3*KK)=H((3*I-1),3*KK)-DH(2,3*K)
    !				    
			            H((3*I),(3*KK-2))=H((3*I),(3*KK-2))-DH(3,(3*K-2))
				        H((3*I),3*KK-1)=H((3*I),3*KK-1)-DH(3,3*K-1)
			            H((3*I),3*KK)=H((3*I),3*KK)-DH(3,3*K)				    
			        ENDDO
    !
		            DEALLOCATE(DG,DH)
                ENDIF
            ENDIF
	    ENDDO 
	ENDDO
!$OMP END DO
!$OMP END PARALLEL 	
!
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   REMOVING THE CHANGEABLE COLLOCATION POINTS CONTRIBUTIONS
!------------------------------------------------------------------------------------------------------------------------------------------------------ 
!   
  	DO I=1,N_COLLOCPOINTS
        IF(CHANGEABLE_COLLOCPOINTS(I).EQ."CHANGEABLE")THEN
             H((3*I-2),:)=0.D0
             H((3*I-1),:)=0.D0
             H((3*I),:)=0.D0  
!
             G((3*I-2),:)=0.D0
             G((3*I-1),:)=0.D0
             G((3*I),:)=0.D0                         
        ENDIF
	ENDDO	
!	   
    END SUBROUTINE