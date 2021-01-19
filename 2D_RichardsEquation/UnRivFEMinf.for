!	NUMERICAL SCHEME TO SOLVE THE RICHARDS EQUATION (2-D)
!	USING INFINITE ELEMENTS	IN THE FINITE ELEMENT ANALYSIS
!
!	**********************************************************************
!	This code simulates the stream seepage to an infinite phreatic aquifer 
!	by implementing the infinite elements in the finite element analysis. 
!	The code has been developed by Ch. Doulgeris (chdoulgeris@gmail.com)
!	and is freely available for downloading and further modification
!	at www.lri.swri.gr and https://github.com/chdoulgeris/InfElem.
!	**********************************************************************

	CHARACTER METHOD*8,METHOD_INF*8
	DOUBLE PRECISION,ALLOCATABLE:: PHI(:),PHIOLD(:)
	DOUBLE PRECISION,ALLOCATABLE:: PHINEW(:),PHIEST(:),PHALF(:)
	DOUBLE PRECISION,ALLOCATABLE:: X(:),Y(:)
	DOUBLE PRECISION,ALLOCATABLE:: AK(:,:),BK(:)
	DOUBLE PRECISION,ALLOCATABLE:: AKQ(:,:),BKQ(:),BKQ1(:),QNFIX(:)
	DOUBLE PRECISION SF,SFM,B,BM
	
	ALLOCATABLE NOD(:,:),AREA(:),NNX(:),YY(:),IFLOW(:),QL(:),Q(:)
	DIMENSION RW(10),RJ(10),RWINF(10),RJINF(10),RW1(10),RJ1(10)
	DIMENSION XY(4,2),PG(4),CKXY(2,2)
	DIMENSION SF(4),B(2,4),SFM(4),BM(2,4)
	DIMENSION B1(4,2),AK1(4,4),AK2(4,4),BK1(4)

	CALL CPU_TIME(TIMEBEGIN)

!	FILES
	OPEN(7,FILE='input.grd')		!Input file - Mesh data and boundary condition
	OPEN(8,FILE='UnRiv1.dat')		!Input file - Hydraulic and other data
	OPEN(9,FILE='UnRiv.res')		!Writing the data

	OPEN(14,FILE='PHI1_UnRiv.res')	!Piezometric head in nodes
	OPEN(10,FILE='PHI2_UnRiv.res')	!Head and pressure in coordinates
	OPEN(11,FILE='VOL_UnRiv.res')	!Water volume in time
	OPEN(13,FILE='Q1_UnRiv.res')	!Discharge in nodes of predefined (known) )head
	OPEN(20,FILE='Q2_UnRiv.res')	!Total discharge at river


	OPEN(12,FILE='INF_UnRiv.res')	!Head at the edge of the area
	OPEN(31,FILE='TEST1_UnRiv.res')	!Dummy file
	OPEN(32,FILE='TEST2_UnRiv.res')	!Dummy file


!	INPUT AND WRITING OF DATA
!	Spatial data
	READ(7,*)
	READ(7,*)NELEM,NNP,NBAND,METHOD,TIMERENUM
	ALLOCATE(PHI(NNP),PHINEW(NNP),PHIEST(NNP),PHALF(NNP),PHIOLD(NNP))
	ALLOCATE(NOD(NELEM,4),X(NNP),Y(NNP),AREA(NELEM))
	ALLOCATE(AK(NNP,NBAND),BK(NNP))
	READ(7,*)
	DO I=1,NELEM
		READ(7,*)NOD(I,1),NOD(I,2),NOD(I,3),NOD(I,4)
	ENDDO
 	READ(7,*)
	DO I=1,NNP
		READ(7,*)X(I),Y(I)
	ENDDO
	READ(7,*)

!	Special nodes
	READ(7,*)NFIX
	ALLOCATE(NNX(NFIX),YY(NFIX))
	ALLOCATE(AKQ(NNP,NNP),BKQ(NNP),BKQ1(NNP),QNFIX(NNP))
	DO I=1,NFIX
		READ(7,*)NNX(I),YY(I)
	ENDDO
	READ(7,*)
	READ(7,*)NQ
	ALLOCATE(IFLOW(NQ),QL(NQ),Q(NQ))
	DO I=1,NQ
		READ(7,*)IFLOW(I),QL(I)
	ENDDO
	READ(7,*)
	READ(7,*)YLAY
	READ(7,*)
	READ(7,*)WOA

!	Initial condition
	READ(8,1000)HSTART,ERROR1,PHIVAR
	DO I=1,NNP
		PHI(I)=HSTART
	ENDDO

!	Initial condition in river
	DO I=1,NFIX
		PHI(NNX(I))=YY(I)
	ENDDO

!	Rain-irrigation
	READ(8,1010)RAIN
	DO I=1,NQ
		Q(I)=QL(I)*RAIN
	ENDDO

!	Simulation time and time step
	READ(8,1010)TSTOP
	READ(8,1010)DT

!	Number of iteration for convergence
	READ(8,1030)ITERMAX

!	Brooks-Corey parameters for each soil layer
	READ(8,1020)TSAT1,TRES1,CKSXX1,CKSYY1
	READ(8,*)HB1,CLM1,CN1
	READ(8,1020)TSAT2,TRES2,CKSXX2,CKSYY2
	READ(8,*)HB2,CLM2,CN2

!	Integration points and weighting factors
	READ(8,1040)NGP
	READ(8,*)
	READ(8,*)(RW(K),K=1,NGP)
	READ(8,*)
	READ(8,*)(RJ(K),K=1,NGP)

	READ(8,1040)NGPINF
	READ(8,*)
	READ(8,*)(RWINF(K),K=1,NGPINF)
	READ(8,*)
	READ(8,*)(RJINF(K),K=1,NGPINF)

!	Aproaches at the edge of the simulation area
	READ(8,1050)INFMETHOD
	READ(8,*)
	READ(8,*)RA
	READ(8,*)
	READ(8,*)XSTART,YSTART
	READ(8,*)
	READ(8,*)
	READ(8,*)Xdirec,Ydirec
 1000 FORMAT(///F10.5,2(///F10.7))
 1010 FORMAT(//F15.8)
 1020 FORMAT(/4(//F10.8)/)
 1030 FORMAT(//I5)
 1040 FORMAT(///I5)
 1050 FORMAT(//////I5)

!	Data writing
	WRITE(9,2000)NELEM,NNP,NBAND,METHOD
	WRITE(9,2010)(I,(NOD(I,J),J=1,4),I=1,NELEM)
	WRITE(9,2020)(I,X(I),Y(I),I=1,NNP)
	WRITE(9,2030)NFIX,(NNX(I),YY(I),I=1,NFIX)
	WRITE(9,2040)NQ,(IFLOW(I),QL(I),I=1,NQ)
	WRITE(9,2050)HSTART,ERROR1,PHIVAR,RAIN,TSTOP,DT,ITERMAX
	WRITE(9,2060)TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1
	WRITE(9,2070)TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2
	WRITE(9,2080)YLAY

	IF (INFMETHOD==0) THEN
	 WRITE(9,2085)NGP
	ELSE
	 WRITE(9,2086)NGP,NGPINF
	 IF (INFMETHOD==1) METHOD_INF='LINEAR'
	 IF (INFMETHOD==2) METHOD_INF='EXPONENTIAL'
	  WRITE(9,2090)METHOD_INF,RA,XSTART,YSTART,Xdirec,Ydirec
	ENDIF

	WRITE(10,2100)
 2000 FORMAT('OUTPUT FILE - RESULTS - UnRiv.res'/33('*')/
     $'ELEMENT NUMBER='I5/
     $'NODE NUMBER='I5/
     $'BANDWIDTH OF MESH='I5/
     $'RENUMBERING METHOD OF MESH = ALGORITHM OF 'A8)
 2010 FORMAT(/'ELEMENT',5X,'NODES OF ELEMENT'/
     $(2X,I5,5X,4(2X,I5)))
 2020 FORMAT(/'COORDINATES OF NODES'/3(4X,'I',8X,'X',8X,'Y',5X)/
     $3(I5,2(1X,F8.2),5X))
 2030 FORMAT(/'NODES OF KNOWN HEAD='1X,I3/
     $2X'NODE',5X,'HEAD'/(3X,I5,3X,F8.3))
 2040 FORMAT(/'NODES OF KNOWN DISCHARGE='1X,I3/
     $2X'NODE',3X,'FLOW LENGHT'/(3X,I5,5X,F8.3))
 2050 FORMAT(
     $/'Initial piezometric head(m)                  = ',F10.2
     $/'Convergence error in units of head = ',F10.7
     $/'Maximun head in each time step = ',F10.7
     $/'Irrigation or rainfall(m/hr)                = ',F10.7
     $/'Simulation time(hours)                      = ',F10.2     
     $/'Initial tiem step(hours)               = ',F10.7     
     $/'Number of iterations for convergence = ',I10)
 2060 FORMAT(//'SOIL PARAMETERS OF UPPER SOIL LAYER'
     $/'SATURATION MOISTURE = ',F8.5
     $/'RESIDENCE MOISTURE = ',F8.5
     $/'HYDRAULIC CONDUCTIVITY IN SATURATION'
     $/'IN HORIZONTAL(m/hr)  = ',F8.5
     $/'HYDRAULIC CONDUCTIVITY IN SATURATION'
     $/'IN VERTICAL(m/hr) = ',F8.5
     $/'PARAMETERS OF SOIL MOISTURE'
     $/'Hb1(m)=',F8.5,2X,'Ln1=',F8.5,2X,'N1=',F8.5)
 2070 FORMAT(//'SOIL PARAMETERS OF LOWER LAYER'
     $/'SATURATION MOISTURE                    = ',F8.5
     $/'RESIDENCE MOISTURE                = ',F8.5
     $/'HYDRAULIC CONDUCTIVITY IN SATURATION'
     $/'IN HORIZONTAL(m/hr)  = ',F8.5
     $/'HYDRAULIC CONDUCTIVITY IN SATURATION'
     $/'IN VERTICAL(m/hr) = ',F8.5
     $/'PARAMETERS OF SOIL MOISTURE'
     $/'Hb2(m)=',F8.5,2X,'Ln2=',F8.5,2X,'N2=',F8.5)	
 2080 FORMAT(//'LOWER LAYER DEPTH'/'YLAY(m)=',F8.5)
 2085 FORMAT(//'NUMBER OF INTEGRATION POINTS'/'NGP=',I2)
 2086 FORMAT(//'NUMBER OF INTEGRATION POINTS'
     $/'Finite elements, NGP    =',I2
     $/'Infinite elements,   NGPINF =',I2)
 2090 FORMAT(//'APROACHING METHOD AT THE EDGE OF THE SIMULATION AREA'
     $/'Selection of method, INFMETHOD = ',A8
     $/'Parameter, a         =',E9.2
     $/'Coordinate in axis -X , XSTART     =',E9.2
     $/'Coordinate in axis -Y , YSTART     =',E9.2
     $/'Direction in axis -X , XX         =',E9.2
     $/'Direction in axis -Y , YY         =',E9.2)
 2100 FORMAT('OUTPUT FILE - RESULTS - PHI2_UnRiv.res'/33('*')
     $/'HEAD IN NODES'
     $/'OF MESH')
!
!
!
!CALCULATION ROUTE
!
!	1.	Loop for convergence
!	2.	Loop for each element
!	3.1	Read of coordinates of element nodes
!	3.2	Calculation of pressure head of element nodes
!	3.3	Searching for the infinite elements
!	3.4	Solving without the infinite elements
!	3.5	Numerical integration depending of the element type
!	4.	Loop for each integration point
!	5.	Matrices calculation B*,Bm*,xizi,N,M,B,Bm,J - SHAPEFUN
!	6.	Calculation of C(h) and K(h) - BROOKS_COREY
!	7.	Loop end for each integration point
!	8.	Creation of element matrices
!	9.	Assemblage of elements
!	10.	Loop end for each element
!	11.1	Boundary condition of known head
!	11.2	Irrigation or rain
!	12.	Solving the equation system - GAUSS
!	13.	Checking of convergence criterion
!	14.	Loop end for convergence
!	14.1	Reduce the time step in case of non-convergence
!		and return to step 1 (GOTO)
!	14.2	Searching for the maximum head difference
!	14.3	Searching for the head difference at the edge of the simulation area
!		and potential exit from the program - βήμα 22 (GOTO)
!	15.	Calculation of discharge of known head
!	16.	Total water volume based on the discharge
!		in the nodes of known head and discharge
!	17.	Calculation of water volume in the flow area
!	18.	Writing of results for the current time step
!	18.1	Preperation for the next time step
!	19.	Check for simulation time
!		and exit from the program - βήμα 22 (GOTO)
!	20.	Increase of time step in case of fast convergence
!		for the next time step
!	21.	Return to step 1 (GOTO)
!	22.	Exit from the program

!	Calculation of element area
	DO NE=1,NELEM
		AREA(NE)=X(NOD(NE,3))-X(NOD(NE,2))+X(NOD(NE,4))-X(NOD(NE,1))
	 	AREA(NE)=0.5*AREA(NE)*(Y(NOD(NE,4))-Y(NOD(NE,3)))
	ENDDO
!	Calculation of the initial water volume 
	CALL VOLUME(PHI,VOLSTART,AREA,NELEM,NNP,
     $NOD,X,Y,WOA,NGP,RJ,RW,NGPINF,RJINF,RWINF,
     $INFMETHOD,RA,XSTART,YSTART,Xdirec,Ydirec,
     $YLAY,TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1,
     $TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2)
	VOLQ=VOLSTART
	VOLPHI=VOLSTART

	WRITE(11,*)' TIME            VOLQ           VOLPHI'
	WRITE(11,*)TIME,VOLQ,VOLPHI
	WRITE(*,*)'   TIME       STEP OF TIME          ITERATION'
	TIME=0.0
 5000	ITERCHECK=0
	XY=0
	PG=0
	AK=0
	BK=0
	AKQ=0
	BKQ=0
	AK1=0
	BK1=0
	PHINEW=0
	PHIEST=0
	PHALF=0

!	1.	Loop for convergence
	DO 999 ITER=1,ITERMAX
	IF (ITERCHECK==0) THEN	!ITERCHECK
	AK=0
	BK=0
	AKQ=0
	BKQ=0
!	2.	Loop for each element
	DO 998 NE=1,NELEM
	INFCHECK=0
	DO 100 I=1,4
!	3.1	Read of coordinates of element nodes
	XY(I,1)=X(NOD(NE,I))
 	XY(I,2)=Y(NOD(NE,I))
!	3.2	Calculation of pressure head of element nodes
	IF (ITER==1) THEN
	 IF (TIME==0.0) THEN
	  PHALF(NOD(NE,I))=PHI(NOD(NE,I))
	  PHIEST(NOD(NE,I))=PHI(NOD(NE,I))
	 ELSE
	  PHALF(NOD(NE,I))=PHI(NOD(NE,I))
     $  +(DT/(2*DTOLD))*(PHI(NOD(NE,I))-PHIOLD(NOD(NE,I)))
	  PHIEST(NOD(NE,I))=PHI(NOD(NE,I))
     $  +(DT/DTOLD)*(PHI(NOD(NE,I))-PHIOLD(NOD(NE,I)))
	 ENDIF
	ENDIF
	PG(I)=PHALF(NOD(NE,I))-Y(NOD(NE,I))
!	3.3 Searching for the infinite elements
	IF (XY(I,1)==WOA) INFCHECK=1
!	3.4 Solving without the infinite elements
	IF (INFMETHOD==0) INFCHECK=0
 100	CONTINUE
	
	AK1=0
	AK2=0
	NGP1=0
	RJ1=0
	RW1=0

!	3.5 Numerical integration depending of the element type
	IF (INFCHECK==0) THEN
	 NGP1=NGP
	 RJ1=RJ
	 RW1=RW
	ELSE
	 NGP1=NGPINF
	 RJ1=RJINF
	 RW1=RWINF
	ENDIF
!	4.	Loop for each integration point
	DO IN=1,NGP1
	DO JN=1,NGP
!	5.	Matrices calculation B*,Bm*,xizi,N,M,B,Bm,J - SHAPEFUN
	CALL SHAPEFUN(IN,JN,RJ1,RJ,XY,DJ,SF,B,
     $INFCHECK,SFM,BM,INFMETHOD,RA,XSTART,YSTART,Xdirec,Ydirec)

!	6.	Calculation of C(h) και K(h) - BROOKS_COREY
	CALL BROOKS_COREY(YLAY,SF,SFM,XY,PG,CKXY,C,WET,
     $TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1,
     $TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2)

	B1=MATMUL(TRANSPOSE(BM),CKXY)	
	IF (INFCHECK==0) AK1=AK1+RW1(IN)*RW(JN)*DJ*DT*MATMUL(B1,BM)
	IF (INFCHECK==1) AK1=AK1+0.5*RW1(IN)*RW(JN)*DJ*DT*MATMUL(B1,BM)
	
	DO I=1,4					
		IF (INFCHECK==0) AK2(I,I)=AK2(I,I)+RW1(IN)*RW(JN)*C*DJ*SFM(I)
		IF (INFCHECK==1) AK2(I,I)=AK2(I,I)+0.5*RW1(IN)*RW(JN)*C*DJ*SFM(I)
	ENDDO

	ENDDO
	ENDDO
!	7.	Loop end for each integration point

!	8.	Creation of element matrices
	AK1=AK1+AK2
	DO I=1,4
		BK1(I)=0
		DO J=1,4
			BK1(I)=BK1(I)+AK2(I,J)*PHI(NOD(NE,J))
		ENDDO
	ENDDO

!	9.	Assemblage of elements
	DO I=1,4
	BK(NOD(NE,I))=BK(NOD(NE,I))+BK1(I)
	DO J=1,4
	IF(NOD(NE,I)<=NOD(NE,J)) THEN		!SYMETRY
	 AK(NOD(NE,I),NOD(NE,J)-NOD(NE,I)+1)=
     $ AK(NOD(NE,I),NOD(NE,J)-NOD(NE,I)+1)+AK1(I,J)
	ENDIF
	ENDDO
	ENDDO



	DO I=1,4
		BKQ(NOD(NE,I))=BKQ(NOD(NE,I))+BK1(I)
		DO J=1,4
			AKQ(NOD(NE,I),NOD(NE,J))=AKQ(NOD(NE,I),NOD(NE,J))+AK1(I,J)
		ENDDO
	ENDDO

 998	CONTINUE
!	10.	Loop end for each element

!	11.1	Boundary condition of known head
	DO I=1,NNP
		DO J=1,NFIX
			IF(NNX(J)==I) THEN
			 AK(I,1)=AK(I,1)*10.**30
			 BK(I)=YY(J)*AK(I,1)
			ENDIF
		ENDDO
	ENDDO

!	11.2	Irrigation or rain
	VOL2=0.
	DO I=1,NNP
		DO J=1,NQ
			IF(IFLOW(J)==I) THEN
			 BK(I)=BK(I)+Q(J)*DT
			 VOL2=VOL2+Q(J)*DT
			ENDIF
		ENDDO
	ENDDO




!	------------------------------------------------
!	12.	Solving the equation system - GAUSS
!	------------------------------------------------
	CALL GAUSS(AK,BK,NNP,NBAND,PHINEW)




!	13.	Checking of convergence criterion
	DMAX1=0.0
	DO I=1,NNP		
		DIFF=ABS(PHINEW(I)-PHIEST(I))
		IF (DIFF>=DMAX1) THEN
		 DMAX1=DIFF
		 IMAX1=I
		ENDIF
	ENDDO
	IF (DMAX1>=ERROR1) THEN
	 PHALF=(PHI+PHINEW)/2.
	 PHIEST=PHINEW
	 ITERCHECK=0
	ELSE
	 ITERCHECK=1
	 ITER1=ITER
	ENDIF

	ENDIF				!ITERCHECK
 999	CONTINUE
!	14.	Loop end for convergence
!	14.1	Reduce the time step in case of non-convergence
!		and return to step 1 (GOTO)
	IF (ITERCHECK==0) THEN
	 DT=DT*0.5
	 WRITE(*,*)'DT=',DT,' DMAX1=',DMAX1,' ITER=',ITER
	 WRITE(14,2105)DT,DMAX1,IMAX1
	 GOTO 5000
	ENDIF
 2105	FORMAT(/39('*')/'TIME STEP REDUCED IN :',E10.3/
     $'DUE TO NO CONVERGENCE IN HEAD:',E10.3/
     $'IN NODE:',I10/39('*'))

!	14.2	Searching for the maximum head difference
!		and return to step 1 (GOTO)
	DMAX=0.0
	DO I=1,NNP		
		DIFF=ABS(PHINEW(I)-PHI(I))
		IF (DIFF>=DMAX) THEN
		 DMAX=DIFF
		 IMAX=I
		ENDIF
	ENDDO

	
	IF (DMAX>=PHIVAR) THEN
	 DT=DT*0.5
	 WRITE(*,*)'DT=',DT,' DMAX=',DMAX
	 WRITE(14,2106)DT,DMAX,IMAX
	 GOTO 5000
	ENDIF
 2106	FORMAT(/41('*')/'TIME STEP REDUCED IN:',E10.3/
     $'DUE TO HEAD DIFFERENCE:',E10.3/
     $'IN NODE:',I10/41('*'))

!	14.3	Searching for the head difference at the edge of the simulation area
!		and potential exit from the program - step 22 (GOTO)
	IF (INFMETHOD==0) THEN
	 DO I=1,NNP		
		IF (X(I)==WOA) THEN	
		 IF (ABS(HSTART-PHINEW(I))>=ERROR1) THEN
		  WRITE(*,*) 'PROGRAMME TERMINATED'
		  WRITE(*,*) 'Increase width of area OR use INFINITE elements'
		  WRITE(10,*) 'PROGRAMME TERMINATED'
		  WRITE(10,*) 'Increase width of area OR use INFINITE elements'
		  GOTO 6000
		 ENDIF
		ENDIF
	 ENDDO
	ENDIF
		
!	15.	Calculation of discharge of known head
	DO I=1,NNP
		BKQ1(I)=0.
		DO J=1,NNP
			BKQ1(I)=BKQ1(I)+AKQ(I,J)*PHINEW(J)
		ENDDO
	ENDDO
	QNFIX=(BKQ1-BKQ)/DT
	QSUM=0.
	DO I=1,NNP
		DO J=1,NFIX
			IF (NNX(J)==I) THEN
			 QSUM=QSUM+QNFIX(I)
			ENDIF
		ENDDO
	ENDDO

!	16.	Total water volume based on the discharge
!		in the nodes of known head and discharge
	VOLUMEQ=VOLQ+QSUM*DT+VOL2


!	17.	Calculation of water volume in the flow area
!		
	CALL VOLUME(PHINEW,VOLUMEPHI,AREA,NELEM,NNP,
     $NOD,X,Y,WOA,NGP,RJ,RW,NGPINF,RJINF,RWINF,
     $INFMETHOD,RA,XSTART,YSTART,Xdirec,Ydirec,
     $YLAY,TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1,
     $TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2)


!	18.	Writing of results for the current time step
	TIME=TIME+DT
	WRITE(*,*)TIME,DT,ITER1,DMAX
	IF (MOD(TIME,1.)==0.) THEN	
	WRITE(10,2110)TIME,DT,ITER1,DMAX,IMAX,DMAX1,IMAX1
	WRITE(12,2110)TIME,DT,ITER1,DMAX,IMAX,DMAX1,IMAX1
	WRITE(13,2110)TIME,DT,ITER1,DMAX,IMAX,DMAX1,IMAX1
	WRITE(14,2110)TIME,DT,ITER1,DMAX,IMAX,DMAX1,IMAX1
	WRITE(10,*)'      X(I)       Y(I)     PHI(I)      PG(I)'
	WRITE(12,*)'      X(I)       Y(I)     PHI(I)'
	WRITE(13,*)'      X(I)       Y(I)     PHI(I)'
	DO I=1,NNP
		WRITE(10,2121)X(I),Y(I),PHINEW(I),PHINEW(I)-Y(I)
		IF (X(I)==WOA) WRITE(12,2120)X(I),Y(I),PHINEW(I)
		DO J=1,NFIX
			IF (NNX(J)==I) WRITE(13,2130)X(I),Y(I),PHINEW(I),QNFIX(I)
		ENDDO
	ENDDO
	WRITE(14,2140)(I,PHINEW(I),I=1,NNP)
	WRITE(13,2150)QSUM,TIME
	ENDIF						
	WRITE(11,*)TIME,VOLUMEQ,VOLUMEPHI
	WRITE(20,*)TIME,QSUM

 2110 FORMAT(///16X,'TIME(hours)=',E10.3/
     $16X,'DT(hours)  =',E10.3/
     $16X,'ITER       =',I10/
     $8X,'Maximum head difference(m)=',E10.3,1X,'in node=',I5/
     $8X,'Maximum convrgence error(m)=',E10.3,1X,'in node=',I5/
     $73('-')/)
 2120 FORMAT(2(1X,F10.3),(1X,F10.3))
 2121 FORMAT(2(1X,F10.3),2(1X,F10.3))
 2130 FORMAT(2(1X,F10.3),2(1X,F10.5))
 2140 FORMAT(5(1X,I5,':',F6.3,2X))
 2150 FORMAT('DISCHARGE='1X,E8.3,2X,'IN TIME= ',F6.3)

!	18.1	Preperation for the next time step
	DTOLD=DT
	PHIOLD=PHI
	PHI=PHINEW
	VOLQ=VOLUMEQ
	VOLPHI=VOLUMEPHI

!	19.	Check for simulation time
!		and exit from the program - step 22 (GOTO)
	IF (TIME>=TSTOP) GOTO 6000

!	20.	Increase of time step in case of fast convergence
	DT1=DTOLD*PHIVAR/DMAX	
	DT2=DTOLD*1.2		
	DT3=AINT(TIME)+1-TIME	
	DT=MIN(DT1,DT2,DT3)	

!	21.	Return to step 1 (GOTO)
	GOTO 5000

!	22.	Exit from the program

 6000	CONTINUE


							
	CALL CPU_TIME(TIMEEND)
	TIMEEQ=TIMEEND-TIMEBEGIN
	TIMETOTAL=TIMERENUM+TIMEEQ
	WRITE(9,2200)TIMERENUM,TIMEEQ/60,TIMETOTAL/60
 2200 FORMAT(//'COMPUTATIONAL TIMES'/
     $'for the renumbering of mesh(seconds) ='F10.3/
     $'for the equation(minutes) ='F10.3/
     $'Total computational time(minutes) ='F10.3)
	WRITE(*,*)
	WRITE(*,*)"RESULTS ARE (mainly) IN FILES UnRiv.res"
	WRITE(*,*)
	WRITE(*,*)'EXECUTABLE TIME WAS',TIMEEQ,' SECONDS'
	WRITE(*,*)
	END


	SUBROUTINE VOLUME(PHI,VOLTOTAL,AREA,NELEM,NNP,
     $NOD,X,Y,WOA,NGP,RJ,RW,NGPINF,RJINF,RWINF,
     $INFMETHOD,RA,XSTART,YSTART,Xdirec,Ydirec,
     $YLAY,TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1,
     $TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2)
	DOUBLE PRECISION PHI,X,Y,SF,SFM,B,BM
	DIMENSION PHI(NNP),X(NNP),Y(NNP)
	DIMENSION SF(4),B(2,4),SFM(4),BM(2,4)
	DIMENSION NOD(NELEM,4),AREA(NELEM),WETELEM(NELEM)
	DIMENSION RW(10),RJ(10),RWINF(10),RJINF(10),RW1(10),RJ1(10)
	DIMENSION XY(4,2),PG(4),CKXY(2,2)

	VOLTOTAL=0.
	DO 5 NE=1,NELEM
	INFCHECK=0
	DO I=1,4
		XY(I,1)=X(NOD(NE,I))
 		XY(I,2)=Y(NOD(NE,I))
		PG(I)=PHI(NOD(NE,I))-Y(NOD(NE,I))
		IF (XY(I,1)==WOA) INFCHECK=1
		IF (INFMETHOD==0) INFCHECK=0
	ENDDO
	
	NGP1=0
	RJ1=0
	RW1=0
	WETELEM(NE)=0
	IF (INFCHECK==0) THEN
	 NGP1=NGP
	 RJ1=RJ
	 RW1=RW
	ELSE
	 NGP1=NGPINF
	 RJ1=RJINF
	 RW1=RWINF
	ENDIF

	
	DO IN=1,NGP1
	  DO JN=1,NGP
	    CALL SHAPEFUN(IN,JN,RJ1,RJ,XY,DJ,SF,B,INFCHECK,SFM,BM,
     $INFMETHOD,RA,XSTART,YSTART,Xdirec,Ydirec)
	    CALL BROOKS_COREY(YLAY,SF,SFM,XY,PG,CKXY,C,WET,TSAT1,TRES1,
     $CKSXX1,CKSYY1,HB1,CLM1,CN1,TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2)
	    IF (INFCHECK==0) WETELEM(NE)=WETELEM(NE)+RW1(IN)*RW(JN)*WET/4.
	    IF (INFCHECK==1) WETELEM(NE)=WETELEM(NE)+0.5*RW1(IN)*RW(JN)*WET/4.
	  ENDDO
	ENDDO
	VOLTOTAL=VOLTOTAL+WETELEM(NE)*AREA(NE)
  5	CONTINUE
	RETURN
	END



	SUBROUTINE SHAPEFUN(IN,JN,RJ,RJ2,XY,DJ,SF,B,
     $INFCHECK,SFM,BM,INFMETHOD,RA,XSTART,YSTART,XX,YY)
	DOUBLE PRECISION SF,SFM,BS,BSM,B,BM
	DIMENSION RJ(10),RJ2(10),SF(4),BS(2,4),B(2,4),XY(4,2),AJ(2,2)
	DIMENSION SFM(4),BSM(2,4),BM(2,4),SR(4)
	DIMENSION SF1(4),BS1(2,4),PHI(4),PG(4)

	IF (INFCHECK==0) THEN	!finite element
	 RIN1=1.-RJ(IN)
	 RIN2=1.+RJ(IN)
	 RJN1=1.-RJ(JN)
	 RJN2=1.+RJ(JN)

	 SF(1)=RIN1*RJN2*0.25
	 SF(2)=RIN1*RJN1*0.25
	 SF(3)=RIN2*RJN1*0.25
	 SF(4)=RIN2*RJN2*0.25

	 BS(1,1)=-RJN2*0.25
	 BS(1,2)=-RJN1*0.25
	 BS(1,3)=RJN1*0.25
	 BS(1,4)=RJN2*0.25
	 BS(2,1)=RIN1*0.25
	 BS(2,2)=-RIN1*0.25
	 BS(2,3)=-RIN2*0.25
	 BS(2,4)=RIN2*0.25

	 AJ=MATMUL(BS,XY)
	 SFM=SF
	 BSM=BS
	ENDIF					



	IF (INFCHECK==1) THEN	!infinite element
	 RIN1=1.-RJ(IN)
	 RIN2=1.+RJ(IN)
	 RJN1=1.-RJ2(JN)
	 RJN2=1.+RJ2(JN)

	 IF (RJ(IN)<=0) THEN		
	  SF(1)=-RJ(IN)*RJN2*0.5
	  SF(2)=-RJ(IN)*RJN1*0.5
	  SF(3)=RIN2*RJN1*0.5
	  SF(4)=RIN2*RJN2*0.5

	  BS(1,1)=-RJN2*0.5
	  BS(1,2)=-RJN1*0.5
	  BS(1,3)=RJN1*0.5
	  BS(1,4)=RJN2*0.5
	  BS(2,1)=-RJ(IN)*0.5
	  BS(2,2)=RJ(IN)*0.5
	  BS(2,3)=-RIN2*0.5
	  BS(2,4)=RIN2*0.5

	  AJ=MATMUL(BS,XY)
	  SFM=SF
	  BSM=BS
	 ENDIF					

	 IF (RJ(IN)>0) THEN		
	  SF(1)=-0.5*RJ(IN)*RJN2/RIN1
	  SF(2)=-0.5*RJ(IN)*RJN1/RIN1
	  SF(3)=0.5*RJN1/RIN1
	  SF(4)=0.5*RJN2/RIN1

	  BS(1,1)=-0.5*RJN2/RIN1**2
	  BS(1,2)=-0.5*RJN1/RIN1**2
	  BS(1,3)=0.5*RJN1/RIN1**2
	  BS(1,4)=0.5*RJN2/RIN1**2
	  BS(2,1)=-0.5*RJ(IN)/RIN1
	  BS(2,2)=0.5*RJ(IN)/RIN1
	  BS(2,3)=-0.5/RIN1
	  BS(2,4)=0.5/RIN1

	  SF1(1)=0.
	  SF1(2)=0.
	  SF1(3)=RJN1*0.5
	  SF1(4)=RJN2*0.5

	  BS1(1,1)=0.
	  BS1(1,2)=0.
	  BS1(1,3)=0.
	  BS1(1,4)=0.
	  BS1(2,1)=0.
	  BS1(2,2)=0.
	  BS1(2,3)=-0.5
	  BS1(2,4)=0.5

	  AJ=MATMUL(BS,XY)
	  SX=0.
	  SY=0.
	  DO I=1,4
		SX=SX+SF(I)*XX*(XY(I,1)-XSTART)
		SY=SY+SF(I)*YY*(XY(I,2)-YSTART)
		SR(I)=SQRT((XY(I,1)-XSTART)**2+(XY(I,2)-YSTART)**2)
	  ENDDO
	  SRS=SQRT(SX**2+SY**2)
	  SXX=SX*AJ(1,1)+SY*AJ(1,2)
	  SYY=SX*AJ(2,1)+SY*AJ(2,2)

	  DO I=1,4
		IF (INFMETHOD==1) THEN
		 SFM(I)=SF1(I)*(SR(I)/SRS)**RA
		 BSM(1,I)=(BS1(1,I)-RA*SXX*SF1(I)/SRS**2)*(SR(I)/SRS)**RA
		 BSM(2,I)=(BS1(2,I)-RA*SYY*SF1(I)/SRS**2)*(SR(I)/SRS)**RA
		ENDIF
		IF (INFMETHOD==2) THEN
		 SFM(I)=SF1(I)*EXP(RA*(SR(I)-SRS))
		 BSM(1,I)=(BS1(1,I)-RA*SXX*SF1(I)/SRS)*EXP(RA*(SR(I)-SRS))
		 BSM(2,I)=(BS1(2,I)-RA*SYY*SF1(I)/SRS)*EXP(RA*(SR(I)-SRS))
		ENDIF
	  ENDDO
	 ENDIF					
	ENDIF					

	DJ=AJ(1,1)*AJ(2,2)-AJ(1,2)*AJ(2,1)
	DUMMY=AJ(1,1)
	AJ(1,1)=AJ(2,2)/DJ
	AJ(2,2)=DUMMY/DJ
	AJ(1,2)=-AJ(1,2)/DJ
	AJ(2,1)=-AJ(2,1)/DJ
	DJ=ABS(DJ)
	B=MATMUL(AJ,BSM)
	BM=MATMUL(AJ,BSM)
	RETURN
	END
	


	SUBROUTINE BROOKS_COREY(YLAY,SF,SFM,XY,PG,CKXY,C,WET,
     $TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1,
     $TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2)
	DOUBLE PRECISION SF,SFM
	DIMENSION SF(4),SFM(4),XY(4,2),PG(4),CKXY(2,2)
	YNGP=0.0
  	PGNGP=0.0
	DO I=1,4
		YNGP1=SF(I)*XY(I,2)
		PGNGP1=SFM(I)*PG(I)
		YNGP=YNGP+YNGP1
  		PGNGP=PGNGP+PGNGP1
	ENDDO

	IF(YNGP>YLAY)THEN
	 IF(PGNGP>-HB1)THEN
	  WET=TSAT1
	  CKXX=CKSXX1
	  CKYY=CKSYY1
	  ! The specific moisture capacity is not considered zero
	  ! in the saturated flow region and rather is taken equal
	  ! to a particularly small value, i.e. C(h)=1e-5 m-1,
	  ! in order to avoid non-convergence and ensure
	  ! a numerical solution in any case. 
	  C=0.00001
	 ELSE
	  SE=ABS(HB1/PGNGP)
	  WET=TRES1+(TSAT1-TRES1)*SE**CLM1
	  C=CLM1*(TSAT1-TRES1)*SE**(CLM1+1.)/HB1
	  CKXX=CKSXX1*SE**CN1
	  CKYY=CKSYY1*SE**CN1
	 ENDIF
	ELSE
	 IF(PGNGP>-HB2)THEN
	  WET=TSAT2
	  CKXX=CKSXX2
	  CKYY=CKSYY2
	  ! The specific moisture capacity is not considered zero
	  ! in the saturated flow region and rather is taken equal
	  ! to a particularly small value, i.e. C(h)=1e-5 m-1,
	  ! in order to avoid non-convergence and ensure
	  ! a numerical solution in any case. 
	  C=0.00001
	 ELSE
	  SE=ABS(HB2/PGNGP)
	  WET=TRES2+(TSAT2-TRES2)*SE**CLM2
	  C=CLM2*(TSAT2-TRES2)*SE**(CLM2+1.)/HB2
	  CKXX=CKSXX2*SE**CN2
	  CKYY=CKSYY2*SE**CN2
	 ENDIF
	ENDIF
	CKXY(1,1)=CKXX
	CKXY(2,2)=CKYY
	CKXY(2,1)=0.
	CKXY(1,2)=0.
	RETURN
	END




	SUBROUTINE GAUSS(A,B,N,M,PHINEW)
	DOUBLE PRECISION A,B,PHINEW
	DIMENSION A(N,M),C(M),B(N),PHINEW(N)

	DO K=1,N-1
		B(K)=B(K)/A(K,1)
		DO J=2,M
			C(J)=A(K,J)
			A(K,J)=A(K,J)/A(K,1)
		ENDDO
		DO L=2,M
			I=K+L-1
			IF(N>=I)THEN
			 J=0
			 DO LL=L,M
			 	J=J+1
				A(I,J)=A(I,J)-C(L)*A(K,LL)
			 ENDDO
			 B(I)=B(I)-C(L)*B(K)
			ENDIF
		ENDDO
	ENDDO

	B(N)=B(N)/A(N,1)
	DO K=N-1,1,-1
		DO J=2,M
			L=K+J-1
			IF(N>=L) B(K)=B(K)-A(K,J)*B(L)
		ENDDO
	ENDDO

	DO I=1,N
		PHINEW(I)=B(I)
	ENDDO
	RETURN
	END
