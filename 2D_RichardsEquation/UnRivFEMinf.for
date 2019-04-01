!	������� ��� ����������� �������� ��� RICHARDS
!	�� �� ������ ��� ������������ ���������
!	��� �� ����� ���������� ���������
!
!	���������� ��� ������������
!NELEM		������� ��� ��������� ��� �������
!NNP			������� ��� ������ ��� �������
!NBAND		�����(BANDWIDTH) ��� �������
!
!PHIOLD(NNP)	������ ���� ����������� ������� ������
!PHI(NNP)		������ ���� ������� ������� ������
!PHINEW(NNP)	������������� ������ ���� ������� ������� ������
!PHIEST(NNP)	����������� ������ ���� ������� ������� ������
!PHALF(NNP)	������ ���� ��������� ������� ������
!AREA(NELEM)	������ ��� ���� ���������
!NOD(NELEM,NNP)	�������� ���� ������� ��� ���� ��������
!X(NNP),Y(NNP)	������������� ��� ������ ��� �������
!
!NFIX			������� ������ ��������������� �������
!NNX(NFIX)	������ ��. �������
!YY(NFIX)		���� ��. �������
!NQ			������� ������ ��������������� �������
!IFLOW(NQ)	������ ��. �������
!QL(NQ)		����� ���� ����� ������� ��. �������
!Q(NQ)		������������� ������ ����� ������� ��. �������
!	
!NGP			������ ����������� ����������� ���� ���� ����������
!RW(10)		������������ ���������� ������� ��. �����������
!RJ(10)		������������� ������� ��. �����������
!XY(4,2)		������������� ��� ������ ���� ���������
!AJ(2,2)		�� ��������� ������
!
!			���������� ����������
!SF(4)		����������� ������ ��� ���������
!BS(2,4)		��������� ��� ����������� ������
!B(2,4)		�� ������ �=AJ*BS
!
!			���������� ��� ����������� ����������
!SFM(4)		����������� ������ ��� ���������
!BSM(2,4)		��������� ��� ����������� ������
!BM(2,4)		�� ������ �M=AJ*BSM
!
!WOA			�� ������ ��� �������� �������
!YLAY			�� ����� ��� ���� �������� �������
!
!PG(4)		���� ������ ��� ������ ���� ���������
!TSAT1,TRES1	������� �������� ��� ������������ �������
!CKSXX1		��������� ����������� �������� ���� ��� ��������� ���������
!CKSYY1		��������� ����������� �������� ���� ��� ���������� ���������
!HB1,CLM1,CN1	�� ���������� (Hb,� ��� n) ���� ������� ��� BROOKS-KOREY
!CKXY(2,2)	������ ���������� ������������
!C			��������� ������������
!
!			��������� ������ ��� ���� ������������ �� ���� ��������
!B1(4,2),AK1(4,4),AK2(4,4),BK1(4)
!
!AK(NNP,NBAND),BK(NNP)	������ ������ ���� ��� ��������� ��� GAUSS
!
!QNFIX(NNP)		������������� ������ ����� ������� ��������������� �������
!				��������� ������ ��� ��� ���������� ��� �������
!AKQ(NNP,NBAND),BKQ(NNP),BKQ1(NNP)
!
!HSTART		� ������ ���� ��� ������� �� ����� ���� �������
!ERROR1		�� �������� ��������� ��� �������
!PHIVAR		� ������� �������� ��� ������� �� ���� ������� ����
!RAIN			������ �������� � �����������
!
!TIME			� ������� ������ ���� ���� ������������
!TSTOP		� ������� ������
!DT			�� ������ ������� ����
!ITERMAX		� ������� ��� ����������� ���������
!ITERCHECK	������� ��� ��� ������ ��� ���������
!
!INFMETHOD	������� ��� ��� ���������� ��� ����� ��� �������� �������
!INFCHECK		������� ��� ��� ������ ���������� ���������



C	������� ���������� ��� �������
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

C	������ ��� ����������������
	OPEN(7,FILE='input.grd')		!�������� ������� ��� ������� ��������
	OPEN(8,FILE='UnRiv1.dat')		!��������� ��� ����� ��������
	OPEN(9,FILE='UnRiv.res')		!�������� ��� ��������� ��� �����������

	OPEN(14,FILE='PHI1_UnRiv.res')	!����� ��� ������� �� ���� �����
	OPEN(10,FILE='PHI2_UnRiv.res')	!������ ��� ���� ������ ���� �������������
	OPEN(11,FILE='VOL_UnRiv.res')	!����� ����� �� �� �����
	OPEN(13,FILE='Q1_UnRiv.res')	!������ ����� ������� ������. �������
	OPEN(20,FILE='Q2_UnRiv.res')	!�������� ������ ��� ���� �� �� �����


	OPEN(12,FILE='INF_UnRiv.res')	!����� ������� ��� ����� ��� ��������
	OPEN(31,FILE='TEST1_UnRiv.res')	!������ ��� �������
	OPEN(32,FILE='TEST2_UnRiv.res')	!������ ��� �������


C	�������� ��� ������� ��� ���������
C	���������� ��������
	READ(7,*)
	READ(7,*)NELEM,NNP,NBAND,METHOD,TIMERENUM
	ALLOCATE(PHI(NNP),PHINEW(NNP),PHIEST(NNP),PHALF(NNP),PHIOLD(NNP))
	ALLOCATE(NOD(NELEM,4),X(NNP),Y(NNP),AREA(NELEM))
	ALLOCATE(AK(NNP,NBAND),BK(NNP))
	READ(7,*)
	DO 10 I=1,NELEM
 10	READ(7,*)NOD(I,1),NOD(I,2),NOD(I,3),NOD(I,4)
 	READ(7,*)
	DO 11 I=1,NNP
 11	READ(7,*)X(I),Y(I)
	READ(7,*)

C	��������������� ������ ��� �������
	READ(7,*)NFIX
	ALLOCATE(NNX(NFIX),YY(NFIX))
	ALLOCATE(AKQ(NNP,NNP),BKQ(NNP),BKQ1(NNP),QNFIX(NNP))
	DO 12 I=1,NFIX
 12	READ(7,*)NNX(I),YY(I)
	READ(7,*)
	READ(7,*)NQ
	ALLOCATE(IFLOW(NQ),QL(NQ),Q(NQ))
	DO 13 I=1,NQ
 13	READ(7,*)IFLOW(I),QL(I)
	READ(7,*)
	READ(7,*)YLAY
	READ(7,*)
	READ(7,*)WOA

C	�������� ������� ��������
	READ(8,1000)HSTART,ERROR1,PHIVAR
	DO 14 I=1,NNP
 14	PHI(I)=HSTART

C	�������� ������� �������� ���� ������� ��� �������
	DO 15 I=1,NFIX
 15	PHI(NNX(I))=YY(I)

C	������ ������ �����������-��������
	READ(8,1010)RAIN
	DO 16 I=1,NQ
 16	Q(I)=QL(I)*RAIN

C	������ ���� ��� ������ ��� ������� ���� ���������
	READ(8,1010)TSTOP
	READ(8,1010)DT

C	������� ����������� ���������
	READ(8,1030)ITERMAX

C	���������� Brooks-Corey ��� ���� ������� ������
	READ(8,1020)TSAT1,TRES1,CKSXX1,CKSYY1
	READ(8,*)HB1,CLM1,CN1
	READ(8,1020)TSAT2,TRES2,CKSXX2,CKSYY2
	READ(8,*)HB2,CLM2,CN2

C	������ ����������� ��� ������������ ����������
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

C	������� ����������� ��� ����� ��� �������� �������
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

C	������� ��� ��������� ��� �����������
	WRITE(9,2000)NELEM,NNP,NBAND,METHOD
	WRITE(9,2010)(I,(NOD(I,J),J=1,4),I=1,NELEM)
	WRITE(9,2020)(I,(X(I),Y(I)),I=1,NNP)
	WRITE(9,2030)NFIX,((NNX(I),YY(I)),I=1,NFIX)
	WRITE(9,2040)NQ,((IFLOW(I),QL(I)),I=1,NQ)
	WRITE(9,2050)HSTART,ERROR1,PHIVAR,RAIN,TSTOP,DT,ITERMAX
	WRITE(9,2060)TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1
	WRITE(9,2070)TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2
	WRITE(9,2080)YLAY

	IF (INFMETHOD==0) THEN
	 WRITE(9,2085)NGP
	ELSE
	 WRITE(9,2086)NGP,NGPINF
	 IF (INFMETHOD==1) METHOD_INF='��������'
	 IF (INFMETHOD==2) METHOD_INF='��������'
	  WRITE(9,2090)METHOD_INF,RA,XSTART,YSTART,Xdirec,Ydirec
	ENDIF

	WRITE(10,2100)
 2000 FORMAT('������ ������������� - UnRiv.res'/33('*')/
     $'������� ��� ��������� ��� ������� ='I5/
     $'������� ��� ������ ��� �������    ='I5/
     $'����� ��� �������                 ='I5/
     $'������� ����������� ��� �������   = ���������� ��� 'A8)
 2010 FORMAT(/'�/� ���������',5X,'������ ��� ���������'/
     $(2X,I5,5X,4(2X,I5)))
 2020 FORMAT(/'������������� ��� ������'/3(4X'I',8X'X',8X,'Y',5X)/
     $3(I5,2(1X,F8.2),5X))
 2030 FORMAT(/'������ ������� ������� ='1X,I3/
     $2X'������',5X,'������'/(3X,I5,3X,F8.3))
 2040 FORMAT(/'������ ������� ������� ='1X,I3/
     $2X'������',3X,'����� ����'/(3X,I5,5X,F8.3))
 2050 FORMAT(
     $/'������ ���� ���������� �������(m)                  = ',F10.2
     $/'������ ��������� �� ������� ��� ���������� ������� = ',F10.7
     $/'������� �������� ��� ������� �� ���� ������� ����  = ',F10.7
     $/'������ �������� � �����������(m/hr)                = ',F10.7
     $/'������ ���� ��� ������(hours)                      = ',F10.2     
     $/'������ ������� ���� ���������(hours)               = ',F10.7     
     $/'������� ����������� ���������                      = ',I10)
 2060 FORMAT(//'�������� ���������� ��� ���� �������'
     $/'������� ��������                    = ',F8.5
     $/'������������ �������                = ',F8.5
     $/'��������� ����������� ��������'
     $/'���� ��� ��������� ���������(m/hr)  = ',F8.5
     $/'��������� ����������� ��������'
     $/'���� ��� ���������� ���������(m/hr) = ',F8.5
     $/'���������� �������� ��������'
     $/'Hb1(m)=',F8.5,2X,'Ln1=',F8.5,2X,'N1=',F8.5)
 2070 FORMAT(//'�������� ���������� ��� ���� �������'
     $/'������� ��������                    = ',F8.5
     $/'������������ �������                = ',F8.5
     $/'��������� ����������� ��������'
     $/'���� ��� ��������� ���������(m/hr)  = ',F8.5
     $/'��������� ����������� ��������'
     $/'���� ��� ���������� ���������(m/hr) = ',F8.5
     $/'���������� �������� ��������'
     $/'Hb2(m)=',F8.5,2X,'Ln2=',F8.5,2X,'N2=',F8.5)	
 2080 FORMAT(//'����� ��� ���� �������'/'YLAY(m)=',F8.5)
 2085 FORMAT(//'������� ������� ����������� �����������'/'NGP=',I2)
 2086 FORMAT(//'������� ������� ����������� �����������'
     $/'������������ ���������, NGP    =',I2
     $/'���������� ���������,   NGPINF =',I2)
 2090 FORMAT(//'������� ����������� ��� ����� ��� �������� �������'
     $/'������� �������    , INFMETHOD = ',A8
     $/'���� ��� ����������, a         =',E9.2
     $/'���� ��� ����� -X , XSTART     =',E9.2
     $/'���� ��� ����� -Y , YSTART     =',E9.2
     $/'���� ��� ����� -X , XX         =',E9.2
     $/'���� ��� ����� -Y , YY         =',E9.2)
 2100 FORMAT('������ ������������� - PHI2_UnRiv.res'/33('*')
     $/'����� ������� ����� �������'
     $/'��� ������� ��� �������������')
!
!
!
!������ ��� �����������
!
C	1.	���� ������� ����������� ���������
C	2.	���� ������� ��� ���� ��������
C	3.1	�������� ��� ������������� ��� ������ ��� ���������
C	3.2	����������� ��� ����� ������ ��� ������ ��� ���������
C	3.3	������ ��� ���������� ���������
C	3.4	������� ����� �� ����� ���������� ���������
C	3.5	���������� ���������� ������� �� �� ����� ��� ���������
C	4.	���� ������� ��� ���� ������ �����������
C	5.	����������� ��� ������� B*,Bm*,xizi,N,M,B,Bm,J - SHAPEFUN
C	6.	����������� ��� C(h) ��� K(h) - BROOKS_COREY
C	7.	�������� ������� ��� ���� ������ �����������
C	8.	���������� ��� ������� ��� ���������	
C	9.	�������� ��� ���������
C	10.	�������� ������� ��� ���� ��������
C	11.1	�������� ������� �������� ��������������� �������
C	11.2	�������� ������� �������� � �����������
C	12.	������� ��� ���������� ��� ��������� - GAUSS
C	13.	������� ��� ��������� ���������
C	14.	�������� ��� ������� ���������
C	14.1	������ ��� �������� ������� �� ��������� �� ���������
C		��� ��������� ��� ���� 1 (GOTO)
C	14.2	������ ��� �������� ��������� ��� �������
C	14.3	������� ��� ��������� ��� ������� ��� ����� ��� ��������
C		��� ������ ������ ��� �� ��������� - ���� 22 (GOTO)
C	15.	����������� ��� ������� ����� ������� ��������������� �������
C	16.	������������� ����� ����� ���� ��� ����� ��� ������� �����
C		������� ��������������� ������� ��� ��������������� �������
C	17.	����������� ��� ��������� ��� ����� ���� ������� ����
C	18.	�������� ��� ������������� ��� ��������� �������� �������
C	18.1	������������ ��� ��� ������� ������� ������
C	19.	������� ��� ������� ����� ��� ������
C		��� ������ ��� �� ��������� - ���� 22 (GOTO)
C	20.	������ ��� �������� ������� �� ��������� �������� ���������
C		��� ��� ������� ������� ������
C	21.	��������� ��� ���� 1 (GOTO)
C	22.	������ ��� �� ���������

!	����������� ��� ������� ���� ���������
	DO 50 NE=1,NELEM
	AREA(NE)=X(NOD(NE,3))-X(NOD(NE,2))+X(NOD(NE,4))-X(NOD(NE,1))
 50 	AREA(NE)=0.5*AREA(NE)*(Y(NOD(NE,4))-Y(NOD(NE,3)))
!	����������� ��� ������ ������������� ����� �����
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

C	1.	���� ������� ����������� ���������
	DO 999 ITER=1,ITERMAX
	IF (ITERCHECK==0) THEN	!ITERCHECK
	AK=0
	BK=0
	AKQ=0
	BKQ=0
C	2.	���� ������� ��� ���� ��������
	DO 998 NE=1,NELEM
	INFCHECK=0
	DO 100 I=1,4
C	3.1	�������� ��� ������������� ��� ������ ��� ���������
	XY(I,1)=X(NOD(NE,I))
 	XY(I,2)=Y(NOD(NE,I))
C	3.2	����������� ��� ����� ������ ��� ������ ��� ���������
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
C	3.3 ������ ��� ���������� ���������
	IF (XY(I,1)==WOA) INFCHECK=1
C	3.4 ������� ����� �� ����� ���������� ���������
	IF (INFMETHOD==0) INFCHECK=0
 100	CONTINUE
	
	AK1=0
	AK2=0
	NGP1=0
	RJ1=0
	RW1=0

C	3.5 ���������� ���������� ������� �� �� ����� ��� ���������
	IF (INFCHECK==0) THEN
	 NGP1=NGP
	 RJ1=RJ
	 RW1=RW
	ELSE
	 NGP1=NGPINF
	 RJ1=RJINF
	 RW1=RWINF
	ENDIF
C	4.	���� ������� ��� ���� ������ �����������
	DO 997 IN=1,NGP1
	DO 997 JN=1,NGP
C	5.	����������� ��� ������� B*,Bm*,xizi,N,M,B,Bm,J - SHAPEFUN
	CALL SHAPEFUN(IN,JN,RJ1,RJ,XY,DJ,SF,B,
     $INFCHECK,SFM,BM,INFMETHOD,RA,XSTART,YSTART,Xdirec,Ydirec)

C	6.	����������� ��� C(h) ��� K(h) - BROOKS_COREY
	CALL BROOKS_COREY(YLAY,SF,SFM,XY,PG,CKXY,C,WET,
     $TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1,
     $TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2)

	B1=MATMUL(TRANSPOSE(BM),CKXY)	!����������� ��� ������� �������
	IF (INFCHECK==0) AK1=AK1+RW1(IN)*RW(JN)*DJ*DT*MATMUL(B1,BM)
	IF (INFCHECK==1) AK1=AK1+0.5*RW1(IN)*RW(JN)*DJ*DT*MATMUL(B1,BM)
	
	DO 150 I=1,4					!����������� ��� �������� �������
	IF (INFCHECK==0) AK2(I,I)=AK2(I,I)+RW1(IN)*RW(JN)*C*DJ*SFM(I)
	IF (INFCHECK==1) AK2(I,I)=AK2(I,I)+0.5*RW1(IN)*RW(JN)*C*DJ*SFM(I)
 150	CONTINUE

 997	CONTINUE
C	7.	�������� ������� ��� ���� ������ �����������

C	8.	���������� ��� ������� ��� ���������	
	AK1=AK1+AK2
	DO 200 I=1,4
	BK1(I)=0
	DO 200 J=1,4
 200	BK1(I)=BK1(I)+AK2(I,J)*PHI(NOD(NE,J))

C	9.	�������� ��� ���������
	DO 250 I=1,4
	BK(NOD(NE,I))=BK(NOD(NE,I))+BK1(I)
	DO 250 J=1,4
	IF(NOD(NE,I)<=NOD(NE,J)) THEN		!���������� ����� � ��������������
	 AK(NOD(NE,I),NOD(NE,J)-NOD(NE,I)+1)=
     $ AK(NOD(NE,I),NOD(NE,J)-NOD(NE,I)+1)+AK1(I,J)
	ENDIF
 250	CONTINUE

C	9.	���������� ��� ����������� ��� ����������� ��� ��� ����������
C		��� ������� ����� ������� ��������������� �������
	DO 251 I=1,4
	BKQ(NOD(NE,I))=BKQ(NOD(NE,I))+BK1(I)
	DO 251 J=1,4
	AKQ(NOD(NE,I),NOD(NE,J))=AKQ(NOD(NE,I),NOD(NE,J))+AK1(I,J)
 251	CONTINUE

 998	CONTINUE
C	10.	�������� ������� ��� ���� ��������

C	11.1	�������� ������� �������� ��������������� �������
	DO 300 I=1,NNP
	DO 300 J=1,NFIX
	IF(NNX(J)==I) THEN
	 AK(I,1)=AK(I,1)*10**30
	 BK(I)=YY(J)*AK(I,1)
	ENDIF
 300	CONTINUE

C	11.2	�������� ������� �������� � �����������
	VOL2=0.
	DO 350 I=1,NNP
	DO 350 J=1,NQ
	IF(IFLOW(J)==I) THEN
	 BK(I)=BK(I)+Q(J)*DT
	 VOL2=VOL2+Q(J)*DT
	ENDIF
 350	CONTINUE




C	------------------------------------------------
C	12.	������� ��� ���������� ��� ��������� - GAUSS
C	------------------------------------------------
	CALL GAUSS(AK,BK,NNP,NBAND,PHINEW)




C	13.	������� ��� ��������� ���������
	DMAX1=0.0
	DO 400 I=1,NNP		!������ ���� ��� ������
	DIFF=ABS(PHINEW(I)-PHIEST(I))
	IF (DIFF>=DMAX1) THEN
	 DMAX1=DIFF
	 IMAX1=I
	ENDIF
 400	CONTINUE
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
C	14.	�������� ��� ������� ���������

C	14.1	������ ��� �������� ������� �� ��������� �� ���������
C		��� ��������� ��� ���� 1 (GOTO)
	IF (ITERCHECK==0) THEN
	 DT=DT*0.5
	 WRITE(*,*)'DT=',DT,' DMAX1=',DMAX1,' ITER=',ITER
	 WRITE(14,2105)DT,DMAX1,IMAX1
	 GOTO 5000
	ENDIF
 2105	FORMAT(/39('*')/'�� ������� ���� �������� ��:',E10.3/
     $'���� �� ��������� ��� ������:',E10.3/
     $'���� �����:',I10/39('*'))

C	14.2	������ ��� �������� ��������� ��� �������
C		��� ������ ��������� ��� ���� 1 (GOTO)
	DMAX=0.0
	DO 410 I=1,NNP		!������ ���� ��� ������
	DIFF=ABS(PHINEW(I)-PHI(I))
	IF (DIFF>=DMAX) THEN
	 DMAX=DIFF
	 IMAX=I
	ENDIF
 410	CONTINUE
	!������� ��� �������� ������������� ��������� ��� �������
	!��� ������ ��� �������� �������
	IF (DMAX>=PHIVAR) THEN
	 DT=DT*0.5
	 WRITE(*,*)'DT=',DT,' DMAX=',DMAX
	 WRITE(14,2106)DT,DMAX,IMAX
	 GOTO 5000
	ENDIF
 2106	FORMAT(/41('*')/'�� ������� ���� �������� ��:',E10.3/
     $'���� ������� ��������� �������:',E10.3/
     $'���� �����:',I10/41('*'))

C	14.3	������� ��� ��������� ��� ������� ��� ����� ��� ��������
C		��� ������ ������ ��� �� ��������� - ���� 22 (GOTO)
	IF (INFMETHOD==0) THEN
	DO 430 I=1,NNP		!������ ���� ��� ������
	IF (X(I)==WOA) THEN	!������� ��� ������ ��� ����� ��� �������
	 IF (ABS(HSTART-PHINEW(I))>=ERROR1) THEN
	  WRITE(*,*) 'PROGRAMME TERMINATED'
	  WRITE(*,*) 'Increase width of area OR use INFINITE elements'
	  WRITE(10,*) 'PROGRAMME TERMINATED'
	  WRITE(10,*) 'Increase width of area OR use INFINITE elements'
	  GOTO 6000
	 ENDIF
	ENDIF
 430	CONTINUE
	ENDIF
		
C	15.	����������� ��� ������� ����� ������� ��������������� �������
	DO 500 I=1,NNP
	BKQ1(I)=0.
	DO 500 J=1,NNP
 500	BKQ1(I)=BKQ1(I)+AKQ(I,J)*PHINEW(J)
	QNFIX=(BKQ1-BKQ)/DT
	QSUM=0.
	DO 501 I=1,NNP
	DO 501 J=1,NFIX
	IF (NNX(J)==I) THEN
	 QSUM=QSUM+QNFIX(I)
	ENDIF
 501	CONTINUE

C	16.	������������� ����� ����� ���� ��� ����� ��� ������� �����
C		������� ��������������� ������� ��� ��������������� �������
	VOLUMEQ=VOLQ+QSUM*DT+VOL2


C	17.	����������� ��� ��������� ��� ����� ���� ������� ����		
C		���� ��� �������������� ����� ��� �������
	CALL VOLUME(PHINEW,VOLUMEPHI,AREA,NELEM,NNP,
     $NOD,X,Y,WOA,NGP,RJ,RW,NGPINF,RJINF,RWINF,
     $INFMETHOD,RA,XSTART,YSTART,Xdirec,Ydirec,
     $YLAY,TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1,
     $TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2)


C	18.	�������� ��� ������������� ��� ��������� �������� �������
	TIME=TIME+DT
	WRITE(*,*)TIME,DT,ITER1,DMAX
	IF (MOD(TIME,1.)==0.) THEN	!�������� ��� ��� ���
	WRITE(10,2110)TIME,DT,ITER1,DMAX,IMAX,DMAX1,IMAX1
	WRITE(12,2110)TIME,DT,ITER1,DMAX,IMAX,DMAX1,IMAX1
	WRITE(13,2110)TIME,DT,ITER1,DMAX,IMAX,DMAX1,IMAX1
	WRITE(14,2110)TIME,DT,ITER1,DMAX,IMAX,DMAX1,IMAX1
	WRITE(10,*)'      X(I)       Y(I)     PHI(I)      PG(I)'
	WRITE(12,*)'      X(I)       Y(I)     PHI(I)'
	WRITE(13,*)'      X(I)       Y(I)     PHI(I)'
	DO 900 I=1,NNP
	WRITE(10,2121)X(I),Y(I),PHINEW(I),PHINEW(I)-Y(I)
	IF (X(I)==WOA) WRITE(12,2120)X(I),Y(I),PHINEW(I)
	DO 900 J=1,NFIX
	IF (NNX(J)==I) WRITE(13,2130)X(I),Y(I),PHINEW(I),QNFIX(I)
 900	CONTINUE
	WRITE(14,2140)(I,PHINEW(I),I=1,NNP)
	WRITE(13,2150)QSUM,TIME
	ENDIF						!�������� ��� ��� ���
	WRITE(11,*)TIME,VOLUMEQ,VOLUMEPHI
	WRITE(20,*)TIME,QSUM

 2110 FORMAT(///16X,'TIME(hours)=',E10.3/
     $16X,'DT(hours)  =',E10.3/
     $16X,'ITER       =',I10/
     $8X,'������� �������� �������(m)=',E10.3,1X,'���� �����=',I5/
     $8X,'������� ������ ���������(m)=',E10.3,1X,'���� �����=',I5/
     $73('-')/)
 2120 FORMAT(2(1X,F10.3),(1X,F10.3))
 2121 FORMAT(2(1X,F10.3),2(1X,F10.3))
 2130 FORMAT(2(1X,F10.3),2(1X,F10.5))
 2140 FORMAT(5(1X,I5,':',F6.3,2X))
 2150 FORMAT('������='1X,E8.3,2X,'��� �����= ',F6.3)

C	18.1	������������ ��� ��� ������� ������� ������
	DTOLD=DT
	PHIOLD=PHI
	PHI=PHINEW
	VOLQ=VOLUMEQ
	VOLPHI=VOLUMEPHI

C	19.	������� ��� ������� ����� ��� ������
C		��� ������ ��� �� ��������� - ���� 22 (GOTO)
	IF (TIME>=TSTOP) GOTO 6000

C	20.	�������� ��� �������� ������� ��� ��� ������� ������� ������
	DT1=DTOLD*PHIVAR/DMAX	!���� ��� ��������� �������
	DT2=DTOLD*1.2			!������ ���� 20%
	DT3=AINT(TIME)+1-TIME	!��������������� ��� ���
	DT=MIN(DT1,DT2,DT3)		!������� �������� ������� ��� �� ��������

C	21.	��������� ��� ���� 1 (GOTO)
	GOTO 5000

C	22.	������ ��� �� ���������
 6000	CONTINUE


							!������������� ������ ��� �����������
	CALL CPU_TIME(TIMEEND)
	TIMEEQ=TIMEEND-TIMEBEGIN
	TIMETOTAL=TIMERENUM+TIMEEQ
	WRITE(9,2200)TIMERENUM,TIMEEQ/60,TIMETOTAL/60
 2200 FORMAT(//'������������� ������ ��� �����������'/
     $'��� ��� ���������� ��� ������� (seconds) ='F10.3/
     $'��� ��� ������� ��� ��������   (minutes) ='F10.3/
     $'��������� ������������� ������ (minutes) ='F10.3)
	WRITE(*,*)
	WRITE(*,*)'RESULTS ARE IN FILES UnRiv1.res & UnRiv2.res'
	WRITE(*,*)
	WRITE(*,*)'EXECUTABLE TIME WAS',TIMEEQ,' SECONDS'
	WRITE(*,*)
	END


	SUBROUTINE VOLUME(PHI,VOLTOTAL,AREA,NELEM,NNP,
     $NOD,X,Y,WOA,NGP,RJ,RW,NGPINF,RJINF,RWINF,
     $INFMETHOD,RA,XSTART,YSTART,Xdirec,Ydirec,
     $YLAY,TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1,
     $TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2)
C	������������ � ������������� ����� �����
C	���� ��� �������������� ����� ��� �������
	DOUBLE PRECISION PHI,X,Y,SF,SFM,B,BM
	DIMENSION PHI(NNP),X(NNP),Y(NNP)
	DIMENSION SF(4),B(2,4),SFM(4),BM(2,4)
	DIMENSION NOD(NELEM,4),AREA(NELEM),WETELEM(NELEM)
	DIMENSION RW(10),RJ(10),RWINF(10),RJINF(10),RW1(10),RJ1(10)
	DIMENSION XY(4,2),PG(4),CKXY(2,2)

	VOLTOTAL=0.
	DO 5 NE=1,NELEM
	INFCHECK=0
	DO 6 I=1,4
	XY(I,1)=X(NOD(NE,I))
 	XY(I,2)=Y(NOD(NE,I))
	PG(I)=PHI(NOD(NE,I))-Y(NOD(NE,I))
	IF (XY(I,1)==WOA) INFCHECK=1
	IF (INFMETHOD==0) INFCHECK=0
  6	CONTINUE
	
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

	
	DO 10 IN=1,NGP1
	DO 10 JN=1,NGP
	CALL SHAPEFUN(IN,JN,RJ1,RJ,XY,DJ,SF,B,
     $INFCHECK,SFM,BM,INFMETHOD,RA,XSTART,YSTART,Xdirec,Ydirec)
	CALL BROOKS_COREY(YLAY,SF,SFM,XY,PG,CKXY,C,WET,
     $TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1,
     $TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2)

	IF (INFCHECK==0) WETELEM(NE)=WETELEM(NE)+RW1(IN)*RW(JN)*WET/4.
	IF (INFCHECK==1) WETELEM(NE)=WETELEM(NE)+0.5*RW1(IN)*RW(JN)*WET/4.
  10	CONTINUE
	VOLTOTAL=VOLTOTAL+WETELEM(NE)*AREA(NE)
  5	CONTINUE
	RETURN
	END



	SUBROUTINE SHAPEFUN(IN,JN,RJ,RJ2,XY,DJ,SF,B,
     $INFCHECK,SFM,BM,INFMETHOD,RA,XSTART,YSTART,XX,YY)
C	������������� ��� ���� ������ �����������
C	�� ����������� ������ (SF)
C	�� ��������� ����� (BS)
C	�� ��������� ������ (AJ) ��� � �������� ����� (DJ)
C	�� ������ �
	DOUBLE PRECISION SF,SFM,BS,BSM,B,BM
	DIMENSION RJ(10),RJ2(10),SF(4),BS(2,4),B(2,4),XY(4,2),AJ(2,2)
	DIMENSION SFM(4),BSM(2,4),BM(2,4),SR(4)
	DIMENSION SF1(4),BS1(2,4),PHI(4),PG(4)

	IF (INFCHECK==0) THEN	!�������� ��������
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
	ENDIF					!�������� ��������



	IF (INFCHECK==1) THEN	!��������� ��������
	RIN1=1.-RJ(IN)
	RIN2=1.+RJ(IN)
	RJN1=1.-RJ2(JN)
	RJN2=1.+RJ2(JN)

	IF (RJ(IN)<=0) THEN		!����������� ����� ��� ���������� ���������
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
	ENDIF					!����������� ����� ��� ���������� ���������

	IF (RJ(IN)>0) THEN		!��������� ����� ��� ���������� ���������
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
	DO 5 I=1,4
	SX=SX+SF(I)*XX*(XY(I,1)-XSTART)
	SY=SY+SF(I)*YY*(XY(I,2)-YSTART)
   5	SR(I)=SQRT((XY(I,1)-XSTART)**2+(XY(I,2)-YSTART)**2)
	SRS=SQRT(SX**2+SY**2)
	SXX=SX*AJ(1,1)+SY*AJ(1,2)
	SYY=SX*AJ(2,1)+SY*AJ(2,2)

	DO 10 I=1,4
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
  10	CONTINUE
	ENDIF					!��������� ����� ��� ���������� ���������
	ENDIF					!��������� ��������

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
C	������������� ��� ���� ������ �����������
C	������� �� ��� ������� ������ ��� ��������� ���
C	������� �� ��� ��������� ��������� ��� (���������-��������)
C	�� ������ ��� ���������� ������������ (CKXY)
C	� ��������� ������������ (C)
	DOUBLE PRECISION SF,SFM
	DIMENSION SF(4),SFM(4),XY(4,2),PG(4),CKXY(2,2)
	YNGP=0.0
  	PGNGP=0.0
	DO 10 I=1,4
	YNGP1=SF(I)*XY(I,2)
	PGNGP1=SFM(I)*PG(I)
	YNGP=YNGP+YNGP1
  	PGNGP=PGNGP+PGNGP1
  10	CONTINUE

	IF(YNGP>YLAY)THEN
	 IF(PGNGP>-HB1)THEN
	  WET=TSAT1
	  CKXX=CKSXX1
	  CKYY=CKSYY1
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
C	��������� ������� �� ��� ��������� ��� GAUSS
C	�� ������� ��� ��������� ����������� ����� �� ����� ��� �������
C	��� ������������� �������������� ����� ��� ������� (PHINEW)
	DOUBLE PRECISION A,B,PHINEW
	DIMENSION A(N,M),C(M),B(N),PHINEW(N)

	DO 10 K=1,N-1
	B(K)=B(K)/A(K,1)
	DO 20 J=2,M
	C(J)=A(K,J)
 20	A(K,J)=A(K,J)/A(K,1)
	DO 10 L=2,M
	I=K+L-1
	IF(N>=I)THEN
	 J=0
	 DO 40 LL=L,M
	 J=J+1
 40	 A(I,J)=A(I,J)-C(L)*A(K,LL)
	 B(I)=B(I)-C(L)*B(K)
	ENDIF
 10	CONTINUE

	B(N)=B(N)/A(N,1)
	DO 50 K=N-1,1,-1
	DO 50 J=2,M
	L=K+J-1
	IF(N>=L) B(K)=B(K)-A(K,J)*B(L)
  50	CONTINUE

	DO 60 I=1,N
  60	PHINEW(I)=B(I)
	RETURN
	END