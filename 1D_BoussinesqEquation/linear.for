C	NUMERICAL SCHEME TO SOLVE THE NON-DIMENSIONAL
C	BOUSSINESQ EQUATION USING
C	THE FINITE ELEMENT METHOD
C	
C	OPTIONS TO APROXIMATE THE FLOW AT THE EDGE OF THE FLOW AREA
C	1. FINITE ELEMENT
C	2. INFINITE ELEMENT, f(ri/r)=(ri/r)**a
C	3. INFINITE ELEMENT, f(ri/r)=exp{a*(ri-r)}
C
C
C	INPUT DATA IN FILE 'linear.in'
C
C
	ALLOCATABLE::H(:),A(:),B(:),C(:),D(:)
	ALLOCATABLE::PC(:),PC1(:),PC2(:),PC3(:)
	ALLOCATABLE::PV(:),PD(:),PZ1(:),PZ2(:),PZ3(:)
	ALLOCATABLE::PV1(:),PV2(:),PV3(:)
	ALLOCATABLE::Z1(:),Z2(:),Z3(:),U(:),V(:)
	DIMENSION RW(10),RJ(10)
	AX(X)=-(N-1)*DX*X/(1-X)+N*DX/(1-X)
	AM2(X)=((N*DX/AX(X))**a4)**2/(1-X)**2
	ADM2(X)=(N*DX**a4/AX(X)**(a4+1))**2/(1-X)**2
	BM2(X)=(EXP(a4*(N*DX-AX(X))))**2/(1-X)**2
	BDM2(X)=(EXP(a4*(N*DX-AX(X))))**2/(1-X)**2

	BX1(X)=(N*PDX*(1-X)/(N*PDX-(N-1)*PDX*X))**a5
	BX2(X)=-a5*(N*PDX)**a5*(1-X)**(a5+1)*PDX/
     $(N*PDX-(N-1)*PDX*X)**(a5+1)/(1-X)**2
	CX1(X)=EXP(-a5*X*PDX/(1-X))
	CX2(X)=-a5*EXP(-a5*X*PDX/(1-X))*PDX/(1-X)**2
	AOLD(X)=BX1(X)*BX2(X)*(1-X)**2/PDX
	AOLV(X)=BX1(X)**2*BX2(X)
	AOLC2(X)=(BX1(X)/(1-X))**2
	BOLD(X)=CX1(X)*CX2(X)*(1-X)**2/PDX
	BOLV(X)=CX1(X)**2*CX2(X)
	BOLC2(X)=(CX1(X)/(1-X))**2

	OPEN(8,FILE='linear.in')
	OPEN(9,FILE='lineara.out')
	OPEN(10,FILE='linearb.out')
	OPEN(11,FILE='linearc.out')
	OPEN(12,FILE='lineard.out')

	OPEN(15,FILE='Q1.out')
	OPEN(16,FILE='Q2.out')
	OPEN(17,FILE='Q3.out')
	OPEN(19,FILE='IX1.out')
	OPEN(20,FILE='IX2.out')
	OPEN(21,FILE='IX3.out')
	OPEN(23,FILE='Z1.out')
	OPEN(24,FILE='Z2.out')
	OPEN(25,FILE='Z3.out')
	OPEN(30,FILE='POLLUTION.out')
	OPEN(31,FILE='POLLUTIONa.out')
	OPEN(32,FILE='POLLUTION1.out')
	OPEN(33,FILE='POLLUTION2.out')
	OPEN(34,FILE='POLLUTION3.out')
	OPEN(35,FILE='VELOCITY1.out')
	OPEN(36,FILE='VELOCITY2.out')
	OPEN(37,FILE='VELOCITY3.out')



	WRITE(*,300)
	READ(8,100)DX,DT,N,NB,NNJ,NT,NN,IX,ICOMPUTa4,
     $a4,a5,RK,RS,RNE,HO,HA,RD,PCO,PCA,AL,PN,PDD
	WRITE(*,*)PCO,PCA,AL,PN,PDD
	READ(8,*)
	READ(8,*)
	READ(8,*)KMAX
	READ(8,*)
	READ(8,*)(RW(K),K=1,KMAX)
	READ(8,*)
	READ(8,*)(RJ(K),K=1,KMAX)
	NK=N;NBK=NB-1
	RA=RK*RD/RS
	PDX=DX*ABS(HO-HA)
	PDT=DT*(HO-HA)**2/RA
	ALLOCATE(H(NBK),A(NBK),B(NBK),C(NBK),D(NBK))
	ALLOCATE(PC(NBK),PC1(NBK),PC2(NBK),PC3(NBK))
	ALLOCATE(PV(0:NBK),PD(0:NBK),PZ1(0:NBK),PZ2(0:NBK),PZ3(0:NBK))
	ALLOCATE(PV1(0:NBK),PV2(0:NBK),PV3(0:NBK))
	ALLOCATE(Z1(NBK),Z2(NBK),Z3(NBK),U(NBK),V(NBK))
	PZ1(0)=HO;PZ2(0)=HO;PZ3(0)=HO

	WRITE(9,200)DX,DT,N,NB,NT,KMAX,a4
	WRITE(10,201)
	WRITE(11,202)DX*ABS(HO-HA),DT*(HO-HA)**2/RA,N,NB,NT,KMAX,a4,
     $RK,RS,HO,HA,RD
	WRITE(30,204)DX*ABS(HO-HA),DT*(HO-HA)**2/RA,N,NB,NT,KMAX,
     $a4,RK,RS,HO,HA,RD
	WRITE(12,203)IX,DX*IX
 300	FORMAT('Numerical experiments in linear Boussinesqs eguation'/
     $'with three differnte methods'/'1. FINITE ELEMENTS'/
     $'2. USING INFINITE ELEMENT FOR f(ri/r)=(ri/r)**a'/
     $'3. USING INFINITE ELEMENT FOR f(ri/r)=exp{a*(ri-r)}'//
     $'Input data in file linear.in'/
     $'WAIT...'/)
 100	FORMAT(/////////////F10.5///F10.5,7(///I10),13(///F10.5))
 200	FORMAT('аявеио енодоу апотекеслатым - lineara.out'//
     $'епикусг тгс цяаллийопоиглемгс(ыс пяос г**2)'/
     $'йаи адиастатгс енисысгс тоу BOUSSINESQ'/
     $'йахыс йаи тгс енисысгс тгс летажояас лафас'/
     $'ле тг леходо тым пепеяаслемым стоивеиым'//
     $'тяопои пяоссеццисгс сто текос тгс пеяиовгс лекетгс'/
     $'1. пепеяаслемо стоивеио - Z1(J)'/
     $'2. глиапеияо стоивеио опоу f(ri/r)=(ri/r)**a - Z2(J)'/
     $'3. глиапеияо стоивеио опоу f(ri/r)=exp{a*(ri-r)} - Z3(J)'//
     $'та апотекеслата еимаи се адиастатг лояжг'//
     $'DX=',F10.5/'DT=',F10.5/'N=',I10/'NB=',I10/
     $'NT=',I10/'KMAX=',I10/'a=',F10.5/
     $'I=вяомийос деийтгс летабокгс'/'J=выяийос деийтгс летабокгс'/
     $'TIME=адиастатос вяомос'/////'апотекеслата')
 201	FORMAT('аявеио енодоу апотекеслатым - linearb.out'//
     $'епикусг тгс цяаллийопоиглемгс(ыс пяос г**2)'/
     $'йаи адиастатгс енисысгс тоу BOUSSINESQ'//
     $'тилес тгс паяовгс се адиастатг лояжг'//
     $4X,'I',5X,'Q1(I)',5X,'Q2(I)',5X,'Q3(I)')
 202	FORMAT('аявеио енодоу апотекеслатым - linf1rc.out'//
     $'епикусг тгс цяаллийопоиглемгс(ыс пяос г**2)'/
     $'йаи адиастатгс енисысгс тоу BOUSSINESQ'/
     $'йахыс йаи тгс енисысгс тгс летажояас лафас'/
     $'ле тг леходо тым пепеяаслемым стоивеиым'//
     $'тяопои пяоссеццисгс сто текос тгс пеяиовгс лекетгс'/
     $'1. пепеяаслемо стоивеио - Z1(J)'/
     $'2. глиапеияо стоивеио опоу f(ri/r)=(ri/r)**a - Z2(J)'/
     $'3. глиапеияо стоивеио опоу f(ri/r)=exp{a*(ri-r)} - Z3(J)'//
     $'та апотекеслата еимаи се ломадес тоу летяийоу сустглатос'//
     $'DX(m)=',F10.5/'DT(day)=',F10.5/'N=',I10/'NB=',I10/'NT=',I10/
     $'KMAX=',I10/'a=',F10.5/'RK(m/day)=',F10.5/'RS=',F10.5/
     $'HO(m)=',F10.5/'HA(m)=',F10.5/'RD(m)=',F10.5/
     $'I=вяомийос деийтгс летабокгс'/'J=выяийос деийтгс летабокгс'/
     $'TIME=адиастатос вяомос'/////'апотекеслата')
 203	FORMAT('аявеио енодоу апотекеслатым - lineard.out'//
     $'епикусг тгс цяаллийопоиглемгс(ыс пяос г**2)'/
     $'йаи адиастатгс енисысгс тоу BOUSSINESQ'//
     $'тилес тгс упоцеиас стахлгс стом йолбо ив'//
     $'йолбос ив=',I5,2X,'хесг тоу йолбоу =',E8.2//
     $4X,'I',4X,'ф1(Iв)',4X,'ф2(Iв)',4X,'ф3(Iв)',5X,'a4-Z2',6X,'a4-Z3')	
 204	FORMAT('аявеио енодоу апотекеслатым - Polution.out'//
     $'епикусг тгс цяаллийопоиглемгс(ыс пяос г**2)'/
     $'йаи адиастатгс енисысгс тоу BOUSSINESQ'/
     $'йахыс йаи тгс енисысгс тгс летажояас лафас'/
     $'ле тг леходо тым пепеяаслемым стоивеиым'//
     $'та апотекеслата еимаи се ломадес тоу летяийоу сустглатос'/
     $'циа тис тилес тоу жоятиоу се йкеисто удяожояеа-PZ1(m)'/
     $'йаи циа тис тилес тгс суцйемтяысгсг-PC1(mg/l)'//
     $'DX(m)=',F10.5/'DT(day)=',F10.5/'N=',I10/
     $'NB=',I10/'NT=',I10/'KMAX=',I10/'a=',F10.5/'RK(m/day)=',F10.5/
     $'RS=',F10.5/'HO(m)=',F10.5/'HA(m)=',F10.5/'RD(m)=',F10.5/
     $'I=вяомийос деийтгс летабокгс'/'J=выяийос деийтгс летабокгс'/
     $'TIME=адиастатос вяомос'/////'апотекеслата')


C	еИСАЦЫЦч ОЯИАЙчР СУМХчЙГР СТГМ ТэЖЯО & АЯВИЙчР СУМХчЙГР
C	*******************************************************
	IF(HO>HA)THEN
	HOD=1
	ELSE
	HOD=0
	ENDIF
	DO 5 J=1,NBK
	Z1(J)=1-HOD
	Z2(J)=1-HOD
	Z3(J)=1-HOD
	PC1(J)=PCA
	PC2(J)=PCA
	PC3(J)=PCA
 5    CONTINUE

	PZZ1=HA*HOD+HO*(1-HOD)
	PZZ2=ABS(HO-HA)
	write(*,*)PZZ1,PZZ2

	DO 10 I=1,NT
	TIME=DT*I

	DO 15 METHOD=1,3
	IF(METHOD==1) H=Z1
	IF(METHOD==2) H=Z2
	IF(METHOD==3) H=Z3

	
	IF(METHOD==2 .AND. HO>HA .AND. ICOMPUTa4==1) THEN
	IF (ABS(H(NK))>=0.0001) THEN
	a4=ABS(LOG((1-H(NK))/(1+H(NK))))/LOG(((N-1)*DX))
	ELSE
	a4=0.
	ENDIF
	ENDIF

	IF(METHOD==3 .AND. HO>HA .AND. ICOMPUTa4==1) THEN
	IF (ABS(H(NK))>=0.0001) THEN
	a4=ABS(LOG((1-H(NK))/(1+H(NK))))/((N-1)*DX)
	ELSE
	a4=0.
	ENDIF
	ENDIF	
	

	OLM1=0;OLM2=0;DOLM1=0;DOLM2=0
	IF (METHOD==2)THEN
	DO 2 K=1,KMAX
	DOLM2=DOLM2+RW(K)*ADM2(RJ(K))/2
	OLM2=OLM2+RW(K)*AM2(RJ(K))/2
 2	CONTINUE
	ENDIF
	IF (METHOD==3)THEN
	DO 3 K=1,KMAX
	DOLM2=DOLM2+RW(K)*BDM2(RJ(K))/2
	OLM2=OLM2+RW(K)*BM2(RJ(K))/2
 3	CONTINUE
	ENDIF



C	сУМТЕКЕСТщР ПЯЧТГР ЕНъСЫСГР
C	***************************
	B(1)=2/DX+2*DX/3/DT
	C(1)=1/DX-DX/6/DT
	D(1)=(1/DX)*HOD+(2*DX/3/DT)*H(1)+(DX/6/DT)*H(2)

C	сУМТЕКЕСТщР ТЫМ УПЭКОИПЫМ ЕНИСЧСЕЫМ
C	***********************************
	DO 20 J=2,NBK-1
	A(J)=1/DX-DX/6/DT
	B(J)=2/DX+2*DX/3/DT
	C(J)=1/DX-DX/6/DT
	D(J)=(DX/6/DT)*H(J-1)+(2*DX/3/DT)*H(J)+(DX/6/DT)*H(J+1)
 20   CONTINUE
	
C	сУМТЕКЕСТщР ТЕКЕУТАъАР ЕНъСЫСГР
C	*******************************
	IF(METHOD==1) THEN
	A(NBK)=1/DX-DX/6/DT
	B(NBK)=2/DX+2*DX/3/DT
	D(NBK)=(DX/6/DT)*H(NBK-1)+(2*DX/3/DT+1/DX)*H(NBK)
	ELSE
	A(NK-1)=1/DX-DX/6/DT
	B(NK-1)=2/DX+2*DX/3/DT
	C(NK-1)=1/DX-DX/6/DT
	D(NK-1)=(DX/6/DT)*H(NK-2)+(2*DX/3/DT)*H(NK-1)+
     $(DX/6/DT)*H(NK)

	A(NK)=1/DX-DX/6/DT
	B(NK)=1/DX+DX/3/DT+(A4**2*DX*DOLM2+OLM2*DX/DT)
	D(NK)=(DX/6/DT)*H(NK-1)+(DX/3/DT)*H(NK)+(OLM2*DX/DT)*H(NK)
	ENDIF
	

C	аКЦЭЯИХЛОР ТОУ THOMAS
C	*********************
	U(1)=C(1)/B(1)
	V(1)=D(1)/B(1)
	IF(METHOD==1) THEN
	DO 30 J=2,NBK
	U(J)=C(J)/(B(J)-A(J)*U(J-1))
	V(J)=(D(J)+A(J)*V(J-1))/(B(J)-A(J)*U(J-1))
 30   CONTINUE
	DO 40 J=NBK-1,1,-1
	V(J)=V(J)+U(J)*V(J+1)
 40   CONTINUE
	ELSE
	DO 31 J=2,NK
	U(J)=C(J)/(B(J)-A(J)*U(J-1))
	V(J)=(D(J)+A(J)*V(J-1))/(B(J)-A(J)*U(J-1))
 31   CONTINUE
	DO 41 J=NK-1,1,-1
	V(J)=V(J)+U(J)*V(J+1)
 41   CONTINUE
	ENDIF

C	пЯОЕТОИЛАСъА ЦИА ТГМ ЕПЭЛЕМГ ВЯОМИЙч СТИЦЛч
C	*******************************************
	IF(METHOD==1) Z1=V
	IF(METHOD==2) Z2=V
	IF(METHOD==3) Z3=V
 15	CONTINUE



C	лЕТАЖОЯэ ЛэФАР
C	**************
	DO METHODC=1,3
	IF(METHODC==1) THEN
	DO J=1,NBK
	PC(J)=PC1(J)
	PZ1(J)=PZZ1+PZZ2*Z1(J)
	PV(J-1)=RK*ABS(PZ1(J-1)-PZ1(J))/PDX/RNE
	PV1(J-1)=PV(J-1)
	PD(J-1)=(PDD+AL*(PV(J-1)*100)**PN)/10000
	ENDDO
	PV(NBK)=PV(NBK-1)
	PV1(NBK)=PV(NBK-1)
	PD(NBK)=PD(NBK-1)
	ENDIF
	IF(METHODC==2) THEN
	DO J=1,NK
	PC(J)=PC2(J)
	PZ2(J)=PZZ1+PZZ2*Z2(J)
	PV(J-1)=RK*ABS(PZ2(J-1)-PZ2(J))/PDX/RNE
	PV2(J-1)=PV(J-1)
	PD(J-1)=(PDD+AL*(PV(J-1)*100)**PN)/10000
	ENDDO
	PV(NK)=PV(NK-1)
	PV2(NK)=PV(NK-1)
	PD(NK)=PD(NK-1)
	ENDIF
	IF(METHODC==3) THEN
	DO J=1,NK
	PC(J)=PC3(J)
	PZ3(J)=PZZ1+PZZ2*Z3(J)
	PV(J-1)=RK*ABS(PZ3(J-1)-PZ3(J))/PDX/RNE
	PV3(J-1)=PV(J-1)
	PD(J-1)=(PDD+AL*(PV(J-1)*100)**PN)/10000
	ENDDO
	PV(NK)=PV(NK-1)
	PV3(NK)=PV(NK-1)
	PD(NK)=PD(NK-1)
	ENDIF

	OLD=0;OLV=0;OLC2=0
	IF (METHODC==2)THEN
	DO K=1,KMAX
	OLD=OLD+RW(K)*AOLD(RJ(K))/2
	OLV=OLV+RW(K)*AOLV(RJ(K))/2
	OLC2=OLC2+RW(K)*AOLC2(RJ(K))/2
	ENDDO
	ENDIF

	IF (METHODC==3)THEN
	DO K=1,KMAX
	OLD=OLD+RW(K)*BOLD(RJ(K))/2
	OLV=OLV+RW(K)*BOLV(RJ(K))/2
	OLC2=OLC2+RW(K)*BOLC2(RJ(K))/2
	ENDDO
	ENDIF



C	сУМТЕКЕСТщР ПЯЧТГР ЕНъСЫСГР
C	************************
	B(1)=(PD(0)+2*PD(1)+PD(2))/2/PDX-(PV(0)-PV(2))/6+2/PDX+2*PDX/3/PDT
	C(1)=(PD(1)+PD(2))/2/PDX-(PV(1)+2*PV(2))/6+1/PDX-PDX/6/PDT
	D(1)=((PD(0)+PD(1))/2/PDX+(2*PV(0)+PV(1))/6+1/PDX)*PCO
     $+(2*PDX/3/PDT)*PC(1)+(PDX/6/PDT)*PC(2)

C	сУМТЕКЕСТщР ТЫМ УПЭКОИПЫМ ЕНИСЧСЕЫМ
C	***********************************
	DO  J=2,NBK-1
	A(J)=(PD(J-1)+PD(J))/2/PDX+(2*PV(J-1)+PV(J))/6+1/PDX-PDX/6/PDT
	B(J)=(PD(J-1)+2*PD(J)+PD(J+1))/2/PDX-(PV(J-1)-PV(J+1))/6
     $+2/PDX+2*PDX/3/PDT
	C(J)=(PD(J)+PD(J+1))/2/PDX-(PV(J)+2*PV(J+1))/6+1/PDX-PDX/6/PDT
	D(J)=(PDX/6/PDT)*PC(J-1)+(2*PDX/3/PDT)*PC(J)+(PDX/6/PDT)*PC(J+1)
	ENDDO
	
C	сУМТЕКЕСТщР ТЕКЕУТАъАР ЕНъСЫСГР
C	*******************************
	IF (METHODC==1) THEN
	A(NBK)=(PD(NBK-1)+PD(NBK))/2/PDX+(2*PV(NBK-1)+PV(NBK))/6
     $+1/PDX-PDX/6/PDT
	B(NBK)=(PD(NBK-1)+PD(NBK))/2/PDX-(PV(NBK-1)-4*PV(NBK))/6
     $+1/PDX+2*PDX/3/PDT
	D(NBK)=(PDX/6/PDT)*PC(NBK-1)+(2*PDX/3/PDT)*PC(NBK)
	ELSE
	A(NK-1)=(PD(NK-2)+PD(NK-1))/2/PDX+(2*PV(NK-2)
     $+PV(NK-1))/6+1/PDX-PDX/6/PDT
	B(NK-1)=(PD(NK-2)+2*PD(NK-1)+PD(NK))/2/PDX-(PV(NK-2)
     $-PV(NK))/6+2/PDX+2*PDX/3/PDT
	C(NK-1)=(PD(NK-1)+PD(NK))/2/PDX-(PV(NK-1)+2*PV(NK))/6
     $+1/PDX-PDX/6/PDT
	D(NK-1)=(PDX/6/PDT)*PC(NK-2)+(2*PDX/3/PDT)*PC(NK-1)
     $+(PDX/6/PDT)*PC(NK)

	A(NK)=(PD(NK-1)+PD(NK))/2/PDX+(2*PV(NK-1)+PV(NK))/6
     $+1/PDX-PDX/6/PDT
	B(NK)=(PD(NK-1)+PD(NK))/2/PDX-(PV(NK-1)+2*PV(NK))/6
     $+1/PDX+PDX/3/PDT+PD(NK)*OLD-PV(NK)*OLV+PDX*OLC2/PDT
	D(NK)=(PDX/6/PDT)*PC(NK-1)+(PDX/3/PDT+PDX*OLC2/PDT)*PC(NK)
	ENDIF

C	аКЦЭЯИХЛОР ТОУ THOMAS
C	*********************
	U(1)=C(1)/B(1)
	V(1)=D(1)/B(1)
	IF(METHODC==1) THEN
	DO J=2,NBK
	U(J)=C(J)/(B(J)-A(J)*U(J-1))
	V(J)=(D(J)+A(J)*V(J-1))/(B(J)-A(J)*U(J-1))
	ENDDO
	DO J=NBK-1,1,-1
	V(J)=V(J)+U(J)*V(J+1)
	ENDDO
	ELSE
	DO J=2,NK
	U(J)=C(J)/(B(J)-A(J)*U(J-1))
	V(J)=(D(J)+A(J)*V(J-1))/(B(J)-A(J)*U(J-1))
	ENDDO
	DO J=NK-1,1,-1
	V(J)=V(J)+U(J)*V(J+1)
	ENDDO
	ENDIF

C	пЯОЕТОИЛАСъА ЦИА ТГМ ЕПЭЛЕМГ ВЯОМИЙч СТИЦЛч
C	*******************************************
	IF(METHODC==1) PC1=V
	IF(METHODC==2) PC2=V
	IF(METHODC==3) PC3=V
	ENDDO


C	еЙТЩПЫСГ ТЫМ АПОТЕКЕСЛэТЫМ
C	*************************
	IF(MOD(I,NN)==0)THEN
c	IF(time==1. .or. time==1.5 .or. time==2. .or. time==3.
c     $ .or. time==5. .or. time==8.
c     $ .or. time==12. .or. time==20. .or. time==30.
c     $ .or. time==50. .or. time==100. .or. time==200.
c     $ .or. time==400. .or. time==1000.)THEN
C	еЙТЩПЫСГ ТЫМ ЖОЯТъЫМ СЕ АДИэСТАТЕР ЙАИ ПЯАЦЛАТИЙщР ТИЛщР
	WRITE(9,210)I,TIME,DX
	WRITE(23,*)I,TIME
	WRITE(24,*)I,TIME
	WRITE(25,*)I,TIME
	WRITE(11,210)I,TIME*(HO-HA)**2/RA,DX*ABS(HO-HA)
	WRITE(30,210)I,TIME*(HO-HA)**2/RA,DX*ABS(HO-HA)
210   FORMAT(//'I=',I10,3X,'TIME=',F10.2,3X,'DX=',F10.2/46('*'))
	WRITE(9,*)'   J     Z1(J)     Z2(J)     Z3(J)'
	WRITE(11,*)'   J     Z1(J)     Z2(J)     Z3(J)'
	WRITE(30,*)'          J     PZ1(J)         PC1(J)         PV1(J)'
	DO 50 J=NNJ,NK,NNJ
	WRITE(23,*)Z1(J)
	WRITE(24,*)Z2(J)
	WRITE(25,*)Z3(J)
	WRITE(9,220)J,Z1(J),Z2(J),Z3(J)
	WRITE(11,220)J,HA**2*HOD+HO**2*(1-HOD)+ABS(HO**2-HA**2)*Z1(J),
     $HA**2*HOD+HO**2*(1-HOD)+ABS(HO**2-HA**2)*Z2(J),
     $HA**2*HOD+HO**2*(1-HOD)+ABS(HO**2-HA**2)*Z3(J)
 50   CONTINUE
	WRITE(9,*)' NBK   Z1(NBK)   NBK*DX'
	WRITE(9,230)NBK,Z1(NBK),NBK*DX
	WRITE(11,*)' NBK   Z1(NBK)   NBK*DX'
	WRITE(11,230)NBK,HA**2*HOD+HO**2*(1-HOD)+ABS(HO**2-HA**2)*Z1(NBK),
     $NBK*DX*ABS(HO-HA)

C	уПОКОЦИСЛЭР ЙАИ ЕЙТЩПЫСГ ТГР ПАЯОВчР СЕ АДИэСТАТЕР ТИЛщР
	Q1=ABS(HOD-Z1(1))/DX
	Q2=ABS(HOD-Z2(1))/DX
	Q3=ABS(HOD-Z3(1))/DX
	WRITE(10,220)I,Q1,Q2,Q3
	WRITE(15,*)Q1
	WRITE(16,*)Q2
	WRITE(17,*)Q3


C	еЙТЩПЫСГ ТГР УПЭЦЕИАР СТэХЛГР СТОМ ЙЭЛБО IX
	WRITE(12,222)I,Z1(IX),Z2(IX),Z3(IX),XI2,XI3
	WRITE(19,*)Z1(IX)
	WRITE(20,*)Z2(IX)
	WRITE(21,*)Z3(IX)

C	еЙТЩПЫСГ ТЫМ АПОТЕКЕСЛэТЫМ ЛЕТАЖОЯэР ЛэФАР
	WRITE(31,*)I,TIME*(HO-HA)**2/RA
	DO J=NNJ,NBK,NNJ
	WRITE(30,*)J,PZ1(J),PC1(J),PV1(J)
	WRITE(31,*)PC1(J)
	ENDDO


	WRITE(32,*)I,TIME*(HO-HA)**2/RA
	WRITE(33,*)I,TIME*(HO-HA)**2/RA
	WRITE(34,*)I,TIME*(HO-HA)**2/RA
	WRITE(35,*)I,TIME*(HO-HA)**2/RA
	WRITE(36,*)I,TIME*(HO-HA)**2/RA
	WRITE(37,*)I,TIME*(HO-HA)**2/RA
	DO J=NNJ,NK,NNJ
	WRITE(32,*)PC1(J)
	WRITE(33,*)PC2(J)
	WRITE(34,*)PC3(J)
	WRITE(35,*)PV1(J)
	WRITE(36,*)PV2(J)
	WRITE(37,*)PV3(J)
	ENDDO


	ENDIF
220	FORMAT(I5,2X,F8.5,2X,F8.5,2X,F8.5)
222	FORMAT(I5,2X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5)
221	FORMAT(I5,2X,E10.3,2X,E10.3,2X,E10.3)
230	FORMAT(I5,2X,F8.5,2X,F8.3)
 10   CONTINUE
	WRITE(*,310)
310	FORMAT('Results are in files linear(a,b,c,d).out')
	END
