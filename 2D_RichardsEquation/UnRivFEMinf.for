!	ΕΠΙΛΥΣΗ ΤΗΣ ΔΙΔΙΑΣΤΑΤΗΣ ΕΞΙΣΩΣΗΣ ΤΟΥ RICHARDS
!	ΜΕ ΤΗ ΜΕΘΟΔΟ ΤΩΝ ΠΕΠΕΡΑΣΜΕΝΩΝ ΣΤΟΙΧΕΙΩΝ
!	ΚΑΙ ΤΗ ΧΡΗΣΗ ΗΜΙΑΠΕΙΡΟΥ ΣΤΟΙΧΕΙΟΥ
!
!	ΜΕΤΑΒΛΗΤΕΣ ΤΟΥ ΠΡΟΓΡΑΜΜΑΤΟΣ
!NELEM		ΑΡΙΘΜΟΣ ΤΩΝ ΣΤΟΙΧΕΙΩΝ ΤΟΥ ΔΙΚΤΥΟΥ
!NNP			ΑΡΙΘΜΟΣ ΤΩΝ ΚΟΜΒΩΝ ΤΟΥ ΔΙΚΤΥΟΥ
!NBAND		ΕΥΡΟΣ(BANDWIDTH) ΤΟΥ ΔΙΚΤΥΟΥ
!
!PHIOLD(NNP)	ΦΟΡΤΙΟ ΣΤΗΝ ΠΡΟΗΓΟΥΜΕΝΗ ΧΡΟΝΙΚΗ ΣΤΙΓΜΗ
!PHI(NNP)		ΦΟΡΤΙΟ ΣΤΗΝ ΠΑΡΟΥΣΑ ΧΡΟΝΙΚΗ ΣΤΙΓΜΗ
!PHINEW(NNP)	ΥΠΟΛΟΓΙΖΟΜΕΝΟ ΦΟΡΤΙΟ ΣΤΗΝ ΕΠΟΜΕΝΗ ΧΡΟΝΙΚΗ ΣΤΙΓΜΗ
!PHIEST(NNP)	ΕΚΤΙΜΟΥΜΕΝΟ ΦΟΡΤΙΟ ΣΤΗΝ ΕΠΟΜΕΝΗ ΧΡΟΝΙΚΗ ΣΤΙΓΜΗ
!PHALF(NNP)	ΦΟΡΤΙΟ ΣΤΗΝ ΕΝΔΙΑΜΕΣΗ ΧΡΟΝΙΚΗ ΣΤΙΓΜΗ
!AREA(NELEM)	ΕΜΒΑΔΟ ΤΟΥ ΚΑΘΕ ΣΤΟΙΧΕΙΟΥ
!NOD(NELEM,NNP)	ΠΕΡΙΕΧΕΙ ΤΟΥΣ ΚΟΜΒΟΥΣ ΓΙΑ ΚΑΘΕ ΣΤΟΙΧΕΙΟ
!X(NNP),Y(NNP)	ΣΥΝΤΕΤΑΓΜΕΝΕΣ ΤΩΝ ΚΟΜΒΩΝ ΤΟΥ ΔΙΚΤΥΟΥ
!
!NFIX			ΑΡΙΘΜΟΣ ΚΟΜΒΩΝ ΠΡΟΚΑΘΟΡΙΣΜΕΝΟΥ ΦΟΡΤΙΟΥ
!NNX(NFIX)	ΚΟΜΒΟΙ ΠΡ. ΦΟΡΤΙΟΥ
!YY(NFIX)		ΤΙΜΗ ΠΡ. ΦΟΡΤΙΟΥ
!NQ			ΑΡΙΘΜΟΣ ΚΟΜΒΩΝ ΠΡΟΚΑΘΟΡΙΣΜΕΝΗΣ ΠΑΡΟΧΗΣ
!IFLOW(NQ)	ΚΟΜΒΟΙ ΠΡ. ΠΑΡΟΧΗΣ
!QL(NQ)		ΜΗΚΟΣ ΡΟΗΣ ΣΤΟΥΣ ΚΟΜΒΟΥΣ ΠΡ. ΠΑΡΟΧΗΣ
!Q(NQ)		ΥΠΟΛΟΓΙΖΟΜΕΝΗ ΠΑΡΟΧΗ ΣΤΟΥΣ ΚΟΜΒΟΥΣ ΠΡ. ΠΑΡΟΧΗΣ
!	
!NGP			ΣΗΜΕΙΑ ΑΡΙΘΜΗΤΙΚΗΣ ΟΛΟΚΛΗΡΩΣΗΣ ΠΡΟΣ ΚΑΘΕ ΚΑΤΕΥΘΥΝΣΗ
!RW(10)		ΣΤΑΘΜΙΣΤΙΚΟΙ ΠΑΡΑΓΟΝΤΕΣ ΣΗΜΕΙΩΝ ΑΡ. ΟΛΟΚΛΗΡΩΣΗΣ
!RJ(10)		ΣΥΝΤΕΤΑΓΜΕΝΕΣ ΣΗΜΕΙΩΝ ΑΡ. ΟΛΟΚΛΗΡΩΣΗΣ
!XY(4,2)		ΣΥΝΤΕΤΑΓΜΕΝΕΣ ΤΩΝ ΚΟΜΒΩΝ ΕΝΟΣ ΣΤΟΙΧΕΙΟΥ
!AJ(2,2)		ΤΟ ΙΑΚΩΒΙΑΝΟ ΜΗΤΡΩΟ
!
!			ΓΕΩΜΕΤΡΙΚΗ ΠΡΟΣΕΓΓΙΣΗ
!SF(4)		ΣΥΝΑΡΤΗΣΕΙΣ ΜΟΡΦΗΣ ΤΩΝ ΣΤΟΙΧΕΙΩΝ
!BS(2,4)		ΠΑΡΑΓΩΓΟΙ ΤΩΝ ΣΥΝΑΡΤΗΣΕΩΝ ΜΟΡΦΗΣ
!B(2,4)		ΤΟ ΜΗΤΡΩΟ Β=AJ*BS
!
!			ΠΡΟΣΕΓΓΙΣΗ ΤΗΣ ΕΞΑΡΤΗΜΕΝΗΣ ΜΕΤΑΒΛΗΤΗΣ
!SFM(4)		ΣΥΝΑΡΤΗΣΕΙΣ ΜΟΡΦΗΣ ΤΩΝ ΣΤΟΙΧΕΙΩΝ
!BSM(2,4)		ΠΑΡΑΓΩΓΟΙ ΤΩΝ ΣΥΝΑΡΤΗΣΕΩΝ ΜΟΡΦΗΣ
!BM(2,4)		ΤΟ ΜΗΤΡΩΟ ΒM=AJ*BSM
!
!WOA			ΤΟ ΠΛΑΤΟΣ ΤΗΣ ΠΕΡΙΟΧΗΣ ΜΕΛΕΤΗΣ
!YLAY			ΤΟ ΠΑΧΟΣ ΤΗΣ ΚΑΤΩ ΕΔΑΦΙΚΗΣ ΣΤΡΩΣΗΣ
!
!PG(4)		ΥΨΟΣ ΠΙΕΣΗΣ ΤΩΝ ΚΟΜΒΩΝ ΕΝΟΣ ΣΤΟΙΧΕΙΟΥ
!TSAT1,TRES1	ΥΓΡΑΣΙΑ ΚΟΡΕΣΜΟΥ ΚΑΙ ΥΠΟΛΕΙΜΑΤΙΚΗ ΥΓΡΑΣΙΑ
!CKSXX1		ΥΔΡΑΥΛΙΚΗ ΑΓΩΓΙΜΟΤΗΤΑ ΚΟΡΕΣΜΟΥ ΚΑΤΑ ΤΗΝ ΟΡΙΖΟΝΤΙΑ ΔΙΕΥΘΥΝΣΗ
!CKSYY1		ΥΔΡΑΥΛΙΚΗ ΑΓΩΓΙΜΟΤΗΤΑ ΚΟΡΕΣΜΟΥ ΚΑΤΑ ΤΗΝ ΚΑΤΑΚΟΡΥΦΗ ΔΙΕΥΘΥΝΣΗ
!HB1,CLM1,CN1	ΟΙ ΠΑΡΑΜΕΤΡΟΙ (Hb,λ και n) ΣΤΙΣ ΣΧΕΣΕΙΣ ΤΩΝ BROOKS-KOREY
!CKXY(2,2)	ΜΗΤΡΩΟ ΥΔΡΑΥΛΙΚΗΣ ΑΓΩΓΙΜΟΤΗΤΑΣ
!C			ΥΔΡΑΥΛΙΚΗ ΧΩΡΗΤΙΚΟΤΗΤΑ
!
!			ΒΟΗΘΗΤΙΚΑ ΜΗΤΡΩΑ ΓΙΑ ΤΟΥΣ ΥΠΟΛΟΓΙΣΜΟΥΣ ΣΕ ΚΑΘΕ ΣΤΟΙΧΕΙΟ
!B1(4,2),AK1(4,4),AK2(4,4),BK1(4)
!
!AK(NNP,NBAND),BK(NNP)	ΤΕΛΙΚΑ ΜΗΤΡΩΑ ΠΡΙΝ ΤΟΝ ΑΛΓΟΡΙΘΜΟ ΤΟΥ GAUSS
!
!QNFIX(NNP)		ΥΠΟΛΟΓΙΖΟΜΕΝΗ ΠΑΡΟΧΗ ΣΤΟΥΣ ΚΟΜΒΟΥΣ ΠΡΟΚΑΘΟΡΙΣΜΕΝΟΥ ΦΟΡΤΙΟΥ
!				ΒΟΗΘΗΤΙΚΑ ΜΗΤΡΩΑ ΓΙΑ ΤΟΝ ΥΠΟΛΟΓΙΣΜΟ ΤΗΣ ΠΑΡΟΧΗΣ
!AKQ(NNP,NBAND),BKQ(NNP),BKQ1(NNP)
!
!HSTART		Η ΑΡΧΙΚΗ ΤΙΜΗ ΤΟΥ ΦΟΡΤΙΟΥ ΣΕ ΟΛΟΥΣ ΤΟΥΣ ΚΟΜΒΟΥΣ
!ERROR1		ΤΟ ΚΡΙΤΗΡΙΟ ΣΥΓΚΛΙΣΗΣ ΤΟΥ ΦΟΡΤΙΟΥ
!PHIVAR		Η ΜΕΓΙΣΤΗ ΜΕΤΑΒΟΛΗ ΤΟΥ ΦΟΡΤΙΟΥ ΣΕ ΚΑΘΕ ΧΡΟΝΙΚΟ ΒΗΜΑ
!RAIN			ΠΑΡΟΧΗ ΑΡΔΕΥΣΗΣ Ή ΒΡΟΧΩΠΤΩΣΗΣ
!
!TIME			Η ΧΡΟΝΙΚΗ ΣΤΙΓΜΗ ΚΑΤΑ ΤΟΥΣ ΥΠΟΛΟΓΙΣΜΟΥΣ
!TSTOP		Ο ΤΕΛΙΚΟΣ ΧΡΟΝΟΣ
!DT			ΤΟ ΑΡΧΙΚΟ ΧΡΟΝΙΚΟ ΒΗΜΑ
!ITERMAX		Ο ΑΡΙΘΜΟΣ ΤΩΝ ΕΠΑΝΑΛΗΨΕΩΝ ΣΥΓΚΛΙΣΗΣ
!ITERCHECK	ΔΕΙΚΤΗΣ ΓΙΑ ΤΟΝ ΕΛΕΓΧΟ ΤΗΣ ΣΥΓΚΛΙΣΗΣ
!
!INFMETHOD	ΔΕΙΚΤΗΣ ΓΙΑ ΤΗΝ ΠΡΟΣΕΓΓΙΣΗ ΣΤΟ ΤΕΛΟΣ ΤΗΣ ΠΕΡΙΟΧΗΣ ΜΕΛΕΤΗΣ
!INFCHECK		ΔΕΙΚΤΗΣ ΓΙΑ ΤΗΝ ΕΥΡΕΣΗ ΗΜΙΑΠΕΙΡΟΥ ΣΤΟΙΧΕΙΟΥ



C	ΟΡΙΣΜΟΙ ΜΕΤΑΒΛΗΤΩΝ ΚΑΙ ΜΗΤΡΩΩΝ
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

C	ΑΡΧΕΙΑ ΠΟΥ ΧΡΗΣΙΜΟΠΟΙΟΥΝΤΑΙ
	OPEN(7,FILE='input.grd')		!Δεδομένα δικτύου και οριακών συνθηκών
	OPEN(8,FILE='UnRiv1.dat')		!Υδραυλικά και Λοιπά δεδομένα
	OPEN(9,FILE='UnRiv.res')		!Εκτύπωση των δεδομένων του προβλήματος

	OPEN(14,FILE='PHI1_UnRiv.res')	!Τιμές του φορτίου σε κάθε κόμβο
	OPEN(10,FILE='PHI2_UnRiv.res')	!Φορτίο και ύψος πίεσης στις συντεταγμένες
	OPEN(11,FILE='VOL_UnRiv.res')	!Όγκος νερού με το χρόνο
	OPEN(13,FILE='Q1_UnRiv.res')	!Παροχή στους κόμβους προκαθ. φορτίου
	OPEN(20,FILE='Q2_UnRiv.res')	!Συνολική παροχή στο ρέμα με το χρόνο


	OPEN(12,FILE='INF_UnRiv.res')	!Τιμές φορτίου στο τέλος της περιοχής
	OPEN(31,FILE='TEST1_UnRiv.res')	!Αρχείο για δοκιμές
	OPEN(32,FILE='TEST2_UnRiv.res')	!Αρχείο για δοκιμές


C	ΕΙΣΑΓΩΓΗ ΚΑΙ ΕΓΓΡΑΦΗ ΤΩΝ ΔΕΔΟΜΕΝΩΝ
C	Γεωμετρικά δεδομένα
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

C	Χαρακτηριστικοί κόμβοι του δικτύου
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

C	Εισαγωγή αρχικής συνθήκης
	READ(8,1000)HSTART,ERROR1,PHIVAR
	DO 14 I=1,NNP
 14	PHI(I)=HSTART

C	Εισαγωγή οριακής συνθήκης στην περιοχή του ποταμού
	DO 15 I=1,NFIX
 15	PHI(NNX(I))=YY(I)

C	Ειδική παροχή βροχόπτωσης-άρδευσης
	READ(8,1010)RAIN
	DO 16 I=1,NQ
 16	Q(I)=QL(I)*RAIN

C	Τελική τιμή του χρόνου και χρονικό βήμα εκτέλεσης
	READ(8,1010)TSTOP
	READ(8,1010)DT

C	Αριθμός επαναλήψεων σύγκλισης
	READ(8,1030)ITERMAX

C	Παράμετροι Brooks-Corey για κάθε εδαφική στρώση
	READ(8,1020)TSAT1,TRES1,CKSXX1,CKSYY1
	READ(8,*)HB1,CLM1,CN1
	READ(8,1020)TSAT2,TRES2,CKSXX2,CKSYY2
	READ(8,*)HB2,CLM2,CN2

C	Σημεία ολοκλήρωσης και σταθμιστικοί παράγοντες
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

C	Μέθοδος προσέγγισης στο τέλος της περιοχής μελέτης
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

C	Εγγραφή των δεδομένων του προβλήματος
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
	 IF (INFMETHOD==1) METHOD_INF='ΓΡΑΜΜΙΚΗ'
	 IF (INFMETHOD==2) METHOD_INF='ΕΚΘΕΤΙΚΗ'
	  WRITE(9,2090)METHOD_INF,RA,XSTART,YSTART,Xdirec,Ydirec
	ENDIF

	WRITE(10,2100)
 2000 FORMAT('ΑΡΧΕΙΟ ΑΠΟΤΕΛΕΣΜΑΤΩΝ - UnRiv.res'/33('*')/
     $'ΑΡΙΘΜΟΣ ΤΩΝ ΣΤΟΙΧΕΙΩΝ ΤΟΥ ΔΙΚΤΥΟΥ ='I5/
     $'ΑΡΙΘΜΟΣ ΤΩΝ ΚΟΜΒΩΝ ΤΟΥ ΔΙΚΤΥΟΥ    ='I5/
     $'ΕΥΡΟΣ ΤΟΥ ΔΙΚΤΥΟΥ                 ='I5/
     $'ΜΕΘΟΔΟΣ ΑΝΑΡΙΘΜΗΣΗΣ ΤΟΥ ΔΙΚΤΥΟΥ   = ΑΛΓΟΡΙΘΜΟΣ ΤΟΥ 'A8)
 2010 FORMAT(/'α/α ΣΤΟΙΧΕΙΩΝ',5X,'ΚΟΜΒΟΙ ΤΟΥ ΣΤΟΙΧΕΙΟΥ'/
     $(2X,I5,5X,4(2X,I5)))
 2020 FORMAT(/'ΣΥΝΤΕΤΑΓΜΕΝΕΣ ΤΩΝ ΚΟΜΒΩΝ'/3(4X'I',8X'X',8X,'Y',5X)/
     $3(I5,2(1X,F8.2),5X))
 2030 FORMAT(/'ΚΟΜΒΟΙ ΓΝΩΣΤΟΥ ΦΟΡΤΙΟΥ ='1X,I3/
     $2X'ΚΟΜΒΟΣ',5X,'ΦΟΡΤΙΟ'/(3X,I5,3X,F8.3))
 2040 FORMAT(/'ΚΟΜΒΟΙ ΓΝΩΣΤΗΣ ΠΑΡΟΧΗΣ ='1X,I3/
     $2X'ΚΟΜΒΟΣ',3X,'ΜΗΚΟΣ ΡΟΗΣ'/(3X,I5,5X,F8.3))
 2050 FORMAT(
     $/'Αρχική τιμή υδραυλικού φορτίου(m)                  = ',F10.2
     $/'Σφάλμα σύγκλισης σε μονάδες του υδραυλικού φορτίου = ',F10.7
     $/'Μεγιστη μεταβολή του φορτίου σε κάθε χρονικό βήμα  = ',F10.7
     $/'Παροχή άρδευσης ή βροχόπτωσης(m/hr)                = ',F10.7
     $/'Τελική τιμή του χρόνου(hours)                      = ',F10.2     
     $/'Αρχικό χρονικό βήμα εκτέλεσης(hours)               = ',F10.7     
     $/'Αριθμός επαναλήψεων σύγκλισης                      = ',I10)
 2060 FORMAT(//'ΕΔΑΦΙΚΕΣ ΠΑΡΑΜΕΤΡΟΙ ΤΗΣ ΠΑΝΩ ΣΤΡΩΣΗΣ'
     $/'Υγρασία κορεσμού                    = ',F8.5
     $/'Υπολειματική υγρασία                = ',F8.5
     $/'Υδραυλική αγωγιμότητα κορεσμού'
     $/'κατά την οριζόντια διεύθυνση(m/hr)  = ',F8.5
     $/'Υδραυλική αγωγιμότητα κορεσμού'
     $/'κατά την κατακόρυφη διεύθυνση(m/hr) = ',F8.5
     $/'Παράμετροι εδαφικής υγρασίας'
     $/'Hb1(m)=',F8.5,2X,'Ln1=',F8.5,2X,'N1=',F8.5)
 2070 FORMAT(//'ΕΔΑΦΙΚΕΣ ΠΑΡΑΜΕΤΡΟΙ ΤΗΣ ΚΑΤΩ ΣΤΡΩΣΗΣ'
     $/'Υγρασία κορεσμού                    = ',F8.5
     $/'Υπολειματική υγρασία                = ',F8.5
     $/'Υδραυλική αγωγιμότητα κορεσμού'
     $/'κατά την οριζόντια διεύθυνση(m/hr)  = ',F8.5
     $/'Υδραυλική αγωγιμότητα κορεσμού'
     $/'κατά την κατακόρυφη διεύθυνση(m/hr) = ',F8.5
     $/'Παράμετροι εδαφικής υγρασίας'
     $/'Hb2(m)=',F8.5,2X,'Ln2=',F8.5,2X,'N2=',F8.5)	
 2080 FORMAT(//'ΠΑΧΟΣ ΤΗΣ ΚΑΤΩ ΣΤΡΩΣΗΣ'/'YLAY(m)=',F8.5)
 2085 FORMAT(//'ΑΡΙΘΜΟΣ ΣΗΜΕΙΩΝ ΑΡΙΘΜΗΤΙΚΗΣ ΟΛΟΚΛΗΡΩΣΗΣ'/'NGP=',I2)
 2086 FORMAT(//'ΑΡΙΘΜΟΣ ΣΗΜΕΙΩΝ ΑΡΙΘΜΗΤΙΚΗΣ ΟΛΟΚΛΗΡΩΣΗΣ'
     $/'Πεπερασμένων στοιχείων, NGP    =',I2
     $/'Ημιάπειρων στοιχείων,   NGPINF =',I2)
 2090 FORMAT(//'ΜΕΘΟΔΟΣ ΠΡΟΣΕΓΓΙΣΗΣ ΣΤΟ ΤΕΛΟΣ ΤΗΣ ΠΕΡΙΟΧΗΣ ΜΕΛΕΤΗΣ'
     $/'Επιλογή μεθόδου    , INFMETHOD = ',A8
     $/'Τιμή της παραμέτρου, a         =',E9.2
     $/'Αρχή του άξονα -X , XSTART     =',E9.2
     $/'Αρχή του άξονα -Y , YSTART     =',E9.2
     $/'Φορά του άξονα -X , XX         =',E9.2
     $/'Φορά του άξονα -Y , YY         =',E9.2)
 2100 FORMAT('ΑΡΧΕΙΟ ΑΠΟΤΕΛΕΣΜΑΤΩΝ - PHI2_UnRiv.res'/33('*')
     $/'ΤΙΜΕΣ ΦΟΡΤΙΟΥ ΣΤΟΥΣ ΚΟΜΒΟΥΣ'
     $/'ΤΟΥ ΔΙΚΤΥΟΥ ΠΟΥ ΥΠΟΛΟΓΙΣΘΗΚΑΝ')
!
!
!
!ΠΟΡΕΙΑ ΤΩΝ ΥΠΟΛΟΓΙΣΜΩΝ
!
C	1.	Αρχή βρόγχου επαναλήψεων σύγκλισης
C	2.	Αρχή βρόγχου για κάθε στοιχείο
C	3.1	Ανάγνωση των συντεταγμένων των κόμβων του στοιχείου
C	3.2	Υπολογισμός του ύψους πίεσης των κόμβων του στοιχείου
C	3.3	Εύρεση των ημιάπειρων στοιχείων
C	3.4	Επίλυση χωρίς τη χρήση ημιάπειρων στοιχείων
C	3.5	Αριθμητική ολοκλήρωση ανάλογα με το είδος του στοιχείου
C	4.	Αρχή βρόγχων για κάθε σημείο ολοκλήρωσης
C	5.	Υπολογισμός των μητρώων B*,Bm*,xizi,N,M,B,Bm,J - SHAPEFUN
C	6.	Υπολογισμός των C(h) και K(h) - BROOKS_COREY
C	7.	Κλείσιμο βρόγχου για κάθε σημείο ολοκλήρωσης
C	8.	Δημιουργία των μητρώων του στοιχείου	
C	9.	Συνένωση των στοιχείων
C	10.	Κλείσιμο βρόγχου για κάθε στοιχείο
C	11.1	Εισαγωγή οριακής συνθήκης προκαθορισμένου φορτίου
C	11.2	Εισαγωγή παροχής άρδευσης ή βροχόπτωσης
C	12.	Επίλυση του συστήματος των εξισώσεων - GAUSS
C	13.	Έλεγχος του κριτηρίου σύγκλισης
C	14.	Κλείσιμο του βρόγχου σύγκλισης
C	14.1	Μείωση του χρονικού βήματος σε περίπτωση μη σύγκλισης
C		και επιστροφή στο βήμα 1 (GOTO)
C	14.2	Εύρεση της μέγιστης μεταβολής του φορτίου
C	14.3	Έλεγχος της μεταβολής του φορτίου στο τέλος της περιοχής
C		και πιθανή έξοδος από το πρόγραμμα - βήμα 22 (GOTO)
C	15.	Υπολογισμός της παροχής στους κόμβους προκαθορισμένου φορτίου
C	16.	Προστιθέμενος όγκος νερού βάση των τιμών της παροχής στους
C		κόμβους προκαθορισμένου φορτίου και προκαθορισμένης παροχής
C	17.	Υπολογισμός της μεταβολής του όγκου στην περιοχή ροής
C	18.	Εκτύπωση των αποτελεσμάτων της δεδομένης χρονικής στιγμής
C	18.1	Προετοιμασία για την επόμενη χρονική στιγμή
C	19.	Έλεγχος της τελικής τιμής του χρόνου
C		και έξοδος από το πρόγραμμα - βήμα 22 (GOTO)
C	20.	Αύξηση του χρονικού βήματος σε περίπτωση γρήγορης σύγκλισης
C		για την επόμενη χρονική στιγμή
C	21.	Επιστροφή στο βήμα 1 (GOTO)
C	22.	Έξοδος από το πρόγραμμα

!	Υπολογισμός του εμβαδού κάθε στοιχείου
	DO 50 NE=1,NELEM
	AREA(NE)=X(NOD(NE,3))-X(NOD(NE,2))+X(NOD(NE,4))-X(NOD(NE,1))
 50 	AREA(NE)=0.5*AREA(NE)*(Y(NOD(NE,4))-Y(NOD(NE,3)))
!	Υπολογισμός του αρχικά αποθηκευμένου όγκου νερού
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

C	1.	Αρχή βρόγχου επαναλήψεων σύγκλισης
	DO 999 ITER=1,ITERMAX
	IF (ITERCHECK==0) THEN	!ITERCHECK
	AK=0
	BK=0
	AKQ=0
	BKQ=0
C	2.	Αρχή βρόγχου για κάθε στοιχείο
	DO 998 NE=1,NELEM
	INFCHECK=0
	DO 100 I=1,4
C	3.1	Ανάγνωση των συντεταγμένων των κόμβων του στοιχείου
	XY(I,1)=X(NOD(NE,I))
 	XY(I,2)=Y(NOD(NE,I))
C	3.2	Υπολογισμός του ύψους πίεσης των κόμβων του στοιχείου
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
C	3.3 Εύρεση των ημιάπειρων στοιχείων
	IF (XY(I,1)==WOA) INFCHECK=1
C	3.4 Επίλυση χωρίς τη χρήση ημιάπειρων στοιχείων
	IF (INFMETHOD==0) INFCHECK=0
 100	CONTINUE
	
	AK1=0
	AK2=0
	NGP1=0
	RJ1=0
	RW1=0

C	3.5 Αριθμητική ολοκλήρωση ανάλογα με το είδος του στοιχείου
	IF (INFCHECK==0) THEN
	 NGP1=NGP
	 RJ1=RJ
	 RW1=RW
	ELSE
	 NGP1=NGPINF
	 RJ1=RJINF
	 RW1=RWINF
	ENDIF
C	4.	Αρχή βρόγχων για κάθε σημείο ολοκλήρωσης
	DO 997 IN=1,NGP1
	DO 997 JN=1,NGP
C	5.	Υπολογισμός των μητρώων B*,Bm*,xizi,N,M,B,Bm,J - SHAPEFUN
	CALL SHAPEFUN(IN,JN,RJ1,RJ,XY,DJ,SF,B,
     $INFCHECK,SFM,BM,INFMETHOD,RA,XSTART,YSTART,Xdirec,Ydirec)

C	6.	Υπολογισμός των C(h) και K(h) - BROOKS_COREY
	CALL BROOKS_COREY(YLAY,SF,SFM,XY,PG,CKXY,C,WET,
     $TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1,
     $TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2)

	B1=MATMUL(TRANSPOSE(BM),CKXY)	!Υπολογισμός του χωρικού μητρώου
	IF (INFCHECK==0) AK1=AK1+RW1(IN)*RW(JN)*DJ*DT*MATMUL(B1,BM)
	IF (INFCHECK==1) AK1=AK1+0.5*RW1(IN)*RW(JN)*DJ*DT*MATMUL(B1,BM)
	
	DO 150 I=1,4					!Υπολογισμός του χρονικού μητρώου
	IF (INFCHECK==0) AK2(I,I)=AK2(I,I)+RW1(IN)*RW(JN)*C*DJ*SFM(I)
	IF (INFCHECK==1) AK2(I,I)=AK2(I,I)+0.5*RW1(IN)*RW(JN)*C*DJ*SFM(I)
 150	CONTINUE

 997	CONTINUE
C	7.	Κλείσιμο βρόγχου για κάθε σημείο ολοκλήρωσης

C	8.	Δημιουργία των μητρώων του στοιχείου	
	AK1=AK1+AK2
	DO 200 I=1,4
	BK1(I)=0
	DO 200 J=1,4
 200	BK1(I)=BK1(I)+AK2(I,J)*PHI(NOD(NE,J))

C	9.	Συνένωση των στοιχείων
	DO 250 I=1,4
	BK(NOD(NE,I))=BK(NOD(NE,I))+BK1(I)
	DO 250 J=1,4
	IF(NOD(NE,I)<=NOD(NE,J)) THEN		!Λαμβάνεται υπόψη η συμμετρικότητα
	 AK(NOD(NE,I),NOD(NE,J)-NOD(NE,I)+1)=
     $ AK(NOD(NE,I),NOD(NE,J)-NOD(NE,I)+1)+AK1(I,J)
	ENDIF
 250	CONTINUE

C	9.	Αποθήκευση των συντελεστών που απαιτούνται για τον υπολογισμό
C		της παροχής στους κόμβους προκαθορισμένου φορτίου
	DO 251 I=1,4
	BKQ(NOD(NE,I))=BKQ(NOD(NE,I))+BK1(I)
	DO 251 J=1,4
	AKQ(NOD(NE,I),NOD(NE,J))=AKQ(NOD(NE,I),NOD(NE,J))+AK1(I,J)
 251	CONTINUE

 998	CONTINUE
C	10.	Κλείσιμο βρόγχου για κάθε στοιχείο

C	11.1	Εισαγωγή οριακής συνθήκης προκαθορισμένου φορτίου
	DO 300 I=1,NNP
	DO 300 J=1,NFIX
	IF(NNX(J)==I) THEN
	 AK(I,1)=AK(I,1)*10**30
	 BK(I)=YY(J)*AK(I,1)
	ENDIF
 300	CONTINUE

C	11.2	Εισαγωγή παροχής άρδευσης ή βροχόπτωσης
	VOL2=0.
	DO 350 I=1,NNP
	DO 350 J=1,NQ
	IF(IFLOW(J)==I) THEN
	 BK(I)=BK(I)+Q(J)*DT
	 VOL2=VOL2+Q(J)*DT
	ENDIF
 350	CONTINUE




C	------------------------------------------------
C	12.	Επίλυση του συστήματος των εξισώσεων - GAUSS
C	------------------------------------------------
	CALL GAUSS(AK,BK,NNP,NBAND,PHINEW)




C	13.	Έλεγχος του κριτηρίου σύγκλισης
	DMAX1=0.0
	DO 400 I=1,NNP		!Σάρωση όλων των κόμβων
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
C	14.	Κλείσιμο του βρόγχου σύγκλισης

C	14.1	Μείωση του χρονικού βήματος σε περίπτωση μη σύγκλισης
C		και επιστροφή στο βήμα 1 (GOTO)
	IF (ITERCHECK==0) THEN
	 DT=DT*0.5
	 WRITE(*,*)'DT=',DT,' DMAX1=',DMAX1,' ITER=',ITER
	 WRITE(14,2105)DT,DMAX1,IMAX1
	 GOTO 5000
	ENDIF
 2105	FORMAT(/39('*')/'ΤΟ ΧΡΟΝΙΚΟ ΒΗΜΑ ΜΕΙΩΘΗΚΕ ΣΕ:',E10.3/
     $'ΛΟΓΩ ΜΗ ΣΥΓΚΛΙΣΗΣ ΣΤΟ ΦΟΡΤΙΟ:',E10.3/
     $'ΣΤΟΝ ΚΟΜΒΟ:',I10/39('*'))

C	14.2	Εύρεση της μέγιστης μεταβολής του φορτίου
C		και πιθανή επιστροφή στο βήμα 1 (GOTO)
	DMAX=0.0
	DO 410 I=1,NNP		!Σάρωση όλων των κόμβων
	DIFF=ABS(PHINEW(I)-PHI(I))
	IF (DIFF>=DMAX) THEN
	 DMAX=DIFF
	 IMAX=I
	ENDIF
 410	CONTINUE
	!Έλεγχος της μέγιστης επιτρεπόμενης μεταβολής του φορτίου
	!και μείωση του χρονικού βήματος
	IF (DMAX>=PHIVAR) THEN
	 DT=DT*0.5
	 WRITE(*,*)'DT=',DT,' DMAX=',DMAX
	 WRITE(14,2106)DT,DMAX,IMAX
	 GOTO 5000
	ENDIF
 2106	FORMAT(/41('*')/'ΤΟ ΧΡΟΝΙΚΟ ΒΗΜΑ ΜΕΙΩΘΗΚΕ ΣΕ:',E10.3/
     $'ΛΟΓΩ ΜΕΓΑΛΗΣ ΜΕΤΑΒΟΛΗΣ ΦΟΡΤΙΟΥ:',E10.3/
     $'ΣΤΟΝ ΚΟΜΒΟ:',I10/41('*'))

C	14.3	Έλεγχος της μεταβολής του φορτίου στο τέλος της περιοχής
C		και πιθανή έξοδος από το πρόγραμμα - βήμα 22 (GOTO)
	IF (INFMETHOD==0) THEN
	DO 430 I=1,NNP		!Σάρωση όλων των κόμβων
	IF (X(I)==WOA) THEN	!Έλεγχος των κόμβων στο τέλος του δικτύου
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
		
C	15.	Υπολογισμός της παροχής στους κόμβους προκαθορισμένου φορτίου
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

C	16.	Προστιθέμενος όγκος νερού βάση των τιμών της παροχής στους
C		κόμβους προκαθορισμένου φορτίου και προκαθορισμένης παροχής
	VOLUMEQ=VOLQ+QSUM*DT+VOL2


C	17.	Υπολογισμός της μεταβολής του όγκου στην περιοχή ροής		
C		βάση των υπολογιζόμενων τιμών του φορτίου
	CALL VOLUME(PHINEW,VOLUMEPHI,AREA,NELEM,NNP,
     $NOD,X,Y,WOA,NGP,RJ,RW,NGPINF,RJINF,RWINF,
     $INFMETHOD,RA,XSTART,YSTART,Xdirec,Ydirec,
     $YLAY,TSAT1,TRES1,CKSXX1,CKSYY1,HB1,CLM1,CN1,
     $TSAT2,TRES2,CKSXX2,CKSYY2,HB2,CLM2,CN2)


C	18.	Εκτύπωση των αποτελεσμάτων της δεδομένης χρονικής στιγμής
	TIME=TIME+DT
	WRITE(*,*)TIME,DT,ITER1,DMAX
	IF (MOD(TIME,1.)==0.) THEN	!Εκτύπωση ανά μία ώρα
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
	ENDIF						!Εκτύπωση ανά μία ώρα
	WRITE(11,*)TIME,VOLUMEQ,VOLUMEPHI
	WRITE(20,*)TIME,QSUM

 2110 FORMAT(///16X,'TIME(hours)=',E10.3/
     $16X,'DT(hours)  =',E10.3/
     $16X,'ITER       =',I10/
     $8X,'Μέγιστη μεταβολή φορτίου(m)=',E10.3,1X,'στον κόμβο=',I5/
     $8X,'Μέγιστο σφάλμα σύγκλισης(m)=',E10.3,1X,'στον κόμβο=',I5/
     $73('-')/)
 2120 FORMAT(2(1X,F10.3),(1X,F10.3))
 2121 FORMAT(2(1X,F10.3),2(1X,F10.3))
 2130 FORMAT(2(1X,F10.3),2(1X,F10.5))
 2140 FORMAT(5(1X,I5,':',F6.3,2X))
 2150 FORMAT('ΠΑΡΟΧΗ='1X,E8.3,2X,'ΣΤΟ ΧΡΟΝΟ= ',F6.3)

C	18.1	Προετοιμασία για την επόμενη χρονική στιγμή
	DTOLD=DT
	PHIOLD=PHI
	PHI=PHINEW
	VOLQ=VOLUMEQ
	VOLPHI=VOLUMEPHI

C	19.	Έλεγχος της τελικής τιμής του χρόνου
C		και έξοδος από το πρόγραμμα - βήμα 22 (GOTO)
	IF (TIME>=TSTOP) GOTO 6000

C	20.	Μεταβολή του χρονικού βήματος για την επόμενη χρονική στιγμή
	DT1=DTOLD*PHIVAR/DMAX	!Βάση της μεταβολής φορτίου
	DT2=DTOLD*1.2			!Αύξηση κατά 20%
	DT3=AINT(TIME)+1-TIME	!Στρογγυλοποίηση ανά ώρα
	DT=MIN(DT1,DT2,DT3)		!Επιλογή χρονικού βήματος από τα παραπάνω

C	21.	Επιστροφή στο βήμα 1 (GOTO)
	GOTO 5000

C	22.	Έξοδος από το πρόγραμμα
 6000	CONTINUE


							!Υπολογιστικοί χρόνοι που απαιτήθηκαν
	CALL CPU_TIME(TIMEEND)
	TIMEEQ=TIMEEND-TIMEBEGIN
	TIMETOTAL=TIMERENUM+TIMEEQ
	WRITE(9,2200)TIMERENUM,TIMEEQ/60,TIMETOTAL/60
 2200 FORMAT(//'ΥΠΟΛΟΓΙΣΤΙΚΟΙ ΧΡΟΝΟΙ ΠΟΥ ΑΠΑΙΤΗΘΗΚΑΝ'/
     $'Για την αναρίθμηση του δικτύου (seconds) ='F10.3/
     $'Για την επίλυση της εξίσωσης   (minutes) ='F10.3/
     $'Συνολικός υπολογιστικός χρόνος (minutes) ='F10.3)
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
C	ΥΠΟΛΟΓΙΖΕΤΑΙ Ο ΑΠΟΘΗΚΕΥΜΕΝΟΣ ΟΓΚΟΣ ΝΕΡΟΥ
C	ΒΑΣΗ ΤΩΝ ΥΠΟΛΟΓΙΖΟΜΕΝΩΝ ΤΙΜΩΝ ΤΟΥ ΦΟΡΤΙΟΥ
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
C	ΥΠΟΛΟΓΙΖΟΝΤΑΙ ΓΙΑ ΚΑΘΕ ΣΗΜΕΙΟ ΟΛΟΚΛΗΡΩΣΗΣ
C	ΟΙ ΣΥΝΑΡΤΗΣΕΙΣ ΜΟΡΦΗΣ (SF)
C	ΟΙ ΠΑΡΑΓΩΓΟΙ ΑΥΤΩΝ (BS)
C	ΤΟ ΙΑΚΩΒΙΑΝΟ ΜΗΤΡΩΟ (AJ) ΚΑΙ Η ΟΡΙΖΟΥΣΑ ΑΥΤΟΥ (DJ)
C	ΤΟ ΜΗΤΡΩΟ Β
	DOUBLE PRECISION SF,SFM,BS,BSM,B,BM
	DIMENSION RJ(10),RJ2(10),SF(4),BS(2,4),B(2,4),XY(4,2),AJ(2,2)
	DIMENSION SFM(4),BSM(2,4),BM(2,4),SR(4)
	DIMENSION SF1(4),BS1(2,4),PHI(4),PG(4)

	IF (INFCHECK==0) THEN	!Κανονικό στοιχείο
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
	ENDIF					!Κανονικό στοιχείο



	IF (INFCHECK==1) THEN	!Ημιάπειρο στοιχείο
	RIN1=1.-RJ(IN)
	RIN2=1.+RJ(IN)
	RJN1=1.-RJ2(JN)
	RJN2=1.+RJ2(JN)

	IF (RJ(IN)<=0) THEN		!Πεπερασμένο τμήμα του ημιάπειρου στοιχείου
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
	ENDIF					!Πεπερασμένο τμήμα του ημιάπειρου στοιχείου

	IF (RJ(IN)>0) THEN		!Απειρίζον τμήμα του ημιάπειρου στοιχείου
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
	ENDIF					!Απειρίζον τμήμα του ημιάπειρου στοιχείου
	ENDIF					!Ημιάπειρο στοιχείο

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
C	ΥΠΟΛΟΓΙΖΟΝΤΑΙ ΓΙΑ ΚΑΘΕ ΣΗΜΕΙΟ ΟΛΟΚΛΗΡΩΣΗΣ
C	ΑΝΑΛΟΓΑ ΜΕ ΤΗΝ ΕΔΑΦΙΚΗ ΣΤΡΩΣΗ ΠΟΥ ΒΡΙΣΚΕΤΑΙ ΚΑΙ
C	ΑΝΑΛΟΓΑ ΜΕ ΤΗΝ ΥΓΡΑΣΙΑΚΗ ΚΑΤΑΣΤΑΣΗ ΤΟΥ (ΚΟΡΕΣΜΕΝΗ-ΑΚΟΡΕΣΤΗ)
C	ΤΟ ΜΗΤΡΩΟ ΤΗΣ ΥΔΡΑΥΛΙΚΗΣ ΑΓΩΓΙΜΟΤΗΤΑΣ (CKXY)
C	Η ΥΔΡΑΥΛΙΚΗ ΧΩΡΗΤΙΚΟΤΗΤΑ (C)
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
C	ΕΠΙΛΥΕΤΑΙ ΣΥΜΦΩΝΑ ΜΕ ΤΟΝ ΑΛΓΟΡΙΘΜΟ ΤΟΥ GAUSS
C	ΤΟ ΣΥΣΤΗΜΑ ΤΩΝ ΕΞΙΣΩΣΕΩΝ ΛΑΜΒΑΝΟΝΤΑΣ ΥΠΟΨΗ ΤΟ ΕΥΡΟΣ ΤΟΥ ΔΙΚΤΥΟΥ
C	ΚΑΙ ΕΠΙΣΤΡΕΦΟΝΤΑΙ ΥΠΟΛΟΓΙΖΟΜΕΝΕΣ ΤΙΜΕΣ ΤΟΥ ΦΟΡΤΙΟΥ (PHINEW)
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