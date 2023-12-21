; Pre-Control stream template
;DDE
;TAU1=75.0
;TAU2=42.0
;TAU3=147.0
;TAU4=114.0
;TAU5=1587.0
;TAU6=1554.0
;TSTOP=200.0
$SIZES PC=80 DIMNEW=3000 DIMTMP=1000 PG=200
$PROB EPO
$ABBR FUNCTION DDEFUNC(VDDE,4)
$ABBR DERIV2=NO DERIV2=NOCOMMON ; DERIV1=NO
$INPUT ID AMT TIME DV EVID MDV  CMT
$DATA EPO.csv IGNORE=C
$SUBROUTINES ADVAN16  TOL=9 ATOL=9
$MODEL NCOMPARTMENTS=6

$PK

KON  =THETA(1)+ETA(1)
KOFF =THETA(2)
KEL  =THETA(3)
KPT  =THETA(4)
KTP  =THETA(5)
VP   =THETA(6)
KINT =THETA(7)
SMAX =THETA(8)
SC50 =THETA(9)
IMAX =THETA(10)
IC50 =THETA(11)
MCH  =THETA(12)
C0   =THETA(13)
RR0   =THETA(14)
KDEG =THETA(15)
RBC0 =THETA(16)
TP1  =THETA(17)
TP2  =THETA(18)
TRET =THETA(19)
TRBC =THETA(20)

;TAUy
TAU1=TP1+TP2
TAU2=TP2
TAU3=TP1+TP2+TRET
TAU4=TP2+TRET
TAU5=TP1+TP2+TRET+TRBC
TAU6=TP2+TRET+TRBC

RET0=TRET*RBC0/(TRET+TRBC)
RBCM0=RBC0-RET0
HB0=MCH*RBC0
AT0=KPT*C0*VP/KTP
RC0=KON*RR0*C0/(KOFF+KINT)
KEPO=KEL*C0*VP+KINT*RC0*VP
KSYN=KDEG*RR0+KINT*RC0
KIN=RET0/(TRET*(1+SMAX*RC0/(SC50+RC0))**2)

; Initial conditions for Base equations
A_0(1)=C0*VP
A_0(2)=AT0
A_0(3)=RR0
A_0(4)=RC0
A_0(5)=RET0
A_0(6)=RBCM0


TSTOP=200.0
; INITIALIZING EQUATIONS FOR DDE COMPARTMENTS
TAU_1=1*TAU1
TAU_14=1*TAU2
TAU_15=1*TAU3
TAU_21=1*TAU4
TAU_22=1*TAU5
TAU_23=1*TAU6

$DES
; AD_x_y is the State value of A(x) delayed for time TAUy.
; AP_x_y is the State value of A(x) in the past, for time delay TAUy.
; PASTS

; DELAY SETUP FOR EQUATION SET 1
 AP_4_1=RC0
 AP_5_1=RET0
 AP_6_1=RBCM0
 DTAU_1=0.0
 IF(T>=TAU_1) DTAU_1=1.0
VDDE(1)=TAU_1
VDDE(2)=T
VDDE(3)=4
VDDE(4)=1
 AZ_4_1=DDEFUNC(VDDE)
 AD_4_1=(1.0-DTAU_1)*AP_4_1+DTAU_1*AZ_4_1
VDDE(1)=TAU_1
VDDE(2)=T
VDDE(3)=5
VDDE(4)=1
 AZ_5_1=DDEFUNC(VDDE)
 AD_5_1=(1.0-DTAU_1)*AP_5_1+DTAU_1*AZ_5_1
VDDE(1)=TAU_1
VDDE(2)=T
VDDE(3)=6
VDDE(4)=1
 AZ_6_1=DDEFUNC(VDDE)
 AD_6_1=(1.0-DTAU_1)*AP_6_1+DTAU_1*AZ_6_1

; DELAY SETUP FOR EQUATION SET 14
 AP_4_2=RC0
 DTAU_14=0.0
 IF(T>=TAU_14) DTAU_14=1.0
VDDE(1)=TAU_14
VDDE(2)=T
VDDE(3)=4
VDDE(4)=2
 AZ_4_2=DDEFUNC(VDDE)
 AD_4_2=(1.0-DTAU_14)*AP_4_2+DTAU_14*AZ_4_2

; DELAY SETUP FOR EQUATION SET 15
 AP_4_3=RC0
 AP_5_3=RET0
 AP_6_3=RBCM0
 DTAU_15=0.0
 IF(T>=TAU_15) DTAU_15=1.0
VDDE(1)=TAU_15
VDDE(2)=T
VDDE(3)=4
VDDE(4)=3
 AZ_4_3=DDEFUNC(VDDE)
 AD_4_3=(1.0-DTAU_15)*AP_4_3+DTAU_15*AZ_4_3
VDDE(1)=TAU_15
VDDE(2)=T
VDDE(3)=5
VDDE(4)=3
 AZ_5_3=DDEFUNC(VDDE)
 AD_5_3=(1.0-DTAU_15)*AP_5_3+DTAU_15*AZ_5_3
VDDE(1)=TAU_15
VDDE(2)=T
VDDE(3)=6
VDDE(4)=3
 AZ_6_3=DDEFUNC(VDDE)
 AD_6_3=(1.0-DTAU_15)*AP_6_3+DTAU_15*AZ_6_3

; DELAY SETUP FOR EQUATION SET 21
 AP_4_4=RC0
 DTAU_21=0.0
 IF(T>=TAU_21) DTAU_21=1.0
VDDE(1)=TAU_21
VDDE(2)=T
VDDE(3)=4
VDDE(4)=4
 AZ_4_4=DDEFUNC(VDDE)
 AD_4_4=(1.0-DTAU_21)*AP_4_4+DTAU_21*AZ_4_4

; DELAY SETUP FOR EQUATION SET 22
 AP_4_5=RC0
 AP_5_5=RET0
 AP_6_5=RBCM0
 DTAU_22=0.0
 IF(T>=TAU_22) DTAU_22=1.0
VDDE(1)=TAU_22
VDDE(2)=T
VDDE(3)=4
VDDE(4)=5
 AZ_4_5=DDEFUNC(VDDE)
 AD_4_5=(1.0-DTAU_22)*AP_4_5+DTAU_22*AZ_4_5
VDDE(1)=TAU_22
VDDE(2)=T
VDDE(3)=5
VDDE(4)=5
 AZ_5_5=DDEFUNC(VDDE)
 AD_5_5=(1.0-DTAU_22)*AP_5_5+DTAU_22*AZ_5_5
VDDE(1)=TAU_22
VDDE(2)=T
VDDE(3)=6
VDDE(4)=5
 AZ_6_5=DDEFUNC(VDDE)
 AD_6_5=(1.0-DTAU_22)*AP_6_5+DTAU_22*AZ_6_5

; DELAY SETUP FOR EQUATION SET 23
 AP_4_6=RC0
 DTAU_23=0.0
 IF(T>=TAU_23) DTAU_23=1.0
VDDE(1)=TAU_23
VDDE(2)=T
VDDE(3)=4
VDDE(4)=6
 AZ_4_6=DDEFUNC(VDDE)
 AD_4_6=(1.0-DTAU_23)*AP_4_6+DTAU_23*AZ_4_6

; DELAY EQUATIONS FOR EQUATION SET 0 (BASE EQUATIONS)

 CC=A(1)/VP
 AT=A(2)
 RR=A(3)
 RC=A(4)
 RET=A(5)
 RBCM=A(6)

 X1=1+SMAX*AD_4_1/(SC50+AD_4_1)
 X2=1+SMAX*AD_4_2/(SC50+AD_4_2)
 X3=1+SMAX*AD_4_3/(SC50+AD_4_3)
 X4=1+SMAX*AD_4_4/(SC50+AD_4_4)
 X5=1+SMAX*AD_4_5/(SC50+AD_4_5)
 X6=1+SMAX*AD_4_6/(SC50+AD_4_6)
 X0=1+SMAX*RC0/(SC50+RC0)

 I1=1-IMAX*(MCH*(AD_5_1+AD_6_1)-HB0)/(IC50+(MCH*(AD_5_1+AD_6_1)-HB0))
 I3=1-IMAX*(MCH*(AD_5_3+AD_6_3)-HB0)/(IC50+(MCH*(AD_5_3+AD_6_3)-HB0))
 I5=1-IMAX*(MCH*(AD_5_5+AD_6_5)-HB0)/(IC50+(MCH*(AD_5_5+AD_6_5)-HB0))



; BASE EQUATIONS.

 DADT(1)=KEPO-KON*CC*VP*RR+KOFF*RC*VP-(KEL+KPT)*CC*VP+KTP*AT
 DADT(2)=KPT*CC*VP-KTP*AT
 DADT(3)=KSYN-KON*CC*RR+KOFF*RC-KDEG*RR
 DADT(4)=KON*CC*RR-(KOFF+KINT)*RC
 DADT(5)=KIN*X1*X2*I1-KIN*X3*X4*I3
 DADT(6)=KIN*X3*X4*I3-KIN*X5*X6*I5



$ERROR

Y1=A(1)/VP
Y2=A(2)
Y3=A(3)
Y4=A(4)
Y5=A(5)
Y6=A(5)+A(6)
Y7=MCH*Y6
IF(CMT==1) IPRED=Y1*(1.0+EPS(1))
IF(CMT==2) IPRED=Y2
IF(CMT==3) IPRED=Y3
IF(CMT==4) IPRED=Y4
IF(CMT==5) IPRED=Y5
IF(CMT==6) IPRED=Y6

Y=IPRED

$THETA
0.01132    ; 1: KON    1/nM/h
1.297      ; 2: KOFF   1/h
0.2256     ; 3: KEL    1/h
0.2092     ; 4: KPT    1/h
0.1721     ; 5: KTP    1/h
0.05694    ; 6: VP     mL/kg
0.8228     ; 7: KINT   1/h
3.48       ; 8: SMAX   1
1.7        ; 9: SC50   pM
1.0        ; 10: IMAX  1
1.79       ; 11: IC50  g/dL
2.0        ; 12: MCH   0.1 g/dL
3.248      ; 13: C0    pM
63.2       ; 14: RR0    pM
0.1133     ; 15: KDEG  1/h
6.128      ; 16: RBC0  10^6 cells/uL
42.97      ; 17: TP1   h
33.6       ; 18: TP2   h
72.33      ; 19: TRET  h
1440       ; 20: TRBC  h

$OMEGA (0.0 FIXED)
$SIGMA (0.0 FIXED)

$SIMULATION (567811 NORMAL) (2933012 UNIFORM) ONLYSIMULATION SUBPROBLEMS=1
$TABLE ID TIME IPRED CMT MDV Y7
NOAPPEND NOPRINT ONEHEADER FILE=EPOd.tab