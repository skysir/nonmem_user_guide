;Model Desc: Receptor Mediated Clearance model with Dynamic Change in Receptors
;Project Name: antibody
;Project ID: NO PROJECT DESCRIPTION

;Project ID: NO PROJECT DESCRIPTION

$PROB RUN# wexample6x (from r2compl), using FAST
$ABBR DERIV2=NO DERIV2=NOCOMMON ; DERIV1=NO
$INPUT C SET ID JID TIME DV AMT RATE EVID MDV CMT
$DATA a_uflg.csv IGNORE=C IGNORE=(AMT.EQN.1.0) IGNORE=(RATE.EQN.1.0E+08)

; The new numerical integration solver is used, although ADVAN=9 is also efficient
; for this problem.
$SUBROUTINES ADVAN13 TRANS1 TOL=6 ATOL=6
$MODEL NCOMPARTMENTS=3

;Initial Omegas
$OMEGA BLOCK(8) VALUES(0.5,0.01)

$PK
MU_1=THETA(1)
MU_2=THETA(2)
MU_3=THETA(3)
MU_4=THETA(4)
MU_5=THETA(5)
MU_6=THETA(6)
MU_7=THETA(7)
MU_8=THETA(8)
VC=EXP(MU_1+ETA(1))
K10=EXP(MU_2+ETA(2))
K12=EXP(MU_3+ETA(3))
K21=EXP(MU_4+ETA(4))
VM=EXP(MU_5+ETA(5))
KMC=EXP(MU_6+ETA(6))
K03=EXP(MU_7+ETA(7))
K30=EXP(MU_8+ETA(8))
S3=VC
S1=VC
KM=KMC*S1
A_0(3)=K03/K30
MTIME(1)=7.1
MTDIFF=1
AZTEST=A_0FLG
IF(TSTATE==MTIME(1).AND.AZTEST==0) A_UFLG=1
A_U(1)=A(1)+1000.0
A_U(2)=A(2)
A_U(3)=A(3)

$DES
DADT(1) = -(K10+K12)*A(1) + K21*A(2) - VM*A(1)*A(3)/(A(1)+KM)
DADT(2) = K12*A(1) - K21*A(2)
DADT(3) =  -(VM-K30)*A(1)*A(3)/(A(1)+KM) - K30*A(3) + K03

$ERROR
ETYPE=1
IF(CMT.NE.1) ETYPE=0
CP=A(1)/S1
CR=A(3)/S3
IPRE=CP
IF(CMT.NE.1) IPRE=CR
Y = IPRE + IPRE*ETYPE*EPS(1) + IPRE*(1.0-ETYPE)*EPS(2)


$THETA 
;Initial Thetas
( 4.0 )  ;[MU_1]
( -3.1 ) ;[MU_2]
( 0.5 )  ;[MU_3]
( -0.2 );[MU_4]      
( 3.2 ) ;[MU_5]
( 0.01 )  ;[MU_6]
( 4.0 )  ;[MU_7]
( -0.1) ;[MU_8]



$SIGMA  
0.1 ;[p]
0.1 ;[p]

$EST METHOD=1 INTERACTION PRINT=1 NOABORT NOPRIOR=1 SIGL=6 MAXEVAL=0 NSIG=2 MCETA=10 
     NOTHETABOUNDTEST NOOMEGABOUNDTEST NOSIGMABOUNDTEST FAST
