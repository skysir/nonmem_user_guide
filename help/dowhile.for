


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           DOWHILE BLOCK                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Coding technique
 CONTEXT: Abbreviated code

 DISCUSSION:

 SAMPLE:
 $ABBR DECLARE DOWHILE ILOOP $PK
 ILOOP=1
 DOWHILE (condition)
  .. statements ..
 ILOOP=ILOOP+1
 ENDDO

 DISCUSSION:

 DOWHILE blocks allow loops in abbreviated code.  Previously, they were
 permitted only in initialization-finalization, simulation, data  aver-
 age,  and  expectation blocks.  See the separate descriptions of these
 blocks.  With NONMEM 7.3, DOWHILE may also be  used  in  data-analytic
 blocks (i.e., in code that is executed at ICALL=2).

 When  it  is  used  in  data-anlytic code, DOWHILE must be used with a
 looping variable.  Such a variable meets three criteria:
 (1) It is declared as such, e.g.
 $ABBR DECLARE DOWHILE ILOOP
 (2) It is given a value in the statement immediately before DOWHILE
 (3) It is given a value in the statement immediately before ENDDO.

 Usually, the looping variable will be tested  in  the  condition,  but
 this is not necessary.

 DOWHILE blocks may be nested.

 DOWHILE  blocks  may  not  include IF/THEN/ENDIF blocks.  Only single-
 statement IF statements are permitted.  E.g.,
 IF (condition) statement

 Subscripted variables may be used within the block.  These  are  vari-
 ables  that  have  been  declared arrays using the $ABBR statement, or
 subscripted reserved variables such as THETA.  Subscripts may be inte-
 ger  expressions  using declared integer variables and integer values.
 E.g.
 $ABBR DECLARE X(10)
 $ABBR DECLARE DOWHILE ILOOP
  ...
 $PK
 ILOOP=1
 DOWHILE (condition)
 X(ILOOP)=THETA(ILOOP+1)
 ILOOP=ILOOP+1
 ENDDO

 (See abbreviated code).

 A DOWHILE loop may compute  recursive  random  variables.   These  are
 variables  that  are modified recursively in a dowhile block and which
 have eta derivatives. For example,

 TERM=THETA(1)*EXP(ETA(1))
 SUM=0
 ILOOP=1
 DO WHILE(ILOOP<=IMAX)
 SUM=SUM+TERM
 ILOOP=ILOOP+1
 ENDDO

 The dowhile recursive variable must be  initialized  to  a  non-random
 variable  outside  the loop and must appear on both sides of the equal
 sign only once within the DOWHILE loop: V= ... V ...

 The syntax is very limited.  The following are not permitted:
 ILOOP=1
 DO WHILE(ILOOP<=IMAX)
 IF (ILOOP == 1) SUM=0  ; initialization within the loop
 SUM=SUM+TERM
 ILOOP=ILOOP+1
 ENDDO

 ILOOP=1
 DO WHILE(ILOOP<=IMAX)
 IF (ILOOP == 1) THEN
 SUM=0                  ; initialization within the loop
 ELSE                   ; else statement
 SUM=SUM+TERM
 ENDIF
 ILOOP=ILOOP+1
 ENDDO

 Although the example given is a summing loop, a product loop such as
   PROD=PROD*TERM
 is also possible, as are other ways the dowhile recursive variable can
 be  used.   The  looping variable ILOOP may also be set and tested and
 incremented with different integers than shown.

 EXAMPLES:

 Several examples are present in the examples directory:

 Auto-correlation (See ar1mod.ctl)
 Auto-correlation with Simulation (See ar1newsim.ctl).
 Dose Superposition (See sumdosetn.ctl)
 Summing elements of THETA (See superid3_6.ctl)

 REFERENCES: None.
