


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           RANDOM MODELS                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Models for random variables
 CONTEXT: Abbreviated code

 DISCUSSION:
 The  following are examples of commonly used models using random vari-
 ables.

 Models for CL in terms of TVCL ("typical value of clearance") and  ETA
 are  examples  of  models expressing population inter-individual vari-
 ability.

 Models for Y involving F ("prediction  based  on  the  pharmacokinetic
 parameters")  and  ERR  are  examples  of  models for intra-individual
 ("residual") variability. They are used with both population  or  sin-
 gle-subject  data,  in  which  case ERR stands for EPS or ETA, respec-
 tively.

 Additive models

      CL=TVCL+ETA(1)
      Y=F+ERR(1)

 Proportional (CCV; Constant Coefficient of Variation) models

      CL=TVCL*(1+ETA(1))
      Y=F*(1+ERR(1))

      An equivalent way of coding the proportional model is:

      CL=TVCL+TVCL*ETA(1)
      Y=F+F*ERR(1)

 Exponential models

      CL=TVCL*EXP(ETA(1))
      Y=F*EXP(ERR(1))

      During estimation by  the  first-order  method,  the  exponential
      model  and proportional models give identical results, i.e., NON-
      MEM cannot distinguish between them.

      During estimation by a conditional estimation method,  the  expo-
      nential  and proportional models for inter-individual variability
      give different results.  During simulation, the two  models  give
      different results, in both the inter- and intra-individual cases.

 Power Function model

      Y=F+F**P*ERR(1)

      The  Power Function model has both the additive and the CCV error
      models as special cases, and smoothly interpolates  between  them
      in other cases.

 Combined Additive and Proportional model (slope-intercept model)

      Y =F*(1+ERR(1))+ERR(2)

      Here  is  an alternative parameterization for the same model when
      there is no covariance between ERR(1) and ERR(2).  Any theta  may
      be used.

      W=(1+THETA(5)*THETA(5)*F*F)**.5
      Y=F+W*ERR(1)

 REFERENCES:  Guide V, section 3 , 4.1 , 7.5 , 8.3 
 REFERENCES: Guide V, section 8 
 REFERENCES: Guide VII, section I , III 
