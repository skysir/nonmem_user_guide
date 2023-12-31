


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $CHAIN                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Supplies initial estimates for an entire problem
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $CHAIN  [FILE=filename]
         [FORMAT=s] [ORDER=xxxf]
         [NOTITLE=[0|1]] [NOLABEL=[0|1]]
         [ISAMPLE=n][NSAMPLE=n]
         [SEED=n]
         [CLOCKSEED=[0|1]]
         [SELECT=n]
         [RANMETHOD=[n|S|m] ]
         [CTYPE=[0|1|2|3|4]]
         [DF=n]
         [DFS=n]
         [IACCEPT=x]
         [TBLN=n]

 SAMPLE:
 $CHAIN FILE=example1_previous.txt NSAMPLE=0 ISAMPLE=-1000000000

 DISCUSSION:
 Option  METHOD=CHAIN  of  the $ESTIMATION record will only set thetas,
 omegas, and sigmas for initial values of the estimation process.   Its
 scope  is  therefore limited in that it will not impact the parameters
 used in simulating data for the Simulation step. To introduce  initial
 THETAs  omegas  and sigmas that will cover the entire scope of a given
 problem,  use the $CHAIN record.

 The options have the same meanings  as  for  the  $ESTIMATION  record.
 Setting  SEED  or  CLOCKSEED  or RANMETHOD in a $CHAIN record does not
 propagate to $EST METHOD=CHAIN or any other $EST record.
 (See $ESTIMATION_record_options).
 (See $ESTIMATION_record).

 REFERENCES: Guide Introduction_7
