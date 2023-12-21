


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $OMIT                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Defines data item types to be excluded from template matching
 when raw data averages are computed
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $OMIT  item1  item2  item3 ...

 SAMPLE:
 $OMIT   TIME

 DISCUSSION:
 Optional.  If a label of a data item type listed in the $INPUT record,
 or  a  synonym for such a data item type, appears in the $OMIT record,
 then data items of this type are excluded from template matching.

 REFERENCES: None
