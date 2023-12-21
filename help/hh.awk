{ x[NR]=$0 } 
  END { h=int(NR/2) ; if (h*2!=NR) h=h+1
  for (i=1; i<=h; i++) 
 {printf("%-39.39s %-39.39s\n", x[i],x[i+h])} }
