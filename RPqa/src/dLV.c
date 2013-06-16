/* Version 0.1
 Date: 10-09-2012
 */

#define INDEX(a, b, nr) ((b)*(nr) + (a))

void  dLV(int  *D, int  *N, int *T, int *histLV)
{
int n=*N;
int t=*T;
int count, i, j;

// IMPORTANT: the LOI must go from [N,1] to [1, N]
// The input matrix must be square!!
// but this can be easily changed

// First calculate L lines
for(i = 0; i< ((n-1)-t); i++)
{
 count = 0;
 for(j= 0; j<=i; j++)
 {
   if(D[INDEX((i-j), j, n)] == 1) count++; else {
      histLV[0]++ ;
	  if (count > 0) {
	         histLV[count]++ ;
	         count = 0;}
	  } // if
 } //for
 if (count > 0) histLV[count]++;
} //for

// Calculate V
for(i = 0; i< ((n-1)-t); i++)
{
 count = 0;
 for(j= 0; j< ((n-1)-i-t); j++)
 {
   if(D[INDEX(i, j, n)] == 1) count++; else {
      histLV[0+n]++ ;
	  if (count > 0) {
	         histLV[count+n]++ ;
	         count = 0;}
	  } // if
 } //for
 if (count > 0) histLV[count+n]++;
} //for
}

