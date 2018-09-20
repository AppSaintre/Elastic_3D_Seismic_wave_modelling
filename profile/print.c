#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>

#define array (Ny*Nz*i)+(Nz*j)+k
#define mat(i,j,k) (Ny*Nz*i)+(Nz*j)+k
int main()
{
int i,j;
int k1 = 4001;
int k2 = 600;
FILE *fp;
double * forward = malloc(245*45*225*sizeof(double));
fp=fopen("output","r");if(!fp){printf("Can't open file \n");exit(1);}
fread(forward,sizeof(double),k1*k2,fp);
int count2  = 0;

for (i = 0; i <  k1; i++)
{
      for (j = 0; j < k2; j++)
         {

printf("%3.50lf	 ",forward[count2]);
count2 = count2 +1;
}
printf(" \n");
}



}
