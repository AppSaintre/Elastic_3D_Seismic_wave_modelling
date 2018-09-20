#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>

#define array (Ny*Nz*i)+(Nz*j)+k
#define mat(i,j,k) (Ny*Nz*i)+(Nz*j)+k
void write_to_file(double *ptr, int size)
{
int i,j;
if(deb == 1) 
   printf("Writing wavefields recorded to disk\n");
 grvy_timer_begin(__func__);
FILE *fp;
fp=fopen("output","wb");if(!fp){printf("Can't open file \n");exit(1);}
fwrite(ptr,sizeof(double),size,fp);
 grvy_timer_end(__func__);
fclose(fp);
}



 





