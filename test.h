#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"omp.h"
#include <string.h>

#define array (Ny*Nz*i)+(Nz*j)+k
#define mat(i,j,k) (Ny*Nz*i)+(Nz*j)+k

void test()
{
int count,count1;
count1 =0;
double t;
a = 0; 
h = 8;
int i,j;
for(ml=1;ml<=3;ml++)  // loop for 6 shots
{
count =0;
for(t=0;t<(totaltime/2);t=t+delt)
{
printf("%d Current time step is %lf\n",ml,t);
forcing(t);   // Enter forcing function 
stressupdate();  // update stresses
velocityupdate(); // update velocities
boundarycns();  // apply boundary conditions
count = count+1;
}
count1 = count1 +1;
printf ("Code succesful for 2nd order. Success count = %d",count1);
h= h/2;
}
a = -1/24.0; 
h = 8;
count1 = 0; 
for(ml=1;ml<=3;ml++)  // loop for 6 shots
{
count =0;
for(t=0;t<(totaltime/2);t=t+delt)
{
printf("%d Current time step is %lf\n",ml,t);
forcing(t);   // Enter forcing function 
stressupdate();  // update stresses
velocityupdate(); // update velocities
boundarycns();  // apply boundary conditions
count1 = count1 +1;
printf ("Code succesful for 4rth order. Success count = %d",count1);
h= h/2;
count = count+1;
}
h= h/2;
}


}















