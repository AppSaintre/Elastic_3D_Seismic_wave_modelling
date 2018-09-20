#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>
#include<unistd.h>
#include<grvy.h>
#include<sys/time.h>
#include<time.h>
#include"omp.h"
#include"globalvar.h"
#include"damping3d.h"
#include"routines.h"
#include"write.h"
#include"test.h"

int main(int argc, char **argv)
{
int count,i,j,k;
double t,time1,time2;
allocations();  //Dynamic allocatins of all variables
grvy_timer_init("GRVY Example Timing");
for(ml=1;ml<=1;ml++)  // loop for 6 shots
{
grvy_timer_reset();
grvy_timer_begin("Main Program");
time1 = omp_get_wtime();
count2 = 0;
count =0;
//FILE *matrix = fopen("check3.txt","w");
for(t=0;t<totaltime;t=t+delt)
{
printf("%d Current time step is %lf\n",ml,t);
forcing(t);   // Enter forcing function 
stressupdate();  // update stresses
velocityupdate(); // update velocities
boundarycns();  // apply boundary conditions
storewave();// storing recorded signals
//for(i=22;i<=221;i++)  //storing velocities at receivers
//{
//fprintf(matrix,"%3.50lf \n",W3[mat(25,24,4)]);
//}
count = count+1;
}
//fclose(matrix);
if (modev == 1) test();
printf("code is good %d \n",count2);
write_to_file(forward,count2); // writing to disk the stored velocities
time2 = omp_get_wtime();
printf("time is %3.50lf \n",time2-time1); //prints time take for one shot
grvy_timer_end("Main Program");
grvy_timer_finalize();
grvy_timer_summarize();
const int num_procs = 1;
grvy_timer_save_hist("C-Example1","My clever comment",num_procs,"hist.h5");
exit(1);
}
}




























