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

int main(int argc, char **argv)
{
int count,i;
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
for(t=0;t<totaltime;t=t+delt)
{
printf("%d Current time step is %lf\n",ml,t);
forcing(t);   // Enter forcing function 
stressupdate();  // update stresses
velocityupdate(); // update velocities
boundarycns();  // apply boundary conditions
storewave();// storing recorded signals
count = count+1;
//printf("%3.50lf \n",U3[mat(23,24,4)]);
}
printf("code is good %d \n",count2);
time2 = omp_get_wtime();
printf("time is %3.50lf \n",time2-time1); //prints time take for one shot
write_to_file(forward,count2); // writing to disk the stored velocities
grvy_timer_end("Main Program");
grvy_timer_finalize();
grvy_timer_summarize();
const int num_procs = 1;

  grvy_timer_save_hist("C-Example1","My clever comment",num_procs,"hist.h5");
exit(1);
}
}




























