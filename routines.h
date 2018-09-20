#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"omp.h"
#include <string.h>

#define array (Ny*Nz*i)+(Nz*j)+k
#define mat(i,j,k) (Ny*Nz*i)+(Nz*j)+k

void allocations()  // dynamic allocations of variables made as needed and defined
{
 grvy_timer_begin(__func__);
int i,j,k,temp,igot,ii;
double lol1,lol2,lol3;
l1 = malloc(1*sizeof(double));
  /* Initialize/read the file */

  igot = grvy_input_fopen("./input2.txt");

  /* Read specific variables and echo locally */

grvy_input_fread_double("x0",&x0);
grvy_input_fread_double("x",&x);
grvy_input_fread_double("yy0",&yy0);
grvy_input_fread_double("y",&y);
grvy_input_fread_double("z0",&z0);
grvy_input_fread_double("z",&z);
grvy_input_fread_double("h",&h);
grvy_input_fread_double("totaltime",&totaltime);
grvy_input_fread_double("delt",&delt);
grvy_input_fread_double("PML",&Pml);
grvy_input_fread_double("fc",&fc);
grvy_input_fread_int("nol",&nol);
rho = (double *)malloc(nol*sizeof(double));
Vp = (double *)malloc(nol*sizeof(double));
Vs = (double *)malloc(nol*sizeof(double));
mu = (double *)malloc(nol*sizeof(double));
lamda = (double *)malloc(nol*sizeof(double));
kappa = (double *)malloc(nol*sizeof(double));
rec = (int *)malloc(6*sizeof(int));
src = (int *)malloc(3*sizeof(int));
grvy_input_fread_double_vec("Du",rho,nol);
grvy_input_fread_double_vec("Vs",Vs,nol);
grvy_input_fread_double_vec("Vp",Vp,nol);
grvy_input_fread_int("order",&order);
if(order == 2) a = 0;
grvy_input_fread_int("modev",&modev);
grvy_input_fread_int_vec("rec",rec,6);
grvy_input_fread_int_vec("src",src,3);
grvy_input_fread_int("debug",&deb);

if(deb == 1)
   printf("reading data from input file and dynamically allocating memory for the global variables \n");

//printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d \n",x0,x,yy0,y,z0,z,h,totaltime,delt,Pml,fc,nol);
//printf("%lf %lf %lf %lf %lf %lf %d %d %d %d %d \n",rho[0],rho[1],Vp[0],Vp[1],Vs[0],Vs[1],order,modev,src[0],src[1],src[2]);
//printf("%d %d %d %d %d %d \n",rec[0],rec[1],rec[2],rec[3],rec[4],rec[5]);
/*printf("\n ------ Full Dump ------\n\n");
  grvy_input_fdump();
  printf(" ---- End Full Dump ----\n\n");

  printf("\n ------ Full Dump (delimited) ------\n\n");
  grvy_input_fdump_delim("# ");
  printf(" ---- End Full Dump ----\n\n");

  // Dump the whole file to a file 

  printf("\n ------ Full Dump to test.out ------\n\n");
  grvy_input_fdump_file("% ","test.out");
  printf(" -------    End Full Dump    -------\n\n");*/

  grvy_input_fclose();

x0 = x0 - Pml;
x  = x + Pml;
yy0 = yy0 - Pml;
y = y + Pml;
z = z + Pml;
hx = h;
hy=h;
hz=h;
Nx = ((x-x0)/hx)+1;
Ny = ((y-yy0)/hy)+1;
Nz = ((z-z0)/hz)+1;
Nx = Nx+4;
Ny = Ny+4;
Nz = Nz+4;
//printf("%d %d %d \n",Nx,Ny,Nz);
for(i=0;i<nol;i++)
{
mu[i]=rho[i]*Vs[i]*Vs[i];
lamda[i]= rho[i]*((Vp[i]*Vp[i])-(2*Vs[i]*Vs[i]));
kappa[i]=(lamda[i])+(2*mu[i]/3);
}
forward = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Xu = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Yu = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Zu = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 U3 = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Du = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Xv = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Yv = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Zv = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 V3 = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Dv = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Xw = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Yw = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Zw = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 W3 = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Dw = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Xxy = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Yxy = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Zxy = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Txy = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 muxy = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Xyz = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Yyz = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Zyz = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Tyz = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 muyz = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Xzx = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Yzx = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Zzx = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Tzx = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 muzx =(double *)malloc(Nx*Ny*Nz*sizeof(double));
 Xxyz = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Yxyz = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Zxyz = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Txx = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Tyy = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Tzz = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 muxyz = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 lamdaxyz = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Fu = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Fv = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 Fw = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 dx = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 dxyx = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 dzxx = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 dyzx = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 dux = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 dvx = (double *)malloc(Nx*Ny*Nz*sizeof(double));
 dwx = (double *)malloc(Nx*Ny*Nz*sizeof(double));
if(!dwx)
{
printf("Failed"); 
exit(1);
}
d0 = (-log10(R)*1.5*Vp[1]/(Pml));
for(i=0;i<Nx;i++)
{
for(j=0;j<Ny;j++)
{
for(k=0;k<Nz;k++)
{
Xu[array]=x0+((i-2)*hx);
Yu[array]=yy0+(0.5*hy)+((j-2)*hy);
Zu[array]=z0+(0.5*hz)+((k-2)*hz);
U3[array]=0.0;
Fu[array]=0.0;
lol1 = Xu[array]; lol2 = Yu[array] ; lol3 = Zu[array];
dux[array] = damping(Xu[array],Yu[array],Zu[array]);
temp = Zu[array];
for(ii=0;ii<nol-1;ii++)
{
if(temp<(z0+((ii+1)*((z-z0)/nol))))
Du[array] = rho[ii];
if(temp>(z0+((ii+1)*((z-z0)/nol))))
Du[array] = rho[ii+1];
if(temp ==(z0+((ii+1)*((z-z0)/nol))))
Du[array]= (rho[ii]+rho[ii+1])*0.5;
}

Xv[array]=x0+(0.5*hx)+((i-2)*hx);
Yv[array]=yy0+((j-2)*hy);
Zv[array]=z0+(0.5*hz)+((k-2)*hz);
V3[array]=0;
Fv[array]=0;
l1[0]=damping(Xv[array],Yv[array],Zv[array]);
dvx[array] = l1[0];
temp = Zv[array];
for(ii=0;ii<nol-1;ii++)
{
if(temp<(z0+((ii+1)*((z-z0)/nol))))
Dv[array] = rho[ii];
if(temp>(z0+((ii+1)*((z-z0)/nol))))
Dv[array] = rho[ii+1];
if(temp ==(z0+((ii+1)*((z-z0)/nol))))
Dv[array]= (rho[ii]+rho[ii+1])*0.5;
}
Xw[array]=x0+(0.5*hx)+((i-2)*hx);
Yw[array]=yy0+(0.5*hy)+((j-2)*hy);
Zw[array]=z0+((k-2)*hz);
W3[array]=0;
Fw[array]=0;
l1[0]=damping(Xw[array],Yw[array],Zw[array]);
dwx[array] = l1[0];
temp = Zw[array];
for(ii=0;ii<nol-1;ii++)
{
if(temp<(z0+((ii+1)*((z-z0)/nol))))
Dw[array] = rho[ii];
if(temp>(z0+((ii+1)*((z-z0)/nol))))
Dw[array] = rho[ii+1];
if(temp ==(z0+((ii+1)*((z-z0)/nol))))
Dw[array]= (rho[ii]+rho[ii+1])*0.5;
}

Xxy[array]=x0+((i-2)*hx);
Yxy[array]=yy0+((j-2)*hy);
Zxy[array]=z0+(0.5*hz)+((k-2)*hz);
Txy[array] = 0;

l1[0]=damping(Xxy[array],Yxy[array],Zxy[array]);
dxyx[array] = l1[0];
temp = Zxy[array];
for(ii=0;ii<nol-1;ii++)
{
if(temp<(z0+((ii+1)*((z-z0)/nol))))
muxy[array] = mu[ii];
if(temp>(z0+((ii+1)*((z-z0)/nol))))
muxy[array] = mu[ii+1];
if(temp ==(z0+((ii+1)*((z-z0)/nol))))
muxy[array]= (2.0/((1/mu[ii])+(1/mu[ii+1])));
}
Xyz[array]=x0+(0.5*hx)+((i-2)*hx);
Yyz[array]=yy0+((j-2)*hy);
Zyz[array]=z0+((k-2)*hz);
Tyz[array] = 0;
l1[0]= damping(Xyz[array],Yyz[array],Zyz[array]);
dyzx[array] = l1[0];
temp = Zyz[array];
for(ii=0;ii<nol-1;ii++)
{
if(temp<(z0+((ii+1)*((z-z0)/nol))))
muyz[array] = mu[ii];
if(temp>(z0+((ii+1)*((z-z0)/nol))))
muyz[array] = mu[ii+1];
if(temp ==(z0+((ii+1)*((z-z0)/nol))))
muyz[array]= (2.0/((1/mu[ii])+(1/mu[ii+1])));
}
Xzx[array]=x0+((i-2)*hx);
Yzx[array]=yy0+(0.5*hy)+((j-2)*hy);
Zzx[array]=z0+((k-2)*hz);
Tzx[array]=0;
l1[0]=damping(Xzx[array],Yzx[array],Zzx[array]);
dzxx[array]= l1[0];
temp = Zzx[array];
for(ii=0;ii<nol-1;ii++)
{
if(temp<(z0+((ii+1)*((z-z0)/nol))))
muzx[array] = mu[ii];
if(temp>(z0+((ii+1)*((z-z0)/nol))))
muzx[array] = mu[ii+1];
if(temp ==(z0+((ii+1)*((z-z0)/nol))))
muzx[array]= (2.0/((1/mu[ii])+(1/mu[ii+1])));
}
Xxyz[array]=x0+(0.5*hx)+((i-2)*hx);
Yxyz[array]=yy0+(0.5*hy)+((j-2)*hy);
Zxyz[array]=z0+(0.5*hz)+((k-2)*hz);
Txx[array]=0;
Tyy[array]=0;
Tzz[array]=0;

l1[0]= damping(Xxyz[array],Yxyz[array],Zxyz[array]);
dx[array] = l1[0];
temp = Zxyz[array];
for(ii=0;ii<nol-1;ii++)
{
if(temp<(z0+((ii+1)*((z-z0)/nol))))
{
muxyz[array] = mu[ii];lamdaxyz[array] = lamda[ii];}
if(temp>(z0+((ii+1)*((z-z0)/nol))))
{
muxyz[array] = mu[ii+1];lamdaxyz[array] = lamda[ii+1];}
if(temp ==(z0+((ii+1)*((z-z0)/nol))))
{
muxyz[array]= (2.0/((1/mu[ii])+(1/mu[ii+1])));lamdaxyz[array] = (2.0/((1/lamda[ii])+(1/lamda[ii+1])));}
}
}}}

  grvy_timer_end(__func__);

}

void storewave()
{
 grvy_timer_begin(__func__);
int i,j,k;
if(deb == 1)
   printf("Storng wavefields recorded\n");

for(i=rec[0];i<=rec[1];i++)  //storing velocities at receivers
{
for(j=rec[2];j<=rec[3];j++)
{
for(k=rec[4];k<=rec[5];k++)
{
forward[count2] =  U3[mat(i,j,k)];
forward[count2+1] = V3[mat(i,j,k)];
forward[count2+2]= W3[mat(i,j,k)];
count2 = count2+3;
}}}

  grvy_timer_end(__func__);
}



