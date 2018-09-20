#include<stdio.h>
#include<stdlib.h>
#include<math.h>


double damping(double Xu,double Yu,double Zu)  // absorbing boundaries 
{
double du;
du=1; //0;
du=1;//0;
du=1;//0;
if(Xu < (x0+Pml))
du = exp(-(0.015*0.015*((x0+Pml-Xu)/2)*((x0+Pml-Xu)/2)));   
if(Xu > (x-Pml))
du =  exp(-(0.015*0.015*((Xu-(x-Pml))/2)*((Xu-(x-Pml))/2)));
if(Yu < (yy0+Pml))
du =  exp(-(0.015*0.015*((yy0+Pml-Yu)/2)*((yy0+Pml-Yu)/2)));//d0*((y0+Pml-Yu)/Pml)*((y0+Pml-Yu)/Pml);
if(Yu > (y-Pml))
du = exp(-(0.015*0.015*((Yu-(y-Pml))/2)*((Yu-(y-Pml))/2)));//d0*((Yu-(y-Pml))/Pml)*((Yu-(y-Pml))/Pml);
if(Zu > (z-Pml))
du = exp(-(0.015*0.015*((Zu-(z-Pml))/2)*((Zu-(z-Pml))/2))); 
return du;
}


void stressupdate()   // function that updates stresses
{
if(deb == 1)
   printf("Updating stresses \n");
 grvy_timer_begin(__func__);
int i,j,k;
double value_x,value_y,value_z,lol1;
#pragma omp parallel for private(value_x,value_y,value_z,lol1,j,k)
for(i=2;i<Nx-2;i++)
{
for(j=2;j<Ny-2;j++)
{
for(k=2;k<Nz-2;k++)
{

value_x = (a*(U3[mat((i+2),j,k)]-U3[mat((i-1),j,k)]))+(b*(U3[mat((i+1),j,k)]-U3[mat(i,j,k)]));
value_y = (a*(V3[mat(i,(j+2),k)]-V3[mat(i,(j-1),k)]))+(b*(V3[mat(i,(j+1),k)]-V3[mat(i,j,k)]));
value_z = (a*(W3[mat(i,j,(k+2))]-W3[mat(i,j,(k-1))]))+(b*(W3[mat(i,j,(k+1))]-W3[mat(i,j,k)]));

lol1 = dx[array];

Txx[array] = (Txx[array]+(((delt*(lamdaxyz[array]+(2*muxyz[array]))*(value_x)/(h))))+(((delt*(lamdaxyz[array])*(value_y)/(h))))+(((delt*(lamdaxyz[array])*(value_z)/(h)))))*lol1;

Tyy[array] = (Tyy[array]+(((delt*(lamdaxyz[array])*(value_x)/(h))))+(((delt*(lamdaxyz[array]+(2*muxyz[array]))*(value_y)/(h))))+(((delt*(lamdaxyz[array])*(value_z)/(h)))))*lol1;
Tzz[array] = (Tzz[array]+(((delt*(lamdaxyz[array])*(value_x)/(h))))+(((delt*(lamdaxyz[array])*(value_y)/(h))))+(((delt*(lamdaxyz[array]+(2*muxyz[array]))*(value_z)/(h)))))*lol1;

value_x = (a*(U3[mat(i,(j+1),k)]-U3[mat(i,(j-2),k)]))+(b*(U3[mat(i,j,k)]-U3[mat(i,(j-1),k)]));
value_y = (a*(V3[mat((i+1),j,k)]-V3[mat((i-2),j,k)]))+(b*(V3[mat(i,j,k)]-V3[mat((i-1),j,k)]));
lol1 = dxyx[array];
Txy[array] = (Txy[array]+ (((delt*(muxy[array])*(value_y)/(h))))+(((delt*(muxy[array])*(value_x)/(h)))))*lol1;

value_x = (a*(U3[mat(i,j,(k+1))]-U3[mat(i,j,(k-2))]))+(b*(U3[mat(i,j,k)]-U3[mat(i,j,(k-1))]));
value_z = (a*(W3[mat((i+1),j,k)]-W3[mat((i-2),j,k)]))+(b*(W3[mat(i,j,k)]-W3[mat((i-1),j,k)]));
lol1 = dzxx[array];

Tzx[array] = (Tzx[array]+ (((delt*(muzx[array])*(value_z)/(h)))) +(((delt*(muzx[array])*(value_x)/(h)))))*lol1 ;
if(Zzx[array]==0)
{
Tzx[array] = 0;
}
value_y = (a*(V3[mat(i,j,(k+1))]-V3[mat(i,j,(k-2))]))+(b*(V3[mat(i,j,k)]-V3[mat(i,j,(k-1))]));
value_z = (a*(W3[mat(i,(j+1),k)]-W3[mat(i,(j-2),k)]))+(b*(W3[mat(i,j,k)]-W3[mat(i,(j-1),k)]));
lol1 = dyzx[array];
Tyz[array] = (Tyz[array]+ (((delt*(muyz[array])*(value_z)/(h)))) + (((delt*(muyz[array])*(value_y)/(h)))))*lol1 ;
if(Zyz[array]==0)
{
Tyz[array] = 0;
}


}
Tzz[mat(i,j,1)] = -Tzz[mat(i,j,2)];
Tzz[mat(i,j,0)] = -Tzz[mat(i,j,3)];
Tzx[mat(i,j,1)]=-Tzx[mat(i,j,3)];

Tyz[mat(i,j,1)]=-Tyz[mat(i,j,3)];

}
}

  grvy_timer_end(__func__);
}


void velocityupdate()// function that updates velocities
{
if(deb == 1)
   printf("Updating velocities \n");
 grvy_timer_begin(__func__);
int i,j,k;
double value_x,value_y,value_z,lol1;
#pragma omp parallel for private(value_x,value_y,value_z,lol1,j,k)
for(i=2;i<Nx-2;i++)
{
for(j=2;j<Ny-2;j++)
{
for(k=2;k<Nz-2;k++)
{

value_x = (a*(Txx[mat((i+1),j,k)]-Txx[mat((i-2),j,k)]))+(b*(Txx[mat(i,j,k)]-Txx[mat((i-1),j,k)]));
value_y = (a*(Txy[mat(i,(j+2),k)]-Txy[mat(i,(j-1),k)]))+(b*(Txy[mat(i,(j+1),k)]-Txy[mat(i,j,k)]));
value_z = (a*(Tzx[mat(i,j,(k+2))]-Tzx[mat(i,j,(k-1))]))+(b*(Tzx[mat(i,j,(k+1))]-Tzx[mat(i,j,k)]));
lol1 = dux[array];

U3[array] = U3[array]+(delt*(value_x+value_y+value_z)/(Du[array]*h))+(delt*Fu[array]/Du[array]);
U3[array] = U3[array]*lol1;

if(Xu[array]==x0 || Xu[array]==x || Yu[array]==yy0 || Yu[array]==y || Zu[array]==z )
{
U3[array] = 0;
}

value_x = (a*(Txy[mat((i+2),j,k)]-Txy[mat((i-1),j,k)]))+(b*(Txy[mat((i+1),j,k)]-Txy[mat(i,j,k)]));
value_y = (a*(Tyy[mat(i,(j+1),k)]-Tyy[mat(i,(j-2),k)]))+(b*(Tyy[mat(i,j,k)]-Tyy[mat(i,(j-1),k)]));
value_z = (a*(Tyz[mat(i,j,(k+2))]-Tyz[mat(i,j,(k-1))]))+(b*(Tyz[mat(i,j,(k+1))]-Tyz[mat(i,j,k)]));
lol1 = dvx[array];

V3[array] = V3[array]+(delt*(value_x+value_y+value_z)/(Dv[array]*h))+(delt*Fv[array]/Dv[array]);
V3[array] = V3[array]*lol1;

if(Xv[array]==x0 || Xv[array]==x || Yv[array]==yy0 || Yv[array]==y || Zv[array]==z)
{
V3[array] = 0;
}

value_x = (a*(Tzx[mat((i+2),j,k)]-Tzx[mat((i-1),j,k)]))+(b*(Tzx[mat((i+1),j,k)]-Tzx[mat(i,j,k)]));
value_y = (a*(Tyz[mat(i,(j+2),k)]-Tyz[mat(i,(j-1),k)]))+(b*(Tyz[mat(i,(j+1),k)]-Tyz[mat(i,j,k)]));
value_z = (a*(Tzz[mat(i,j,(k+1))]-Tzz[mat(i,j,(k-2))]))+(b*(Tzz[mat(i,j,k)]-Tzz[mat(i,j,(k-1))]));
lol1 = dwx[array];
W3[array] = W3[array]+(delt*(value_x+value_y+value_z)/(Dw[array]*h))+(delt*Fw[array]/Dw[array]);
W3[array] = W3[array]*lol1;


if(Xw[array]==x0 || Xw[array]==x || Yw[array]==yy0 || Yw[array]==y || Zw[array]==z )
{
W3[array] = 0;
}

}

}}

  grvy_timer_end(__func__);
}

void boundarycns() // function to apply free surface conditions
{
if(deb == 1)
   printf("Applying free surface conditions \n");
 grvy_timer_begin(__func__);
int i,j,k;
double temp;
#pragma omp parallel for private(j,k)
for(i=2;i<Nx-2;i++)
{
for(j=2;j<Ny-2;j++)
{
for(k=2;k<Nz-2;k++)
{
U3[mat(i,j,1)] = U3[mat(i,j,2)] + W3[mat(i,j,2)] - W3[mat((i-1),j,2)];
if(Xu[array]==x0 || Xu[array]==x || Yu[array]==yy0 || Yu[array]==y || Zu[array]==z)
{
U3[array] = 0;
}

V3[mat(i,j,1)] = V3[mat(i,j,2)] + W3[mat(i,j,2)] - W3[mat(i,(j-1),2)];
if(Xv[array]==x0 || Xv[array]==x || Yv[array]==yy0 || Yv[array]==y || Zv[array]==z)
{
V3[array] = 0;
}

temp = lamdaxyz[array]/(lamdaxyz[array]+(2.0*muxyz[array]));
W3[mat(i,j,1)] = W3[mat(i,j,3)] + (temp*(U3[mat((i+1),j,1)]-U3[mat(i,j,1)]+V3[mat(i,(j+1),1)]-V3[mat(i,j,1)]+U3[mat((i+1),j,2)]-U3[mat(i,j,2)]+V3[mat(i,(j+1),2)]-V3[mat(i,j,2)]));
if(Xw[array]==x0 || Xw[array]==x || Yw[array]==yy0 || Yw[array]==y || Zw[array]==z)
{
W3[array] = 0;
}
if(U3[array] > 1e-8 || V3[array] > 1e-8 || W3[array] > 1e-8)
{
printf("outof u limits");
exit(-1);
}
}}}

  grvy_timer_end(__func__);
}

void forcing(double t)  // point source function 
{
if(deb == 1)
   printf("Injecting the source at source location \n");
 grvy_timer_begin(__func__);
double ts,tp, gamma,temp;
tp = 1/fc;
ts = 1.4*tp;
if(t<=2*ts)
{
gamma = (3.14159*fc*(t-ts));
temp = ((gamma*gamma)-0.5)*exp(-gamma*gamma);
Fw[mat(src[0],src[1],src[2])] = temp;

}
else
{
Fw[mat(src[0],src[1],src[2])] = 0; 
 }

  grvy_timer_end(__func__);
}


