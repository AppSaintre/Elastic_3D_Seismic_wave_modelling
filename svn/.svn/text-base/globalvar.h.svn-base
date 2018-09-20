#define array (Ny*Nz*i)+(Nz*j)+k
#define mat(i,j,k) (Ny*Nz*i)+(Nz*j)+k     // declaring all the global variables needed

int shotloc[6] = {23,62,102,142,182,219};
int k1 = 4001;
int k2 = 600;
int Nx,Ny,Nz,ml,count2,nol,order,modev,*rec,*src,deb;
char *griddens,*gridVp,*gridVs;
double x0,x,y,z0,z,h,hx,hy,hz,totaltime,delt,Pml,yy0,d0;
double *rho,*Vp,*Vs,*mu,*lamda,*kappa;
double *l1,*forward,*Xu,*Yu,*Zu,*Du,*U3,*Xv,*Yv,*Zv,*Dv,*V3,*Xw,*Yw,*Zw,*Dw,*W3,*Fu,*Fv,*Fw;
double  *dx,*dxyx,*dzxx,*dyzx,*dux,*dvx,*dwx;
double *Xxy,*Yxy,*Zxy,*Txy,*muxy,*Xyz,*Yyz,*Zyz,*Tyz,*muyz,*Xzx,*Yzx,*Zzx,*Tzx,*muzx,*Xxyz,*Yxyz,*Zxyz,*lamdaxyz,*muxyz,*Txx,*Tyy,*Tzz;
double a = -1.0/24.0;
double b = 9.0/8.0;
double mid = 200.0;
double R = 0.000001;
double fc;


