
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

void arenstorf(double* yt, double* y);
void step(double* r4, double* r5, double* error, double& emax1, double& emax2, double& errormax);


int main(void){
  
  ofstream out("solution");
  
  double t=0.0, T=17.065216560157, t1=17.1;
  double dt=0.001, tol=1e-5;
  double k1[4], k2[4], k3[4], k4[4], k5[4], k6[4], k7[4];
  double y0[4];
  double ytemp[4];
  double r4[4];
  double r5[4];
  double error[4], emax1, emax2, errormax;
  
  y0[0]=0.994;
  y0[1]=0.0;
  y0[2]=0.0;
  y0[3]=-2.00158510637908;
  
  double a21=1.0/5.0;
  double a31=3.0/40.0, a32=9.0/40.0;
  double a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0;
  double a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0;
  double a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0;
  double a71=35.0/384.0, a73=500.0/1113.0, a74=125.0/192.0, a75=-2187.0/6784.0, a76=11.0/84.0;
  
  double b11=a71, b13=a73, b14=a74, b15=a75, b16=a76;
  double b21=5179.0/57600.0, b23=7571.0/16695.0, b24=393.0/640.0, b25=-92097.0/339200.0, b26=187.0/2100.0, b27=1.0/40.0;
  
  
  out << t << "\t" << y0[0] << "\t" << y0[1] << "\t" << y0[2] << "\t" << y0[3] << "\t" << dt << endl;

while(t<t1){
   
   //K-Berechnung
   arenstorf(k1, y0);
   
   for(int j=0; j<4; j++){
    ytemp[j]=y0[j]+a21*k1[j]*dt; 
   }
   
   arenstorf(k2, ytemp);
   
   for(int j=0; j<4; j++){
    ytemp[j]=y0[j]+dt*(a31*k1[j]+a32*k2[j]); 
   }
   
   arenstorf(k3, ytemp);
   
   for(int j=0; j<4; j++){
    ytemp[j]=y0[j]+dt*(a41*k1[j]+a42*k2[j]+a43*k3[j]); 
   }
   
   arenstorf(k4, ytemp);
   
   for(int j=0; j<4; j++){
    ytemp[j]=y0[j]+dt*(a51*k1[j]+a52*k2[j]+a53*k3[j]+a54*k4[j]); 
   }
   
   arenstorf(k5, ytemp);
   
   for(int j=0; j<4; j++){
    ytemp[j]=y0[j]+dt*(a61*k1[j]+a62*k2[j]+a63*k3[j]+a64*k4[j]+a65*k5[j]); 
   }
   
   arenstorf(k6, ytemp);
   
   for(int j=0; j<4; j++){
    ytemp[j]=y0[j]+dt*(a71*k1[j]+a73*k3[j]+a74*k4[j]+a75*k5[j]+a76*k6[j]); 
   }
   
   arenstorf(k7, ytemp);
   
   //RK 5. Ordnung
   for(int j=0; j<4; j++){
    r5[j]=y0[j]+dt*(b11*k1[j]+b13*k3[j]+b14*k4[j]+b15*k5[j]+b16*k6[j]); 
   }
   
   arenstorf(k7, ytemp);
   
   //RK 5. Ordnung
   for(int j=0; j<4; j++){
    r5[j]=y0[j]+dt*(b11*k1[j]+b13*k3[j]+b14*k4[j]+b15*k5[j]+b16*k6[j]); 
   }
   
   
   //RK 4. Ordnung
   for(int j=0; j<4; j++){
    r4[j]=y0[j]+dt*(b21*k1[j]+b23*k3[j]+b24*k4[j]+b25*k5[j]+b26*k6[j]+b27*k7[j]); 
   }
   
   step(r4, r5, error, emax1, emax2, errormax);
   
   dt=dt*pow(tol/errormax, 0.2);
   
   for(int i=0; i<4; i++){
    y0[i]=r4[i]; 
   }
   
   t += dt;
   out << t << "\t" << y0[0] << "\t" << y0[1] << "\t" << y0[2] << "\t" << y0[3] << "\t" << dt << endl;
   
  }
 
 
 
 out.close(); 
 return(0); 
}

void arenstorf(double* yt, double* x){
  const double mu1 = 0.012277471;
  const double mu2 = 1-mu1;
 
 double r = sqrt(pow(x[0]+mu1,2.0)+pow(x[1],2.0));
 double s = sqrt(pow(x[0]-1+mu1,2.0)+pow(x[1],2.0));
 
 yt[0]=x[2];
 yt[1]=x[3];
 yt[2]=x[0]+2.0*x[3]-(mu2*(x[0]+mu1)/pow(r,3.0))-(mu1*(x[0]-1+mu1)/pow(s,3.0));
 yt[3]=x[1]-2.0*x[2]-(mu2*x[1]/pow(r,3.0))-(mu1*x[1]/pow(s,3.0));
}

//step-size step
void step(double* r4, double* r5, double* error, double& emax1, double& emax2, double& errormax){
  for(int k=0; k<4; k++){
   error[k] = abs(r4[k]-r5[k]);
  emax1= max(error[0],error[1]);
   emax2= max(error[2],error[3]);
   errormax= max(emax1,emax2);
  }
}
