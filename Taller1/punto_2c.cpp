#include <iostream>
#include <cmath>

const double a=1;
double Lambda=1;

//Ecuacion de Bessel con alpha=0. t=r, x1=dR/dr, x2=R(r)
double f1(double t, double x1, double x2){
  return -x1/t-Lambda*x2;
}

double f2(double t, double x1, double x2){
  return x1;
}

void runge_kutta_step(double &t0, double dt, double &x10, double &x20){
  double dx11,dx21,dx31,dx41;
  double dx12,dx22,dx32,dx42;
  dx11 = dt*f1(t0,x10, x20);                     dx12 = dt*f2(t0,x10,x20);
  dx21 = dt*f1(t0+dt/2,x10+dx11/2,x20+dx12/2);   dx22 = dt*f2(t0+dt/2,x10+dx11/2,x20+dx12/2);
  dx31 = dt*f1(t0+dt/2,x10+dx21/2, x20+dx22/2);  dx32 = dt*f2(t0+dt/2,x10+dx21/2, x20+dx22/2);
  dx41 = dt*f1(t0+dt,x10+dx31, x20+dx32);        dx42 = dt*f2(t0+dt,x10 + dx31, x20+dx32);
  x10+=(dx11 + 2*(dx21+dx31)+dx41)/6;            x20+=(dx12 + 2*(dx22+dx32)+dx42)/6;
  t0+=dt;
}

double f(double t0, double dt, double x10, double x20, double tf){
  double t, x1, x2;
  for(t=t0, x1=x10, x2=x20; t<tf+dt/2; ){
    runge_kutta_step(t,dt,x1,x2);
  }
  return x2;
}

int main(){
  double t=0.01,x1=0, x2=1, tf=1., f_lambda; double dt=0.005;
  std::cout << "Lambda" << "\t" << "f(Lambda)" << std::endl;
  for(Lambda=0.1; Lambda<=15.; Lambda+=0.1){
    f_lambda=f(t,dt,x1,x2,tf);
    std::cout << Lambda << "\t" << f_lambda << std::endl;
  }
  return 0;
}
