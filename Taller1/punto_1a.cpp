#include <iostream>
#include <cmath>

const double Beta = 0.35;
const double Gamma = 0.08;

//En el modelo SIR s=x1, i=x2
double f1(double t, double x1, double x2){
  return -Beta*x1*x2;
}

double f2(double t, double x1, double x2){
  return Beta*x1*x2 - Gamma*x2;
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

int main(){
  double t,x1, x2, r, tf=300; double dt=0.01;
  std::cout << "t" << "\t" << "s" << "\t" << "i" << "\t" << "r" << std::endl;
  for(t=0, x1=0.999, x2=0.001; t<tf+dt/2; ){
    r=1.-(x1+x2);
    std::cout << t << "\t" << x1 << "\t" << x2 << "\t" << r << std::endl;
    runge_kutta_step(t,dt,x1, x2);
  }
  return 0;
}
