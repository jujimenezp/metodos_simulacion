#include <iostream>
#include <cmath>

const double omega = 3.0;

double f1(double t, double x1, double x2){
  return -omega*omega*x2;
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

int main(){
  double t,x1, x2, tf=50; double dt=0.01;

  for(t=0, x=1, x2=0; t<tf+dt/2; ){
    std::cout << t << "\t" << x << std::endl;
    runge_kutta_step(t,dt,x1, x2);
  }
  return 0;
}
