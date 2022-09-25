#include <iostream>
#include <cmath>

double Beta = 0.35;
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

void S_inf(double &t0, double dt, double &x10, double &x20, double tf){
  for( ;t0<tf+dt/2; ){
    runge_kutta_step(t0,dt,x10, x20);
  }
}

int main(){
  double t=0,x1=0.999, x2=0.001, tf=2000; double dt=0.01;
  double R0=Beta/Gamma;
  std::cout << "R_0" << "\t" << "S_inf" << std::endl;
  for(Beta=0.001; Beta<1.; Beta+=0.01){
    t=0; x1=0.999; x2=0.001;
    S_inf(t, dt, x1, x2, tf);
    R0=Beta/Gamma;
    std::cout << R0 << "\t" << x1 << std::endl;
  }
  return 0;
}
