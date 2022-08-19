#include <iostream>
#include <cmath>

double f(double t, double x){
  return sin(x);
}

void runge_kutta_step(double &t, double dt, double &x){
  double dx1,dx2,dx3,dx4;
  dx1 = dt*f(t,x);
  dx2 = dt*f(t+dt/2,x + dx1/2);
  dx3 = dt*f(t+dt/2,x + dx2/2);
  dx4 = dt*f(t+dt,x + dx3);
  x+=(dx1 + 2*(dx2+dx3)+dx4)/6; t+=dt;
}

int main(){
  double t,x, tf=50; double dt=0.01;

  for(t=0, x=1; t<tf+dt/2; ){
    std::cout << t << "\t" << x << std::endl;
    runge_kutta_step(t,dt,x);
  }
  return 0;
}
