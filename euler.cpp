#include <iostream>
#include <cmath>

double f(double t, double x){
  return sin(x);
}

void euler_step(double &t, double dt, double &x){
  double dx;
  dx = dt*f(t,x);
  x+=dx; t+=dt;
}

int main(){
  double t,x; double dt=0.01;

  for(t=0, x=1; t<2+dt/2; ){
    std::cout<<t<<"\t"<<x<<std::endl;
    euler_step(t,dt,x);
  }
  return 0;
}
