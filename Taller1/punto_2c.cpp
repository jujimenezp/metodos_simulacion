#include <iostream>
#include <vector>
#include <cmath>

double Lambda=1;

//Ecuacion de Bessel con alpha=0. t=r, x1=dR/dr, x2=R(r)
double f1(double t, double x1, double x2);
double f2(double t, double x1, double x2);

//Métod de Runge-Kutta 4to oren
void runge_kutta_step(double &t0, double dt, double &x10, double &x20);

//funcion que calcula f(lambda)=R(r=x2, lambda)
double f(double t0, double dt, double x10, double x20, double tf);

//Metodo de bisección
double biseccion(double Lambda_a,double Lambda_b,double res,double tol,int maxiter,
                 double t0, double dt, double x10, double x20, double tf);

int main(){
  double t=0.01,x1=0, x2=1, tf=1.; double dt=0.01, Lambda_f=0;
  std::vector<double> Lambdas;
  for(double ii=0; ii<=14.5;ii+=0.5){
    Lambda_f=biseccion(ii,ii+0.5,0.,1e-16,50,t,dt,x1,x2,tf);
    if(Lambda_f!=0){
      Lambdas.push_back(Lambda_f);
      std::clog << "Lambda: " << Lambda_f << std::endl;
    }
  }

  Lambda=Lambdas[4];
  for(; t<tf+dt/2; t+=dt){
    runge_kutta_step(t,dt,x1,x2);
    std::cout << t << "\t" << x2 << std::endl;
  }
  return 0;
}

double f1(double t, double x1, double x2){
  return -x1/t-pow(Lambda,2)*x2;
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
}

double f(double t0, double dt, double x10, double x20, double tf){
  double t, x1, x2;
  for(t=t0, x1=x10, x2=x20; t<tf+dt/2; t+=dt){
    runge_kutta_step(t,dt,x1,x2);
  }
  return x2;
}

double biseccion(double Lambda_a,double Lambda_b,double res,double tol,int maxiter,
                 double t0, double dt, double x10, double x20, double tf){
  double fa, fb, fc, Lambda_c=0;
  Lambda=Lambda_a; fa=f(t0, dt, x10, x20, tf);
  Lambda=Lambda_b; fb=f(t0, dt, x10, x20, tf);
  if(std::signbit(fa)!=std::signbit(fb)){
    for(int ii=0;ii<maxiter;ii++){
      Lambda_c=(Lambda_b+Lambda_a)/2;
      Lambda=Lambda_c; fc=f(t0, dt, x10, x20, tf);
      if(std::signbit(fc)==std::signbit(fa)){
        Lambda_a=Lambda_c;
        fa=fc;
      }
      else{
        Lambda_b=Lambda_c;
        fb=fc;
      }
      if(fabs(Lambda_b-Lambda_a)/2<=tol || fabs(fc)<=res) break;
    }
  }
return Lambda_c;
}
