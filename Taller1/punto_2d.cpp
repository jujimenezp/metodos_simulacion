#include <iostream>
#include <cmath>

const double Lambda=5.67129;
double f(double alpha, double x, double t){
  return cos(alpha*t - Lambda*x*sin(t));
}

double Simpson_Integral(double alpha, double x, double a, double b, int n){
  double t,h,suma; int i;
  n*=2; h=(b-a)/n;
  for(suma=0, i=0;i<=n;i++){
    t=a+i*h;
    if(i==0 || i==n)
      suma+=f(alpha,x,t);
    else if(i%2==0)
      suma+=2*f(alpha,x,t);
    else
      suma+=4*f(alpha,x,t);
    }
  return suma*h/3;
}

double Bessel(double alpha, double x){
  double a=0,b=M_PI; int n=100;
  return Simpson_Integral(alpha,x,a,b,n)/M_PI;
}

int main(){
  double alpha=0,x;
  for(x=0;x<=1;x+=0.05){
    std::cout << x << "\t" << Bessel(alpha, x) << std::endl;
  }
  return 0;
}
