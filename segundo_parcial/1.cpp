#include "ondas.hpp"

int main(){
  LatticeBoltzmann Ondas;
  int t,tmax=6000;
  double rho0=0,Jx0=0,Jy0=0;

  Ondas.Start(rho0,Jx0,Jy0);
  for(t=0;t<tmax;t++){
    Ondas.Collision();
    Ondas.ImposeFields(t);
    Ondas.Advection();
    if(t==2000) Ondas.Print("data/1_"+std::to_string(t)+".dat");
  }
  Ondas.Print("data/1_6000.dat");

  return 0;
}
