#include "ondas.hpp"

int main(){
  LatticeBoltzmann Ondas;
  int t,tmax=80000;
  double rho0=0,Jx0=0,Jy0=0;
  std::vector<double> pmin(Lx-2,0.), pmax(Lx-2,0.);

  Ondas.Start(rho0,Jx0,Jy0);
  for(t=0;t<tmax;t++){
    Ondas.Collision();
    Ondas.ImposeFields(t);
    Ondas.Advection();
  }
  for(t=tmax;t < tmax+1.1*lambda/C;t++){
    Ondas.Collision();
    Ondas.ImposeFields(t);
    Ondas.Advection();
    Ondas.envolventes(pmin, pmax);
    if(t%20==0) Ondas.Print("data/2_"+std::to_string(t)+".dat");
  }

  std::ofstream envo("data/2_envolventes.dat");
  for(int ix=1; ix < Lx-1; ix++){
    envo << ix <<"\t"<< pmin[ix-1] <<"\t"<< pmax[ix-1] << std::endl;
  }
  envo.close();
  Ondas.Print("data/2_final.dat");

  return 0;
}
