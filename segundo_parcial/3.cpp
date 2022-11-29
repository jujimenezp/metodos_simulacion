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
  }
  auto Amax= std::max_element(std::begin(pmax), std::end(pmax));
  auto Amin= std::min_element(std::begin(pmax), std::end(pmax));
  std::clog << "Punto máximo en envolvente superior: " <<*Amax << std::endl
            << "Punto mínimo en envolvente superior: " << *Amin << std::endl;

  double SWR = *Amax/(*Amin);
  double C_abs = 4*SWR/pow(SWR+1,2.);
  std::clog << "Standing Wave Ratio: " << SWR << std::endl
            << "Coeficiente de absorcion: " << C_abs << std::endl;

  Ondas.Print("data/2_final.dat");

  return 0;
}
