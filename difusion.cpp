#include <cmath>
#include "Random64.h"

const int Lx=16; //1024
const double p=0.5;

const int Q=2;

//Clase LatticeGas
class LatticeGas{
private:
  int V[Q]; //V[i] i=0 derecha, Q=1 izquierda
  int n[Lx][Q], nnew[Lx][Q]; //n[ix][i]
public:
  LatticeGas(void);
  void Borrar();
  void Inicie(int N,double mu,double sigma,Crandom & ran64);
  void Show();
  void Colisione(Crandom & ran64);
  void Adveccione();
  double rho(int ix){return n[ix][0]+n[ix][1];}
  double Varianza();
};

LatticeGas::LatticeGas(){
  V[0]=1; V[1]=-1;
}

void LatticeGas::Borrar(){
  for(int ix=0; ix<Lx;ix++){
    for(int i=0;i<Q;i++){
      n[ix][i]=nnew[ix][i]=0;
    }
  }
}

void LatticeGas::Inicie(int N,double mu,double sigma,Crandom & ran64){
  int ix,i;
  while(N>0){
    ix=(int) ran64.gauss(mu,sigma);//Escoger una celda al azar;
    if(ix<0) ix=0;
    if(ix>Lx-1) ix=Lx-1;
    i=(int) Q*ran64.r();
    if(n[ix][i]==0){
      n[ix][i]=1; N--;
    }
  }
}


void LatticeGas::Show(){
  for (int i=0; i < Q; i++) {
    for (int ix = 0; ix < Lx; ix++) {
      std::cout << n[ix][i]<< std::endl;
    }
  }
}

void LatticeGas::Colisione(Crandom & ran64){
  for(int ix=0;ix<Lx;ix++){
    if(ran64.r() >p){
      nnew[ix][0]=n[ix][1]; nnew[ix][1]=n[ix][0];
    }
    else{
      nnew[ix][0]=n[ix][0]; nnew[ix][1]=n[ix][1];
    }
  }
}

void LatticeGas::Adveccione(){
  for(int ix=0; ix <Lx;ix++){
    for(int i=0;i <Q;i++){
      n[ix+V[i]][i]=nnew[ix][i];
    }
  }
}

double LatticeGas::Varianza(){
  int ix; double N,Xprom,Sigma2;
  for(N=0, ix=0; ix <Lx;ix++){
    N+=rho(ix);
  }
  for(Xprom=0,ix=0;ix <Lx;ix++) Xprom+=ix*rho(ix);
  Xprom/=N;
  for(Sigma2=0,ix=0;ix <Lx;ix++) Sigma2+=pow(ix-Xprom,2.0)*rho(ix);
  Sigma2/=(N-1);
  return Sigma2;
}

int main(){
  LatticeGas Difusion;
  Crandom ran64(1);
  int N=50; double mu=Lx/2, sigma=Lx/8;
  int t, tmax=400;

  Difusion.Borrar();
  Difusion.Inicie(N, mu, sigma, ran64);
  for(t=0;t <tmax;t++){
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
  }
  Difusion.Show();

  return 0;
}
