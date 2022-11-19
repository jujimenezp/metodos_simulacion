#include <cmath>
#include <iostream>

const int Lx=1024;
const double p=0.5;

const int Q=2;

//Clase LatticeGas
class LatticeGas{
private:
  int V[Q]; //V[i] i=0 derecha, Q=1 izquierda
  double f[Lx][Q], fnew[Lx][Q]; //f[ix][i]
public:
  LatticeGas(void);
  void Inicie(int N,double mu,double sigma);
  void Colisione();
  void Adveccione();
  double rho(int ix);
  double Varianza();
  void GrafiqueRho();
};

LatticeGas::LatticeGas(){
  V[0]=1; V[1]=-1;
}

void LatticeGas::Inicie(int N,double mu,double sigma){
  for (int ix=0;ix <Lx;ix++){
    double rho=N/(sigma*sqrt(2*M_PI))*exp(-0.5*pow((ix-mu)/sigma,2.0));
    for(int i=0;i <Q;i++){
      f[ix][i]=rho/Q;
    }
  }
}

double LatticeGas::rho(int ix){
  double suma=0; int i;
  for (i = 0; i < Q; i++) {
    suma+=f[ix][i];
  }
  return suma;
}

void LatticeGas::GrafiqueRho(){
  for(int ix=0;ix < Lx;ix++){
    std::cout << ix << "\t"<<rho(ix) <<std::endl;
  }
}

void LatticeGas::Colisione(){
  int ix,i,j;
  for(ix=0;ix<Lx;ix++){
    for(i=0;i <Q;i++){
      j=(i+1)%Q;
      fnew[ix][i]=f[ix][i]+(1-p)*(f[ix][j]-f[ix][i]);
    }
  }
}

void LatticeGas::Adveccione(){
  for(int ix=0; ix <Lx;ix++){
    for(int i=0;i <Q;i++){
      f[(ix+V[i]+Lx)%Lx][i]=fnew[ix][i];
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
  Sigma2/=N;
  return Sigma2;
}

int main(){
  LatticeGas Difusion;
  int N=50; double mu=Lx/2, sigma=Lx/8;
  int t, tmax=400;

  Difusion.Inicie(N, mu, sigma);
  for(t=0;t <tmax;t++){
    std::cout <<t <<"\t"<< Difusion.Varianza()<< std::endl;
    Difusion.Colisione();
    Difusion.Adveccione();
  }

  return 0;
}
