#include <cmath>

const int Lx=1024;
const double p=0.5;

const int Q=2;

//Clase LatticeGas
class LatticeGas{
private:
  int V[Q]; //V[i] i=0 derecha, Q=1 izquierda
  double f[Lx][Q], fnew[Lx][Q]; //n[ix][i]
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
    double rho=N*1.0/(sigma*sqrt(2*M_PI))*exp(-0.5*pow((ix-mu)/sigma,2.0));
    for(int i=0;i <Q;i++){
      f[ix][i]=rho/Q;
    }
  }
}

double LatticeGas::rho(int ix){
  double suma; int i;
  for (suma=0, i = 0; i < Q; i++) {
    suma+=f[ix][i];
  }
  return suma;
}

void LatticeGas::Colisione(){
  int ix,i,j;
  for(ix=0;ix<Lx;ix++){
    for(i=0;i <Q;i++){
      j=(i+1)%Q;
      fnew[ix][i]=f[ix][i]+(1-p)*(f[ix][i]);
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
