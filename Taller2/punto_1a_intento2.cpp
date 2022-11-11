#include <cmath>
#include <fstream>
#include "../Random64.h"

const int Lx=256;
const int Ly=256;
const double p0=0.25;
const double p=0.25;

const int Q=4;

//Clase LatticeGas
class LatticeGas{
private:
  int V[Q]; //V[i] i=0 derecha, i=1 izquierda
  int n[Lx][Ly][Q], nnew[Lx][Ly][Q]; //n[ix][i]
public:
  LatticeGas(void);
  void Borrar();
  void Inicie(int N,double mux, double muy, double sigmax,double sigmay,Crandom & ran64);
  void Show();
  void GrafiqueRho();
  void Colisione(Crandom & ran64);
  void Adveccione();
  double rho(int ix, int iy){return n[ix][iy][0]+n[ix][iy][1]+n[ix][iy][2]+n[ix][iy][3];}
  double Varianza();
  void GrafiqueDistribucion(std::string filename);
};

LatticeGas::LatticeGas(){
  V[0]=V[1]=1; V[2]=V[3]=-1;
}

void LatticeGas::Borrar(){
  for(int ix=0; ix<Lx;ix++){
    for(int iy=0; iy<Ly;iy++){
      for(int i=0;i<Q;i++){
        n[ix][iy][i]=nnew[ix][iy][i]=0;
      }
    }
  }
}

void LatticeGas::Inicie(int N,double mux, double muy, double sigmax, double sigmay,Crandom & ran64){
  int ix,iy,i;
  while(N>0){
    ix=(int) ran64.gauss(mux,sigmax); //Escoger una celda al azar;
    iy=(int) ran64.gauss(muy,sigmay); //Escoger una celda al azar;
    if(ix<0) {ix=0;} if(ix>Lx-1) {ix=Lx-1;}
    if(iy<0) {iy=0;} if(iy>Ly-1) {iy=Ly-1;}
    i=(int) Q*ran64.r();
    if(n[ix][iy][i]==0){
      n[ix][iy][i]=1; N--;
    }
  }
}



void LatticeGas::Show(){
  for (int i=0; i < Q; i++) {
    for (int ix = 0; ix < Lx; ix++) {
      for(int iy=0; iy < Ly; iy++){
      std::cout << n[ix][iy][i];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void LatticeGas::Colisione(Crandom & ran64){
  double ran=ran64.r();
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      //No rotar
      if(ran64.r() <p0){
        nnew[ix][iy][0]=n[ix][iy][0];  nnew[ix][iy][1]=n[ix][iy][1];
        nnew[ix][iy][2]=n[ix][iy][2];  nnew[ix][iy][3]=n[ix][iy][3];
      }
      //Rotar a la izquierda
      else if(ran64.r() >=p0 && ran64.r() <p+p0){
        nnew[ix][iy][0]=n[ix][iy][3];  nnew[ix][iy][1]=n[ix][iy][0];
        nnew[ix][iy][2]=n[ix][iy][1];  nnew[ix][iy][3]=n[ix][iy][2];
      }
      //Rotar a la derecha
      else if(ran64.r() >=p0+p && ran64.r() < 2*p+p0){
        nnew[ix][iy][0]=n[ix][iy][1];  nnew[ix][iy][1]=n[ix][iy][2];
        nnew[ix][iy][2]=n[ix][iy][3];  nnew[ix][iy][3]=n[ix][iy][0];
      }
      //rotar 180
      else{
        nnew[ix][iy][0]=n[ix][iy][2];  nnew[ix][iy][1]=n[ix][iy][3];
        nnew[ix][iy][2]=n[ix][iy][0];  nnew[ix][iy][3]=n[ix][iy][1];
      }
    }
  }
}

void LatticeGas::Adveccione(){
  for(int ix=0; ix <Lx;ix++){
    for(int iy=0; iy <Ly; iy++){
      for(int i=0;i <Q;i++){
        if(i%2==0){ n[(ix+V[i]+Lx)%Lx][iy][i]=nnew[ix][iy][i];}
        else{ n[ix][(iy+V[i]+Ly)%Ly][i]=nnew[ix][iy][i];}
        //n[(ix+V[i]+Lx)%Lx][i]=nnew[ix][i];
      }
    }
  }
}

double LatticeGas::Varianza(){
  int ix,iy; double N,Xprom=0,Yprom=0,Sigma2,Rho;
  for(N=0, ix=0; ix <Lx;ix++){
    for(iy=0; iy <Ly;iy++){
      N+=rho(ix,iy);
    }
  }
  for(ix=0;ix < Lx;ix++){
    for(iy=0;iy < Ly;iy++){
      Rho=rho(ix,iy);
      Xprom+=ix*Rho; Yprom+=iy*Rho;
    }
  }

  Xprom/=N; Yprom/=N;
  for(Sigma2=0,ix=0;ix <Lx;ix++){
    for(iy=0;iy < Ly;iy++){
      Sigma2+=(pow(iy-Yprom,2.0)+pow(ix-Xprom,2.0))*rho(ix,iy);
    }
  }
  Sigma2/=(N-1);
  return Sigma2;
}

void LatticeGas::GrafiqueRho(){
  for(int ix=0;ix <Lx;ix++)
    std::cout <<ix <<"\t" <<std::endl;
}

void LatticeGas::GrafiqueDistribucion(std::string filename){
  int ix,iy;
  std::ofstream file(filename);
  for(ix=0; ix <Lx;ix++){
    for(iy=0;iy < Ly;iy++){
      for(int i=0;i <Q; i++){
        if(n[ix][iy][i] > 0){
          file << ix <<"\t" <<iy <<std::endl;
        }
      }
    }
  }
  file.close();
}

int main(){
  LatticeGas Difusion;
  Crandom ran64(1);
  int N=2400; double mu=Lx/2, sigma=16;
  int t, tmax=3500;

  Difusion.Borrar();
  Difusion.Inicie(N, mu, mu, sigma, sigma, ran64);
  for(t=0;t <tmax;t++){
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
    std::cout << t <<"\t"<< Difusion.Varianza() << std::endl;
    if(t%50==0) Difusion.GrafiqueDistribucion("data/Taller2/punto_1aintento2_dist"+std::to_string(t)+".dat");
  }
  //Difusion.Show();

  return 0;
}
