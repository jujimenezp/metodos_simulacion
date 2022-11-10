//Automata celular de difusion 2D
#include <cmath>
#include "../Random64.h"
#include "fstream"

const int Lx=256;
const int Ly=256;
const int L=Lx*Ly;
const double p0=0.25;
const double p=0.25;

const int Q=4;

//Clase LatticeGas
class LatticeGas{
private:
  int V[Q]; //V[i] i=0 derecha, i=1 arriba, i=2 izquierda, i=3 abajo
  int n[L][Q], nnew[L][Q]; //n[ix*iy][i]
public:
  LatticeGas(void);
  void Borrar();
  void Inicie(int N,double mux, double muy, double sigmax, double sigmay,Crandom & ran64);
  void Show();
  void GrafiqueRho();
  void Colisione(Crandom & ran64);
  void Adveccione();
  double rho(int ix){return n[ix][0]+n[ix][1]+n[ix][2]+n[ix][3];}
  double Varianza();
  void GrafiqueDistribucion(std::string filename);
};

LatticeGas::LatticeGas(){
  V[0]=V[3]=1; V[2]=V[1]=-1;
}

void LatticeGas::Borrar(){
  for(int ir=0; ir<L;ir++){
    for(int i=0;i<Q;i++){
      n[ir][i]=nnew[ir][i]=0;
    }
  }
}

void LatticeGas::Inicie(int N,double mux, double muy, double sigmax, double sigmay,Crandom & ran64){
  int ix, iy, i;
  while(N>0){
    ix=(int) ran64.gauss(mux,sigmax);//Escoger una celda al azar;
    iy=(int) ran64.gauss(muy,sigmay);//Escoger una celda al azar;

    if(ix<0) ix=0;
    if(iy<0) iy=0;
    if(ix>Lx-1) ix=Lx-1;
    if(iy>Ly-1) iy=Ly-1;
    i=(int) Q*ran64.r();
    if(n[Lx*iy+ix][i]==0){
      n[Lx*iy+ix][i]=1; N--;
    }
  }
}

void LatticeGas::Show(){
  for (int i=0; i < Q; i++) {
    for (int ir = 0; ir < L; ir++) {
      std::cout << n[ir][i];
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void LatticeGas::Colisione(Crandom & ran64){
  double ran=ran64.r();
  for(int ir=0;ir<L;ir++){
    //No rotar
    if(ran < p0){
      nnew[ir][0]=n[ir][0]; nnew[ir][1]=n[ir][1];
      nnew[ir][2]=n[ir][2]; nnew[ir][3]=n[ir][3];
    }
    //Rotar a la izquierda
    else if(ran >= p0 && ran <p+p0){
      nnew[ir][0]=n[ir][3]; nnew[ir][1]=n[ir][0];
      nnew[ir][2]=n[ir][1], nnew[ir][3]=n[ir][2];
    }
    //Rotar a la derecha
    else if(ran >= p+p0 && ran <2*p+p0){
      nnew[ir][0]=n[ir][1]; nnew[ir][1]=n[ir][2];
      nnew[ir][2]=n[ir][3], nnew[ir][3]=n[ir][0];
    }
    //Rotar 180 grados
    else {
      nnew[ir][0]=n[ir][2]; nnew[ir][1]=n[ir][3];
      nnew[ir][2]=n[ir][0], nnew[ir][3]=n[ir][1];
    }
  }
}

void LatticeGas::Adveccione(){
  int ix, iy;
  for(int ir=0; ir <L;ir++){
    ix=ir%Lx; iy=floor(ir/Lx);
    for(int i=0;i <Q;i++){
      //Cuando i=0 o i=2 se mueve en x
      if(i%2==0) {n[((ix+V[i]+Lx)%Lx)+Lx*iy][i]=nnew[ir][i];}
      else {n[(ir+Lx*V[i]+L)%L][i]=nnew[ir][i];}
    }
  }
}

double LatticeGas::Varianza(){
  int ir,ix,iy; double N,Xprom,Yprom,Sigma2,Rho;
  for(N=0, ir=0; ir <L;ir++){
    N+=rho(ir);
  }
  for(Xprom=0,Yprom=0,ir=0;ir <L;ir++){
    Rho=rho(ir);
    ix=ir%Lx; iy=floor(ir/Lx);
    Xprom+=ix*Rho; Yprom+=iy*Rho;
  }
  Xprom/=N; Yprom/=N;
  for(Sigma2=0,ir=0;ir <L;ir++){
    ix=ir%Lx; iy=floor(ir/Lx);
    Sigma2+=(pow(ix-Xprom,2.0)+pow(iy-Yprom,2.0))*rho(ir);
  }
  Sigma2/=(N-1);
  return Sigma2;
}

void LatticeGas::GrafiqueRho(){
  for(int ir=0;ir <L;ir++)
    std::cout <<ir <<"\t" <<std::endl;
}

void LatticeGas::GrafiqueDistribucion(std::string filename){
  int ix,iy;
  std::ofstream file(filename);
  for(int ir=0; ir <L;ir++){
    for(int i=0;i <Q; i++){
      if(n[ir][i] > 0){
        ix=ir%Lx; iy=floor(ir/Lx);
        file << ix <<"\t" <<iy <<std::endl;
      }
    }
  }
  file.close();
}

int main(){
  LatticeGas Difusion;
  Crandom ran64(1);
  int N=2400; double mu=Lx/2, sigma=16;
  int t, tmax=350;

  Difusion.Borrar();
  Difusion.Inicie(N, mu, mu, sigma, sigma, ran64);
  for(t=0;t <tmax;t++){
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
    std::cout << t << "\t" << Difusion.Varianza() <<std::endl;
    if(t%50==0) Difusion.GrafiqueDistribucion("data/Taller2/punto_1a_dist"+std::to_string(t)+".dat");
  }
  //Difusion.Show();

  return 0;
}
