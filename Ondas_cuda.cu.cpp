#include <iostream>
#include <fstream>
#include <cmath>

#define Lx 128
#define Ly 128
#define N 32 //Threads per block
const int M=(Lx*Ly+N-1)/N; //Blocks per Grid
#define Q 5
const int ArraySize=Lx*Ly*Q;

const double W0=1.0/3;

const double C=0.5; // C<0.707 cells/click
const double C2=C*C;
const double AUX0=1-3*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//------------ PROGRAMMING ON THE DEVICE ----------------
//---------------Constants (Symbols)----------------
__constant__ float d_w[5];
__constant__ int d_Vx[5];
__constant__ int d_Vy[5];
__constant__ float d_C[3];   // d_C[0]=C,  d_C[1]=C2,  d_C[2]=AUX,
__constant__ float d_tau[3]; // d_tau[0]=tau,  d_tau[1]=Utau,  d_tau[2]=UmUtau,

//Clase LatticeGas
class LatticeBoltzmann{
private:
  double w[Q]; //Weights
  int Vx[Q],Vy[Q]; //Velocity vectors
  double *f, *fnew; //Distribution functions

public:
  LatticeBoltzmann();
  ~LatticeBoltzmann();
  int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;}
  double rho(int ix, int iy, bool UseNew);
  double Jx(int ix,int iy, bool UseNew);
  double Jy(int ix,int iy,bool UseNew);
  double feq(double rho0,double Jx0,double Jy0,int i);
  void Start(double rho0,double Jx0,double Jy0);
  void Collision();
  void ImposeFields(int t);
  void Advection();
  void Print(const char * NameFile);
};

LatticeBoltzmann::LatticeBoltzmann(){
  //CONSTANTS(d_Symbols)
  //---Charge constantes on the Host-----------------
  //running constants
  h_C[0]=C;  h_C[1]=C2;  h_C[2]=AUX0;
  h_tau[0]=tau;  h_tau[1]=Utau;  h_tau[2]=UmUtau;
  //Set the weights
  h_w[0]=W0; h_w[1]=h_w[2]=h_w[3]=h_w[4]=(1.0-W0)/4;
  //Set the velocity vectors
  h_Vx[0]=0;  h_Vx[1]=1;  h_Vx[2]=0;  h_Vx[3]=-1; h_Vx[4]=0;
  h_Vy[0]=0;  h_Vy[1]=0;  h_Vy[2]=1;  h_Vy[3]=0;  h_Vy[4]=-1;
  //------Send to the Device-----------------
  cudaMemcpyToSymbol(d_w,h_w,Q*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vx,h_Vx,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy,h_Vy,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_C,h_C,3*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_tau,h_tau,3*sizeof(float),0,cudaMemcpyHostToDevice);

  //Create the dynamic arrays
  int ArraySize=Lx*Ly*Q;
  f= new double [ArraySize]; fnew= new double [ArraySize];
}

LatticeBoltzmann::~LatticeBoltzmann(){
  delete[] f; delete[] fnew;
}

double LatticeBoltzmann::rho(int ix, int iy, bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=fnew[n0];
    else sum+=f[n0];
  }
  return sum;
}

double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vx[i]*fnew[n0];
    else sum+=Vx[i]*f[n0];
  }
  return sum;
}

double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vy[i]*fnew[n0];
    else sum+=Vy[i]*f[n0];
  }
  return sum;
}

double LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,int i){
  if(i>0)
    return 3*w[i]*(C2*rho0+Vx[i]*Jx0+Vy[i]*Jy0);
  else
    return rho0*AUX0;
}

void LatticeBoltzmann::Start(double rho0,double Jx0,double Jy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(i=0;i<Q;i++){
        n0=n(ix,iy,i);
        f[n0]=feq(rho0,Jx0,Jy0,i);
      }
    }
  }
}

void LatticeBoltzmann::Collision(){
  int ix,iy,i,n0; double rho0,Jx0,Jy0;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      Jx0=Jx(ix,iy,false);
      Jy0=Jy(ix,iy,false);
      for(i=0;i<Q;i++){
        n0=n(ix,iy,i);
        fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,i);
      }
    }
  }
}

void LatticeBoltzmann::ImposeFields(int t){
  int i,ix,iy,n0;
  double lambda,omega,rho0,Jx0,Jy0;
  lambda=10;
  omega=2*M_PI/lambda*C;
  ix=Lx/2;iy=Ly/2;
  rho0=10*sin(omega*t);
  Jx0=Jx(ix,iy,false);
  Jy0=Jy(ix,iy,false);
  for(i=0;i<Q;i++){
    n0=n(ix,iy,i);
    fnew[n0]=feq(rho0,Jx0,Jy0,i);
  }
}

void LatticeBoltzmann::Advection(){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(i=0;i<Q;i++){
        ixnext=(ix+Vx[i]+Lx)%Lx;
        iynext=(iy+Vy[i]+Ly)%Ly;
        n0=n(ix,iy,i);
        n0next=n(ixnext,iynext,i);
        f[n0next]=fnew[n0];
      }
    }
  }
}

void LatticeBoltzmann::Print(const char * NameFile){
  std::ofstream MyFile(NameFile); double rho0; int ix, iy;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true);
      MyFile <<ix <<"\t" <<iy <<"\t"<<rho0 <<std::endl;
    }
    MyFile <<std::endl;
  }
  MyFile.close();
}

int main(){
  LatticeBoltzmann Ondas;
  int t,tmax=4000;
  double rho0=0,Jx0=0,Jy0=0;

  Ondas.Start(rho0,Jx0,Jy0);
  for(t=0;t<tmax;t++){
    Ondas.Collision();
    Ondas.ImposeFields(t);
    Ondas.Advection();
  }
  Ondas.Print("data/ondas.dat");

  return 0;
}
