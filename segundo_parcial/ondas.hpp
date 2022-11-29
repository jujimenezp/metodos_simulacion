//Lattice-Boltzmann para ondas
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

const int Lx=2000;
const int Ly=1;

const int Q=5;
const double W0=1.0/3;

const double C=0.5; // C<0.707 cells/click
const double C2=C*C;
const double AUX0=1-3*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

const double lambda=1000; //Longitud de la onda entrante
const double A=1.; //Amplitud de la onda entrante
double D=0.6; //Factor de atenuacion de la pared reflectante, no se hace constante para el punto 4

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
  void Print(std::string NameFile);
  void envolventes(std::vector<double> &pmin, std::vector<double> &pmax);
};

LatticeBoltzmann::LatticeBoltzmann(){
  //Set the weights
  w[0]=W0;w[1]=w[2]=w[3]=w[4]=(1.0-W0)/4;
  //set the velocity vectors
  Vx[0]=0; Vx[1]=1; Vx[2]=0; Vx[3]=-1; Vx[4]=0;
  Vy[0]=0; Vy[1]=0; Vy[2]=1; Vy[3]=0; Vy[4]=-1;
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

  for(iy=0;iy<Ly;iy++){
    for(ix=0;ix<Lx-1;ix++){ //Se itera hasta la penultima celda en ix
      rho0=rho(ix,iy,false);
      Jx0=Jx(ix,iy,false);
      Jy0=Jy(ix,iy,false);
      for(i=0;i<Q;i++){
        n0=n(ix,iy,i);
        fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,i);
      }
    }
    ix=Lx-1; //Se implementa condicion de reflexion para las celdas en Lx-1
    n0=n(ix,iy,0); fnew[n0] = D*f[n0];
    n0=n(ix,iy,1); fnew[n0] = D*f[n(ix,iy,3)];
    n0=n(ix,iy,2); fnew[n0] = D*f[n(ix,iy,4)];
    n0=n(ix,iy,3); fnew[n0] = D*f[n(ix,iy,1)];
    n0=n(ix,iy,4); fnew[n0] = D*f[n(ix,iy,2)];
  }
}

void LatticeBoltzmann::ImposeFields(int t){
  int i,ix,iy,n0;
  double omega,rho0,Jx0,Jy0;
  omega=2*M_PI*C/lambda;
  ix=0;
  rho0=A*sin(omega*t);
  for(iy=0;iy < Ly;iy++){
    Jx0=Jx(ix,iy,false);
    Jy0=Jy(ix,iy,false);
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=feq(rho0,Jx0,Jy0,i);
    }
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

void LatticeBoltzmann::Print(std::string NameFile){
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

void LatticeBoltzmann::envolventes(std::vector<double> &pmin, std::vector<double> &pmax){
  double rho0;
  double iy=0;
  for(int ix=1;ix <Lx-1;ix++){
    rho0=rho(ix,iy,true);
    if(rho0 > pmax[ix-1]) pmax[ix-1]=rho0;
    if(rho0 < pmin[ix-1]) pmin[ix-1]=rho0;
  }
}
