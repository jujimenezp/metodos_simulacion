//Fuerza de arrastre
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

const int Lx=512;
const int Ly=64;

const int Q=9;

const double tau=1.5; //Preguntar
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

const double viscosity=(tau-0.5)/3;

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
  double feq(double rho0,double Ux0,double Uy0,int i);
  void Start(double rho0,double Ux0,double Uy0);
  void Collision();
  void ImposeFields(double Ufan, int ixc, int iyc, int R);
  void Advection();
  void Print(const char * NameFile, double Ufan);
  double sigmaxx(int ix, int iy);
  double sigmayy(int ix, int iy);
  double sigmaxy(int ix, int iy);
  double dFx(double x, double y, double dAx, double dAy);
  double dFy(double x, double y, double dAx, double dAy);
  std::vector<double> Fcilindro(int N, int  ixc, int  iyc, int R);
};

LatticeBoltzmann::LatticeBoltzmann(){
  //Set the weights
  w[0]=4.0/9;  w[1]=w[2]=w[3]=w[4]=1.0/9;  w[5]=w[6]=w[7]=w[8]=1.0/36;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;

  Vx[5]=1;  Vx[6]=-1; Vx[7]=-1; Vx[8]=1;
  Vy[5]=1;  Vy[6]=1;  Vy[7]=-1; Vy[8]=-1;

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

double  LatticeBoltzmann::feq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*Vx[i]+Uy0*Vy[i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}

void LatticeBoltzmann::Start(double rho0,double Ux0,double Uy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(i=0;i<Q;i++){
        n0=n(ix,iy,i);
        f[n0]=feq(rho0,Ux0,Uy0,i);
      }
    }
  }
}

void LatticeBoltzmann::Collision(){
  int ix,iy,i,n0; double rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      //compute the macroscopic fields on the cell
      rho0=rho(ix,iy,false);
      Ux0=Jx(ix,iy,false)/rho0;
      Uy0=Jy(ix,iy,false)/rho0;
      for(i=0;i<Q;i++){
        n0=n(ix,iy,i);
        fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i);
      }
    }
  }
}

void LatticeBoltzmann::ImposeFields(double Ufan, int ixc, int iyc, int R){
  int i,ix,iy,n0;
  double rho0, R2=R*R;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      //fan
      if(ix==0)
        for(i=0;i<Q;i++) {n0=n(ix,iy,i); fnew[n0]=feq(rho0,Ufan,0,i);}
      //obstacle
      else if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R2)
        for(i=0;i<Q;i++) {n0=n(ix,iy,i); fnew[n0]=feq(rho0,0,0,i);}
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
        f[n0next]=fnew[n0]; //periodic boundaries
      }
    }
  }
}

void LatticeBoltzmann::Print(const char * NameFile,double Ufan){
  std::ofstream MyFile(NameFile); double rho0,Ux0,Uy0; int ix,iy;
  for(ix=0;ix<Lx;ix+=4){
    for(iy=0;iy<Ly;iy+=4){
      rho0=rho(ix,iy,true); Ux0=Jx(ix,iy,true)/rho0; Uy0=Jy(ix,iy,true)/rho0;
      MyFile<<ix<< "\t"<< iy<< "\t"<< Ux0/Ufan*4<<"\t"<< Uy0/Ufan*4<<std::endl;
    }
    MyFile<<std::endl;
  }
  MyFile.close();
}

double LatticeBoltzmann::sigmaxx(int ix, int iy){
  double rho0, sigma=0;
  rho0=rho(ix,iy,true);
  for(int i=0;i < Q;i++){
    sigma+=w[i]*Vx[i]*Jx(ix+Vx[i],iy+Vy[i],true)/rho0;
  }
  sigma*=6*viscosity; sigma-=rho0/3;
  return sigma;
}
double LatticeBoltzmann::sigmayy(int ix, int iy){
  double rho0, sigma=0;
  rho0=rho(ix,iy,true);
  for(int i=0;i < Q;i++){
    sigma+=w[i]*Vy[i]*Jy(ix+Vx[i],iy+Vy[i],true)/rho0;
  }
  sigma*=6*viscosity; sigma-=rho0/3;
  return sigma;
}
double LatticeBoltzmann::sigmaxy(int ix, int iy){
  double rho0=rho(ix,iy,true);
  double sigma=0;
  for(int i=0;i < Q;i++){
    sigma+=w[i]*(Vy[i]*Jx(ix+Vx[i],iy+Vy[i],true)/rho0+Vx[i]*Jy(ix+Vx[i],iy+Vy[i],true)/rho0);
  }
  sigma*=3*viscosity;
  return sigma;
}

double LatticeBoltzmann::dFx(double x, double y, double dAx, double dAy){
  double Sxx, Sxy, u, v; int ix, iy;
  //Celda para x,y
  ix=floor(x); iy=floor(y);
  u=x-ix; v=y-iy;
  //Interpolacion del tensor de esfuerzos
  Sxx=sigmaxx(ix,iy)*(1-u)*(1-v)+sigmaxx(ix+1,iy)*u*(1-v)+sigmaxx(ix,iy+1)*(1-u)*v+sigmaxx(ix+1,iy+1)*u*v;
  Sxy=sigmaxy(ix,iy)*(1-u)*(1-v)+sigmaxy(ix+1,iy)*u*(1-v)+sigmaxy(ix,iy+1)*(1-u)*v+sigmaxy(ix+1,iy+1)*u*v;
  //dF=sigma(x,y)*dA
  return Sxx*dAx+Sxy*dAy;
}

double LatticeBoltzmann::dFy(double x, double y, double dAx, double dAy){
  double Syy, Sxy, u, v; int ix, iy;
  //Celda para x,y
  ix=floor(x); iy=floor(y);
  u=x-ix; v=y-iy;
  //Interpolacion del tensor de esfuerzos
  Syy=sigmayy(ix,iy)*(1-u)*(1-v)+sigmayy(ix+1,iy)*u*(1-v)+sigmayy(ix,iy+1)*(1-u)*v+sigmayy(ix+1,iy+1)*u*v;
  Sxy=sigmaxy(ix,iy)*(1-u)*(1-v)+sigmaxy(ix+1,iy)*u*(1-v)+sigmaxy(ix,iy+1)*(1-u)*v+sigmaxy(ix+1,iy+1)*u*v;
  //dF=sigma(x,y)*dA
  return Sxy*dAx+Syy*dAy;
}

std::vector<double> LatticeBoltzmann::Fcilindro(int N, int  ixc, int  iyc, int R){ //Revisar todo
  std::vector<double> F(2,0); double x,y,dAx,dAy,theta;
  for(int i=0; i < N;i++){
    theta=2*M_PI*i/N;
    x=ixc+R*cos(theta);
    y=iyc+R*sin(theta);
    dAx=pow(2*M_PI*R/N,2.0)*cos(theta);
    dAy=pow(2*M_PI*R/N,2.0)*sin(theta);
    F[0]+=dFx(x,y,dAx,dAy);
    F[1]+=dFy(x,y,dAx,dAy);
  }
  return F;
}

int main(){
  LatticeBoltzmann Aire;
  int t,tmax=10000;
  double rho0=1., Ufan=0.1;
  int ixc=128, iyc=32, R=8;
  int Ncilindro=25;
  std::vector<double> F(2,0);

  Aire.Start(rho0,Ufan,0);
  for(t=0;t<tmax;t++){
    Aire.Collision();
    Aire.ImposeFields(Ufan,ixc,iyc,R);
    Aire.Advection();
    if(t==1) Aire.Print("data/Taller2/punto_3a_distribución0.dat",Ufan);
    F=Aire.Fcilindro(Ncilindro,ixc,iyc,R);
    std::cout << t << "\t" << F[0] << "\t" <<F[1] << std::endl;
  }
  Aire.Print("data/Taller2/punto_3a_distribución.dat",Ufan);

  return 0;
}
