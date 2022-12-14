//Conductores en 3D
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "../vector.h"

#ifndef LB_EM_H_
#define LB_EM_H_
//------------------------CONSTANTS-------------------------------
const int Lx = 200;   //
const int Ly = 200;   //
const int Lz = 200; //
const int Qr = 2, Qp = 3, Qi = 4, Qj = 2;
//-------------------
const double Tau = 0.5;
const double UTau = 1/Tau;
const double UmUTau=1-1/Tau;
//-------------------
const double Epsilon0=1, Mu0=2;
const double Sigma0=0;
const double C=1.0/sqrt(2.0);

const double E00=0.001,B00=E00/C,J00=0.0001;
const double Z0=sqrt(Mu0/Epsilon0);

const int ix_ant=Lx/2, iy_ant=Ly/2, iz_ant=Lz/2;
const double alpha = 0.5;
//Constantes lambda/2
const double lambda=20., T=lambda/C, omega=2*M_PI/T;
const double r=3*lambda;

//------------------Electromagnetic Constants for the Media------------------------------
double mur(int ix,int iy,int iz){
  return 1.0;
}
double epsilonr(int ix,int iy,int iz){
  return 1.0;
}
double parasito(int ix, int iy, int iz,
                int ix_par, int iy_par, double l, double alpha_par){
  double sigma;
  sigma=exp(-alpha_par*(pow(ix-(ix_ant+ix_par),2.)+pow(iy-(iy_ant+iy_par),2.)))*0.5*(tanh(100*(iz-(iz_ant-l/2)))-tanh(100*(iz-(iz_ant+l/2))));
  return sigma;
}

double sigma(int ix,int iy,int iz){
  double director1,director2,director3,reflector1,magnitud=16;
  director1=parasito(ix,iy,iz,0.125*lambda,0,0.55*lambda,0.5);
  director2=parasito(ix,iy,iz,(0.125+0.2)*lambda,0,0.4*lambda,0.5);
  director3=parasito(ix,iy,iz,(0.125+0.4)*lambda,0,0.35*lambda,0.5);

  reflector1=parasito(ix,iy,iz,-0.35*lambda,0,0.56*lambda,0.5);
  return magnitud*(director1+director2+director3+reflector1);
}

//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
  private:
    int V[3][3][4], V0[3]; /*V[xyz][p][i]*/  vector3D v[3][4],v0; //v[p][i]
    vector3D e[3][4][2], e0; //e[p][i][j]
    vector3D b[3][4][2], b0; //b[p][i][j]
    double *f=nullptr,*fnew=nullptr;//f[ix][iy][iz][r][p][i][j]
    double *f0 = nullptr,*f0new = nullptr;//f0[ix][iy][iz] (r=0)
  public:
  LatticeBoltzmann();
  ~LatticeBoltzmann();
  int index(int ix,int iy,int iz,int r,int p,int i,int j);
  int index0(int ix,int iy,int iz);
  //Auxiliary variables
  double prefactor(double epsilonr0,double sigma0){return sigma0/(1+(sigma0*Mu0)/(4*epsilonr0));};
  //Fields from direct sums
  double rhoc(int ix,int iy,int iz,bool UseNew);
  vector3D D(int ix,int iy,int iz,bool UseNew);
  vector3D B(int ix,int iy,int iz,bool UseNew);
  //Fields deduced from the first ones through electromagnetic constants
  vector3D E(vector3D & D0,double Epsilonr);
  vector3D H(vector3D & B0,double Mur);
  //Forced (Actual) Macroscopic Fields (inline)
  vector3D Jprima(vector3D & E0,double prefactor0){return E0*prefactor0;};
  vector3D Eprima(vector3D & E0,vector3D & Jprima0,double epsilonr0){return E0-Jprima0*(Mu0/(4*epsilonr0));};
  vector3D Dprima(vector3D & Eprima0,double epsilonr0){return Eprima0/epsilonr0;};
  double Jz(int ix, int iy, int iz, double t);
  vector3D Poynting(int ix, int iy, int iz, double t);
  //Equilibrium Functions
  double feq(vector3D & Jprima0,vector3D & Eprima0,vector3D & B0,
             double Epsilonr,double Mur,
             int r,int p,int i,int j);
  double feq0(double rhoc0);
  //Simulation Functions
  void Start(void);
  void Collision(double t);
  void ImposeFields(int t);
  void Advection(void);
  void Print(std::string filename, std::string filename_perfil);
  double inter_bilineal(double x, double y, vector3D Sxy, vector3D Sxy1, vector3D Sx1y,vector3D Sx1y1);
  void Patron_PlanoH(int N, std::vector<double> &S, double t);
  void Patron_PlanoE(int N, std::vector<double> &S, double t);
  void Print_Patron(int N,std::vector<double> S_maxH,std::vector<double> S_maxE,std::string filenameH, std::string filenameE);
  friend class Dipole;
  friend class Lambda2;
  friend class Yagi_Uda;
};


LatticeBoltzmann::LatticeBoltzmann(){
  int ix,iy,iz,alpha,r,p,i,j;
  //Velocity vectors V[p][i]=V^p_i (in components)
  V0[0]=V0[1]=V0[2]=0;

  V[0][0][0]=V[0][1][0]=V[1][2][0]=1;
  V[1][0][0]=V[2][1][0]=V[2][2][0]=1;
  V[2][0][0]=V[1][1][0]=V[0][2][0]=0;

  V[0][0][1]=V[0][1][1]=V[1][2][1]=-1;
  V[1][0][1]=V[2][1][1]=V[2][2][1]=1;
  V[2][0][1]=V[1][1][1]=V[0][2][1]=0;

  V[0][0][2]=V[0][1][2]=V[1][2][2]=-1;
  V[1][0][2]=V[2][1][2]=V[2][2][2]=-1;
  V[2][0][2]=V[1][1][2]=V[0][2][2]=0;

  V[0][0][3]=V[0][1][3]=V[1][2][3]=1;
  V[1][0][3]=V[2][1][3]=V[2][2][3]=-1;
  V[2][0][3]=V[1][1][3]=V[0][2][3]=0;
  //Velocity vectors V[p][i]=V^p_i (as vectors)
  v0.load(V0[0],V0[1],V0[2]); //cargue= load (in Spanish)
  for(p=0;p<3;p++)
    for(i=0;i<4;i++){
      v[p][i].load(V[0][p][i],V[1][p][i],V[2][p][i]);
    }
  //Electric vectors e[p][i][j]=e^p_{ij}
  e0.load(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++){
      e[p][i][0]=v[p][(i+1)%4]*0.5;
      e[p][i][1]=v[p][(i+3)%4]*0.5;
    }
  //Magnetic vectors b[p][i][j]=b^p_{ij}=v^p_i x e^p_{ij}
  b0.load(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
        b[p][i][j]=(v[p][i]^e[p][i][j]);

  f = new double[Lx*Ly*Lz*Qr*Qp*Qi*Qj]; fnew = new double[Lx*Ly*Lz*Qr*Qp*Qi*Qj];
  f0 = new double[Lx*Ly*Lz]; f0new=new double[Lx*Ly*Lz];
}

LatticeBoltzmann::~LatticeBoltzmann(void){
  delete[] f;  delete[] fnew;
  delete[] f0; delete[] f0new;
}

int LatticeBoltzmann::index(int ix,int iy,int iz,int r,int p,int i,int j){
  return (iz*Lx*Ly+iy*Lx+ix)*Qr*Qp*Qi*Qj + (r*Qp*Qi*Qj + p*Qi*Qj + i*Qj +j);
}

int LatticeBoltzmann::index0(int ix,int iy,int iz){
  return (iz*Lx*Ly+iy*Lx+ix);
}

//-----------------MACROSCOPIC FIELDS------------------
double LatticeBoltzmann::Jz(int ix, int iy, int iz, double t){
  double J, J0;
  J0=J00/2*(tanh(100*(iz-iz_ant+lambda/4))-tanh(100*(iz-iz_ant-lambda/4)));
  J=J0*cos(2*M_PI/lambda*fabs(iz-iz_ant))*sin(omega*t)*exp(-alpha*(pow(ix-ix_ant,2.)+pow(iy-iy_ant,2.)));
  return J;
}
vector3D LatticeBoltzmann::Poynting(int ix, int iy, int iz, double t){
  vector3D S;
  double sigma0, mur0, epsilonr0, prefactor0, rhoc0;
  vector3D D0, B0, E0, H0, Jprima0, Eprima0;
  //Compute the constants
  sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
  prefactor0=prefactor(epsilonr0,sigma0);
  //Compute the fields
  rhoc0=rhoc(ix,iy,iz,false); D0=D(ix,iy,iz,false); B0=B(ix,iy,iz,false);
  E0=E(D0,epsilonr0); H0=H(B0,mur0);
  Jprima0.load(0,0,Jz(ix,iy,iz,t)); Eprima0=Eprima(E0,Jprima0,epsilonr0);
  S = (Eprima0^H0);
  return S;
}

//Fields from direct sums
double LatticeBoltzmann::rhoc(int ix,int iy,int iz,bool UseNew){
  int p,i,j; double sum;
  int id0,id;
  id0 = index0(ix,iy,iz);
  //Start for the distribution for the central (zero) vector
  if(UseNew)
    sum=f0new[id0];
  else
    sum=f0[id0];
  //Add all the others
  for(p=0;p<2;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++){
        id = index(ix,iy,iz,0,p,i,j);
        if(UseNew)
          sum+=fnew[id];
        else
          sum+=f[id];
      }
  return sum;
}
vector3D LatticeBoltzmann::D(int ix,int iy,int iz,bool UseNew){
  int p,i,j,id; vector3D sum; sum.load(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++){
        id = index(ix,iy,iz,0,p,i,j);
        if(UseNew)
          sum+=e[p][i][j]*fnew[id];
        else
          sum+=e[p][i][j]*f[id];
      }
  return sum;
}
vector3D LatticeBoltzmann::B(int ix,int iy,int iz,bool UseNew){
  int p,i,j,id; vector3D sum; sum.load(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++){
        id = index(ix,iy,iz,1,p,i,j);
        if(UseNew)
          sum+=b[p][i][j]*fnew[id];
        else
          sum+=b[p][i][j]*f[id];
      }
  return sum;
}
//Fields deduced from the first ones through electromagnetic constants
vector3D LatticeBoltzmann::E(vector3D & D0,double Epsilonr){
  return D0*(1.0/Epsilonr);
}
vector3D LatticeBoltzmann::H(vector3D & B0,double Mur){
  return B0*(1.0/Mur);
}
//---------------EQUILIBRIUM FUNCTIONS-------------
double LatticeBoltzmann::feq(vector3D & Jprima0,vector3D & Eprima0,vector3D & B0,
                             double epsilonr0,double mur0,
                             int r,int p,int i,int j){
  double VdotJp=(v[p][i]*Jprima0),Epdote=(e[p][i][j]*Eprima0),Bdotb=(b[p][i][j]*B0),aux;
  if(r==0)
    aux=0.25*(0.25*VdotJp+epsilonr0*Epdote+0.5/mur0*Bdotb);
  if(r==1)
    aux=0.25*(0.25*VdotJp+Epdote+0.5*Bdotb);
  return aux;
}
double LatticeBoltzmann::feq0(double rhoc0){
  return rhoc0;
}

//-------------------SIMULATION FUNCTIONS ----------------------------
void LatticeBoltzmann::Start(void){
  int ix,iy,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  int id0,id;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        prefactor0=prefactor(epsilonr0,sigma0);
        //Impose the fields
        rhoc0=0; D0.load(0,0,0); B0.load(0,0,0);
        E0=E(D0,epsilonr0); H0=H(B0,mur0);
        Jprima0=Jprima(E0,prefactor0); Eprima0=Eprima(E0,Jprima0,epsilonr0);
        //Impose f=fnew=feq with the desired fields
        id0 = index0(ix,iy,iz);
        f0new[id0]=f0[id0]=feq0(rhoc0);
        for(r=0;r<2;r++)
          for(p=0;p<3;p++)
            for(i=0;i<4;i++)
              for(j=0;j<2;j++){
                id = index(ix,iy,iz,r,p,i,j);
                fnew[id]=f[id]=feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
              }
      }
}

void LatticeBoltzmann::Collision(double t){
  int ix,iy,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  int id0,id;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  for(ix=0;ix<Lx;ix++){ //para cada celda
    for(iy=0;iy<Ly;iy++){
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        prefactor0=prefactor(epsilonr0,sigma0);
        //Compute the fields
        rhoc0=rhoc(ix,iy,iz,false); D0=D(ix,iy,iz,false); B0=B(ix,iy,iz,false);
        E0=E(D0,epsilonr0); H0=H(B0,mur0);
        Jprima0.load(0,0,Jz(ix,iy,iz,t)); Eprima0=Eprima(E0,Jprima0,epsilonr0);
        //BGK evolution rule
        id0 = index0(ix,iy,iz);
        f0new[id0]=UmUTau*f0[id0]+UTau*feq0(rhoc0);
        for(r=0;r<2;r++)
          for(p=0;p<3;p++)
            for(i=0;i<4;i++)
              for(j=0;j<2;j++){
                id = index(ix,iy,iz,r,p,i,j);
                fnew[id]=UmUTau*f[id]+UTau*feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
              }
      }
    }
  }
}
// void LatticeBoltzmann::ImposeFields(int t){
//   int ix,iy,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0,J0,J;
//   int ix_ant, iy_ant, iz_ant;
//   int id0,id;
//   double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
//   double T=25.,omega=2*M_PI/T;//25<T<50
//   ix_ant=iy_ant=iz_ant=Lz/2;
//   for(iz=0;iz<Lz;iz++){ //for each cell
//     for(iy=0;iy<Ly;iy++){
//       for(ix=0;ix<Lx;ix++){
//         //Compute the constants
//         J0=J00*exp(-0.25*(pow(ix-ix_ant,2.)+pow(iy-iy_ant,2.)+pow(iz-iz_ant,2.)));
//         sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
//         prefactor0=prefactor(epsilonr0,sigma0);
//         J=J0*sin(omega*t);
//         //Compute the fields
//         //Primary fields (i.e. those from the sums)
//         rhoc0=rhoc(ix,iy,iz,false);
//         //D0.load(0,0,0);
//         //B0.load(0,0,0);
//         //Secundary fields (i.e. computed from the primary fields)
//         E0=E(D0,epsilonr0); H0=H(B0,mur0);
//         Jprima0.load(0,0,J);
//         Eprima0=Eprima(E0,Jprima0,epsilonr0);
//         //Impose fnew=feq with the desired fields
//         for(r=0;r<2;r++)
//           for(p=0;p<3;p++)
//             for(i=0;i<4;i++)
//               for(j=0;j<2;j++){
//                 id0 = index0(ix,iy,iz); id = index(ix,iy,iz,r,p,i,j);
//                 fnew[id]=feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
//                 f0new[id0]=feq0(rhoc0);
//               }
//       }
//     }
//   }
// }

void LatticeBoltzmann::Advection(void){
  int ix,iy,iz,r,p,i,j,ixnew,iynew,iznew;
  int id0,id,idnew;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){//for each cell
        for(r=0;r<2;r++)
          for(p=0;p<3;p++)
            for(i=0;i<4;i++)
              for(j=0;j<2;j++){
                ixnew=(ix+V[0][p][i]+Lx)%Lx; iynew=(iy+V[1][p][i]+Ly)%Ly; iznew=(iz+V[2][p][i]+Lz)%Lz;
                id0=index0(ix,iy,iz);
                id = index(ix,iy,iz,r,p,i,j);idnew=index(ixnew,iynew,iznew,r,p,i,j);
                f[idnew]=fnew[id];
                f0new[id0]=f0[id0];
              }
      }
}
void LatticeBoltzmann::Print(std::string filename, std::string filename_perfil){
  int ix,iy=Ly/2,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0; double E2,B2;
  std::ofstream file(filename);
  std::ofstream file2(filename_perfil);
  for(iz=0;iz < Lz; iz++)
    for(ix=0;ix<Lx;ix++){
      //Compute the electromagnetic constants
      sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
      prefactor0=prefactor(epsilonr0,sigma0);
      //Compute the Fields
      rhoc0=rhoc(ix,iy,iz,true); D0=D(ix,iy,iz,true); B0=B(ix,iy,iz,true);
      //Print
      B2 = B0.y()/J00;
      file<<ix<<"\t"<<iz<<"\t"<<B2<<std::endl;
      if(iz==Lz/2)
        file2 << ix << "\t" << B2 << std::endl;
    }
  file.close();
  file2.close();
}

double LatticeBoltzmann::inter_bilineal(double x, double y, vector3D Sxy, vector3D Sxy1, vector3D Sx1y,vector3D Sx1y1){
  int ix, iy, iz=Lz/2;
  double u, v, S_max;
  ix=floor(x); iy=floor(y);
  u=x-ix; v=y-iy;
  S_max=(1-u)*(1-v)*Sxy.norm()+(1-u)*v*Sxy1.norm()+u*(1-v)*Sx1y.norm()+u*v*Sx1y1.norm();
  return S_max;
}

void LatticeBoltzmann::Patron_PlanoH(int N, std::vector<double> &S, double t){
  double x,y,z=Lz/2,theta=0,S_max;
  int ix,iy,iz;
  double dtheta=2*M_PI/N;
  vector3D Sxy,Sxy1,Sx1y,Sx1y1;
  for(int i=0;i < N;i++,theta+=dtheta){
    x=ix_ant+r*cos(theta); y=iy_ant+r*sin(theta);
    ix=floor(x); iy=floor(y); iz=floor(z);
    Sxy=Poynting(ix,iy,iz,t);
    Sxy1=Poynting(ix,iy+1,iz,t);
    Sx1y=Poynting(ix+1,iy,iz,t);
    Sx1y1=Poynting(ix+1,iy+1,iz,t);;
    S[i]=inter_bilineal(x,y,Sxy,Sxy1,Sx1y,Sx1y1);
  }
}

void LatticeBoltzmann::Patron_PlanoE(int N, std::vector<double> &S, double t){
  double x,y=Ly/2,z,theta=0,S_max;
  int ix,iy,iz;
  double dtheta=2*M_PI/N;
  vector3D Sxz,Sxz1,Sx1z,Sx1z1;
  for(int i=0;i < N;i++,theta+=dtheta){
    x=ix_ant+r*cos(theta); z=iz_ant+r*sin(theta);
    ix=floor(x); iy=floor(y); iz=floor(z);
    Sxz=Poynting(ix,iy,iz,t);
    Sxz1=Poynting(ix,iy,iz+1,t);
    Sx1z=Poynting(ix+1,iy,iz,t);
    Sx1z1=Poynting(ix+1,iy,iz+1,t);
    S[i]=inter_bilineal(x,y,Sxz,Sxz1,Sx1z,Sx1z1);
  }
}

void LatticeBoltzmann::Print_Patron(int N, std::vector<double> S_maxH, std::vector<double> S_maxE,std::string filenameH, std::string filenameE){
  double x,y,theta=0;
  double dtheta=2*M_PI/N;
  std::ofstream file(filenameH);
  std::ofstream file2(filenameE);
  for(int i=0;i < N;i++,theta+=dtheta){
    x=ix_ant+r*cos(theta); y=iy_ant +r*sin(theta);
    file <<theta <<"\t"<< log10(S_maxH[i]/(Z0*J00*J00))<<std::endl;
    file2 <<theta <<"\t"<< log10(S_maxE[i]/(Z0*J00*J00))<<std::endl;
  }
  file.close();
  file2.close();
}

#endif // LB_EM_H_
