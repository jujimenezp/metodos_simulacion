//Skin Effect in 1D
#include<iostream>
#include<cmath>
#include "Vector.h"
using namespace std;

//------------------------CONSTANTS-------------------------------
const int Lx = 1;   //
const int Ly = 1;   //
const int Lz = 1000; //
//-------------------
const double Tau = 0.5;
const double UTau = 1/Tau;
const double UmUTau=1-1/Tau;
//-------------------
const double Epsilon0=1, Mu0=2;
const double Sigma0=0.0125;
const double C=1.0/sqrt(2.0);

const double E00=0.001,B00=E00/C;

//------------------Electromagnetic Constants for the Media------------------------------
double mur(int ix,int iy,int iz){
  return 1.0;
}
double epsilonr(int ix,int iy,int iz){
  return 1.0;
}
double sigma(int ix,int iy,int iz){
  return Sigma0*(tanh(iz-100)-tanh(iz-900))/2;
}

//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
private:
  int V[3][3][4], V0[3]; /*V[xyz][p][i]*/  vector3D v[3][4],v0; //v[p][i]
  vector3D e[3][4][2], e0; //e[p][i][j]
  vector3D b[3][4][2], b0; //b[p][i][j]
  double f[Lx][Ly][Lz][2][3][4][2],fnew[Lx][Ly][Lz][2][3][4][2];//f[ix][iy][iz][r][p][i][j]
  double f0[Lx][Ly][Lz],f0new[Lx][Ly][Lz];//f0[ix][iy][iz] (r=0)
public:
  LatticeBoltzmann(void);
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
  //Equilibrium Functions
  double feq(vector3D & Jprima0,vector3D & Eprima0,vector3D & B0,
	     double Epsilonr,double Mur,
	     int r,int p,int i,int j);
  double feq0(double rhoc0);
  //Simulation Functions
  void Start(void);
  void Collision(void);
  void ImposeFields(int t);
  void Advection(void);
  void Print(void);
};

LatticeBoltzmann::LatticeBoltzmann(void){
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
  v0.cargue(V0[0],V0[1],V0[2]); //cargue= load (in Spanish)
  for(p=0;p<3;p++)
    for(i=0;i<4;i++){
      v[p][i].cargue(V[0][p][i],V[1][p][i],V[2][p][i]);
  }
  //Electric vectors e[p][i][j]=e^p_{ij}
  e0.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++){
      e[p][i][0]=v[p][(i+1)%4]*0.5;
      e[p][i][1]=v[p][(i+3)%4]*0.5;
  }
  //Magnetic vectors b[p][i][j]=b^p_{ij}=v^p_i x e^p_{ij}
  b0.cargue(0,0,0);  
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	b[p][i][j]=(v[p][i]^e[p][i][j]);
}
//-----------------MACROSCOPIC FIELDS------------------
//Fields from direct sums
double LatticeBoltzmann::rhoc(int ix,int iy,int iz,bool UseNew){
  int p,i,j; double sum;
  //Start for the distribution for the central (zero) vector
  if(UseNew) 
    sum=f0new[ix][iy][iz]; 
  else 
    sum=f0[ix][iy][iz];
  //Add all the others
  for(p=0;p<2;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	if(UseNew) 
	  sum+=fnew[ix][iy][iz][0][p][i][j]; 
	else 
	  sum+=f[ix][iy][iz][0][p][i][j];
  return sum;
}
vector3D LatticeBoltzmann::D(int ix,int iy,int iz,bool UseNew){
  int p,i,j; vector3D sum; sum.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	if(UseNew) 
	    sum+=e[p][i][j]*fnew[ix][iy][iz][0][p][i][j];
	else
	    sum+=e[p][i][j]*f[ix][iy][iz][0][p][i][j];
  return sum;
}
vector3D LatticeBoltzmann::B(int ix,int iy,int iz,bool UseNew){
  int p,i,j; vector3D sum; sum.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	if(UseNew) 
	  sum+=b[p][i][j]*fnew[ix][iy][iz][1][p][i][j];
	else
	  sum+=b[p][i][j]*f[ix][iy][iz][1][p][i][j];
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
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
	//Compute the constants
	sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
	prefactor0=prefactor(epsilonr0,sigma0);
	//Impose the fields
	rhoc0=0; D0.cargue(0,0,0); B0.cargue(0,0,0);
	E0=E(D0,epsilonr0); H0=H(B0,mur0);
	Jprima0=Jprima(E0,prefactor0); Eprima0=Eprima(E0,Jprima0,epsilonr0); 
	//Impose f=fnew=feq with the desired fields
	f0new[ix][iy][iz]=f0[ix][iy][iz]=feq0(rhoc0);
	for(r=0;r<2;r++)
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++)
		fnew[ix][iy][iz][r][p][i][j]=f[ix][iy][iz][r][p][i][j]=
		  feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
      }
}
void LatticeBoltzmann::Collision(void){
  int ix,iy,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
	//Compute the constants
	sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
	prefactor0=prefactor(epsilonr0,sigma0);
	//Compute the fields
	rhoc0=rhoc(ix,iy,iz,false); D0=D(ix,iy,iz,false); B0=B(ix,iy,iz,false);
	E0=E(D0,epsilonr0); H0=H(B0,mur0);
	Jprima0=Jprima(E0,prefactor0); Eprima0=Eprima(E0,Jprima0,epsilonr0); 
	//BGK evolution rule
	f0new[ix][iy][iz]=UmUTau*f0[ix][iy][iz]+UTau*feq0(rhoc0);
	for(r=0;r<2;r++)
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++)
		fnew[ix][iy][iz][r][p][i][j]=UmUTau*f[ix][iy][iz][r][p][i][j]
		  +UTau*feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
      }
}
void LatticeBoltzmann::ImposeFields(int t){
  int ix,iy,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  double T=25.0,omega=2*M_PI/T;//25<T<50
  iz=0; //On the whole incident plane
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      //Compute the constants
      sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
      prefactor0=prefactor(epsilonr0,sigma0);
      //Compute the fields
      //Primary fields (i.e. those from the sums)
      rhoc0=rhoc(ix,iy,iz,false); 
      D0.cargue(E00*sin(omega*t)*epsilonr0,0,0);
      B0.cargue(0,B00*sin(omega*t),0);
      //Secundary fields (i.e. computed from the primary fields)
      E0=E(D0,epsilonr0); H0=H(B0,mur0);
      Jprima0=Jprima(E0,prefactor0); Eprima0=Eprima(E0,Jprima0,epsilonr0); 
      //Impose fnew=feq with the desired fields
      for(r=0;r<2;r++)
	for(p=0;p<3;p++)
	  for(i=0;i<4;i++)
	    for(j=0;j<2;j++){
	      fnew[ix][iy][iz][r][p][i][j]=feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
	      f0new[ix][iy][iz]=feq0(rhoc0);
	    }
    }
}
void LatticeBoltzmann::Advection(void){
  int ix,iy,iz,r,p,i,j,ixnew,iynew,iznew;
  for(ix=0;ix<Lx;ix++) 
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){//for each cell
	for(r=0;r<2;r++)
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++){
		ixnew=(ix+V[0][p][i]+Lx)%Lx; iynew=(iy+V[1][p][i]+Ly)%Ly; iznew=(iz+V[2][p][i]+Lz)%Lz;
		f[ixnew][iynew][iznew][r][p][i][j]=fnew[ix][iy][iz][r][p][i][j];
		f0new[ix][iy][iz]=f0[ix][iy][iz];
	      }
      }
}
void LatticeBoltzmann::Print(void){
  int ix=0,iy=0,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0; double E2,B2;
  for(iz=0;iz<Lz/2;iz++){
    //Compute the electromagnetic constants
    sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
    prefactor0=prefactor(epsilonr0,sigma0);
    //Compute the Fields
    rhoc0=rhoc(ix,iy,iz,true); D0=D(ix,iy,iz,true); B0=B(ix,iy,iz,true);
    E0=E(D0,epsilonr0); H0=H(B0,mur0);
    Jprima0=Jprima(E0,prefactor0); Eprima0=Eprima(E0,Jprima0,epsilonr0); 
    //Print
    cout<<iz<<" "<<E0.x()/E00<<endl;
  }
}


int main(){
  LatticeBoltzmann OndaSkin;
  int t, tmax=700;
  
  OndaSkin.Start();
  OndaSkin.ImposeFields(0);
  
  for(t=0;t<tmax;t++){
    OndaSkin.Collision();
    OndaSkin.ImposeFields(t);
    OndaSkin.Advection();
  }
  
  OndaSkin.Print();
  
  return 0;
}
