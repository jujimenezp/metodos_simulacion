#include "LB_EM.hpp"

class Lambda2: public LatticeBoltzmann{
  public:
  double Jz_lambda2(int ix,int iy,int iz,double t);
  void Collision_lambda2(double t);
  vector3D Poynting_lambda2(int ix, int iy, int iz, double t);
  void Patron_PlanoH_lambda2(int N, std::vector<double> &S, double t);
  void Patron_PlanoE_lambda2(int N, std::vector<double> &S, double t);
};

void max(std::vector<double> &V1, std::vector<double> &V2);

int main(){
  Lambda2 Conductor;
  int t, tmax=100;
  int N=150;
  std::vector<double> S_maxH(N,0), S_maxH_aux(N,0), S_maxE(N,0), S_maxE_aux(N,0);
  
  Conductor.Start();
  
  for(t=0;t<r/C-1;t++){
    Conductor.Collision_lambda2(t);
    Conductor.Advection();
    if(t%5==0) Conductor.Print("data/lambda_"+std::to_string(t)+".dat", "data/lambda_perfil_"+std::to_string(t)+".dat");
  }
  for(t=r/C;t<r/C+T;t++){
    Conductor.Collision_lambda2(t);
    Conductor.Advection();
    Conductor.Patron_PlanoH_lambda2(N,S_maxH,t);
    Conductor.Patron_PlanoE_lambda2(N,S_maxE,t);
    max(S_maxH, S_maxH_aux);
    max(S_maxE, S_maxE_aux);
    if(t%5==0) Conductor.Print("data/lambda_"+std::to_string(t)+".dat", "data/lambda_perfil_"+std::to_string(t)+".dat");

  }
  for(t=r/C+T+1;t<tmax;t++){
    Conductor.Collision_lambda2(t);
    Conductor.Advection();
    if(t%5==0) Conductor.Print("data/lambda_"+std::to_string(t)+".dat", "data/lambda_perfil_"+std::to_string(t)+".dat");
  }

  Conductor.Print_Patron(N,S_maxH_aux,S_maxE_aux,"data/patronH.dat","data/patronE.dat");

  return 0;
}


void max(std::vector<double> &V1, std::vector<double> &V2){
  double size=V1.size();
  for(int i=0; i < size;i++){
    if(V1[i] > V2[i]){
      V2[i]=V1[i];
    }
  }
}

void Lambda2::Collision_lambda2(double t){
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
        Jprima0.load(0,0,Jz_lambda2(ix,iy,iz,t)); Eprima0=Eprima(E0,Jprima0,epsilonr0);
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

double Lambda2::Jz_lambda2(int ix, int iy, int iz, double t){
  double J, J0;
  J0=J00/2*(tanh(100*(iz-iz_ant+lambda/4))-tanh(100*(iz-iz_ant-lambda/4)));
  J=J0*cos(2*M_PI/lambda*fabs(iz-iz_ant))*sin(omega*t)*exp(-alpha*(pow(ix-ix_ant,2.)+pow(iy-iy_ant,2.)));
  return J;
}

vector3D Lambda2::Poynting_lambda2(int ix, int iy, int iz, double t){
  vector3D S;
  double sigma0, mur0, epsilonr0, prefactor0, rhoc0;
  vector3D D0, B0, E0, H0, Jprima0, Eprima0;
  //Compute the constants
  sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
  prefactor0=prefactor(epsilonr0,sigma0);
  //Compute the fields
  rhoc0=rhoc(ix,iy,iz,false); D0=D(ix,iy,iz,false); B0=B(ix,iy,iz,false);
  E0=E(D0,epsilonr0); H0=H(B0,mur0);
  Jprima0.load(0,0,Jz_lambda2(ix,iy,iz,t)); Eprima0=Eprima(E0,Jprima0,epsilonr0);
  S = (Eprima0^H0);
  return S;
}
void Lambda2::Patron_PlanoH_lambda2(int N, std::vector<double> &S, double t){
  double x,y,z=Lz/2,theta=0,S_max;
  int ix,iy,iz;
  double dtheta=2*M_PI/N;
  vector3D Sxy,Sxy1,Sx1y,Sx1y1;
  for(int i=0;i < N;i++,theta+=dtheta){
    x=ix_ant+r*cos(theta); y=iy_ant+r*sin(theta);
    ix=floor(x); iy=floor(y); iz=floor(z);
    Sxy=Poynting_lambda2(ix,iy,iz,t);
    Sxy1=Poynting_lambda2(ix,iy+1,iz,t);
    Sx1y=Poynting_lambda2(ix+1,iy,iz,t);
    Sx1y1=Poynting_lambda2(ix+1,iy+1,iz,t);;
    S[i]=inter_bilineal(x,y,Sxy,Sxy1,Sx1y,Sx1y1);
  }
}
void Lambda2::Patron_PlanoE_lambda2(int N, std::vector<double> &S, double t){
  double x,y=Ly/2,z,theta=0,S_max;
  int ix,iy,iz;
  double dtheta=2*M_PI/N;
  vector3D Sxz,Sxz1,Sx1z,Sx1z1;
  for(int i=0;i < N;i++,theta+=dtheta){
    x=ix_ant+r*cos(theta); z=iz_ant+r*sin(theta);
    ix=floor(x); iy=floor(y); iz=floor(z);
    Sxz=Poynting_lambda2(ix,iy,iz,t);
    Sxz1=Poynting_lambda2(ix,iy,iz+1,t);
    Sx1z=Poynting_lambda2(ix+1,iy,iz,t);
    Sx1z1=Poynting_lambda2(ix+1,iy,iz+1,t);
    S[i]=inter_bilineal(x,y,Sxz,Sxz1,Sx1z,Sx1z1);
  }
}
