#include "LB_EM.hpp"

void max(std::vector<double> &V1, std::vector<double> &V2);

int main(){
  LatticeBoltzmann Conductor;
  int t, tmax=100;
  int N=150;
  std::vector<double> S_maxH(N,0), S_maxH_aux(N,0), S_maxE(N,0), S_maxE_aux(N,0);
  
  Conductor.Start();
  
  for(t=0;t<r/C-1;t++){
    Conductor.Collision(t);
    Conductor.Advection();
    if(t%5==0) Conductor.Print("data/lambda_"+std::to_string(t)+".dat", "data/lambda_perfil_"+std::to_string(t)+".dat");
  }
  for(t=r/C;t<r/C+T;t++){
    Conductor.Collision(t);
    Conductor.Advection();
    Conductor.Patron_PlanoH(N,S_maxH,t);
    Conductor.Patron_PlanoE(N,S_maxE,t);
    max(S_maxH, S_maxH_aux);
    max(S_maxE, S_maxE_aux);
    if(t%5==0) Conductor.Print("data/lambda_"+std::to_string(t)+".dat", "data/lambda_perfil_"+std::to_string(t)+".dat");

  }
  for(t=r/C+T+1;t<tmax;t++){
    Conductor.Collision(t);
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
