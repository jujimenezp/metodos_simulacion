#include "LB_EM.hpp"

void max(std::vector<double> &V1, std::vector<double> &V2);

int main(){
  LatticeBoltzmann Conductor;
  int t, tmax=100;
  int N=100;
  std::vector<double> S_max(N,0), S_max_aux(N,0);
  
  Conductor.Start();
  
  for(t=0;t<r/C-1;t++){
    Conductor.Collision(t);
    Conductor.Advection();
    if(t%5==0) Conductor.Print("data/lambda_"+std::to_string(t)+".dat", "data/lambda_perfil_"+std::to_string(t)+".dat");
  }
  for(t=r/C;t<r/C+T;t++){
    Conductor.Collision(t);
    Conductor.Advection();
    Conductor.Patron_Radio(N,S_max,t);
    max(S_max, S_max_aux);
    if(t%5==0) Conductor.Print("data/lambda_"+std::to_string(t)+".dat", "data/lambda_perfil_"+std::to_string(t)+".dat");

  }
  for(t=r/C+T+1;t<tmax;t++){
    Conductor.Collision(t);
    Conductor.Advection();
    if(t%5==0) Conductor.Print("data/lambda_"+std::to_string(t)+".dat", "data/lambda_perfil_"+std::to_string(t)+".dat");
  }

  Conductor.Print_Patron(N,S_max_aux,"data/Patron.dat");

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
