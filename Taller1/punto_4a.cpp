#include "../pendulo.hpp"

int main(){
  Cuerpo Pendulo[N];
  Colisionador Newton;
  double m0=100, l0=12, R=1.5;
  //double M=m0+m1, x0=-m1*r/M, x1=m0*r/M;
  //double omega=sqrt(G*M/(r*r*r)), T=2*M_PI/omega, V0=omega*x0, V1=omega*x1;
  double T=2*M_PI*sqrt(l0/G);
  double t,tf=0.6*T,dt=1e-6;
  double tdibujo,tcuadro=T/1000;
  int i;
  double k_values[]={0.1e10, 0.2e10, 0.5e10, 1e10, 2.5e10, 5e10, 10e10};
  int k_size = sizeof(k_values)/sizeof(k_values[0]);
  double Torques[k_size] = {0},Tiempos[k_size]={0}, Tau_min, tmin, tmax;

  //-------Inicie( Theta0, omega0, T0, m0, R0, l0, x00)
  Pendulo[0].Inicie(-M_PI/12, 0, 0, m0, R, l0, 0);
  for(i=1;i<N;i++){
    Pendulo[i].Inicie(0, 0, 0, m0, R, l0, 2*R*i);
  }

 std::stringstream ss;
 std::string filename;

 InicieAnimacion();
 std::ofstream file_maximos("data/Taller1/punto_4b.dat");
 for(int ii=0; ii<k_size; ii++){
   Newton.Actualizar_K(k_values[ii]);
   Tau_min=0;

   //-------Inicie( Theta0, omega0, T0, m0, R0, l0, x00)
   Pendulo[0].Inicie(-M_PI/12, 0, 0, m0, R, l0, 0);
   for(i=1;i<N;i++){
     Pendulo[i].Inicie(0, 0, 0, m0, R, l0, 2*R*i);
   }

   //Archivo de salida
   ss <<"data/Taller1/punto_4a_"<< Newton.GetK() <<".dat";
   filename=ss.str();
   filename.erase(std::remove(filename.begin(), filename.end(), '+'), filename.end());
   std::ofstream file(filename);
   std::clog <<"K: " <<Newton.GetK() <<std::endl;


   for(t=0,tdibujo=0; t<tf; t+=dt,tdibujo+=dt){
     //Dibujar
     if(tdibujo>tcuadro){
       InicieCuadro();
       for(i=0;i<N;i++) Pendulo[i].Dibujese();
       TermineCuadro();
       tdibujo=0;
     }

     step(Pendulo, Newton, dt);
     Tiempos[ii]=Maximos(Torques[ii], t, tmax, tmin, Pendulo[1].GetT(), Tau_min);
     file << t << "\t" << Pendulo[1].GetT() <<std::endl;
   }
   file_maximos << Newton.GetK() <<"\t"<< Torques[ii] << "\t" <<Tiempos[ii] <<std::endl;
   file.close();
   ss.str("");
 }
 file_maximos.close();

 return 0;
}
