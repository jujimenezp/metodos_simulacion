// Simular el movimiento de N moleculas en un gas 2D de Lennard-Jones
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include "vector.h"
#include "Random64.h"

//---- declarar constantes ---
const double K=1.0e4;
const double Rpared=50;
const int Nx=5, Ny=5, N=Nx*Ny;

const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--- declarar clases -----
class Cuerpo;
class Colisionador;

//---- interface e implementacion de clases ----
//---- clase cuerpo ---
class Cuerpo{
private:
  vector3D r,V,F; double m,R; double theta,omega,tau,I;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0,
              double theta0, double omega0);
  void BorreFuerza(){F.load(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0, double theta0, double omega0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0, theta=theta0; omega=omega0; I=2.0/5*m*R*R;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt); theta+=omega*(Coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m); omega+=tau*(Coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      << r.x() << "+" << R*cos(theta)/7.0<<"*t, " <<r.y() <<"+"<< R*sin(theta)/7.0 <<"*t";
}

//--- clase Colisionador ----
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Molecula);
  void CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2);
  vector3D r_cm(Cuerpo * Molecula);
  double r_giro(Cuerpo * Molecula);
};

void Colisionador::CalculeFuerzas(Cuerpo * Molecula){
  int i,j;
  vector3D F0;
  //--- Borrar todas las fuerzas ---
  for(i=0;i<N;i++){
    Molecula[i].BorreFuerza();
    //Fuerzas de pared con una circunferencia centrada en el origen
    double s=Molecula[i].r.norm()+Molecula[i].R-Rpared;
    vector3D n=Molecula[i].r*(1/Molecula[i].r.norm());
    if(s>0){
      F0=n*(-K*s);
      Molecula[i].AdicioneFuerza(F0);
    }
  }
  //--- Calcular Fuerzas entre pares de cuerpos ---
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      CalculeFuerzaEntre(Molecula[i], Molecula[j]);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2){
  vector3D r21=Molecula2.r-Molecula1.r;
  double d=r21.norm();
  double epsilon=1, r0=10;
  vector3D n=r21*(1.0/d);

  //Fuerza de Lennard-Jones
  vector3D F2=n*((pow(r0/d,12)-pow(r0/d,6))*12*epsilon/d);
  Molecula2.AdicioneFuerza(F2); Molecula1.AdicioneFuerza(F2*(-1));

}

//Calculo del centro de masa
vector3D Colisionador::r_cm(Cuerpo * Molecula){
  vector3D r_cm(0,0,0);
  for(int ii=0;ii<N; ii++){
    r_cm+=Molecula[ii].r;
  }
  r_cm=r_cm*(1/N);
  return r_cm;
}

double Colisionador::r_giro(Cuerpo * Molecula){
  vector3D r, aux;
  double R_giro=0;
  r=r_cm(Molecula);
  for(int ii=0;ii<N;ii++){
    aux=Molecula[ii].r-r;
    R_giro+=aux.norm2();
  }
  R_giro=sqrt(R_giro/N);
  return R_giro;
}

//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  std::cout<<"set terminal gif animate"<<std::endl;
  std::cout<<"set output 'data/2.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set xrange[-"<<Rpared+10<<":"<<Rpared+10<<"]"<<std::endl;
  std::cout<<"set yrange[-"<<Rpared+10<<":"<<Rpared+10<<"]"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange [0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(void){
    std::cout<<"plot " <<Rpared <<"*cos(2*pi*t/7),"
             <<Rpared <<"*sin(2*pi*t/7)";
}
void TermineCuadro(void){
    std::cout<<std::endl;
}

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Molecula[N];
  Colisionador Lennard_Jones;
  Crandom ran64(1);
  double m=1, R0=3, kT=10, V0=sqrt(2*kT/m);
  int i;
  double t,tdibujo,tf=100,tcuadro=tf/2000,dt=5e-4;
  double dx=10, dy=10;
  double Theta, OmegaMax=1.0;
  double R_giro, R_masa, giro_prom=0;
  int prom_N=0;
  const double kT_values[]={0.05,0.5,10};
  std::string filename;

  for(int j=0; j<3; j++){
    //Archivo de salida
    std::stringstream ss;
    ss <<"data/3_"<<kT_values[j]<<".dat";
    filename=ss.str();
    std::ofstream resultados(filename);

    prom_N=0; giro_prom=0;
    kT=kT_values[j]; V0=sqrt(2*kT/m);
    //Inicializar las moléculas
    for(int ix=0;ix<Nx;ix++)
      for(int iy=0;iy<Ny;iy++){
        Theta=2*M_PI*ran64.r();
        //-----------------------(   x0,   y0,      Vx0, Vy0, m0,R0, theta0,omega0)
        Molecula[Nx*iy+ix].Inicie((ix-2)*dx,(iy-2)*dy, V0*cos(Theta),V0*sin(Theta), m,R0,   0, OmegaMax);
      }

    for(t=0,tdibujo=0; t<tf ; t+=dt,tdibujo+=dt){

      R_giro=Lennard_Jones.r_giro(Molecula);
      resultados <<t <<"\t" <<R_giro <<std::endl;
      if(t>20){
        giro_prom+=R_giro;
        prom_N++;
      }

      //--- Muevase por PEFRL ---
      for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,epsilon);
      Lennard_Jones.CalculeFuerzas(Molecula);
      for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda2);
      for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,chi);
      Lennard_Jones.CalculeFuerzas(Molecula);
      for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda);
      for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,chiepsilon);
      Lennard_Jones.CalculeFuerzas(Molecula);
      for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda);
      for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,chi);
      Lennard_Jones.CalculeFuerzas(Molecula);
      for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda2);
      for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,epsilon);

    }
    giro_prom*=1./prom_N;
    std::clog <<"El radio de giro promedio es: " << giro_prom <<std::endl;
    resultados.close();
  }
  
  return 0;
}

  
