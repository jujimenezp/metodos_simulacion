// Simular el movimiento de N moleculas en un gas 2D
#include <iostream>
#include <cmath>
#include "../../Vector3D/vector.h"
#include "../Random64.h"
#include <vector>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_histogram.h>
using namespace std;

//---- declarar constantes ---
const double K=1.0e4;
const double Lx=60, Ly=120;
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
  double GetVx(void){return V.x();};
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
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      << r.x() << "+" << R*cos(theta)/7.0<<"*t, " <<r.y() <<"+"<< R*sin(theta)/7.0 <<"*t";
}

//--- clase Colisionador ----
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Molecula);
  void CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2);
};

void Colisionador::CalculeFuerzas(Cuerpo * Molecula){
  int i,j;
  vector3D F0;
  //--- Borrar todas las fuerzas ---
  for(i=0;i<N;i++){
    Molecula[i].BorreFuerza();
    //Fuerzas de pared en x=0, x=Lx, y=0, y=Ly
    double s_Lx = Molecula[i].Getx()+Molecula[i].R -Lx;
    double s_x0= Molecula[i].R-Molecula[i].Getx();
    double s_Ly = Molecula[i].Gety()+Molecula[i].R -Ly;
    double s_y0= Molecula[i].R-Molecula[i].Gety();
    if(s_Lx>0){
      F0.load(-K*pow(s_Lx,1.5),0,0);
      Molecula[i].AdicioneFuerza(F0);
    }
    if(s_x0>0){
      F0.load(K*pow(s_x0,1.5),0,0);
      Molecula[i].AdicioneFuerza(F0);
      }
    if(s_Ly>0){
      F0.load(0,-K*pow(s_Ly,1.5),0);
      Molecula[i].AdicioneFuerza(F0);
    }
    if(s_y0>0){
      F0.load(0,K*pow(s_y0,1.5),0);
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
  Molecula2.AdicioneFuerza(F2);   Molecula1.AdicioneFuerza(F2*(-1));

}

//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'data/Taller1/punto_5f_K10.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}

void Graficar_yprom(double t, Cuerpo *Molecula){
  double S=0;
  for(int ii=0; ii<N; ii++){
    S+=Molecula[ii].Gety();
  }
  S/=N;
  std::cout << t <<"\t"<< S <<std::endl;
}

void Graficar_Vx(double t, Cuerpo *molecula){
  std::cout <<t <<"\t";
  for(int ii=0; ii<N; ii++){
    std::cout << molecula[ii].GetVx() << "\t";
  }
  std::cout << std::endl;
}

void Save_Vx(std::vector<double> &vec, Cuerpo *Molecula){
  for(int ii=0; ii<N; ii++){
    vec.push_back(Molecula[ii].GetVx());
  }
}

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Molecula[N];
  Colisionador Lennard_Jones;
  Crandom ran64(1);
  double m=1, R0=2.5, kT=10, V0=sqrt(2*kT/m);
  int i;
  double t,tdibujo,tf=200,tcuadro=tf/1000,dt=1e-3;
  double dx=10, dy=10;
  double Theta, OmegaMax=1.0;
  bool histo=true;
  std::vector<double> Vx;


  //Inicializar las moléculas
  for(int ix=0;ix<Nx;ix++)
    for(int iy=0;iy<Ny;iy++){
      Theta=2*M_PI*ran64.r();
      //-----------------------(   x0,   y0,      Vx0, Vy0, m0,R0, theta0,omega0)
      Molecula[Nx*iy+ix].Inicie((ix+1)*dx,(iy+1)*dy, V0*cos(Theta),V0*sin(Theta), m,R0,   0, OmegaMax);//OJO
    }

  //InicieAnimacion(); //Dibujar
  for(t=0,tdibujo=0  ; t<tf ; t+=dt,tdibujo+=dt){
    // //Dibujar
    // if(tdibujo>tcuadro){

    //   InicieCuadro();
    //   for(i=0;i<N;i++) Molecula[i].Dibujese();
    //   TermineCuadro();

    //   tdibujo=0;
    // }

    //Graficar_yprom(t,Molecula);

    if(t>60){
      Graficar_Vx(t,Molecula);
      if(histo==true) Save_Vx(Vx, Molecula);
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

  if(histo==true){
    double Vx_dev = gsl_stats_sd(&Vx[0],1,Vx.size());
    double Vx_std = gsl_stats_mean(&Vx[0],1,Vx.size());
    std::clog << "Promedio: " << Vx_std << std::endl
              << "Desviacion: " << Vx_dev << std::endl;
  }
  
  return 0;
}

  
