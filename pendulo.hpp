#ifndef PENDULO_H_
#define PENDULO_H_

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "../Vector3D/vector.h"

//Constantes globales
const int N=3;
const double G=980;

//constantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Declaración de las clases
class Cuerpo;
class Colisionador;

//---------- Clase Cuerpo --------------
class Cuerpo{
private:
  double Theta,Omega,T,m,R,l,I,x0;
public:
  void Inicie(double Theta0,double omega0,double T0,double m0,double R0,double l0,double x00);
  void BorreTorque(void){T=0;};
  void SumeTorque(double T0){T+=T0;};
  void Mueva_Theta(double dt,double coeficiente);
  void Mueva_Omega(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return x0+l*sin(Theta);}; //Inline
  double Gety(void){return -l*cos(Theta);}; //Inline
  double GetT(void){return T;}
  friend class Colisionador;
};
void Cuerpo::Inicie(double Theta0,double Omega0, double T0,
  double m0,double R0,double l0,double x00){
  m=m0; R=R0; Theta=Theta0; Omega=Omega0; T=T0; l=l0; x0=x00; I=m*l*l;
}
void Cuerpo::Mueva_Theta(double dt,double coeficiente){
  Theta+=Omega*(dt*coeficiente);
}
void Cuerpo::Mueva_Omega(double dt,double coeficiente){
  Omega+=T*(dt*coeficiente/I);
}
void Cuerpo::Dibujese(void){
  std::cout<<" , "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety() <<"+"<<R<<"*sin(t)";
  std::cout<<" , "<<x0<<"+"<<l/7<<"*t*sin("<<Theta<<"),-"<<l/7<<"*t*cos("<<Theta<<")";
}

//---------- Clase Colisionador --------------
class Colisionador{
private:
  double K;
public:
  double GetK(){return K;}
  void Actualizar_K(double kv){K=kv;}
  void CalculeTorques(Cuerpo * Pendulo);
  void CalculeTorqueEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2);
};
void Colisionador::CalculeTorques(Cuerpo * Pendulo){
  int i;
  //Borrar fuerzas
  double T0;
  for(i=0;i<N;i++){
    Pendulo[i].BorreTorque();
    T0=-Pendulo[i].l*Pendulo[i].m*G*sin(Pendulo[i].Theta);
    Pendulo[i].SumeTorque(T0);
  }
  //Calcular las fuerzas entre todas las parejas de planetas
  for(i=N-1;i>0;i--)
    CalculeTorqueEntre(Pendulo[i],Pendulo[i-1]);
}
 void Colisionador::CalculeTorqueEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2){
   double s=(Pendulo2.Getx()+Pendulo2.R)-(Pendulo1.Getx()-Pendulo1.R); double F=0;
   if(s>0) F=K*pow(s,1.5); //Fuerza de Hertz
   Pendulo1.SumeTorque(F*Pendulo1.l); Pendulo2.SumeTorque(-F*Pendulo2.l);
 }

//----------- Funciones para graficar -----------

void InicieAnimacion(void){
  //  std::cout<<"set terminal gif animate"<<std::endl;
  //  std::cout<<"set output 'DosPendulos.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set xrange[-5:10]"<<std::endl;
  std::cout<<"set yrange[-14:1]"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange [0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(void){
    std::cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    std::cout<<std::endl;
}

//Paso de evolucion
void step(Cuerpo *Pendulo, Colisionador Newton, double dt){
  // Mover por PEFRL
  int i=0;
  for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Zeta);
  Newton.CalculeTorques(Pendulo);
  for(i=0;i<N;i++) Pendulo[i].Mueva_Omega(dt,Coeficiente1);
  for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Chi);
  Newton.CalculeTorques(Pendulo);
  for(i=0;i<N;i++) Pendulo[i].Mueva_Omega(dt,Lambda);
  for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Coeficiente2);
  Newton.CalculeTorques(Pendulo);
  for(i=0;i<N;i++) Pendulo[i].Mueva_Omega(dt,Lambda);
  for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Chi);
  Newton.CalculeTorques(Pendulo);
  for(i=0;i<N;i++) Pendulo[i].Mueva_Omega(dt,Coeficiente1);
  for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Zeta);
}

//calculo de torque maximo y tiempo de media oscilacion
double Maximos(double &Torques, double t, double &tmax, double &tmin, double Tau, double &Tau_min){
     if(Tau > Torques){
       Torques=Tau;
       tmax=t;
     }
     else if(Tau < Tau_min){
       Tau_min=Tau;
       tmin=t;
     }
     return fabs(tmax-tmin);
}
#endif // PENDULO_H_
