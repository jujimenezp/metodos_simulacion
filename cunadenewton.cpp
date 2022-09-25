#include <iostream>
#include <cmath>
#include "../Vector3D/vector.h"
using namespace std;

//Constantes globales

const int N=3;
const double G=980;
const double K=1e9;

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
  cout<<" , "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety() <<"+"<<R<<"*sin(t)";
  cout<<" , "<<x0<<"+"<<l/7<<"*t*sin("<<Theta<<"),-"<<l/7<<"*t*cos("<<Theta<<")";
}
//---------- Clase Colisionador --------------
class Colisionador{
private:
public:
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

//----------- Funciones Globales -----------

void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'DosPendulos.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-5:10]"<<endl;
  cout<<"set yrange[-14:1]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}


int main(){
  Cuerpo Pendulo[N];
  Colisionador Newton;
  double m0=100, l0=12, R=1.5;
  //double M=m0+m1, x0=-m1*r/M, x1=m0*r/M;
  //double omega=sqrt(G*M/(r*r*r)), T=2*M_PI/omega, V0=omega*x0, V1=omega*x1;
  double T=2*M_PI*sqrt(l0/G);
  double t,tmax=3*T,dt=1e-6;
  double tdibujo,tcuadro=T/1000;
  int i;
  
  //-------Inicie( Theta0, omega0, T0, m0, R0, l0, x00)
  Pendulo[0].Inicie(-M_PI/12, 0, 0, m0, R, l0, 0);
  for(i=1;i<N;i++){
    Pendulo[i].Inicie(0, 0, 0, m0, R, l0, 2*R*i);
  }
  
  InicieAnimacion();
  
  for(t=0,tdibujo=0; t<tmax; t+=dt,tdibujo+=dt){
    //Dibujar
    if(tdibujo>tcuadro){
      InicieCuadro();
      for(i=0;i<N;i++) Pendulo[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }         
    
    //cout<<Pendulo[1].Getx()<<" "<<Pendulo[1].Gety()<<endl;
    // Mover por PEFRL
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
  
  return 0;
}
