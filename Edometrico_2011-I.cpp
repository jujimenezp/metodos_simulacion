//UnPlaneta por Velocidad Verlet Optimizado
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;

const double Deltat=0.01;
const double chi=0.193183325037836;

const double g=9.8,VEL0=10;
const double K=500, Gamma=0.3, Kcundall=500, mu=0.4;

const double Lx=100,Ly=100;
const int Nx=5,Ny=5,N=Nx*Ny;

enum Restriccion{ninguna,posicion,velocidad,fuerza};

class Cuerpo;
class Colisionador;
//--------------------class Cuerpo -----------------------
class Cuerpo{
private:
  double m,R; vector3D r,V,F;  double theta,omega,tau,I;
  Restriccion RestriccionX,RestriccionY,RestriccionTheta;  double ValorX,ValorY,ValorTheta;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double theta0,double omega0,double m0,double R0);
  void Mueva_r1(double dt);
  void Mueva_V(double dt);
  void Mueva_r2(double dt);
  void Dibujese(void);

  void Fijex(double x0){RestriccionX=posicion; ValorX=x0;};
  void Fijey(double y0){RestriccionY=posicion; ValorY=y0;};
  void Fijetheta(double theta0){RestriccionTheta=posicion; ValorTheta=theta0;};
  void FijeVx(double Vx0){RestriccionX=velocidad; ValorX=Vx0;};
  void FijeVy(double Vy0){RestriccionY=velocidad; ValorY=Vy0;};
  void Fijeomega(double omega0){RestriccionTheta=velocidad; ValorTheta=omega0;};
  void FijeFx(double Fx0){RestriccionX=fuerza; ValorX=Fx0;};
  void FijeFy(double Fy0){RestriccionY=fuerza; ValorY=Fy0;};
  void Fijetau(double tau0){RestriccionTheta=fuerza; ValorTheta=tau0;};

  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  double Getm(void){return m;};
  double GetEc(void){return m*norma2(V)/2+I*omega*omega/2;};
  double GetV(void){return norma(V);};
  double GetVx(void){return V.x();};
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double theta0,double omega0,double m0,double R0){
  r.cargue(x0,y0,0); V.cargue(Vx0,Vy0,0); m=m0; R=R0; theta=theta0; omega=omega0; I=2.0/5*m*R*R;
  RestriccionX=RestriccionY=RestriccionTheta=ninguna;       ValorX=ValorY=ValorTheta=0;
}
void Cuerpo::Mueva_r1(double dt){
  r+=V*chi*dt;            theta+=omega*chi*dt;
  if(RestriccionX==posicion) r.v[0]=ValorX;  
  if(RestriccionY==posicion) r.v[1]=ValorY;  
  if(RestriccionTheta==posicion) theta=ValorTheta;
}
void Cuerpo::Mueva_V(double dt){
  V+=F*(dt/(2*m));      omega+=tau/I*dt/2;
  if(RestriccionX==velocidad) V.v[0]=ValorX;  
  if(RestriccionY==velocidad) V.v[1]=ValorY;  
  if(RestriccionTheta==velocidad) omega=ValorTheta;
}
void Cuerpo::Mueva_r2(double dt){
  r+=V*((1-2*chi)*dt);  theta+=omega*(1-2*chi)*dt;
  if(RestriccionX==posicion) r.v[0]=ValorX;  
  if(RestriccionY==posicion) r.v[1]=ValorY;  
  if(RestriccionTheta==posicion) theta=ValorTheta;
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t"; 
}
//--------------------class Colisionador -----------------------
class Colisionador{
private:
  double dcontacto[N+4][N+4],hold[N+4][N+4];
public:
  void Inicie(void);
  void CalculeFuerzas(Cuerpo *Grano,double dt);
  void Choque(Cuerpo &Cuerpo1,Cuerpo &Cuerpo2,double &dcontacto,double &hold,double dt);
  double GetEp(void);
};
double Colisionador::GetEp(void){
  double Ep=0,h;
  for(int i=0;i<N;i++)  
    for(int j=i+1;j<N+4;j++){
      Ep+=Kcundall*pow(dcontacto[i][j],2)/2;
      h=hold[i][j]; if(h>0) Ep+K*pow(h,3.5)/3.5;
    }
  return Ep;
}
void Colisionador::Inicie(void){
  for(int i=0;i<N+4;i++)  
    for(int j=0;j<N+4;j++)
      dcontacto[i][j]=hold[i][j]=0;
}
void Colisionador::CalculeFuerzas(Cuerpo *Grano,double dt){
  int i,j;
  for(i=N+1;i<N+3;i++)
      {Grano[i].F.cargue(0,0,0); Grano[i].tau=0;}
  Grano[N+3].F.cargue(0,-Grano[N+3].m*g,0); Grano[i].tau=0;
  for(i=0;i<N;i++)  
    {Grano[i].F.cargue(0,-Grano[i].m*g,0); Grano[i].tau=0;}
  for(i=0;i<N;i++)  
    for(j=i+1;j<N+4;j++)
      Choque(Grano[i],Grano[j],dcontacto[i][j],hold[i][j],dt);
}
void Colisionador::Choque(Cuerpo &Cuerpo1,Cuerpo &Cuerpo2,double &dcontacto,double &hold,double dt){
  vector3D r12,V12,Vc,n,Rw,t; double h,R1,R2,m1,m2,m12,Vn,Vt,Fn,Ft,Ftmax;
  R1=Cuerpo1.R;  R2=Cuerpo2.R;
  r12=Cuerpo1.r-Cuerpo2.r;
  h=R1+R2-norma(r12);

  if(h>0){//si hay contacto
    V12=Cuerpo1.V-Cuerpo2.V;
    n=r12/norma(r12);  Rw.cargue(0,0,R1*Cuerpo1.omega+R2*Cuerpo2.omega);
    Vc=V12-(Rw^n); Vn=Vc*n; t.cargue(n.y(),-n.x(),0); Vt=Vc*t;

    m1=Cuerpo1.m;    m2=Cuerpo2.m; m12=m1*m2/(m1+m2);
    Fn=K*pow(h,1.5)-m12*sqrt(h)*Gamma*Vn;

    dcontacto+=Vt*dt;
    Ft=-Kcundall*dcontacto;
    if(fabs(Ft)>(Ftmax=mu*fabs(Fn))) Ft=Ft/fabs(Ft)*Ftmax;
    
    Cuerpo1.F+=n*Fn;    Cuerpo1.F+=t*Ft;    Cuerpo1.tau+=R1*Ft;
    Cuerpo2.F-=n*Fn;    Cuerpo2.F-=t*Ft;    Cuerpo2.tau+=R2*Ft;

    if(Cuerpo1.RestriccionX==fuerza) Cuerpo1.F.v[0]=Cuerpo1.ValorX;
    if(Cuerpo1.RestriccionY==fuerza) Cuerpo1.F.v[1]=Cuerpo1.ValorY;
    if(Cuerpo1.RestriccionTheta==fuerza) Cuerpo1.tau=Cuerpo1.ValorTheta;
    if(Cuerpo2.RestriccionX==fuerza) Cuerpo2.F.v[0]=Cuerpo1.ValorX;
    if(Cuerpo2.RestriccionY==fuerza) Cuerpo2.F.v[1]=Cuerpo1.ValorY;
    if(Cuerpo2.RestriccionTheta==fuerza) Cuerpo2.tau=Cuerpo1.ValorTheta;
  }

  if(hold>0 && h<0){
    dcontacto=0;
  }

  hold=h;
  
}

//--------------------Funciones de Animacion -----------------------

void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'Relajacion.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-10:110]"<<endl;
  cout<<"set yrange [-10:110]"<<endl;
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
  Cuerpo Grano[N+4];
  Colisionador Cundall;
  Crandom Ran2(0);
  double t,tdibujo,dx,dy,Rmax,R,m,Alpha,Ec,Ep;  
  int i;

  InicieAnimacion();  
  Cundall.Inicie();
  //-------------(x0 ,y0,Vx0,Vy0,theta0,omega0,m0,R0);
  //PAREDES
  //Pared izquierda
  Grano[N].Inicie(-10000,Ly/2,0,0,0,0,100000,10000); 
  Grano[N].Fijex(-10000);  Grano[N].Fijey(Ly/2);  Grano[N].Fijetheta(0);
  //Pared derecha
  Grano[N+1].Inicie(Lx+10000,Ly/2,0,0,0,0,100000,10000); 
  Grano[N+1].Fijex(Lx+10000);  Grano[N+1].Fijey(Ly/2);  Grano[N+1].Fijetheta(0);
  //Pared abajo
  Grano[N+2].Inicie(Lx/2,-10000,0,0,0,0,100000,10000); 
  Grano[N+2].Fijex(Lx/2);  Grano[N+2].Fijey(-10000);  Grano[N+2].Fijetheta(0);
  //Pared arriba
  Grano[N+3].Inicie(Lx/2,Ly+10000,0,0,0,0,1,10000); 
  Grano[N+3].Fijex(Lx/2);  /*Grano[N+3].Fijey(Ly+10000);*/   Grano[N+3].Fijetheta(0);
  //GRANOS
  dx=Lx/(Nx+1); dy=Ly/(Ny+1); Rmax=dx/3; if(dx/3 > dy/3) Rmax=dy/3; 
  for(i=0;i<N;i++){
    Alpha=2*M_PI*Ran2.r();   R=Rmax*(1+Ran2.r())/2; m=pow(R/Rmax,3);
    Grano[i].Inicie(((i%Nx)+1)*dx,((i/Ny)+1)*dy,VEL0*cos(Alpha),VEL0*sin(Alpha),0,0,m,R);  
  }
  //RELAJACION
  for(t=tdibujo=0;t<60;t+=Deltat,tdibujo+=Deltat){
    if(tdibujo>20/120.0){
      InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
      }
    //Velocidad Verlet Optimizado
    for(i=0;i<N+4;i++) Grano[i].Mueva_r1(Deltat);
    Cundall.CalculeFuerzas(Grano,0.5*Deltat);
    for(i=0;i<N+4;i++) {Grano[i].Mueva_V(Deltat);Grano[i].Mueva_r2(Deltat);}
    Cundall.CalculeFuerzas(Grano,0.5*Deltat);
    for(i=0;i<N+4;i++) {Grano[i].Mueva_V(Deltat);Grano[i].Mueva_r1(Deltat);}
  }
  //COMPRESION EDOMETRICA
  //Pared arriba
  Grano[N+3].Inicie(Lx/2,Ly+10000,0,0,0,0,100,10000); 
  Grano[N+3].Fijex(Lx/2);  /*Grano[N+3].Fijey(Ly+10000);*/   Grano[N+3].Fijetheta(0);
  for(t=tdibujo=0;t<60;t+=Deltat,tdibujo+=Deltat){
    if(tdibujo>20/120.0){
      InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
      }
    clog<<t<<Grano[N+3].Gety()-10000<<endl;
    //Velocidad Verlet Optimizado
    for(i=0;i<N+4;i++) Grano[i].Mueva_r1(Deltat);
    Cundall.CalculeFuerzas(Grano,0.5*Deltat);
    for(i=0;i<N+4;i++) {Grano[i].Mueva_V(Deltat);Grano[i].Mueva_r2(Deltat);}
    Cundall.CalculeFuerzas(Grano,0.5*Deltat);
    for(i=0;i<N+4;i++) {Grano[i].Mueva_V(Deltat);Grano[i].Mueva_r1(Deltat);}
  }

  return 0;
}
