#include <iostream>
#include <cmath>
#include "../../Vector3D/vector.h"

//Constantes globales
const int N=1;

//consantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Clases
class Cuerpo;
class Colisionador;

class Cuerpo{
private:
	vector3D r,V,F;
	double m,R;
public:
	void Init(double x0,double y0, double z0,
			  double Vx0,double Vy0, double Vz0,
			  double m0,double R0);
	void Force_erase();
	void Force_add(vector3D F0);
	void Move_r(double dt, double coeff);
	void Move_V(double dt, double coeff);
	double Getx(){return r.x();}
	double Gety(){return r.y();}
	double Getz(){return r.z();}
	double Get_rnorm(){return r.norm();}
	void Dibujese(void);
	friend class Colisionador;
};
void Cuerpo::Init(double x0,double y0, double z0,double Vx0,double Vy0, double Vz0,double m0,double R0){
	r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0);
	m=m0; R=R0;
}

void Cuerpo::Force_erase(){
	F.load(0,0,0);
}

void Cuerpo::Force_add(vector3D F0){
	F+=F0;
}

void Cuerpo::Move_r(double dt, double coeff){
	r+=V*(dt*coeff);
}

void Cuerpo::Move_V(double dt, double coeff){
	V+=F*(dt*coeff/m);
}

void Cuerpo::Dibujese(void){
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

class Colisionador{
private:
public:
	void Force(Cuerpo *Planeta);
	void Force_betw(Cuerpo & Planeta1, Cuerpo & Planeta2);
};

void Colisionador::Force(Cuerpo *Planeta){
	//fuerza de Lennard-Jones central
	vector3D F0, n;
	double epsilon=1,r0=10,d;
	for(int i=0;i<N;i++){
		Planeta[i].Force_erase();
		d=Planeta[i].r.norm(); n=Planeta[i].r*(1/d);
		F0=n*((pow(r0/d,12)-pow(r0/d,6))*12*epsilon/d);
		Planeta[i].Force_add(F0);
	}
	for(int i=0;i<N;i++){
		for(int j=i+1;j<N;j++){
			Force_betw(Planeta[i], Planeta[j]);
		}
	}
}

void Colisionador::Force_betw(Cuerpo & Planeta1, Cuerpo & Planeta2){

}

//----------- Funciones Globales -----------

void InicieAnimacion(void){
  std::cout<<"set terminal gif animate"<<std::endl;
  std::cout<<"set output 'data/Taller1/punto_5a.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set xrange[0:15]"<<std::endl;
  std::cout<<"set yrange[-5:5]"<<std::endl;
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

// void Graficar_rotado(Cuerpo *planetas){
// 	double x_jupiter,y_jupiter,xpos,ypos;
// 	x_jupiter=planeta[1].Getx();
// 	y_jupiter=(planeta[1]
// 	for(int i=0; i<N; i++){
// 		xpos=(planeta)
// 	}
// }


int main(){
	Cuerpo planeta[N];
	Colisionador Lennard_Jones;
	double m=1;
	double x0=10;
	double V0=sqrt(2*0.5/m);
	double t,tmax=100, dt=0.01;
	double tdibujo,tcuadro=tmax/500;
	int i;

	planeta[0].Init(x0,0,0,V0,0,0,m,2.5);

	InicieAnimacion();
	for(t=0,tdibujo=0; t<tmax; t+=dt,tdibujo+=dt){
		//Dibujar animación
		if(tdibujo>tcuadro){
			InicieCuadro();
			for(i=0;i<N;i++) planeta[i].Dibujese();
			TermineCuadro();
			tdibujo=0;
		}

		//Dibujar trayectoria
		//std::cout << t <<"\t"<< planeta[0].Getx() << std::endl;

		//Move by PEFRL
		for(i=0;i<N;i++) planeta[i].Move_r(dt,Zeta);
		Lennard_Jones.Force(planeta);
		for(i=0;i<N;i++) planeta[i].Move_V(dt,Coeficiente1);
		for(i=0;i<N;i++) planeta[i].Move_r(dt,Chi);
		Lennard_Jones.Force(planeta);
		for(i=0;i<N;i++) planeta[i].Move_V(dt,Lambda);
		for(i=0;i<N;i++) planeta[i].Move_r(dt,Coeficiente2);
		Lennard_Jones.Force(planeta);
		for(i=0;i<N;i++) planeta[i].Move_V(dt,Lambda);
		for(i=0;i<N;i++) planeta[i].Move_r(dt,Chi);
		Lennard_Jones.Force(planeta);
		for(i=0;i<N;i++) planeta[i].Move_V(dt,Coeficiente1);
		for(i=0;i<N;i++) planeta[i].Move_r(dt,Zeta);
	}
	return 0;
}
