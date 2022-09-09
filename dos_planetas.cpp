#include <iostream>
#include <cmath>
#include "vector.h"

//Constantes globales
const int N=2;
const double g=9.8;

//consantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;

const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Lambda);

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

class Colisionador{
private:

public:
	void Force(Cuerpo *Planeta);
	void Force_betw(Cuerpo & Planeta1, Cuerpo & Planeta2);
};

void Colisionador::Force(Cuerpo *Planeta){
	for(int i=0;i<N;i++){
		Planeta[i].Force_erase();
		for(int j=i+1;j<N;j++)
			Force_betw(Planeta[i], Planeta[j]);
	}
}

void Colisionador::Force_betw(Cuerpo & Planeta1, Cuerpo & Planeta2){
	vector3D r21,n,F1; double d,F;
	r21=Planeta2.r-Planeta1.r; d=r21.norm(); n=r21/d;
	F=g*Planeta1.m*Planeta2.m*pow(d, -3);
	F1=F*n; Planeta1.Force_add(F1); Planeta2.Force_add(F1*(-1));


}

int main(){
	Cuerpo planeta[N];
	double m=1,r_i=10;
	double t, dt=0.01;
	double w=sqrt(g*m*pow(r_i,-3)), v_i=w*r_i, T=2*M_PI/w;

	planeta.Init(r_i,0,0,0,v_i,0,m,r_i);
	planeta.Force();
	for(t=0;t<1.1*T;t+=dt){
		std::cout << planeta.Getx() << "\t" << planeta.Gety() << std::endl;
		planeta.Move_r(dt,Zeta);
		planeta.Force();
		planeta.Move_V(dt,Coeficiente1);
		planeta.Move_r(dt,Chi);
		planeta.Force();
		planeta.Move_V(dt,Lambda);
		planeta.Move_r(dt,Coeficiente2);
		planeta.Force();
		planeta.Move_V(dt,Lambda);
		planeta.Move_r(dt,Chi);
		planeta.Force();
		planeta.Move_V(dt,Coeficiente1);
		planeta.Move_r(dt,Zeta);
	}
	return 0;
}