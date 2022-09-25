#include <iostream>
#include <cmath>
#include "../Vector3D/vector.h"

double g=9.8;

//consantes de forest Roth
const double Theta = 1/(2-pow(2., 1/3));
const double Theta_half = Theta/2;
const double UmTheta_half = (1-Theta)/2;
const double Um2Theta=1-2*Theta;

class Cuerpo{
private:
	vector3D r,V,F;
	double m,R;
public:
	void Init(double x0,double y0, double z0,
			  double Vx0,double Vy0, double Vz0,
			  double m0,double R0);
	void Force();
	void Move_r(double dt, double coeff);
	void Move_V(double dt, double coeff);
	double Getx(){return r.x();}
	double Gety(){return r.y();}
	double Getz(){return r.z();}
};
void Cuerpo::Init(double x0,double y0, double z0,double Vx0,double Vy0, double Vz0,double m0,double R0){
	r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0);
	m=m0; R=R0;
}

void Cuerpo::Force(){
	double aux = -g*m*pow(r.norm2(),-1.5);
	F = aux*r;
}

void Cuerpo::Move_r(double dt, double coeff){
	r+=V*(dt*coeff);
}

void Cuerpo::Move_V(double dt, double coeff){
	V+=F*(dt*coeff/m);
}

int main(){
	Cuerpo planeta;
	double m=1,r_i=10;
	double t, dt=0.1;
	double w=sqrt(g*m*pow(r_i,-3)), v_i=w*r_i, T=2*M_PI/w;

	planeta.Init(r_i,0,0,0,v_i,0,m,r_i);
	planeta.Force();
	for(t=0;t<1.1*T;t+=dt){
		std::cout << planeta.Getx() << "\t" << planeta.Gety() << std::endl;
		planeta.Move_r(dt,Theta_half);
		planeta.Force();
		planeta.Move_V(dt,Theta);
		planeta.Move_r(dt,UmTheta_half);
		planeta.Force();
		planeta.Move_V(dt,Um2Theta);
		planeta.Move_r(dt,UmTheta_half);
		planeta.Force();
		planeta.Move_V(dt,Theta);
		planeta.Move_r(dt,Theta_half);
}
	return 0;
}
