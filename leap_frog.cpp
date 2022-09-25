#include <iostream>
#include <cmath>
#include "../Vector3D/vector.h"

double g=9.8;

class Cuerpo{
private:
	vector3D r,V,F;
	double m,R;
public:
	void Init(double x0,double y0, double z0,
			  double Vx0,double Vy0, double Vz0,
			  double m0,double R0);
	void start(double dt);
	void Force();
	void Move(double dt);
	double Getx(){return r.x();}
	double Gety(){return r.y();}
	double Getz(){return r.z();}
};
void Cuerpo::Init(double x0,double y0, double z0,double Vx0,double Vy0, double Vz0,double m0,double R0){
	r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0);
	m=m0; R=R0;
}

void Cuerpo::start(double dt){
	V-=F*(dt*(2*m));
}

void Cuerpo::Force(){
	double aux = -g*m*pow(r.norm2(),-1.5);
	F = aux*r;
}

void Cuerpo::Move(double dt){
	V+=F*(dt/m);
	r+=V*dt;
}

int main(){
	Cuerpo planeta;
	double m=1,r_i=10;
	double t, dt=0.01;
	double w=sqrt(g*m*pow(r_i,-3)), v_i=w*r_i, T=2*M_PI/w;

	planeta.Init(r_i,0,0,0,v_i,0,m,r_i);
	planeta.Force();
	planeta.start(dt);
	for(t=0;t<1.1*T;t+=dt){
		std::cout << planeta.Getx() << "\t" << planeta.Gety() << std::endl;
		planeta.Force();
		planeta.Move(dt);
}

	return 0;
}
