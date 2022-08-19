#include <iostream>
#include <cmath>

double g=9.8;

class Cuerpo{
private:
	double x,y,Vx,Vy,Fx,Fy,m,R;
public:
	void Init(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
	void Force();
	void Move(double dt);
	double Getx(){return x;}
	double Gety(){return y;}
};
void Cuerpo::Init(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
	x=x0; y=y0; Vx=Vx0; Vy=Vy0; m=m0; R=R0;
}

void Cuerpo::Force(){
	double aux = -g*m*pow(x*x + y*y,-1.5);
	Fx=aux*x; Fy=aux*y;
}

void Cuerpo::Move(double dt){
	x+=Vx*dt; y+=Vy*dt;
	Vx+=Fx/m*dt; Vy+=Fy/m*dt;
}

int main(){
	Cuerpo planeta;
	double m=1,r_i=10;
	double t, dt=0.0001;
	double w=sqrt(g*m*pow(r_i,-3)), v_i=w*r_i, T=2*M_PI/w;

	planeta.Init(r_i,0,0,v_i,m,r_i);
	for(t=0;t<1.1*T;t+=dt){
		std::cout << planeta.Getx() << "\t" << planeta.Gety() << std::endl;
		planeta.Force();
		planeta.Move(dt);
}

	return 0;
}
