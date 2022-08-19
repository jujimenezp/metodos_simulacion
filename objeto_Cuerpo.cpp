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
	Fx=0; Fy=-m*g;
}

void Cuerpo::Move(double dt){
	x+=Vx*dt; y+=Vy*dt;
	Vx+=Fx/m*dt; Vy+=Fy/m*dt;
}

int main(){
	Cuerpo balon;
	double t, dt=0.1;
	balon.Init(0,0,4,3,0.453,0.15);
	for(t=0;t<3;t+=dt){
		std::cout << balon.Getx() << "\t" << balon.Gety();
		balon.Force();
		balon.Move(dt);
}

	return 0;
}
