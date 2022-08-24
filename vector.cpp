#include <iostream>
#include <cmath>
#include "vector.h"

int main(int argc, char* argv[]){
  vector3D a, b, c;

  a.load(1,0,0);
  b.load(0,1,0);
  c=a^b;
  c.show();

  return 0;
}
