#include <iostream>
#include <fstream>
#include <cmath>

#define Lx 16
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx;

//--------- CODIGO DEL DEVICE  --------
//--------- Kernels ------------
__global__ void AddTwoVectors(float *d_a, float * d_b, float * d_c){
  int ix; ix=blockIdx.x*blockDim.x+threadIdx.x;
  d_c[ix]=d_a[ix]+d_b[ix];
}
//-------- CODIGO DEL HOST -------
int main(void){
  int ix;
  //Declarar todas las variables por duplicado
  //--------en el host
  float h_a[Lx],h_b[Lx],h_c[Lx];
  //--------en el device
  float *d_a; cudaMalloc((void**) &d_a, Lx*sizeof(float));
  float *d_b; cudaMalloc((void**) &d_b, Lx*sizeof(float));
  float *d_c; cudaMalloc((void**) &d_c, Lx*sizeof(float));

  //LLenar los datos que vamos a procesar
  for(ix=0; ix<Lx;ix++){
    h_a[ix]=ix; h_b[ix]=2*ix;}

    //Enviarlos al Device
    cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_b,h_b,Lx*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_c,h_c,Lx*sizeof(float),cudaMemcpyHostToDevice);

  //correr en el Device
  dim3 ThreadsPerBlock(Nx,0,0);
  dim3 BlocksPerGrid(Mx,0,0);
  AddTwoVectors<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,d_b,d_c);

  //Devolver el resultado al host
  cudaMemcpy(h_c,d_c,Lx*sizeof(float),cudaMemcpyDevicetoHost);

  //Imprimir los resultados
  for(ix=0;ix<Lx;ix++){
    std::cout <<h_c[ix] <<std::endl;
  }

  //Liberar la memoria dinamica
  cudaFree(d_a);
  cudaFree(d_b);
  cudafree(d_c);

  return 0;
}
