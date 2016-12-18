// ATTENTION version bugge F. Hecht 17 oct 2016 
//  a corriger erreur:
//  MBP-FH3:C-5 hecht$ ./lap2d
// Assertion failed: (!(gCg != gCg)), function GradientConjugue, file GC.c, line 112.
//  et   GradienConjugue => GradientConjugue
// 
#include "GC.hpp"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include "gnuplot.hpp"
#include <iostream>
using namespace std; 
void PrintVect(const char * str, int n, double *x,double h)
{
    cout << str << endl; ;
    for (int i=0; i<n; ++i)
        cout << i*h << " "  <<x[i]<<endl;;
    cout << endl;
}
class MatLap2d : public MatVirt {
public:
    int N;
    int M;
	MatLap2d(int NN,int MM) : MatVirt((NN+1)*(MM+1)) ,N(NN), M(MM) {}
	double * MatMul (double *au, double *u);
   };
   
double *MatLap2d::MatMul(double *u,double *au)
   { // a matrice est pleine
       assert( ((N+1)*(M+1)==n ) && (n == m) ); 
    	for( int i=0; i<n ; ++i)
   		au[i]=0.; 
		int N1= N+1;
       double  hx = 1./N;
       double  hy = 1./M;
       double  c = 2./(hx*hx)+2./(hy*hy) ;// I,J: i-1,j;i,j  i+1,j;i,j
       double  a = -1./(hx*hx);//I,J  i,j-1;i,j   i,j+1;i,j
       double  b = -1./(hy*hy); //I,J:  ij,ij
       // I = numerotion(i,j)
       for( int i=1; i< N; ++i)
           for( int j=1; j< M; ++j)
           {
               int k = i+j*N1;
               au[k] =   u[k]*c
                      + (u[k-1]+u[k+1])*a
                      + (u[k-N1]+u[k+N1])*b;
           }
        //    au[i] =  u[i]*a + b *(u[i-1]+u[i+1]) + c ;
               ;// afaire
       return au;
   }

double ue(double x,double y)
{
    return sin(M_PI*y)*(1-x)*(x)*4.;
}
double f(double x,double y)
{
    return 4*M_PI*M_PI*sin(M_PI*y)*(1-x)*(x)+8*sin(M_PI*y);
}

int main(int argc, const char ** argv)
{
    int N=10,M=10;
    if(argc>1) N=atoi(argv[1]);
    if(argc>2) M=atoi(argv[2]);
    int n = (N+1)*(M+1);
    double hx = 1./N;
    double hy = 1./M;

    MatLap2d A(N,M);
	    MatId Id(n );
    double * u = new double[n];
    double * b = new double[n];
    double * ux = new double[n];
    
    for(int i=0; i< n; ++i)
        u[i]=0;
    
    for(int j=0,k=0; j<=M ; ++j)
        for(int i=0; i<=N ; ++i)
          {
              int kij=k++; //numerotation(i,j) row major
              ux[kij] = ue(hx*i,hy*j);
              if( i==0 || j==0 || i == N || j == M)
              {
                  b[kij]=0;
                  u[kij]=ux[kij];
              }
              else
                  b[kij]= f(hx*i,hy*j);
          }
    GradientConjugue(&A,&Id, b, u, n,1e-10,10);
   // PrintVect(" u : ",n,u,h);

    //char buf[256];
    //sprintf(buf,"gnuplot");
    //system(buf);
    //sprintf(buf,"set pm3d");
    //system(buf);

    gnuplot2d2("u-ue.txt",N,M,u,ux);

    

    for(int i=0; i< n; ++i)
        b[i]=u[i]-ux[i];
    gnuplot2d("err.txt",N,M,b);

 //   PrintVect(" b : ",n,b,h);
    delete [] b;
    delete [] u;
    delete [] ux;
    return 0;
}
