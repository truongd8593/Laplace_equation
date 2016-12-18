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
class MatLap3d : public MatVirt {
public:
    int N;
    int M;
	int K;
	MatLap3d(int NN,int MM,int KK) : MatVirt((NN+1)*(MM+1)*(KK+1)) ,N(NN), M(MM), K(KK) {}
	double * MatMul (double *au, double *u);
   };
   
double *MatLap3d::MatMul(double *u,double *au)
   { // a matrice est pleine
       assert( ((N+1)*(M+1)*(K+1)==n ) && (n == m) ); 
    	for( int i=0; i<n ; ++i)
   		au[i]=0.; 
	   int N1= N+1;
	   int NM1 = N1*(M+1);
       double  hx = 1./N;
       double  hy = 1./M;
       double  hz = 1./K;
       double  d = 2./(hx*hx)+2./(hy*hy)+2./(hz*hz) ;// I,J: i-1,j;i,j  i+1,j;i,j
       double  a = -1./(hx*hx);//I,J  i,j-1;i,j   i,j+1;i,j
       double  b = -1./(hy*hy); //I,J:  ij,ij
       double  c = -1./(hz*hz); 
       // I = numerotion(i,j)
	    // ijk = i + j*(N+1)+ k *( N+1)*(M+1)
       for( int k=1; k< K; ++k)
       for( int i=1; i< N; ++i)
           for( int j=1; j< M; ++j)
           {
               int l = i+j*N1+ k*NM1;
               au[l] =   u[l]*d
                      + (u[l-1]+u[l+1])*a
                      + (u[l-N1]+u[l+N1])*b
					  + (u[l-NM1]+u[l+NM1])*c;  
           }
        //    au[i] =  u[i]*a + b *(u[i-1]+u[i+1]) + c ;
               ;// afaire
       return au;
   }

double ue(double x,double y,double z)
{
    return sin(x)*sin(2*y)*sin(3*z);
}
double f(double x,double y,double z)
{

    return 14.*ue(x,y,z);
}
double NormeSup(int n, double * x, double *y)
	{
		double s=0;
		for(int i=0; i<n; ++i)
			s = std::max(s, std::abs(x[i]-y[i]));
		return s; 
	}
	
int main(int argc, const char ** argv)
{
    int N=10,M=10,K=10;
    if(argc>1) N=atoi(argv[1]);
    if(argc>2) M=atoi(argv[2]);
    if(argc>3) K=atoi(argv[3]);
    int n = (N+1)*(M+1)*(K+1);
    double hx = 1./N;
    double hy = 1./M;
    double hz = 1./K;

    MatLap3d A(N,M,K);
	    MatId Id(n );
    double * u = new double[n];
    double * b = new double[n];
    double * ux = new double[n];
    
    for(int i=0; i< n; ++i)
        u[i]=0;
    
    for(int k=0,l=0; k<=K ; ++k)
     for(int j=0; j<=M ; ++j)
        for(int i=0; i<=N ; ++i)
          {
              int ijk=l++; //numerotation(i,j) row major
			  // ijk = i + j*(N+1)+ k *(N+1)*(M+1)
              ux[ijk] = ue(hx*i,hy*j,hz*k);
              if( i==0 || j==0 || i == N || j == M || k==0 || k==K)
              {
                  b[ijk]=0;
                  u[ijk]=ux[ijk];
              }
              else
                  b[ijk]= f(hx*i,hy*j,hz*k);
          }
    GradientConjugue(&A,&Id, b, u, n,1e-10,10);
	cout << " Norme Sup" << NormeSup(n,u,ux) << endl; 
    delete [] b;
    delete [] u;
    delete [] ux;
    return 0;
}
