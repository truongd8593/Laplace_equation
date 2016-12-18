#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
using namespace std;
void gnuplot1d(const char * file, int N, double *x)
{
    double h = 1./N;
    int n=N+1;
	{
    ofstream  f(file);
    assert(f);
    for(int i=0; i<n ; ++i)
        f << i*h << " " <<x[i] << endl; 

    }// ici f est ferme 
    char buf[256];
    sprintf(buf,"echo 'plot \"%s\" w l'|  gnuplot -p",file);
    system(buf);
}
void gnuplot2d2(const char * file, int N,int M, double *x,double *y)
{
    double hx = 1./N;
    double hy = 1./M;
	{
    ofstream  f(file);
    assert(f);
    for (int j=0,k=0; j<=M; ++j)
    {
	  for (int i=0; i<=N; ++i,++k)
        f << i*hx<< " " << j*hy << " " <<x[k] << " " << y[k] << endl; 
        f << "\n\n"; 
    }
    
    }
    char buf[256];
    sprintf(buf,"echo 'splot \"%s\" w l , \"%s\" u 1:2:4 w l '|  gnuplot -p",file,file);
    system(buf);
}

void gnuplot2d(const char * file, int N,int M, double *x)
{
    double hx = 1./N;
    double hy = 1./M;
    FILE *f=fopen(file,"w");
    assert(f);
    for (int j=0,k=0; j<=M; ++j)
    {for (int i=0; i<=N; ++i,++k)
        fprintf(f,"  %f %f %e \n",i*hx,j*hy,x[k]);
        fprintf(f,"\n\n");
    }
    
    fclose(f);
    char buf[256];
    sprintf(buf,"echo 'splot \"%s\" w l  '|  gnuplot -p",file);
    system(buf);
}

