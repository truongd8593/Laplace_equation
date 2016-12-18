
struct  MatVirt {
	int n,m;
	MatVirt(int nn,int mm): n(nn),m(mm) {}
	MatVirt(int nn): n(nn),m(nn) {}
    virtual double *MatMul(double *u,double *au)=0; 
	virtual ~MatVirt() {}
}  ;


int GradientConjugue(MatVirt *A, // fonction et pointeur data pour A
		     MatVirt *C, // fonction et pointeur data pour C
		     double * b, // second membre
		     double * x, // solution qui contient une initialisation
		     int nbitermax,double eps,
                     int niveauimpression)
;


double * ProduitMatVec(MatVirt *A,double *x, double *Ax);
double mysdot(int n,double *x,double *y);
double * mysaxpy(int n,double a,double *x,double *y);
double * myscal(int n,double a,double *x);

class matPleineColMajor  : public MatVirt 
{
public:
	double *A; 
	matPleineColMajor(int nn,int mm,double *AA): MatVirt(nn,mm),A(AA) {}
	
 double *MatMul(double *x,double *ax)
{ // a matrice est pleine
    for(int i=0;i< n ;++i)
        ax[i]=0.;
    for(int i=0;i< n ;++i)
        for(int j=0;j< m ;++j)
            ax[i] +=A[i*m+j]*x[j] ;
    return ax;
}

};

class matPleineRowMajor  : public MatVirt 
{
public:
	double *A; 
	matPleineRowMajor(int nn,int mm,double *AA): MatVirt(nn,mm),A(AA) {}
	
    double *MatMul( double *x,double *ax)
    { // a matrice est pleine
          for(int i=0;i< n ;++i)
           ax[i]=0.;
       for(int j=0;j< m ;++j)
           for(int i=0;i< n ;++i)
               ax[i] +=A[j*n+i]*x[j] ;
       return ax;
   }
};

class MatId  : public MatVirt 
{
public:
	
	MatId(int nn): MatVirt(nn) {}
	
double *MatMul( double *x,double *ax)
{ 
 for(int i=0; i<n; ++i)
	 ax[i]=x[i];
   return ax;
}

};