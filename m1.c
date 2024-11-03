#define USE_FC_LEN_T

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include "R_ext/RS.h"
#include "R_ext/Applic.h"
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>        /* for Lapack (dpotrf, etc.) and BLAS */
#include <R_ext/Linpack.h>

//#ifndef FCONE
//# define FCONE
//#endif

double **create_matrixd (int row, int col)
{
    register int i;
    double **matrix;
    
    matrix = Calloc(row, double *);
    assert(matrix);
    for(i = 0; i < row; ++i) {
        matrix[i] = Calloc(col, double);
        assert(matrix[i]);
    }
    return matrix;
}

double *create_vectord (int dim)
{
    double *vector;
    
    vector = Calloc (dim, double);
    assert(vector);
    
    return vector;
}


int *create_vector (int dim)
{
    int *vector;
    
    vector = Calloc (dim, int);
    assert(vector);
    
    return vector;
}


int **create_matrix (int row, int col)
{
    register int i;
    int **matrix;
    
    matrix = Calloc(row, int *);
    assert(matrix);
    for(i = 0; i < row; ++i) {
        matrix[i] = Calloc(col,int);
        assert(matrix[i]);
    }
    return matrix;
    
}

void matvector (double *vector, double **matrix,int nrow, int ncol)
{
    int i, j;
    
    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) vector[j*nrow+i]=matrix[i][j];
    }
    
}

void matsclpr (double **matrix, double a, double **b, int rows, int cols)
{
    int i, j;
    
   	for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            matrix[i][j] = a * b[i][j];
        }
   	}
}

void matsum (double **matrix, double **a, double **b, int rows, int cols)
{
    int i, j;
    
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) matrix[i][j] = a[i][j] + b[i][j];
    }
}


void matdiff (double **matrix, double **a, double **b, int rows, int cols)
{
    int i, j;
    
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) matrix[i][j] = a[i][j] - b[i][j];
    }
}


void freematd (double **a, int row, int col)
{
    int i;
    
    for (i = 0; i < row; i++) Free(a[i]);
    Free(a);
}

void zero_mat(double **matrix, int rows, int cols)
{
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++) 	matrix[i][j]=0;
    }
    
}
void ident_mat(double **matrix, int dim)
{
    int i, j;
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        { if (i!=j) matrix[i][j]=0;
	       else matrix[i][j]=1;
        }
    }
    
}

void vectormat (double **matrix, double *vector, int nrow, int ncol)
{
    int i, j;
    
    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) matrix[i][j] = vector[j*nrow+i];
    }
    
}


double *inverse (double *invertmat, int dim, double **matrix)
{
    
    double *qraux, *b,  *y, **ymat, *temp1, *temp2;
    double *tol, tol1;
    int *pivot, *n, *k, rank, *inf, l;
    n=&dim;
    pivot=create_vector(dim);
    qraux=create_vectord(dim);
    b=create_vectord(dim*2);
    k=&rank;
    tol1=1.0e-07;
    tol=&tol1;
    inf=&l;
    
    temp1=create_vectord(dim*dim);
    matvector (temp1, matrix,dim, dim);
    
    ymat=create_matrixd(dim, dim);
    ident_mat(ymat, dim);
    y=create_vectord(dim*dim);
    matvector (y, ymat,dim, dim);
    
    temp2=create_vectord(dim*dim);
    
    F77_CALL(dqrdc2)(temp1,n,n,n,tol,k,qraux,pivot,b);
    
   
    
    F77_CALL(dqrcf)(temp1,n,k,qraux,y, n, temp2,inf);
    Free(pivot);
    Free(qraux);
    Free(b);
    Free(temp1);
    Free(y);
   
    return temp2;
    
}


void transpose (double **matrix, double **x, int rows, int cols)
{
    int  i, j;
    for (i = 0; i < rows; i++) {
        for ( j = 0; j < cols; j++) {
            matrix[j][i] = x[i][j];
        }
    }
}

void matvecpr (double *vector, double **a, int rows, int cols, double *x)
{
    int i, j;
    double sum;
    
    for (i = 0; i < rows; i++) {
        sum = 0;
        for (j = 0; j < cols; j++) {
            sum += a[i][j] * x[j];
        }
        vector[i] = sum;
    }
}

void matpr (double **c, double **a, double **b, int m, int p, int n)
{
    int i, j, k;
    double sum;
    
    for (i = 0; i < m ; i++) {
        for(j = 0; j < n; j++) {
            sum = 0 ;
            for (k = 0; k < p; k++) {
                sum += a[i][k] * b[k][j];
            }
            c[i][j]  = sum;
        }
    }
}

void multinormal (double *vector, double *mean, double **disp, int dim)
{
    int i, *inf, l, *p;
    double *x, *z, **A, *temp1, *temp2, **tA;
    
    temp1=create_vectord(dim*dim);
    temp2=create_vectord(dim*dim);
    p=&dim; inf=&l;
    matvector(temp1, disp,dim, dim);
    F77_CALL(dpotrf)('U',p,temp1,p,inf FCONE);
    //LAPACK_dpotrf("U",p,temp1,p,0);
    //F77_CALL(chol)(temp1, p, p,temp2,inf);
    if(inf[0]!=0)
    {
        REprintf("This is not P.d matrix\n");
        Free (temp1);
        Free (temp2);
        return;
    }
    x = create_vectord (dim);
    z = create_vectord (dim);
    A=create_matrixd (dim, dim);
    tA=create_matrixd (dim, dim);
    vectormat (A, temp1, dim,dim);
    transpose (tA, A, dim, dim);
    for (i = 0; i < dim; i++) z[i] = norm_rand( );
    matvecpr (x, tA, dim, dim, z);
    for (i = 0; i < dim; i++) vector[i] = mean[i] + x[i];
    Free (x);
    Free (z);
    Free (temp1);
    Free (temp2);
    freematd (A, dim, dim);
    freematd (tA, dim, dim);
}

double *std_rWishart_factor(double nu, int p, int upper, double ans[])
{
    int pp1 = p + 1;
    memset(ans, 0, p * p * sizeof(double));
    
    for (int j = 0; j < p; j++) {	/* jth column */
        ans[j * pp1] = sqrt(rchisq(nu - (double) j));
        for (int i = 0; i < j; i++) {
            int uind = i + j * p, /* upper triangle index */
            lind = j + i * p; /* lower triangle index */
            ans[(upper ? uind : lind)] = norm_rand();
            ans[(upper ? lind : uind)] = 0;
        }
    }
    return ans;
}

void rwishart(double *ansj, double nuP, double *scCp, const int *dims){
    int info,psqr;
    double *tmp,nu=nuP, one = 1.0,zero = 0.0;
    // allocate early to avoid memory leaks in Callocs below.
    psqr = dims[0]*dims[0];
    tmp = Calloc(psqr, double);
    scCp = Calloc(psqr, double);
    F77_CALL(dpotrf)('U',&(dims[0]),scCp,&(dims[0]),&info FCONE);
    std_rWishart_factor(nu, dims, 1, tmp);
    F77_CALL(dtrmm)('R', 'U', 'N', 'N',&(dims[0]),&(dims[0]),&one, 
             scCp, &(dims[0]), tmp,&(dims[0]) FCONE FCONE FCONE FCONE);
    F77_CALL(dsyrk)('U', 'T', &(dims[1]), &(dims[1]),
                   &one, tmp, &(dims[1]),
                 &zero, ansj, &(dims[1]) FCONE FCONE);
    
    for (int i = 1; i < dims[0]; i++){
        for (int k = 0; k < i; k++){
            ansj[i + k * dims[0]] = ansj[k + i * dims[0]];
        }
    }
    Free(tmp);
}


/*********beta0 *****************/

double target_b0 (double *offset,double x, double b1, double b2,  double *phi1, 
                  double sig, int *Cores, int *Stratums,double *year, double *x1, 
                  double *x2, double *y, int *s, int N)
{
    double sum=0.0,mu=0.0,lam=0.0,a= 0.0;
    for (int i = 0; i <N; i++)
    {
        mu=offset[Stratums[i]]+x + b1*x1[Stratums[i]]+b2*x2[Stratums[i]]+phi1[Cores[i]];
        lam=exp(mu-offset[Stratums[i]]);
        sum=sum+y[i]*(mu-offset[Stratums[i]])-log(lam);
    }
       return sum;
}




double target_b1 (double *offset,double b0, double x, double b2,  double *phi1, 
                  double sig, int *Cores, int *Stratums,double *year, double *x1, 
                  double *x2, double *y, int *s, int N)
{
    double sum=0.0,mu,lam,a= 0.0;
    for (int i = 0; i < N; i++)
    {
      mu=offset[Stratums[i]]+b0 + x*x1[Stratums[i]]+b2*x2[Stratums[i]]+phi1[Cores[i]];
      lam=exp(mu-offset[Stratums[i]]);
      sum=sum+y[i]*(mu-offset[Stratums[i]])-log(lam);
    }
        return sum;
}


double target_b2 (double *offset,double b0, double b1, double x,  double *phi1, 
                  double sig, int *Cores, int *Stratums,double *year, double *x1, 
                  double *x2, double *y, int *s, int N)
{
    double sum=0.0,mu,lam,a= 0.0;
    for (int i = 0; i < N; i++){

      mu=offset[Stratums[i]]+b0 + b1*x1[Stratums[i]]+x*x2[Stratums[i]]+phi1[Cores[i]];
      lam=exp(mu-offset[Stratums[i]]);
      sum=sum+y[i]*(mu-offset[Stratums[i]])-log(lam);
    }
    return sum;
}



/*********for phi *****************/

double target_phi(double *offset,double b0, double b1, double b2, double *phi1, 
                  double sig,int *Cores,int *Stratums, double *year, double *x1, 
                  double *x2, double *y, int *s, int N, 
                  int no_regions, double **RTR11P1, 
                  double *lamda, double x1_new, int k)
{
    double sum = 0.0, func=0.0,mu=0.0,a=0.0,lam=0.0;
    int i ,j, m;
    
    double *temp1;
    temp1=create_vectord(no_regions);
    
    for(m=0; m<no_regions; m++)
    {
      if(m==k) {
        temp1[m]=x1_new;
       }
    else
    {
        temp1[m]=phi1[m];
     }
    }

    
    for (i = 0; i < N; i++)
    {
        mu=offset[Stratums[i]]+b0 + b1*x1[Stratums[i]]+b2*x2[Stratums[i]]+temp1[Cores[i]];
      lam=exp(mu-offset[Stratums[i]]);
      sum=sum+y[i]*(mu-offset[Stratums[i]])-log(lam);
    }
    
    
    func=func+lamda[0]* RTR11P1[k][k]*temp1[k]*temp1[k];
    
    for (j=0; j<no_regions; j++)
    { if (j != k)
		  {
              func=func+ 2*RTR11P1[j][k]*temp1[j]*temp1[k]*lamda[0];
          }
    }
    sum=sum-0.5*func;
    Free(temp1);
    return (sum);
}

double target_sig (double *offset,double b0, double b1, double b2, double *phi1, 
                   double x, int *Cores,int *Stratums, double *year,double *x1, 
                   double *x2,double *y, int *s, int N)
{
    double sum=0.0,mu,lam,a= 0.0;
     a=1/(x);
  
    for (int i = 0; i < N; i++)
    {
      mu=offset[Stratums[i]]+b0 + b1*x1[Stratums[i]]+b2*x2[Stratums[i]]+phi1[Cores[i]];
      lam=exp(mu-offset[Stratums[i]]);
      sum=sum+y[i]*(mu-offset[Stratums[i]])-log(lam);
    }
    
    sum=sum+dgamma(x,0.001,1000,1);
    
    return sum;
}


/***** for lamda matrix ***********************/
double target_lamda( double *phi1, int no_regions, double **RTR11P1)
{
    double sum=0.0;
    for(int i=0; i<no_regions; i++)
    {
        sum=sum+phi1[i]*phi1[i];
        for (int j=0; j<i; j++)
        {
        sum=sum+ 2*phi1[i]*phi1[j]*RTR11P1[i][j];
        }
        
    }
    return(sum);
}


/********* Metroplis-Hasting for rho and beta0  ***********************/

double metroplis_b0(double *offset,double x1_prev, double b1,double b2,
                    double *phi1, double sig,  int *Cores,int *Stratums, 
                    double *year,double *x1, double *x2,double *y, 
                    int *s, int N )
{
    double lgratio=0, u ,y1, temp1=0, temp2=0,newx;
    u=runif(0,1);

    y1=x1_prev+0.05*rnorm(0,1);
  
    temp1=target_b0(offset,y1,b1,b2,phi1, sig,Cores,Stratums,year,x1,x2,y,s,N);
    temp2=target_b0(offset,x1_prev,b1,b2,phi1,sig,Cores,Stratums,year,x1,x2,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  { newx=y1; }
        else { newx=x1_prev;}
    }
    else {newx=x1_prev;}
    return (newx);
}


double metroplis_b1(double *offset,double b0, double x_prev,double b2,double *phi1, 
                    double sig, int *Cores, int *Stratums,double *year,double *x1, 
                    double *x2,double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.05*rnorm(0,1);
  
    temp1=target_b1(offset,b0,y_new,b2,phi1,sig,Cores,Stratums,year,x1,x2,y,s,N);
    temp2=target_b1(offset,b0,x_prev,b2,phi1,sig,Cores,Stratums,year,x1,x2,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}


double metroplis_b2(double *offset, double b0, double b1, double x_prev,double *phi1,
                     double sig, int *Cores,int *Stratums, double *year,double *x1, 
                     double *x2,double *y, int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+0.05*rnorm(0,1);
    temp1=target_b2(offset,b0,b1,y_new,phi1,sig,Cores,Stratums, year,x1,x2,y,s,N);
    temp2=target_b2(offset,b0,b1,x_prev,phi1,sig,Cores,Stratums, year,x1,x2,y,s,N);
    lgratio=temp1-temp2;
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    return(x);
}


double *center(double *x,int n_region){
    double sum=0.0, avg;
    for (int i=0;i<n_region;i++){
        sum=sum+x[i];
    }
    avg=sum/n_region;
    for (int j=0;j<n_region;j++){
        x[j]=x[j]-avg;
    }
    return (x);
}

/********* Metroplis-Hasting for random effects phi ***********************/

void metroplis_phi(double *offset,double b0, double b1,double b2,double *phi_curr1, 
                   double sig, double *lamda, int *Cores, int *Stratums,
                   double *year, double *x1, double *x2, 
                   double *y,int *s, int N, int no_regions, 
                   double **RTR11P1,double x1_prev,  int j)
{
    double lgratio=0, u ,temp1=0, temp2=0;
    double x1_curr;
   
    
    u=runif(0,1);
    //if(j==0) printf("x1_prev is %f \n",x1_prev);
   // if(j==0) printf("x2_prev is %f \n",x2_prev);
    
    x1_curr=x1_prev+rnorm(0,0.003);
   
    //if(j==0) printf("x1_curr is %f \n", x1_curr);
   // if(j==0) printf("x2_curr is %f \n", x2_curr);
    
    temp1=target_phi(offset,b0,b1,b2,phi_curr1,sig, Cores,Stratums, year, x1,x2,y,s,N, no_regions,RTR11P1,lamda, x1_curr,j);
    temp2=target_phi(offset,b0,b1,b2,phi_curr1,sig, Cores, Stratums,year, x1,x2,y,s,N, no_regions,RTR11P1,lamda, x1_prev,j);
   lgratio=temp1-temp2;
    // if(j==0) printf("lgratio is %f \n", lgratio);
    if (u>0)
    { if (lgratio>=log(u))
    {
        phi_curr1[j]=x1_curr;
    }
    else
    {
        phi_curr1[j]=x1_prev;
    }
    }
    else
        
    {
        phi_curr1[j]=x1_prev;
    }
    
}



double metroplis_sig(double *offset, double b0, double b1, double b2, double *phi_curr1, 
                      double x_prev, int *Cores,int *Stratums, double *year,
                      double *x1, double *x2, double *y,int *s, int N)
{
    double lgratio=0, u ,y_new, x=0, temp1=0, temp2=0;
    u=runif(0,1);
    y_new=x_prev+rnorm(0,0.1);
    if (y_new>0){
    temp1=target_sig(offset,b0,b1,b2,phi_curr1,y_new,Cores,Stratums,year,x1,x2,y,s,N);
   
    temp2=target_sig(offset,b0,b1,b2,phi_curr1,x_prev,Cores,Stratums,year,x1,x2,y,s,N);
    
    lgratio=temp1-temp2;
    
    if (u>0)
    {
        if (lgratio>=log(u))  x=y_new ;
        else x=x_prev;
    }
    else x=x_prev;
    }
    else x=x_prev;
    
    
    
    return(x);
}


double LL (double *offset,double b0, double b1, double b2, double *phi1,  
           double sig, int *Cores,int *Stratums, double *year, double *x1, 
           double *x2, double *y, int *s, int N)
{
    double sum=0.0,mu,lam,a= 0.0;
    for (int i = 0; i < N; i++){
      mu=offset[Stratums[i]]+b0 + b1*x1[Stratums[i]]+b2*x2[Stratums[i]]+phi1[Cores[i]];
      lam=exp(mu-offset[Stratums[i]]);
      sum=sum+y[i]*(mu-offset[Stratums[i]])-log(lam);
    }
    return -2*sum;
}


void weibul(double *offset,double *b0, double *b1, double *b2, double *sig, 
            double *lamda, double *phi1,  double *vecD, 
            double *vecC,double *y, int *s, int *N_p, double *x1, 
            double *x2, double *year, int *Cores,int *Stratums, int *nend, 
            int *no_regions_p,int *dims,double *alpha_gamma,
            double *beta_gamma,double *sumLnF,double *phi_sample,int *n_burn){
	int i, j;
    int p=dims[0];
    int N=N_p[0];
    int no_regions=no_regions_p[0];
    double **D,**C,**RTR11P1;
    double x1_prev;
    double lamda_region=0.0;
	
    D=create_matrixd(no_regions,no_regions);
    C=create_matrixd(no_regions,no_regions);
    RTR11P1=create_matrixd(no_regions,no_regions);
	

    vectormat(D,vecD,no_regions,no_regions);//transfer diagnol matrix
    vectormat(C,vecC,no_regions,no_regions);//transfer adjacent matrix
    
    matdiff(RTR11P1, D, C, no_regions, no_regions);//D_w-C

	 GetRNGstate( );
    
   for (i = 1; i < nend[0]; i++)
     //for (i = 1; i < 2; i++)
	{
     //printf("the i th sample is %d \n",i);
	 b0[i]=metroplis_b0(offset,b0[i-1],b1[i-1],b2[i-1],phi1,sig[i-1],Cores,Stratums,year,x1,x2,y,s,N);
     b1[i]=metroplis_b1(offset,b0[i],b1[i-1],b2[i-1],phi1,sig[i-1],Cores,Stratums,year,x1,x2,y,s,N);
     b2[i]=metroplis_b2(offset,b0[i],b1[i],b2[i-1],phi1,sig[i-1],Cores,Stratums,year,x1,x2,y,s,N);
        for (j = 0; j < no_regions; j++)
	    {
            x1_prev=phi1[j];
           
            metroplis_phi(offset,b0[i],b1[i],b2[i],phi1,sig[i-1],lamda, Cores,Stratums,year,x1,x2,y,s,N,no_regions,RTR11P1,x1_prev,j);
            phi1=center(phi1,no_regions);
            
        }
    //printf("Before center %d is %f %f %f \n",i+1,phi1[0],phi1[1],phi1[2]);
    //printf("Before center %d is %f %f %f \n",i+1,phi2[0],phi2[1],phi2[2]);
     
      //sig[i]=metroplis_sig(offset,b0[i],b1[i],b2[i],phi1,sig[i-1],Cores,Stratums,year,x1,x2,y,s,N);

      lamda_region=target_lamda(phi1,no_regions,RTR11P1);
      lamda[0]=rgamma(0.001+no_regions/2,1/(0.001+1/2*lamda_region));
      
      if(i>=n_burn[0]) sumLnF[0]=sumLnF[0]+LL(offset,b0[i],b1[i],b2[i],phi1,sig[i],
         Cores,Stratums,year,x1,x2,y,s,N);
      for(int k=0;k<no_regions;k++){
          phi_sample[k*nend[0]+i]=phi1[k];
          
        }
    }
    
   
   
    freematd(RTR11P1,no_regions,no_regions);
    freematd(C, no_regions, no_regions);
    freematd(D, no_regions, no_regions);
   
    PutRNGstate( );

}






