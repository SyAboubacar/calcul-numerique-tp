/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

double pi = acos(-1.0);

void eig_poisson1D(double* eigval, int *la){
  // TODO: Compute all eigenvalues for the 1D Poisson operator
   
  int n = *la;

  for(int k=0;k<n;k++){
    eigval[k] = 2 - 2*cos((k* pi)/(n+1));
  }
}

double eigmax_poisson1D(int *la){
  // TODO: Compute and return the maximum eigenvalue for the 1D Poisson operator
  int n = *la;
  int k = n;
  double eigmax = 2 - 2*cos((k* pi)/(n+1));
  return eigmax;
}

double eigmin_poisson1D(int *la){
  // TODO: Compute and return the minimum eigenvalue for the 1D Poisson operator
  int n = *la;
  int k = 1;

  double eigmin = 2 - 2*cos((k* pi)/(n+1));
  return eigmin;
}

double richardson_alpha_opt(int *la){
  // TODO: Compute alpha_opt
  double opt = 2/(eigmin_poisson1D(la) + eigmax_poisson1D(la));
  return opt;
}

/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iteration
  // 1. Compute residual r = b - A*x (use dgbmv for matrix-vector product)
  // 2. Update x = x + alpha*r (use daxpy)
  // 3. Check convergence: ||r||_2 < tol (use dnrm2)
  // 4. Store residual norm in resvec and repeat
  
  //SUBROUTINE DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
  //y= alpha*A*x + beta*y

  double alpha = -1.;
  double beta = 1.;
  int inc = 1;
  char trans = 'N';

  
  double norm_r,norm_b,rel;

  double *residus = (double*)malloc(*la * sizeof(double));
  
  // dcopy_(la,RHS,&inc,residus,&inc);
  // dcopy_(la,X,&inc,old_x,&inc);

  // dgbmv_(&trans,la,lab,kl,ku,&alpha,AB,lab,X,&inc,&beta,residus,&inc); // r = RHS 
  // daxpy_(la,alpha_rich,residus,&inc,X,&inc); // r = alpha*r + x


  norm_b = cblas_dnrm2(*la,RHS,inc);
  // norm_r = cblas_dnrm2(*la,residus,inc);
  // rel = norm_r/norm_b;

  // resvec[0] = rel;

  *nbite = 0;
  for(int k=0;k < *maxit;k++){
    dcopy_(la, RHS, &inc, residus, &inc);

    dgbmv_(&trans,la,la,kl,ku,&alpha,AB,lab,X,&inc,&beta,residus,&inc);

    norm_r = cblas_dnrm2(*la,residus,inc);
    rel = norm_r/norm_b;
    resvec[k] = rel;

    if(rel < *tol){
      break;
    }

    daxpy_(la,alpha_rich,residus,&inc,X,&inc);

    
    printf("resvec[%d] : %lf \n",k+1,rel);
    *nbite = k+1;
  }

  free(residus);
}

/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal elements from AB and store in MB
  // MB should contain only the diagonal of A
}

/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal and lower diagonal from AB
  // MB should contain the lower triangular part (including diagonal) of A
}

/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iterative method
}

