/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the tridiagonal Poisson operator

  int m_la = *la;
  int m_lab = *lab;
  int m_kv = *kv;
  
  for(size_t i=0;i<m_la*m_lab;i++){
    AB[i] = 0.;
  }

  for(size_t i=0;i<m_la*m_lab-2;i+=4){
    AB[i] = 0.;
    AB[i+1] = -1.;
    AB[i+2] = 2.;
    AB[i+3] = -1.;
  }

  AB[1] = 0.;
  AB[m_la*m_lab-1] = 0.;

  for(size_t i=0;i<m_la*m_lab;i++){
      //printf("%f ", AB[i]);
  }
  //printf("\n");


}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the identity matrix
  // Only the main diagonal should have 1, all other entries are 0

  int m_la = *la;
  int m_lab = *lab;
  int m_kv = *kv;
  
  for(size_t i=1;i<(m_lab*m_la)-1;i+=4){
    AB[i] = 0.;
    AB[i+1] = 1.;
    AB[i+2] = 0.;
  }


}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // TODO: Compute RHS vector
  RHS[0] += *BC0;
  RHS[(*la)-1] += *BC1;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D

  set_dense_RHS_DBC_1D(EX_SOL,la,BC0,BC1); // initialise the source vector
  for(int i =0;i<(*la);i++){
    EX_SOL[i] = (*BC0) +X[i]*((*BC1)-(*BC0)); // T(x) = T0 + x*(T1 − T0)
  }

}  

void set_grid_points_1D(double* x, int* la){
  // TODO: Generate uniformly spaced grid points in [0,1]

  float k = 1./((*la)+1); // calcule le séparement des points sur le grid
  float x0 = k;
  //x[0] = 0.;
  for(size_t i=0; i<(*la); i++){
    
    x[i] = x0;
    x0 += k;
  }
  //x[(*la)-1] = 1.;
}



double relative_forward_error(double* x, double* y, int* la){
  // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)
  int incx = 1;
  double nx = dnrm2_(la,x,&incx); // ||x||
  double alpha = -1.;
  daxpy_(la,&alpha,y,&incx,x,&incx); // x = alpha*y + x
  double n_x_minus_y = dnrm2_(la,x,&incx); // ||x-y|| 

  double relative_forward = n_x_minus_y / nx;
  return relative_forward;
}

int indexABCol(int i, int j, int *lab){
  // TODO: Return the correct index formula for column-major band storage
  //printf("%d\n",(i-j)+((*lab))*j +1);
  return (i-j)+(*lab)*j +2;
  //return (i-j)+((*lab))*j +1;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices

  double e = 1e-12;
  int _n = *n;
  for(size_t i =1;i<_n;i++){

    double mi = AB[indexABCol(i,i-1,lab)] / AB[indexABCol(i-1,i-1,lab)];
    AB[indexABCol(i,i-1,lab)] = mi; // mi -> multiplicateur

    AB[indexABCol(i,i,lab)] = AB[indexABCol(i,i,lab)] - AB[indexABCol(i-1,i,lab)]*mi;


  }
  //int* _info = 0;
  //info = *_info;
  return *info;
}
