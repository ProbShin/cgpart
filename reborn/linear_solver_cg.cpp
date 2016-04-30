
/*
 A = [2 1 0 0 0
      1 2 1 0 0 
      0 1 2 1 0
      0 0 1 2 1
      0 0 0 1 2]
 A_ia     = {0 2 5 8 11 13}
 A_ja     = {0 1 0 1 2 1 2 3 2 3 4 3 4}
 A_values = {2 1 1 2 1 1 2 1 1 2 1 1 2}
*/

#include "linear_solver.h" // parameters
///#include <omp.h>

#include <cmath>
#include <cstdio>

typedef int MKL_INT;
#define CG_MAX_ITER 1

void imqqtax(MKL_INT*, MKL_INT*, double*, double*, double*, double*);
void cg_core(MKL_INT*, MKL_INT*, double*, double*, double*,  double*, double*, double*, int, int&);


/* implements CG linear solver with OpenMp */

void linear_solver(
                   MKL_INT* A_ia,
                   MKL_INT* A_ja,
                   double*  A_values, // CSR format of full matri A
                   double*  Q1,          // n * s double, column major
                   double*  RHS,         // n * s double, column major
                   double*  solution,    // n * s double, column major
		           MKL_INT n,
		           MKL_INT s) {
  int its;
  double  *RHSCOPY = new double [n*s];   //workspace of p in cg
  double  *WKSPACE1 = new double [n*s];   //workspace of Ap in cg      each col need size== p
  double  *WKSPACE2 = new double [n*s];    //workspace of I-QQ'Ax in cg  each col need size== sizex	  
  //omp parallal for
  for(int i=0; i<n*s; i++){
    RHSCOPY[i] = RHS[i];
  }

  //omp parallal for
  for(int j=0; j<1; j++){
    cg_core(A_ia, A_ja, A_values, Q1, RHS+n*j, RHSCOPY+n*j, WKSPACE1+n*j, WKSPACE2+n*j, solution+n*j, n, its);
  }


  delete[] RHSCOPY;
  return;
}



int main(int argc, char** argv){
  MKL_INT is[]  = {0, 2, 5, 8, 11, 13};
  MKL_INT js[]  = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4 };
  double  vs[]  = {2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2 };
  double  q1[]  = {1, 1, 1, 1, 1, 1};
  double  rhs[] = {1, 1, 1, 1, 1, 1};
  double  xs[]  = {1, 1, 1, 1, 1, 1};
  MKL_INT n=3;
  MKL_INT s=2;

  linear_solver(is,js,vs,q1,rhs,xs,n,s);	
}


/* implement (I-QQ')Ax */
void imqqtax(MKL_INT* ia,
		MKL_INT* ja,
		double*  va,
		double*  Q,
		double*  X,
		double*  QQTAX;
		double*  Y,
		MKL_INT  n,
		MKL_INT  s) {

  printf("X:");
  for(int i=0; i<n; i++)
	  printf("%f ", X[i]);
  printf("\n");
  
  for(int i=0; i<n; i++){
    Y[i]=0.0;
    for(int j, jt=ia[i]; jt<ia[i+1]; jt++){
        Y[i]+=X[ ja[jt] ]*va[jt];
	}
  }

  dgemv('T',n,s,1, Q1,n, Y,  1,  0,QtAX, 1);
  
printf("wocao:");
for(int i=0; i<1; i++){
	printf("%f\n", QtAX[i]);
}

  dgemv('N',n,s,-1,Q1,n, QtAX,1, 1,Y, 1);



  printf("X:");
  for(int i=0; i<n; i++)
	  printf("%f ", X[i]);
  printf("\n");

  printf("Y:");
  for(int i=0; i<n; i++)
	  printf("%f ", Y[i]);
  printf("\n");
}


double norm2(double *x,  MKL_INT n){
	double rlt=0.0;
	//todo omp parallel
	for(int i=0; i<n; i++){
		rlt+=x[i]*x[i];
	}
	return sqrt(rlt);
}



void cg_core(MKL_INT* A_ia, 
             MKL_INT* A_ja,
             double*  A_v,
             double*  Q1,
             double*  rhs,
             double*  rhscopy,
			 double*  wkspace1,
			 double*  wkspace2,
             double*  solution,
             MKL_INT  n,
             int      &itr_out
		) {
	double *&r =rhs; 
	double *&p =rhscopy;
	double rnorm=0.0;
	double *&Ap = wkspace1;
	double *&qqtax = wkspace2;
	double alpha=0.0;
	double beta =0.0;
	double *&x = solution;

	for(int iter=1;  iter<=CG_MAX_ITER; iter++) {
      rnorm = norm2(r,n);
	  imqqtax(A_ia, A_ja, A_v, Q1, p, qqtax, Ap, n);
	

      printf("\n%d: rnorm %f\n", iter, rnorm);
	}

}


