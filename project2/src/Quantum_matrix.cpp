//  Solves linear equations for simple tridiagonal matrix using LU decomposition
//  This is armadillo version that calls the function solve. 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include <cstdlib>
#include "time.h"
// use namespace for output and input
using namespace std;
using namespace arma;

// the offdiag function, using Armadillo
double offdiag(mat &A, int &p, int &q, int n)
{
double max;
for (int i = 0; i < n; ++i)
{
for ( int j = i+1; j < n; ++j)
{
double aij = fabs(A(i,j));
if ( aij > max)
{
max = aij; p = i; q = j;
}
}
}
return max;
}
// more statements

void Jacobi_rotate ( mat &A, mat &R, int k, int l, int n )
{
double s, c;
if ( A(k,l) != 0.0 ) {
double t, tau;
tau = (A(l,l) - A(k,k))/(2*A(k,l));
if ( tau >= 0 ) {
t = 1.0/(tau + sqrt(1.0 + tau*tau));
} else {
t = -1.0/(-tau +sqrt(1.0 + tau*tau));
}
c = 1/sqrt(1+t*t);
s = c*t;
} else {
c = 1.0;
s = 0.0;
}
double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
a_kk = A(k,k);
a_ll = A(l,l);
A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
A(k,l) = 0.0; // hard-coding non-diagonal elements by hand
A(l,k) = 0.0; // same here
for ( int i = 0; i < n; i++ ) {
if ( i != k && i != l ) {
a_ik = A(i,k);
a_il = A(i,l);
A(i,k) = c*a_ik - s*a_il;
A(k,i) = A(i,k);
A(i,l) = c*a_il + s*a_ik;
A(l,i) = A(i,l);
}
// And finally the new eigenvectors
r_ik = R(i,k);
r_il = R(i,l);
R(i,k) = c*r_ik - s*r_il;
R(i,l) = c*r_il + s*r_ik;
}
return;
} // end of function jacobi_rotate


int main(int argc, char *argv[]){
      int n = atoi(argv[1]);


      //variables
      double rho_min = 0.0;
      double rho_max = 7.0;
      double h = (rho_max-rho_min)/(n);
      double hh = h*h;
      double hh_inv = 1.0/hh;


      vec rho(n); vec V(n);
      for(int i = 0; i<n ; i++){

          rho[i] = rho_min + (i+1)*h;
          V[i] = rho[i]*rho[i];

      }

      mat A = zeros<mat>(n,n);
      for (int i=0; i< n-1; i++){
            A(i+1,i+1)= 2.0*hh_inv + V[i+1];
            A(i+1,i)  = -1.0*hh_inv;
            A(i,i)    = 2.0*hh_inv + V[i];
            A(i,i+1)  = -1.0*hh_inv;

      }

      vec eigval;
      mat eigvec;
      eig_sym(eigval,eigvec,A);

      cout << A << endl;
      cout << "Eigenvalues:  " <<endl << eigval << endl;

      mat R = eye (n,n);
      double maxnondiag = 1;
      double maxiter = 100;
      double tolerance = 1.0E-10;
      int iterations = 0;
      while ( maxnondiag > tolerance && iterations <= maxiter)
{
    int p, q;
    maxnondiag = offdiag(A, p, q, n);
    Jacobi_rotate(A, R, p, q, n);
    iterations++;

}
    cout << A << endl;

      return 0;
}
