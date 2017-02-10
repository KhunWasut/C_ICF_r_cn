#include "matrix.h"

#include <atlas/cblas.h>
#include <atlas/clapack.h>
#include <stdlib.h>


// Matrix inverse calculation
// BLAS + LAPACK implementation
// NOTE: The value of A got replaced by the inverse, so deep copy before calling this function!
void matrix_invert(double* A, int N)
{
   int *IPIV, info1, info2;

   // Allocations (don't understand why an example use N+1 for IPIV
   IPIV = (int*)malloc(sizeof(int)*N);

   // We call dgetrf first to perform LU factorization
   clapack_dgetrf(CblasRowMajor, N, N, A, N, IPIV);
   clapack_dgetri(CblasRowMajor, N, A, N, IPIV);

   free(IPIV);
}


double* zeros(int size)
{
   double *arr;
   int i;
   arr = (double*)malloc(sizeof(double)*size);

   for(i = 0; i < size; i++)
      arr[i] = 0.0;

   return arr;
}
