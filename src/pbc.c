#include "pbc.h"

#include <atlas/cblas.h>


double pair_dx(const double* x1, double* x2, const double L)
{
   int i;
   // void cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);
   // Calculates y <- ax + y
   cblas_daxpy(3, -1.0, x1, 1, x2, 1);

   // Loop over each element.
   for(i = 0; i < 3; i++)
   {
      if(x2[i] < (-L/2.0))
      {
         x2[i] += L;
      }
      else if(x1[i] > (L/2.0))
      {
         x2[i] -= L;
      }
   }

   return cblas_dnrm2(3, x2, 1);
}
