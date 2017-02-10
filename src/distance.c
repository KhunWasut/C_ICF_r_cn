#include "distance.h"
#include "chemmatrixaux.h"
#include "pbc.h"
#include "utils.h"

#include <stdlib.h>


double dr_dxi(const double* x_a, const double* x_b, const int a, const int b, const int i, const double L, const double r_ab)
{
   return ((d_dx_pbc(a, b, i, i)*dx_pbc(x_a, x_b, i, L))/r_ab);
}


double d2r_dxjxi(const double* x_a, const double* x_b, const int a, const int b, const int i, const int j, const double L, const double r_ab)
{
   double factor, hess_firstterm, hess_secondterm;

   factor = d_dx_pbc(a, b, i, i)/(r_ab*r_ab);
   hess_firstterm = r_ab * d_dx_pbc(a, b, i, j);
   hess_secondterm = dx_pbc(x_a, x_b, i, L) * dr_dxi(x_a, x_b, a, b, j, L, r_ab);

   return factor*(hess_firstterm - hess_secondterm);
}


double* r_grad_x(const double* X_m, const int size_X_m, const int atom_a_index, const int atom_b_index, const double L)
{
   double *grad_vec, *x_a, *x_b, r_ab;
   double *tmp_x_a, *tmp_x_b;
   int i;
   
   // Allocation
   grad_vec = (double*)malloc(sizeof(double)*size_X_m);
   x_a = (double*)malloc(sizeof(double)*3);
   x_b = (double*)malloc(sizeof(double)*3);
  
   // Populate x_a and x_b
   for(i = 0; i < 3; i ++)
   {
      x_a[i] = X_m[atomindex_to_vecindex(atom_a_index, i)];
      x_b[i] = X_m[atomindex_to_vecindex(atom_b_index, i)];
   }

   // pair_dx implemented, but x_b value changes! Use a deep copy of x_b to calculate instead
   tmp_x_a = deep_copy(x_a, 3);
   tmp_x_b = deep_copy(x_b, 3);
   r_ab = pair_dx(tmp_x_a, tmp_x_b, L);
   free(tmp_x_a);
   free(tmp_x_b);

   for(i = 0; i < size_X_m; i++)
      grad_vec[i] = dr_dxi(x_a, x_b, atom_a_index, atom_b_index, i, L, r_ab);

   free(x_a);
   free(x_b);
   return grad_vec;
}


double* r_hess_x_j(const double* X_m, const int size_X_m, const int atom_a_index, const int atom_b_index, const int vec_index_j, const double L)
{
   double *hess_j_vec, *x_a, *x_b, *tmp_x_a, *tmp_x_b, r_ab;
   int i;

   // Allocation
   hess_j_vec = (double*)malloc(sizeof(double)*size_X_m);
   x_a = (double*)malloc(sizeof(double)*3);
   x_b = (double*)malloc(sizeof(double)*3);

   // Populate x_a and x_b
   for(i = 0; i < 3; i++)
   {
      x_a[i] = X_m[atomindex_to_vecindex(atom_a_index, i)];
      x_b[i] = X_m[atomindex_to_vecindex(atom_b_index, i)];
   }

   tmp_x_a = deep_copy(x_a, 3);
   tmp_x_b = deep_copy(x_b, 3);
   r_ab = pair_dx(tmp_x_a, tmp_x_b, L);
   free(tmp_x_a);
   free(tmp_x_b);

   for(i = 0; i < size_X_m; i++)
      hess_j_vec[i] = d2r_dxjxi(x_a, x_b, atom_a_index, atom_b_index, i, vec_index_j, L, r_ab);

   free(x_a);
   free(x_b);
   return hess_j_vec;
}
