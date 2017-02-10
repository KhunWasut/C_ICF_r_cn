#include "coordnum.h"
#include "chemmatrixaux.h"
#include "distance.h"
#include "pbc.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>


double n_ab_k(const double r_ab_k, const double r0, const double numer_pow)
{
   return (1.0-pow((r_ab_k/r0), numer_pow));
}


double d_ab_k(const double r_ab_k, const double r0, const double denom_pow)
{
   return (1.0-pow((r_ab_k/r0), denom_pow));
}


double dN_xi(const double n, const double r0, const double* x_a, const double* x_b, const int a, const int b, const int i, const double L, 
      const double r_ab)
{
   double N;
   N = n_ab_k(r_ab, r0, n);
   return (n*(N-1.0)*(dr_dxi(x_a, x_b, a, b, i, L, r_ab))/r_ab);
}


double dD_xi(const double m, const double r0, const double* x_a, const double* x_b, const int a, const int b, const int i, const double L, 
      const double r_ab)
{
   double D;
   D = d_ab_k(r_ab, r0, m);
   return (m*(D-1.0)*(dr_dxi(x_a, x_b, a, b, i, L, r_ab))/r_ab);
}


double* cn_grad_x(const double* X_m, const int size_X_m, const int a_index, const int* b_indices, const int size_b_indices, const double n, 
      const double m, const double r0, const double L)
{
   double *x_a, *x_b, *grad_x_vec, *tmp_x_a, *tmp_x_b;
   double sum_dcn_dxi, r_ab, D_val, N_val;
   double factor, grad_firstterm, grad_secondterm;
   int i, k, x_a_i, x_b_i;

   // Allocation of x_a and grad_x_vec
   x_a = (double*)malloc(sizeof(double)*3);
   grad_x_vec = (double*)malloc(sizeof(double)*size_X_m);

   // Value assignment for x_a
   for(x_a_i = 0; x_a_i < 3; x_a_i++)
      x_a[x_a_i] = X_m[atomindex_to_vecindex(a_index, x_a_i)];

   for(i = 0; i < size_X_m; i++)
   {
      sum_dcn_dxi = 0.0;

      for(k = 0; k < size_b_indices; k++)
      {
         // Allocation of x_b
         x_b = (double*)malloc(sizeof(double)*3);

         for(x_b_i = 0; x_b_i < 3; x_b_i++)
            x_b[x_b_i] = X_m[atomindex_to_vecindex(b_indices[k], x_b_i)];

         tmp_x_a = deep_copy(x_a, 3);
         tmp_x_b = deep_copy(x_b, 3);
         r_ab = pair_dx(tmp_x_a, tmp_x_b, L);
         free(tmp_x_a);
         free(tmp_x_b);

         D_val = d_ab_k(r_ab, r0, m);
         N_val = n_ab_k(r_ab, r0, n);

         factor = pow(D_val, -2.0);
         grad_firstterm = D_val * dN_xi(n, r0, x_a, x_b, a_index, b_indices[k], i, L, r_ab);
         grad_secondterm = N_val * dD_xi(m, r0, x_a, x_b, a_index, b_indices[k], i, L, r_ab);

         free(x_b);
         sum_dcn_dxi += (factor*(grad_firstterm - grad_secondterm));
      }
      grad_x_vec[i] = sum_dcn_dxi;
   }

   free(x_a);
   return grad_x_vec;
}


double* cn_hess_x_j(const double* X_m, const int size_X_m, const int vec_index_j, const int a_index, const int* b_indices, const int size_b_indices, 
      const double n, const double m, const double r0, const double L)
{
   double *x_a, *x_b, *hess_x_j, *tmp_x_a, *tmp_x_b;
   double sum_d2cn_dxjxi, r_ab, D_val, N_val;
   double dDj, dNj, dr_i, dr_j, d2r_ji;
   double factor;
   double hess_t1_1, hess_t1_2, hess_t1_3, hess_t1_4, hess_t1;
   double hess_t2_1, hess_t2_2, hess_t2_3, hess_t2_4, hess_t2;
   double hess_t3_1, hess_t3_2, hess_t3;
   double MAINFAC_II, MAINFAC_IJ, MAINFAC_JJ;
   int i, x_a_i, j, k, x_b_i;

   j = vec_index_j;

   // Allocation of x_a and hess_x_j
   x_a = (double*)malloc(sizeof(double)*3); 
   hess_x_j = (double*)malloc(sizeof(double)*size_X_m);

   // Value assignment for x_a
   for(x_a_i = 0; x_a_i < 3; x_a_i++)
      x_a[x_a_i] = X_m[atomindex_to_vecindex(a_index, x_a_i)];

   for(i = 0; i < size_X_m; i++)
   {
      sum_d2cn_dxjxi = 0.0;
      for(k = 0; k < size_b_indices; k++)
      {
         // Allocation of x_b
         x_b = (double*)malloc(sizeof(double)*3);

         for(x_b_i = 0; x_b_i < 3; x_b_i++)
            x_b[x_b_i] = X_m[atomindex_to_vecindex(b_indices[k], x_b_i)];

         tmp_x_a = deep_copy(x_a, 3);
         tmp_x_b = deep_copy(x_b, 3);
         r_ab = pair_dx(tmp_x_a, tmp_x_b, L);
         free(tmp_x_a);
         free(tmp_x_b);

         // CALCULATION OF MAINFACS: REDUCING FUNTION CALL OVERHEAD
         MAINFAC_II = d_dx_pbc(a_index, b_indices[k], i, i);
         MAINFAC_IJ = d_dx_pbc(a_index, b_indices[k], i, j);
         MAINFAC_JJ = d_dx_pbc(a_index, b_indices[k], j, j);

         // Exploting zeros for fast calculations!
         if ((MAINFAC_II == 0.0) && (MAINFAC_IJ == 0.0) && (MAINFAC_JJ == 0.0))
            sum_d2cn_dxjxi = 0.0;
         else
         {
            D_val = d_ab_k(r_ab, r0, m);
            N_val = n_ab_k(r_ab, r0, n);

            if (MAINFAC_JJ == 0.0)
            {
               dDj = 0.0;
               dNj = 0.0;
               dr_j = 0.0;
            }
            else
            {
               dDj = dD_xi(m, r0, x_a, x_b, a_index, b_indices[k], j, L, r_ab);
               dNj = dN_xi(n, r0, x_a, x_b, a_index, b_indices[k], j, L, r_ab);
               dr_j = dr_dxi(x_a, x_b, a_index, b_indices[k], j, L, r_ab);
            }

            if (MAINFAC_II == 0.0)
               dr_i = 0.0;
            else
               dr_i = dr_dxi(x_a, x_b, a_index, b_indices[k], i, L, r_ab);

            if (((MAINFAC_II == 0.0) || (MAINFAC_IJ == 0.0)) && (MAINFAC_JJ == 0.0))
               d2r_ji = 0.0;
            else
               d2r_ji = d2r_dxjxi(x_a, x_b, a_index, b_indices[k], i, j, L, r_ab);

            factor = pow(D_val, -4.0);

            hess_t1_1 = dr_i*(n*D_val*dNj)/r_ab;
            hess_t1_2 = dr_i*(n*(N_val-1.0)*dDj)/r_ab;
            hess_t1_3 = n*D_val*(N_val-1.0)*d2r_ji/r_ab;
            hess_t1_4 = (-1.0)*n*D_val*(N_val-1.0)*dr_i*(dr_j/(r_ab*r_ab));
            hess_t1 = (D_val * D_val) * (hess_t1_1 + hess_t1_2 + hess_t1_3 + hess_t1_4);

            hess_t2_1 = dr_i*(m*N_val*dDj)/r_ab;
            hess_t2_2 = dr_i*(m*(D_val-1.0)*dNj)/r_ab;
            hess_t2_3 = m*N_val*(D_val-1.0)*d2r_ji/r_ab;
            hess_t2_4 = (-1.0)*m*N_val*(D_val-1.0)*dr_i*(dr_j/(r_ab*r_ab));
            hess_t2 = (D_val * D_val) * (hess_t2_1 + hess_t2_2 + hess_t2_3 + hess_t2_4);

            hess_t3_1 = n*D_val*(N_val - 1.0)*dr_i*dDj/r_ab;
            hess_t3_2 = m*N_val*(D_val - 1.0)*dr_i*dDj/r_ab;
            hess_t3 = 2.0*D_val*(hess_t3_1 - hess_t3_2);

            sum_d2cn_dxjxi += (factor*(hess_t1 - hess_t2 - hess_t3));
         }
         free(x_b);
      }
      hess_x_j[i] = sum_d2cn_dxjxi;
   }
   free(x_a);
   return hess_x_j;
}
