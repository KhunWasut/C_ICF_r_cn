#include "icf2d.h"
#include "utils.h"
#include "matrix.h"
#include "distance.h"
#include "coordnum.h"
#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <atlas/cblas.h>


double* icf_firstterm(const double* W, const double* G_w_inv, const double* grad_V, const int size_grad_V, const int num_cv)
{
   // Final output is 2 x 1 vector
   double *tmp1, *ans;
   tmp1 = zeros(num_cv*size_grad_V);    // D x 3N
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_cv, size_grad_V, num_cv, 1.0, G_w_inv, num_cv, W, size_grad_V, 0.0, tmp1, size_grad_V);
   ans = zeros(num_cv);
   // Producing a vector, so using 'dgemv' instead of 'dgemm'
   cblas_dgemv(CblasRowMajor, CblasNoTrans, num_cv, size_grad_V, -1.0, tmp1, size_grad_V, grad_V, 1, 0.0, ans, 1);
   free(tmp1);
   return ans;
}


double* icf_secondterm(const double* X_m, const int size_X_m, const int r_a, const int r_b, const int cn_a, 
      const int* cn_b, const int size_cn_b, const double n, const double m, const double r0, const double L, const double* G_w_inv, 
      const double* W, const double* mu_inv, const double* grad_cv_t, const int num_cv, const double kT)
{
   double *sum_div, *hess_x_i_r, *hess_x_i_cn, *hess_cv_t;
   double *tmp1, *tmp2;
   double *first_subterm;
   double *second_subterm_innermost, *second_subterm_front, *second_subterm;
   double *d_GwinvW_dxi;
   double eval_time = 0.0;
   clock_t t1, t2;
   int i, j;

   sum_div = zeros(num_cv);
   for(i = 0; i < size_X_m; i++)
   {
      t1 = clock();
      hess_x_i_r = r_hess_x_j(X_m, size_X_m, r_a, r_b, i, L);
      hess_x_i_cn = cn_hess_x_j(X_m, size_X_m, i, cn_a, cn_b, size_cn_b, n, m, r0, L);

      // Allocation for hess_cv_t
      hess_cv_t = (double*)malloc(sizeof(double)*num_cv*size_X_m);      // D x 3N

      for(j = 0; j < size_X_m; j++)
      {
         hess_cv_t[j] = hess_x_i_r[j];
         hess_cv_t[size_X_m+j] = hess_x_i_cn[j];
      }

      free(hess_x_i_r);
      free(hess_x_i_cn);

      tmp1 = zeros(num_cv*size_X_m);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_cv, size_X_m, num_cv, 1.0, G_w_inv, num_cv, hess_cv_t, size_X_m, 0.0, tmp1, size_X_m);
      first_subterm = zeros(num_cv*size_X_m);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_cv, size_X_m, size_X_m, 1.0, tmp1, size_X_m, mu_inv, size_X_m, 0.0, first_subterm, size_X_m);
      free(tmp1);

      tmp1 = zeros(num_cv*num_cv);
      tmp2 = zeros(num_cv*size_X_m);
      second_subterm_innermost = zeros(num_cv*num_cv);

      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, num_cv, num_cv, size_X_m, 1.0, W, size_X_m, hess_cv_t, size_X_m, 0.0, tmp1, num_cv);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_cv, size_X_m, size_X_m, 1.0, hess_cv_t, size_X_m, mu_inv, 
            size_X_m, 0.0, tmp2, size_X_m);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, num_cv, num_cv, size_X_m, 1.0, tmp2, size_X_m, grad_cv_t, 
            size_X_m, 0.0, second_subterm_innermost, num_cv);
      free(tmp2);
      // second_subterm_innermost <- (1.0)*tmp1 + ssi
      // In reality, this daxpy adds vector but in BLAS all matrices are 1D anyway.
      cblas_daxpy(num_cv*num_cv, 1.0, tmp1, 1, second_subterm_innermost, 1);
      free(tmp1);

      tmp1 = zeros(num_cv*num_cv);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_cv, num_cv, num_cv, -1.0, G_w_inv, num_cv, second_subterm_innermost, 
            num_cv, 0.0, tmp1, num_cv);
      free(second_subterm_innermost);
      second_subterm_front = zeros(num_cv*num_cv);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_cv, num_cv, num_cv, 1.0, tmp1, num_cv, G_w_inv, num_cv, 0.0, second_subterm_front, num_cv);
      free(tmp1);
      second_subterm = zeros(num_cv*size_X_m);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_cv, size_X_m, num_cv, 1.0, second_subterm_front, num_cv, W, 
            size_X_m, 0.0, second_subterm, size_X_m);
      free(second_subterm_front);

      d_GwinvW_dxi = deep_copy(second_subterm, num_cv*size_X_m);
      cblas_daxpy(num_cv*size_X_m, 1.0, first_subterm, 1, d_GwinvW_dxi, 1);

      sum_div[0] += d_GwinvW_dxi[i];
      sum_div[1] += d_GwinvW_dxi[i+size_X_m];

      t2 = clock();

      printf("The time it took to evaluate hessians of element %d is %.2lf\r", i+1, (double)(t2-t1)/CLOCKS_PER_SEC);

      eval_time = eval_time + ((double)(t2-t1)/CLOCKS_PER_SEC);

      free(first_subterm);
      free(second_subterm);
      free(hess_cv_t);

      free(d_GwinvW_dxi);
   }

   printf("\nTotal time to evaluate %d elements is %.2lf minutes\n", size_X_m, eval_time/60.0);
   printf("Average evaluation time per step is %.4lf seconds\n", eval_time/((double)(size_X_m)));

   cblas_dscal(num_cv, kT, sum_div, 1);

   return sum_div;
}


/* NOTE: Matrix inversion function PERMANENTLY changes the value! */
/* DEEP COPY A MATRIX BEFORE EVALUATING AN INVERSE!! */

/* icf_construct returns a D x 1 array!! */
double* icf_construct_2d(const double* X_m, const int size_X_m, const double* mu, const double* grad_V, const double kT, const int ra_ind, 
      const int rb_ind, const int cn_a_ind, const int* cn_b_inds, const int size_cn_b_inds, const double n, const double m, 
      const double r0, const double L, const int num_cv)
{
   double *grad_r, *grad_cn, *grad_cv_t;
   double *W, *G_w_inv, *mu_inv;
   double *first, *icf;
   int i;

   // grad_r and grad_cn are allocated in the functions! Free them after use!
   grad_r = r_grad_x(X_m, size_X_m, ra_ind, rb_ind, L);
   grad_cn = cn_grad_x(X_m, size_X_m, cn_a_ind, cn_b_inds, size_cn_b_inds, n, m, r0, L);

   // Allocation for grad_cv_t
   grad_cv_t = (double*)malloc(sizeof(double)*2*size_X_m);  // D x 3N

   // Populating grad_cv_t as a 1-dimensional D x 3N matrix in C representation
   for(i = 0; i < size_X_m; i++)
   {
      grad_cv_t[i] = grad_r[i];
      grad_cv_t[size_X_m+i] = grad_cn[i];
   }

   free(grad_r);
   free(grad_cn);

   mu_inv = deep_copy(mu, size_X_m*size_X_m);
   matrix_invert(mu_inv, size_X_m);     // mu_inv: 3N x 3N

   // Matrix multiplication with cblas_dgemm()
   // Calculating C <- aA * B + bC (A: M x K, B: K x N, C: M x N)
   // In this case, we are multiplying grad_cv_t with mu_inv (D x 3N) x (3N x 3N)
   // So, M = D, K = N = 3N
   W = zeros(num_cv*size_X_m);  // D x 3N
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_cv, size_X_m, size_X_m, 1.0, grad_cv_t, size_X_m, mu_inv, size_X_m, 0.0, W, size_X_m);

   G_w_inv = zeros(num_cv*num_cv);
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, num_cv, num_cv, size_X_m, 1.0, W, size_X_m, grad_cv_t, size_X_m, 0.0, G_w_inv, num_cv);
   matrix_invert(G_w_inv, 2);          // G_w_inv: D x D

   first = icf_firstterm(W, G_w_inv, grad_V, size_X_m, num_cv);
   icf = icf_secondterm(X_m, size_X_m, ra_ind, rb_ind, cn_a_ind, cn_b_inds, size_cn_b_inds, n, m, r0, L, G_w_inv, W, mu_inv, grad_cv_t, num_cv, kT);

   cblas_daxpy(num_cv, 1.0, first, 1, icf, 1);

   free(first);
   free(grad_cv_t);
   free(W);
   free(G_w_inv);
   free(mu_inv);

   return icf;
}
