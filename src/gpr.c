#include "gpr.h"
#include "matrix.h"

#include <math.h>
#include <atlas/cblas.h>


/* GPR.C -- Gaussian Process Regression Code
 *
 * Description:
 *
 *    Predict function from derivative data in 1D and 2D
 *
 */

double kexp(const double* CV_diff, const double* theta_vec, const int num_cv)
{
   // For num_cv > 1
   int i;
   double exp_sum = 0.0;

   for (i = 0; i < num_cv; i++)
      exp_sum += (-0.5*(pow(CV_diff[i], 2.0))/(pow(theta_vec[i], 2.0)));

   return exp(exp_sum);
}


double kexp_1d(const double CV_diff, const double theta)
{
   return exp((-0.5*(pow(CV_diff, 2.0)))/(pow(theta, 2.0)));
}


// Constructs K for a relatively small matrix (< 20,000 x 20,000)
// For larger matrices, will create another code 'gprlarge.c'
// Big data implementation will need to split rows into files!!
double* K_construct(const double* CV_mat, const double* theta, const double chi, const int num_data, const int num_cv)
{
   // Based on the code on KLNX at path '/home/wpornpat/KTempWorkSpace/free-energy-modules/1-test-part1/python-script/fe_fitting.py'
   // For num_cv > 1
   // CV_mat is a D x n_data matrix
   // theta is D x 1 vector
   // For 2D, the order is r, cn
   int i, j, ii, jj, iii;
   double *K, *K_submat, *CV_diff;
   double frontterm, factor;

   K = (double*)malloc(sizeof(double)*num_cv*num_data*num_cv*num_data);

   for (i = 0; i < num_data; i++)
   {
      for (j = 0; j < num_data; j++)
      {
         // Populate CV_i, CV_j vectors, each of which D x 1 sized vector
         if (num_cv > 1)
         {
            CV_diff = (double*)malloc(sizeof(double)*num_cv);
            K_submat = zeros(num_cv*num*cv);

            for(iii = 0; iii < num_cv; iii++)
               CV_diff[iii] = CV_mat[num_data*iii + i] - CV_mat[num_data*iii + j];

            // Populate the submatrix!!
            for (ii = 0; ii < num_cv; ii++)
            {
               for (jj = 0; jj < num_cv; jj++)
               {
                  if (ii == jj)
                     frontterm = pow(theta[jj], 2.0);
                  else
                     frontterm = 0.0;

                  factor = frontterm - (((CV_diff[ii])/(pow(theta[ii], 2.0))) * ((CV_diff[jj])/(pow(theta[jj], 2.0))));

                  if (factor != 0.0)
                     K_submat[num_cv*ii + jj] = pow(chi, 2.0)*kexp(CV_diff, theta, num_cv)*factor;

                  // Populating K from the mapping scheme
                  // Carefully check if this mapping is correct!!
                  K[(num_data*num_cv)*(num_data*i + ii) + (num_data*j + jj)] = K_submat[num_cv*ii + jj];
               }
            }

            free(CV_diff);
            free(K_submat);
         }

         // Implement for 1D case later!!
         // Don't forget to think about cases for big data!!
      }
   }

   return K;
}


double* k_star(const double* CV_mat, const double* theta, const double* xi_star, const double chi, const int num_data, const int num_cv)
{
   int i, ii;
   double *CV_diff, *k_star;
   double frontterm;

   k_star = zeros(num_cv*num_data);

   for (i = 0; i < num_data; i++)
   {
      if (num_cv > 1)
      {
         CV_diff = (double*)malloc(sizeof(double)*num_cv);
         for(ii = 0; ii < num_cv; ii++)
            CV_diff[ii] = CV_mat[num_data*ii + i] - xi_star[ii];

         for(ii = 0; ii < num_cv; ii++)
         {
            frontterm = (-1.0)*(pow(chi, 2.0)*kexp(CV_diff, theta, num_cv));

            // Check this mapping very carefully!!
            k_star[num_cv*i + ii] = frontterm*((CV_diff[ii])/(pow(theta[ii], 2.0)));
         }

         free(CV_diff);
      }

      // Implement for 1D case later!!
      // Don't forget to think about cases for big data!!
   }

   return k_star;
}


double cond_expected(const double* CV_mat, const double* theta, const double* xi_star, const double* icf, const double* icf_var, const double chi,
      const int num_data, const int num_cv)
{
   double *KPlusVar_inv, *T, *k_star_vec;
   double ystar_expected;
   int i;

   KPlusVar_inv = K_construct(CV_mat, theta, chi, num_data, num_cv);

   // Add variance to K matrix to form (K + \sigma^2 I)
   // Carefully check the mapping!!
   for (i = 0; i < num_cv*num_data; i++)
      KPlusVar_inv[num_data*num_cv*i + i] += icf_var[i];

   // Inverting K + \sigma^2 I
   matrix_invert(KPlusVar_inv, num_data*num_cv);

   T = zeros(num_cv*num_data);
   cblas_dgemv(CblasRowMajor, CblasNoTrans, num_cv*num_data, num_cv*num_data, 1.0, KPlusVar_inv, icf, 1, 0.0, T, 1);

   // Free the matrix to save memory
   free(KPlusVar_inv);

   k_star_vec = k_star(CV_mat, theta, xi_star, chi, num_data, num_cv);

   ystar_expected = cblas_ddot(num_cv*num_data, T, 1, k_star_vec, 1);

   free(k_star_vec);
   free(T);

   return ystar_expected;
}
