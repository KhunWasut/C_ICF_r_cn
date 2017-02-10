#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <atlas/cblas.h>

#include "utils.h"
#include "matrix.h"
#include "icf2d.h"


int main(int argc, char* argv[])
{
   // The main program of 2D ICF calculation in C, implement to make it work on boom as well!
   // Command line options: -xdir, -fdir, -idir, -first, -last, --help
   // DIR prefixes require trailing slashes!!

   // Cycle through arguments
   int i, arglen, first_icf, last_icf, N;
   double r0, n, m, L;
   char *xdir, *fdir, *idir, *mpath;

   if(argc < 23)
   {
      // 11 arguments. If less than that, then print error message and a guideline, as well as exiting the program.
      printf("There should be 11 command line options, each with a specific argument!!");
      printf("\n");
      printf("-xdir (str): Directory storing position vectors.\n");
      printf("-fdir (str): Directory storing gradient vectors.\n");
      printf("-idir (str): Directory to store ICF values (in atomic units after this calculation).\n");
      printf("-mpath (str): The full path to the atomic mass vector formerly generated with Python.\n");
      printf("-first (int): The 1-based index of the first ICF file to be calculated.\n");
      printf("-last (int): The 1-based index of the last ICF file to be calculated. Cannot be equal to -first.\n");
      printf("-numpart (int): Number of atoms in the system.\n");
      printf("-r0 (double): The 'r0' parameter for C.N. calculation in Angstroms. (For LiCl sim, it's 2.59.)\n");
      printf("-n (double): The 'n' parameter for C.N. calculation. (For LiCl sim, it's 20.0)\n");
      printf("-m (double): The 'm' parameter for C.N. calculation. (For LiCl sim, it's 34.0)\n");
      printf("-L (double): PBC cell size in angstroms. (For LiCl sim, it's 19.033)");
      printf("\n");
      printf("ERROR: PLEASE INPUT ALL ABOVE OPTIONS!! EXITING...\n");
      exit(0);
   }

   for(i = 1; i < argc; i+=2)
   {
      if (strcmp(argv[i], "-xdir") == 0)
      {
         arglen = strlen(argv[i+1]);
         xdir = (char*)malloc(sizeof(char)*arglen + 1);  
         strcpy(xdir, argv[i+1]);
      }
      else if (strcmp(argv[i], "-fdir") == 0)
      {
         arglen = strlen(argv[i+1]);
         fdir = (char*)malloc(sizeof(char)*arglen + 1);  
         strcpy(fdir, argv[i+1]);
      }
      else if (strcmp(argv[i], "-idir") == 0)
      {
         arglen = strlen(argv[i+1]);
         idir = (char*)malloc(sizeof(char)*arglen + 1);  
         strcpy(idir, argv[i+1]);
      }
      else if (strcmp(argv[i], "-mpath") == 0)
      {
         arglen = strlen(argv[i+1]);
         mpath = (char*)malloc(sizeof(char)*arglen + 1); 
         strcpy(mpath, argv[i+1]);
      }
      else if (strcmp(argv[i], "-first") == 0)
         first_icf = atoi(argv[i+1]);
      else if (strcmp(argv[i], "-last") == 0)
         last_icf = atoi(argv[i+1]);
      else if (strcmp(argv[i], "-numpart") == 0)
         N = atoi(argv[i+1]);
      else if (strcmp(argv[i], "-r0") == 0)
         r0 = (double)(atof(argv[i+1]));
      else if (strcmp(argv[i], "-n") == 0)
         n = (double)(atof(argv[i+1]));
      else if (strcmp(argv[i], "-m") == 0)
         m = (double)(atof(argv[i+1]));
      else if (strcmp(argv[i], "-L") == 0)
         L = (double)(atof(argv[i+1]));
   }

   // Processing input files and store position and gradient data
   // Note that initially, these vectors will have NAMD units! Will convert to atomic unit for calculations.
   // If this code works, can test using NAMD unit and see if we have same result! Then if so, we can eliminate the conversion.
   double **positions, **gradients;
   char *x_filename, *f_filename;
   char numbuff[6];
   char *path_x, *path_f;
   int vector_counter = 0;

   // Allocate 2D arrays for positions and gradients. Size of outer dimension is number of ICFs to be calculated.
   positions = (double**)malloc(sizeof(double*)*(last_icf - first_icf + 1));
   gradients = (double**)malloc(sizeof(double*)*(last_icf - first_icf + 1));

   for(i = first_icf; i <= last_icf; i++)
   {
      // Path processing
      sprintf(numbuff, "%d", i);
      x_filename = (char*)malloc(sizeof(char)*(5 + strlen(numbuff)) + 1);
      f_filename = (char*)malloc(sizeof(char)*(10 + strlen(numbuff)) + 1);

      sprintf(x_filename, "x%s.vec", numbuff);
      sprintf(f_filename, "grad_V%s.vec", numbuff);

      path_x = (char*)malloc(sizeof(char)*(strlen(xdir) + strlen(x_filename)) + 1);
      path_f = (char*)malloc(sizeof(char)*(strlen(fdir) + strlen(f_filename)) + 1);

      sprintf(path_x, "%s%s", xdir, x_filename);
      sprintf(path_f, "%s%s", fdir, f_filename);

      free(x_filename);
      free(f_filename);

      // Actually read the file. Allocate the second dimension of positions and gradients
      positions[vector_counter] = readvec(path_x, 3*N);    
      gradients[vector_counter] = readvec(path_f, 3*N);    

      vector_counter++;

      free(path_x);
      free(path_f);
   }

   // Free xdir and fdir
   free(xdir);
   free(fdir);

   // Reading mass info
   double *masses, *mass_matrix;

   masses = readvec(mpath, N);

   free(mpath);

   // Unit conversion of all important vectors to atomic units. Use BLAS routines for scalar multiplication!!
   // BE CAREFUL!! BLAS operations do change the arrays!! This pointer business in C makes debugging a much harder business!!

   // Converting mass to atomic units
   cblas_dscal(N, (1.66054e-27/9.10938e-31), masses, 1);

   // Converting positions and gradients to atomic units
   for(i = 0; i < (last_icf - first_icf + 1); i++)
   {
      cblas_dscal(3*N, (1.0/0.529177), positions[i], 1);
      cblas_dscal(3*N, (0.529177/627.509), gradients[i], 1);
   }

   // Defining atomic unit constants and parameters
   const double L_atomic = L/0.529177;
   const double kT_atomic = 1.38064852e-23*300.0/4.35974e-18;
   const double r0_atomic = r0/0.529177;

   // Diagonalize mass vector to a 3N x 3N matrix (1D rep, row major, LDA is column size)
   // This is 'mu' in Gabor's paper
   mass_matrix = zeros(3*N*3*N);
   for(i = 0; i < 3*N; i+=3)
   {
      mass_matrix[3*N*i + i] = masses[i/3];
      mass_matrix[3*N*(i+1) + (i+1)] = masses[i/3];
      mass_matrix[3*N*(i+2) + (i+2)] = masses[i/3];
   }

   free(masses);

   // Defining indices
   const int li_index = 0;
   const int cl_index = 1;
   int *o_indices;
   const int size_o_indices = 230;

   o_indices = (int*)malloc(sizeof(int)*size_o_indices);

   for(i = 0; i < size_o_indices; i++)
      o_indices[i] = 3*i+2;

   // Actually calculate the 2D ICFs!!
   double **icf_atomic;
   char *path_icf, *icf_filename;

   icf_atomic = (double**)malloc(sizeof(double)*(last_icf - first_icf + 1));

   int k;
   k = 0;
   for(i = first_icf; i <= last_icf; i++)
   {
      printf("Evaluating ICF %d...\n", i);
      icf_atomic[k] = icf_construct_2d(positions[k], 3*N, mass_matrix, gradients[k], kT_atomic, li_index, cl_index, li_index, 
            o_indices, 230, n, m, r0_atomic, L_atomic, 2);

      // Write the k-th ICF into a file!!
      sprintf(numbuff, "%d", i);

      icf_filename = (char*)malloc(sizeof(char)*(14 + strlen(numbuff)) + 1);
      sprintf(icf_filename, "icf%s_atomic.vec", numbuff);

      path_icf = (char*)malloc(sizeof(char)*(strlen(idir) + strlen(icf_filename)) + 1);
      sprintf(path_icf, "%s%s", idir, icf_filename);

      writevec(path_icf, icf_atomic[k], 2);

      free(icf_filename);
      free(path_icf);

      k++;
   }

   // Free all the data arrays!!

   for(i = 0; i < (last_icf - first_icf + 1); i++)
   {
      free(icf_atomic[i]);
      free(positions[i]);
      free(gradients[i]);
   }

   free(idir);
   free(o_indices);
   free(icf_atomic);
   free(positions);
   free(gradients);
   free(mass_matrix);

   return 0;
}
