#include "utils.h"

#include <stdio.h>
#include <stdlib.h>

// Read the number from .vec file!
double* readvec(const char* path, const int size)
{
   FILE *fp;
   double *arr;
   double val;
   int counter;

   // Allocate a 1D array for the vector
   arr = (double *)malloc(sizeof(double)*size);

   // Open file
   fp = fopen(path, "r");

   // Read until EOF
   counter = 0;
   while (!feof(fp))
   {
      // Use fscanf
      fscanf(fp, "%lf", &val);
      arr[counter] = val;
      counter++;
   }

   fclose(fp);

   return arr;
}


void writevec(const char* path, const double* v, const int size)
{
   FILE *fp;
   int counter;

   fp = fopen(path, "w");

   for(counter = 0; counter < size; counter++)
      fprintf(fp, "%.12lf\n", v[counter]);

   fclose(fp);
}


double* deep_copy(const double* origin, const int origin_size)
{
   int i;
   double *arr;
   arr = (double*)malloc(sizeof(double)*origin_size);

   for(i = 0; i < origin_size; i++)
      arr[i] = origin[i];

   return arr;
}
