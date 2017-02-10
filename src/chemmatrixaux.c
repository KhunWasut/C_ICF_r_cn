#include "chemmatrixaux.h"

#include <stdlib.h>
#include <string.h>


int atomindex_to_vecindex(const int atom_index, const int axis_index)
{
   int vec_index;
   vec_index = (3*atom_index + axis_index);
   return vec_index;
}


int vecindex_to_atomindex(const int vec_index_i)
{
   int atom_index;
   atom_index = vec_index_i / 3;
   return atom_index;
}


int vecindex_to_axisindex(const int vec_index_i)
{
   int axis_index;
   axis_index = vec_index_i % 3;
   return axis_index;
}


char* index_verify_1d(const int vec_index_i, const int atom_a_index, const int atom_b_index)
{
   // C implementation of index_verify_1d
   char* ans;
   int atom_index;

   atom_index = vecindex_to_atomindex(vec_index_i);

   if (atom_index == atom_a_index)
   {
      ans = "FT";
      return ans;
   }
   else if (atom_index == atom_b_index)
   {
      ans = "LT";
      return ans;
   }
   else
   {
      ans = "NA";
      return ans;
   }
}


double dx_pbc(const double* x_a, const double* x_b, const int vec_index_i, const double L)
{
   int axis_index;
   double dx;

   // x_a and x_b are 3x1 vectors
   axis_index = vecindex_to_axisindex(vec_index_i);
   dx = x_a[axis_index] - x_b[axis_index];

   if ((x_a[axis_index] - x_b[axis_index]) > (L/2))
   {
      dx -= L;
   }
   else if ((x_a[axis_index] - x_b[axis_index]) < (-L/2))
   {
      dx += L;
   }

   return dx;
}


double d_dx_pbc(const int atom_a_index, const int atom_b_index, const int vec_index_k, const int vec_index_l)
{
   int domain_vec_indices[6] = {
      atomindex_to_vecindex(atom_a_index, 0),
      atomindex_to_vecindex(atom_a_index, 1),
      atomindex_to_vecindex(atom_a_index, 2),
      atomindex_to_vecindex(atom_b_index, 0),
      atomindex_to_vecindex(atom_b_index, 1),
      atomindex_to_vecindex(atom_b_index, 2)
   };

   int l_not_in_domain = 1;
   int i;

   // Check if l is in the domain
   for(i = 0; i < 6; i++)
   {
      if (vec_index_l == domain_vec_indices[i])
      {
         l_not_in_domain = 0;
         break;
      }
   }

   if (l_not_in_domain)
   {
      return 0.0;
   }
   else
   {
      if (vecindex_to_axisindex(vec_index_l) == vecindex_to_axisindex(vec_index_k))
      {
         if (strcmp("FT", index_verify_1d(vec_index_l, atom_a_index, atom_b_index)) == 0)
         {
            return 1.0;
         }
         else if (strcmp("LT", index_verify_1d(vec_index_l, atom_a_index, atom_b_index)) == 0)
         {
            return -1.0;
         }
         else
         {
            return 0.0;
         }
      }
      else
      {
         return 0.0;
      }
   }
}
