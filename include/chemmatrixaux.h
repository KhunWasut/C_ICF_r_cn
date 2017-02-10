#ifndef CHEMMATRIXAUX_H
#define CHEMMATRIXAUX_H

extern int
atomindex_to_vecindex(const int, const int);
extern int
vecindex_to_atomindex(const int);
extern int
vecindex_to_axisindex(const int);
extern char*
index_verify_1d(const int, const int, const int);
extern double
dx_pbc(const double*, const double*, const int, const double);
extern double
d_dx_pbc(const int, const int, const int, const int);

#endif
