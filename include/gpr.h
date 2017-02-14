#ifndef GPR_H
#define GPR_H

extern double
kexp(const double*, const double*, const int);

extern double
kexp_1d(const double, const double);

extern double*
K_construct(const double*, const double*, const double, const int, const int);

extern double*
k_star(const double*, const double*, const double*, const double, const int, const int);

extern double
cond_expected(const double*, const double*, const double*, const double*, const double*, const double, const int, const int);

#endif
