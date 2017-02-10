#ifndef COORDNUM_H
#define COORDNUM_H

extern double
n_ab_k(const double, const double, const double);

extern double
d_ab_k(const double, const double, const double);

extern double
dN_xi(const double, const double, const double*, const double*, const int, const int, const int, const double, const double);

extern double
dD_xi(const double, const double, const double*, const double*, const int, const int, const int, const double, const double);

extern double*
cn_grad_x(const double*, const int, const int, const int*, const int, const double, const double, const double, const double);

extern double*
cn_hess_x_j(const double*, const int, const int, const int, const int*, const int, const double, const double, const double, const double);

#endif
