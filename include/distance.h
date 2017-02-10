#ifndef DISTANCE_H
#define DISTANCE_H

extern double
dr_dxi(const double*, const double*, const int, const int, const int, const double, const double);
extern double
d2r_dxjxi(const double*, const double*, const int, const int, const int, const int, const double, const double);
extern double*
r_grad_x(const double*, const int, const int, const int, const double);
extern double*
r_hess_x_j(const double*, const int, const int, const int, const int, const double);

#endif
