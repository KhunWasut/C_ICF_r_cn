#ifndef ICF2D_H
#define ICF2D_H

extern double*
icf_firstterm(const double*, const double*, const double*, const int, const int);

extern double*
icf_secondterm(const double*, const int, const int, const int, const int, const int*, const int, const double, const double,
      const double, const double, const double*, const double*, const double*, const double*, const int, const double);

extern double*
icf_construct_2d(const double*, const int, const double*, const double*, const double, const int, const int, const int, const int*, const int,
      const double, const double, const double, const double, const int);

#endif
