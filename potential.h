#ifndef POTENTIAL_H
#define POTENTIAL_H
#include "romberg.h"
double b_coeff(double n, double l, double np, double lp, double p);
double compute_potential(double n, double np, double l, double lp, int l_op, double k, int r_pow);
double a_coeff(double n, double l, double k);
double talmi(double p, int l, double k, int r_pow);
double talmi_integrand(double q, double p, int l, double k, int r_pow);
#endif
