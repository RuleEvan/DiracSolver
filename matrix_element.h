#ifndef MATRIX_ELEMENT_H
#define MATRIX_ELEMENT_H
#include "potential.h"
double dpot1(int np, int lp, int n, int l, int l_op, double k);
double dpot2(int np, int lp, int n, int l, int l_op, double k);
double em_mat(int np, int lp, double jp, int n, int l, double j, int j_op, double k);
double em_dot_sig_mat(int np, int lp, double jp, int n, int l, double j, int j_op, int l_op, double k);
double em_dot_grad_mat(int np, int lp, double jp, int n, int l, double j, int j_op, int l_op, double k);
double t_el_c(int np, int lp, double jp, int n, int l, double j, int j_op, double k);
double t_el_mu(int np, int lp, double jp, int n, int l, double j, int j_op, double k);
double t_mag_c(int np, int lp, double jp, int n, int l, double j, int j_op, double k);
double t_mag_mu(int np, int lp, double jp, int n, int l, double j, int j_op, double k);
double decay_rate();
#endif
