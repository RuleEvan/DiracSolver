#include "matrix_element.h"
void four_pt_b_derivative(double (*df) (double, double, double, double, double), double (*dg) (double, double, double, double, double), double epsilon, double x_max, double x_match, double k, double b_tilde, int l, double a0, double a1, double *delta_d_f_d_b, double *delta_f_g_d_b);

void four_pt_a0_derivative(double (*df) (double, double, double, double, double), double (*dg) (double, double, double, double, double), double epsilon, double x_max, double x_match, double k, double e_tilde, int l, double a0, double a1, double *delta_d_f_d_a0, double *delta_f_g_d_a0);


void two_pt_a0_derivative(double (*df) (double, double, double, double, double), double (*dg) (double, double, double, double, double), double epsilon, double x_max, double x_match, double k, double e_tilde, int l, double a0, double a1, double *delta_d_f_d_a0, double *delta_d_g_d_a0);


void two_pt_b_derivative(double (*df) (double, double, double, double, double), double (*dg) (double, double, double, double, double), double epsilon, double x_max, double x_match, double k, double b_tilde, int l, double a0, double a1, double *delta_d_f_d_b, double *delta_d_g_d_b);

void delta_left_right(double (*df) (double, double, double, double, double), double (*dg) (double, double, double, double, double), double epsilon, double x_max, double x_match, double f_min, double f_max, double g_min, double g_max, double k, double e_tilde, double* delta_f, double* delta_g);

void delta_left_right_derivs(double (*df) (double, double, double, double, double), double (*dg) (double, double, double, double, double), double epsilon, double x_max, double x_match, double f_min, double f_max, double g_min, double g_max, double k, double e_tilde, double* delta_f, double* delta_g, double* delta_df, double* delta_dg);


void runge_kutta4_right(double (*df)(double, double, double, double, double), double (*dg)(double, double, double, double, double), double t0, double f0, double g0, double tmax, double k, double e_tilde, double *f_right, double *g_right);

void runge_kutta4_left(double (*df)(double, double, double, double, double), double (*dg)(double, double, double, double, double), double t0, double f0, double g0, double tmin, double k, double e_tilde, double *f_left, double *g_left);

void runge_kutta4_right_d(double (*df)(double, double, double, double, double), double (*dg)(double, double, double, double, double), double t0, double f0, double g0, double tmax, double k, double e_tilde, double *f_right, double *g_right, double *df_right, double *dg_right);

void runge_kutta4_left_d(double (*df)(double, double, double, double, double), double (*dg)(double, double, double, double, double), double t0, double f0, double g0, double tmin, double k, double e_tilde, double *f_left, double *g_left, double *df_left, double *dg_left);

double df_hydrogen(double x, double f, double g, double k, double e_tilde);

double dg_hydrogen(double t, double f, double g, double k, double e_tilde);

