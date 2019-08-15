#include "runge_kutta.h"

void two_pt_b_derivative(double (*df) (double, double, double, double, double), double (*dg) (double, double, double, double, double), double epsilon, double x_max, double x_match, double k, double b_tilde, int l, double a0, double a1, double *d_delta_f_d_b, double *d_delta_g_d_b) {

  double deriv_tol = 0.0001;
  double b_tilde_up = b_tilde*(1.0 + deriv_tol);
  double b_tilde_down = b_tilde*(1.0 - deriv_tol);

  double e_tilde_up = 1.0 - b_tilde_up;
  double e_tilde_down = 1.0 - b_tilde_down;
  
  double mu_e_up = sqrt(1.0 - pow(e_tilde_up, 2.0));
  double g_max_e_up = a1*exp(-mu_e_up*x_max);
  double f_max_e_up = -a1*pow((1.0 - e_tilde_up)/(1.0 + e_tilde_up), 0.5)*exp(-mu_e_up*x_max);

  double mu_e_down = sqrt(1.0 - pow(e_tilde_down, 2.0));
  double g_max_e_down = a1*exp(-mu_e_down*x_max);
  double f_max_e_down = -a1*pow((1.0 - e_tilde_down)/(1.0 + e_tilde_down), 0.5)*exp(-mu_e_down*x_max);

  double delta_f_b_up = 0.0;
  double delta_g_b_up = 0.0;

  double f_min, g_min;
  if (k < 0) {
    g_min = a0*pow(epsilon, l + 1);
    f_min = -a0*Z_ATOM*ALPHA_FS/(2.0*l + 3.0)*pow(epsilon, l + 1);
  } else if (k >= 0) {
    f_min = a0*pow(epsilon, l);
    g_min = a0*Z_ATOM*ALPHA_FS/(2.0*l + 1.0)*pow(epsilon, l);
  }

  delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min, f_max_e_up, g_min, g_max_e_up, k, e_tilde_up, &delta_f_b_up, &delta_g_b_up); 
  double delta_f_b_down = 0.0;
  double delta_g_b_down = 0.0;
  delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min, f_max_e_down, g_min, g_max_e_down, k, e_tilde_down, &delta_f_b_down, &delta_g_b_down); 

  *d_delta_f_d_b = (delta_f_b_up - delta_f_b_down)/(2.0*deriv_tol*b_tilde);
  *d_delta_g_d_b = (delta_g_b_up - delta_g_b_down)/(2.0*deriv_tol*b_tilde);
  printf("2pt: df/dB: %g dg/dB: %g\n", *d_delta_g_d_b, *d_delta_f_d_b);
 
  return;
}


void two_pt_a0_derivative(double (*df) (double, double, double, double, double), double (*dg) (double, double, double, double, double), double epsilon, double x_max, double x_match, double k, double e_tilde, int l, double a0, double a1, double *d_delta_f_d_a0, double *d_delta_g_d_a0) {

  double deriv_tol = 0.0001;
  double a0_up = a0 + fabs(a0)*deriv_tol;
  double a0_down = a0 - fabs(a0)*deriv_tol;
  double mu = sqrt(1.0 - pow(e_tilde, 2.0));
  double g_max = a1*exp(-mu*x_max);
  double f_max = -a1*pow((1.0 - e_tilde)/(1.0 + e_tilde), 0.5)*exp(-mu*x_max);

  double f_min_a0_up, g_min_a0_up, f_min_a0_down, g_min_a0_down;
  if (k < 0) {
    g_min_a0_up = a0_up*pow(epsilon, l + 1);
    f_min_a0_up = -a0_up*Z_ATOM*ALPHA_FS/(2.0*l + 3.0)*pow(epsilon, l + 1);
    g_min_a0_down = a0_down*pow(epsilon, l + 1);
    f_min_a0_down = -a0_down*Z_ATOM*ALPHA_FS/(2.0*l + 3.0)*pow(epsilon, l + 1);

  } else if (k >= 0) {
    f_min_a0_up = a0_up*pow(epsilon, l);
    g_min_a0_up = a0_up*Z_ATOM*ALPHA_FS/(2.0*l + 1.0)*pow(epsilon, l);
    f_min_a0_down = a0_down*pow(epsilon, l);
    g_min_a0_down = a0_down*Z_ATOM*ALPHA_FS/(2.0*l + 1.0)*pow(epsilon, l);

  }
  double delta_f_a0_up = 0.0;
  double delta_g_a0_up = 0.0;

  delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min_a0_up, f_max, g_min_a0_up, g_max, k, e_tilde, &delta_f_a0_up, &delta_g_a0_up); 
  double delta_f_a0_down = 0.0;
  double delta_g_a0_down = 0.0;
  delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min_a0_down, f_max, g_min_a0_down, g_max, k, e_tilde, &delta_f_a0_down, &delta_g_a0_down); 

  *d_delta_f_d_a0 = (delta_f_a0_up - delta_f_a0_down)/(2.0*deriv_tol*fabs(a0));
  *d_delta_g_d_a0 = (delta_g_a0_up - delta_g_a0_down)/(2.0*deriv_tol*fabs(a0));
  printf("2pt: df/da0: %g dg/da0: %g\n", *d_delta_g_d_a0, *d_delta_f_d_a0);
 
  return;
}

void four_pt_a0_derivative(double (*df) (double, double, double, double, double), double (*dg) (double, double, double, double, double), double epsilon, double x_max, double x_match, double k, double e_tilde, int l, double a0, double a1, double *d_delta_f_d_a0, double *d_delta_g_d_a0) {
  double deriv_tol = 0.0001;
  double a0_up1 = a0 + deriv_tol;
  double a0_up2 = a0 + 2.0*deriv_tol;
  double a0_down1 = a0 - deriv_tol;
  double a0_down2 = a0 - 2.0*deriv_tol;

  double mu = sqrt(1.0 - pow(e_tilde, 2.0));
  double g_max = a1*exp(-mu*x_max);
  double f_max = -a1*pow((1.0 - e_tilde)/(1.0 + e_tilde), 0.5)*exp(-mu*x_max);
 
  double delta_f_a0_up1 = 0.0;
  double delta_g_a0_up1 = 0.0;
  double delta_f_a0_up2 = 0.0;
  double delta_g_a0_up2 = 0.0;

  double f_min_a0_up1, f_min_a0_up2, f_min_a0_down1, f_min_a0_down2, g_min_a0_up1, g_min_a0_up2, g_min_a0_down1, g_min_a0_down2;

  if (k < 0) {
    g_min_a0_up1 = a0_up1*pow(epsilon, l + 1);
    f_min_a0_up1 = -a0_up1*Z_ATOM*ALPHA_FS/(2.0*l + 3.0)*pow(epsilon, l + 1);
    g_min_a0_up2 = a0_up2*pow(epsilon, l + 1);
    f_min_a0_up2 = -a0_up2*Z_ATOM*ALPHA_FS/(2.0*l + 3.0)*pow(epsilon, l + 1);
    g_min_a0_down1 = a0_down1*pow(epsilon, l + 1);
    f_min_a0_down1 = -a0_down1*Z_ATOM*ALPHA_FS/(2.0*l + 3.0)*pow(epsilon, l + 1);
    g_min_a0_down2 = a0_down2*pow(epsilon, l + 1);
    f_min_a0_down2 = -a0_down2*Z_ATOM*ALPHA_FS/(2.0*l + 3.0)*pow(epsilon, l + 1);

  } else if (k >= 0) {
    f_min_a0_up1 = a0_up1*pow(epsilon, l);
    g_min_a0_up1 = a0_up1*Z_ATOM*ALPHA_FS/(2.0*l + 1.0)*pow(epsilon, l);
    f_min_a0_up2 = a0_up2*pow(epsilon, l);
    g_min_a0_up2 = a0_up2*Z_ATOM*ALPHA_FS/(2.0*l + 1.0)*pow(epsilon, l);
    f_min_a0_down1 = a0_down1*pow(epsilon, l);
    g_min_a0_down1 = a0_down1*Z_ATOM*ALPHA_FS/(2.0*l + 1.0)*pow(epsilon, l);
    f_min_a0_down2 = a0_down2*pow(epsilon, l);
    g_min_a0_down2 = a0_down2*Z_ATOM*ALPHA_FS/(2.0*l + 1.0)*pow(epsilon, l);
  }

  delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min_a0_up1, f_max, g_min_a0_up1, g_max, k, e_tilde, &delta_f_a0_up1, &delta_g_a0_up1); 
  delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min_a0_up2, f_max, g_min_a0_up2, g_max, k, e_tilde, &delta_f_a0_up2, &delta_g_a0_up2); 

  double delta_f_a0_down1 = 0.0;
  double delta_g_a0_down1 = 0.0;
  double delta_f_a0_down2 = 0.0;
  double delta_g_a0_down2 = 0.0;

  delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min_a0_down1, f_max, g_min_a0_down1, g_max, k, e_tilde, &delta_f_a0_down1, &delta_g_a0_down1); 
  delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min_a0_down2, f_max, g_min_a0_down2, g_max, k, e_tilde, &delta_f_a0_down2, &delta_g_a0_down2); 

  *d_delta_f_d_a0 = (-delta_f_a0_up2 + 8.0*delta_f_a0_up1 - 8.0*delta_f_a0_down1 + delta_f_a0_down2)/(12.0*deriv_tol);
  *d_delta_g_d_a0 = (-delta_g_a0_up2 + 8.0*delta_g_a0_up1 - 8.0*delta_g_a0_down1 + delta_g_a0_down2)/(12.0*deriv_tol);

  printf("4pt: df/da0: %g dg/da0: %g\n", *d_delta_g_d_a0, *d_delta_f_d_a0);
 
  return;
}

void four_pt_b_derivative(double (*df) (double, double, double, double, double), double (*dg) (double, double, double, double, double), double epsilon, double x_max, double x_match, double k, double b_tilde, int l, double a0, double a1, double *d_delta_f_d_b, double *d_delta_g_d_b) {
  double deriv_tol = 0.000001;
  double b_tilde_up1 = b_tilde*(1.0 + deriv_tol);
  double b_tilde_up2 = b_tilde*(1.0 + 2.0*deriv_tol);
  double b_tilde_down1 = b_tilde*(1.0 - deriv_tol);
  double b_tilde_down2 = b_tilde*(1.0 - 2.0*deriv_tol);

  double e_tilde_up1 = 1.0 - b_tilde_up1;
  double e_tilde_up2 = 1.0 - b_tilde_up2;
  double e_tilde_down1 = 1.0 - b_tilde_down1;
  double e_tilde_down2 = 1.0 - b_tilde_down2;

  double mu_e_up1 = sqrt(1.0 - pow(e_tilde_up1, 2.0));
  double g_max_e_up1 = a1*exp(-mu_e_up1*x_max);
  double f_max_e_up1 = -a1*pow((1.0 - e_tilde_up1)/(1.0 + e_tilde_up1), 0.5)*exp(-mu_e_up1*x_max);
 
  double mu_e_up2 = sqrt(1.0 - pow(e_tilde_up2, 2.0));
  double g_max_e_up2 = a1*exp(-mu_e_up2*x_max);
  double f_max_e_up2 = -a1*pow((1.0 - e_tilde_up2)/(1.0 + e_tilde_up2), 0.5)*exp(-mu_e_up2*x_max);


  double mu_e_down1 = sqrt(1.0 - pow(e_tilde_down1, 2.0));
  double g_max_e_down1 = a1*exp(-mu_e_down1*x_max);
  double f_max_e_down1 = -a1*pow((1.0 - e_tilde_down1)/(1.0 + e_tilde_down1), 0.5)*exp(-mu_e_down1*x_max);

  double mu_e_down2 = sqrt(1.0 - pow(e_tilde_down2, 2.0));
  double g_max_e_down2 = a1*exp(-mu_e_down2*x_max);
  double f_max_e_down2 = -a1*pow((1.0 - e_tilde_down2)/(1.0 + e_tilde_down2), 0.5)*exp(-mu_e_down2*x_max);

  double delta_f_b_up1 = 0.0;
  double delta_g_b_up1 = 0.0;
  double delta_f_b_up2 = 0.0;
  double delta_g_b_up2 = 0.0;

  double f_min, g_min;
  if (k < 0) {
    g_min = a0*pow(epsilon, l + 1);
    f_min = -a0*Z_ATOM*ALPHA_FS/(2.0*l + 3.0)*pow(epsilon, l + 1);
  } else if (k >= 0) {
    k = l;
    f_min = a0*pow(epsilon, l);
    g_min = a0*Z_ATOM*ALPHA_FS/(2.0*l + 1.0)*pow(epsilon, l);
  }

  delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min, f_max_e_up1, g_min, g_max_e_up1, k, e_tilde_up1, &delta_f_b_up1, &delta_g_b_up1); 
  delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min, f_max_e_up2, g_min, g_max_e_up2, k, e_tilde_up2, &delta_f_b_up2, &delta_g_b_up2); 

  double delta_f_b_down1 = 0.0;
  double delta_g_b_down1 = 0.0;
  double delta_f_b_down2 = 0.0;
  double delta_g_b_down2 = 0.0;

  delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min, f_max_e_down1, g_min, g_max_e_down1, k, e_tilde_down1, &delta_f_b_down1, &delta_g_b_down1); 
  delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min, f_max_e_down2, g_min, g_max_e_down2, k, e_tilde_down2, &delta_f_b_down2, &delta_g_b_down2); 

  *d_delta_f_d_b = (-delta_f_b_up2 + 8.0*delta_f_b_up1 - 8.0*delta_f_b_down1 + delta_f_b_down2)/(12.0*deriv_tol*b_tilde);
  *d_delta_g_d_b = (-delta_g_b_up2 + 8.0*delta_g_b_up1 - 8.0*delta_g_b_down1 + delta_g_b_down2)/(12.0*deriv_tol*b_tilde);

  printf("4pt: df/dB: %g dg/dB: %g\n", *d_delta_g_d_b, *d_delta_f_d_b);
 
  return;
}


void delta_left_right(double (*df) (double, double, double, double, double), double (*dg) (double, double, double, double, double), double epsilon, double x_max, double x_match, double f_min, double f_max, double g_min, double g_max, double k, double e_tilde, double* delta_f, double* delta_g) {

  double g_right = 0.0;
  double g_left = 0.0;
  double f_right = 0.0;
  double f_left = 0.0;
  runge_kutta4_right(&df_hydrogen, &dg_hydrogen, epsilon, f_min, g_min, x_match, k, e_tilde, &f_right, &g_right);
  runge_kutta4_left(&df_hydrogen, &dg_hydrogen, x_max, f_max, g_max, x_match, k, e_tilde, &f_left, &g_left);
 // printf("f_right: %g, g_right: %g, f_left: %g, g_left: %g\n", f_right, g_right, f_left, g_left);
  *delta_g = 2.0*(g_right - g_left)/(g_right + g_left);
  *delta_f = 2.0*(f_right - f_left)/(f_right + f_left);

  return;
}

void delta_left_right_derivs(double (*df) (double, double, double, double, double), double (*dg) (double, double, double, double, double), double epsilon, double x_max, double x_match, double f_min, double f_max, double g_min, double g_max, double k, double e_tilde, double* delta_f, double* delta_g, double* delta_df, double* delta_dg) {

  double g_right = 0.0;
  double g_left = 0.0;
  double f_right = 0.0;
  double f_left = 0.0;
  double dg_right = 0.0;
  double dg_left = 0.0;
  double df_right = 0.0;
  double df_left = 0.0;
  printf("fmin: %g fmax: %g gmin: %g gmax: %g\n", f_min, f_max, g_min, g_max);
  runge_kutta4_right_d(&df_hydrogen, &dg_hydrogen, epsilon, f_min, g_min, x_match, k, e_tilde, &f_right, &g_right, &df_right, &dg_right);
  runge_kutta4_left_d(&df_hydrogen, &dg_hydrogen, x_max, f_max, g_max, x_match, k, e_tilde, &f_left, &g_left, &df_left, &dg_left);
  printf("f_right: %g, g_right: %g, f_left: %g, g_left: %g\n", f_right, g_right, f_left, g_left);

  *delta_g = 2.0*(g_right - g_left)/(g_right + g_left);
  *delta_f = 2.0*(f_right - f_left)/(f_right + f_left);
  *delta_dg = 2.0*(dg_right - dg_left)/(dg_right + dg_left);
  *delta_df = 2.0*(df_right - df_left)/(df_right + df_left);

  return;
}


void runge_kutta4_right(double (*df)(double, double, double, double, double), double (*dg)(double, double, double, double, double), double t0, double f0, double g0, double tmax, double k, double e_tilde, double* f_final, double* g_final) {
  double tn = t0;
  double fn = f0;
  double gn = g0;
  double h = 0.001;
//  printf("t0: %g f0: %g g0: %g tmax: %g\n", tn, fn, gn, tmax);
//  printf("k: %g ET: %g\n", k, e_tilde);
  while (tn <= tmax - h) {
    double k1 = h*df(tn, fn, gn, k, e_tilde);
    double h1 = h*dg(tn, fn, gn, k, e_tilde);
    double k2 = h*df(tn + h/2.0, fn + k1/2.0, gn + h1/2.0, k, e_tilde);
    double h2 = h*dg(tn + h/2.0, fn + k1/2.0, gn + h1/2.0, k, e_tilde);
    double k3 = h*df(tn + h/2.0, fn + k2/2.0, gn + h2/2.0, k, e_tilde);
    double h3 = h*dg(tn + h/2.0, fn + k2/2.0, gn + h2/2.0, k, e_tilde);
    double k4 = h*df(tn + h, fn + k3, gn + h3, k, e_tilde);
    double h4 = h*dg(tn + h, fn + k3, gn + h3, k, e_tilde);
    fn += 1.0/6.0*(k1 + 2*k2 + 2*k3 + k4);
    gn += 1.0/6.0*(h1 + 2*h2 + 2*h3 + h4);
    tn += h;
  }
  double dhf = tmax - tn;
  double k1f = dhf*df(tn, fn, gn, k, e_tilde);
  double h1f = dhf*dg(tn, fn, gn, k, e_tilde);
  double k2f = dhf*df(tn + dhf/2.0, fn + k1f/2.0, gn + h1f/2.0, k, e_tilde);
  double h2f = dhf*dg(tn + dhf/2.0, fn + k1f/2.0, gn + h1f/2.0, k, e_tilde);
  double k3f = dhf*df(tn + dhf/2.0, fn + k2f/2.0, gn + h2f/2.0, k, e_tilde);
  double h3f = dhf*dg(tn + dhf/2.0, fn + k2f/2.0, gn + h2f/2.0, k, e_tilde);
  double k4f = dhf*df(tn + dhf, fn + k3f, gn + h3f, k, e_tilde);
  double h4f = dhf*dg(tn + dhf, fn + k3f, gn + h3f, k, e_tilde);
  fn += 1.0/6.0*(k1f + 2*k2f + 2*k3f + k4f);
  gn += 1.0/6.0*(h1f + 2*h2f + 2*h3f + h4f);

  *f_final = fn;
  *g_final = gn;

  return;
}

void runge_kutta4_right_d(double (*df)(double, double, double, double, double), double (*dg)(double, double, double, double, double), double t0, double f0, double g0, double tmax, double k, double e_tilde, double* f_final, double* g_final, double* df_final, double* dg_final) {
  double tn = t0;
  double fn = f0;
  double gn = g0;
  double h = 0.001;
//  printf("t0: %g f0: %g g0: %g tmax: %g\n", tn, fn, gn, tmax);
//  printf("k: %g ET: %g\n", k, e_tilde);
  while (tn < tmax - h) {
    double k1 = h*df(tn, fn, gn, k, e_tilde);
    double h1 = h*dg(tn, fn, gn, k, e_tilde);
    double k2 = h*df(tn + h/2.0, fn + k1/2.0, gn + h1/2.0, k, e_tilde);
    double h2 = h*dg(tn + h/2.0, fn + k1/2.0, gn + h1/2.0, k, e_tilde);
    double k3 = h*df(tn + h/2.0, fn + k2/2.0, gn + h2/2.0, k, e_tilde);
    double h3 = h*dg(tn + h/2.0, fn + k2/2.0, gn + h2/2.0, k, e_tilde);
    double k4 = h*df(tn + h, fn + k3, gn + h3, k, e_tilde);
    double h4 = h*dg(tn + h, fn + k3, gn + h3, k, e_tilde);
    fn += 1.0/6.0*(k1 + 2*k2 + 2*k3 + k4);
    gn += 1.0/6.0*(h1 + 2*h2 + 2*h3 + h4);
    tn += h;
  }
  double dhf = tmax - tn;
  double k1f = dhf*df(tn, fn, gn, k, e_tilde);
  double h1f = dhf*dg(tn, fn, gn, k, e_tilde);
  double k2f = dhf*df(tn + dhf/2.0, fn + k1f/2.0, gn + h1f/2.0, k, e_tilde);
  double h2f = dhf*dg(tn + dhf/2.0, fn + k1f/2.0, gn + h1f/2.0, k, e_tilde);
  double k3f = dhf*df(tn + dhf/2.0, fn + k2f/2.0, gn + h2f/2.0, k, e_tilde);
  double h3f = dhf*dg(tn + dhf/2.0, fn + k2f/2.0, gn + h2f/2.0, k, e_tilde);
  double k4f = dhf*df(tn + dhf, fn + k3f, gn + h3f, k, e_tilde);
  double h4f = dhf*dg(tn + dhf, fn + k3f, gn + h3f, k, e_tilde);
  fn += 1.0/6.0*(k1f + 2*k2f + 2*k3f + k4f);
  gn += 1.0/6.0*(h1f + 2*h2f + 2*h3f + h4f);

  *f_final = fn;
  *g_final = gn;
  *df_final = df(tmax, fn, gn, k, e_tilde);
  *dg_final = dg(tmax, fn, gn, k, e_tilde);

  return;
}

 
void runge_kutta4_left(double (*df)(double, double, double, double, double), double (*dg)(double, double, double, double, double), double t0, double f0, double g0, double tmin, double k, double e_tilde, double* f_left, double *g_left) {
  double tn = t0;
  double fn = f0;
  double gn = g0;
  double h = 0.001;
  while (tn >= tmin + h) {
    double k1 = -h*df(tn, fn, gn, k, e_tilde);
    double h1 = -h*dg(tn, fn, gn, k, e_tilde);
    double k2 = -h*df(tn - h/2.0, fn + k1/2.0, gn + h1/2.0, k, e_tilde);
    double h2 = -h*dg(tn - h/2.0, fn + k1/2.0, gn + h1/2.0, k, e_tilde);
    double k3 = -h*df(tn - h/2.0, fn + k2/2.0, gn + h2/2.0, k, e_tilde);
    double h3 = -h*dg(tn - h/2.0, fn + k2/2.0, gn + h2/2.0, k, e_tilde);
    double k4 = -h*df(tn - h, fn + k3, gn + h3, k, e_tilde);
    double h4 = -h*dg(tn - h, fn + k3, gn + h3, k, e_tilde);
    fn += 1.0/6.0*(k1 + 2*k2 + 2*k3 + k4);
    gn += 1.0/6.0*(h1 + 2*h2 + 2*h3 + h4);
    tn -= h;
  }
  double dhf = tn - tmin;
  double k1f = -dhf*df(tn, fn, gn, k, e_tilde);
  double h1f = -dhf*dg(tn, fn, gn, k, e_tilde);
  double k2f = -dhf*df(tn - dhf/2.0, fn + k1f/2.0, gn + h1f/2.0, k, e_tilde);
  double h2f = -dhf*dg(tn - dhf/2.0, fn + k1f/2.0, gn + h1f/2.0, k, e_tilde);
  double k3f = -dhf*df(tn - dhf/2.0, fn + k2f/2.0, gn + h2f/2.0, k, e_tilde);
  double h3f = -dhf*dg(tn - dhf/2.0, fn + k2f/2.0, gn + h2f/2.0, k, e_tilde);
  double k4f = -dhf*df(tn - dhf, fn + k3f, gn + h3f, k, e_tilde);
  double h4f = -dhf*dg(tn - dhf, fn + k3f, gn + h3f, k, e_tilde);
  fn += 1.0/6.0*(k1f + 2*k2f + 2*k3f + k4f);
  gn += 1.0/6.0*(h1f + 2*h2f + 2*h3f + h4f);

  *f_left = fn;
  *g_left = gn;

  return;
}

void runge_kutta4_left_d(double (*df)(double, double, double, double, double), double (*dg)(double, double, double, double, double), double t0, double f0, double g0, double tmin, double k, double e_tilde, double* f_left, double *g_left, double *df_left, double *dg_left) {
  double tn = t0;
  double fn = f0;
  double gn = g0;
  double h = 0.001;
  while (tn > tmin + h) {
    double k1 = -h*df(tn, fn, gn, k, e_tilde);
    double h1 = -h*dg(tn, fn, gn, k, e_tilde);
    double k2 = -h*df(tn - h/2.0, fn + k1/2.0, gn + h1/2.0, k, e_tilde);
    double h2 = -h*dg(tn - h/2.0, fn + k1/2.0, gn + h1/2.0, k, e_tilde);
    double k3 = -h*df(tn - h/2.0, fn + k2/2.0, gn + h2/2.0, k, e_tilde);
    double h3 = -h*dg(tn - h/2.0, fn + k2/2.0, gn + h2/2.0, k, e_tilde);
    double k4 = -h*df(tn - h, fn + k3, gn + h3, k, e_tilde);
    double h4 = -h*dg(tn - h, fn + k3, gn + h3, k, e_tilde);
    fn += 1.0/6.0*(k1 + 2*k2 + 2*k3 + k4);
    gn += 1.0/6.0*(h1 + 2*h2 + 2*h3 + h4);
    tn -= h;
  }
  double dhf = tn - tmin;
  double k1f = -dhf*df(tn, fn, gn, k, e_tilde);
  double h1f = -dhf*dg(tn, fn, gn, k, e_tilde);
  double k2f = -dhf*df(tn - dhf/2.0, fn + k1f/2.0, gn + h1f/2.0, k, e_tilde);
  double h2f = -dhf*dg(tn - dhf/2.0, fn + k1f/2.0, gn + h1f/2.0, k, e_tilde);
  double k3f = -dhf*df(tn - dhf/2.0, fn + k2f/2.0, gn + h2f/2.0, k, e_tilde);
  double h3f = -dhf*dg(tn - dhf/2.0, fn + k2f/2.0, gn + h2f/2.0, k, e_tilde);
  double k4f = -dhf*df(tn - dhf, fn + k3f, gn + h3f, k, e_tilde);
  double h4f = -dhf*dg(tn - dhf, fn + k3f, gn + h3f, k, e_tilde);
  fn += 1.0/6.0*(k1f + 2*k2f + 2*k3f + k4f);
  gn += 1.0/6.0*(h1f + 2*h2f + 2*h3f + h4f);

  *f_left = fn;
  *g_left = gn;
  *df_left = df(tmin, fn, gn, k, e_tilde);
  *dg_left = dg(tmin, fn, gn, k, e_tilde);

  return;
}



double df_hydrogen(double x, double f, double g, double k, double e_tilde) {
  double df = k/x*f + (1.0 - e_tilde - Z_ATOM*ALPHA_FS/x)*g;

  return df;
}

double dg_hydrogen(double x, double f, double g, double k, double e_tilde) {
  double dg = -k/x*g + (1.0 + e_tilde + Z_ATOM*ALPHA_FS/x)*f;

  return dg;
}
