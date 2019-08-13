#include "runge_kutta.h"

int main(int argc, char *argv[]) {
  double epsilon = 0.0001;
  double x_bohr = 1.0/ALPHA_FS;
  double x_max = 7.0*x_bohr;
  double x_match = 0.5*x_bohr;
  double a1 = 10.0;
  double a0 = 0.015;
  double e_tilde = 0.999977;
  double b = 12.0*pow(10, -6);
  int l = 0;
  double j = 0.5;
  double mu = sqrt(1.0 - pow(e_tilde, 2.0));
  double g_max = a1*exp(-mu*x_max);
  double f_max = -a1*pow((1.0 - e_tilde)/(1.0 + e_tilde), 0.5)*exp(-mu*x_max);
  double g_min, f_min;
  double k;
  if (j == l + 0.5) {
    k = -(l + 1);
    g_min = a0*pow(epsilon, l + 1);
    f_min = -a0*Z_ATOM*ALPHA_FS/(2.0*l + 3.0)*pow(epsilon, l + 1);
  } else if (j == l - 0.5) {
    k = l;
    f_min = a0*pow(epsilon, l);
    g_min = a0*Z_ATOM*ALPHA_FS/(2.0*l + 1.0)*pow(epsilon, l);
  }
  double delta_f = 1.0, delta_g = 1.0;
//  printf("k: %g\n", k);
 // printf("x_min: %g x_max: %g x_match: %g\n", epsilon, x_max, x_match);
 // printf("fmin: %g fmax: %g gmin: %g gmax: %g\n", f_min, f_max, g_min, g_max);
  double tol = pow(10, -3);
  while ((fabs(delta_f) > tol) && (fabs(delta_g) > tol)) {
    delta_left_right(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min, f_max, g_min, g_max, k, e_tilde, &delta_f, &delta_g); 
    printf("delta_f: %g delta_g: %g\n", delta_f, delta_g);
    double d_delta_f_d_e, d_delta_g_d_e, d_delta_f_d_a0, d_delta_g_d_a0;
    two_pt_e_derivative(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, k, e_tilde, l, a0, a1, &d_delta_f_d_e, &d_delta_g_d_e);
//  four_pt_e_derivative(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, k, e_tilde, l, a0, a1);
    two_pt_a0_derivative(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, k, e_tilde, l, a0, a1, &d_delta_f_d_a0, &d_delta_g_d_a0);

    double c1 = d_delta_g_d_e*e_tilde + d_delta_g_d_a0*a0 - delta_g;
    double c2 = d_delta_f_d_e*e_tilde + d_delta_f_d_a0*a0 - delta_f;

    double det_m = d_delta_g_d_e*d_delta_f_d_a0 - d_delta_f_d_e*d_delta_g_d_a0;
  
    e_tilde = 1.0/det_m*(c1*d_delta_f_d_a0 - c2*d_delta_g_d_a0);
    a0 = 1.0/det_m*(-c1*d_delta_f_d_e + c2*d_delta_g_d_e);
  }
  double binding = (1.0 - e_tilde)*M_ELECTRON*pow(10, 6);
  printf("e_tilde: %g B: %g a0: %g\n", e_tilde, binding, a0);
 
  return 0;
}
