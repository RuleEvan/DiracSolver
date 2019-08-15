#include "runge_kutta.h"

int main(int argc, char *argv[]) {
  double epsilon = 0.0001;
  double x_bohr = 1.0/ALPHA_FS;
  double x_max = 16.0*x_bohr/Z_ATOM;
  double x_match = 2.0*x_bohr/Z_ATOM;
  double a1 = 1000.0;
  double a0 = 0.0002;
  double b = 3.2*pow(10, -6);
  double b_tilde = b/M_ELECTRON;
  double e_tilde = 1.0 - b_tilde;
  int l = 1;
  double j = 0.5;
  double mu = sqrt(1.0 - pow(e_tilde, 2.0));
  double g_max = a1*exp(-mu*x_max);
  double f_max = -a1*pow((1.0 - e_tilde)/(1.0 + e_tilde), 0.5)*exp(-mu*x_max);
  printf("mu*xmax: %g f_max: %g g_max: %g\n", mu*x_max, f_max, g_max);
  printf("ep: %g x_max: %g x_match: %g\n", epsilon, x_max, x_match);
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
  double delta_f = 1.0;
  double delta_g = 1.0;
  double delta_df = 1.0; 
  double delta_dg = 1.0;
  printf("k: %g\n", k);
 // printf("x_min: %g x_max: %g x_match: %g\n", epsilon, x_max, x_match);
 // printf("fmin: %g fmax: %g gmin: %g gmax: %g\n", f_min, f_max, g_min, g_max);
  double tol = pow(10, -6); // Tolerance required for convergence
  while ((fabs(delta_f) > tol) || (fabs(delta_g) > tol) || (fabs(delta_df) > tol) || (fabs(delta_dg) > tol)) {
    printf("E: %g mu: %g B: %g a0: %g\n", e_tilde, mu, b_tilde, a0);
    delta_left_right_derivs(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, f_min, f_max, g_min, g_max, k, e_tilde, &delta_f, &delta_g, &delta_df, &delta_dg);
 
    printf("delta_f: %g delta_g: %g delta_df: %g delta_dg: %g\n", delta_f, delta_g, delta_df, delta_dg);
    double d_delta_f_d_b, d_delta_g_d_b, d_delta_f_d_a0, d_delta_g_d_a0;

    two_pt_b_derivative(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, k, b_tilde, l, a0, a1, &d_delta_f_d_b, &d_delta_g_d_b);
    two_pt_a0_derivative(df_hydrogen, dg_hydrogen, epsilon, x_max, x_match, k, e_tilde, l, a0, a1, &d_delta_f_d_a0, &d_delta_g_d_a0);

//    double c1 = d_delta_g_d_b*b_tilde + d_delta_g_d_a0*a0 - delta_g;
  //  double c2 = d_delta_f_d_b*b_tilde + d_delta_f_d_a0*a0 - delta_f;

    double det_m = d_delta_g_d_b*d_delta_f_d_a0 - d_delta_f_d_b*d_delta_g_d_a0;
//    b_tilde = 1.0/det_m*(c1*d_delta_f_d_a0 - c2*d_delta_g_d_a0);
    b_tilde -= 1.0/det_m*(d_delta_f_d_a0*delta_g - d_delta_g_d_a0*delta_f);
//    a0 = 1.0/det_m*(-c1*d_delta_f_d_b + c2*d_delta_g_d_b);
    a0 -= 1.0/det_m*(-d_delta_f_d_b*delta_g + d_delta_g_d_b*delta_f);
    e_tilde = 1.0 - b_tilde;
    mu = sqrt(1.0 - pow(e_tilde, 2.0));
    g_max = a1*exp(-mu*x_max);
    f_max = -a1*pow((1.0 - e_tilde)/(1.0 + e_tilde), 0.5)*exp(-mu*x_max);
    if (j == l + 0.5) {
      g_min = a0*pow(epsilon, l + 1);
      f_min = -a0*Z_ATOM*ALPHA_FS/(2.0*l + 3.0)*pow(epsilon, l + 1);
    } else if (j == l - 0.5) {
      f_min = a0*pow(epsilon, l);
      g_min = a0*Z_ATOM*ALPHA_FS/(2.0*l + 1.0)*pow(epsilon, l);
    }
  }
  double binding = b_tilde*M_ELECTRON*pow(10, 6);
  printf("B: %10.9lf a0: %g\n", binding, a0);
 
  return 0;
}
