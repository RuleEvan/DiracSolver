#include "potential.h"
double a_coeff(double n, double l, double k) {
  double a = pow(-1.0, k)/gsl_sf_gamma(k + 1.0);
  a *= sqrt(2.0*gsl_sf_gamma(n + 1.0)/gsl_sf_gamma(n + l + 1.5))*gsl_sf_gamma(n + l + 1.5)/(gsl_sf_gamma(n - k + 1.0)*gsl_sf_gamma(k + l + 1.5));
  
  return a;
}

double b_coeff(double n, double l, double np, double lp, double p) {
  double b = 0.5*gsl_sf_gamma(p + 1.5);
  int k, kp;
  double b_sum = 0.0;
  for (k = 0; k <= n; k++) {
    for (kp = 0; kp <= np; kp++) {
      if (0.5*(2*k + 2*kp + l + lp) != p) {continue;}
      b_sum += a_coeff(n,l,k)*a_coeff(np, lp, kp);
    }
  }
  b *= b_sum;

  return b;
}

double compute_potential(double n, double np, double l, double lp, int l_bessel, double k, int r_pow) {
  // Sum the required Talmi integrals with the corresponding B coefficients
  double v = 0;
  for (int ip = l + lp; ip <= l + lp + 2*n + 2*np; ip += 2) {
    double p = ip/2.0;
    v += b_coeff(n, l, np, lp, p)*talmi(p, l_bessel, k, r_pow);
  }
  return v;
}

double talmi(double p, int l, double k, int r_pow) {
  // Compute the order p Talmi integral
  // Set the limits of the integral and the error tolerance
  double r_min = 0.0001;
  double r_max = 10.0;
  double tol = pow(10, -6);
  double I_p = Romberg5Vars(&talmi_integrand, r_min, r_max, p, l, k, r_pow, tol);
  I_p *= 2.0/gsl_sf_gamma(p + 1.5);
  
  return I_p;
}
    
double talmi_integrand(double q, double p, int l, double k, int r_pow) {
  // The integrand of the Talmi integral
  // Plug in the required potential here
  double v = pow(q, 2.0*p + 2.0)*exp(-q*q);
  if (COR_FAC == 1) {
    double beta = exp(-1.1*pow(B_OSC*q, 2))*(1.0 - 0.68*pow(B_OSC*q,2.0));
    v *= pow(1.0 - beta, 2.0);
  }
  v *= gsl_sf_bessel_jl(l, k*q*B_OSC)*pow(q, r_pow);

  return v;
}
