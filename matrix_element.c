#include "matrix_element.h"

double decay_rate() {

  FILE* density_file;
  density_file = fopen("al27-al27_J2_T1_1_0.dens", "r");
  double t_tot = 0.0;
  double delta_E = 0.88241; // [MeV]
  double k = delta_E/(197.3); // [fm]^-1 
  double ji = 0.5;
  double omega = 8.0*M_PI*ALPHA_FS*k*3.0*pow(10, 8)*pow(10, 12)/(2*ji + 1);
  double lambda_p = 2.793;
  double lambda_n = -1.913;
  printf("k: %g b: %g\n", k, B_OSC);
  while (0 == 0) {
    int j_op, n_entries;
    int ret = fscanf(density_file, "%d, %d\n", &j_op, &n_entries);
    printf("J_op: %d\n", j_op);
    if (ret == EOF) {break;}
    double t_el = 0;
    double t_mag = 0;
    for (int i = 0; i < n_entries; i++) {
      int np, lp, ijp, n, l, ij;
      double jp, j, coeff_t0, coeff_t1;
      fscanf(density_file, "%d, %d, %d, %d, %d, %d, %lf, %lf\n", &np, &lp, &ijp, &n, &l, &ij, &coeff_t1, &coeff_t0);
  //    printf("%g, %g\n", coeff_t0, coeff_t1);
      jp = ijp/2.0;
      j = ij/2.0;
  //    printf("%d %d %g %d %d %g\n", np, lp, jp, n, l, j);
      t_el += 0.5*(coeff_t0*sqrt(2) + coeff_t1*sqrt(6))*t_el_c(np, lp, jp, n, l, j, j_op, k);
      t_el += 0.5*(coeff_t0*(lambda_p + lambda_n) + coeff_t1*(lambda_p - lambda_n))*t_el_mu(np, lp, jp, n, l, j, j_op, k);
      t_mag += 0.5*(coeff_t0*sqrt(2) + coeff_t1*sqrt(6))*t_mag_c(np, lp, jp, n, l, j, j_op, k);
      t_mag += 0.5*(coeff_t0*(lambda_p + lambda_n) + coeff_t1*(lambda_p - lambda_n))*t_mag_mu(np, lp, jp, n, l, j, j_op, k);
  //    printf("El: %g Mag: %g\n", t_el, t_mag);
    }
    t_mag *= t_mag;
    t_el *= t_el;
    t_tot += t_el + t_mag;
    printf("T: %g\n", t_tot);
  }
  omega *= t_tot;

  return omega;
}

double t_el_c(int np, int lp, double jp, int n, int l, double j, int j_op, double k) {

  double t_el_1 = em_dot_grad_mat(np, lp, jp, n, l, j, j_op, j_op + 1, k);  
  double t_el_2 = em_dot_grad_mat(np, lp, jp, n, l, j, j_op, j_op - 1, k);
  double jx = j_op;
//  printf("j_op: %d tel1: %g tel2: %g\n", j_op, t_el_1, t_el_2);
  double t_tot = 197.3/(M_NEUTRON)*(-sqrt(jx/(2*jx + 1))*t_el_1 + sqrt((jx + 1)/(2*jx + 1))*t_el_2);
//  double t_tot = 1.0/(k)*(-sqrt(jx/(2*jx + 1))*t_el_1 + sqrt((jx + 1)/(2*jx + 1))*t_el_2);

  printf("T_el_c: %d %d %g %d %d %g %g\n", np, lp, jp, n, l, j, t_tot);
  return t_tot;
}

double t_el_mu(int np, int lp, double jp, int n, int l, double j, int j_op, double k) {
//  double t_el = em_dot_sig_mat(np, lp, jp, n, l, j, j_op, j_op, k);

  double t_el = k*197.3/(2.0*M_NEUTRON)*em_dot_sig_mat(np, lp, jp, n, l, j, j_op, j_op, k);
  printf("T_el_mu: %d %d %g %d %d %g %g\n", np, lp, jp, n, l, j, t_el); 
  return t_el;
}

double t_mag_c(int np, int lp, double jp, int n, int l, double j, int j_op, double k) {
  
  double t_mag = -197.3/(M_NEUTRON)*em_dot_grad_mat(np, lp, jp, n, l, j, j_op, j_op, k);
  return t_mag;
}

double t_mag_mu(int np, int lp, double jp, int n, int l, double j, int j_op, double k) {

  double t_mag_1 = em_dot_sig_mat(np, lp, jp, n, l, j, j_op, j_op + 1, k);
  double t_mag_2 = em_dot_sig_mat(np, lp, jp, n, l, j, j_op, j_op - 1, k);
  double jx = j_op;
  double t_tot = 197.3*k/(2.0*M_NEUTRON)*(-sqrt(jx/(2*jx + 1))*t_mag_1 + sqrt((jx + 1)/(2*jx + 1))*t_mag_2);

  return t_tot;
}

double dpot1(int np, int lp, int n, int l, int l_op, double k) {
  double mat = -compute_potential(n, np, l, lp, l_op, k, 1);
  if (n == 0) {return mat;}
  mat += - 2.0*sqrt(n)*compute_potential(n - 1, np, l + 1, lp, l_op, k, 0);
//  printf("Dpot1: %g\n", mat);
  return mat;
}

double dpot2(int np, int lp, int n, int l, int l_op, double k) {
  double mat = (2*l + 1)*compute_potential(n, np, l, lp, l_op, k, -1) - compute_potential(n, np, l, lp, l_op, k, 1);
  if (n == 0) {return mat;}
  mat += - 2.0*sqrt(n)*compute_potential(n - 1, np, l + 1, lp, l_op, k, 0);
//  printf("Dpot 2: %g\n", mat);
  return mat;
}

double em_mat(int np, int lp, double jp, int n, int l, double j, int j_op, double k) {
  double mat = pow(-1.0, j + j_op + 0.5)/sqrt(4*M_PI)*sqrt((2*lp + 1)*(2*l + 1)*(2*jp + 1)*(2*j + 1))*sqrt(2*j_op + 1)*six_j(lp, jp, 0.5, j, l, j_op)*three_j(lp, j_op, l, 0, 0, 0)*compute_potential(n, np, l, lp, j_op, k, 0);

  return mat;
}

double em_dot_sig_mat(int np, int lp, double jp, int n, int l, double j, int j_op, int l_op, double k) {
  double mat = pow(-1.0, lp)*sqrt(6/(4*M_PI))*sqrt((2*lp + 1)*(2*l + 1)*(2*jp + 1)*(2*j + 1))*sqrt((2*j_op + 1)*(2*l_op + 1))*nine_j(lp, l, l_op, 0.5, 0.5, 1, jp, j, j_op)*three_j(lp, l_op, l, 0, 0, 0)*compute_potential(n, np, l, lp, l_op, k, 0);
//  printf("emdot sig: %d %d %g %d %d %g %d %g\n", np, lp, jp, n, l, j, l_op, mat);
  return mat;
}

double em_dot_grad_mat(int np, int lp, double jp, int n, int l, double j, int j_op, int l_op, double k) {
  double mat = 1.0/B_OSC*pow(-1.0, lp + j - 0.5)/sqrt(4*M_PI)*sqrt((2*lp + 1)*(2*l + 1)*(2*jp + 1)*(2*j + 1))*sqrt((2*l_op + 1)*(2*j_op + 1))*six_j(lp, jp, 0.5, j, l, j_op);
  double m1 = six_j(l_op, 1, j_op, l, lp, l + 1)*three_j(lp, l_op, l + 1, 0, 0, 0)/three_j(l + 1, 1, l, 0, 0, 0)*(l + 1)/(2*l + 1)*dpot1(np, lp, n, l, l_op, k);
  if (l == 0) {return mat*m1;}
  double m2 = six_j(l_op, 1, j_op, l, lp, l - 1)*three_j(lp, l_op, l - 1, 0, 0, 0)/three_j(l - 1, 1, l, 0, 0, 0)*l/(2*l + 1)*dpot2(np, lp, n, l, l_op, k);
  mat = mat*(m1 + m2);
//  printf("emdot grad: %d %d %g %d %d %g %d %d %g\n", np, lp, jp, n, l, j,j_op, l_op, m1); 
  return mat;
}
