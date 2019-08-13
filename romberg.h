#ifndef ROMBERG_H
#define ROMBERG_H
#include "angular.h"
double RombergIntegrator(double (*f)(double), double a, double b, double tol); 
double Romberg2Vars(double (*f)(double, double), double a, double b, double r, double tol);
double Romberg3Vars(double (*f)(double, double, int), double a, double b, double p, int iv, double tol);
double Romberg5Vars(double (*f)(double, double, int, double, int), double a, double b, double p, int l, double k, int r_pow, double tol);

#endif
