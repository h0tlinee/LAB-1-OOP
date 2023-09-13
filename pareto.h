#pragma once
#include <iostream>
#include <time.h>
#include <math.h>
#include <fstream>

double randomizer();
double std_normal_distribution_plotnost(double x);
double std_normal_distribution_function(double x);
double param_normal_distribution_plotnost(double x, double mu, double lambda);
double std_pareto_normal_distribution(double x, double form_param);
double param_pareto_normal_distribution(double x, double mu, double lambda, double form_param);
double param_math_expect(double mu);
double std_pareto_normal_dispersion(double form_param);
double std_excess_koef(double form_param);
double param_pareto_normal_dispersion(double lambda, double form_param);
double central_interavl_P(double form_param);
double x1_realization(double form_param);
double pareto_normal_modulation(double form_param);
void out_to_file(double form_param);
