#pragma once
#include <iostream>
#include <time.h>
#include <math.h>
#include <fstream>
#include "pareto.h"

double mixed_pareto_distribution(double p, double x, double mu1, double mu2, double l1, double l2, double form1, double form2);
double mixed_math_expect(double p, double mu1, double mu2);
double mixed_pareto_dispersion(double p, double mu1, double mu2, double lambda1, double lambda2, double form1, double form2);
double mixed_assymetry(double p, double mu1, double mu2, double lambda1, double lambda2, double form1, double form2);
double mixed_excess(double p, double mu1, double mu2, double lambda1, double lambda2, double form1, double form2);
double mixed_pareto_modulation(double p, double form1, double form2);
void mix_out_to_file(double p, double form1, double form2, double mu1, double mu2, double lambda1, double lambda2);
