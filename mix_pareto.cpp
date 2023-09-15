#include <iostream>
#include <time.h>
#include <math.h>
#include "pareto.h"
#include <fstream>
#include "mix_pareto.h"

const double PI = 3.141592653589793;

double mixed_pareto_distribution(double p, double x, double mu1, double mu2, double l1, double l2, double form1, double form2) {
	return ((1 - p) * param_pareto_normal_distribution(x, mu1, l1, form1)) + (p* param_pareto_normal_distribution(x, mu2, l2, form2));
}

double mixed_math_expect(double p, double mu1, double mu2 ) {
	return ((1 - p) * param_math_expect(mu1)) + (p * param_math_expect(mu2));
}

double mixed_pareto_dispersion(double p, double mu1, double mu2,double lambda1, double lambda2, double form1, double form2) {
	return (1 - p) * (pow(param_math_expect(mu1), 2) + param_pareto_normal_dispersion(lambda1, form1)) + p * (pow(param_math_expect(mu2), 2) + param_pareto_normal_dispersion(lambda2, form2))
		- pow(mixed_math_expect(p, mu1, mu2), 2);
}

double mixed_assymetry(double p,double mu1, double mu2, double lambda1, double lambda2, double form1, double form2){
	return ((1 - p) * (pow((param_math_expect(mu1) - mixed_math_expect(p,mu1,mu2)), 3) + 3 * (param_math_expect(mu1) - mixed_math_expect(p,mu1,mu2)) * param_pareto_normal_dispersion(lambda1,form1) + pow(param_pareto_normal_dispersion(lambda1,form1), 3 / 2) * assymetry()) +
		p * (pow((param_math_expect(mu2) - mixed_math_expect(p, mu1, mu2)), 3) + 3 * (param_math_expect(mu2) - mixed_math_expect(p, mu1, mu2)) * param_pareto_normal_dispersion(lambda2, form2) + pow(param_pareto_normal_dispersion(lambda2, form2), 3 / 2) * assymetry())) /
		pow(mixed_pareto_dispersion(p,mu1,mu2,lambda1,lambda2,form1,form2), 3 / 2);
}

double mixed_excess(double p, double mu1, double mu2, double lambda1, double lambda2, double form1, double form2) {
	return ((1 - p) * (pow((param_math_expect(mu1) - mixed_math_expect(p, mu1, mu2)), 4) + 6 * param_pareto_normal_dispersion(lambda1, form1) * pow((param_math_expect(mu1) - mixed_math_expect(p, mu1, mu2)), 2) +
		4 * (param_math_expect(mu1) - mixed_math_expect(p, mu1, mu2)) * pow(param_pareto_normal_dispersion(lambda1, form1), 3 / 2) * assymetry() + pow(param_pareto_normal_dispersion(lambda1, form1), 2) * std_excess_koef(form1)) +
		p * (pow((param_math_expect(mu2) - mixed_math_expect(p, mu1, mu2)), 4) + 6 * param_pareto_normal_dispersion(lambda2, form2) * pow((param_math_expect(mu2) - mixed_math_expect(p, mu1, mu2)), 2) +
			4 * (param_math_expect(mu2) - mixed_math_expect(p, mu1, mu2)) * pow(param_pareto_normal_dispersion(lambda2, form2), 3 / 2) * assymetry() + pow(param_pareto_normal_dispersion(lambda2, form2), 2) * std_excess_koef(form2)) - 3) /
		pow(mixed_pareto_dispersion(p, mu1, mu2, lambda1, lambda2, form1, form2), 2);
}

double mixed_pareto_modulation(double p,double form1, double form2) {
	double r = randomizer();
	if (r > p) {
		return pareto_normal_modulation(form1);
	}
	else {
		return pareto_normal_modulation(form2);
	}
}
void mix_out_to_file(double p, double form1, double form2,double mu1, double mu2, double lambda1,double lambda2) {
	std::ofstream file("2.txt");
	for (int i = 0; i < 1000000; i++) {
		double t = mixed_pareto_modulation(p,form1,form2);
		file << t << "\t" << mixed_pareto_distribution(p,t,mu1,mu2,lambda1,lambda2,form1,form2) << std::endl;
	};
}