#pragma once
#include <iostream>
#include <time.h>
#include <math.h>
#include "pareto.h"
#include <fstream>

const double PI = 3.141592653589793;

double randomizer() {
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}

double std_normal_distribution_plotnost(double x) {//вычисление плотности нормального распределения
	double plotnost = (1 / sqrt(2 * PI)) * (exp(-((pow(x, 2) / 2))));
	return plotnost;
}
double std_normal_distribution_function(double x) {//вычисление функции нормального распределения
	double distrib_f = (0.5 * (1 + erf(x / sqrt(2))));
	return distrib_f;
}
double param_normal_distribution_plotnost(double x, double mu, double lambda) {//вычисление плотности нормального распределения со сдвигом и масштабом
	double new_arg = (x - mu) / lambda;
	double plotnost = (1 / lambda) * std_normal_distribution_plotnost(x);
	return plotnost;
}
double std_pareto_normal_distribution(double x, double form_param) {//стандартное парето-нормальное распределение
	double k = 2 * std_normal_distribution_function(form_param) - 1 + (((2 * form_param) / (pow(form_param, 2) - 1)) * std_normal_distribution_plotnost(form_param));
	double plotnost;
	if (abs(x) <= form_param) {
		plotnost = (1 / k) * std_normal_distribution_plotnost(x);
	}
	else {
		plotnost = (1 / k) * std_normal_distribution_plotnost(x) * (pow((form_param / abs(x)), pow(form_param, 2)));
	}
	return plotnost;
}

double param_pareto_normal_distribution(double x, double mu, double lambda, double form_param) {//парето-нормальное распределение с учетом сдвига и масштаба
	double new_arg = (x - mu) / lambda;
	double plotnost = (1 / lambda) * std_pareto_normal_distribution(new_arg, form_param);
	return plotnost;
}

double std_math_expect = 0;//стд.матожидание=0, при сдвиге-увеличивается на мю
double param_math_expect(double mu) {
	return (std_math_expect + mu);
}
//коэф асимметрии всегда 0

double assymetry() {
	return 0;
}
double std_pareto_normal_dispersion(double form_param) {//стандартная дисперсия
	double k = 2 * std_normal_distribution_function(form_param) - 1 + (((2 * form_param) / (pow(form_param, 2) - 1)) * std_normal_distribution_plotnost(form_param));
	double dispersion;
	dispersion = 1 + (4 * pow(form_param, 3) * std_normal_distribution_plotnost(form_param)) / (k * (pow(form_param, 2) - 1) * (pow(form_param, 2) - 3));
	return dispersion;
}
double std_excess_koef(double form_param) {//стнадратный коэф экцесса , не меняется при сдвиге и масштабе
	double k = 2 * std_normal_distribution_function(form_param) - 1 + (((2 * form_param) / (pow(form_param, 2) - 1)) * std_normal_distribution_plotnost(form_param));
	double excess;
	double a = std_pareto_normal_dispersion(form_param);
	excess = (1 / (k * (pow(a, 2)))) * (3 * (2 * std_normal_distribution_function(form_param) - 1) - 2 * form_param * std_normal_distribution_plotnost(form_param) * ((3 + pow(form_param, 2) - (pow(form_param, 4)) / (pow(form_param, 2) - 5)))) - 3;
	return excess;
}

double param_pareto_normal_dispersion(double lambda, double form_param) {
	return pow(lambda, 2) * (std_pareto_normal_dispersion(form_param));
}

double central_interavl_P(double form_param) {
	double k = 2 * std_normal_distribution_function(form_param) - 1 + (((2 * form_param) / (pow(form_param, 2) - 1)) * std_normal_distribution_plotnost(form_param));
	double p = (2 * std_normal_distribution_function(form_param) - 1) / k;
	return p;
}

double x1_realization(double form_param) {//реализация x1
	double r2 = randomizer();
	double r3 = randomizer();
	//std::cout << r2 << " " << r3;
	double x1 = (sqrt(-2 * std::log(r2))) * cos(2 * PI * r3);
	return x1;
}

double pareto_normal_modulation(double form_param) {//моделирование случайной величины
	//шаг 1
	double x1;
	double x2;
	double r1 = randomizer();
	double p = central_interavl_P(form_param);
	if (r1 <= p) {
		//шаг 2,3
		do {
			x1=x1_realization(form_param);
		} while (x1 <= (-(form_param)) || (x1 >= form_param));
		return x1;
	}
	else {      //шаг 4
		double r4=randomizer();
		x2 = form_param / (pow(r4, (1 / (pow(form_param, 2) - 1))));
		if (r1 < (1 + p) / 2) {
			return x2;
		}
		else {
			return ( - (x2));
		}
	}
}

void out_to_file(double form_param) {
	std::ofstream file("1.txt");
	for (int i = 0; i < 1000000; i++) {
		double t = pareto_normal_modulation(form_param);
		file << t << "\t" << param_pareto_normal_distribution(t,1,0.3, form_param) << std::endl;
	};
}
