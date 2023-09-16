#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
const double PI = 3.141592653589793;
using namespace std;

struct pareto {
	double v = 3.3;
	double mu = 0;
	double lambda = 0.3;
};

struct mixed_pareto {
	pareto params1;
	pareto params2;
	double p = 0.5;
};

double new_argument(const double x, const double lambda, const double mu);
double randomizer();
double normal_distribution_function(double v);
double normal_distribution_plotnost(double x, const pareto& params);
double std_normal_distribution_plotnost(double x);
double k_koef(const pareto& params);
double pareto_normal_distribution(double x, const pareto& params);
double math_expect(const pareto& params);
double asymmetry_k(const pareto& params);
double dispersion(const pareto& params);
double excess_koef(const pareto& params);
double central_interval(const pareto& params);
double pareto_modeling(const pareto& params);
vector<double> generate_x_std(int n, const pareto& params);
vector<pair<double, double>> std_to_graph(int n, const pareto& params, vector<double> x);


double mixed_plotnost(const double x, const mixed_pareto& mix);
double mix_math_expect(const mixed_pareto& mix);
double mixed_dispersion(const mixed_pareto& mix);
double mixed_asymmetry(const mixed_pareto& mix);
double mixed_excess(const mixed_pareto& mix);
double mixed_pareto_modeling(const mixed_pareto& mix);
vector<double> generate_x_mixed(const mixed_pareto& mix, int n);
vector<pair<double, double>> mix_to_graph(int n, const mixed_pareto& mix, vector<double> x);

double delta_calc(const double n, const double min, const double max);
vector<double> create_intervals(const double delta, const double min, const double max);
int get_interval_index(const vector<double>& intervals, const double x);
double empirc_density(const vector<double>& sample, const double x);
double empirc_expected_value(const vector<double>& sample);
double empirc_dispersion(const vector<double>& sample);
double empirc_asymmetry(const vector<double>& sample);
double empirc_kurtosis(const vector<double>& sample);
vector<pair<double, double>> generate_empric_graph_selection(vector<double> selection, const int n);