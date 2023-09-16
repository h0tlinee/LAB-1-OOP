#include "structures_and_functions.h"

double mixed_plotnost(const double x,const mixed_pareto& mix) {
	return (1 - mix.p) * pareto_normal_distribution(x,mix.params1) + mix.p * pareto_normal_distribution(x,mix.params2);
}

double mix_math_expect(const mixed_pareto& mix ) {
	return (1 - mix.p) * math_expect(mix.params1) + mix.p * math_expect(mix.params2);
}

double mixed_dispersion(const mixed_pareto& mix) {
	return (1 - mix.p) * (pow(math_expect(mix.params1), 2) + dispersion(mix.params1)) +
		mix.p * (pow(math_expect(mix.params2), 2) + dispersion(mix.params2)) -
		pow(mix_math_expect(mix), 2);
}

double mixed_asymmetry(const mixed_pareto& mix) {
	return ((1 - mix.p) * (pow((math_expect(mix.params1) - mix_math_expect(mix)), 3) + 3 * (math_expect(mix.params1) - mix_math_expect(mix)) * dispersion(mix.params1) + pow(dispersion(mix.params1), 3 / 2) * asymmetry_k(mix.params1)) +
		mix.p * (pow((math_expect(mix.params2) - mix_math_expect(mix)), 3) + 3 * (math_expect(mix.params2) - mix_math_expect(mix)) * dispersion(mix.params2) + pow(dispersion(mix.params2), 3 / 2) * asymmetry_k(mix.params2))) /
		pow(mixed_dispersion(mix), 3 / 2);
}

double mixed_excess(const mixed_pareto& mix) {
	return ((1 - mix.p) * (pow((math_expect(mix.params1) - mix_math_expect(mix)), 4) + 6 * dispersion(mix.params1) * pow((math_expect(mix.params1) - mix_math_expect(mix)), 2) +
		4 * (math_expect(mix.params1) - mix_math_expect(mix)) * pow(dispersion(mix.params1), 3 / 2) * asymmetry_k(mix.params1) + pow(dispersion(mix.params1), 2) * excess_koef(mix.params1)) +
		mix.p * (pow((math_expect(mix.params2) - mix_math_expect(mix)), 4) + 6 * dispersion(mix.params2) * pow((math_expect(mix.params2) - mix_math_expect(mix)), 2) +
			4 * (math_expect(mix.params2) - mix_math_expect(mix)) * pow(dispersion(mix.params2), 3 / 2) * asymmetry_k(mix.params2) + pow(dispersion(mix.params2), 2) * excess_koef(mix.params2)) - 3) /
		pow(mixed_dispersion(mix), 2);
}

double mixed_pareto_modeling(const mixed_pareto& mix) {
	double r = randomizer();
	if (r > mix.p) {
		return pareto_modeling(mix.params1);
	}
	else {
		return pareto_modeling(mix.params2);
	}
}

vector<double> generate_x_mixed(const mixed_pareto& mix, int n) {
	vector<double> x;
	for (int i = 0; i < n; ++i) {
		x.push_back(mixed_pareto_modeling(mix));
	}
	sort(x.begin(), x.end());
	return x;
}

vector<pair<double, double>> mix_to_graph(int n, const mixed_pareto& mix, vector<double> x) {
	vector<pair<double, double>> result;
	for (int i = 0; i < n; i++) {
		result.push_back(make_pair(x[i], mixed_plotnost(x[i], mix)));
	}
	return result;
}