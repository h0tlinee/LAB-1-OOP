#include "structures_and_functions.h"

double new_argument(const double x, const double lambda, const double mu) {//пересчет аргумента для сдвига
	return (x - mu) / lambda;
}

double randomizer() {//гнератор случайных чисел
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}


double normal_distribution_function(double v) {//вычисление функции нормального распределения, в нашем варианте зависима только от v
	double distrib_f = (0.5 * (1 + erf(v / sqrt(2))));
	return distrib_f;
}

double normal_distribution_plotnost(double x, const pareto& params) {//вычисление плотности нормального распределения со сдвигом и масштабом, и от v и от x
	double new_arg = new_argument(x, params.lambda, params.mu);
	double plotnost = ((1 / sqrt(2 * PI)) * (exp(-((pow(new_arg, 2) / 2)))))/params.lambda;
	return plotnost;
}

double std_normal_distribution_plotnost(double x) {//вычисление плотности нормального распределения
	double plotnost = (1 / sqrt(2 * PI)) * (exp(-((pow(x, 2) / 2))));
	return plotnost;
}

double k_koef(const pareto& params) {
	return 2 * normal_distribution_function(params.v) - 1 + (2 * params.v * std_normal_distribution_plotnost(params.v) / (pow(params.v, 2) - 1));
}
double pareto_normal_distribution(double x, const pareto& params) {//плотность парето-нормального распределения
	double new_arg = new_argument(x, params.lambda, params.mu);
	if (abs(x) <= params.v) {
		return ((1 / k_koef(params)) * normal_distribution_plotnost(x, params)) / params.lambda;
	}
	else {
		return ((1 / k_koef(params))* std_normal_distribution_plotnost(params.v)* (pow((params.v / abs(x)), pow(params.v, 2))))/params.lambda;
	}
}

double math_expect(const pareto& params) { //мат ожидание для станд. распред. равно сдвигу
	return params.mu;
}

double asymmetry_k(const pareto& params) {//для стд. распред. коэф. асимметрии 0
	return 0;
}

double dispersion(const pareto& params) { //увеличивается в lambda^2 при сдвиге
	return (1 + (4 * pow(params.v, 3) * std_normal_distribution_plotnost(params.v)) / (k_koef(params) * (pow(params.v, 2) - 1) * (pow(params.v, 2) - 3))) * pow(params.lambda,2);
}

double excess_koef(const pareto& params) {//не меняется
	return (1 / (k_koef(params) * (pow(dispersion(params), 2)))) * (3 * (2 * normal_distribution_function(params.v) - 1) - 2 * params.v * std_normal_distribution_plotnost(params.v) * ((3 + pow(params.v, 2) - (pow(params.v, 4)) / (pow(params.v, 2) - 5)))) - 3;
}

double central_interval(const pareto& params) {
	return (2 * normal_distribution_function(params.v) - 1) / k_koef(params);
}

double pareto_modeling(const pareto& params) {
	double x1;
	double r1 = randomizer();
	//cout << "r1   " << r1 << endl;
	double p = central_interval(params);
	if (r1 <= p) {//шаг1
		do {//шаг 2,3
			double r2 = randomizer();
			double r3 = randomizer();
			//cout << "r2   " << r2 << endl;
			//cout << "r3   " << r3 << endl;
			x1 = sqrt((-2) * log(r2)) * cos(2 * PI * r3);
		} while (x1 <= (-(params.v)) || (x1 >= params.v));
		return x1;
	}
	else {//шаг 4 
		double r4 = randomizer();
		//cout << "r4   " << r4 << endl;
		double x2 = params.v / (pow(r4, (1 / (pow(params.v, 2) - 1))));
		if (r1 < ((1 + central_interval(params)) / 2)) {
			return x2;
		}
		else {
			return ( - (x2));
		}
	}
}



vector<double> generate_x_std(int n,const pareto& params){
	vector <double> x;
	for (int i = 0; i < n; i++) {
		x.push_back(pareto_modeling(params));
		//cout << x[i] << endl;
	}
	sort(x.begin(), x.end());
	return x;
}

vector<pair<double, double>> std_to_graph(int n, const pareto& params, vector<double> x) {
	vector<pair<double, double>> result;
	for (int i = 0; i < n; i++) {
		result.push_back(make_pair(x[i], pareto_normal_distribution(x[i], params)));
	}
	return result;
}