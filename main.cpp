#include "structures_and_functions.h"
#include <time.h>

int main() {
	setlocale(LC_ALL, "ru");
	srand((unsigned)time(0));
	mixed_pareto mix;
	mix.params1.v = 3.3;
	mix.params2.v = 2.3;
	mix.params1.lambda = 1;
	mix.params2.lambda = 1;
	mix.params1.mu = -3;
	mix.params2.mu = 3;
	mix.p = 0.7;

	/*auto x = generate_x_mixed(mix, 10000);
	auto result = mix_to_graph(10000, mix, x);
	ofstream file2("test2.txt");
	for (int i = 0; i < result.size(); i++) {
		file2 << result[i].first << "\t" << result[i].second << endl;
	}*/

	

	pareto params;
	params.v = 3.3;
	params.mu = 0;
	params.lambda = 1;
	auto x1 = generate_x_std(100000,params);
	auto result1 = std_to_graph(100000, params, x1);
	ofstream file1("test1.txt");
	for (int i = 0; i < result1.size(); i++) {
		file1 << result1[i].first << "\t" << result1[i].second << endl;
	}

	auto func4 = generate_empric_graph_selection(x1, 100000);
	ofstream file4("test3.txt");
	for (int i = 0; i < func4.size(); ++i) {
		file4 << func4[i].first << "\t" << func4[i].second << endl;
	}
}