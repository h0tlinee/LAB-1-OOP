#include "structures_and_functions.h"


double delta_calc(const double n, const double min, const double max) {
	return (max - min) / ((int)log2(n) + 1);
}

vector<double> create_intervals(const double delta, const double min, const double max) {
	vector<double> intervals;
	double slider = min;
	intervals.push_back(slider);
	while (slider < max) {
		intervals.push_back(slider + delta);
		slider += delta;
	}
	return intervals;
}

int get_interval_index(const vector<double>& intervals, const double x) {
	if (intervals.size() < 2) {
		return 0;
	}
	if (x >= intervals[intervals.size() - 2] && x <= intervals[intervals.size() - 1]) {
		return intervals.size() - 2;
	}
	for (int i = 0; i < intervals.size() - 1; ++i) {
		if (x >= intervals[i] && x < intervals[i + 1]) {
			return i;
		}
	}
}

double empirc_density(const vector<double>& sample, const double x) {
	double delta = delta_calc(sample.size(), sample[0], sample[sample.size() - 1]);
	auto intervals = create_intervals(delta, sample[0], sample[sample.size() - 1]);
	int indx = get_interval_index(intervals, x);
	int left = 0;
	int right = sample.size() - 1;
	if (indx == 0) {
		right = 0;
		while (sample[right] < intervals[1]) {
			++right;
		}
	}
	else if (indx == intervals.size() - 2) {
		left = sample.size() - 1;
		while (sample[left] > intervals[intervals.size() - 2]) {
			--left;
		}
	}
	else {
		while (sample[left] < intervals[indx]) {
			++left;
		}
		right = left;
		while (sample[right] < intervals[indx + 1]) {
			++right;
		}
	}
	return (right - left + 1) / (sample.size() * delta);
}

double empirc_expected_value(const vector<double>& sample) {
	double sum = 0;
	for (int i = 0; i < sample.size(); ++i) {
		sum += sample[i];
	}
	return sum / sample.size();
}

double empirc_dispersion(const vector<double>& sample) {
	double sum = 0;
	double exp_val = empirc_expected_value(sample);
	for (int i = 0; i < sample.size(); ++i) {
		sum += pow(sample[i] - exp_val, 2);
	}
	return sum / sample.size();
}

double empirc_asymmetry(const vector<double>& sample) {
	double sum = 0;
	double exp_val = empirc_expected_value(sample);
	for (int i = 0; i < sample.size(); ++i) {
		sum += pow(sample[i] - exp_val, 3);
	}
	return sum / (sample.size() * pow(empirc_dispersion(sample), 3 / 2));
}

double empirc_kurtosis(const vector<double>& sample) {
	double sum = 0;
	double exp_val = empirc_expected_value(sample);
	for (int i = 0; i < sample.size(); ++i) {
		sum += pow(sample[i] - exp_val, 4);
	}
	return (sum / (sample.size() * pow(empirc_dispersion(sample), 2))) - 3;
}

vector<pair<double, double>> generate_empric_graph_selection(vector<double> selection, const int n) {
	vector<pair<double, double>> result;
	for (int i = 0; i < n; ++i) {
		result.push_back(make_pair(selection[i], empirc_density(selection, selection[i])));
	}
	return result;
}