#include "MemTimeTester.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "obs_calc_fast.cuh"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/count.h>
#include <thrust/inner_product.h>
#include <thrust/functional.h>
#include <string>


#include <stdio.h>


int main()
{

	MemTimeTester timer;
	thrust::host_vector<int> state;
	std::ifstream file;
	std::string line;
	int Lx, Ly;
/*
	//check single site energy, correlation function
	std::cout << "Calculating Single Site energy and correlation\n\n";
	file.open("C:/Users/Matthew/Dropbox/Code/class_mc/windows/ClassMCTest/x64/Debug/fast_calc_testing/single_site/dump/state0.csv");
	if (file.is_open()) {
		file >> Lx >> Ly;
		state.resize(Lx * Ly);
		for (int i = 0; i < Lx; ++i) {
			for (int j = 0; j < Ly; ++j) {
				file >> state[i*Ly + j];
			}
		}
	}
	else {
		std::cout << "Error: input file not opened\n";
	}
	timer.flag_start_time("single site correlation measurement");
	thrust::host_vector<int>& state_ref = state;
	thrust::host_vector<double> corr(Lx*Ly);
	thrust::host_vector<double>& corr_ref = corr;
	calc_corr_fast_1site(corr_ref, state_ref, Ly);
	timer.flag_end_time("single site correlation measurement");
	std::cout << "\nSingle Site Correlation function:\n";
	for (int i = 0; i < corr.size(); ++i) {
		std::cout << corr[i] << ",";
	}
	std::cout << "\n";
	file.close();

	thrust::host_vector<double> interactions(Lx*Ly);
	thrust::host_vector<double>& int_ref = interactions;
	file.open("C:/Users/Matthew/Dropbox/Code/class_mc/windows/ClassMCTest/x64/Debug/fast_calc_testing/single_site/interactions.csv");
	if (file.is_open()) {
		for (int i = 0; i < Lx*Ly; ++i) {
			file >> interactions[i];
		}
	}
	else {
		std::cout << "Interactions file failed to open.\n";
	}
	file.close();
	std::cout << "\nInteractions:\n";
	for (int i = 0; i < interactions.size(); ++i) {
		std::cout << interactions[i] << ",";
	}
	std::cout << "\n";
	timer.flag_start_time("single site energy calculation");
	std::cout << "Energy calculated: " << calc_action_fast(corr_ref, int_ref) << "\n";
	timer.flag_end_time("single site energy calculation");

	std::cout << "Sx: " << 0.5*(1 - corr[1]) << "\n";

	//check double site energy, correlation function
	std::cout << "Calculating Double Site energy and correlation\n\n";
	file.open("C:/Users/Matthew/Dropbox/Code/class_mc/windows/ClassMCTest/x64/Debug/fast_calc_testing/double_site/dump/state0.csv");
	if (file.is_open()) {
		file >> Lx >> Ly;
		state.resize(Lx * Ly);
		for (int i = 0; i < Lx; ++i) {
			for (int j = 0; j < Ly; ++j) {
				file >> state[i*Ly + j];
			}
		}
	}
	else {
		std::cout << "Error: input file not opened\n";
	}
	file.close();
	timer.flag_start_time("double site correlation measurement");
	state_ref = state;
	corr.resize(Lx * Ly);
	corr_ref = corr;
	std::cout << "Lx, Ly: " << Lx << ", " << Ly << "\n";
	std::cout << "state size: " << state.size() << "\n";
	std::cout << "correlation size: " << corr.size() << "\n";
	calc_corr_fast_2site(corr_ref, state_ref, Ly);
	timer.flag_end_time("double site correlation measurement");
	std::cout << "Double Site Correlation function:\n";
	for (int i = 0; i < corr.size(); ++i) {
		std::cout << corr[i] << ",";
	}
	std::cout << "\n";
	file.close();

	file.open("C:/Users/Matthew/Dropbox/Code/class_mc/windows/ClassMCTest/x64/Debug/fast_calc_testing/double_site/interactions.csv");
	interactions.resize(Lx*Ly);
	int_ref = interactions;
	if (file.is_open()) {
		for (int i = 0; i < Lx; ++i) {
			for (int j = 0; j < Ly; ++j) {
				file >> interactions[i*Ly + j];
			}
			std::getline(file, line);
		}
	}
	else {
		std::cout << "Interactions file failed to open.\n";
	}
	file.close();
	std::cout << "\nInteractions:\n";
	for (int i = 0; i < interactions.size(); ++i) {
		std::cout << interactions[i] << ",";
	}
	std::cout << "\n";
	timer.flag_start_time("double site energy calculation");
	std::cout << "Energy calculated: " << calc_action_fast(corr_ref, int_ref) << "\n";
	timer.flag_end_time("double site energy calculation");

	std::cout << "Sx: " << 0.5*(1 - corr[1]) << "\n";

	//check average energy, correlation function
	std::cout << "Calculating average of 3 runs energy and correlation\n\n";
	int nruns = 3;
	char filename[300];
	Lx = 1, Ly = 100;
	state.resize(Lx * Ly);
	state_ref = state;
	corr.resize(Lx * Ly);
	thrust::device_vector<double> corr_temp = corr;
	thrust::host_vector<double>& corr_temp_ref = corr;
	thrust::device_vector<double> corr_total(Lx*Ly, 0.0);
	for (int run = 0; run < nruns; ++run) {
		sprintf(filename, "C:/Users/Matthew/Dropbox/Code/class_mc/windows/ClassMCTest/x64/Debug/fast_calc_testing/avg_states/dump/state%d.csv", run);
		file.open(filename);
		if (file.is_open()) {
			file >> Lx >> Ly;
			for (int i = 0; i < Lx; ++i) {
				for (int j = 0; j < Ly; ++j) {
					file >> state[i*Ly + j];
				}
			}
		}
		else {
			std::cout << "Error: input file not opened\n";
		}
		file.close();
		timer.flag_start_time("double site correlation measurement");
		calc_corr_fast_1site(corr_temp_ref, state_ref, Ly);
		corr_temp = corr_temp_ref;
		thrust::transform(corr_temp.begin(), corr_temp.end(), corr_total.begin(), corr_total.begin(), thrust::plus<double>());
		timer.flag_end_time("double site correlation measurement");
	}
	thrust::constant_iterator<double> factor(1.0 / 3.0);
	thrust::transform(corr_total.begin(), corr_total.end(), factor, corr_temp.begin(), thrust::multiplies<double>());
	corr = corr_temp;
	std::cout << "Double Site Correlation function:\n";
	for (int i = 0; i < corr.size(); ++i) {
		std::cout << corr[i] << ",";
	}
	std::cout << "\n";

	/*file.open("C:/Users/Matthew/Dropbox/Code/class_mc/windows/ClassMCTest/x64/Debug/fast_calc_testing/double_site/interactions.csv");
	interactions.resize(Lx*Ly);
	int_ref = interactions;
	if (file.is_open()) {
		for (int i = 0; i < Lx; ++i) {
			for (int j = 0; j < Ly; ++j) {
				file >> interactions[i*Ly + j];
			}
			std::getline(file, line);
		}
	}
	else {
		std::cout << "Interactions file failed to open.\n";
	}
	file.close();
	std::cout << "\nInteractions:\n";
	for (int i = 0; i < interactions.size(); ++i) {
		std::cout << interactions[i] << ",";
	}
	std::cout << "\n";
	timer.flag_start_time("double site energy calculation");
	std::cout << "Energy calculated: " << calc_action_fast(corr_ref, int_ref) << "\n";
	timer.flag_end_time("double site energy calculation");

	std::cout << "Sx: " << 0.5*(1 - corr[1]) << "\n";
*/
	//Perform the resampling procedure using the states saved in classmctest
	std::cout << "Creating resampling data\n";
	double base_alpha = 1.1, temp_mag;
	int num_alphas = 11;
	Ly = 100;
	int num_runs = 10000;
	char filename[300];
	std::vector<double> alphas(num_alphas);
	std::vector<double> avg_mag2s(num_alphas);
	std::vector<double> state_sums(num_runs);
	std::vector<double> mag2s(num_runs);
	thrust::host_vector<double> resample_interactions(Ly);
	thrust::host_vector<double> temp_corr(Ly);
	thrust::host_vector<double>& corr_ref = temp_corr;
	thrust::host_vector<double>& int_ref = resample_interactions;
	thrust::host_vector<int>& state_ref = state;
	file.open("C:/Users/Matthew/Dropbox/Code/class_mc/windows/ClassMCTest/x64/Debug/fast_calc_testing/resampling/alphas.txt");
	if (file.is_open()) {
		for (int i = 0; i < num_alphas; ++i) {
			file >> alphas[i];
		}
	}
	else {
		std::cout << "Error: alpha file not opened\n";
	}
	file.close();

	file.open("C:/Users/Matthew/Dropbox/Code/class_mc/windows/ClassMCTest/x64/Debug/fast_calc_testing/resampling/resample_interactions.csv");
	if (file.is_open()) {
		for (int i = 0; i < Ly; ++i) {
			file >> resample_interactions[i];
		}
	}
	else {
		std::cout << "Error: interactions file not opened\n";
	}
	file.close();

	for (int run = 0; run < num_runs; ++run) {
		sprintf(filename, "C:/Users/Matthew/Dropbox/Code/class_mc/windows/ClassMCTest/x64/Debug/fast_calc_testing/resampling/state%d.csv", run);
		file.open(filename);
		if (file.is_open()) {
			file >> Lx >> Ly;
			state.resize(Lx * Ly);
			for (int i = 0; i < Lx; ++i) {
				for (int j = 0; j < Ly; ++j) {
					file >> state[i*Ly + j];
				}
			}
		}
		else {
			std::cout << "Error: state file not opened\n";
		}

		calc_corr_fast_1site(corr_ref, state_ref, Ly);
		state_sums[run] = calc_action_fast(corr_ref, int_ref);
		temp_mag = thrust::inner_product(state.begin(), state.end(), thrust::make_constant_iterator(1.0), 0.0);
		mag2s[run] = temp_mag*temp_mag/Ly/Ly;
		file.close();
	}
	//redo averages with different parameters
	double m2avg = 0;
	for (int i = 0; i < mag2s.size(); ++i) {
		m2avg += mag2s[i];
	}
	std::cout << "avg mag2 for alpha = 1: " << m2avg/mag2s.size() << "\n";

	//find the initial histogram - assume relatively gaussian
	//use the max, min, and mean
	double max = state_sums[0], min = state_sums[0], sum_temp = 0;
	for (int i = 0; i < state_sums.size(); ++i) {
		sum_temp = state_sums[i];
		if (max < sum_temp) { max = sum_temp; }
		if (min > sum_temp) { min = sum_temp; }
	}
	int num_bins = state_sums.size() / 10;
	double bin_step = (max - min) / (num_bins - 1);
	std::vector<int> histogram(num_bins);
	int bin;
	for (int i = 0; i < state_sums.size(); ++i) {
		bin = (int)((state_sums[i] - min) / bin_step);
		if (bin >= 0 && bin < num_bins) {
			histogram[bin] += 1;
		}
		else {
			std::cout << "Error: bin index out of bounds\n";
		}
	}
	for (int i = 0; i < num_bins; ++i) {
		std::cout << histogram[i] << ",";
	}
	std::cout << "\n";

	//resample using other alphas
	double part_func, temp_prob;
	for (int i = 0; i < alphas.size(); ++i) {
		temp_mag = 0.0;
		part_func = 0.0;
		//find the approximate partition function for the new value of alpha
		for (int s = 0; s < num_runs; ++s) {
			bin = (int)((state_sums[s] - min) / bin_step);
			temp_prob = (1.0*histogram[bin]) / num_runs;
			part_func += temp_prob*exp(-(alphas[i] - base_alpha) * state_sums[s]);
		}
		//find the average mag^2 with the new probability distribution
		for (int s = 0; s < num_runs; ++s) {
			bin = (int)((state_sums[s] - min) / bin_step);
			temp_prob = (1.0*histogram[bin]) / num_runs;
			temp_mag += mag2s[s] * temp_prob * exp(-(alphas[i] - base_alpha) * state_sums[s]) / part_func;
		}
		std::cout << "partition function: " << part_func << " and mag2: " << temp_mag << "\n";
		avg_mag2s[i] = temp_mag;
	}
	std::ofstream outfile;
	outfile.open("C:/Users/Matthew/Dropbox/Code/class_mc/windows/ClassMCTest/x64/Debug/fast_calc_testing/resampling/resample_results.csv");
	for (int i = 0; i < alphas.size(); ++i) {
		outfile << alphas[i] << ",";
	}
	outfile << "\n";
	for (int i = 0; i < alphas.size(); ++i) {
		outfile << avg_mag2s[i] << ",";
	}
	outfile.close();

	timer.print_timers();
    return 0;
}
