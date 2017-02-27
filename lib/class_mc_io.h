#pragma once
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "IsingLattice2D.h"

//classical monte carlo params for a number of different models and algorithms
//Models: n-vector models (ising, heisenberg), spinboson
//Lattices: chain, square, triangular, kagome
//Algorithms: wolff, generalizedwolff

const double PI = 3.1415897535;

struct spin_boson_params {
	double g, A0, delta, v;
};

struct class_mc_params {
	int dim, rand_seed, eq_time, steps_per_measure, measures_per_dump, max_dumps;
	double beta, h, kT;
	std::vector<int> lengths;
	std::vector<double> spacings;
	std::vector<double> Js;
	std::string model, lattice, alg;
	spin_boson_params sbparams;
	std::string to_string() {
		std::stringstream ss;
		ss << "Dimension: " << dim << "\n";
		ss << "Lattice type: " << lattice << "\n";
		ss << "Dimension lengths: ";
		for (auto const& value : lengths) {
			ss << value << " ";
		}
		ss << "\n";
		ss << "Dimension spacings: ";
		for (auto const& value : spacings) {
			ss << value << " ";
		}
		ss << "\n";
		ss << "Model: " << model << "\n";
		ss << "J couplings: ";
		for (auto const& value : Js) {
			ss << value << " ";
		}
		ss << "\n";
		ss << "kT: " << kT << "\n";
		ss << "Beta: " << beta << "\n";
		ss << "H: " << h << "\n";
		ss << "Algorithm: " << alg << "\n";
		ss << "Random Seed: " << rand_seed << "\n";
		ss << "Equilibration time: " << eq_time << "\n";
		ss << "Steps Per Measurement: " << steps_per_measure << "\n";
		ss << "Measures Per Dump: " << measures_per_dump << "\n";
		ss << "Max Dumps: " << max_dumps << "\n";

		if (model.compare("spin_boson") == 0) {
			ss << "Spin Boson Params:\n";
			ss << "g: " << sbparams.g << "\n";
			ss << "A0: " << sbparams.A0 << "\n";
			ss << "Delta: " << sbparams.delta << "\n";
			ss << "V: " << sbparams.v << "\n";
			ss << "Spatial Antiferro Coupling: " << 2.0 * sbparams.A0 * PI * PI / lengths[1] / lengths[1] / sinh(PI * spacings[0] / beta / sbparams.v) / sinh(PI * spacings[0] / beta / sbparams.v) << "\n";
		}

		return ss.str();
	}

	void set_q_ising_dummy(int length_in, double beta_in, double delta_in) {
		dim = 2;
		rand_seed = length_in;
		eq_time = 1;
		steps_per_measure = 0;
		measures_per_dump = 0;
		max_dumps = 0;
		beta = beta_in;
		h = 0.01;
		kT = 1 / beta;
		lengths = { length_in, length_in };
		spacings = { 1.0, 1.0 };
		Js = { 0.0, 1.0 };
		model = "spin_boson";
		lattice = "square";
		alg = "long_range_cluster";
		sbparams.A0 = 0.1;
		sbparams.delta = delta_in;
		sbparams.g = 1.0;
		sbparams.v = 1.0;
	}
};



struct class_mc_measurements {
	//vector for the names (steps, energies, mags, etc.)
	//parallel vector for the recorded values
	std::vector<std::string> names;
	std::vector<std::vector<double>> values;

	std::vector<std::string> func_names;
	std::vector<std::vector<double>> functions;
	std::vector<int> function_num_measures;

	std::vector<double> get_vals(std::string identifier) {
		int name_ind = -1;
		for (int i = 0; i < names.size(); ++i) {
			if (identifier.compare(names[i]) == 0) {
				name_ind = i;
			}
		}
		if (name_ind == -1) {
			std::cout << identifier << " is not a valid observable value\n";
			return{};
		}
		else {
			return values[name_ind];
		}
	}

	std::vector<double> get_func(std::string identifier) {
		int name_ind = -1;
		for (int i = 0; i < func_names.size(); ++i) {
			if (identifier.compare(func_names[i]) == 0) {
				name_ind = i;
			}
		}
		if (name_ind == -1) {
			std::cout << identifier << " is not a valid function observable\n";
			return{};
		}
		else {
			std::vector<double> result(functions[name_ind].size());
			for (int i = 0; i < result.size(); ++i) {
				result[i] = functions[name_ind][i] / function_num_measures[name_ind];
			}
			return result;
		}
	}

	void record(std::string name, double val) {
		bool found = false;
		for (int i = 0; i < names.size(); ++i) {
			if (name.compare(names[i]) == 0) {
				found = true;
				values[i].push_back(val);
			}
		}
		if (!found) {
			std::cout << name << " is not a valid observable\n";
		}
	}

	void record(std::string name, std::vector<double> function) {
		bool found = false;
		for (int i = 0; i < func_names.size(); ++i) {
			if (name.compare(func_names[i]) == 0) {
				found = true;
				if (function.size() == functions[i].size()) {
					for (int j = 0; j < functions[i].size(); ++j) {
						functions[i][j] += function[j];
					}
					function_num_measures[i] += 1;
				}
				else if (functions[i].size() == 0) {
					functions[i] = function;
					function_num_measures[i] = 1;
				}
				else { found = false; }
			}
		}
		if (!found) {
			std::cout << name << " is not a valid observable function, or there was a size mismatch between recordings\n";
		}
	}


	void record(std::vector<double> vals) {
		if (vals.size() != names.size()) {
			std::cout << "Attempt to record failed: too many/few values recorded\n";
		}
		else {
			for (int i = 0; i < vals.size(); ++i) {
				values[i].push_back(vals[i]);
			}
		}
	}
};

void create_input();

bool str_is_equal(std::string, std::string);

template<typename T>
std::string vec2str(std::vector<T> vec);

void read_input_ising(std::ifstream*, class_mc_params*);

void read_input_spin_boson(std::ifstream*, spin_boson_params*);

void apply_spin_boson_params(class_mc_params*);

void write_outputs(int, std::vector<int>, std::vector<double>, std::vector<double>, std::vector<double>);

void write_outputs_var(int, class_mc_measurements);

void write_final_outputs(class_mc_measurements results, int bins);

void write_params(class_mc_params params);

void write_state(int, IsingLattice2D&);

bool isDirExist(const std::string& path);

bool makePath(const std::string& path);

double mean(std::vector<double> vals);

double error(std::vector<double> vals, double mean, int bins);
