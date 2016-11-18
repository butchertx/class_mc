#pragma once
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include "IsingLattice2D.h"

//classical monte carlo params for a number of different models and algorithms
//Models: n-vector models (ising, heisenberg), spinboson
//Lattices: chain, square, triangular, kagome
//Algorithms: wolff, generalizedwolff
struct class_mc_params {
	int dim, rand_seed, eq_time, steps_per_measure, measures_per_dump, max_dumps;
	double beta, h, kT;
	std::vector<int> lengths;
	std::vector<double> spacings;
	std::vector<double> Js;
	std::string model, lattice, alg;
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
		return ss.str();
	}
};

struct spin_boson_params {
	double g, A0, delta;
};

struct class_mc_measurements {
	//vector for the names (steps, energies, mags, etc.)
	//parallel vector for the recorded values
	std::vector<std::string> names;
	std::vector<std::vector<double>> values;

	std::vector<double> get_vals(std::string identifier) {
		int name_ind = -1;
		for (int i = 0; i < names.size(); ++i) {
			if (identifier.compare(names[i]) == 0) {
				name_ind = i;
			}
		}
		if (name_ind == -1) {
			std::cout << identifier << " is not a valid observable\n";
			return{};
		}
		else {
			return values[name_ind];
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

void write_outputs(int, std::vector<int>, std::vector<double>, std::vector<double>, std::vector<double>);

void write_state(int, IsingLattice2D, double);

bool isDirExist(const std::string& path);

bool makePath(const std::string& path);