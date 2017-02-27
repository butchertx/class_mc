#pragma once
#include "IsingLattice2D.h"
#include "class_mc_io.h"
#include "Matrix.h"
#include "MemTimeTester.h"
extern "C" {
#include "random.h"
}


class LongRangeWolff2D {
	IsingLattice2D lat;//underlying lattice
	class_mc_params params;//problem parameters
	Matrix interactions;//interaction matrix
	std::vector<double> cumulative_probs;//cumulative bond probabilities
	double E;//current "energy" (action divided by beta, probably corresponds to free energy)
	double dS;//double the running action due to cluster spin interactions
	int mag;//current magnetization


	std::vector<spin> buffer;//holds the indices of spins to check
	IntMatrix cluster;//0 at indices that are not in cluster, 1 at indices that are
	std::vector<spin> not_cluster;//vector of spins that have not been added to the cluster
	std::vector<std::vector<int>> cluster_alt;//alternative cluster that may or may not be faster
	int cluster_size;//size of the cluster during a given MC step
	int cluster_mag;//total magnetization of the cluster
	void fill_not_cluster();
	
	
	void set_model();//set the model based on the model parameter in "params"
	void set_mean_field_model();//create interaction matrix for the mean field model
	void set_spin_boson_model();//create interaction matrix for the spin boson model
	void set_cumulative_probs(int, int);//set cumulative probability vector for long range cluster forming


	bool LONG_RANGE_CLUSTER;//if true, use long range cluster forming
	bool NEAREST_NEIGHBOR_CLUSTER;//if true, use short range cluster forming
	bool TIMERS;//if true, timers will be used throughout
	MemTimeTester timers;
	bool VERIFY_TESTING;//use to print steps to check if algorithm is working properly

public:
	LongRangeWolff2D(IsingLattice2D*, class_mc_params*);

	void step();

	int get_cluster_size() { return cluster_size; }

	IsingLattice2D get_lat() { return lat; }

	double calc_mag();

	double calc_E_mean_field_model();

	double calc_E();

	std::vector<double> calc_corr(int dimension);

	std::vector<double> calc_corr_slow();

	double calc_sx();

	double calc_sz();

	void set_beta(double new_beta) {
		//only use for parallel tempering
		params.beta = new_beta;
		set_model();
	}

	double get_beta() {
		return params.beta;
	}

	double get_E() {
		return E;
	}

	double get_mag() {
		return mag < 0 ? -((double)mag)/lat.get_N() : ((double)mag) / lat.get_N();
	}

	void print_lat() {
		lat.print_lattice();
	}

	void print_interactions() {
		std::cout << interactions.to_string() << "\n";
	}

	void output_state();

	void test_spins(spin);

	void set_alg_long_range_cluster() { LONG_RANGE_CLUSTER = true; NEAREST_NEIGHBOR_CLUSTER = false; }
	void set_alg_short_range_cluster() { NEAREST_NEIGHBOR_CLUSTER = true; LONG_RANGE_CLUSTER = false; }
	void set_alg_naive_cluster() { LONG_RANGE_CLUSTER = false; NEAREST_NEIGHBOR_CLUSTER = false; }
	void set_timers() { TIMERS = true;}
	void set_testing() { VERIFY_TESTING = true; }

};
