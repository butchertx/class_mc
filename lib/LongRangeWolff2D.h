#pragma once
#include "IsingLattice2D.h"
#include "class_mc_io.h"
#include "Matrix.h"
#include "MemTimeTester.h"
#include <cmath>
extern "C" {
#include "random.h"
}


class LongRangeWolff2D {
	IsingLattice2D lat;
	class_mc_params params;
	Matrix interactions;
	IntMatrix cluster;//0 at indices that are not in cluster, 1 at indices that are
	void set_mean_field_model();
	void set_spin_boson_model();
	std::vector<spin> buffer;//holds the indices of spins to check
	bool LONG_RANGE_CLUSTER;

public:
	LongRangeWolff2D(IsingLattice2D*, class_mc_params*, std::string);

	void step();

	double calc_mag();

	double calc_E_mean_field_model();

	void set_beta(double new_beta) {
		//only use for parallel tempering
		params.beta = new_beta;
		set_mean_field_model();
	}

	double get_beta() {
		return params.beta;
	}

	void print_lat() {
		lat.print_lattice();
	}

	void test_spins(spin);

	void set_alg_long_range_cluster() { LONG_RANGE_CLUSTER = true; }

};
