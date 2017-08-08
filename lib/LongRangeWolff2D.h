#pragma once
#include <thrust/host_vector.h>
#include <sstream>
#include "IsingLattice2D.h"
#include "class_mc_io.h"
#include "Matrix.h"
#include "MemTimeTester.h"
extern "C" {
#include "random.h"
}
#include <gsl/gsl_sf_zeta.h>
#include <cmath>
#include "obs_calc_fast.cuh"

double real_trigamma_cmplxarg(double x, double y, int m, int l);


class LongRangeWolff2D {
	IsingLattice2D lat;//underlying lattice
	class_mc_params params;//problem parameters
	Matrix interactions;//interaction matrix
	//std::vector<double> cumulative_probs;//cumulative bond probabilities
	Matrix interaction_sum;//values from interaction sum integral
	double insum;
	double outsum;//helper sums for the interaction sum and for finding lattice interactions
	double E;//current "energy" (action divided by beta, probably corresponds to free energy)
	double dS;//double the running action due to cluster spin interactions
	int mag;//current magnetization
	int tau_x;//time index of last antiferromagnetic value in interactions matrix


	std::vector<spin> buffer;//holds the indices of spins to check
	IntMatrix cluster;//0 at indices that are not in cluster, 1 at indices that are
	std::vector<spin> not_cluster;//vector of spins that have not been added to the cluster
	std::vector<std::vector<int>> cluster_alt;//alternative cluster that may or may not be faster
	int cluster_size;//size of the cluster during a given MC step
	int cluster_mag;//total magnetization of the cluster
	int cutoff;//cutoff distance for checking individual spins in the long range algorithm
	void fill_not_cluster();
	
	
	void set_model();//set the model based on the model parameter in "params"
	void set_mean_field_model();//create interaction matrix for the mean field model
	void set_spin_boson_model();//create interaction matrix for the spin boson model
	//void set_cumulative_probs(int, int);//set cumulative probability vector for long range cluster forming


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
	
	double calc_sz_stagger();

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

	void print_interaction_sum() {
		std::cout << interaction_sum.to_string() << "\n";
	}

	void write_interactions() {
		std::ofstream file;
		file.open("interactions.csv");
		file << interactions.to_string() << "\n";
		file.close();
	}

	void write_resample_interactions();

	void output_state();

	void test_spins(spin);

	void set_alg_long_range_cluster() { LONG_RANGE_CLUSTER = true; NEAREST_NEIGHBOR_CLUSTER = false; }
	void set_alg_short_range_cluster() { NEAREST_NEIGHBOR_CLUSTER = true; LONG_RANGE_CLUSTER = false; }
	void set_alg_naive_cluster() { LONG_RANGE_CLUSTER = false; NEAREST_NEIGHBOR_CLUSTER = false; }
	void set_timers() { TIMERS = true;}
	void set_testing() { VERIFY_TESTING = true; }

};

class GeneralLRW{
	//This is a more lightweight version of the long range wolff class, designed to be
	//more general - hopefully to admit the use of user-defined "test_spins" functions
	//so that it can be used for many different MC simulations.  It does not hold the
	//lattice in order to keep the lattice easily available for GPU computations
	std::vector<std::vector<double>> interactions;//interaction matrix
	std::vector<std::vector<double>> interaction_sum;//values from interaction sum integral
	std::vector<double> site_sums;//absolute sum of interactions at each site
	double mag;//current magnetization
	double h; //applied field
	std::vector<spin> buffer;//holds the indices of spins to check
	std::vector<std::vector<int>> cluster;//0 at indices that are not in cluster, 1 at indices that are
	int cluster_size;//size of the cluster during a given MC step
	int cluster_mag;//total magnetization of the cluster
	void set_spin_boson_model_wc(class_mc_params);//version for finite cutoff frequency
	void set_spin_boson_model(class_mc_params);//create interaction matrix for the spin boson model
	void set_interaction_sum();
	void set_site_sums();
	MemTimeTester timer;

public:
	GeneralLRW(class_mc_params);

	void step(IsingLattice2D& lat);

    bool step_one_site(IsingLattice2D& lat, double *prev_action,
		thrust::host_vector<double>& thrustlat_ref, thrust::host_vector<double>& thrustint_ref, thrust::host_vector<double>& corr_ref);

	void test_spins(spin, IsingLattice2D& lat);

    void test_spins_one_site(spin, IsingLattice2D& lat);

	void set_mag(IsingLattice2D& lat){
		mag = 0;
		for (int i = 0; i < lat.get_Lx(); ++i){
			for (int j = 0; j < lat.get_Ly(); ++j){
				mag += lat.get_spin(i, j);
			}
		}
		mag = ((double) mag) / lat.get_N();
	}

	double get_mag(){
		return ((mag < 0) ? -mag : mag);
	}

	int get_cluster_size() { return cluster_size; }

    void get_thrust_interactions(thrust::host_vector<double>& result){
		if(result.size() != interactions.size()*interactions[0].size()){
			result.resize(interactions.size()*interactions[0].size());
		}
		for (int i = 0; i < interactions.size(); ++i){
			for (int j = 0; j < interactions[i].size(); ++j){
				result[i*interactions[i].size() + j] = interactions[i][j];
			}
		}
    }
    
	double calc_sx(IsingLattice2D&);

	double calc_sz_stagger(IsingLattice2D&);

	double calc_loc(IsingLattice2D&);

	double calc_s1s2(IsingLattice2D&);

	double calc_action_slow(IsingLattice2D&);

	void print_site_sums(){
		std::cout << "Site sums:\n" << vec2str(site_sums) << "\n";
	}

	void print_interactions() {
		std::stringstream outstring;
		outstring << "Interactions:\n";
		for (int i = 0; i < interactions.size(); ++i){
			outstring << vec2str(interactions[i]) << "\n";
		}
		std::cout << outstring.str() << "\n";
	}


	void print_interaction_sum() {
		std::stringstream outstring;
		outstring << "Interaction sum:\n";
		for (int i = 0; i < interaction_sum.size(); ++i){
			outstring << vec2str(interaction_sum[i]) << "\n";
		}
		std::cout << outstring.str() << "\n";
	}

	std::string get_int_string() {
		std::stringstream ss;
		ss << "Interactions:\n" ;
		for (int i = 0; i < interactions.size(); ++i){
			ss << vec2str(interactions[i]) << "\n";
		}
		return ss.str();
	}
/*
	void write_interactions() {
		std::ofstream file;
		file.open("interactions.csv");
		file << interactions.to_string() << "\n";
		file.close();
	}

	void write_resample_interactions();
*/
};
