#include <iostream>
#include <fstream>
#include <omp.h>
#include "MemTimeTester.h"
#include "IsingLattice2D.h"
#include "LongRangeWolff2D.h"
#include "class_mc_io.h"
#include "Matrix.h"
extern "C" {
#include "random.h"
}

LongRangeWolff2D* parallel_temp(LongRangeWolff2D*, LongRangeWolff2D*, LongRangeWolff2D*, LongRangeWolff2D*, LongRangeWolff2D*, LongRangeWolff2D*);

int main(int argc, char *argv[]) {
	IsingLattice2D lat = IsingLattice2D(1000, 1000, 4);
	MemTimeTester timer;

	/*
	Steps for Monte Carlo that need to be tested
	1. set up everything - lattice and parameters, set Wolff alg to begin (also prepare measurement struct
	2. begin simulation
		2.a. equilibrate - parallel tempering should be considered here
		2.b. take a measurement and record all relevent info
	3. compute any averages (if needed)
	4. output results
	*/

	//1. Test set up
	timer.flag_start_time("mc_setup");
	class_mc_params params;
	std::ifstream infile;
	std::cout << "Program name: " << argv[0] << "\n";
	std::cout << "Input file: " << argv[1] << "\n\n";
	infile.open(argv[1]);
	read_input_ising(&infile, &params);
	std::cout << params.to_string();

	IsingLattice2D test_lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
	LongRangeWolff2D wolff = LongRangeWolff2D(&test_lat, &params, params.model);
	class_mc_measurements results;
	results.names = { "steps", "times", "mags" , "energies"};
	results.values = { {}, {}, {} , {} };

	timer.flag_end_time("mc_setup");

	//2. Test Simulation
	//lattice and wolff are set up, just need to run and record measurements
	//include a parallel tempering step. use 4 threads for testing

	//set up parallel tempering lattices - goal is to have some higher-temperature lattices (lower beta)
	//stepping at the same time to allow for a better traversal of the state space
	double beta1 = 1.1*params.beta, beta2 = 1.05*params.beta, beta3 = .95*params.beta, beta4 = .9*params.beta;
	class_mc_params params1 = params, params2 = params, params3 = params, params4 = params;
	IsingLattice2D test_lat1 = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
	IsingLattice2D test_lat2 = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
	IsingLattice2D test_lat3 = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
	IsingLattice2D test_lat4 = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
	LongRangeWolff2D wolff1 = LongRangeWolff2D(&test_lat1, &params1, params.model);
	LongRangeWolff2D wolff2 = LongRangeWolff2D(&test_lat2, &params2, params.model);
	LongRangeWolff2D wolff3 = LongRangeWolff2D(&test_lat3, &params3, params.model);
	LongRangeWolff2D wolff4 = LongRangeWolff2D(&test_lat4, &params4, params.model);
	wolff1.set_beta(beta1);
	wolff2.set_beta(beta2);
	wolff3.set_beta(beta3);
	wolff4.set_beta(beta4);
	LongRangeWolff2D* primary_wolff_p = &wolff;

	//start stepping
	for (int eq = 0; eq < params.eq_time; ++eq) {
		wolff.step();
		//wolff1.step();
		//wolff2.step();
		//wolff3.step();
		//wolff4.step();
	}

	timer.flag_start_time("mc_run");
	for (int dump = 0; dump < params.max_dumps; ++dump) {
		for (int measure = 0; measure < params.measures_per_dump; ++measure) {
			for (int steps = 0; steps < params.steps_per_measure; ++steps) {
					wolff.step();
					//wolff1.step();
					//wolff2.step();
					//wolff3.step();
					//wolff4.step();
			}
			results.record("steps", (double)(dump*params.measures_per_dump*params.steps_per_measure + measure*params.steps_per_measure + params.eq_time));
			results.record("times", timer.get_running_time("mc_run"));
			results.record("mags", primary_wolff_p->calc_mag());
			results.record("energies", primary_wolff_p->calc_E_mean_field_model());
			//primary_wolff_p = parallel_temp(&wolff, &wolff1, &wolff2, &wolff3, &wolff4, primary_wolff_p);
		}
		std::vector<int> real_steps;
		std::vector<double> bad_steps = results.get_vals("steps");
		for (int i = 0; i < bad_steps.size(); ++i) {
			real_steps.push_back((int)bad_steps[i]);
		}
		write_outputs(dump, real_steps, results.get_vals("times"), results.get_vals("energies"), results.get_vals("mags"));
	}
	timer.flag_end_time("mc_run");

	//3. Compute averages
	//mag, mag^2, mag^4 (for binder cumulant)

	std::vector<double> mags = results.get_vals("mags");
	double mag2 = 0;
	double mag2_avg = 0;
	double mag4_avg = 0;
	for (int i = 0; i < mags.size(); ++i) {
		mag2 = mags[i] * mags[i];
		mag2_avg += mag2;
		mag4_avg += mag2*mag2;
	}
	mag2_avg = mag2_avg / mags.size();
	mag4_avg = mag4_avg / mags.size();
	double Q = mag2_avg * mag2_avg / mag4_avg;

	std::vector<double> energies = results.get_vals("energies");
	double e = 0;
	double e2 = 0;
	for (int i = 0; i < energies.size(); ++i) {
		e += energies[i];
		e2 += energies[i] * energies[i];
	}
	e = e / energies.size();
	e2 = e2 / energies.size();
	std::cout << "Spec heat: " << params.beta * params.beta / test_lat.get_N() * (e2 - e*e) << "\n";

	//4. output results

	std::cout << "Q = " << Q << "\n";

	timer.print_timers();

	return 0;
}

LongRangeWolff2D* parallel_temp(LongRangeWolff2D* set1, LongRangeWolff2D* set2, LongRangeWolff2D* set3, LongRangeWolff2D* set4, LongRangeWolff2D* set5, LongRangeWolff2D* prim_p) {
	double swap;

	double delta_beta_E = (set2->get_beta() - set1->get_beta())*(set2->calc_E_mean_field_model() - set1->calc_E_mean_field_model());
	//std::cout << "Temp swap: " << set2->get_beta() << ", " << set1->get_beta() << ", Probability: " << (exp(delta_beta_E) < 1 ? exp(delta_beta_E) : 1) << "\n";
	if (drand1_() < exp(delta_beta_E)) {
		swap = set1->get_beta();
		if (swap == prim_p->get_beta()) { prim_p = set2; }
		set1->set_beta(set2->get_beta());
		set2->set_beta(swap);
	}

	delta_beta_E = (set3->get_beta() - set2->get_beta())*(set3->calc_E_mean_field_model() - set2->calc_E_mean_field_model());
	//std::cout << "Temp swap: " << set3->get_beta() << ", " << set2->get_beta() << ", Probability: " << (exp(delta_beta_E) < 1 ? exp(delta_beta_E) : 1) << "\n";
	if (drand1_() < exp(delta_beta_E)) {
		swap = set2->get_beta();
		if (swap == prim_p->get_beta()) { prim_p = set3; }
		set2->set_beta(set3->get_beta());
		set3->set_beta(swap);
	}


	delta_beta_E = (set4->get_beta() - set3->get_beta())*(set4->calc_E_mean_field_model() - set3->calc_E_mean_field_model());
	//std::cout << "Temp swap: " << set4->get_beta() << ", " << set3->get_beta() << ", Probability: " << (exp(delta_beta_E) < 1 ? exp(delta_beta_E) : 1) << "\n";
	if (drand1_() < exp(delta_beta_E)) {
		swap = set3->get_beta();
		if (swap == prim_p->get_beta()) { prim_p = set4; }
		set3->set_beta(set4->get_beta());
		set4->set_beta(swap);
	}

	delta_beta_E = (set5->get_beta() - set4->get_beta())*(set5->calc_E_mean_field_model() - set4->calc_E_mean_field_model());
	//std::cout << "Temp swap: " << set5->get_beta() << ", " << set4->get_beta() << ", Probability: " << (exp(delta_beta_E) < 1 ? exp(delta_beta_E) : 1) << "\n";
	if (drand1_() < exp(delta_beta_E)) {
		swap = set4->get_beta();
		if (swap == prim_p->get_beta()) { prim_p = set5; }
		set4->set_beta(set5->get_beta());
		set5->set_beta(swap);
	}

	return prim_p;

}