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

void process_timer_results(timer_list, std::vector<class_mc_params>);
void process_autocorrelation_results(std::vector<std::vector<double>>, std::vector<class_mc_params>);
void parallel_tempering(std::vector<LongRangeWolff2D>::iterator, std::vector<LongRangeWolff2D>::iterator);

int main(int argc, char *argv[]) {

	/*
	Benchmarking the Spin-Boson Monte Carlo code
	Hard-code inputs, can maybe have input file later
	0. Set up an array of different wolff engines with pre-defined parameters
	1. Set timers and track them for different parts of the algorithm, then output the results in
		an easy-to-read format.  Should be used to find the average step time for different
		parameter sets.
	2. Calculate the autocorrelation function for a given set of parameters - this can be used
		to find the best choice for number of steps between mc measurements
	3. Calculate the average acceptance ratio for a set of parallel tempering temperatures
		- this will probably depend on other parameters of the hamiltonian as well so this
			should be possible to check for a variety of different input situations
	4. Parallelization: want to do a test run of some small lattices including omp and mpi to ensure
		that these can be used for large simulations
	5. Data analysis: want to be able to do some basic scaling and fitting calculations within the c++
		code rather than just outputting every mc step.
	6. Test the long-range-wolff for the spin-boson model.  Should be able to do some outputs and check
		by hand that the algorithm is working.
	*/

	/**
	--------------------------------------
	0. Set up range of wolff engines
	--------------------------------------
	**/
	const int NUM_TESTS = 100;
	std::vector<int> lattice_lengths = { 8, 16, 32, 128 };//go from small lattice to large lattice.  Runs will be on square lattices.
	std::vector<double> betas = { 0.1, 0.5, 1.0, 5.0, 10.0 }, deltas = { 0.1, 0.5, 1.0, 5.0, 10.0 };//set J = 1.  Then test a range of the other params.
	MemTimeTester timer;
	std::vector<class_mc_params> params_list;
	std::cout << "Benchmarking Quantum Ising for param list: \n" << "Lattice Sizes: " << vec2str(lattice_lengths) << std::endl
		<< "Betas: " << vec2str(betas) << std::endl << "Deltas: " << vec2str(deltas) << std::endl;
	std::cout << "Setting J = 1.0, spacial lattice spacing = 1.0" << std::endl << std::endl << std::endl;
	//loop through parameters and make a lattice, class_mc_params, and wolff for each
	std::vector<IsingLattice2D> lattices;
	std::vector<LongRangeWolff2D> wolffs;
	for (int lengths : lattice_lengths) {
		for (double beta : betas) {
			for (double delta : deltas) {
				class_mc_params params;
				params.set_q_ising_dummy(lengths, beta, delta);
				apply_spin_boson_params(&params);
				//std::cout << params.to_string() << std::endl;
				params_list.push_back(params);
				IsingLattice2D lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
				LongRangeWolff2D wolff = LongRangeWolff2D(&lat, &params);
				wolff.set_alg_short_range_cluster();
				wolff.set_timers();
				lattices.push_back(lat);
				wolffs.push_back(wolff);
			}
		}
	}

	if (argc > 1) {

		/**
		---------------------------------------------------------------------------------------------------------
		1.  Run some equilibrations and use timers to find average step times and normalized average cluster size
		---------------------------------------------------------------------------------------------------------
		**/
		//also calculate variances to see if some might need more extensive testing
		//(for example, some might start fast and get longer as they equilibrate or vice versa
	//#pragma omp parallel for num_threads(4)
		for (int i = 0; i < wolffs.size(); ++i) {
			std::string timer_name = "wolff_equil_timer" + std::to_string(i);
			std::cout << "Running wolff equilibration timer #" << i << std::endl;
			timer.flag_start_time(timer_name);
			for (int j = 0; j < params_list[i].eq_time; ++j) {
				//for(int j = 0; j < 10; ++j){
				wolffs[i].step();
			}
			timer.flag_end_time(timer_name);
		}
		process_timer_results(timer.get_timers(), params_list);


		/**
		--------------------------------------------------------------------------------------------
		2.  Autocorrelation functions - find how long each wolff takes to change to a "random" state
		--------------------------------------------------------------------------------------------
		**/
		std::vector<std::vector<double>> autocorr;
		//#pragma omp parallel for num_threads(4)
		for (int i = 0; i < wolffs.size(); ++i) {
			std::cout << "Running wolff autocorrelation calc #" << i << std::endl;
			autocorr.push_back({});
			for (int j = 0; j < params_list[i].eq_time; ++j) {
				for (int k = 0; k < 1; ++k) {
					wolffs[i].step();
				}
				autocorr[i].push_back(abs(wolffs[i].calc_mag()));
			}
		}
		process_autocorrelation_results(autocorr, params_list);
	}

	/**
	----------------------------------------------------
	3.  Parallel Tempering and average acceptance ratios
	----------------------------------------------------
	**/
	//to perform parallel tempering on a q-c mapped lattice, we must choose
	//lattices with different spacings in the time direction to maintain the same size
	//for easy switching during the parallel tempering step
	std::cout << std::endl << "Testing parallel tempering and acceptance ratios" << std::endl << std::endl;

	//do around J = Delta = Beta = 1 for 32x32 lattice
	betas = { 0.9, 0.95, 1.0, 1.05, 1.1 };
	params_list = {};
	lattices = {};
	wolffs = {};
	for (double beta : betas) {
		class_mc_params params;
		params.set_q_ising_dummy(32, beta, 1.0);
		apply_spin_boson_params(&params);
		//std::cout << params.to_string() << std::endl;
		params_list.push_back(params);
		IsingLattice2D lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
		LongRangeWolff2D wolff = LongRangeWolff2D(&lat, &params);
		wolff.set_alg_short_range_cluster();
		wolff.set_timers();
		lattices.push_back(lat);
		wolffs.push_back(wolff);
	}
	//attempt to do this in parallel
	for (int i = 0; i < 100; ++i) {
		for (int j = 0; j < 10; ++j) {
#pragma omp parallel for num_threads(4)
			for (int k = 0; k < wolffs.size(); ++k) {
				wolffs[k].step();
			}
		}
		parallel_tempering(wolffs.begin(), wolffs.end());
	}




	//class_mc_measurements results;
	//results.names = { "steps", "times", "sx" , "sz"};
	//results.values = { {}, {}, {} , {} };


	//results.func_names = { "corr_t" };
	//results.functions = { {} };
	//results.function_num_measures = { 0 };


	//timer.flag_end_time("mc_setup");

	////2. Test Simulation
	////lattice and wolff are set up, just need to run and record measurements

	////set up parallel tempering lattices - goal is to have some higher-temperature lattices (lower beta)
	////stepping at the same time to allow for a better traversal of the state space
	////double beta1 = 1.1*params.beta, beta2 = 1.05*params.beta, beta3 = .95*params.beta, beta4 = .9*params.beta;
	////class_mc_params params1 = params, params2 = params, params3 = params, params4 = params;
	////IsingLattice2D test_lat1 = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
	////IsingLattice2D test_lat2 = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
	////IsingLattice2D test_lat3 = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
	////IsingLattice2D test_lat4 = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
	////LongRangeWolff2D wolff1 = LongRangeWolff2D(&test_lat1, &params1, params.model);
	////LongRangeWolff2D wolff2 = LongRangeWolff2D(&test_lat2, &params2, params.model);
	////LongRangeWolff2D wolff3 = LongRangeWolff2D(&test_lat3, &params3, params.model);
	////LongRangeWolff2D wolff4 = LongRangeWolff2D(&test_lat4, &params4, params.model);
	////wolff1.set_beta(beta1);
	////wolff2.set_beta(beta2);
	////wolff3.set_beta(beta3);
	////wolff4.set_beta(beta4);
	//LongRangeWolff2D* primary_wolff_p = &wolff;

	////start stepping
	//for (int eq = 0; eq < params.eq_time; ++eq) {
	//	wolff.step();
	//	//wolff1.step();
	//	//wolff2.step();
	//	//wolff3.step();
	//	//wolff4.step();
	//}

	//timer.flag_start_time("mc_run");
	//for (int dump = 0; dump < params.max_dumps; ++dump) {
	//	for (int measure = 0; measure < params.measures_per_dump; ++measure) {
	//		for (int steps = 0; steps < params.steps_per_measure; ++steps) {
	//				wolff.step();
	//				//wolff1.step();
	//				//wolff2.step();
	//				//wolff3.step();
	//				//wolff4.step();
	//		}
	//		results.record("steps", (double)(dump*params.measures_per_dump*params.steps_per_measure + measure*params.steps_per_measure + params.eq_time));
	//		results.record("times", timer.get_running_time("mc_run"));
	//		results.record("sx", primary_wolff_p->calc_sx());
	//		results.record("sz", primary_wolff_p->calc_sz());
	//		results.record("corr_t", primary_wolff_p->calc_corr(1));
	//		//primary_wolff_p = parallel_temp(&wolff, &wolff1, &wolff2, &wolff3, &wolff4, primary_wolff_p);
	//		//write_state(dump*params.measures_per_dump + measure, lat, 5.0);
	//	}
	//	std::vector<int> real_steps;
	//	std::vector<double> bad_steps = results.get_vals("steps");
	//	for (int i = 0; i < bad_steps.size(); ++i) {
	//		real_steps.push_back((int)bad_steps[i]);
	//	}
	//	write_outputs(dump, real_steps, results.get_vals("times"), results.get_vals("sz"), results.get_vals("sx"), results.get_func("corr_t"));
	//	std::cout << "Dump " << dump + 1 << " out of " << params.max_dumps << "\n";
	//}
	//timer.flag_end_time("mc_run");

	////3. Compute averages
	////mag, mag^2, mag^4 (for binder cumulant)

	//std::vector<double> mags = results.get_vals("sz");
	//double mag2 = 0;
	//double mag2_avg = 0;
	//double mag4_avg = 0;
	//for (int i = 0; i < mags.size(); ++i) {
	//	mag2 = mags[i] * mags[i];
	//	mag2_avg += mag2;
	//	mag4_avg += mag2*mag2;
	//}
	//mag2_avg = mag2_avg / mags.size();
	//mag4_avg = mag4_avg / mags.size();
	//double Q = mag2_avg * mag2_avg / mag4_avg;

	////std::vector<double> energies = results.get_vals("sz");
	////double e = 0;
	////double e2 = 0;
	////for (int i = 0; i < energies.size(); ++i) {
	////	e += energies[i];
	////	e2 += energies[i] * energies[i];
	////}
	////e = e / energies.size();
	////e2 = e2 / energies.size();
	////std::cout << "Spec heat: " << params.beta * params.beta / test_lat.get_N() * (e2 - e*e) << "\n";

	////4. output results

	//std::cout << "Q = " << Q << "\n";

	//timer.print_timers();

	return 0;
}



void process_timer_results(timer_list timers, std::vector<class_mc_params> params_list) {
	//calculate step times as a function of different simulation parameters
	//first, output a file that just lists params and times
	std::ofstream file;
	file.open("bare_timer_results.txt");
	for (int i = 0; i < timers.size; ++i) {
		file << "Size: " << params_list[i].lengths[0] << "x" << params_list[i].lengths[0] << std::endl;
		file << "Beta: " << params_list[i].beta << std::endl;
		file << "Delta: " << params_list[i].sbparams.delta << std::endl;
		file << "Eq time: " << timers.times[i] << std::endl << std::endl << std::endl;
	}

}

void process_autocorrelation_results(std::vector<std::vector<double>> mags, std::vector<class_mc_params> params_list) {
	//take a list of autocorrelation functions, match them with their input params, and output to csv for easy plotting
	std::ofstream file;
	file.open("autocorrelation_results.csv");
	std::vector<std::vector<double>> func;
	//calculate autocorrelation function
	for (int i = 0; i < mags.size(); ++i) {
		int tmax = mags[i].size();
		func.push_back({});
		for (int t = 1; t < tmax; ++t) {
			double m0 = 0, mt = 0, both = 0;
			for (int tprime = 0; tprime < tmax - t; ++tprime) {
				both += mags[i][tprime] * mags[i][tprime + t];
				m0 += mags[i][tprime];
				mt += mags[i][tprime + t];
			}
			func[i].push_back(0.1 / (tmax - t) * both - 0.01 / (tmax - t) / (tmax - t) * m0 * mt);
		}
		//normalize f
		double factor = func[i][0];
		for (int j = 0; j < func[i].size(); ++j) {
			func[i][j] = func[i][j] / factor;
		}
	}
	for (int i = 0; i < params_list.size(); ++i) {
		file << "Size," << params_list[i].lengths[0] << "," << params_list[i].lengths[0] << std::endl;
		file << "Beta," << params_list[i].beta << std::endl;
		file << "Delta," << params_list[i].sbparams.delta << std::endl;
		file << "autocorr," << vec2str(func[i]) << std::endl << std::endl << std::endl;
	}
}

void parallel_tempering(std::vector<LongRangeWolff2D>::iterator  wolffs_in_begin, std::vector<LongRangeWolff2D>::iterator wolffs_in_end) {

}