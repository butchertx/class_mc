#include <iostream>
#include <cmath>
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

void process_timer_cluster_results(timer_list, std::vector<double>, std::vector<class_mc_params>);
void process_autocorrelation_results(std::vector<std::vector<double>>, std::vector<class_mc_params>);
std::vector<double> parallel_tempering(std::vector<LongRangeWolff2D>::iterator, std::vector<LongRangeWolff2D>::iterator);



int main(int argc, char *argv[]) {

	/*
	Benchmarking the Spin-Boson Monte Carlo code
	Hard-code inputs, can maybe have input file later
	0. Set up small wolff sim and use it to check individual steps
	1. Set timers and track them for different parts of the algorithm, then output the results in
		an easy-to-read format.  Should be used to find the average step time for different
		parameter sets. Also find average cluster sizes.
	2. Calculate the autocorrelation function for a given set of parameters - this can be used
		to find the best choice for number of steps between mc measurements
	3. Calculate the average acceptance ratio for a set of parallel tempering temperatures
		- this will probably depend on other parameters of the hamiltonian as well so this
			should be possible to check for a variety of different input situations
	4. Parallelization: want to do a test run of some small lattices including omp and mpi to ensure
		that these can be used for large simulations
	5. Data analysis: want to be able to do some basic scaling and fitting calculations within the c++
		code rather than just outputting every mc step.
	6. Create images for gif creation of algorithm/equilibration
	7. Create test data for the Swendson resampling procedure
	*/


	std::string mystr;

	/**
	---------------------------------------------------------------
	-1. Test fast calc of observables
	---------------------------------------------------------------
	**/
	std::cout << "Test observable calculation? Y or N\n";
	getline(std::cin, mystr);
	if (mystr.compare("Y") == 0 || mystr.compare("y") == 0) {
		class_mc_params params;
		std::ifstream infile;
		std::cout << "Program name: " << argv[0] << "\n";
		std::cout << "Input file: " << argv[1] << "\n\n";
		infile.open(argv[1]);
		read_input_ising(&infile, &params);
		read_input_spin_boson(&infile, &(params.sbparams));
		IsingLattice2D lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
		LongRangeWolff2D wolff = LongRangeWolff2D(&lat, &params);
		MemTimeTester timer;
		apply_spin_boson_params(&params);
		std::cout << "Base parameters: \n" << params.to_string();
		//wolff.print_interactions();
		//wolff.set_testing();
		wolff.set_alg_long_range_cluster();
		class_mc_measurements results;
		results.names = { "mag", "mag2"};
		results.values = { {},{} };
		timer.flag_start_time("steps");
		for (int i = 0; i < params.measures_per_dump; ++i) {
			for (int j = 0; j < params.steps_per_measure; ++j) {
				wolff.step();
			}
			results.record("mag", wolff.get_mag());
			results.record("mag2", wolff.get_mag()*wolff.get_mag());
		}
		timer.flag_end_time("steps");
		//write the state



		timer.print_timers();
		write_final_outputs(results, 1);
		wolff.write_interactions();
		double sus_err = bootstrap(results.get_vals("mag"), 500, "susceptibility");
		std::cout << "bootstrap susceptibility error: " << sus_err << "\n";
		std::cout << "naive susceptiblity: " << mean(results.get_vals("mag2")) - mean(results.get_vals("mag"))*mean(results.get_vals("mag")) << " \n";
	}

	/**
	---------------------------------------------------------------
	0. Set up small simulation and use it to check individual steps
	---------------------------------------------------------------
	**/

	std::cout << "Test algorithm? Y or N\n";
	getline(std::cin, mystr);
	if (mystr.compare("Y") == 0 || mystr.compare("y") == 0) {
		class_mc_params params;
		std::ifstream infile;
		std::cout << "Program name: " << argv[0] << "\n";
		std::cout << "Input file: " << argv[1] << "\n\n";
		infile.open(argv[1]);
		read_input_ising(&infile, &params);
		read_input_spin_boson(&infile, &(params.sbparams));
		IsingLattice2D lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
		LongRangeWolff2D wolff = LongRangeWolff2D(&lat, &params);
		MemTimeTester timer;
		apply_spin_boson_params(&params);
		std::cout << "Base parameters: \n" << params.to_string();
		//wolff.set_testing();
		//wolff.print_interaction_sum();


/*
		wolff.set_alg_naive_cluster();
		timer.flag_start_time("naive cluster");
		for (int i = 0; i < params.steps_per_measure; ++i) {
			wolff.step();
			timer.flag_start_time("naive measurements");
			wolff.calc_mag();
			wolff.calc_corr(1);
			wolff.calc_sz();
			wolff.calc_sx();
			timer.flag_end_time("naive measurements");
		}
		timer.flag_end_time("naive cluster");
*/
		wolff.set_alg_long_range_cluster();
		timer.flag_start_time("long range cluster");
		for (int i = 0; i < params.steps_per_measure; ++i) {
			wolff.step();
		}
		timer.flag_end_time("long range cluster");
		timer.print_timers();
	}
	/**
	----------------------------------------------------------
	1, 2. Set up range of wolff engines and run step timers
	----------------------------------------------------------
	**/

	std::cout << "Test step time? Y or N\n";
	getline(std::cin, mystr);
	if (mystr.compare("Y") == 0 || mystr.compare("y") == 0) {
		//get base parameters
		class_mc_params base_params;
		std::ifstream infile;
		std::cout << "Program name: " << argv[0] << "\n";
		std::cout << "Input file: " << argv[1] << "\n\n";
		infile.open(argv[1]);
		read_input_ising(&infile, &base_params);
		read_input_spin_boson(&infile, &(base_params.sbparams));
		std::cout << "Base parameters: \n" << base_params.to_string();

		//set up list of different params
		const int NUM_TESTS = 1;
		std::vector<double> alphas = {0.5};
		MemTimeTester timer;
		std::vector<class_mc_params> params_list;
		std::cout << "Benchmarking Spin Boson for param list: \n" << "Alpha = " << vec2str(alphas) << "\n";
		//loop through parameters and make a lattice, class_mc_params, and wolff for each
		std::vector<IsingLattice2D> lattices;
		std::vector<LongRangeWolff2D> wolffs;
		std::vector<class_mc_measurements> results_list;
		for (double alpha : alphas) {
			class_mc_params params = base_params;
			params.sbparams.A0 = alpha;
			apply_spin_boson_params(&params);
			params_list.push_back(params);
			IsingLattice2D lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
			LongRangeWolff2D wolff = LongRangeWolff2D(&lat, &params);
			class_mc_measurements results;
			results.names = { "mag", "mag2", "mag4" };
			results.values = { {}, {}, {} };
			wolff.set_timers();
			wolff.set_alg_long_range_cluster();
			//wolff.set_alg_naive_cluster();
			lattices.push_back(lat);
			wolffs.push_back(wolff);
			results_list.push_back(results);
		}

		/**
		---------------------------------------------------------------------------------------------------------
		1.  Run some equilibrations and use timers to find average step times and normalized average cluster size
		---------------------------------------------------------------------------------------------------------
		**/
		std::cout << "Test eq timers? Y or N\n";
		getline(std::cin, mystr);
		if (mystr.compare("Y") == 0 || mystr.compare("y") == 0) {
			//also calculate variances to see if some might need more extensive testing
			//(for example, some might start fast and get longer as they equilibrate or vice versa
			std::cout << "Running Equilibration Timers...\n\n";
			std::vector<double> cluster_sizes;
			//std::vector<int> seeds = { 1, 2, 3, 4, 5, 6, 7, 8 };
	//#pragma omp parallel for shared(seeds, wolffs, params_list, results_list, timer) num_threads(4)
			for (int i = 0; i < wolffs.size(); ++i) {
				double mag;
				//rand_init_(&seeds[i]);
				double cluster_size = 0;
				std::string timer_name = "wolff_equil_timer" + std::to_string(i);
				std::cout << "Running wolff equilibration timer #" << i << std::endl;
				std::cout << "Parameters for A0 = " << params_list[i].sbparams.A0 << ": \n" << params_list[i].to_string() << "\n\n";
				timer.flag_start_time(timer_name);
				for (int j = 0; j < params_list[i].eq_time; ++j) {
					for (int k = 0; k < params_list[i].measures_per_dump; ++k) {
						wolffs[i].step();
					}
					cluster_size += wolffs[i].get_cluster_size();
					mag = wolffs[i].calc_mag();
					results_list[i].record("mag", mag);
					results_list[i].record("mag2", mag*mag);
					results_list[i].record("mag4", mag*mag*mag*mag);
				}
				timer.flag_end_time(timer_name);
				cluster_sizes.push_back(cluster_size / params_list[i].eq_time / wolffs[i].get_lat().get_N());
			}
			process_timer_cluster_results(timer.get_timers(), cluster_sizes, params_list);
			//output alphas, then binder cumulants
			std::ofstream file;
			file.open("results.csv");
			for (double alpha : alphas) {
				file << alpha << ",";
			}
			file << "\n";
			for (class_mc_measurements results : results_list) {
				file << mean(results.get_vals("mag2")) << ",";
			}
			file.close();
		}


		/**
		--------------------------------------------------------------------------------------------
		2.  Autocorrelation functions - find how long each wolff takes to change to a "random" state
		--------------------------------------------------------------------------------------------
		**/
		std::cout << "Test autocorrelation? Y or N\n";
		getline(std::cin, mystr);
		if (mystr.compare("Y") == 0 || mystr.compare("y") == 0){
			std::cout << "Running Autocorrelation Tests for m...\n\n";
			std::vector<std::vector<double>> autocorr;
			double avg_cluster = 0;
			//		#pragma omp parallel for num_threads(4)
			for (int i = 0; i < wolffs.size(); ++i) {
				std::cout << "Running wolff autocorrelation calc #" << i << std::endl;
				autocorr.push_back({});
				for (int j = 0; j < params_list[i].eq_time; ++j) {
					for (int k = 0; k < 1; ++k) {
						wolffs[i].step();
					}
					autocorr[i].push_back(wolffs[i].calc_sz());
					avg_cluster += wolffs[i].get_cluster_size();
				}
				std::cout << "Average cluster fraction for wolff #" << i << ": " << avg_cluster / params_list[i].eq_time / wolffs[i].get_lat().get_N() << "\n";
			}
			process_autocorrelation_results(autocorr, params_list);
		}
	}


	/**
	----------------------------------------------------
	3.  Parallel Tempering and average acceptance ratios
	----------------------------------------------------
	**/
	//to perform parallel tempering on a q-c mapped lattice, we must choose
	//lattices with different spacings in the time direction to maintain the same size
	//for easy switching during the parallel tempering step

	std::cout << "Test parallel tempering acceptance ratios? Y or N\n";
	getline(std::cin, mystr);
	if (mystr.compare("Y") == 0 || mystr.compare("y") == 0) {
		std::cout << std::endl << "Testing parallel tempering and acceptance ratios" << std::endl << std::endl;

		//do around J = Delta = Beta = 1 for 32x32 lattice
		std::vector<double> betas = { 0.9, 0.95, 1.0, 1.05, 1.1 };
		std::vector<class_mc_params> params_list = {};
		std::vector<IsingLattice2D> lattices = {};
		std::vector<LongRangeWolff2D> wolffs = {};
		std::cout << "Successfully allocated parallel tempering lattices\n";
		for (double beta : betas) {
			class_mc_params params;
			params.set_q_ising_dummy(32, beta, 1.0);
			apply_spin_boson_params(&params);
			//std::cout << params.to_string() << std::endl;
			params_list.push_back(params);
			IsingLattice2D lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
			LongRangeWolff2D wolff = LongRangeWolff2D(&lat, &params);
			//wolff.set_alg_short_range_cluster();
			wolff.set_timers();
			lattices.push_back(lat);
			wolffs.push_back(wolff);
		}
		//attempt to do this in parallel
		std::vector<double> acc_ratios(wolffs.size() - 1);
		std::vector<double> avg_acc_ratios(wolffs.size() - 1, 0.0);
		for (int i = 0; i < 100; ++i) {
			for (int j = 0; j < 10; ++j) {
				for (int k = 0; k < wolffs.size(); ++k) {
					wolffs[k].step();
				}
			}
			acc_ratios = parallel_tempering(wolffs.begin(), wolffs.end());
			for (int m = 0; m < acc_ratios.size(); ++m) {
				avg_acc_ratios[m] += acc_ratios[m];
			}
		}
		std::cout << "Average acceptance ratios: \n";
		for (int i = 0; i < avg_acc_ratios.size(); ++i) {
			std::cout << avg_acc_ratios[i] / 100.0 << "\n";
		}
	}

	/*
	-----------------------------------------------------
	4. Test Parallelization
	-----------------------------------------------------
	*/
	//set up some parallel lattices and see if they can be iterated in a parallelized manner
	//specifically, check time and ability to choose independent random numbers
	std::cout << "Test parallelization? Y or N\n";
	getline(std::cin, mystr);
	if (mystr.compare("Y") == 0 || mystr.compare("y") == 0) {

		std::cout << "First, check parallel random numbers.  Start a loop with 4 threads and output random numbers from each" << "\n";
		std::vector<double> rands;
		std::vector<int> seeds = { 1, 2, 3, 4 };
		MemTimeTester timer;

		timer.flag_start_time("parallel_num_test");
#pragma omp parallel for private(rands) shared(seeds) num_threads(4)
		for (int i = 0; i < 4; ++i) {
			rand_init_(&seeds[i]);
			for (int j = 0; j < 30; ++j) {
				rands.push_back(drand1_());
			}
			std::cout << vec2str(rands) << "\n";
		}
		timer.flag_end_time("parallel_num_test");


		std::cout << "\nNow check serial random numbers\n";
		timer.flag_start_time("serial_num_test");
		rands = {};
		rand_init_(&seeds[0]);
		for (int i = 0; i < 120; ++i) {
			rands.push_back(drand1_());
		}
		std::cout << vec2str(rands) << "\n";
		timer.flag_end_time("serial_num_test");
		timer.print_timers();


	}

	/*
	-----------------------------------------------------
	6.  Make a gif of state equilibration
	-----------------------------------------------------
	*/
	std::cout << "Output images of states? Y or N\n";
	getline(std::cin, mystr);
	if (mystr.compare("Y") == 0 || mystr.compare("y") == 0){
		class_mc_params base_params;
		std::ifstream infile;
		std::cout << "Program name: " << argv[0] << "\n";
		std::cout << "Input file: " << argv[1] << "\n\n";
		infile.open(argv[1]);
		read_input_ising(&infile, &base_params);
		read_input_spin_boson(&infile, &(base_params.sbparams));
		std::cout << "Base parameters: \n" << base_params.to_string();

		IsingLattice2D lat = IsingLattice2D(base_params.lengths[0], base_params.lengths[1], base_params.rand_seed);
		lat.set_ferro();
		LongRangeWolff2D wolff = LongRangeWolff2D(&lat, &base_params);
		MemTimeTester timer;
		apply_spin_boson_params(&base_params);

		IsingLattice2D lat_copy = wolff.get_lat();
		IsingLattice2D& lat_ref = lat_copy;
		write_state_pbm(0, lat_ref);
		//wolff.set_alg_short_range_cluster();
		wolff.set_alg_short_range_cluster();

		std::cout << "Equilibrating and outputting each step...\n";
		for (int i = 1; i < base_params.eq_time; ++i){
			wolff.step();
			lat_copy = wolff.get_lat();
			write_state_pbm(i, lat_ref);
		}

		

	}

	/*
	-----------------------------------------------------
	7.  Make test data for Swendson resampling procedure
	-----------------------------------------------------
	*/
	std::cout << "Output states for resampling? Y or N\n";
	getline(std::cin, mystr);
	if (mystr.compare("Y") == 0 || mystr.compare("y") == 0) {
		//get base parameters
		class_mc_params base_params;
		std::ifstream infile;
		std::cout << "Program name: " << argv[0] << "\n";
		std::cout << "Input file: " << argv[1] << "\n\n";
		infile.open(argv[1]);
		read_input_ising(&infile, &base_params);
		read_input_spin_boson(&infile, &(base_params.sbparams));
		std::cout << "Base parameters: \n" << base_params.to_string();

		//set up list of different params.  Will calculate properties for these but only output states for one of them
		const int NUM_TESTS = 11;
		std::vector<double> alphas = {1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2};
		MemTimeTester timer;
		std::vector<class_mc_params> params_list;
		std::cout << "Generating resampling test data for param list: \n" << "Alpha = " << vec2str(alphas) << "\n";
		//loop through parameters and make a lattice, class_mc_params, and wolff for each
		std::vector<IsingLattice2D> lattices;
		std::vector<LongRangeWolff2D> wolffs;
		std::vector<class_mc_measurements> results_list;
		for (double alpha : alphas) {
			class_mc_params params = base_params;
			params.sbparams.A0 = alpha;
			apply_spin_boson_params(&params);
			params_list.push_back(params);
			IsingLattice2D lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
			LongRangeWolff2D wolff = LongRangeWolff2D(&lat, &params);
			class_mc_measurements results;
			results.names = { "mag", "mag2", "mag4" };
			results.values = { {},{},{} };
			wolff.set_timers();
			wolff.set_alg_long_range_cluster();
			lattices.push_back(lat);
			wolffs.push_back(wolff);
			results_list.push_back(results);
		}
		IsingLattice2D latcopy = wolffs[5].get_lat();
		wolffs[5].write_resample_interactions();
		IsingLattice2D& latref = latcopy;
		double mag_temp = 0;
		for (int i = 0; i < wolffs.size(); ++i) {
			std::cout << "Calculating wolff #" << i + 1 << " out of " << wolffs.size() << "...\n";
			for (int eq = 0; eq < params_list[i].eq_time; ++eq) {
				wolffs[i].step();
			}
			for (int meas = 0; meas < params_list[i].measures_per_dump; ++meas) {
				for (int step = 0; step < params_list[i].steps_per_measure; ++step) {
					wolffs[i].step();
				}
				if (params_list[i].sbparams.A0 == 1.1) {
					latcopy = wolffs[i].get_lat();
					write_state(meas, latref);
				}
				mag_temp = wolffs[i].get_mag();
				results_list[i].record("mag", mag_temp);
				results_list[i].record("mag2", mag_temp*mag_temp);
				results_list[i].record("mag4", mag_temp*mag_temp*mag_temp*mag_temp);
			}
		}

		//write results in an easy-to-plot manner
		std::ofstream file;
		file.open("results.csv");
		for (double alpha : alphas) {
			file << alpha << ",";
		}
		file << "\n";
		for (class_mc_measurements results : results_list) {
			file << mean(results.get_vals("mag2")) << ",";
		}
		file.close();
	}

	return 0;
}



void process_timer_cluster_results(timer_list timers, std::vector<double> avg_clusters, std::vector<class_mc_params> params_list) {
	//calculate step times as a function of different simulation parameters
	//first, output a file that just lists params and times
	std::ofstream file;
	file.open("bare_timer_results.txt");
	for (int i = 0; i < timers.size; ++i) {
		file << "J: " << params_list[i].Js[0] << "\n";
		file << "Gamma: " << params_list[i].Js[1] << "\n";
		file << "A0: " << params_list[i].sbparams.A0 << "\n";
		file << "a / beta / v: " << params_list[i].spacings[0] << "\n";
		file << "average cluster size (fraction of lattice): " << avg_clusters[i] << "\n";
		file << "Eq time: " << timers.times[i] << std::endl << std::endl << std::endl;
	}
	////output average step time as a function of lattice size
	//int lattice_size, prev_lattice_size, num_times = 0;
	//std::vector<double> avg_eq_times = {0.0};
	//prev_lattice_size = params_list[0].lengths[0] * params_list[0].lengths[1];
	//for (int i = 0; i < timers.size; ++i){
	//	lattice_size = params_list[i].lengths[0] * params_list[i].lengths[1];
	//	if(lattice_size == prev_lattice_size){
	//		avg_eq_times[avg_eq_times.size() - 1] += timers.times[i];
	//		++num_times;
	//	}
	//	else {
	//		std::cout << "lattice size: " << prev_lattice_size << std::endl;
	//		avg_eq_times[avg_eq_times.size() - 1] /= (num_times * prev_lattice_size);
	//		avg_eq_times.push_back(0.0);
	//		prev_lattice_size = lattice_size;
	//	}
	//}
	//std::cout << "lattice size: " << lattice_size << std::endl;
	//avg_eq_times[avg_eq_times.size() - 1] /= (num_times * lattice_size);
	//for (int i = 0; i < avg_eq_times.size(); ++i){
	//	std::cout << avg_eq_times[i] << "\n";
	//}
}

void process_autocorrelation_results(std::vector<std::vector<double>> mags, std::vector<class_mc_params> params_list) {
	//take a list of autocorrelation functions, match them with their input params, and output to csv for easy plotting
	std::ofstream file;
	file.open("autocorrelation_results.csv");
	std::vector<std::vector<double>> func;
	double f0;
	//calculate autocorrelation function
	for (int i = 0; i < mags.size(); ++i) {
		int tmax = mags[i].size()/2;
		func.push_back({});
		for (int dt = 0; dt < tmax; ++dt){
			func[i].push_back(0.0);
			for (int t0 = 0; t0 < tmax; ++t0){
				func[i][dt] += mags[i][t0 + dt]*mags[i][t0];
			}
			if (dt == 0) {
				f0 = func[i][0];
				func[i][0] = 1;
			}
			else {
				func[i][dt] /= f0;
			}
		}
/*
		for (int t = 1; t < tmax; ++t) {
			double m0 = 0, mt = 0, both = 0;
			for (int tprime = 0; tprime < tmax - t; ++tprime) {
				both += mags[i][tprime] * mags[i][tprime + t];
				m0 += mags[i][tprime];
				mt += mags[i][tprime + t];
			}
			func[i].push_back(0.1 / (tmax - t) * both - 0.01 / (tmax - t) / (tmax - t) * m0 * mt);
		}
*/
	}
	for (int i = 0; i < params_list.size(); ++i) {
		file << "Size," << params_list[i].lengths[0] << "," << params_list[i].lengths[0] << std::endl;
		file << "Beta," << params_list[i].beta << std::endl;
		file << "Delta," << params_list[i].sbparams.delta << std::endl;
		file << "autocorr," << vec2str(func[i]) << std::endl << std::endl << std::endl;
	}
}

std::vector<double> parallel_tempering(std::vector<LongRangeWolff2D>::iterator  wolffs_in_begin, std::vector<LongRangeWolff2D>::iterator wolffs_in_end) {
	int num_switches = std::distance(wolffs_in_begin, wolffs_in_end) - 1;
	std::vector<double> ratios(num_switches);
	std::vector<double>::iterator space = ratios.begin();
	double temp_beta, factor;
	for (;wolffs_in_begin != wolffs_in_end - 1; ++wolffs_in_begin, ++space){
		factor = exp(wolffs_in_begin->calc_E()*wolffs_in_begin->get_beta() - (wolffs_in_begin+1)->calc_E()*(wolffs_in_begin+1)->get_beta());
		if (drand1_() < factor){
			
			temp_beta= wolffs_in_begin->get_beta();
			wolffs_in_begin->set_beta((wolffs_in_begin+1)->get_beta());
			(wolffs_in_begin+1)->set_beta(temp_beta);
			*space = factor > 1 ? 1.0 : factor;
		}
	}
	return ratios;
}











