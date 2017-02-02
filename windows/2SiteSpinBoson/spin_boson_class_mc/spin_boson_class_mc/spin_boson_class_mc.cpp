#include <iostream>
#include <fstream>
#include "MemTimeTester.h"
#include "IsingLattice2D.h"
#include "LongRangeWolff2D.h"
#include "class_mc_io.h"
#include "Matrix.h"
extern "C" {
#include "random.h"
}

int main(int argc, char* argv[]) {

	//Get input parameters
	class_mc_params params;
	std::ifstream infile;
	std::cout << "Program name: " << argv[0] << "\n";
	std::cout << "Input file: " << argv[1] << "\n\n";
	infile.open(argv[1]);
	read_input_ising(&infile, &params);
	read_input_spin_boson(&infile, &(params.sbparams));
	apply_spin_boson_params(&params);
	std::cout << "Parameters: \n" << params.to_string();
	write_params(params);



	//Set up lattice, wolff, and measurement objects
	IsingLattice2D test_lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
	LongRangeWolff2D wolff = LongRangeWolff2D(&test_lat, &params);
	std::cout << "Interactions: \n";
	wolff.print_interactions();
	class_mc_measurements results;
	MemTimeTester timer;
	results.names = { "steps", "times", "sx", "sz", "mag", "mag2", "mag4", "ising_order" };
	results.values = { {}, {}, {}, {}, {}, {}, {}, {} };

	results.func_names = { "corr_x" };
	results.functions = { {} };
	results.function_num_measures = { 0 };

	
	//Set flags within wolff algorithm
	if (params.alg.compare("long_range_cluster") == 0) { wolff.set_alg_long_range_cluster(); }
	else if (params.alg.compare("nearest_neighbor_cluster") == 0) { wolff.set_alg_short_range_cluster(); }

	
	//Set up parallel tempering lattices
	int num_ptemp = 1;
	std::cout << "Number of parallel tempering lattices: " << num_ptemp << "\n";
	std::vector<LongRangeWolff2D> wolffs = { wolff };


	//Run simulation
	timer.flag_start_time("simulation");
		//equilibrate
	for (int i = 0; i < params.eq_time; ++i) {
		for (int ptemp = 0; ptemp < num_ptemp; ++ptemp) {
			wolffs[ptemp].step();
		}
	}
		//measure
	double mag;
	for (int dump = 0; dump < params.max_dumps; ++dump) {
		for (int measure = 0; measure < params.measures_per_dump; ++measure) {
			for (int step = 0; step < params.steps_per_measure; ++step) {
				for (int ptemp = 0; ptemp < num_ptemp; ++ptemp) {
					wolffs[ptemp].step();
				}
			}
			//make measurements
			results.record("steps", (double)(dump*params.measures_per_dump*params.steps_per_measure + measure*params.steps_per_measure + params.eq_time));
			results.record("times", timer.get_running_time("simulation"));
			results.record("sx", wolffs[0].calc_sx());
			results.record("sz", wolffs[0].calc_sz());
			mag = wolffs[0].calc_mag();
			results.record("mag", mag);
			results.record("mag2", mag * mag);
			results.record("mag4", mag*mag*mag*mag);
			results.record("ising_order", (wolffs[0].calc_corr(0))[1]);
			results.record("corr_x", (wolffs[0].calc_corr(0)));
			
		}
		write_state(dump, wolffs[0].get_lat(), 1.0);
		std::vector<int> real_steps;
		std::vector<double> bad_steps = results.get_vals("steps");
		for (int i = 0; i < bad_steps.size(); ++i) {
			real_steps.push_back((int)bad_steps[i]);
		}
		write_outputs_var(dump, results);
		std::cout << "Dump " << dump + 1 << " out of " << params.max_dumps << "\n";
	}
	write_final_outputs(results, params.max_dumps);
	timer.flag_end_time("simulation");

	timer.print_timers();
}