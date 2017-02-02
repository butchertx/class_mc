#include <IsingLattice2D.h>
#include <cmath>
#include <stdio.h>
#include <MemTimeTester.h>
#include <LongRangeWolff2D.h>
#include <class_mc_io.h>

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
	std::cout << params.to_string();


	//Set up lattice, wolff, and measurement objects
	IsingLattice2D test_lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
	LongRangeWolff2D wolff = LongRangeWolff2D(&test_lat, &params);
	class_mc_measurements results;

	
	//Set flags within wolff algorithm
	if (params.alg.compare("long_range_cluster") == 0) { wolff.set_alg_long_range_cluster(); }
	else if (params.alg.compare("nearest_neighbor_cluster") == 0) { wolff.set_alg_short_range_cluster(); }
}