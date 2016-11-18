#include "LongRangeWolff2D.h"

LongRangeWolff2D::LongRangeWolff2D(IsingLattice2D* lat_in, class_mc_params* params_in, std::string model_type) {
	lat = *lat_in;
	params = *params_in;
	interactions = Matrix(lat.get_Lx(), lat.get_Ly());
	cluster = IntMatrix(lat.get_Lx(), lat.get_Ly());
	cluster.fill(0);

	//generate interactions matrix
	if (model_type.compare("mean_field") == 0) {
		set_mean_field_model();
	}
	else {
		std::cout << model_type << " is not a valid lattice model\n";
	}

	//set flags
	LONG_RANGE_CLUSTER = false;
}

void LongRangeWolff2D::set_mean_field_model() {
	//matrix giving the action contributed by a bond of dimension (i,j)
	double int_val = - params.Js[0] * params.beta / lat.get_N();
	interactions.fill(int_val);
}

void LongRangeWolff2D::set_spin_boson_model() {

}

void LongRangeWolff2D::step() {
	//1. determine a "seed" spin
	//2. go to all interacting spins and test to add to buffer and cluster
	//		- use long range procedure in Int. J. Mod. Phys. C 6, 359-370 (1995).
	//3. when buffer is empty, flip cluster

	//1. get seed
	spin seed;
	seed.x = (int) (lat.get_Lx() * drand1_());
	seed.y = (int) (lat.get_Ly() * drand1_());
	seed.proj = lat.get_spin(seed.x, seed.y);
	cluster.setval(seed.x, seed.y, 1);
	//std::cout << "seed spin: (" << seed.x << "," << seed.y << ") " << seed.proj << "\n";


	//2. add to buffer, then set loop to go through buffer
	test_spins(seed);
	while (buffer.size() != 0) {
		seed = buffer.back();
		buffer.pop_back();
		test_spins(seed);
	}

	//std::cout << "cluster: \n" << cluster.to_string() << "\n";

	//3. flip cluster
	for (int x = 0; x < lat.get_Lx(); ++x) {
		for (int y = 0; y < lat.get_Ly(); ++y) {
			if (cluster.getval(x, y) == 1) {
				lat.flip_spin(x, y);
			}
		}
	}
	//clear cluster
	cluster.fill(0);
}

double LongRangeWolff2D::calc_mag() {
	double mag = 0;
	for (int i = 0; i < lat.get_Lx(); ++i) {
		for (int j = 0; j < lat.get_Ly(); ++j) {
			mag += lat.get_spin(i, j);
		}
	}
	return mag / lat.get_N();
}

double LongRangeWolff2D::calc_E_mean_field_model() {
	int N = lat.get_N();
	double M = calc_mag() * N;
	double E = 0;
	for (int i = 0; i < lat.get_Lx(); ++i) {
		for (int j = 0; j < lat.get_Ly(); ++j) {
			E += - params.Js[0] * .5 / lat.get_N() * (M - lat.get_spin(i, j)) * lat.get_spin(i,j);
		}
	}
	return E;
}

void LongRangeWolff2D::test_spins(spin seed){
	//choose a random number.  Then determine how far away the next added spin is
	//test next added spin for validity
	//continue to next candidate

	int N = lat.get_N(), Lx = lat.get_Lx(), Ly = lat.get_Ly(), i = 0, j = 0, dist = 0;
	double prob, rando, K = interactions.getval(0, 0);
	spin newspin;
	int x = 0, y = 0;
	//current version is specific to the mean field model
	if (LONG_RANGE_CLUSTER) {

		while (i < Lx || j < Ly) {
			rando = drand1_();
			dist = ceil(log(1 - rando) / K);
			i = (i + dist);
			j = (j + (dist / Lx));
			x = (seed.x + i) % Lx;
			y = (seed.y + j) % Ly;
			if (cluster.getval(x, y) == 0 && lat.get_spin(x, y) == seed.proj) {
				newspin.x = x;
				newspin.y = y;
				newspin.proj = seed.proj;
				cluster.setval(x, y, 1);
				buffer.push_back(newspin);
			}
		}
	}
	else {

		//do it the naive way - choose a random number then go through
		//lattice until you hit the right probability range and correct spin
		//this will make things easier for adding spin boson model
		for (int i = 0; i < Lx; ++i) {
			for (int j = 0; j < Ly; ++j) {
				x = (seed.x + i) % Lx;
				y = (seed.y + j) % Ly;
				prob = 1 - exp(2 * interactions.getval(i, j) * seed.proj * lat.get_spin(x, y));
				if (drand1_() < prob && cluster.getval(x, y) == 0) {
					newspin.x = x;
					newspin.y = y;
					newspin.proj = lat.get_spin(x, y);
					cluster.setval(x, y, 1);
					buffer.push_back(newspin);
				}
			}
		}
	}
}
