#define _USE_MATH_DEFINES
#include <cmath>
#include <omp.h>
#include <algorithm>
#include "LongRangeWolff2D.h"

LongRangeWolff2D::LongRangeWolff2D(IsingLattice2D* lat_in, class_mc_params* params_in) {
	lat = *lat_in;
	params = *params_in;


	//generate member objects
	interactions = Matrix(lat.get_Lx(), lat.get_Ly());
	cluster = IntMatrix(lat.get_Lx(), lat.get_Ly());
	cluster.fill(0);
	cluster_size = 0;
	set_model();
	set_cumulative_probs(lat.get_Lx(), lat.get_Ly());
	mag = calc_sz()*lat.get_N();
	E = calc_E();

	//set flags
	LONG_RANGE_CLUSTER = false;
	NEAREST_NEIGHBOR_CLUSTER = false;
	TIMERS = false;
	VERIFY_TESTING = false;
}

//output state of members for checking
//should only be used for testing
void LongRangeWolff2D::output_state() {
	std::cout << "\n\nOUTPUTTING CURRENT WOLFF STATE\n\n";

	std::cout << "Flags:\n";
	std::cout << "Long range cluster building: " << LONG_RANGE_CLUSTER << "\n";
	std::cout << "Nearest neighbor cluster building: " << NEAREST_NEIGHBOR_CLUSTER << "\n";

	std::cout << "\nState\n" << lat.to_string() << "\n";

	std::cout << "\nCluster and buffer info\n";
	std::cout << "Cluster: \n" << cluster.to_string();
	std::cout << "Buffer: \n";
	for (int i = 0; i < buffer.size(); ++i) {
		std::cout << buffer[i].x << "," << buffer[i].y << "  ";
	}
	std::cout << "\nNot Cluster: \n";
	for (int i = 0; i < not_cluster.size(); ++i) {
		std::cout << not_cluster[i].x << "," << not_cluster[i].y << "  ";
	}
	std::cout << "\nNumber in Cluster: " << cluster_size << "\n";

}

void LongRangeWolff2D::set_model() {
	if (params.model.compare("mean_field") == 0) {
		set_mean_field_model();
	}
	else if (params.model.compare("spin_boson") == 0) {
		set_spin_boson_model();
	}
	else {
		std::cout << params.model << " is not a valid lattice model\n";
	}
}

void LongRangeWolff2D::set_mean_field_model() {
	//matrix giving the action contributed by a bond of dimension (i,j)
	double int_val = - params.Js[0] * params.beta / lat.get_N();
	interactions.fill(int_val);
}

void LongRangeWolff2D::set_spin_boson_model() {
	//model is translationally invariant, so (i,j) will represent a distance in the (spacial, imaginary time) direction
	interactions.fill(0.0);
	double g = params.sbparams.g, A = params.sbparams.A0, gamma = params.Js[1], J = params.Js[0], v = params.sbparams.v, a = params.spacings[0], tc = params.spacings[1], x, y;
	int Nt = params.lengths[1];
	for (int i = 0; i < interactions.get_dimx(); ++i) {
		for (int j = 0; j < interactions.get_dimy(); ++j) {
			if (i == 0) {
				//temporal nearest neighbor interaction
				if (j == 1 || j == interactions.get_dimy() - 1) {
					interactions.setval(i, j, -gamma - 0.5*A*M_PI*M_PI / Nt / Nt / sin(M_PI / Nt) / sin(M_PI / Nt));
				}
				//temporal self-interaction
				if (j > 1 && j < interactions.get_dimy() - 1) {
					interactions.setval(i, j, -0.5*A*M_PI*M_PI / Nt / Nt / sin(M_PI * j / Nt) / sin(M_PI * j / Nt));
				}
			}
			else if (j == 0) {
				//spacial nearest neighbor interactions
				if (i == 1) {
					interactions.setval(i, j, -J + 0.5*A*M_PI*M_PI / Nt / Nt / sinh(M_PI * a / Nt / v / tc) / sinh(M_PI * a / Nt / v / tc));
				}
				//spacial same-time long range interactions
				if (i > 1) {
					interactions.setval(i, j, 0.5*A*M_PI*M_PI / Nt / Nt / sinh(M_PI * a * i / Nt / v / tc) / sinh(M_PI * a * i / Nt / v / tc));
				}
			}
			else {
				//different site, different time interactions
				x = M_PI*j / Nt;
				y = M_PI*i*a / v / Nt / tc;
				interactions.setval(i, j, -M_PI*M_PI*0.5*A / Nt / Nt*(sin(x)*sin(x)*cosh(y)*cosh(y) - cos(x)*cos(x)*sinh(y)*sinh(y)) /
					((sin(x)*sin(x)*cosh(y)*cosh(y) + cos(x)*cos(x)*sinh(y)*sinh(y))*(sin(x)*sin(x)*cosh(y)*cosh(y) + cos(x)*cos(x)*sinh(y)*sinh(y))));
			}
		}
	}
}

//use interaction matrix to set up cumulative bond probabilitites
//inputs allow incorporation of a cutoff
//simple cluster forming: start at seed and go across in the y direction, then jump to next x neighbor
void LongRangeWolff2D::set_cumulative_probs(int dimx, int dimy) {
	double sum = 0.0;
	for (int i = 0; i < dimx; ++i) {
		for (int j = 0; j < dimy; ++j) {
			sum += abs(interactions.getval(i, j));
			cumulative_probs.push_back(1 - exp(-2*sum));
		}
	}
}

void LongRangeWolff2D::fill_not_cluster() {
	spin temp_spin;
	not_cluster = {};
	for (int i = 0; i < lat.get_Lx(); ++i) {
		for (int j = 0; j < lat.get_Ly(); ++j) {
			temp_spin.x = i;
			temp_spin.y = j;
			temp_spin.proj = lat.get_spin(i, j);
			if (cluster.getval(i, j) == 0) {
				not_cluster.push_back(temp_spin);
			}
		}
	}
}

void LongRangeWolff2D::step() {
	if (VERIFY_TESTING) {
		std::cout << "\n\n*********************STARTING STEP*************************\n\n";
		output_state();
	}

	//1. determine a "seed" spin
	//2. go to all interacting spins and test to add to buffer and cluster
	//		- use long range procedure in Int. J. Mod. Phys. C 6, 359-370 (1995).
	//3. when buffer is empty, flip cluster
	cluster_size = 0;
	cluster_mag = 0;
	dS = 0;
	//1. get seed
	spin seed;
    seed.x = (int) (lat.get_Lx() * drand1_());
    seed.y = (int) (lat.get_Ly() * drand1_());
    seed.proj = lat.get_spin(seed.x, seed.y);
    cluster.setval(seed.x, seed.y, 1);
	++cluster_size;
	cluster_mag += seed.proj;
	if (!LONG_RANGE_CLUSTER && !NEAREST_NEIGHBOR_CLUSTER) {
		fill_not_cluster();
	}


	if (VERIFY_TESTING) {
		std::cout << "\n\n*********************SEED SPIN CHOSEN*************************\n\n";
		std::cout << "seed spin: (" << seed.x << "," << seed.y << ") " << seed.proj << "\n";
		output_state();
	}


	//2. add to buffer, then set loop to go through buffer
	test_spins(seed);
	while (buffer.size() != 0) {
		seed = buffer.back();
		cluster_mag += seed.proj;
		++cluster_size;
		buffer.pop_back();
		test_spins(seed);
	}


	if (VERIFY_TESTING) {
		std::cout << "\n\n*********************CLUSTER CHOSEN*************************\n\n";
		output_state();
	}

	//3. flip cluster according to a probability check based on the magnetization change
	//prob is 1 if h*mag is negative, exp(-2*h*mag) if h*mag is positive
	//satisfies detailed balance: this is the "acceptance probability", rest was "selection probability"
	if (drand1_() < exp(-2 * params.h * cluster_mag)) {
		for (int x = 0; x < lat.get_Lx(); ++x) {
			for (int y = 0; y < lat.get_Ly(); ++y) {
				if (cluster.getval(x, y) == 1) {
					lat.flip_spin(x, y);
					cluster.setval(x, y, 0);
				}
			}
		}
		if (VERIFY_TESTING) {
			std::cout << "\n\n*********************CLUSTER FLIPPED WITH PROBABILITY " << exp(-2 * params.h * cluster_mag) << "*************************\n\n";
		}
	}
	else {
		if (VERIFY_TESTING) {
			std::cout << "\n\n*********************CLUSTER NOT FLIPPED WITH PROBABILITY " << exp(-2 * params.h * cluster_mag) << "*************************\n\n";
		}
	}

	mag -= 2*cluster_mag;

}

double LongRangeWolff2D::calc_mag() {
	int m = 0;
	for (int i = 0; i < lat.get_Lx(); ++i) {
		for (int j = 0; j < lat.get_Ly(); ++j) {
			m += lat.get_spin(i, j);
		}
	}
	return m < 0 ? -((double)m) / lat.get_N() : ((double)m) / lat.get_N();
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

double LongRangeWolff2D::calc_E() {
	//essentially, the time-averaged energy of a given set of quantum fluctuations.  Calculate the action then divide by beta
	double S = 0;
	int s1;
	for (int i = 0; i < lat.get_Lx(); ++i) {
		for (int j = 0; j < lat.get_Ly(); ++j) {
			s1 = lat.get_spin(i, j);
			for (int m = i; m < lat.get_Lx(); ++m) {
				for (int n = j; n < lat.get_Ly(); ++n) {
					S += s1 * lat.get_spin(m, n) * interactions.getval(m - i, n - j);
				}
			}
		}
	}

	return S / params.beta;
}

std::vector<double> LongRangeWolff2D::calc_corr(int dimension) {
	//calculate the correlation function along x (if dim = 0) or y (if dim = 1)
	if (dimension != 0 && dimension != 1) {
		std::cout << "Error: invalid correlation measurement attempted\n";
		return{};
	}
	else {
		int length_avg = dimension == 1 ? lat.get_Lx() : lat.get_Ly();
		int length_corr = dimension == 1 ? lat.get_Ly() : lat.get_Lx();
		std::vector<double> corr(length_corr, 0.0);
		for (int i = 0; i < length_avg; ++i) {
			for (int j = 0; j < length_corr; ++j) {
				corr[j] += dimension == 1 ? lat.get_spin(i, 0)*lat.get_spin(i, j) : lat.get_spin(0, i)*lat.get_spin(j, i);
			}
		}
		for (int i = 0; i < corr.size(); ++i) {
			corr[i] = corr[i] / length_avg;
		}
		return corr;
	}
}

std::vector<double> LongRangeWolff2D::calc_corr_slow() {
	//calculate the correlation function in both dimensions
	std::vector<double> corr(lat.get_N(), 0.0);
	for (int i = 0; i < lat.get_Lx(); ++i) {
		for (int j = 0; j < lat.get_Ly(); ++j) {
			for (int r = 0; r < lat.get_Lx(); ++r) {
				for (int s = 0; s < lat.get_Ly(); ++s) {
					corr[i*lat.get_Ly() + j] += lat.get_spin(r, s)*lat.get_spin((r + i) % lat.get_Lx(), (s + j) % lat.get_Ly());
				}
			}
		}
	}
	for (int i = 0; i < corr.size(); ++i) {
		corr[i] = corr[i] / lat.get_N();
	}
	return corr;
}

double LongRangeWolff2D::calc_sx() {
	//assume a quantum-classical mapping.  Then, <flux>_i = 1/beta * sum_t{1/2(1 - <s(t + dt)s(t)>)}
	double final_result = 0;
	for (int i = 0; i < lat.get_Lx(); ++i) {
		double site_mag = 0;
		for (int j = 0; j < lat.get_Ly() - 1; ++j) {
			site_mag += 0.5*(1 - lat.get_spin(i, j)*lat.get_spin(i, j + 1));
		}
		site_mag += 0.5*(1 - lat.get_spin(i, 0)*lat.get_spin(i, lat.get_Ly() - 1));
		final_result += site_mag / lat.get_Ly();
	}
	return final_result / lat.get_Lx();
}

double LongRangeWolff2D::calc_sz() {
	//similar to <flux>, but instead find the magnetization for each point in time, then average and return the absolute value
	double final_result = 0.0;
	for (int j = 0; j < lat.get_Ly(); ++j) {
		double time_mag = 0.0;
		for (int i = 0; i < lat.get_Lx(); ++i) {
			time_mag += lat.get_spin(i, j);
		}
		final_result += time_mag / lat.get_Lx();
	}
	return final_result / lat.get_Ly();
}

void LongRangeWolff2D::test_spins(spin seed){
	//choose a random number.  Then determine how far away the next added spin is
	//test next added spin for validity
	//continue to next candidate

	int N = lat.get_N(), Lx = lat.get_Lx(), Ly = lat.get_Ly(), i = 0, j = 0, k = 0, dist = 0;
	double prob, rando;
	spin newspin;
	int x = 0, y = 0;
	//use look-up table and long-range wolff cluster building
	if (LONG_RANGE_CLUSTER) {
		//std::vector<double>::iterator start = cumulative_probs.begin(), site = cumulative_probs.begin(), finish = cumulative_probs.end();
		////std::cout << "Testing long range cluster building: \n";
		////std::cout << "lattice start: \n" << lat.to_string() << "\n";
		////std::cout << "seed site: {" << seed.x << ", " << seed.y << "}, seed projection: " << seed.proj << "\n";
		//do {
		//	rando = drand1_()*(1 - *site) + *site;
		//	//std::cout << "current add site: " << std::distance(start, site) << "\n";
		//	//std::cout << "scaled random number: " << rando << "\n";
		//	site = std::lower_bound(start, finish, rando);
		//	dist = std::distance(start, site);
		//	//std::cout << "new add site: " << dist << " away from seed location\nSpin added?";
		//	i = (seed.x + dist / Ly) % Lx;
		//	j = (seed.y + dist % Ly) % Ly;
		//	if (site != finish && cluster.getval(i, j) == 0 && interactions.getval(i, j)*seed.proj*lat.get_spin(i, j) < 0) {
		//		newspin.x = i;
		//		newspin.y = j;
		//		newspin.proj = lat.get_spin(i, j);
		//		buffer.push_back(newspin);
		//		cluster.setval(i, j, 1);
		//		++cluster_size;
		//		//std::cout << " yes\n\n";
		//	}
		//	//else std::cout << " no\n\n";
		//} while (site != finish);


		//implement Blote's long range cluster forming
		//here the spin-boson action is used; would need to alter for different forms of the action
		//also, only looks in time dimension from seed spin
		//algorithm: choose a random number, calculate distance from seed spin and add corresponding spin to the cluster until length is exhausted
		//start with a cutoff at a distance of 2*log_10(Nt) or 10 spins, whichever is larger and consider each spin in turn
		int cutoff = 5 > 2 * log10(Ly) ? 5 : 2 * log10(Ly);
		if (cutoff > Ly/2) { cutoff = Ly/2; }
		if (VERIFY_TESTING) {
			std::cout << "Cutoff distance: " << cutoff << "\n";

			//check spins within cutoff distance
			std::cout << "Seed coordinates: (" << seed.x << ", " << seed.y << ")\n";
		}
		for (int i = Ly - cutoff; i < Ly + cutoff; ++i) {
			y = (seed.y + i) % Ly;
			if (cluster.getval(seed.x, y) == 0) {
				if (VERIFY_TESTING) {
					std::cout << "Checking time coordinate: " << y << "\n";
				}
				if (drand1_() < 1 - exp(2 * interactions.getval(0, i % Ly)*seed.proj*lat.get_spin(seed.x, y))) {
					newspin.x = seed.x;
					newspin.y = y;
					newspin.proj = lat.get_spin(seed.x, y);
					buffer.push_back(newspin);
					cluster.setval(newspin.x, newspin.y, 1);
					if (VERIFY_TESTING) {
						std::cout << "Spin at tau bar = " << y << " Added\n";
					}
				}
			}
		}

		//now do the long-range cluster building
		j = cutoff;
		while (j < Ly - cutoff) {
			rando = atan(1 / (tan(M_PI*(.5*(Ly + 1) - j) / Ly) + Ly * log(1 - drand1_()) / params.sbparams.A0 / M_PI));
			k = rando < 0 ? ceil(Ly / M_PI * (rando + M_PI) - .5) : ceil(Ly / M_PI * rando - .5);
			y = (seed.y + k + Ly) % Ly;
			if (k < Ly - cutoff && cluster.getval(seed.x, y) == 0) {
				//std::cout << "Inside if statement: (j, k, y) = (" << j << ", " << k << ", " << y << ")\n";
				newspin.x = seed.x;
				newspin.y = y;
				newspin.proj = lat.get_spin(seed.x, y);
				buffer.push_back(newspin);
				cluster.setval(newspin.x, newspin.y, 1);
			}
			j = k;
		}

	}
	else if (NEAREST_NEIGHBOR_CLUSTER) {
		//Only look at nearest neighbors, and only care about two possible couplings
		std::cout << "Dm and DE are not implemented\n";
		double Jx = params.Js[0], Jy = params.Js[1];
		x = seed.x;
		y = seed.y;

		//test x neighbors
		if (seed.proj == lat.get_spin((x + Lx - 1) % Lx, y) && cluster.getval((x + Lx - 1) % Lx, y) == 0 && drand1_() < (1 - exp(-2 * Jx))) {
			newspin.proj = seed.proj;
			newspin.x = (x + Lx - 1) % Lx;
			newspin.y = y;
			cluster.setval((x + Lx - 1) % Lx, y, 1);
			++cluster_size;
			buffer.push_back(newspin);
		}
		if (seed.proj == lat.get_spin((x + 1) % Lx, y) && cluster.getval((x + 1) % Lx, y) == 0 && drand1_() < (1 - exp(-2 * Jx))) {
			newspin.proj = seed.proj;
			newspin.x = (x + 1) % Lx;
			newspin.y = y;
			cluster.setval((x + 1) % Lx, y, 1);
			++cluster_size;
			buffer.push_back(newspin);
		}

		//test y neighbors
		if (seed.proj == lat.get_spin(x, (y + Ly - 1)%Ly) && cluster.getval(x, (y + Ly - 1) % Ly) == 0 && drand1_() < (1 - exp(-2 * Jy))) {
			newspin.proj = seed.proj;
			newspin.x = x;
			newspin.y = (y + Ly - 1) % Ly;
			cluster.setval(x, (y + Ly - 1) % Ly, 1);
			++cluster_size;
			buffer.push_back(newspin);
		}
		if (seed.proj == lat.get_spin(x, (y + 1) % Ly) && cluster.getval(x, (y + 1) % Ly) == 0 && drand1_() < (1 - exp(-2 * Jy))) {
			newspin.proj = seed.proj;
			newspin.x = x;
			newspin.y = (y + 1) % Ly;
			cluster.setval(x, (y + 1) % Ly, 1);
			++cluster_size;
			buffer.push_back(newspin);
		}
	}
	else {

		//do it the naive way - choose a random number then go through
		//lattice until you hit the right probability range and correct spin
		//this will make things easier for adding spin boson model
		//includes the sz-field term h
		std::vector<spin>::iterator newspin_it = not_cluster.begin();
		spin temp_spin;
		while (newspin_it != not_cluster.end()) {
			prob = 1 - exp(2 * (interactions.getval(abs(newspin_it->x - seed.x), abs(newspin_it->y - seed.y)) * (seed.proj)) * (newspin_it->proj));
			if (drand1_() < prob) {
				buffer.push_back(*newspin_it);
				cluster.setval(newspin_it->x, newspin_it->y, 1);
				//if (VERIFY_TESTING) {
				//	std::cout << "probability: " << prob << " of adding new spin to cluster (negative means 0 chance)\n";
				//	std::cout << "seed spin: (" << seed.x << "," << seed.y << "," << seed.proj << ")\n";
				//	std::cout << "new spin: (" << newspin_it->x << "," << newspin_it->y << "," << newspin_it->proj << ")\n";
				//}
				*newspin_it = not_cluster.back();
				if (newspin_it == not_cluster.end() - 1) {
					not_cluster.pop_back();
					newspin_it = not_cluster.end();
				}
				else { not_cluster.pop_back(); }
			}
			else { ++newspin_it; }
		}
		//fill_not_cluster();
	}
}
