#define _USE_MATH_DEFINES
#include <cmath>
#include <omp.h>
#include <algorithm>
#include "LongRangeWolff2D.h"

double real_trigamma_cmplxarg(double x, double y, int m, int l){
    //compute the real part of the trigamma function at z = x + iy using
    //my accelerated series formula
    //argument should have positive real part (x > 0)
    //m gives the number of zeta function terms to use
    //l gives the number of terms to use in the residual series
    double result = 0;

    //if x is small, the result will be large and thus, inaccurate.  Use the
    //polygamma shift formula to give a larger value of x for the computation
    if(x < 1000){
        return ((x*x - y*y)/(x*x + y*y)/(x*x + y*y)) + real_trigamma_cmplxarg(x + 1, y, m, l);
    }
    else{
        double phase, ypow = 1.0, xpow, denom;
        //compute finite sum of Hurwitz zeta functions
        for (int i = 1; i < m; ++i){
            phase = 2*(i/2) == i ? -1.0 : 1.0;
            result += phase*(2*i - 1)*ypow*gsl_sf_hzeta(2*i, x);
            ypow = ypow*y*y;
        }
        //compute the infinite sum of residuals
        phase = 2*(m/2) == m ? 1.0 : -1.0;
        for (int n = 0; n < l; ++n){
            xpow = pow(x + n, 2*m);
            denom = ((x + n)*(x + n) + y*y)*((x + n)*(x + n) + y*y);
            result += phase*ypow*((2*m - 1)*y*y + (2*m + 1)*(x + n)*(x + n)) / xpow / denom;
        }
        return result;
    }
}

LongRangeWolff2D::LongRangeWolff2D(IsingLattice2D* lat_in, class_mc_params* params_in) {
	lat = *lat_in;
	params = *params_in;


	//generate member objects
	interactions = Matrix(lat.get_Lx(), lat.get_Ly());
	cluster = IntMatrix(lat.get_Lx(), lat.get_Ly());
	cluster.fill(0);
	cluster_size = 0;
	set_model();
	//set_cumulative_probs(lat.get_Lx(), lat.get_Ly());
	mag = calc_sz()*lat.get_N();
	E = calc_E();
	cutoff = 5 > log10(lat.get_Ly()) ? 5 : log10(lat.get_Ly());
	tau_x = floor(lat.get_Ly() / M_PI * atan(tanh(M_PI * params.lengths[0] / params.beta / params.sbparams.v)));
	//std::cout << "tau_x: " << tau_x << "\n";
	if (cutoff > lat.get_Ly() / 2) { cutoff = lat.get_Ly() / 2; }

	//set flags
	LONG_RANGE_CLUSTER = false;
	NEAREST_NEIGHBOR_CLUSTER = false;
	TIMERS = false;
	VERIFY_TESTING = false;

	//set up interaction sum
	if (lat.get_Lx() > 1) {
		//interaction_sum(i, j) = integral(abs(chi)) at i, j + 1/2
		double xc = M_PI * params.lengths[0] / params.beta / params.sbparams.v;
		double real_taux = lat.get_Ly() * atan(tanh(xc)) / M_PI;
		outsum = -sin(2 * atan(tanh(xc))) / (cos(2 * atan(tanh(xc))) - cosh(2 * xc));//integral of abs(chi) from 0 to real_taux
		//std::cout << "outsum: " << outsum << "\n";
		//std::cout << "taux: " << real_taux << "\n";
		interaction_sum = Matrix(lat.get_Lx(), lat.get_Ly());
		double coshx = cosh(2 * xc);
		for (int j = 0; j < lat.get_Ly(); ++j) {
			interaction_sum.setval(0, j, sin(2 * M_PI*(j + 0.5) / lat.get_Ly()) / (cos(2 * M_PI*(j + 0.5) / lat.get_Ly()) - 1));
		}
		for (int j = 0; j < lat.get_Ly(); ++j) {
			if (j + 0.5 < real_taux) {
				interaction_sum.setval(1, j, -sin(2 * M_PI*(j + 0.5) / lat.get_Ly()) / (cos(2 * M_PI*(j + 0.5) / lat.get_Ly()) - coshx));
			}
			else if (j + 0.5 >= real_taux && j + 0.5 < params.beta - real_taux) {
				interaction_sum.setval(1, j, 2 * outsum + sin(2 * M_PI*(j + 0.5) / lat.get_Ly()) / (cos(2 * M_PI*(j + 0.5) / lat.get_Ly()) - coshx));
			}
			else {
				interaction_sum.setval(1, j, 4 * outsum - sin(2 * M_PI*(j + 0.5) / lat.get_Ly()) / (cos(2 * M_PI*(j + 0.5) / lat.get_Ly()) - coshx));
			}
		}
	}
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
	std::cout << "\nCluster Mag: " << cluster_mag << "\n";

}

void LongRangeWolff2D::set_model() {
	set_spin_boson_model();
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

void LongRangeWolff2D::write_resample_interactions() {
	//only does the long-range terms (proportional to alpha)
	Matrix resample_interactions = Matrix(lat.get_Lx(), lat.get_Ly());

	resample_interactions.fill(0.0);
	double g = params.sbparams.g, A = 1.0, v = params.sbparams.v, a = params.spacings[0], tc = params.spacings[1], x, y;
	int Nt = params.lengths[1];
	for (int i = 0; i < resample_interactions.get_dimx(); ++i) {
		for (int j = 0; j < resample_interactions.get_dimy(); ++j) {
			if (i == 0) {
				//temporal self-interaction
				if (j >= 1 && j <= interactions.get_dimy() - 1) {
					resample_interactions.setval(i, j, -0.5*A*M_PI*M_PI / Nt / Nt / sin(M_PI * j / Nt) / sin(M_PI * j / Nt));
				}
			}
			else if (j == 0) {
				//spacial same-time long range interactions
				if (i >= 1) {
					resample_interactions.setval(i, j, 0.5*A*M_PI*M_PI / Nt / Nt / sinh(M_PI * a * i / Nt / v / tc) / sinh(M_PI * a * i / Nt / v / tc));
				}
			}
			else {
				//different site, different time interactions
				x = M_PI*j / Nt;
				y = M_PI*i*a / v / Nt / tc;
				resample_interactions.setval(i, j, -M_PI*M_PI*0.5*A / Nt / Nt*(sin(x)*sin(x)*cosh(y)*cosh(y) - cos(x)*cos(x)*sinh(y)*sinh(y)) /
					((sin(x)*sin(x)*cosh(y)*cosh(y) + cos(x)*cos(x)*sinh(y)*sinh(y))*(sin(x)*sin(x)*cosh(y)*cosh(y) + cos(x)*cos(x)*sinh(y)*sinh(y))));
			}
		}
	}


	std::ofstream file;
	file.open("resample_interactions.csv");
	file << resample_interactions.to_string() << "\n";
}

////use interaction matrix to set up cumulative bond probabilitites
////inputs allow incorporation of a cutoff
////simple cluster forming: start at seed and go across in the y direction, then jump to next x neighbor
//void LongRangeWolff2D::set_cumulative_probs(int dimx, int dimy) {
//	double sum = 0.0;
//	for (int i = 0; i < dimx; ++i) {
//		for (int j = 0; j < dimy; ++j) {
//			sum += abs(interactions.getval(i, j));
//			cumulative_probs.push_back(1 - exp(-2*sum));
//		}
//	}
//}

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

double LongRangeWolff2D::calc_sz_stagger() {
	//antiferromagnetic analogue of sz().  calc mag for each site, then multiply by (-1)^i and average
	double final_result = 0.0;
	for (int j = 0; j < lat.get_Ly(); ++j) {
		double time_mag = 0.0;
		for (int i = 0; i < lat.get_Lx(); ++i) {
			time_mag += (2*(i/2) == i) ? lat.get_spin(i, j) : -lat.get_spin(i, j);
		}
		final_result += time_mag / lat.get_Lx();
	}
	return (final_result < 0) ? -final_result / lat.get_Ly() : final_result / lat.get_Ly();
}

void LongRangeWolff2D::test_spins(spin seed){
	//choose a random number.  Then determine how far away the next added spin is
	//test next added spin for validity
	//continue to next candidate

	int N = lat.get_N(), Lx = lat.get_Lx(), Ly = lat.get_Ly(), i = 0, j = 0, k = 0, dist = 0;
	double prob, rando, q, shift, xc = M_PI*params.lengths[0]/params.beta/params.sbparams.v, alpha = params.sbparams.A0;
	spin newspin;
	int x = 0, y = 0;
	//use look-up table and long-range wolff cluster building
	if (LONG_RANGE_CLUSTER) {
		//implement Blote's long range cluster forming
		//here the spin-boson action is used; would need to alter for different forms of the action
		//also, only looks in time dimension from seed spin
		//algorithm: choose a random number, calculate distance from seed spin and add corresponding spin to the cluster until length is exhausted
		//start with a cutoff at a distance of 2*log_10(Nt) or 10 spins, whichever is larger and consider each spin in turn

		//if (VERIFY_TESTING) {
		//	std::cout << "Cutoff distance: " << cutoff << "\n";
		//	//check spins within cutoff distance
		//	std::cout << "Seed coordinates: (" << seed.x << ", " << seed.y << ")\n";
		//}
		if (Lx == 1) {
			for (int i = Ly - cutoff; i < Ly + cutoff; ++i) {
				y = (seed.y + i) % Ly;
				if (cluster.getval(seed.x, y) == 0) {
					//if (VERIFY_TESTING) {
					//	std::cout << "Checking time coordinate: " << y << "\n";
					//}
					if (drand1_() < 1 - exp(2 * interactions.getval(0, i % Ly)*seed.proj*lat.get_spin(seed.x, y))) {
						newspin.x = seed.x;
						newspin.y = y;
						newspin.proj = lat.get_spin(seed.x, y);
						buffer.push_back(newspin);
						cluster.setval(newspin.x, newspin.y, 1);
						//if (VERIFY_TESTING) {
						//	std::cout << "Spin at tau bar = " << y << " Added\n";
						//}
					}
				}
			}

			//now do the long-range cluster building
			j = cutoff;
			while (j < Ly - cutoff) {
				rando = atan(1 / (tan(M_PI*(.5*(Ly + 1) - j) / Ly) + Ly * log(1 - drand1_()) / alpha / M_PI));
				k = rando < 0 ? ceil(Ly / M_PI * (rando + M_PI) - .5) : ceil(Ly / M_PI * rando - .5);
				y = (seed.y + k + Ly) % Ly;
				//if (VERIFY_TESTING) {
				//	std::cout << "next long range activated lattice site: " << y << "\n";
				//}
				if (k < Ly - cutoff && cluster.getval(seed.x, y) == 0 && interactions.getval(0, k)*seed.proj*lat.get_spin(seed.x, y) < 0) {
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
		else if (Lx == 2) {
			for (int i = 0; i < 2; ++i) {
				for (int j = Ly - cutoff; j < Ly + cutoff; ++j) {
					y = (seed.y + j) % Ly;
					x = (seed.x + i) % 2;
					if (cluster.getval(x, y) == 0) {
						if (drand1_() < 1 - exp(2 * interactions.getval(i, j % Ly)*seed.proj*lat.get_spin(x, y))) {
							newspin.x = x;
							newspin.y = y;
							newspin.proj = lat.get_spin(x, y);
							buffer.push_back(newspin);
							cluster.setval(newspin.x, newspin.y, 1);
						}
					}
				}
			}

			//long-range cluster building
			i = 0;
			j = cutoff;
			double pi_Nt = M_PI / Ly;
			while (j < Ly - cutoff) {
				rando = drand1_();
				if (i == 0) {
					//try to place in row 0
					shift = atan(1 / (tan(pi_Nt*(.5*(Ly + 1) - j)) + log(1 - rando) / alpha / pi_Nt));
					k = shift < 0 ? ceil((shift + M_PI)/pi_Nt - .5) : ceil(shift/pi_Nt - .5);
					if (k > Ly - cutoff) {
						//not placed in row 0
						i = 1;
						shift = -log(1 - rando) / alpha / pi_Nt - 0.5*alpha*pi_Nt*(interaction_sum.getval(0, Ly - cutoff) - interaction_sum.getval(0, j - 1) + interaction_sum.getval(1, cutoff - 1));
						//std::cout << "shift after jumping row 0: " << shift << "\n";
						j = cutoff;
						//check for possible solution
						if (1 + shift*shift*(1 - cosh(2 * xc)*cosh(2 * xc)) > 0 && sqrt(0.5 / (1 + shift * shift) * (1 + shift*shift*(1 - cosh(2 * xc)) + sqrt(1 + shift*shift*(1 - cosh(2 * xc)*cosh(2 * xc))))) < 1) {
							k = ceil(0.5*(Ly + 1) + Ly / M_PI * acos(sqrt(0.5 / (1 + shift * shift) * (1 + shift*shift*(1 - cosh(2 * xc)) + sqrt(1 + shift*shift*(1 - cosh(2 * xc)*cosh(2 * xc)))))));
							//std::cout << "k: " << k << "\n";
						}
						//if no solution exit the checker
						else { k = Ly; }
					}
				}
				else if (i == 1) {
					//check for placement into row 1
					shift = -log(1 - rando) / alpha / pi_Nt - 0.5*alpha*pi_Nt*(interaction_sum.getval(1, j));
					//check for possible solution
					if (1 + shift*shift*(1 - cosh(2 * xc)*cosh(2 * xc)) > 0 && sqrt(0.5 / (1 + shift * shift) * (1 + shift*shift*(1 - cosh(2 * xc)) + sqrt(1 + shift*shift*(1 - cosh(2 * xc)*cosh(2 * xc))))) < 1) {
						//std::cout << "shift from row 1: " << shift << "\n";
						k = ceil(0.5*(Ly + 1) + Ly / M_PI * acos(sqrt(0.5 / (1 + shift * shift) * (1 + shift*shift*(1 - cosh(2 * xc)) + sqrt(1 + shift*shift*(1 - cosh(2 * xc)*cosh(2 * xc)))))));
						//if (k < Ly - cutoff) { std::cout << "k: " << k << "\n"; }
					}
					//if no solution exit the checker
					else { k = Ly; }
				}
				y = (seed.y + k + Ly) % Ly;
				x = (seed.x + i) % 2;
				if (k < Ly - cutoff && cluster.getval(x, y) == 0 && interactions.getval(i, k)*seed.proj*lat.get_spin(x, y) < 0) {
					newspin.x = x;
					newspin.y = y;
					newspin.proj = lat.get_spin(x, y);
					buffer.push_back(newspin);
					cluster.setval(newspin.x, newspin.y, 1);
				}
				j = k;
			}
		}
		else { std::cout << "multiple sites not yet implemented\n"; }

	}
	else if (NEAREST_NEIGHBOR_CLUSTER) {
		//Only look at nearest neighbors, and only care about two possible couplings
		double Jx = params.Js[0], Jy = params.Js[1];
		x = seed.x;
		y = seed.y;

		//test x neighbors
		if (seed.proj == lat.get_spin((x + Lx - 1) % Lx, y) && cluster.getval((x + Lx - 1) % Lx, y) == 0 && drand1_() < (1 - exp(-2 * Jx))) {
			newspin.proj = seed.proj;
			newspin.x = (x + Lx - 1) % Lx;
			newspin.y = y;
			cluster.setval((x + Lx - 1) % Lx, y, 1);
			buffer.push_back(newspin);
		}
		if (seed.proj == lat.get_spin((x + 1) % Lx, y) && cluster.getval((x + 1) % Lx, y) == 0 && drand1_() < (1 - exp(-2 * Jx))) {
			newspin.proj = seed.proj;
			newspin.x = (x + 1) % Lx;
			newspin.y = y;
			cluster.setval((x + 1) % Lx, y, 1);
			buffer.push_back(newspin);
		}

		//test y neighbors
		if (seed.proj == lat.get_spin(x, (y + Ly - 1)%Ly) && cluster.getval(x, (y + Ly - 1) % Ly) == 0 && drand1_() < (1 - exp(-2 * Jy))) {
			newspin.proj = seed.proj;
			newspin.x = x;
			newspin.y = (y + Ly - 1) % Ly;
			cluster.setval(x, (y + Ly - 1) % Ly, 1);
			buffer.push_back(newspin);
		}
		if (seed.proj == lat.get_spin(x, (y + 1) % Ly) && cluster.getval(x, (y + 1) % Ly) == 0 && drand1_() < (1 - exp(-2 * Jy))) {
			newspin.proj = seed.proj;
			newspin.x = x;
			newspin.y = (y + 1) % Ly;
			cluster.setval(x, (y + 1) % Ly, 1);
			buffer.push_back(newspin);
		}
	}
	else {

		//do it the naive way - choose a random number then go through
		//lattice until you hit the right probability range and correct spin
		//this will make things easier for adding spin boson model
		//includes the sz-field term h
		std::vector<spin>::iterator newspin_it = not_cluster.begin();
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

/**
----------------------------------------------------------------------------------------------------
General LRW
----------------------------------------------------------------------------------------------------
**/

GeneralLRW::GeneralLRW(class_mc_params params_in) {
	//generate member objects
	for(int i = 0; i < params_in.lengths[0]; ++i){
		interactions.push_back({});
		interaction_sum.push_back({});
		cluster.push_back({});
		for (int j = 0; j < params_in.lengths[1]; ++j){
			interactions[i].push_back(0);
			interaction_sum[i].push_back(0);
			cluster[i].push_back(0);
		}
	}
	cluster_size = 0;
	cluster_mag = 0;

	if (params_in.cutoff_type.compare("exp") == 0){
		std::cout << "Using Finite Exponential Cutoff\n";
		set_spin_boson_model_wc(params_in);
	}
	else if (params_in.cutoff_type.compare("inf") == 0){
		std::cout << "Using Infinite Cutoff\n";
		set_spin_boson_model(params_in);
	}
	else{
		std::cout << "Using Finite Hard Cutoff\n";
		set_spin_boson_model(params_in);
	}
	set_interaction_sum();
	set_site_sums();
	mag = params_in.lengths[0]*params_in.lengths[1] + 1;
	h = params_in.h;
}

void GeneralLRW::set_spin_boson_model(class_mc_params params){
	//STILL NEED TO INCLUDE OMEGA_C AND CHECK DIFFERENT SITE LEADING COEFFICIENT
	//model is translationally invariant, so (i,j) will represent a distance in the (spacial, imaginary time) direction
	double g = params.sbparams.g, A = params.sbparams.A0, gamma = params.Js[1], J = params.Js[0], v = params.sbparams.v, a = params.spacings[0], tc = params.spacings[1], x, y;
	int Nt = params.lengths[1], Lx = params.lengths[0];
	for (int i = 0; i < interactions.size(); ++i) {
		for (int j = 0; j < interactions[i].size(); ++j) {
			if (i == 0) {
				//temporal nearest neighbor interaction
				if (j == 1 || j == Nt - 1) {
					interactions[i][j] = -gamma - 0.5*A*M_PI*M_PI / Nt / Nt / sin(M_PI / Nt) / sin(M_PI / Nt);
				}
				//temporal self-interaction
				if (j > 1 && j < Nt - 1) {
					interactions[i][j] = -0.5*A*M_PI*M_PI / Nt / Nt / sin(M_PI * j / Nt) / sin(M_PI * j / Nt);
				}
			}
			else if (j == 0) {
				//spacial nearest neighbor interactions
				if (i == 1 || i == Lx - 1) {
					interactions[i][j] = -J + A*M_PI*M_PI / Nt / Nt / sinh(M_PI * a / Nt / v / tc) / sinh(M_PI * a / Nt / v / tc);
				}
				//spacial same-time long range interactions
				if (i > 1 && i < Lx - 1) {
					if (i > Lx / 2){
						interactions[i][j] = A*M_PI*M_PI / Nt / Nt / sinh(M_PI * a * (Lx - i) / Nt / v / tc) / sinh(M_PI * a * (Lx - i) / Nt / v / tc);
					}
					else{
						interactions[i][j] = A*M_PI*M_PI / Nt / Nt / sinh(M_PI * a * i / Nt / v / tc) / sinh(M_PI * a * i / Nt / v / tc);
					}
				}
			}
			else {
				//different site, different time interactions
				x = M_PI*j / Nt;
				if (i > Lx / 2){
					y = M_PI*(Lx - i)*a / v / Nt / tc;;
				}
				else{
					y = M_PI*i*a / v / Nt / tc;
				}
				interactions[i][j] = -M_PI*M_PI*A / Nt / Nt*(sin(x)*sin(x)*cosh(y)*cosh(y) - cos(x)*cos(x)*sinh(y)*sinh(y)) /
					((sin(x)*sin(x)*cosh(y)*cosh(y) + cos(x)*cos(x)*sinh(y)*sinh(y))*(sin(x)*sin(x)*cosh(y)*cosh(y) + cos(x)*cos(x)*sinh(y)*sinh(y)));
			}
		}
	}
}

void GeneralLRW::set_spin_boson_model_wc(class_mc_params params){
	//model is translationally invariant, so (i,j) will represent a distance in the (spacial, imaginary time) direction
	double g = params.sbparams.g, A = params.sbparams.A0, gamma = params.Js[1], J = params.Js[0], v = params.sbparams.v, a = params.spacings[0], tc = params.spacings[1], wc = params.sbparams.omega_c, x, y;
	int Nt = params.lengths[1], Lx = params.lengths[0];
	double beta = tc*Nt;
	double omcb = 1.0 / wc / beta;//1/beta/omega_c
	double prefactor = A / Nt / Nt;
	for (int i = 0; i < interactions.size(); ++i) {
		y = (i > Lx / 2) ? a * (Lx - i) / v / beta : y = a * i / v / beta;//periodic boundary condition in spatial direction
		for (int j = 0; j < interactions[i].size(); ++j) {
			x = ((double)j) / Nt;//scaled time coordinate
			if (i == 0) {
				//temporal nearest neighbor interaction
				if (j == 1 || j == Nt - 1) {
					interactions[i][j] = -gamma - prefactor * (real_trigamma_cmplxarg(x + omcb, y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, y, 10, 10));
				}
				//temporal self-interaction
				if (j > 1 && j < Nt - 1) {
					interactions[i][j] = - prefactor * (real_trigamma_cmplxarg(x + omcb, y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, y, 10, 10));
				}
			}
			else if (j == 0) {
				//spacial nearest neighbor interactions
				if (i == 1 || i == Lx - 1) {
					interactions[i][j] = - J - prefactor * (real_trigamma_cmplxarg(x + omcb, y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, y, 10, 10)
																+ real_trigamma_cmplxarg(x + omcb, -y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, -y, 10, 10));
				}
				//spacial same-time long range interactions
				if (i > 1 && i < Lx - 1) {
					interactions[i][j] = - prefactor * (real_trigamma_cmplxarg(x + omcb, y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, y, 10, 10)
																+ real_trigamma_cmplxarg(x + omcb, -y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, -y, 10, 10));
				}
			}
			else {
				//different site, different time interactions
				interactions[i][j] = - prefactor * (real_trigamma_cmplxarg(x + omcb, y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, y, 10, 10)
																+ real_trigamma_cmplxarg(x + omcb, -y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, -y, 10, 10));
			}
		}
	}
}

void GeneralLRW::set_interaction_sum(){
	double total = 0;
	for (int i = 0; i < interaction_sum.size(); ++i){
		for (int j = 0; j < interaction_sum[i].size(); ++j){
			total += interactions[i][j] > 0 ? interactions[i][j] : -interactions[i][j];
			interaction_sum[i][j] = total;
		}
	}
}

void GeneralLRW::set_site_sums(){
	//THIS ASSUMES INTERACTION SUM HAS ALREADY BEEN SET
	//CHECK THIS
	for (int i = 0; i < interaction_sum.size(); ++i){
		site_sums.push_back(interaction_sum[i][interaction_sum[i].size() - 1]);
	}
}

double GeneralLRW::calc_sx(IsingLattice2D& lat) {
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

double GeneralLRW::calc_loc(IsingLattice2D& lat) {
	//calculate the localization - absolute value of sz at each site, averaged
	double final_result = 0.0;
	for (int i = 0; i < lat.get_Lx(); ++i) {
		double site_mag = 0.0;
		for (int j = 0; j < lat.get_Ly(); ++j) {
			site_mag += lat.get_spin(i, j);
		}
		site_mag = site_mag > 0 ? site_mag : -site_mag;
		final_result += site_mag / lat.get_Ly();
	}
	return final_result / lat.get_Lx();
}

double GeneralLRW::calc_s1s2(IsingLattice2D& lat) {
	//calculate the spatial correlation function <S(i)S(i + 1)> (nearest neighbors only)
	double final_result = 0;
	double site_corr;
	int Lx = lat.get_Lx();
	for (int i = 0; i < Lx; ++i) {
		site_corr = 0;
		for (int j = 0; j < lat.get_Ly(); ++j) {
			site_corr += lat.get_spin(i, j)*lat.get_spin((i + 1)%Lx, j);
		}
		final_result += site_corr / lat.get_Ly();
	}
	return final_result / Lx;
}

double GeneralLRW::calc_sz_stagger(IsingLattice2D& lat) {
	//antiferromagnetic analogue of sz().  calc mag for each site, then multiply by (-1)^i and average
	double final_result = 0.0;
	for (int j = 0; j < lat.get_Ly(); ++j) {
		double time_mag = 0.0;
		for (int i = 0; i < lat.get_Lx(); ++i) {
			time_mag += (2*(i/2) == i) ? lat.get_spin(i, j) : -lat.get_spin(i, j);
		}
		final_result += time_mag / lat.get_Lx();
	}
	return (final_result < 0) ? -final_result / lat.get_Ly() : final_result / lat.get_Ly();
}

double GeneralLRW::calc_action_slow(IsingLattice2D& lat) {
	//essentially, the time-averaged energy of a given set of quantum fluctuations
	double S = 0;
	int Lx = lat.get_Lx();
	int Ly = lat.get_Ly();
	int s1;
	for (int i = 0; i < Lx; ++i) {
		for (int j = 0; j < Ly; ++j) {
			s1 = lat.get_spin(i, j);
			for (int m = 0; m < Lx; ++m) {
				for (int n = 0; n < Ly; ++n) {
					S += s1 * lat.get_spin(m, n) * interactions[(m - i + Lx)%Lx][(n - j + Ly)%Ly];
				}
			}
		}
	}

	return 0.5*S;
}

void GeneralLRW::step(IsingLattice2D& lat) {
	//1. determine a "seed" spin
	//2. go to all interacting spins and test to add to buffer and cluster
	//		- use long range procedure in Int. J. Mod. Phys. C 6, 359-370 (1995).
	//3. when buffer is empty, flip cluster
	cluster_size = 0;
	cluster_mag = 0;
	//1. get seed
	spin seed;
    seed.x = (int) (lat.get_Lx() * drand1_());
    seed.y = (int) (lat.get_Ly() * drand1_());
    seed.proj = lat.get_spin(seed.x, seed.y);
    cluster[seed.x][seed.y] = 1;
	++cluster_size;
	cluster_mag += seed.proj;

	//2. add to buffer, then set loop to go through buffer
	test_spins(seed, lat);
	while (buffer.size() != 0) {
		seed = buffer.back();
		cluster_mag += seed.proj;
		++cluster_size;
		buffer.pop_back();
		test_spins(seed, lat);
	}

	//3. flip cluster according to a probability check based on the magnetization change
	//prob is 1 if h*mag is negative, exp(-2*h*mag) if h*mag is positive
	//satisfies detailed balance: this is the "acceptance probability", rest was "selection probability"
	if (drand1_() < exp(-2 * h * cluster_mag)) {
		for (int x = 0; x < lat.get_Lx(); ++x) {
			for (int y = 0; y < lat.get_Ly(); ++y) {
				if (cluster[x][y] == 1) {
					lat.flip_spin(x, y);
					cluster[x][y] = 0;
				}
			}
		}
	}
	mag -= (2.0*cluster_mag)/lat.get_N();

}

bool GeneralLRW::step_one_site(IsingLattice2D& lat, double * prev_action, 
		thrust::host_vector<double>& thrustlat_ref, thrust::host_vector<double>& thrustint_ref, thrust::host_vector<double>& corr_ref){
    //Build a cluster in only the time direction of one spatial site.  Then attempt to flip with Metropolis probability according to previous action
    //and newly calculated action.  Require GPU.  Mostly copied from step() function for full lattice clusters
    //return true if flipped, false if not
	cluster_size = 0;
	cluster_mag = 0;
	//1. get seed
	spin seed;
    seed.x = (int) (lat.get_Lx() * drand1_());
    seed.y = (int) (lat.get_Ly() * drand1_());
    seed.proj = lat.get_spin(seed.x, seed.y);
    cluster[seed.x][seed.y] = 1;
	++cluster_size;
	cluster_mag += seed.proj;

	//2. add to buffer, then set loop to go through buffer
	test_spins_one_site(seed, lat);
	while (buffer.size() != 0) {
		seed = buffer.back();
		cluster_mag += seed.proj;
		++cluster_size;
		buffer.pop_back();
		test_spins_one_site(seed, lat);
	}

	//3. flip cluster according to a probability check based on the magnetization change and interaction change
	//prob is 1 if( new_action - prev_action + 2*h*mag) is negative, exp(prev_action - new_action - 2*h*mag) if h*mag is positive
	//satisfies detailed balance: this is the "acceptance probability", rest was "selection probability"
    //
    //need to subtract off the action contribution from the seed site, because this portion of
    //detailed balance is already satisfied by the selection algorithm.
    double new_action = 0, new_t_cont = 0, prev_t_cont = 0;

    thrust::host_vector<double> thrustsite(lat.get_Ly());
    thrust::host_vector<double>& thrustsite_ref = thrustsite;
    thrust::host_vector<double> thrustint_site(lat.get_Ly());
    thrust::host_vector<double>& thrustint_site_ref = thrustint_site;
	thrust::host_vector<double> site_corr_trash(lat.get_Ly());
	thrust::host_vector<double>& site_corr_trash_ref = site_corr_trash;
	thrust::host_vector<double> full_corr_trash(lat.get_Lx()*lat.get_Ly());
	thrust::host_vector<double>& full_corr_trash_ref = full_corr_trash;
    for(int y = 0; y < lat.get_Ly(); ++y){
        thrustsite[y] = lat.get_spin(seed.x, y);
        thrustint_site[y] = interactions[0][y];
    }
	prev_t_cont = cufft_calc_action(thrustsite_ref, thrustint_site_ref, site_corr_trash, 1, lat.get_Ly());
	//flip cluster to get proposed state
	for (int y = 0; y < lat.get_Ly(); ++y) {
		if (cluster[seed.x][y] == 1) {
			thrustlat_ref[seed.x*lat.get_Ly() + y] *= -1;
			thrustsite_ref[y] *= -1;
		}
	}
    new_action = cufft_calc_action(thrustlat_ref, thrustint_ref, full_corr_trash_ref, lat.get_Lx(), lat.get_Ly());
    new_t_cont = cufft_calc_action(thrustsite_ref, thrustint_site_ref, site_corr_trash, 1, lat.get_Ly());
                                
	if (drand1_() < exp((*prev_action - prev_t_cont) - (new_action - new_t_cont) - 2 * h * cluster_mag)) {
		for (int x = 0; x < lat.get_Lx(); ++x) {
			for (int y = 0; y < lat.get_Ly(); ++y) {
				if (cluster[x][y] == 1) {
					lat.flip_spin(x, y);
					cluster[x][y] = 0;
				}
				//copy new correlation function
				corr_ref[x*lat.get_Ly() + y] = full_corr_trash[x*lat.get_Ly() + y];
			}
		}
	    mag -= (2.0*cluster_mag)/lat.get_N();
        *prev_action = new_action;
        return true;
	}
    else{

		for (int x = 0; x < lat.get_Lx(); ++x) {
			for (int y = 0; y < lat.get_Ly(); ++y) {
                cluster[x][y] = 0;
			}
		}
		//return thrustlat_ref to its original state
		lat.get_thrust_vector(thrustlat_ref);
        return false;
    }
}

void GeneralLRW::test_spins(spin seed, IsingLattice2D& lat){
	int i,j, k, l; //i and j is searched site distance, k and l are the corresponding lattice indices
	double randnum, totalint = 0; //random number, and running interaction sum
	std::vector<double>::iterator x_it_temp = site_sums.begin(), 
			x_it = site_sums.begin(), x_it_end = site_sums.end(), y_it, y_it_end; //x_it is 0 site separation, x_it_end is the end of the site vector
			//y_it will be the beginning of the interaction sum determined by the site search, and y_it_end will be the corresponding ending iterator
	spin add;
	//
	//	Outline of algorithm:
	//		1. get random number.  Find corresponding interaction sum by sum = -1/2 log(rand)
	//		2. search for totalint + randsum; first site_sums, then full interaction sum at the given site
	//		3. add corresponding spin if not in cluster and interaction with seed is favorable
	//		4. update totalint
	//
	//		5. end when totalint + randsum exceeds all interaction sums
	//
	while(x_it_temp != x_it_end){
		randnum = drand1_();
		totalint = totalint - 0.5*log(randnum);
		x_it_temp = std::lower_bound(x_it, x_it_end, totalint);
		if(x_it_temp != x_it_end){
			i = x_it_temp - x_it;
			y_it = interaction_sum[i].begin();
			y_it_end = interaction_sum[i].end();
			j = (std::lower_bound(y_it, y_it_end, totalint) - y_it);
			k = (seed.x + i) % lat.get_Lx();
			l = (seed.y + j) % lat.get_Ly();
			if (j < lat.get_Ly() && cluster[k][l] == 0 && interactions[i][j]*seed.proj*lat.get_spin(k, l) < 0){
				add.x = k;
				add.y = l;
				add.proj = lat.get_spin(k, l);
				buffer.push_back(add);
				cluster[k][l] = 1;
			}
		}
	}
}


void GeneralLRW::test_spins_one_site(spin seed, IsingLattice2D& lat){
    //same as test_spins(), but cluster is restricted to the original x position of the seed
	int i,j, k, l; //i and j is searched site distance, k and l are the corresponding lattice indices
	double randnum, totalint = 0; //random number, and running interaction sum
	std::vector<double>::iterator y_it, y_it_temp, y_it_end;//y_it will be the beginning of the interaction sum determined by the site search, and y_it_end will be the corresponding ending iterator
	spin add;
	//
	//	Outline of algorithm:
	//		1. get random number.  Find corresponding interaction sum by sum = -1/2 log(rand)
	//		2. search for totalint + randsum; first site_sums, then full interaction sum at the given site
	//		3. add corresponding spin if not in cluster and interaction with seed is favorable
	//		4. update totalint
	//
	//		5. end when totalint + randsum exceeds all interaction sums
	//

    i = seed.x;
    y_it = interaction_sum[0].begin();
	y_it_temp = y_it;
    y_it_end = interaction_sum[0].end();
    while(y_it_temp != y_it_end){
		randnum = drand1_();
		totalint = totalint - 0.5*log(randnum);
		y_it_temp = std::lower_bound(y_it, y_it_end, totalint);
		if(y_it_temp != y_it_end){
			j = (y_it_temp - y_it);
			l = (seed.y + j) % lat.get_Ly();
			if (j < lat.get_Ly() && cluster[i][l] == 0 && interactions[0][j]*seed.proj*lat.get_spin(i, l) < 0){
				add.x = i;
				add.y = l;
				add.proj = lat.get_spin(i, l);
				buffer.push_back(add);
				cluster[i][l] = 1;
			}
		}
    }
}
