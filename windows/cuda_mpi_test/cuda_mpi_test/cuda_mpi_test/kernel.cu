
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
#include <thrust/host_vector.h>
#include "MemTimeTester.h"
#include "obs_calc_fast.cuh"
extern "C" {
#include "random.h"
}
/**
For my GPU, the maximum threads per multiprocessor is 2048, max per block is 1024, and the max dimensions are (1024, 1024, 64)
Maximum shared memory per block is 49152 bytes, L2 Cache is 2097152 bytes, 15 multiprocessors, 128 CUDA cores/MP
Cuda capability: 6.1, CUDA driver version 8.0
**/

double calc_action_slow(thrust::host_vector<double>& lat, thrust::host_vector<double>& interactions, int Lx, int Ly) {
	//essentially, the time-averaged energy of a given set of quantum fluctuations.  Calculate the action then divide by beta
	double S = 0;
	double s1;
	for (int i = 0; i < Lx; ++i) {
		for (int j = 0; j < Ly; ++j) {
			s1 = lat[i*Ly + j];
			for (int m = 0; m < Lx; ++m) {
				for (int n = 0; n < Ly; ++n) {
					S += s1 * lat[m*Ly + n] * interactions[((m - i + Lx) % Lx)*Ly + ((n - j + Ly) % Ly)];
				}
			}
		}
	}

	return 0.5*S;
}

thrust::host_vector<double> rand_vector(int length) {
	thrust::host_vector<double> result(length);
	for (int i = 0; i < length; ++i) {
		result[i] = drand1_();
	}
	return result;
}

thrust::host_vector<double> transpose(thrust::host_vector<double> s, int Lx, int Ly) {
	//take a vector indexed like s(x, y) = s[x*Ly + y] and make it s(x,y) = s[y*Lx + x]
	thrust::host_vector<double> new_s(Lx*Ly);
	for (int x = 0; x < Lx; ++x) {
		for (int y = 0; y < Ly; ++y) {
			new_s[y*Lx + x] = s[x*Ly + y];
		}
	}
	return new_s;
}

int main() {
	std::cout << "Timing different partitions of fast calc kernel\n";
	MemTimeTester timer;
	int seed = 1892347;
	rand_init_(&seed);
	double fast_action;
	double cufft_action;
	double slow_action;


	//create different sizes of lattice to calculate, then populate with random numbers
	//sizes to test
	//2x128, 2x1024, 2x16384
	//32x32, 32x64, 32x1024
	//1024x1024

	thrust::host_vector<double> vec2x128 = rand_vector(256),
		int2x128 = rand_vector(256),
		vec32x32 = rand_vector(32 * 32),
		int32x32 = rand_vector(32 * 32),
		vec32x64 = rand_vector(32 * 64),
		int32x64 = rand_vector(32 * 64),
		vec32x1024 = rand_vector(32 * 1024),
		int32x1024 = rand_vector(32 * 1024);

	//2x128
	dim3 threads(128, 1, 1);
	thrust::host_vector<double>& sref = vec2x128;
	thrust::host_vector<double>& iref = int2x128;
	timer.flag_start_time("2x128 {2,0,0} {128,0,0}");
	fast_action = thrust_calc_action_general(sref, iref, 2, 128, threads);
	timer.flag_end_time("2x128 {2,0,0} {128,0,0}");
	slow_action = calc_action_slow(sref, iref, 2, 128);
	cufft_action = cufft_calc_action(sref, iref, 2, 128) / 256.0;
	std::cout << "Fast action = " << fast_action << ", slow action = " << slow_action << ", cufft action: " << cufft_action << "\n";
	show_memory();

	//32x32
	dim3 threads2(4, 4, 2);
	sref = vec32x32;
	iref = int32x32;
	timer.flag_start_time("32x32 {32,0,0} {32,0,0}");
	fast_action = thrust_calc_action_general(sref, iref, 32, 32, threads2);
	timer.flag_end_time("32x32 {32,0,0} {32,0,0}");
	slow_action = calc_action_slow(sref, iref, 32, 32);
	cufft_action = cufft_calc_action(sref, iref, 32, 32) / (32.0*32);
	std::cout << "Fast action = " << fast_action << ", slow action = " << slow_action << ", cufft action: " << cufft_action << "\n";
	show_memory();


	//32x64
	dim3 threads4(4, 4, 4);
	sref = vec32x64;
	iref = int32x64;
	timer.flag_start_time("32x64 {32,0,0} {32,0,0}");
	fast_action = thrust_calc_action_general(sref, iref, 32, 64, threads4);
	timer.flag_end_time("32x64 {32,0,0} {32,0,0}");
	slow_action = calc_action_slow(sref, iref, 32, 64);
	cufft_action = cufft_calc_action(sref, iref, 32, 64) / (32.0 * 64);
	std::cout << "Fast action = " << fast_action << ", slow action = " << slow_action << ", cufft action: " << cufft_action << "\n";
	show_memory();

	//32x1024
	dim3 threads3(32, 32, 1);
	sref = vec32x1024;
	iref = int32x1024;
	timer.flag_start_time("big");
	fast_action = thrust_calc_action_general(sref, iref, 32, 1024, threads3);
	timer.flag_end_time("big");
	timer.flag_start_time("CUFFT");
	cufft_action = cufft_calc_action(sref, iref, 32, 1024)/(32.0*1024);
	timer.flag_end_time("CUFFT");
	std::cout << "Fast action = " << fast_action << ", cufft action: " << cufft_action << "\n";
	show_memory();

	timer.print_timers();

	return 0;
}