#pragma once
#include <fstream>
#include <vector>
#include <thrust/host_vector.h>
#include "cufft.h"

typedef float2 Complex;


template <typename InputIterator1,
	typename InputIterator2,
	typename OutputIterator>
	OutputIterator expand(InputIterator1 first1,
		InputIterator1 last1,
		InputIterator2 first2,
		OutputIterator output);

void calc_corr_fast_1site(thrust::host_vector<double>&, thrust::host_vector<int>&, int);

void calc_corr_fast_2site(thrust::host_vector<double>&, thrust::host_vector<int>&, int);

double calc_action_fast(thrust::host_vector<double>& corr, thrust::host_vector<double>& interactions);

__global__ void calc_shift_state(double* state, double* shift_state, int i, int j);

double thrust_calc_action_general(thrust::host_vector<double>& state, thrust::host_vector<double>& interactions, int Lx, int Ly, dim3 threads);

__global__ void elementwise_product_cmplx(cufftComplex *source);

double cufft_calc_action(thrust::host_vector<double>& state, thrust::host_vector<double>& interactions, thrust::host_vector<double>& corr, int Lx, int Ly);

void show_memory();

