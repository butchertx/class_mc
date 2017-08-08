#include "obs_calc_fast.cuh"
#include <sstream>
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include <thrust/copy.h>
#include <thrust/reduce.h>
#include <thrust/gather.h>
#include <thrust/scan.h>
#include <thrust/fill.h>
#include <thrust/copy.h>

#include <iterator>
#include <iostream>

template <typename InputIterator1,
	typename InputIterator2,
	typename OutputIterator>
	OutputIterator expand(InputIterator1 first1,
		InputIterator1 last1,
		InputIterator2 first2,
		OutputIterator output)
{
	typedef typename thrust::iterator_difference<InputIterator1>::type difference_type;

	difference_type input_size = thrust::distance(first1, last1);
	difference_type output_size = thrust::reduce(first1, last1);

	// scan the counts to obtain output offsets for each input element
	thrust::device_vector<difference_type> output_offsets(input_size, 0);
	thrust::exclusive_scan(first1, last1, output_offsets.begin());

	// scatter the nonzero counts into their corresponding output positions
	thrust::device_vector<difference_type> output_indices(output_size, 0);
	thrust::scatter_if
	(thrust::counting_iterator<difference_type>(0),
		thrust::counting_iterator<difference_type>(input_size),
		output_offsets.begin(),
		first1,
		output_indices.begin());

	// compute max-scan over the output indices, filling in the holes
	thrust::inclusive_scan
	(output_indices.begin(),
		output_indices.end(),
		output_indices.begin(),
		thrust::maximum<difference_type>());

	// gather input values according to index array (output = first2[output_indices])
	OutputIterator output_end = output; thrust::advance(output_end, output_size);
	thrust::gather(output_indices.begin(),
		output_indices.end(),
		first2,
		output);

	// return output + output_size
	thrust::advance(output, output_size);
	return output;
}


void calc_corr_fast_1site(thrust::host_vector<double>& corr, thrust::host_vector<int>& state, int Ly) {

	//begin calculation
	//copy to device
	thrust::device_vector<int> d_state = state;

	//calculate correlation
	//only works for Lx = 1 right now
	for (int i = 0; i < Ly; ++i) {
		corr[i] = 1.0 / (Ly) * (thrust::inner_product(d_state.begin(), d_state.end() - i, d_state.begin() + i, 0) + thrust::inner_product(d_state.begin(), d_state.begin() + i, d_state.end() - i, 0));
	}

}

void calc_corr_fast_2site(thrust::host_vector<double>& corr, thrust::host_vector<int>& state, int Ly) {
	thrust::device_vector<int> d_state = state;

	for (int i = 0; i < Ly; ++i) {
		corr[i] = 1.0 / (2 * Ly) * (thrust::inner_product(d_state.begin(), d_state.begin() + Ly - i, d_state.begin() + i, 0)
						+ thrust::inner_product(d_state.begin(), d_state.begin() + i, d_state.begin() + Ly - i, 0)
						+ thrust::inner_product(d_state.begin() + Ly, d_state.end() - i, d_state.begin() + Ly + i, 0)
						+ thrust::inner_product(d_state.begin() + Ly, d_state.begin() + Ly + i, d_state.end() - i, 0));
		corr[i + Ly] = 1.0 / (2 * Ly) * (thrust::inner_product(d_state.begin(), d_state.begin() + Ly - i, d_state.begin() + Ly + i, 0) 
						+ thrust::inner_product(d_state.begin(), d_state.begin() + i, d_state.end() - i, 0)
						+ thrust::inner_product(d_state.begin() + Ly, d_state.end() - i, d_state.begin() + i, 0) 
						+ thrust::inner_product(d_state.begin() + Ly, d_state.begin() + Ly + i, d_state.begin() + Ly - i, 0));
	}
}

double calc_action_fast(thrust::host_vector<double>& corr, thrust::host_vector<double>& interactions) {
	//interactions matrix at (i) should give interaction between sites x1 - x2 = i and interactions(0) = 0
	//this is actually the inner product of the correlation function with the interactions matrix
	return 0.5*thrust::inner_product(corr.begin(), corr.end(), interactions.begin(), 0.0);
}

__global__ void calc_shift_state(double* state, double* shift_state, int i, int j){
	//i and j are the "starting point" (0,0) of the new lattice state, held in shift_state
	//the lattice is shifted using pbc's to start at the new value of i,j
	//assume it is launched with gridDim.x = Lx and blockDim.x = Ly/32, blockDim.y = 32
	//position in shift_state is determined by thread variables, position in original state is relative to i and j
	//int Lx = gridDim.x, Ly = blockDim.x;
	//int shift_ind = blockIdx.x*Ly + threadIdx.x;//position in shift_state
	//int state_ind_x = (i + blockIdx.x) % Lx;//x coordinate of (i,j) + (m,n)
	//int state_ind_y = (j + threadIdx.x) % Ly;//y coordinate of (i,j) + (m,n)
	//int state_ind = state_ind_x*Ly + state_ind_y;
	//shift_state[shift_ind] = state[state_ind];
	//the above works when grids and blocks are 1D

	int Lx = gridDim.x, Ly = blockDim.x*blockDim.y*blockDim.z;
	int shift_ind = blockIdx.x*Ly + threadIdx.z*blockDim.y*blockDim.x + threadIdx.y*blockDim.x + threadIdx.x;//position in shift_state
	int state_ind_x = (i + blockIdx.x) % Lx;//x coordinate of (i,j) + (m,n)
	int state_ind_y = (j + threadIdx.z*blockDim.y*blockDim.x + threadIdx.y*blockDim.x + threadIdx.x) % Ly;//y coordinate of (i,j) + (m,n)
	int state_ind = state_ind_x*Ly + state_ind_y;
	shift_state[shift_ind] = state[state_ind];
}

double thrust_calc_action_general(thrust::host_vector<double>& state, thrust::host_vector<double>& interactions, int Lx, int Ly, dim3 threads){
	//
	//	Algorithm outline
	//	Def: L = state, A = interactions, Q_ij = sum_(m,n) L_(i + m, j + n) * A_(m,n), action = S = sum_(i,j) L_ij Q_ij
	//	2 steps: calculate Q_ij, then do L.Q (dot product)
	//	1. Q_ij calculation
	//
	if(state.size() != interactions.size() || state.size() != Lx*Ly){
		std::cout << "Error: state/interactions size mismatch\n";
		return 0;
	}
	thrust::device_vector<double> qij(Lx*Ly);
	thrust::device_vector<double> d_state = state;
	thrust::device_vector<double> shift_state = state;
	thrust::device_vector<double> d_int = interactions;
	for(int i = 0; i < Lx; ++i){
		for (int j = 0; j < Ly; ++j){
			calc_shift_state<<<Lx , threads>>>(thrust::raw_pointer_cast(&d_state[0]), thrust::raw_pointer_cast(&shift_state[0]), i, j);
			qij[i*Ly + j] = thrust::inner_product(shift_state.begin(), shift_state.end(), d_int.begin(), 0.0);
		}
	}
	return 0.5*thrust::inner_product(d_state.begin(), d_state.end(), qij.begin(), 0.0);
}

__global__ void elementwise_product_cmplx(cufftDoubleComplex *source) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	source[i].x = source[i].x * source[i].x + source[i].y *source[i].y;
	source[i].y = 0;
}

double cufft_calc_action(thrust::host_vector<double>& state, thrust::host_vector<double>& interactions, thrust::host_vector<double>& corr, int Lx, int Ly) {
	//	Algorithm: Calculate the correlation function and then take the inner product with the interactions
	//	Correlation Calculation: correlation function is the inverse FT of |S_pr . S_pr|, calculated elementwise, where S_pr is the FT of the state vector

	cufftHandle forward_plan, backward_plan;
	cudaError cuda_status;

	cufftDoubleReal *state_rs;//real space state

	cufftDoubleComplex *state_ft;//fourier space state

	int n[2] = { Lx, Ly };

	//allocate real state
	cuda_status = cudaMalloc((void**)&state_rs, sizeof(cufftDoubleReal)*Lx*Ly);
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "Cuda error: Failed to allocate. Code: %s\n", cudaGetErrorName(cuda_status));
		show_memory();
		return 0;
	}

	//allocate transform state
	cuda_status = cudaMalloc((void**)&state_ft, sizeof(cufftDoubleComplex)*Lx*(Ly/2 + 1));
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "Cuda error: Failed to allocate. Code: %s\n", cudaGetErrorName(cuda_status));
		show_memory();
		return 0;
	}

	//copy vals for real state
	if (cudaMemcpy(state_rs, thrust::raw_pointer_cast(&state[0]), sizeof(cufftDoubleReal)*Lx*Ly, cudaMemcpyHostToDevice) != cudaSuccess) {
		fprintf(stderr, "Cuda error: Failed to copy\n");
		return 0;
	}

	//Create 2D R2C FFT plan
	if (cufftPlanMany(&forward_plan, 2, n, NULL, 1, 0, NULL, 1, 0, CUFFT_D2Z, 1) != CUFFT_SUCCESS) {
		fprintf(stderr, "CUFFT Error: Unable to create plan\n");
		return 0;
	}

	//Transform state
	if (cufftExecD2Z(forward_plan, state_rs, state_ft) != CUFFT_SUCCESS) {
		fprintf(stderr, "CUFFT Error: Transform did not work\n");
		return 0;
	}

	//Multiply state_ft with itself.  If Ly/2 >= 1024, this needs to be modified
	if (Ly / 2 >= 1024) { std::cout << "Error: need to modify block structure to make correlation calculation correct\n"; return 0; }
	elementwise_product_cmplx<<<Lx, (Ly/2 + 1)>>>(state_ft);

	//Inverse fourier transform plan
	if (cufftPlanMany(&backward_plan, 2, n, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2D, 1) != CUFFT_SUCCESS) {
		fprintf(stderr, "CUFFT Error: Unable to create plan\n");
		return 0;
	}

	//Find correlation by taking the inverse fourier transform of state_ft
	if (cufftExecZ2D(backward_plan, state_ft, state_rs) != CUFFT_SUCCESS) {
		fprintf(stderr, "CUFFT Error: Transform did not work\n");
		return 0;
	}

	//copy correlation back to a host vector
	if (cudaMemcpy(thrust::raw_pointer_cast(&corr[0]), state_rs, sizeof(cufftDoubleReal)*Lx*Ly, cudaMemcpyDeviceToHost) != cudaSuccess) {
		fprintf(stderr, "Cuda error: Failed to copy\n");
		return 0;
	}

	//Deallocate memory
	cufftDestroy(forward_plan);
	cufftDestroy(backward_plan);
	cudaFree(state_rs);
	cudaFree(state_ft);

	return 0.5*thrust::inner_product(interactions.begin(), interactions.end(), corr.begin(), 0.0) / ((double)Lx*Ly);
}

void show_memory() {

	// show memory usage of GPU

	size_t free_byte;

	size_t total_byte;

	cudaError cuda_status;

	cuda_status = cudaMemGetInfo(&free_byte, &total_byte);

	if (cudaSuccess != cuda_status) {

		printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status));

		exit(1);

	}



	double free_db = (double)free_byte;

	double total_db = (double)total_byte;

	double used_db = total_db - free_db;

	printf("GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n",

		used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0, total_db / 1024.0 / 1024.0);
}