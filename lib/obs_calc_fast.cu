#include "obs_calc_fast.cuh"
#include <sstream>
#include <thrust/host_vector.h>
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


std::vector<int> calc_corr_fast(std::ifstream* file_p) {
	int Lx, Ly;
	thrust::host_vector<int> state;

	if (file_p->is_open()) {
		*file_p >> Lx >> Ly;
		state.resize(Lx * Ly);
		for (int i = 0; i < Lx; ++i) {
			for (int j = 0; j < Ly; ++j) {
				*file_p >> state[i*Ly + j];
			}
		}
	}
	else {
		std::cout << "Error: input file not opened\n";
	}
	////output state vector to check
	//for (int i = 0; i < Lx; ++i) {
	//	for (int j = 0; j < Ly; ++j) {
	//		std::cout << state[i*Ly + j] << ", ";
	//	}
	//	std::cout << "\n";
	//}

	//begin calculation
	//copy to device
	thrust::device_vector<int> d_state = state;
	thrust::device_vector<int> corr(Lx*Ly, 0);


	////output state vector to check
	//for (int i = 0; i < Lx; ++i) {
	//	for (int j = 0; j < Ly; ++j) {
	//		std::cout << d_state[i*Ly + j] << ", ";
	//	}
	//	std::cout << "\n";
	//}

	//calculate correlation
	//only works for Lx = 1 right now
	for (int i = 0; i < Ly; ++i) {
		corr[i] = thrust::inner_product(d_state.begin(), d_state.end() - i, d_state.begin() + i, 0) + thrust::inner_product(d_state.begin(), d_state.begin() + i, d_state.end() - i, 0);
	}

	//for (int i = 0; i < Lx; ++i) {
	//	for (int j = 0; j < Ly; ++j) {
	//		std::cout << corr[i*Ly + j] << ", ";
	//	}
	//	std::cout << "\n";
	//}

	return{ 0 };
}