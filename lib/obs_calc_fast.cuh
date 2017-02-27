#pragma once
#include <fstream>
#include <vector>

template <typename InputIterator1,
	typename InputIterator2,
	typename OutputIterator>
	OutputIterator expand(InputIterator1 first1,
		InputIterator1 last1,
		InputIterator2 first2,
		OutputIterator output);

std::vector<int> calc_corr_fast(std::ifstream* file_p);
