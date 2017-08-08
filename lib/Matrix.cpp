#include "Matrix.h"
#include <sstream>

Matrix::Matrix() {
	vals = std::vector<std::vector<double>>(10,std::vector<double>(10, 0.0));
	rows = 10;
	cols = 10;
}

Matrix::Matrix(int M, int N) {
	//MxN matrix: M rows and N columns
	vals = std::vector<std::vector<double>>(M, std::vector<double>(N, 0.0));
	rows = M;
	cols = N;
}

std::string Matrix::to_string() {
	std::stringstream outstring;
	for (int i = 0; i < vals.size(); ++i) {
		for (int j = 0; j < vals[i].size(); ++j) {
			outstring << vals[i][j] << " ";
		}
		outstring << "\n";
	}
	return outstring.str();
}

IntMatrix::IntMatrix() {
	vals = std::vector<std::vector<int>>(10,std::vector<int>(10, 0));
	rows = 10;
	cols = 10;
}

IntMatrix::IntMatrix(int M, int N) {
	//MxN matrix: M rows and N columns
	vals = std::vector<std::vector<int>>(M, std::vector<int>(N, 0));
	rows = M;
	cols = N;
}



std::string IntMatrix::to_string() {
	std::stringstream outstring;
	for (int i = 0; i < vals.size(); ++i) {
		for (int j = 0; j < vals[i].size(); ++j) {
			outstring << vals[i][j] << " ";
		}
		outstring << "\n";
	}
	return outstring.str();
}
