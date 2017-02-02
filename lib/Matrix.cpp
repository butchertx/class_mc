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

void Matrix::setval(int row, int col, double val) {
	vals[row][col] = val;
}

void Matrix::fill(double val) {
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			vals[i][j] = val;
		}
	}
}

double Matrix::getval(int row, int col) {
	return vals[row][col];
}

double Matrix::getval(int n) {
	//if 0 <= n < dimx * dimy, return val at (i, j) counting across rows
	//n = i*cols + j: n/cols = i + j/cols
	if (n < 0 || n >= rows*cols) {
		std::cout << "Invalid choice of n in Matrix recall\n";
		return 0.0;
	}
	else {
		return vals[n / cols][n % cols];
	}
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

void IntMatrix::setval(int row, int col, int val) {
	vals[row][col] = val;
}

void IntMatrix::fill(int val) {
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			vals[i][j] = val;
		}
	}
}

int IntMatrix::getval(int row, int col) {
	return vals[row][col];
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
