#pragma once
#include <vector>
#include <string>
#include <iostream>

class Matrix
{
	//uses row-first ordering: (x, y) refers to row x, col y
	//store as a vector of row vectors
	std::vector<std::vector<double>> vals;
	int rows, cols;
public:
	Matrix();
	Matrix(int, int);

	int get_dimx() { return rows; }
	int get_dimy() { return cols; }

	void setval(int row, int col, double val) {
		vals[row][col] = val;
	}

	void fill(double val) {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				vals[i][j] = val;
			}
		}
	}

	double getval(int row, int col) {
		return vals[row][col];
	}

	double getval(int n) {
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

	std::string to_string();
};

class IntMatrix
{
	//just like a matrix but with ints instead of doubles
	std::vector<std::vector<int>> vals;
	int rows, cols;
public:
	IntMatrix();
	IntMatrix(int, int);

	int get_dimx() { return rows; }
	int get_dimy() { return cols; }

	void setval(int row, int col, int val) {
		vals[row][col] = val;
	}

	void fill(int val) {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				vals[i][j] = val;
			}
		}
	}

	int getval(int row, int col) {
		return vals[row][col];
	}

	std::string to_string();
};
