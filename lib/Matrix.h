#pragma once
#include <vector>
#include <string>

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

	void setval(int, int, double);

	void fill(double);

	double getval(int, int);

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

	void setval(int, int, int);
	void fill(int);
	int getval(int, int);
	std::string to_string();
};
