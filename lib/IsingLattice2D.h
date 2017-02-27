#pragma once
#include <vector>
#include <sstream>
#include "MemTimeTester.h"
extern "C" {
#include "random.h"
}

class IsingLattice2D {
	int Lx, Ly;
	std::vector<std::vector<int>> spins;

public:
	IsingLattice2D(){}
	IsingLattice2D(int lx, int ly, int randseed) {
		Lx = lx;
		Ly = ly;
		//randomly set spins
		rand_init_(&randseed);
		for (int i = 0; i < lx; ++i) {
			spins.push_back({});
			for (int j = 0; j < ly; ++j) {
				spins[i].push_back(drand1_() > .5 ? 1 : -1);
			}
		}
	}

	//returns the value that the spin changes to
	void flip_spin(int x, int y) {
		spins[x][y] *= -1;
	}

	int get_spin(int x, int y) {
		return spins[x][y];
	}

	std::vector<int>::iterator get_col_it(int x) {
		return spins[x].begin();
	}

	std::vector<std::vector<int>>::iterator get_lat_it_begin() {
		return spins.begin();
	}

	std::vector<std::vector<int>>::iterator get_lat_it_end() {
		return spins.end();
	}

	int get_Lx() { return Lx; }
	int get_Ly() { return Ly; }
	int get_N() { return Lx * Ly; }

	void print_lattice() {
		std::stringstream outstring1;
		for (std::vector<std::vector<int>>::iterator itx = spins.begin(); itx != spins.end(); ++itx) {
			for (std::vector<int>::iterator ity = (*itx).begin(); ity != (*itx).end(); ++ity) {
				outstring1 << *ity << " ";
			}
			outstring1 << "\n";
		}
		std::cout << outstring1.str();

	}

	std::string to_string() {
		std::stringstream outstring1;
		for (std::vector<std::vector<int>>::iterator itx = spins.begin(); itx != spins.end(); ++itx) {
			for (std::vector<int>::iterator ity = (*itx).begin(); ity != (*itx).end(); ++ity) {
				outstring1 << *ity << " ";
			}
			outstring1 << "\n";
		}
		return outstring1.str();
	}
};

struct spin {
	//holds the value (+-1) and the position indices of a given spin
	int proj;
	int x;
	int y;
};