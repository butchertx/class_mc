#pragma once
#include <vector>
#include <sstream>
#include <cmath>
#include <thrust/host_vector.h>
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

	void set_ferro(){
		for (int i = 0; i < Lx; ++i){
			for (int j = 0; j < Ly; ++j){
				spins[i][j] = 1;
			}
		}
	}

	void set_antiferro(){
		for (int i = 0; i < Lx; ++i){
			for (int j = 0; j < Ly; ++j){
				spins[i][j] = pow(-1, i + j);
			}
		}
	}

	//returns the value that the spin changes to
	int flip_spin_return(int x, int y) {
		spins[x][y] *= -1;
		return spins[x][y];
	}

	void flip_spin(int x, int y) {
		spins[x][y] *= -1;
	}

	int get_spin(int x, int y) {
		return spins[(x + Lx)%Lx][(y + Ly)%Ly];
	}

	std::vector<int> get_bool_spins(){
		std::vector<int> bool_spins(Lx*Ly);
		for(int i = 0; i < Lx; ++i){
			for(int j = 0; j < Ly; ++j){
				bool_spins[i*Ly + j] = spins[i][j];
			}
		}
		return bool_spins;
	}

	void copy_bool_spins(std::vector<int> bool_spins){
		for(int i = 0; i < Lx; ++i){
			for (int j = 0; j < Ly; ++j){
				if(bool_spins[i*Ly + j] == 1){
					spins[i][j] = 1;
				}
				else{
					spins[i][j] = -1;
				}
			}
		}
	}
	
    void get_thrust_vector(thrust::host_vector<double>& result){
		if (result.size() != Lx*Ly){
			result.resize(Lx*Ly);
		}
		for (int i = 0; i < Lx; ++i){
			for (int j = 0; j < Ly; ++j){
				result[i*Ly + j] = (double) spins[i][j];
			}
		}
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

	std::string to_bitmap(){
		std::stringstream outstring1;
		for (int i = 0; i < Ly; ++i){
			for (int j = 0; j < Lx; ++j){
				outstring1 << (int) (0.5*(spins[j][i] + 1)) << " ";
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
