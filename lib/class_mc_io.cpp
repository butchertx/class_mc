#include "class_mc_io.h"
#include <cmath>
#include <iostream>
#include <string>
#include <sys/stat.h> // stat
#include <errno.h>    // errno, ENOENT, EEXIST
#if defined(_WIN32)
#include <direct.h>   // _mkdir
#endif

bool isDirExist(const std::string& path)
{
#if defined(_WIN32)
	struct _stat info;
	if (_stat(path.c_str(), &info) != 0)
	{
		return false;
	}
	return (info.st_mode & _S_IFDIR) != 0;
#else 
	struct stat info;
	if (stat(path.c_str(), &info) != 0)
	{
		return false;
	}
	return (info.st_mode & S_IFDIR) != 0;
#endif
}

bool makePath(const std::string& path)
{
#if defined(_WIN32)
	int ret = _mkdir(path.c_str());
#else
	mode_t mode = 0755;
	int ret = mkdir(path.c_str(), mode);
#endif
	if (ret == 0)
		return true;

	switch (errno)
	{
	case ENOENT:
		// parent didn't exist, try to create it
	{
		int pos = path.find_last_of('/');
		if (pos == std::string::npos)
#if defined(_WIN32)
			pos = path.find_last_of('\\');
		if (pos == std::string::npos)
#endif
			return false;
		if (!makePath(path.substr(0, pos)))
			return false;
	}
	// now, try to create again
#if defined(_WIN32)
	return 0 == _mkdir(path.c_str());
#else 
	return 0 == mkdir(path.c_str(), mode);
#endif

	case EEXIST:
		// done!
		return isDirExist(path);

	default:
		return false;
	}
}

void create_input() {

}

bool str_is_equal(std::string str1, std::string str2) {
	return str1.compare(str2) == 0;
}

/*Need to read many params
list order here
*/
void read_input_ising(std::ifstream* file_p, class_mc_params* params) {
	std::string line;
	int dummyi;
	double dummyd;
	std::stringstream iss;

	if (file_p->is_open()) {
		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#lattice parameters

		std::getline(*file_p, line);
		iss << line;//#dimension
		iss >> params->dim;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#type
		iss >> params->lattice;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#lengths
		for (int i = 0; i < params->dim; ++i) {
			iss >> dummyi;
			params->lengths.push_back(dummyi);
		}
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#spacings
		for (int i = 0; i < params->dim; ++i) {
			iss >> dummyd;
			params->spacings.push_back(dummyd);
		}
		iss.str("");

		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#blank

		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#Model Parameters

		std::getline(*file_p, line);
		iss << line;//#model
		iss >> params->model;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#Jcouples
		for (int i = 0; i < params->dim; ++i) {
			iss >> dummyd;
			params->Js.push_back(dummyd);
		}
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#beta
		iss >> params->beta;
		iss.str("");

		params->kT = 1 / params->beta;

		std::getline(*file_p, line);
		iss << line;//#h
		iss >> params->h;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#blank

		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#Algorithm Parameters

		std::getline(*file_p, line);
		iss << line;//#algorithm
		iss >> params->alg;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#random seed
		iss >> params->rand_seed;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#blank

		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#Simulation Parameters

		std::getline(*file_p, line);
		iss << line;//#eq time
		iss >> params->eq_time;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#steps per measure
		iss >> params->steps_per_measure;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#measures per dump
		iss >> params->measures_per_dump;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#max dumps
		iss >> params->max_dumps;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#blank
		iss.str("");
	}
	else {
		std::cout << "Error: input file not opened\n";
	}
}

void read_input_spin_boson(std::ifstream* file_p, spin_boson_params* params) {
	std::string line;
	int dummyi;
	double dummyd;
	std::stringstream iss;

	if (file_p->is_open()) {
		std::getline(*file_p, line);
		iss << line;//#Spin Boson Params
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#g
		iss >> params->g;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#A0
		iss >> params->A0;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#delta
		iss >> params->delta;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#v
		iss >> params->v;
		iss.str("");
	}
	else {
		std::cout << "Error: input file not opened\n";
	}
}

void apply_spin_boson_params(class_mc_params* params) {
	double tc = params->beta / ((double)params->lengths[1]);
	params->spacings[1] = tc;
	params->Js[0] = tc*params->Js[0];
	params->Js[1] = -0.5 * log(tc * params->sbparams.delta);
}

template<typename T>
std::string vec2str(std::vector<T> vec) {
	std::stringstream ss;
	for (int i = 0; i < vec.size() - 1; ++i) {
		ss << vec[i] << ", ";
	}
	ss << vec.back();
	return ss.str();
}

void write_outputs(int dump_num, std::vector<int> steps, std::vector<double> times, std::vector<double> record1, std::vector<double> record2, std::vector<double> record3) {
	//char* dump_path = "./dump";
	makePath("./dump");
	char dump_name[100];
	sprintf(dump_name, "dump/dump%d.csv", dump_num);
	std::ofstream file;
	file.open(dump_name);
	file << "Steps," << vec2str(steps) << "\n";
	file << "Times," << vec2str(times) << "\n";
	file << "<Sz>," << vec2str(record1) << "\n";
	file << "<Sx>," << vec2str(record2) << "\n";
	file << "corr_t," << vec2str(record3) << "\n";
	file.close();
}

void write_state(int state_num, IsingLattice2D lat, double action) {
	//char* dump_path = "./dump";
	makePath("./dump");
	char dump_name[100];
	sprintf(dump_name, "dump/state%d.csv", state_num);
	std::ofstream file;
	file.open(dump_name);
	file << "Action," << action << "\n";
	file << lat.to_string();
	file.close();

}
