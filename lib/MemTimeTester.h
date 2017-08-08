#pragma once
#include <ctime>
#include <chrono>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

struct timer_list {
	int size;
	std::vector<std::string> names;
	std::vector<double> times;
};

class MemTimeTester {

	std::vector<std::chrono::time_point<std::chrono::system_clock>> start_times;
	std::vector<std::chrono::time_point<std::chrono::system_clock>> end_times;
	std::vector<std::string> names;
	std::vector<double> times;

public:

	MemTimeTester() {
		//don't initialize the names, just allow the user to add them as they go
	}

	void flag_start_time(std::string identifier) {
		bool found = false;
		for (int i = 0; i < names.size(); ++i) {
			if (identifier.compare(names[i]) == 0) {
				//if the given identifier already has an entry, start counting from its current elapsed time
				start_times[i] = std::chrono::system_clock::now();
				found = true;
			}
		}
		if (!found) {
			//if the given identifier is not already listed, start a new timer
			names.push_back(identifier);
			times.push_back(0.0);
			start_times.push_back(std::chrono::system_clock::now());
			end_times.push_back(std::chrono::system_clock::now());
		}
	}

	void flag_end_time(std::string identifier) {
		int name_index = -1;
		for (int i = 0; i < names.size(); ++i) {
			if (identifier.compare(names[i]) == 0) {
				name_index = i;
			}
		}
		if (name_index == -1) {
			std::cout << "\n" << identifier << " is not a valid timer identifier\n";
		}
		else {
			end_times[name_index] = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed;
			elapsed = end_times[name_index] - start_times[name_index];
			times[name_index] += elapsed.count();
		}
	}

	double get_running_time(std::string identifier) {
		int name_index = -1;
		for (int i = 0; i < names.size(); ++i) {
			if (identifier.compare(names[i]) == 0) {
				name_index = i;
			}
		}
		if (name_index == -1) {
			std::cout << "\n" << identifier << " is not a valid timer identifier\n";
			return 0;
		}
		else {
			std::chrono::duration<double> elapsed;
			elapsed = std::chrono::system_clock::now() - start_times[name_index];
			return elapsed.count();
		}
	}

	void print_timers() {
		std::stringstream outstring;
		for (int i = 0; i < names.size(); ++i) {
			outstring << names[i] << ": " << times[i] << "s\n";
		}
		std::cout << outstring.str();
	}

	timer_list get_timers() {
		timer_list result;
		std::chrono::duration<double> elapsed;
		for (int i = 0; i < names.size(); ++i) {
			elapsed = end_times[i] - start_times[i];
			result.names.push_back(names[i]);
			result.times.push_back(elapsed.count());
		}
		result.size = result.names.size();
		return result;
	}

	void clear_timers() {
		names = {};
		start_times = {};
		end_times = {};
	}

};