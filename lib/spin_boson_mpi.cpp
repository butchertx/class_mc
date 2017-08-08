/**
This code is to run the full version of the spin-boson simulation, using MPI.
Some of the data that was previously held in classes have been moved so that
they are easier to transfer around and manipulate.  Now, we will run parallel
simulations on each MPI unit (be it processor or node), with a GPU (or two?) 
working in tandem to calculate the action for parallel tempering moves. Each
unit will hold the lattice itself, its parameters (which will be exchanged in 
the parallel tempering), and the corresponding interactions.  The requirement
will be that the lattices are all the same size.  The main difference between
this and the version that uses the "LongRangeWolff" class is that it will allow
the user to define their own Monte Carlo step function, which makes it much more
flexible.
**/

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include "IsingLattice2D.h"
#include "class_mc_io.h"
#include "LongRangeWolff2D.h"
#include "Matrix.h"
#include "MemTimeTester.h"
#include "obs_calc_fast.cuh"
extern "C" {
#include "random.h"
}
#include <mpi.h>

using namespace std;

thrust::host_vector<double> IsingLattice2D::get_thrust_vector(){
	thrust::host_vector<double> result(Lx*Ly);
	for (int i = 0; i < Lx; ++i){
		for (int j = 0; j < Ly; ++j){
			result[i*Ly + j] = (double) spins[i][j];
		}
	}
	return result;
}

thrust::host_vector<double> GeneralLRW::get_thrust_interactions(){
	thrust::host_vector<double> result(interactions.size()*interactions[0].size());
	for (int i = 0; i < interactions.size(); ++i){
		for (int j = 0; j < interactions[i].size(); ++j){
			result[i*interactions[i].size() + j] = interactions[i][j];
		}
	}
	return result;
}

void ptemp(int num_procs, int id, double action, double alpha, IsingLattice2D& lat){
    //use MPI to perform a parallel tempering step
    //when it comes time to send the lattices, master process will send first
    //and child processes will all receive first.  This might not be the quickest way
    //to do it but it ensures nothing gets locked up waiting to send/receive
    //
    //lattice message tags will be the receiving id
    MPI_Status Stat[2];
	MPI_Request req[2];
	std::vector<int> lat_buffer_out = lat.get_bool_spins();
	std::vector<int> lat_buffer_in(lat.get_N());

    //master process
    //1. receive action from other processes
    //2. test for switches, storing the send-to id in elements of an array
    //3. send the send-to id to each other process
    //4. send lattice if send-to[0] != 0
    if(id == 0){
        std::vector<double> actions(num_procs);
        std::vector<double> alphas(num_procs);
        std::vector<int> receive_from(num_procs);
        std::vector<int> send_to(num_procs);
        actions[0] = action;
        alphas[0] = alpha;
        for (int i = 1; i < num_procs; ++i){
            MPI_Recv(&(actions[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            MPI_Recv(&(alphas[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            receive_from[i] = i;
        }

        //test for switches
        //assume for now that action = alpha*C; this will need to be adjust later
        //start with lowest id and travel up
        double prob;
        for (int i = 0; i < num_procs - 1; ++i){
            prob = exp((alphas[receive_from[i]] - alphas[receive_from[i + 1]]) * (actions[receive_from[i]] - actions[receive_from[i + 1]]));
            //cout << "Parallel Tempering compare alphas " << alphas[receive_from[i]] << " and " << alphas[receive_from[i + 1]] << "; and actions " <<
            //            actions[receive_from[i]] << " and " << actions[receive_from[i + 1]] << ", with probability " << prob << "\n";
            if (drand1_() < prob){
                receive_from[i + 1] = receive_from[i];
                receive_from[i] = i + 1;
            }
        }
        //invert receive_from to get send_to
        for (int i = 0; i < num_procs; ++i){
            send_to[receive_from[i]] = i;
        }
        //receive_from[i] now gives the id of the lattice that process i will receive
        //message each process and tell them which process they will send their lattice to and which process will receive their lattice
        //tags for send-to id will be the process id, tags for receive-from id will be num_procs + id
        for (int i = 1; i < num_procs; ++i){
            //send-to
            MPI_Send(&send_to[i], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
            //receive-from
            MPI_Send(&receive_from[i], 1, MPI_DOUBLE, i, num_procs + i, MPI_COMM_WORLD);
        }
        //cout << "MPI master process: send_to = " << send_to[0] << ", receive_from = " << receive_from[0] << "\n";
	
	//send lat_buffer to send_to[0] and receive from receive_from[0]
	if(send_to[0] != 0){
	MPI_Isend(&lat_buffer_out.front(), lat_buffer_out.size(), MPI_INT, send_to[0], send_to[0], MPI_COMM_WORLD, &req[0]);
	MPI_Irecv(&lat_buffer_in.front(), lat_buffer_in.size(), MPI_INT, receive_from[0], 0, MPI_COMM_WORLD, &req[1]);
	MPI_Waitall(2, req, Stat);
	lat.copy_bool_spins(lat_buffer_in);
	}
	}
    else{
    //child processes
    //1. send action to master process
    //2. receive send-to process and receive-from process
    //3. if send-to[id] != id, send lattice and receive new lattice
        action /= alpha;
        MPI_Send(&action, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&alpha, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);

        //cout << "MPI child process sending message action = " << action << ", alpha = " << alpha << "\n";

        int send_to, receive_from;
        MPI_Recv(&send_to, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD, &Stat[0]);
        MPI_Recv(&receive_from, 1, MPI_DOUBLE, 0, num_procs + id, MPI_COMM_WORLD, &Stat[0]);
        //cout << "MPI child process: send_to = " << send_to << ", receive_from = " << receive_from << "\n";

	//receive lattice and send lattice
	if(send_to != id){
        	MPI_Irecv(&lat_buffer_in.front(), lat_buffer_in.size(), MPI_INT, receive_from, id, MPI_COMM_WORLD, &req[0]);
		//cout << "Process id " << id << " received lattice copy\n";
	        MPI_Isend(&lat_buffer_out.front(), lat_buffer_out.size(), MPI_INT, send_to, send_to, MPI_COMM_WORLD, &req[1]);
		MPI_Waitall(2, req, Stat);
		lat.copy_bool_spins(lat_buffer_in);
	}
    }
}

int main(int argc, char* argv[]){
//
//  define global variables
//
    int p;//number of processes
    int id;//ID of this process
    int seed;//RNG seed for each process
    double random_value;//test random value

//
//  Initialize MPI.
//
    MPI_Init ( NULL, NULL );
//
//  Get the number of processes.
//
    MPI_Comm_size ( MPI_COMM_WORLD, &p );
//
//  Get the ID of this process.
//
    MPI_Comm_rank ( MPI_COMM_WORLD, &id );
//
//  After MPI_Init is set, master node should do the following:
//  The master process prints a message.
//  Need a base parameter file: just specify the name of an input file with the parameters
//  Need to know which parameter to vary and how to vary it.  This should be compatible with # of processes
//
    class_mc_params params;
    ifstream infile;
    string vary_param;
    double delta_p;
    infile.open(argv[1]);
    read_input_ising(&infile, &params);
    read_input_spin_boson(&infile, &(params.sbparams));
    read_input_mpi(&infile, &vary_param, &delta_p);
    //change parameter based on id, delta, and given param name.  First, just assume alpha is the one to change
    params.sbparams.A0 += id*delta_p;
    params.rand_seed += id*100;
    apply_spin_boson_params(&params);
    //cout << "Params for process " << id << ":\n" << params.to_string() << "\n\n\n";

//
//  Simulation should be ready to set up for each process with different params and random seed
//  
    IsingLattice2D lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
    IsingLattice2D& latref = lat; //this will be the main workable item, staying local to its original process in the node
    GeneralLRW wolff = GeneralLRW(params);
    wolff.set_mag(latref);
    //cout << "Interactions and Interaction sum for process " << id << ":\n\n";
    //wolff.print_interactions();
    //wolff.print_interaction_sum();
    //wolff.print_site_sums();

//
//  Set up measurement apparatus
//
    class_mc_measurements results;
    results.names = {"loc", "cluster_size"};
    results.values = {{} , {}};
    MemTimeTester timer;
    timer.flag_start_time("full simulation");
//
//  Equilibrate
//
    for (int i = 0; i < params.eq_time; ++i){
        wolff.step(latref);
    }
//
//  Define parallel tempering thrust variables
//
	thrust::host_vector<double> thrustlat = latref.get_thrust_vector();
	thrust::host_vector<double>& thrustlat_ref = thrustlat;
	thrust::host_vector<double> thrustint = wolff.get_thrust_interactions();
	thrust::host_vector<double>& thrustint_ref = thrustint;

//
//  Run the measurement portion of the simulation
//  Use each "dump" as the parallel tempering step
//
	int fast_action;
    for (int dump = 0; dump < params.max_dumps; ++dump){
        for (int measure = 0; measure < params.measures_per_dump; ++measure){
		    timer.flag_start_time("steps");
                for (int n = 0; n < params.steps_per_measure; ++n){
                    wolff.step(latref);
                }
                results.record("loc", wolff.calc_loc(latref));
                results.record("cluster_size", wolff.get_cluster_size());
            timer.flag_end_time("steps");
        }
	    thrustlat = latref.get_thrust_vector();
        //parallel tempering
        timer.flag_start_time("action calculation");
	    fast_action = thrust_calc_action_general(thrustlat_ref, thrustint_ref, latref.get_Lx(), latref.get_Ly());
	    timer.flag_end_time("action calculation");
        timer.flag_start_time("parallel tempering");
        ptemp(p, id, fast_action, params.sbparams.A0, latref);
	    timer.flag_end_time("parallel tempering");
    }
//
//  test fast action calculation
//
    double slow_action = wolff.calc_action_slow(latref);
    fast_action = thrust_calc_action_general(thrustlat_ref, thrustint_ref, latref.get_Lx(), latref.get_Ly());
    std::stringstream outstring;
    outstring << "Printing action test for process " << id << "\nLattice:\n" << latref.to_string() << "\n" << wolff.get_int_string() << "Slow action: " << slow_action << ", Fast action: " << fast_action << "\n";
    cout << outstring.str();


//
//  Write the results
//
    cout << "Results for process #" << id << ": localization = " << mean(results.get_vals("loc")) << ", cluster avg: " << mean(results.get_vals("cluster_size")) << "\n\n";
    timer.flag_end_time("full simulation");
	if(id == 0){
	    timer.print_timers();
	    cout << "total steps = " << params.steps_per_measure*params.measures_per_dump*params.max_dumps << ", parallel tempering steps = " << params.max_dumps << "\n";
	}
	
//
//  End the MPI process
//
    MPI_Finalize();

return 0;
}
