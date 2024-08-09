#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <tuple>
#include <random>
#include "pe_i.h"
#include "generate_positions.h"
#include "random.h"
#include "output_func.h"
#include "mc_eq.h"
#include "pe_total.h"
#include "mc_move.h"

using namespace std;

int main(int argc, char* argv[]) {
    //(void)argc;
    if (argc < 7) {
        cerr << "Usage: " << argv[0] << " <output_name> <total_n_atoms> <T> <rho> <seed> <num_moves>" << endl;
        return 1;  // Exit with error code indicating incorrect usage
    }

    string output_name = argv[1];
    int total_n_atoms = stoi(argv[2]);
    double T = stod(argv[3]);
    double rho = stod(argv[4]);
    int seed = stoi(argv[5]);
    int num_moves = stoi(argv[6]);
    

    // Declare variables
    vector<vector<double> > box_dim(3, vector<double>(3));
    box_dim[0][0] = 20;
    box_dim[0][1] = 0;
    box_dim[0][2] = 0;
    box_dim[1][0] = 0;
    box_dim[1][1] = 20;
    box_dim[1][2] = 0;
    box_dim[2][0] = 0;
    box_dim[2][1] = 0;
    box_dim[2][2] = 20;  
    
    vector<vector<double> > positions(total_n_atoms, vector<double>(3));
    
    const double k = 8.617e-5; // Boltzmann constant in eV/K
    const double beta = 1/(k*T); // Inverse temperature
    const double e = (k*T)/1.2; // epsilon in eV
    const double V = box_dim[0][0]*box_dim[1][1]*box_dim[2][2]; 
    const double s_3 = rho*V/(total_n_atoms-1); 
    const double s = pow(s_3, 1.0/3.0); //sigma in Angstrom

    mt19937 gen(seed);

    // generate positions
    //positions.clear();
    generate_positions(positions, total_n_atoms, rho, box_dim, gen);

    //initial energy
    double E = total_e(e, box_dim, s, positions, total_n_atoms);
    double sum_E = E;

    //print the initial positions
    print_CONTCAR(output_name, total_n_atoms, box_dim, positions);

    //within the loop
    int trials = 0;
    double step_size = 1.0;
    int accepted_moves = 0;
    
    ofstream PE("PE_eq.csv");
    PE << "PE, avg_PE, time" << endl;

    // Record the start time
    // Start time
    auto start_time = chrono::system_clock::now();
    //auto start_time_str = chrono::system_clock::to_time_t(start_time);

    //equilibration
    int num_eq = 10000;
    int total_accepted_moves = 0;
    for(int j = 1; j < num_eq+1; j++) {
        mc_eq(e, beta, s, box_dim, positions, total_n_atoms,step_size, trials, accepted_moves, total_accepted_moves, E, sum_E, gen);  
        // Record the current time and calculate elapsed time
        auto current_time = chrono::system_clock::now();
        chrono::duration<double> elapsed = current_time - start_time;
        if (total_accepted_moves!=0){
            PE << E << " ," << sum_E/total_accepted_moves << " ," << elapsed.count() << endl;
        }
    }
    
    // End time
    auto end_time = chrono::system_clock::now();
    //auto end_time_str = chrono::system_clock::to_time_t(end_time);
    
    // Print processing time
    chrono::duration<double> elapsed_time = end_time - start_time;
    cout << "Eq Processing time: " << fixed << setprecision(6) << elapsed_time.count() << " seconds" << endl;
    print_CONTCAR("EQUBM", total_n_atoms, box_dim, positions);

    //insertion
    trials = 0;
    step_size = 1.0;
    accepted_moves = 0;
    total_accepted_moves = 0;
    double w = 0;

    ofstream PE_wp("PE_wp.csv");
    PE_wp << "PE, avg_PE, time" << endl;

    // Record the start time
    // Start time
    auto start_time_wp = chrono::system_clock::now();
    //auto start_time_str_wp = chrono::system_clock::to_time_t(start_time_wp);

    for(int j = 1; j < num_moves+1; j++) {
        mc_move(e, beta, s, box_dim, positions, total_n_atoms,step_size, trials, accepted_moves, total_accepted_moves, E, w, gen);  
        // Record the current time and calculate elapsed time
        auto current_time_wp = chrono::system_clock::now();
        chrono::duration<double> elapsed_wp = current_time_wp - start_time_wp;
        if (total_accepted_moves!=0){
            PE_wp << E << " ," << sum_E/total_accepted_moves << " ," << elapsed_wp.count() << endl;
        }
    }
    cout << "Num of moves for insertion: " << total_accepted_moves << endl;

    double acceptance_ratio = static_cast<double>(total_accepted_moves) / num_moves;
    cout << "avg Acceptance Ratio: " << acceptance_ratio << endl;
    double avg_w =  w/total_accepted_moves;

    double mu_excess = -log(avg_w/total_n_atoms)/beta;
    cout << "mu_excess: " << mu_excess << endl;

    // End time
    auto end_time_wp = chrono::system_clock::now();
    //auto end_time_str_wp = chrono::system_clock::to_time_t(end_time_wp);
    //cout << "Wp End time: " << put_time(localtime(&end_time_str_wp), "%Y-%m-%d %X") << endl;
    
    // Print processing time
    chrono::duration<double> elapsed_time_wp = end_time_wp - start_time_wp;
    cout << "Wp Processing time: " << fixed << setprecision(6) << elapsed_time_wp.count() << " seconds" << endl;
    
    return 0;
}

//analytical vs the code PE
//equation of state vs the test particle insertion method
//the pressure of the system must change to accomodate the inserted atom

//case study 1: equation of state of LJ fluid 719
//case study 15: widom particle insertion method
//3.3.2 Technical details 63
//5.1.5 Pressure 130




