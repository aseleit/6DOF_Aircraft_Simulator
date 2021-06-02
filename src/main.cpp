#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "lib/numericalIntegrator/ERK.h"
#include "Aircraft.h"

#define PI 3.14159265

// Construct an aircraft object and load data from text file
Aircraft delta_ac("../data/DELTA_Aircraft.txt");
Mat A_long = delta_ac.GetAlong();
Mat B_long = delta_ac.GetBlong();
Mat A_lat = delta_ac.GetAlat();
Mat B_lat = delta_ac.GetBlat();

// Define simulation parameters
double t_init = 0.0;
double t_sim = 10.0;
double time_step = 0.1;
int n_eqs_long = 4;
int n_eqs_lat = 5;

// Define the longitudinal linear system
Mat LongDynXdot(const double t, Mat x, Mat u)
{	
    Mat xdot = A_long * x + B_long * u;
    return xdot;	
};

// Define the lateral linear system
Mat LatDynXdot(const double t, Mat x, Mat u)
{	
    Mat xdot = A_lat * x + B_lat * u;
    return xdot;	
};

int main()
{   
    // Define initial state
    Mat x0_long(4, 0.0);
    Mat x0_lat(5, 0.0);

    // Define the control inputs
    Mat elev_input(100, 0.0);
    double elev_amplitude = 5 * PI / 180;
    double elev_freq_rps = PI;
    for (int i = 0; i < t_sim/time_step; i++) {
        elev_input(i) = elev_amplitude * sin(elev_freq_rps * time_step * i);
    }
    Mat ctrl_input_long = elev_input.Transpose();

    Mat ail_input(100, 0.0);
    double ail_amplitude = 5 * PI / 180;
    double ail_freq_rps = PI;
    for (int i = 0; i < t_sim/time_step; i++) {
        ail_input(i) = ail_amplitude * sin(ail_freq_rps * time_step * i);
    }
    Mat ctrl_input_lat = ail_input.Transpose();

    // Perform numerical integration
    RK4 Long_Dyn_Integrate(LongDynXdot, n_eqs_long, t_init, t_sim, x0_long, time_step, ctrl_input_long);
    Long_Dyn_Integrate.Integrate();
    std::cout << "Long dynamics integration finished." << std::endl;

    RK4 Lat_Dyn_Integrate(LatDynXdot, n_eqs_lat, t_init, t_sim, x0_lat, time_step, ctrl_input_lat);
    Lat_Dyn_Integrate.Integrate();
    std::cout << "Lat dynamics integration finished." << std::endl;

    // Display results
    std::cout << "Final long state is: " << std::endl;
    Long_Dyn_Integrate.y_out_.Print("Output");

    std::cout << "Final lat state is: " << std::endl;
    Lat_Dyn_Integrate.y_out_.Print("Output");

    // Create a csv file and save output data
    std::ofstream outputs_file;
    std::string filename = "sim_out.csv";
    outputs_file.open(filename);
    outputs_file << "time" << "," << "u_mps" << "," << "w_mps" << "," 
                 << "q_dps" << "," << "theta_deg" << ","
                 << "v" << "," << "p_dps" << "," << "r_dps" << ","
                 << "phi_deg" << "," << "psi_deg" << std::endl;
    for (int i = 0; i <= t_sim/time_step; i++) {
        if (i == 0) {
            outputs_file << 0.0 << "," << x0_long(0) << "," << x0_long(1) << "," 
                         << x0_long(2) << "," << x0_long(3) << ","
                         << x0_lat(0) << "," << x0_lat(1) << "," << x0_lat(2) << ","
                         << x0_lat(3) << "," << x0_lat(4) << std::endl;
        } else {
            outputs_file << time_step * i << "," 
                         << Long_Dyn_Integrate.Y_Out_(0, i) << ","
                         << Long_Dyn_Integrate.Y_Out_(1, i) << ","
                         << Long_Dyn_Integrate.Y_Out_(2, i) * 180.0 / PI << ","
                         << Long_Dyn_Integrate.Y_Out_(3, i) * 180.0 / PI << ","
                         << Lat_Dyn_Integrate.Y_Out_(0, i) << ","
                         << Lat_Dyn_Integrate.Y_Out_(1, i) * 180.0 / PI << ","
                         << Lat_Dyn_Integrate.Y_Out_(2, i) * 180.0 / PI << ","
                         << Lat_Dyn_Integrate.Y_Out_(3, i) * 180.0 / PI << ","
                         << Lat_Dyn_Integrate.Y_Out_(4, i) * 180.0 / PI << std::endl;
        }
        
    }

    // Close the output csv file
    outputs_file.close();
    
    return 0; 
}