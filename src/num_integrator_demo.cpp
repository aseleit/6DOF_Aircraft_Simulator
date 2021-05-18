#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "lib/numericalIntegrator/ERK.h"
#include "Aircraft.h"

#define PI 3.14159265

//Create an aircraft object from data file
Aircraft delta_ac("../data/DELTA_Aircraft.txt");
Mat A_long = delta_ac.GetAlong();
Mat B_long = delta_ac.GetBlong();

//Define simulation parameters
double t_sim = 10.0;
double time_step = 0.1;

//Define the system ODEs
Mat Long_Dyn_xdot(const double t, Mat x, Mat u)
{	
    Mat xdot = A_long * x + B_long * u;
    return xdot;	
};

int main()
{   
    //Define initial state
    Mat x0(4, 0.0);
    x0.Print("Initial State");

    //Define control input
    Mat u_vec(100, 0.0);
    double input_amplitude = 5 * PI / 180;
    double input_freq_rps = PI;
    for (int i = 0; i < t_sim/time_step; i++) {
        u_vec(i) = input_amplitude * sin(input_freq_rps * time_step * i);
    }
    u_vec.Print("U Vector");
    Mat ctrl_input = u_vec.Transpose();

    //Perform numerical integration
    RK4 Long_Dyn_Integrate(Long_Dyn_xdot, 4, 0.0, t_sim, x0, time_step, ctrl_input);
    Long_Dyn_Integrate.Integrate();

    //Display results
    Long_Dyn_Integrate.Y_Out_.Print("Output Integration");
    Long_Dyn_Integrate.y_out_.Print("Output");

    //Create a csv file and save output data
    std::ofstream outputs_file;
    std::string filename = "sim_out.csv";
    outputs_file.open(filename);
    outputs_file << "time" << "," << "pitch_rate" << std::endl;
    for (int i = 0; i <= t_sim/time_step; i++) {
        if (i == 0) {
            outputs_file << 0.0 << "," << x0(2) << std::endl;
        } else {
            outputs_file << time_step * i << "," << Long_Dyn_Integrate.Y_Out_(2, i) * 180.0 / PI << std::endl;
        }
        
    }

    //Close the output csv file
    outputs_file.close();
    
    return 0; 
}

