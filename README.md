# Six DOF Aircraft Simulator

The goal of this project is to implement a six degrees of freedom (DOF), nonlinear simulation for fixed-wing aircraft. 
The first iteration of this project is using a linear aircraft model. The second iteration will be using a nonlinear model.
Full documentation is available in `doc/6DOF_Aircraft_Simulator.pdf`.

## Dependencies for Running Locally
* CMake: v3.7 or higher
* Make: v4.1 or higher
* gcc/g++: v5.4 or higher
* gnuplot (optional -  if you want to plot the results)
  * [Click here for installation instructions](http://www.gnuplot.info/)

## Basic Build Instructions
1. Clone this repo
2. Make a build directory in the top level directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./6DOF_Sim_test`
5. To plot the results (assuming gnuplot is installed): `gnuplot ../src/plot_results.gnuplot`

## Results Snippet
![Translational Velocities](./doc/sim_results_uvw.png?raw=true)
![Rotational Velocities](./doc/sim_results_pqr.png?raw=true)
![Attitude](./doc/sim_results_phi_theta_psi.png?raw=true)