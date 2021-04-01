#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

typedef boost::array< double , 4 > state_type;
typedef boost::array< double, 1 > control_type;

//control_type u;
double u[250];
void lorenz( const state_type &x , state_type &dxdt, double t )
{
    // redefining variables for readability
    double &x1dot = dxdt[0];
    double &x2dot = dxdt[1];
    double &x3dot = dxdt[2];
    const double &x1 = x[0];
    const double &x2 = x[1];
    const double &x3 = x[2];
    const double &usys = x[3];

    // System of equations
    x1dot = sigma * ( x2 - x1 ) ;
    x2dot = R * x1 - x2 - x1 * x3; 
    x3dot = -b * x3 + x1 * x2;

    ofstream data_log;
    data_log.open("Data Log.txt", std::ios_base::app);
    data_log << "Time: " << t << endl;
    //for(int i = 0; i<250; i++){
    data_log << "inside system control: " << usys << endl;
   // }
}

void write_lorenz( const state_type &x , const double t )
{
    ofstream data_log;
    data_log.open("Data Log.txt", std::ios_base::app);
    data_log << "t: " << t << '\t' << "x1: " << x[0] << '\t' << "x2: " << x[1] << '\t' << "x3: " << x[2] << endl;
}

int main(int argc, char **argv)
{       
    ofstream data_log;
    data_log.open("Data Log.txt");
    data_log << "Integration Output" << endl;
    
    
    
    data_log << "--- Control History ----" << endl;
    for(int i = 0; i<250; i++){
    u[i] = i;
    data_log << u[i] << endl;
    }
    data_log << " ---------------------" << endl;
    double *uptr = u;   
    cout << "u: " << uptr << "\t" << "*u:" << *uptr << endl;
    state_type x = { 10.0 , 1.0 , 1.0, *uptr}; // initial conditions
    //integrate( lorenz , x, 0.0 , 25.0 , 0.1 , write_lorenz ); //Variable time-step integration


    runge_kutta4< state_type > stepper;
    const double dt = 0.1;
    int i = 0;
    for(double t = 0; t<25; t += dt){
        x[3] = u[i];i++;
        stepper.do_step(lorenz,x,t,dt);
        write_lorenz(x,t);
    }
    data_log.close();
}