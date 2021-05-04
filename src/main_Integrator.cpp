#include <iostream>
#include "lib/numericalIntegrator/ERK.h"
#include <cmath>


Mat twoBodyNumStates(const double t, Mat y,  Mat u)
{	
    double GM_Earth =1.0;
	Mat r(3),v(3);
    Mat dydt(6);
	for( int i =0 ; i < 3 ; i++)
	{
		r(i) = y(i);
		v(i) = y(i+3) ;
	}    

    Mat acc = -r*GM_Earth/pow(r.Norm(),3);
    for( int i =0 ; i < 3 ; i++)
	{
		dydt(i) = v(i);
		dydt(i+3) = acc(i) ;
    } 
    return dydt;	
};

void exponential_function(const double t, const Mat y, Mat& dydt, Mat u)
{	
	dydt=y;	
};

Mat exponential_function2(const double t, Mat y,  Mat u)
{	
	Mat dydt=y;	
    return dydt;
};

int main()
{   
    Mat y0(1,0,0,0,1,0);
    cout<< "d";
    y0.Print(" Initial U");
    Mat u_vec(0.0,10.0,11);
    u_vec.Print(" U Vector ");
    Mat uu = u_vec.Transpose();
    RK4 Orbit_Integrate(twoBodyNumStates,6,0.0,2*3.1,y0,0.2,uu);
    Orbit_Integrate.Integrate();

    Orbit_Integrate.Y_Out_.Print("Output Integration");
    Orbit_Integrate.y_out_.Print("Output");
    return 0; 
}



/*
int main()
{   
    Mat y0(1);
    y0(0)=1.0 ;
    Mat u_vec(0.0,10.0,101);
    u_vec.Print(" U Vector ");
    Mat uu = u_vec.Transpose();
    RK4 Orbit_Integrate(exponential_function2,1,0.0,10.0,y0,0.1,uu);
    Orbit_Integrate.Integrate();

    Orbit_Integrate.Y_Out_.Print("Output Integration");
    return 0; 
}
*/