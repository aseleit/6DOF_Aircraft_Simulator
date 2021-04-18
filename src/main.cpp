#include <iostream>
#include "Aircraft.h"

int main()
{   
    std::cout << "6 DOF Aircraft Simulator" << std::endl;

    Aircraft delta_ac("../data/DELTA_Aircraft.txt");

    // Test that aircraft object has been constructed and initialized sucessfully
    Mat A_long = delta_ac.GetAlong();
    Mat B_long = delta_ac.GetBlong();
    Mat A_lat = delta_ac.GetAlat();
    Mat B_lat = delta_ac.GetBlat();

    std::cout << "A_long(1, 1) = " << A_long(0, 0) << std::endl;
    std::cout << "A_long(1, 4) = " << A_long(0, 3) << std::endl;
    std::cout << "A_long(2, 2) = " << A_long(1, 1) << std::endl;
    std::cout << "B_long(1, 2) = " << B_long(0, 1) << std::endl;
    std::cout << "B_long(3, 1) = " << B_long(2, 0) << std::endl;
    std::cout << "A_lat(1, 1) = " << A_lat(0, 0) << std::endl;
    std::cout << "A_lat(2, 3) = " << A_lat(1, 2) << std::endl;
    std::cout << "A_lat(5, 5) = " << A_lat(4, 4) << std::endl;
    std::cout << "B_lat(1, 2) = " << B_lat(0, 1) << std::endl;
    std::cout << "A_lat(3, 1) = " << B_lat(2, 0) << std::endl;

    return 0; 
}