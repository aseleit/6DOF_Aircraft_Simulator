#include "Aircraft.h"
#include <fstream>
#include <iostream>
#include <math.h>

Aircraft::Aircraft(const std::string &aircraft_data_path) {
    // Load the aircraft parameters from the text file
    std::string line;
    std::string parameter;
    float value;
    std::ifstream file_stream(aircraft_data_path);
    if (file_stream.is_open()) 
    {
        while (file_stream >> parameter >> value)
        {
            if (parameter.compare("mass") == 0) {
                mass_ = value;
            } else if (parameter.compare("Xu") == 0) {
                Xu_ = value;
            } else if (parameter.compare("Xw") == 0) {
                Xw_ = value;
            } else if (parameter.compare("Zu") == 0) {
                Zu_ = value;
            } else if (parameter.compare("Zw") == 0) {
                Zw_ = value;
            } else if (parameter.compare("Mu") == 0) {
                Mu_ = value;
            } else if (parameter.compare("Mw") == 0) {
                Mw_ = value;
            } else if (parameter.compare("Mq") == 0) {
                Mq_ = value;
            } else if (parameter.compare("Yv") == 0) {
                Yv_ = value;
            } else if (parameter.compare("Yp") == 0) {
                Yp_ = value;
            } else if (parameter.compare("Yr") == 0) {
                Yr_ = value;
            } else if (parameter.compare("Lv") == 0) {
                Lv_ = value;
            } else if (parameter.compare("Lp") == 0) {
                Lp_ = value;
            } else if (parameter.compare("Lr") == 0) {
                Lr_ = value;
            } else if (parameter.compare("Nv") == 0) {
                Nv_ = value;
            } else if (parameter.compare("Np") == 0) {
                Np_ = value;
            } else if (parameter.compare("Nr") == 0) {
                Nr_ = value;
            } else if (parameter.compare("Xde") == 0) {
                Xde_ = value;
            } else if (parameter.compare("Zde") == 0) {
                Zde_ = value;
            } else if (parameter.compare("Mde") == 0) {
                Mde_ = value;
            } else if (parameter.compare("Xdt") == 0) {
                Xdt_ = value;
            } else if (parameter.compare("Mdt") == 0) {
                Mdt_ = value;
            } else if (parameter.compare("Lda") == 0) {
                Lda_ = value;
            } else if (parameter.compare("Nda") == 0) {
                Nda_ = value;
            } else if (parameter.compare("Ydr") == 0) {
                Ydr_ = value;
            } else if (parameter.compare("Ldr") == 0) {
                Ldr_ = value;
            } else if (parameter.compare("Ndr") == 0) {
                Ndr_ = value;
            } else if (parameter.compare("theta0") == 0) {
                theta0_ = value;
            } else if (parameter.compare("alpha0") == 0) {
                alpha0_ = value;
            } else if (parameter.compare("U0") == 0) {
                U0_ = value;
            } else if (parameter.compare("W0") == 0) {
                W0_ = value;
            } else {
                std::cout << "Parsing Error: No matching parameter" << std::endl;
            }      
        }
    }
    
    // Initialize system matrices
    double A_long_param_list[16] = {Xu_, Xw_, 0, -1 * GRAVITY * cos(theta0_),
                                    Zu_, Zw_, U0_, -1* GRAVITY * sin(theta0_),
                                    Mu_, Mw_, Mq_, 0,
                                    0, 0, 1, 0};
    double B_long_param_list[8] = {Xde_, Xdt_,
                                   Zde_, 0,
                                   Mde_, Mdt_,
                                   0, 0};
    double A_lat_param_list[25] = {Yv_, Yp_, -1 * (U0_ - Yr_), GRAVITY * cos(theta0_), 0,
                                   Lv_, Lp_, Lr_, 0, 0,
                                   Nv_, Np_, Nr_, 0, 0,
                                   0, 1, 0, 0, 0,
                                   0, 0, 1 / cos(theta0_), 0, 0};
    double B_lat_param_list[10] = {0, Ydr_,
                                   Lda_, Ldr_,
                                   Nda_, Ndr_,
                                   0, 0,
                                   0, 0};

    A_long_ = Mat(4, 4, A_long_param_list);
    B_long_ = Mat(4, 2, B_long_param_list);
    A_lat_ = Mat(5, 5, A_lat_param_list);
    B_lat_ = Mat(5, 2, B_lat_param_list);

    std::cout << "Aicraft object has been constructed." << std::endl;
}

Mat Aircraft::GetAlong() {
    return A_long_;
}

Mat Aircraft::GetBlong() {
    return B_long_;
}

Mat Aircraft::GetAlat() {
    return A_lat_;
}

Mat Aircraft::GetBlat() {
    return B_lat_;
}