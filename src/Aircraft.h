#ifndef AIRCRAFT_H
#define AIRCRAFT_H

#include <string>
#include "lib/matrixLibrary/Mat.h"

#define GRAVITY 9.81
#define PI 3.14159265

class Aircraft {
    public:
        Aircraft(const std::string &aircraft_data_path);
        
        Mat GetAlong();
        Mat GetBlong();
        Mat GetAlat();
        Mat GetBlat();

    private:
        // Aircraft parameters
        double mass_;
        double Xu_;
        double Xw_;
        double Zu_;
        double Zw_;
        double Mu_;
        double Mw_;
        double Mq_;
        double Yv_;
        double Yp_;
        double Yr_;
        double Lv_;
        double Lp_;
        double Lr_;
        double Nv_;
        double Np_;
        double Nr_;
        double Xde_;
        double Zde_;
        double Mde_;
        double Xdt_;
        double Mdt_;
        double Lda_;
        double Nda_;
        double Ydr_;
        double Ldr_;
        double Ndr_;
        double theta0_;
        double alpha0_;
        double U0_;
        double W0_;

        // Longitudinal and Lateral linear systems matrices 
        Mat A_long_;
        Mat A_lat_;
        Mat B_long_;
        Mat B_lat_;
};

#endif