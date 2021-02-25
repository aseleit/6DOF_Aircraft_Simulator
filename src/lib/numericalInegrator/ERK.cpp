
#include "ERK.h"

void RK4::Step (double& t, Mat& y, double& h) {		
	// Elementary RK4 step
	double hh=h*0.5;
	double h6=h/6.0;
	double th=t+hh;
 	f_( t   , y       , k_1_ , pAux_ );
	f_( t+hh, y+hh*k_1_, k_2_, pAux_ );
	f_( t+hh, y+hh*k_2_, k_3_, pAux_ );
	f_( t+h , y+h*k_3_ , k_4_, pAux_ );  
	y = y + h6*( k_1_ + 2.0*k_2_ + 2.0*k_3_ + k_4_ );
	n_fevals_ +=4;
	n_steps_ +=1;  
	 t = t+h; // Update independent variable
}

void RK4::Integrate()
{
	clock_t t1=clock();	
	t_out_ = t_initial_ ;
	y_out_ = y_0_;
	for (int i = 0 ; i < (t_final_-t_initial_)/h_ ;i++)
		Step(t_out_,y_out_,h);
	time_ += (clock()-t1);
}

void RK8::Step (double& t, Mat& y, double& h) {	
	f_( t           , y                                                                                                        , k_1_,   pAux );
	f_( t+h*4.0/27.0, y+(k_1_*4.0/27.0)*h                                                                                      , k_2_,   pAux );
	f_( t+h*2.0/9.0 , y+(k_1_/18.0   +k_2_/6.0)*h                                                                              , k_3_,   pAux );
	f_( t+h/3.0     , y+(k_1_/12.0   +k_3_*0.25)*h                                                                             , k_4_,   pAux );
	f_( t+h*0.5     , y+(k_1_/8.0    +k_4_*3.0*0.125)*h                                                                        , k_5_,   pAux );
	f_( t+h*2.0/3.0 , y+(k_1_*13.0   -k_3_*27.0 +k_4_*42.0  +k_5_*8.0)*h/54.0                                                  , k_6_,   pAux );
	f_( t+h*1.0/6.0 , y+(k_1_*389.0  -k_3_*54.0 +k_4_*966.0 -k_5_*824.0  +k_6_*243.0)*h/4320.0                                 , k_7_,   pAux );
	f_( t+h         , y+(k_1_*-231.0 +k_3_*81.0 -k_4_*1164.0+k_5_*656.0  -k_6_*122.0+k_7_*800.0)*h*0.05                        , k_8_,   pAux );
	f_( t+h*5.0/6.0 , y+(k_1_*-127.0 +k_3_*18.0 -k_4_*678.0 +k_5_*456.0  -k_6_*9.0  +k_7_*576.0 +k_8_*4.0)*h/288.0             , k_9_,   pAux );
	f_( t+h         , y+(k_1_*1481.0 -k_3_*81.0 +k_4_*7104.0-k_5_*3376.0 +k_6_*72.0 -k_7_*5040.0-k_8_*60.0+k_9_*720.0)*h/820.0 , k_10_,  pAux );
	y = y + ( k_1_*41.0 + k_4_*27 + k_5_*272.0 + k_6_*27.0+k_7_*216.0+ k_9_*216.0+k_10_*41.0 )*(h/840.0);
	t = t + h; // Update independent variable
	n_fevals_ +=10;
	n_steps_ +=1;
 }

void RK8::Integrate()
{
	clock_t t1=clock();	
	t_out_ = t_initial_ ;
	y_out_ = y_0_;
	for (int i = 0 ; i < (t_final_-t_initial_)/h_ ;i++)
		Step(t_out_,y_out_,h_);
	time_ += (clock()-t1);
}
