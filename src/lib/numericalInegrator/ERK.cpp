
#include "ERK.h"

void RK4::Step (double& t, Mat& y, double& h) {
		
	// Elementary RK4 step
	double hh=h*0.5;
	double h6=h/6.0;
	double th=t+hh;


 	f( t   , y       , k_1, pAux );
	f( t+hh, y+hh*k_1, k_2, pAux );
	f( t+hh, y+hh*k_2, k_3, pAux );
	f( t+h , y+h*k_3 , k_4, pAux );
  
	y = y + h6*( k_1 + 2.0*k_2 + 2.0*k_3 + k_4 );

	nFevals +=4;
	nSteps +=1;
  // Update independent variable

	 t = t+h;
}

void RK4::Integrate()
{
	clock_t t1=clock();
	
	tOut = tInitial ;
	yOut = y0;
	for (int i = 0 ; i < (tFinal-tInitial)/h ;i++)
	{
		Step(tOut,yOut,h);
	}

	time += (clock()-t1);
}



void RK8::Step (double& t, Mat& y, double& h) {
		
	// Elementary RK4 step
	
  // Elementary RK4 step
  
	
	f( t           , y                                                                                                , k_1,   pAux );
	f( t+h*4.0/27.0, y+(k_1*4.0/27.0)*h                                                                               , k_2,   pAux );
	f( t+h*2.0/9.0 , y+(k_1/18.0   +k_2/6.0)*h                                                                        , k_3,   pAux );
	f( t+h/3.0     , y+(k_1/12.0   +k_3*0.25)*h                                                                       , k_4,   pAux );
	f( t+h*0.5     , y+(k_1/8.0    +k_4*3.0*0.125)*h                                                                  , k_5,   pAux );
	f( t+h*2.0/3.0 , y+(k_1*13.0   -k_3*27.0 +k_4*42.0  +k_5*8.0)*h/54.0                                              , k_6,   pAux );
	f( t+h*1.0/6.0 , y+(k_1*389.0  -k_3*54.0 +k_4*966.0 -k_5*824.0  +k_6*243.0)*h/4320.0                              , k_7,   pAux );
	f( t+h         , y+(k_1*-231.0 +k_3*81.0 -k_4*1164.0+k_5*656.0  -k_6*122.0+k_7*800.0)*h*0.05                      , k_8,   pAux );
	f( t+h*5.0/6.0 , y+(k_1*-127.0 +k_3*18.0 -k_4*678.0 +k_5*456.0  -k_6*9.0  +k_7*576.0 +k_8*4.0)*h/288.0            , k_9,   pAux );
	f( t+h         , y+(k_1*1481.0 -k_3*81.0 +k_4*7104.0-k_5*3376.0 +k_6*72.0 -k_7*5040.0-k_8*60.0+k_9*720.0)*h/820.0 , k_10,  pAux );

	y = y + ( k_1*41.0 + k_4*27 + k_5*272.0 + k_6*27.0+k_7*216.0+ k_9*216.0+k_10*41.0 )*(h/840.0);

  // Update independent variable

	t = t + h;


	nFevals +=10;
	nSteps +=1;
  // Update independent variable

	 t = t+h;
}

void RK8::Integrate()
{
	clock_t t1=clock();
	
	tOut = tInitial ;
	yOut = y0;
	for (int i = 0 ; i < (tFinal-tInitial)/h ;i++)
	{
		Step(tOut,yOut,h);
	}

	time += (clock()-t1);
}
