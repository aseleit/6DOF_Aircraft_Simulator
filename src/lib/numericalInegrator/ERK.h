#ifndef _ERK_H_
#define _ERK_H_


#include "Ode.h"

class RK4 {

	 public:
	// Constructor
    RK4 (

	  Derivs   f_,        // Differential equation
	  int      n_eqn_,    // Dimension
	  double   tInitial_,
	  double   tFinal_,
	  Mat      y0_,
	  double   h_,
	  void*    pAux_      // Pointer to auxiliary data	  		
      ) 
    : f(f_), n_eqn(n_eqn_),tInitial(tInitial_),tFinal(tFinal_),y0(y0_),h(h_), pAux(pAux_)
    {
		clock_t t1 = clock();
		k_1=Mat(n_eqn) , k_2=k_1 , k_3=k_2 , k_4=k_3;
	
		yOut=Mat(n_eqn);
				
		nFevals = 0;
		nSteps  = 0;
		nFails  = 0;

		time = clock()-t1;
	};



    // Integration step
    void Step (         
      double&  t,         // Value of the independent variable; updated by t+h
      Mat&     y,         // Value of y(t); updated by y(t+h)
      double   &h          // Step size
    );

	void Integrate();

	// Elements for Output
	int nFevals;
	int nSteps;
	int nFails;

	Mat yOut;

	double tOut;
	double time;
	// Elements
	Derivs		 f;
    int			 n_eqn;
	double		 tInitial,tFinal;
	Mat			 y0;
	
	double		 h;
	void*		 pAux;

    Mat			 k_1,k_2,k_3,k_4;


//private:
};



class RK8 {

	 public:
	// Constructor
    RK8 (

	  Derivs   f_,        // Differential equation
	  int      n_eqn_,    // Dimension
	  double   tInitial_,
	  double   tFinal_,
	  Mat      y0_,
	  double   h_,
	  void*    pAux_      // Pointer to auxiliary data	  		
      ) 
    : f(f_), n_eqn(n_eqn_),tInitial(tInitial_),tFinal(tFinal_),y0(y0_),h(h_), pAux(pAux_)
    {
		clock_t t1 = clock();
		k_1=Mat(n_eqn) , k_2=k_1 , k_3=k_2 , k_4=k_3 , k_5=k_4 , k_6=k_5 , k_7=k_6  , k_8=k_7 , k_9=k_8 , k_10=k_9;
	
		yOut=Mat(n_eqn);
				
		nFevals = 0;
		nSteps  = 0;
		nFails  = 0;

		time = clock()-t1;
	};



    // Integration step
    void Step (         
      double&  t,         // Value of the independent variable; updated by t+h
      Mat&     y,         // Value of y(t); updated by y(t+h)
      double   &h          // Step size
    );

	void Integrate();

	// Elements for Output
	int nFevals;
	int nSteps;
	int nFails;

	Mat yOut;

	double time;
	double tOut;
	// Elements
	Derivs		 f;
    int			 n_eqn;
	double		 tInitial,tFinal;
	Mat			 y0;
	
	double		 h;
	void*		 pAux;

    Mat			 k_1,k_2,k_3,k_4,k_5,k_6,k_7,k_8,k_9,k_10;


//private:
};


#endif