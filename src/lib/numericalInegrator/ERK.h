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
    : f_(f_), n_eqn_(n_eqn_),t_initial_(tInitial_),t_final_(tFinal_),y_0_(y_0_),h_(h_), pAux_(pAux_)
    {
		clock_t t1_ = clock();
		k_1_=Mat(n_eqn_) , k_2_=k_1_ , k_3_=k_2_ , k_4_=k_3_;	
		y_out_=Mat(n_eqn_);				
		n_fevals_ = 0;
		n_steps_  = 0;
		n_fails_  = 0;
		time_ = clock()-t1_;
	};
    // Integration step
    void Step (         
      double&  t_,         // Value of the independent variable; updated by t+h
      Mat&     y_,         // Value of y(t); updated by y(t+h)
      double   &h_          // Step size
    );
	void Integrate();
	// Elements for Output
	int n_fevals_;
	int n_steps_;
	int n_fails_;
	Mat y_out_;
	double t_out_;
	double time_;
	// Elements
	Derivs		 f_;
    int			 n_eqn_;
	double		 t_initial_,t_final_;
	Mat			 y_0_;	
	double		 h_;
	void*		 pAux_;

	private:
	 Mat k_1_,k_2_,k_3_,k_4_;
};

class RK8 {
	public:
	// Constructor
    RK8 (
	  Derivs   f_,        // Differential equation
	  int      n_eqn_,    // Dimension
	  double   t_initial_,
	  double   t_final_,
	  Mat      y_0_,
	  double   h_,
	  void*    pAux_      // Pointer to auxiliary data	  		
      ) 
    : f_(f_), n_eqn_(n_eqn_),t_initial_(t_initial_),t_final_(t_final_),y_0_(y_0_),h_(h_), pAux_(pAux_)
    {
		clock_t t_1_ = clock();
		k_1_=Mat(n_eqn_) , k_2_=k_1_ , k_3_=k_2_ , k_4_=k_3_ , k_5_=k_4_ , k_6_=k_5_;
		k_7_=k_6_  , k_8_=k_7_ , k_9_=k_8_ , k_10_=k_9_;	
		y_out_=Mat(n_eqn_);				
		n_fevals_ = 0;
		n_steps_  = 0;
		n_fails_  = 0;
		time_ = clock()-t_1_;
	};
    // Integration step
    void Step (         
      double&  t_,         // Value of the independent variable; updated by t+h
      Mat&     y_,         // Value of y(t); updated by y(t+h)
      double   &h_          // Step size
    );
	void Integrate();
	// Elements for Output
	int n_fevals_;
	int n_steps_;
	int n_fails_;
	Mat y_out_;
	double time_;
	double t_out_;
	// Elements
	Derivs		 f_;
    int			 n_eqn_;
	double		 t_initial_,t_final_;
	Mat			 y_0_;	
	double		 h_;
	void*		 pAux_;	
	private:
    Mat			 k_1_,k_2_,k_3_,k_4_,k_5_,k_6_,k_7_,k_8_,k_9_,k_10_;
};
#endif