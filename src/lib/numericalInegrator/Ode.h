#ifndef _ODE_H_
#define _ODE_H_


#include "Mat.h"
#include <iomanip>

inline const double Max(const double &a, const double &b) { return b > a ? b : a; }
inline const double Min(const double &a, const double &b) { return b < a ? b : a; }
inline double SQR(const double a) { return a*a;}

struct AuxDataRecord {
  int     n_step;
  double  t;
};

typedef void (*Derivs)(
  const double  x,     // Independent variable
  const Mat&    y,     // State vector 
  Mat&          yp,    // Derivative y'=f(x,y)
  void*         pAux   // Pointer to auxiliary data used within f
);

typedef void (*Derivs2Order)(
  const double  x,     // Independent variable
  Mat           r,     // State vector
  Mat           v,     // State vector 
  Mat&          yp,    // Derivative y'=f(x,y)
  void*         pAux   // Pointer to auxiliary data used within f
);


typedef void (*Method)(
  int&       S,     
  int&       P,     
  Mat&       A,    
  Mat&       B,   
  Mat&       C
);




#endif