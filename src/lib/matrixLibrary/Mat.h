#ifndef _MAT_H_
#define _MAT_H_

// all the system #include's we'll ever need
#include <iostream>



using namespace std;


class Mat {

public:


	// Cosntructors
	Mat();
	Mat(int n);		        // Zero-based array
	Mat(int m, int n);		// Zero-based array

	Mat(int m, const double &a);	        //initialize to constant value

	Mat(int m, int n , const double &a);	//initialize to constant value

	Mat(int m, const double *a);        	// Initialize to array
	Mat(int m, const int *a);        	// Initialize to array
	Mat(int m, int n , const double *a);	// Initialize to array

	Mat(double x, double y, double z);   // 3dim-Vector
	Mat (double x, double y, double z,   // 6dim-Vector
                double X, double Y, double Z);
    // Assignement
	Mat(const Mat &rhs);	// Copy constructor

	Mat( double Inital, double Final, int N);

	Mat operator=(const Mat& rhs);	//assignment
	
    Mat(const Mat &A, const Mat &B);   // Assemble two matrices

	
    friend Mat Eye(const int n);      // Identity Matrix

	typedef double value_type; // make T available externally

	// Accessing Elements
	//inline double & operator()(const int i);	//i'th element
   // inline double & operator()(const int i, const int j);	//i'th and j'th element
	 double  operator () (int i) const { return v_[i]; };
     double& operator () (int i)       { return v_[i]; };

	 double  operator () (int i,int j) const { return v_[i*n_cols_+j]; };
     double& operator () (int i,int j)       { return v_[i*n_cols_+j]; };
//	 double & Mat::double()(const int i, const int j)
	//inline const double & operator()(const int i) const;1	
	//inline const double & operator()(const int i , const int j) const;
	 

	 Mat operator ~ ();

	 Mat Transpose();

	// Scalar multiplication and division of a vector
	friend Mat operator * (double value, const Mat& Mat);  
	friend Mat operator * (const Mat& V, double value);
	friend Mat operator / (const Mat& V, double value);

	
	

	// Matrix addition and subtraction
    friend Mat operator + (const Mat& left, const Mat& right);
    friend Mat operator - (const Mat& left, const Mat& right); 

	 // Matrix product
    friend Mat operator * (const Mat& left, const Mat& right);
   
	// Unary negative
	friend Mat operator - (const Mat& V);

	// Mat & operator*(const double & value);
	//friend Mat operator*(const double& dFactor, const Mat& oMatrix) ;
	
	// Size
	inline int Size() const;
	//void resize(int newn); // resize (contents not preserved)
	//void assign(int newn, double &a); // resize and assign a constant value
	//void assign(int newn, double &a);

	// Minimum and Maximum Values
	double Min();
	double Max();

	// // Dot product, norm, cross product
	double Dot(const Mat &rhs);
	double Norm();
	Mat    Abs();
	double NormInf();
	double Rms();
	Mat Cross(const Mat &rhs) const;

	// Getting set of elements of the matrix and Assmpling vector in a vector
	Mat Get(const int min,const int max) const;
	Mat Get(const int min_m,const int max_m,const int min_n,const int max_n) const;
	//void Mat::Assemble(const Mat &A, const Mat &B);

	Mat GetCol(int iCol);
	Mat GetRow(int iRow);

	void SetCol(int iCol, Mat Col);
	void SetRow(int iRow,Mat Row);
	// Concatenation 
	friend Mat Assemble(const Mat& a, const Mat& b);


	friend Mat Get(const Mat &A, const double min , const double max);
	// Elementary rotations
    friend Mat R_x(double Angle);
    friend Mat R_y(double Angle);
    friend Mat R_z(double Angle);

	// Sorting
	
	void Sort();
	void Sort(Mat & rhs);
	//Get floati

	float* Getf();
	 // Output
    friend std::ostream& operator << (std::ostream& os, const Mat& Matrix);

	 // Solve Linear System Using Gauss
     
	 
     Mat GaussLin(const Mat rhs,int n) const;

	 Mat Inverse() const;
	// Direct Product 

     Mat DirectProduct(const Mat &rhs) const;

	 // resize
	 void Resize(int newn);

	 void Resize(int nColsNew ,int nRowsNew);
	// Destructor
	~Mat();



	int n_rows_;	// size of array. upper index is nn-1
	int n_cols_;
	double *v_;
};


void vander(Mat x, Mat &w, Mat q) ;


#endif