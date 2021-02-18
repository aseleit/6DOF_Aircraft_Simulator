#include <cmath>
#include <iostream>
#include <iomanip>
#include "mat.h"

Mat::Mat() : n_rows_(0),n_cols_(0), v_(NULL) {}

Mat::Mat(int m) 
	: n_rows_(m), n_cols_(1)
{
  v_ = new double [n_rows_];
  for (int i=0; i<m; i++) v_[i]=0.0;
}

Mat::Mat(int m, int n) : n_rows_(m), n_cols_(n), v_(m*n>0 ? new double[m*n] : NULL) {}

Mat::Mat(int m, const double& a) : n_rows_(m), n_cols_(1),v_(m>0 ? new double[m] : NULL)
{
	for(int i=0; i<m; i++) v_[i] = a;
}

Mat::Mat(int m, int n, const double& a) : n_rows_(m), n_cols_(n),v_(m*n>0 ? new double[m*n] : NULL)
{
	for(int i=0; i<m*n; i++) v_[i] = a;
}

Mat::Mat(int m, const double *a) : n_rows_(m), n_cols_(1), v_(m>0 ? new double[m] : NULL)
{
	for(int i=0; i<m; i++) v_[i] = *a++;
}

Mat::Mat(int m, const int *a) : n_rows_(m), n_cols_(1), v_(m>0 ? new double[m] : NULL)
{
	for(int i=0; i<m; i++) v_[i] = *a++;
}

Mat::Mat(int m, int n, const double *a) : n_rows_(m), n_cols_(n), v_(m*n>0 ? new double[m*n] : NULL)
{
	for(int i=0; i<m*n; i++) v_[i] = *a++;
}

Mat::Mat (double Initial , double Final , int N)
{
	v_= new double [N];
	n_rows_ = N; n_cols_ = 1;
	double h = (Final - Initial)/(N-1);	
	for ( int i=0; i<N; i++)
		v_[i] = Initial + i*h;
}

Mat::Mat (const double x, const double y, const double z)   // 3dim-v_ector             
{
	v_ = new double [3];
	n_rows_=3; n_cols_=1;
	v_[0]=x; v_[1]=y; v_[2]=z;  
}

Mat::Mat (double x, double y, double z,   // 6dim-v_ector
                double X, double Y, double Z)  
{
	v_= new double [6];
	n_rows_=6; n_cols_=1;
	v_[0]=x; v_[1]=y; v_[2]=z;
	v_[3]=X; v_[4]=Y; v_[5]=Z;  
}

Mat::Mat(const Mat &rhs) : n_rows_(rhs.n_rows_),n_cols_(rhs.n_cols_), v_(n_rows_>0 ? new double[n_rows_*n_cols_] : NULL)
{
	for(int i=0; i<n_rows_*n_cols_; i++) v_[i]=rhs.v_[i];
}

Mat::Mat(const Mat &A, const Mat &B)
{
	n_rows_=A.n_rows_+B.n_rows_;
	n_cols_=A.n_cols_;
	v_ = new double [n_rows_];
	int i;
	for (i=0; i<A.n_rows_; i++) v_[i]=A.v_[i];
	for (i=0; i<B.n_rows_; i++) (*this)(i+B.n_rows_)=B.v_[i];
}

Mat Mat::operator=(const Mat& rhs)
{
	delete v_;	
	n_rows_ = rhs.n_rows_; n_cols_=rhs.n_cols_;
	int nEl = n_rows_*n_cols_;
	v_ = new double [nEl];
	for (int i=0; i<nEl; i++)
			v_[i]=rhs.v_[i];
	return *this;
}

Mat Eye(const int n)
{
	Mat eyeMat(n,n);
	for (int j=0; j<n*n; j++)
		eyeMat.v_[j] = 0.0;
	for (int i=0; i<n; i++)
		eyeMat.v_[i*n + i] = 1.0;
	return eyeMat;
}

Mat Mat::operator~()
{
	Mat transpose(n_cols_,n_rows_);
	for (int i=0; i<n_rows_; i++)
		for(int j=0; j<n_cols_; j++)
			transpose.v_[j*n_cols_+i] = v_[i*n_cols_+j];
	return transpose;
}

Mat Mat::Transpose()
{
	Mat trans(n_cols_,n_rows_);

	for (int i=0; i<n_rows_; i++)
		for(int j=0; j<n_cols_; j++)
			trans.v_[j*n_cols_+i] = v_[i*n_cols_+j];	
	return trans;
}

Mat operator*(double value, const Mat& Matrix)
{
	int nEl= Matrix.n_rows_*Matrix.n_cols_;
	Mat Aux(Matrix.n_rows_,Matrix.n_cols_);
	for (int i=0; i<nEl; i++) 	
		Aux.v_[i]=value*Matrix.v_[i];
	return Aux;
}

Mat operator*(const Mat& Matrix, double value)
{
	return value*Matrix;
}

Mat operator/(const Mat& Matrix,double value)
{
	int nEl= Matrix.n_rows_*Matrix.n_cols_;
	Mat Aux(Matrix.n_rows_,Matrix.n_cols_);
	for (int i=0; i<nEl; i++) 	
		Aux.v_[i]=Matrix.v_[i]/value;

	return Aux;
}

Mat operator + (const Mat& left, const Mat& right)
{
	if ( (left.n_rows_!=right.n_rows_) || (left.n_cols_!=right.n_cols_) ) {
		cerr << "ERROR: Incompatible shape in +(Matrix,Matrix)" << endl;
		exit(1);
	};
	Mat Aux(left.n_rows_,left.n_cols_);
	for (int i=0; i<left.n_rows_*left.n_cols_; i++) 	
		Aux.v_[i] = left.v_[i] + right.v_[i];
	return Aux;
}

Mat operator - (const Mat& left, const Mat& right)
{
	if ( (left.n_rows_!=right.n_rows_) || (left.n_cols_!=right.n_cols_) ) {
		cerr << "ERROR: Incompatible shape in +(Matrix,Matrix)" << endl;
		exit(1);
	};
	Mat Aux(left.n_rows_,left.n_cols_);
	for (int i=0; i<left.n_rows_*left.n_cols_; i++) 	
		Aux.v_[i] = left.v_[i] - right.v_[i];
	return Aux;
}

Mat operator * (const Mat& left, const Mat& right)
{
	Mat Aux(left.n_rows_,right.n_cols_);
	double Sum;
	for (int i=0; i<left.n_rows_; i++) 
		for (int j=0; j<right.n_cols_; j++)
		{
			Sum = 0.0;
			for (int k=0; k<left.n_cols_; k++) 
			Sum += left.v_[i*left.n_cols_+k] * right.v_[k*right.n_cols_+j];
			Aux.v_[i*Aux.n_cols_+j] = Sum;
		}
	return Aux;
}

Mat operator - (const Mat& M)
{
	Mat Aux(M.n_rows_,M.n_cols_);
	for (int i=0; i<M.n_rows_*M.n_cols_; i++)  
		Aux.v_[i]=-M.v_[i];		
	return Aux;
}

inline int Mat::Size() const
{
	return n_rows_*n_cols_;
}

double Mat::Min() 
{
	double minValue=v_[0];
	for( int i=0 ; i<Size() ; i++)
		minValue = v_[i] < minValue? v_[i] : minValue;
	return minValue;	
}

double Mat::Max() 
{
	double maxValue=v_[0];
	for( int i=0 ; i<Size() ; i++)
		maxValue = v_[i]> maxValue? v_[i] : maxValue;
	return maxValue;
}

double Mat::Dot(const Mat &rhs)
{
	if ((*this).Size()!=rhs.Size()) {
		cerr << "ERROR: Incompatible shape in Dot(Vector,Vector)" << endl;
		exit(1);
	};
	double Sum = 0.0;
	for (int i=0; i<(*this).Size(); i++) Sum+=(*this).v_[i]*rhs.v_[i];
	return Sum;
}

double Mat::Norm()
{
	return sqrt((*this).Dot((*this)));
}

Mat Mat::Abs() 
{
	Mat absMatrix(n_rows_ , n_cols_);
	for (int i=0; i<Size(); i++)
		absMatrix.v_[i]=abs(v_[i]);
	return absMatrix;
}
double Mat::NormInf()
{	
	double output=abs(v_[0]);	
	for (int i = 1; i < Size(); i++)
		output = output<abs(v_[i]) ? output=abs (v_[i]):output;
	return output;
}

double Mat::Rms(){
	double Result;
	double Diff=0.0;
	int N=n_rows_*n_cols_;
	for (int i=0; i<N; i++) Diff=Diff+pow((*this)(i),2);
	Result = pow(Diff/(N),0.5);
	return Result;
}

Mat Mat::Cross(const Mat &rhs) const
{
	if ( ((*this).Size()!=3 || (rhs.Size()!=3) )) {
	cerr << "ERROR: Invalid dimension in Cross(Vector,Vector)" << endl;
	exit(1);
	};
	Mat Result(3);
	Result.v_[0] = (*this).v_[1]*rhs.v_[2] - (*this).v_[2]*rhs.v_[1];
	Result.v_[1] = (*this).v_[2]*rhs.v_[0] - (*this).v_[0]*rhs.v_[2];
	Result.v_[2] = (*this).v_[0]*rhs.v_[1] - (*this).v_[1]*rhs.v_[0];
	return Result;
}

Mat Mat::Get(const int min,const int max) const
{
	Mat getmat(max-min+1);
	int im;
	for(int i=0, im=min; im<max+1; ++i , ++im)
	{		
		getmat.v_[i]=v_[im]; 	
	}
	return getmat;
}

Mat Mat::Get(const int min_m,const int max_m,const int min_n,const int max_n) const
{
	Mat getmat(max_m-min_m+1,max_n-min_n+1);
	int i, im , in;  
	for(i = 0 , im=min_m; im<max_m+1;  ++im)
		for(in=min_n; in<max_n+1; ++i, ++in)		
			getmat.v_[i]=v_[im*n_cols_+in]; 	
	return getmat;
}

Mat Get(const Mat &A, const double min , const double max)
{
	Mat Aux(max-min+1);
	int i,j;
	for (i=min , j=0 ; i<max+1; i++ , j++)
		Aux.v_[j]=A.v_[i];
	return Aux;
}

Mat Mat::GetCol(int iCol)
{
	Mat Col(n_rows_);
	for ( int i = 0;  i< n_rows_; i++)
		Col.v_[i] = v_[n_cols_ * i+iCol];
	return Col;
}

Mat Mat::GetRow(int iRow)
{

	Mat Row(n_cols_);
	for ( int i = 0; i < n_cols_; i++)
		Row.v_[i] = v_[n_cols_ * iRow+i];
	return Row;
}

void Mat::SetCol(int iCol, Mat Col)
{
	for ( int i = 0; i < n_rows_; i++)
		v_[n_cols_ * i+iCol] = Col.v_[i];
}

void Mat::SetRow(int iRow, Mat Row)
{
	for ( int i = 0; i < n_cols_; i++)
		v_[n_cols_ * iRow+i] = Row.v_[i];
}

Mat R_x(double Angle)
{
	const double C = cos(Angle);
	const double S = sin(Angle);
	Mat U(3,3);
	U.v_[0] = 1.0;  U.v_[1] = 0.0;  U.v_[2] = 0.0;
	U.v_[3] = 0.0;  U.v_[4] =  +C;  U.v_[5] =  +S;
	U.v_[6] = 0.0;  U.v_[7] =  -S;  U.v_[8] =  +C;
	return U;
}

Mat R_y(double Angle)
{
	const double C = cos(Angle);
	const double S = sin(Angle);
	Mat U(3,3);
	U.v_[0] =  +C;  U.v_[1] = 0.0;  U.v_[2] =  -S;
	U.v_[3] = 0.0;  U.v_[4] = 1.0;  U.v_[5] = 0.0;
	U.v_[6] =  +S;  U.v_[7] = 0.0;  U.v_[8] =  +C;
	return U;
}

Mat R_z(double Angle)
{
	const double C = cos(Angle);
	const double S = sin(Angle);
	Mat U(3,3);
	U.v_[0] =  +C;  U.v_[1] =  +S;  U.v_[2] = 0.0;
	U.v_[3] =  -S;  U.v_[4] =  +C;  U.v_[5] = 0.0;
	U.v_[6] = 0.0;  U.v_[7] = 0.0;  U.v_[8] = 1.0;
	return U;
}

void Mat::Sort()
{
	int n = n_rows_*n_cols_;
	double temp;
	for(int i=0; i<n; i++)
	{
		for(int j=i ; j<n ; j++ )
		{
			if(v_[j] < v_[i])
			{
				temp = v_[j] ;
				v_[j] = v_[i];
				v_[i] = temp;
			}			
		}
	}
}

void Mat::Sort(Mat &rhs)
{
	int n = n_rows_*n_cols_;
	double temp ,tempRhs;
	for(int i=0; i<n; i++)
	{
		for(int j=i; j<n; j++ )
		{
			if(v_[j] < v_[i])
			{
				temp = v_[j] ;
				v_[j] = v_[i];
				v_[i] = temp;
				tempRhs = rhs.v_[j] ;
				rhs.v_[j] = rhs.v_[i];
				rhs.v_[i] = tempRhs;
			}			
		}
	}
}

float* Mat::Getf() {
	float * fdata;
	int size =n_rows_*n_cols_;
	fdata=new float[size];
	for(int i=0; i<size; i++) 
		fdata[i]=(float) v_[i];
	return fdata;
}

ostream& operator << (ostream& os, const Mat& Matrix)
{
	int w = os.width();
	os.precision(3);
	for (int i=0; i<Matrix.n_rows_; i++) {
	for (int j=0; j<Matrix.n_cols_; j++)
		os <<"   "<<  Matrix(i,j);
	os << endl<<endl;
	}
	return os;
}

Mat Mat::GaussLin(const Mat rhs,int n) const
{
	Mat Aux = (*this);
	if ( (Aux.n_rows_!=rhs.Size() )) {
    cerr << "ERROR: Invalid dimension in gaussLin" << endl;
    exit(1);
	};
	double m;
	for (int k=0; k<n-1; k++ )
	{
		int l=k;
		// Rearrangement
		for (int i=k+1; i<n; i++)
		{
			if (abs(Aux(i,k)) > abs(Aux(l,k))) 
				l=i;
		}
		if (l!=k)
		{
			for (int j=k; j<n; j++)
			{
				m = Aux(k,j);
				Aux.v_[k*n_cols_+j] = Aux.v_[l*n_cols_+j];
				Aux.v_[l*n_cols_+j] = m;
			}
			m = rhs.v_[k];
			rhs.v_[k] = rhs.v_[l];
		    rhs.v_[l] = m;
		}
		for (int i = k+1; i<n; i++)
		{
			m = Aux.v_[i*n_cols_+k] / Aux.v_[k*n_cols_+k];
			Aux.v_[i*n_cols_+k] = 0.0;
			for ( int j = k+1; j < n; j++)
			{
				Aux.v_[i*n_cols_ + j] = Aux.v_[i*n_cols_ + j] - m * Aux.v_[k*n_cols_+j];
			}
			rhs.v_[i] = rhs.v_[i] - m * rhs.v_[k] ; 
		}
	}
	Mat X(n);
	X.v_[n-1] = rhs.v_[n-1] / Aux.v_[(n-1)*n_cols_+n-1];
	for (int i = n-2 ;i >-1 ;i--)
	{
		double s = 0.0;		
		for (int j=i+1 ; j<n; j++)
			s = s - Aux.v_[i*n_cols_ + j] * X.v_[j];
		X.v_[i] = (rhs.v_[i] + s)/ Aux.v_[i*n_cols_ + i];
	}
	~Aux;
	return X;
}

Mat Mat::Inverse() const
{
	int i,j;
	int n = n_rows_;
	Mat ainv(n,n);
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) ainv(i,j) = 0.;
		ainv(i,i) = 1.;
	}
	int m=n_cols_;	
	Mat xx(n),yy(n);
	for (j=0; j<m; j++) {
		for (i=0; i<n; i++) xx(i) = ainv(i,j);
		yy=GaussLin(xx,n);
		for (i=0; i<n; i++) ainv(i,j) = yy(i);
	}
	return ainv;
}

Mat Mat::DirectProduct(const Mat &rhs) const
{
	Mat product(n_rows_*rhs.n_rows_,n_cols_*rhs.n_cols_);
	for (int im=0; im<n_rows_ ; im++ )
	{
		for (int jm=0; jm<n_cols_; jm++)
		{
			for (int in=0; in<rhs.n_rows_ ; in++ )
			{
				for (int jn=0; jn<rhs.n_cols_; jn++)
				{
					int m = im*rhs.n_rows_ +in;
					int n = jm*rhs.n_cols_+ jn;
					product( m , n) = v_[im*n_cols_+jm] *  rhs.v_[in*rhs.n_cols_+jn];
				}
			}
		} 
	}
	return product;
}

void Mat::Resize(int newn)
{
	Mat ytemp = (*this);
	if (newn != n_cols_*n_rows_) {
		if (v_ != NULL) delete[] (v_);	
		n_rows_ = newn; n_cols_=1;
		v_ = n_rows_*n_cols_ > 0 ? new double[n_rows_] : NULL;		
		for ( int i = 0; i < min(n_rows_,ytemp.n_rows_); i++)
			v_[i] = ytemp.v_[i];
	}
}

void Mat::Resize(int newn, int newm)
{
	int i,nel;
	Mat ytemp = (*this);

	if (newn != n_cols_ || newm != n_rows_) {
		if (v_ != NULL) 
			delete[] (v_);		
		n_rows_ = newn;
		n_cols_ = newm;
		int nEl = n_rows_*n_cols_;
		v_ = n_rows_*n_cols_>0 ? new double[nEl] : NULL;		
		for ( int i = 0; i < min(n_rows_,ytemp.n_rows_) ; i++)
			for ( int j = 0; j < min(n_cols_,ytemp.n_cols_) ; j++)
				v_[i*n_cols_+j] = ytemp.v_[i*ytemp.n_cols_+j];
	}
}

Mat::~Mat()
{
	if (v_ != NULL) delete[] (v_);
};

