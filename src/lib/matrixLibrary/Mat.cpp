
#include <cmath>
#include <iostream>
#include <iomanip>


#include "Mat.h"
// Mat definitions



Mat::Mat() : nRows(0),nCols(0), v(NULL) {}


Mat::Mat(int m) 
	: nRows(m), nCols(1)
{
  v = new double [nRows];
  for (int i=0; i<m; i++) v[i]=0.0;
}


Mat::Mat(int m, int n) : nRows(m), nCols(n), v(m*n>0 ? new double[m*n] : NULL) {}


Mat::Mat(int m, const double& a) : nRows(m), nCols(1),v(m>0 ? new double[m] : NULL)
{
	for(int i=0; i<m; i++) v[i] = a;
}


Mat::Mat(int m, int n, const double& a) : nRows(m), nCols(n),v(m*n>0 ? new double[m*n] : NULL)
{
	for(int i=0; i<m*n; i++) v[i] = a;
}


Mat::Mat(int m, const double *a) : nRows(m), nCols(1), v(m>0 ? new double[m] : NULL)
{
	for(int i=0; i<m; i++) v[i] = *a++;
}

Mat::Mat(int m, const int *a) : nRows(m), nCols(1), v(m>0 ? new double[m] : NULL)
{
	for(int i=0; i<m; i++) v[i] = *a++;
}

Mat::Mat(int m, int n, const double *a) : nRows(m), nCols(n), v(m*n>0 ? new double[m*n] : NULL)
{

	for(int i=0; i<m*n; i++) v[i] = *a++;
}

Mat::Mat (double Initial , double Final , int N)
{
	v= new double [N];
	nRows = N; nCols = 1;
	double h = (Final - Initial)/(N-1);
	
	for ( int i = 0; i<N; i++)
		v[i] = Initial + i*h;

//	return Linspace;
}
Mat::Mat (const double x, const double y, const double z)   // 3dim-Vector             
{
	//init(6,1);
	v = new double [3];

	nRows=3; nCols=1;
	v[0]=x; v[1]=y; v[2]=z;  
}


Mat::Mat (double x, double y, double z,   // 6dim-Vector
                double X, double Y, double Z)  
{
	v= new double [6];
	nRows=6; nCols=1;
	v[0]=x; v[1]=y; v[2]=z;
	v[3]=X; v[4]=Y; v[5]=Z;  
}


Mat::Mat(const Mat &rhs) : nRows(rhs.nRows),nCols(rhs.nCols), v(nRows>0 ? new double[nRows*nCols] : NULL)
{
	for(int i=0; i<nRows*nCols; i++) v[i] = rhs.v[i];
}



Mat::Mat(const Mat &A, const Mat &B)
{
	nRows=A.nRows+B.nRows;
	nCols=A.nCols;

	v = new double [nRows];
	int i;
	for (i=0;i<A.nRows;i++) v[i]=A.v[i];
	for (i=0;i<B.nRows;i++) (*this)(i+B.nRows)=B.v[i];
}

Mat Mat::operator=(const Mat& rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	/*int nEl = rhs.nRows*rhs.nCols;
	if (this != &rhs)
	{
		if (nRows != rhs.nRows || nCols != rhs.nCols) {
			if (v != NULL) delete [] (v);
			nRows = rhs.nRows;
			nCols = rhs.nCols;
			v= nEl>0 ? new double[nEl] : NULL;
		}
		for (int i=0; i<nEl; i++)
			v[i]=rhs.v[i];
	}
	return *this;
	*/
	delete v;
	
	nRows = rhs.nRows; nCols=rhs.nCols;
	int nEl = nRows*nCols;
	v = new double [nEl];
	for (int i=0; i<nEl; i++)
			v[i]=rhs.v[i];

	return *this;

}

Mat eye(const int n)
{
	Mat eyeMat(n,n);

	for (int j = 0 ; j < n*n ; j++)
		eyeMat.v[j] = 0.0;

	for (int i = 0 ; i < n ; i++)
		eyeMat.v[i*n + i] = 1.0;

	return eyeMat;

}

//Mat::Mat(const int nRows_, double& a): nRows(nRows_) , nCols(1) ,  v(nRows>0 ? new double[nRows] : NULL)


/*
inline double & Mat::operator()(const int i)	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nRows) {
	throw("Mat subscript out of bounds");
}
#endif
	return v[i];
}


inline double& Mat::operator()(const int i, const int j)	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nRows || j<0 || j>=nCols) {
	throw("Mat subscript out of bounds");
}
#endif
	return v[i*nCols+j];
}



inline const double & Mat::operator()(const int i) const	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nRows) {
	throw("Mat subscript out of bounds");
}
#endif
	return v[i];
}


inline const double& Mat::operator()(const int i, const int j) const	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nRows || j<0 || j>=nCols) {
	throw("Mat subscript out of bounds");
}
#endif
	return v[i*nCols+j];
}
*/

// Scalar multiplication and division of a vector
/*
Matrix operator * (double value, const Matrix& Mat)
{
  Matrix Aux(Mat.n,Mat.m);
  for (int i=0; i<Mat.n; i++) 
    for (int j=0; j<Mat.m; j++) 
      Aux.M[i][j]=value*Mat.M[i][j];
  return Aux;
}
*/

Mat Mat::operator~()
{
	Mat transpose(nCols,nRows);

	for (int i=0; i<nRows; i++)
	for(int j=0; j<nCols; j++)
	transpose.v[j*nCols+i] = v[i*nCols+j];

	return transpose;


}

Mat Mat::transpose()
{
	Mat trans(nCols,nRows);

	for (int i=0; i<nRows; i++)
	for(int j=0; j<nCols; j++)
	trans.v[j*nCols+i] = v[i*nCols+j];
	
	return trans;
}

Mat operator*(double value, const Mat& Matrix)
{
	int nEl= Matrix.nRows*Matrix.nCols;
	Mat Aux(Matrix.nRows,Matrix.nCols);
	for (int i=0; i<nEl; i++) 	
		Aux.v[i]=value*Matrix.v[i];

	return Aux;
}

Mat operator*(const Mat& Matrix, double value)
{
	return value*Matrix;
}

Mat operator/(const Mat& Matrix,double value)
{
	int nEl= Matrix.nRows*Matrix.nCols;
	Mat Aux(Matrix.nRows,Matrix.nCols);
	for (int i=0; i<nEl; i++) 	
		Aux.v[i]=Matrix.v[i]/value;

	return Aux;
}



// Matrix addition and subtraction

Mat operator + (const Mat& left, const Mat& right)
{
	if ( (left.nRows!=right.nRows) || (left.nCols!=right.nCols) ) {
	cerr << "ERROR: Incompatible shape in +(Matrix,Matrix)" << endl;
	exit(1);
  };

	Mat Aux(left.nRows,left.nCols);
	for (int i=0; i<left.nRows*left.nCols; i++) 	
		Aux.v[i] = left.v[i] + right.v[i];

	return Aux;
}

Mat operator - (const Mat& left, const Mat& right)
{
	if ( (left.nRows!=right.nRows) || (left.nCols!=right.nCols) ) {
	cerr << "ERROR: Incompatible shape in +(Matrix,Matrix)" << endl;
	exit(1);
  };

	Mat Aux(left.nRows,left.nCols);
	for (int i=0; i<left.nRows*left.nCols; i++) 	
		Aux.v[i] = left.v[i] - right.v[i];
	return Aux;
}


 // Matrix product
Mat operator * (const Mat& left, const Mat& right)
{
/*
  if (left.m!=right.n) {
    cerr << "ERROR: Incompatible shape in *(Matrix,Matrix)" << endl;
    exit(1);
  };
  */
  Mat Aux(left.nRows,right.nCols);
  double Sum;
  for (int i=0; i<left.nRows; i++) 
    for (int j=0; j<right.nCols; j++)
	{
      Sum = 0.0;
      for (int k=0; k<left.nCols; k++) 
        Sum += left.v[i*left.nCols+k] * right.v[k*right.nCols+j];
      Aux.v[i*Aux.nCols+j] = Sum;
    }
  return Aux;

}

// Unary minus
Mat operator - (const Mat& M)
{
  Mat Aux(M.nRows,M.nCols);
  for (int i=0; i<M.nRows*M.nCols; i++)  Aux.v[i]=-M.v[i];
     
  return Aux;
}

/*
 
Mat operator * (const Mat& Mat, double value)
{
  return value*Mat;
}

 
Mat operator / (const Mat& Mat, double value)
{
  Matrix Aux(Mat.n,Mat.m);
  for (int i=0; i<Mat.n; i++) 
    for (int j=0; j<Mat.m; j++) 
      Aux.M[i][j]=Mat.M[i][j]/value;
  return Aux;
}
*/


inline int Mat::size() const
{
	return nRows*nCols;
}


double Mat::Min() 
{
	double minValue=v[0];
	for( int i=0 ; i<size() ; i++)
		minValue = v[i] < minValue? v[i] : minValue;
	return minValue;	
}

double Mat::Max() 
{
	double maxValue=v[0];
	for( int i=0 ; i<size() ; i++)
		maxValue = v[i]> maxValue? v[i] : maxValue;
	return maxValue;
}
 /*
void Mat::resize(int newn)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
}

*/



/*
void Mat::assign(int newn, const T& a)
{
	if (newn != nRows) {
		if (v != NULL) delete[] (v);
		nRows = newn;
		v = nRows > 0 ? new T[nRows] : NULL;
	}
	for (int i=0;i<nRows;i++) v[i] = a;
}
*/
double Mat::Dot(const Mat &rhs)
{

	if ((*this).size()!=rhs.size()) {
	cerr << "ERROR: Incompatible shape in Dot(Vector,Vector)" << endl;
	exit(1);
	};
	double Sum = 0.0;
	for (int i=0; i<(*this).size(); i++) Sum+=(*this).v[i]*rhs.v[i];
	return Sum;

}


double Mat::Norm()
{
	return sqrt((*this).Dot((*this)));
}



Mat Mat::Abs() 
{
	Mat absMatrix(nRows , nCols);

	for (int i = 0 ; i < size() ; i++)
		absMatrix.v[i] = abs(v[i]);

	return absMatrix;

}
double Mat::NormInf()
{	double output=abs(v[0]);
	
	for (int i = 1 ; i < size() ; i++)
		output = output < abs(v[i]) ? output =  abs (v[i]) : output;

	return output;
}

double Mat::Rms(){
	double Result;
	double Diff=0.0;
	int N=nRows*nCols;
	for (int i=0;i<N;i++) Diff=Diff+pow((*this)(i),2);
	Result = pow(Diff/(N),0.5);
	return Result;
}


Mat Mat::Cross(const Mat &rhs) const
{
	if ( ((*this).size()!=3 || (rhs.size()!=3) )) {
    cerr << "ERROR: Invalid dimension in Cross(Vector,Vector)" << endl;
    exit(1);
	};
	Mat Result(3);
	Result.v[0] = (*this).v[1]*rhs.v[2] - (*this).v[2]*rhs.v[1];
	Result.v[1] = (*this).v[2]*rhs.v[0] - (*this).v[0]*rhs.v[2];
	Result.v[2] = (*this).v[0]*rhs.v[1] - (*this).v[1]*rhs.v[0];
	return Result;
}

Mat Mat::get(const int min,const int max) const
{
	Mat getmat(max-min+1);
  
	int im;
	for(int i=0, im=min ; im<max+1 ; ++i , ++im)
	{		
	  getmat.v[i]=v[im]; 	
	}
	return getmat;
}

Mat Mat::get(const int min_m,const int max_m,const int min_n,const int max_n) const
{
	Mat getmat(max_m-min_m+1,max_n-min_n+1);

	int i, im , in;
  
	for(i = 0 , im=min_m ; im<max_m+1 ;  ++im)
		for(in=min_n ; in<max_n+1 ; ++i , ++in)		
			getmat.v[i]=v[im*nCols+in]; 	
	
	return getmat;
}

//void Mat::Assemble(const Mat &A, const Mat &B)
//{
//	nRows=A.nRows+B.nRows;
//	nCols=A.nCols;

//	delete v;
//	v = new double [nRows];
//	int i;
//	for (i=0;i<A.nRows;i++) v[i]=A.v[i];
//	for (i=0;i<B.nRows;i++) v[i+A.nRows]=B.v[i];

//}

//Mat Assemble(const Mat &A, const Mat &B)
//{
//	int nRows=A.nRows+B.nRows;
//	int nCols=A.nCols;

//	Mat Aux(nRows,nCols);

//	int i;
//	for (i=0;i<A.nRows;i++) Aux.v[i]=A.v[i];
//	for (i=0;i<B.nRows;i++) Aux.v[i+A.nRows]=B.v[i];

//	return Aux;
//}

Mat Get(const Mat &A, const double min , const double max)
{
	
	Mat Aux(max-min+1);

	int i,j;
	for (i=min , j=0 ;i<max+1;i++ , j++)
		Aux.v[j]=A.v[i];
	

	return Aux;
}


Mat Mat::getCol(int iCol)
{

	Mat Col(nRows);

	//for ( int i = 0 ; i < nRows ; i++)
	//	Col.v[i] = v[nRows * i+iCol];
	for ( int i = 0 ;  i< nRows ; i++)
		Col.v[i] = v[nCols * i+iCol];

	return Col;
}

Mat Mat::getRow(int iRow)
{

	Mat Row(nCols);

	for ( int i = 0 ; i < nCols ; i++)
		Row.v[i] = v[nCols * iRow+i];

	return Row;

}

void Mat::setCol(int iCol, Mat Col)
{

	for ( int i = 0 ; i < nRows ; i++)
		 v[nCols * i+iCol] = Col.v[i];
}

void Mat::setRow(int iRow, Mat Row)
{
	for ( int i = 0 ; i < nCols ; i++)
		v[nCols * iRow+i] = Row.v[i];
}


// Elementary rotations

Mat R_x(double Angle)
{
  const double C = cos(Angle);
  const double S = sin(Angle);
  Mat U(3,3);
  U.v[0] = 1.0;  U.v[1] = 0.0;  U.v[2] = 0.0;
  U.v[3] = 0.0;  U.v[4] =  +C;  U.v[5] =  +S;
  U.v[6] = 0.0;  U.v[7] =  -S;  U.v[8] =  +C;
  return U;
}

Mat R_y(double Angle)
{
  const double C = cos(Angle);
  const double S = sin(Angle);
  Mat U(3,3);
  U.v[0] =  +C;  U.v[1] = 0.0;  U.v[2] =  -S;
  U.v[3] = 0.0;  U.v[4] = 1.0;  U.v[5] = 0.0;
  U.v[6] =  +S;  U.v[7] = 0.0;  U.v[8] =  +C;
  return U;
}

Mat R_z(double Angle)
{
  const double C = cos(Angle);
  const double S = sin(Angle);
  Mat U(3,3);
  U.v[0] =  +C;  U.v[1] =  +S;  U.v[2] = 0.0;
  U.v[3] =  -S;  U.v[4] =  +C;  U.v[5] = 0.0;
  U.v[6] = 0.0;  U.v[7] = 0.0;  U.v[8] = 1.0;
  return U;
}


void Mat::Sort()
{
	int n = nRows*nCols;
	double temp;
	for(int i=0 ; i<n ; i++)
	{
		for(int j=i ; j<n ; j++ )
		{
			if(v[j] < v[i])
			{
				temp = v[j] ;
				v[j] = v[i];
				v[i] = temp;
			}			
		}
	}
}

void Mat::Sort(Mat &rhs)
{
	int n = nRows*nCols;
	double temp ,tempRhs;
	for(int i=0 ; i<n ; i++)
	{
		for(int j=i ; j<n ; j++ )
		{
			if(v[j] < v[i])
			{
				temp = v[j] ;
				v[j] = v[i];
				v[i] = temp;

				tempRhs = rhs.v[j] ;
				rhs.v[j] = rhs.v[i];
				rhs.v[i] = tempRhs;
			}			
		}
	}
}

float* Mat::getf() {
	  //double *fdata;
	float * fdata;
	int size =nRows*nCols;
	 fdata=new float[size];
	 for(int i=0; i<size; i++) 
		fdata[i]=(float)v[i];

	return fdata;
}
//
ostream& operator << (ostream& os, const Mat& Matrix)
{
   int w = os.width();

  os.precision(3);
  for (int i=0; i<Matrix.nRows; i++) {
	for (int j=0; j<Matrix.nCols; j++)
		os <<"   "<<  Matrix(i,j);
    os << endl<<endl;
  }
  return os;
}


Mat Mat::gaussLin(const Mat rhs,int n) const
{
	Mat Aux = (*this);

//	cout << " m A  Inside1 "<<endl <<(*this)<<endl;
	if ( (Aux.nRows!=rhs.size() )) {
    cerr << "ERROR: Invalid dimension in gaussLin" << endl;
    exit(1);
	};

	double m;
	for (int k=0; k<n-1 ;k++ )
	{
		int l=k;
		// Rearrangement

		for (int i=k+1; i<n ; i++)
		{
			if (abs(Aux(i,k)) > abs(Aux(l,k))) 
				l=i;
		}

		if (l!=k)
		{
			for (int j = k ; j < n ; j++)
			{
				m = Aux(k,j);
				Aux.v[k*nCols+j] = Aux.v[l*nCols+j];
				Aux.v[l*nCols+j] = m;
			}
			m = rhs.v[k];
			rhs.v[k] = rhs.v[l];
		    rhs.v[l] = m;
		}
	//	cout << endl << " Ma before Exculsion "<<endl<<(*this)<<endl<< "  rhs  "<<endl <<rhs<<endl; 
		// Exclusion

		for (int i = k+1 ; i<n ; i++)
		{
			m = Aux.v[i*nCols+k] / Aux.v[k*nCols+k];
		//	cout << endl << "   m =  " <<m<<endl;

			Aux.v[i*nCols+k] = 0.0;

			//cout <<endl << " mA  "<<endl<<(*this) <<endl;
			for ( int j = k+1 ; j < n ; j++)
			{
				Aux.v[i*nCols + j] = Aux.v[i*nCols + j] - m * Aux.v[k*nCols+j];

			}

			rhs.v[i] = rhs.v[i] - m * rhs.v[k] ; 
			//cout <<endl << "   rhs.v[i]  "<<rhs.v[i]<<endl;
		}
	}

	// Backward Substitution

	Mat X(n);

//	cout << " m A  after exculsion "<<endl <<(*this)<<endl;
	X.v[n-1] = rhs.v[n-1] / Aux.v[(n-1)*nCols+n-1];

	for (int i = n-2 ; i > -1 ; i--)
	{
		double s = 0.0;
		
		for (int j = i+1 ; j < n ; j++)
			s = s - Aux.v[i*nCols + j] * X.v[j];

		X.v[i] = (rhs.v[i] + s)/ Aux.v[i*nCols + i];
	}
	~Aux;
	return X;
}

	

Mat Mat::inverse() const
{
	int i,j;
	int n = nRows;
	Mat ainv(n,n);
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) ainv(i,j) = 0.;
		ainv(i,i) = 1.;
	}

	int m=nCols;
	
	Mat xx(n),yy(n);
	for (j=0;j<m;j++) {
	//	cout << endl << " this " <<endl<<(*this)<<endl;
		for (i=0;i<n;i++) xx(i) = ainv(i,j);
	//	cout << endl << " this " <<endl<<(*this)<<endl;
	//	cout << endl <<" xx  " <<endl<<xx; 
		yy=gaussLin(xx,n);
	//	cout << endl << " this " <<endl<<(*this)<<endl;
	//	cout << endl <<" yy  " <<endl<<yy; 
	//	cout << endl << " this " <<endl<<(*this)<<endl;
		for (i=0;i<n;i++) ainv(i,j) = yy(i);
	}
	return ainv;

}
Mat Mat::directProduct(const Mat &rhs) const
{
	Mat product(nRows*rhs.nRows,nCols*rhs.nCols);


	
	for (int im=0; im<nRows ;im++ )
	{
		for (int jm=0; jm<nCols ; jm++)
		{
			for (int in=0; in<rhs.nRows ;in++ )
			{
				for (int jn=0; jn<rhs.nCols ; jn++)
				{

					int m = im*rhs.nRows +in;
					int n = jm*rhs.nCols+ jn;

					product( m , n) = v[im*nCols+jm] *  rhs.v[in*rhs.nCols+jn];
				}
			}
		} 
	}
	return product;
}



void vander(Mat x, Mat &w, Mat q) 

{
	int i,j,k,n=q.size();
	double b,s,t,xx;
	Mat c(n);
	if (n == 1) w(0)=q(0);
	else {
		for (i=0;i<n;i++) c(i)=0.0; 
		c(n-1) = -x(0); 
		for (i=1;i<n;i++) 
		{ 
			xx = -x(i);
			for (j=(n-1-i);j<(n-1);j++) c(j) += xx*c(j+1);
			c(n-1) += xx;
		}
		for (i=0;i<n;i++)
		{ 		
			xx=x(i);
			t=b=1.0;
			s=q(n-1);
			for (k=n-1;k>0;k--) 
			{ 
				b=c(k)+xx*b;
				s += q(k-1)*b; 
				t=xx*t+b;
			}
			w(i)=s/t; 
		}
	}
}

void Mat::resize(int newn)
{
	Mat ytemp = (*this);
	if (newn != nCols*nRows) {
		if (v != NULL) delete[] (v);	
		nRows = newn; nCols=1;
		v = nRows*nCols > 0 ? new double[nRows] : NULL;

		
		for ( int i = 0 ; i < min(nRows,ytemp.nRows) ; i++)
			v[i] = ytemp.v[i];
	}
}

void Mat::resize(int newn, int newm)
{
	int i,nel;
	Mat ytemp = (*this);

	if (newn != nCols || newm != nRows) {
		if (v != NULL) 
			delete[] (v);
		
		nRows = newn;
		nCols = newm;
		int nEl = nRows*nCols;
		v = nRows*nCols>0 ? new double[nEl] : NULL;
		
		for ( int i = 0 ; i < min(nRows,ytemp.nRows) ; i++)
			for ( int j = 0 ; j < min(nCols,ytemp.nCols) ; j++)
				v[i*nCols+j] = ytemp.v[i*ytemp.nCols+j];
	}
}

Mat::~Mat()
{
	if (v != NULL) delete[] (v);
};

