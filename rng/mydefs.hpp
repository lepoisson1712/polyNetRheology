/*
 * mydefs.hpp
 *
 *  Created on: Jun 26, 2012
 *      Author: amuasi
 */

/*
 * mydefs.hpp
 *
 *  Created on: Nov 12, 2009
 *      Author: henry
 */

#ifndef MYDEFS_HPP_
#define MYDEFS_HPP_

#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <set>
#include <list>

#define DO_NOT_CHECK_INTERACTIONS0
#define DO_NOT_CHECK_INTERACTIONS
//#define CHECK_ENERGY
#define BOND_CENTRED_ROTATION
//#define CHECK_DELTAE
//#define TEST_TIMING
//#define MY_OWN_RG
//#define TEST_REPETITIVE_BIAS
//#define CHECK_OVERALLOVERLAP
//#define CHECK_NEW_JAC
//#define CHECK_JAC_COMPUTATION
//#define CHECK_TOPO
#define TESTING1
//#define TESTING2
//#define CHECK_RB
//#define CHECK_RG
//#define CHECK_BSR
#define SELECTION_PROBABILITY 1.0
//#define DEBUGGING
//#define DEBUG_JJ 7785
//#define DEBUG_STATEMENT globaljj==7785
// CHOP_TOLERANCE is used in Prismatoid.hpp & IntersectingCylinders.hpp
#define CHOP_TOLERANCE 1e-5
#define GIBBS_ENSEMBLE
#define HELMHOLTZ_ENSEMBLE
#define WRMC
#define SOLVE_PRECISION 1e-4
#define BOND_PRECISION 1e-3
#define ROT_EPSILON 1e-18
#define MY_EPSILON 1e-3
#define MY_EPSILON4CHOP 1e-3
#define kB 1	//1.380650424e-23 	// The real Value
#define TRACTRIX_STEP_LENGTH 0	//Steps for Tractrix in units of bond length
//#define TRACTRIX_FOLD 10
#define SLOW_MOTION
#define CHECK_REVERSIBILITY
//#define FREE_POLYMER
//#define TRACK_ACCEPTANCE_PROBABILITY
//#define TRACK_POLYMER_STRUCTURE
#define USE_NONJACOBIAN_SOLVER
//#define USE_FJC
//#define JUST_R
//#define CHOLESKY_JACOBIAN
//#define LEONTIDIS_JAC
#define MAGGS_JAC
//#define ROBOTICS_JACOBIAN
//#define TRACTRIX_JACOBIAN
//#define TRACTRIX_JACOBIAN2
//#define USE_NewPullOnString2
#define USE_ONLY_TRACTRIX
//#define USE_EULERANGLES
//#define TRIMER_MOVES_ONLY
//#define USE_RAMP
//#define REPORTING_ON
//#define REPORTING_TBM_ON
//#define REPORT_STEPCOUNT
//#define REPORT_JACOBIANS
#define TRACK_END
#define IGNORE_BOND_LENGTH_CHECKS
//#define DIM2 2
#define DIM3 3
//#define TEST4PROPOSALSYMMETRY

#define STREAM_LENGTH 100
#define PI 3.1415926535897932384626433832795028841968
#define Degree PI/180

#define OUTPUT_PRECISION 5

#define JUST_DISPLAY(x) std::cout << x

#define DISPLAY(x) std::cout << x << std::endl

#define DISPLAY2(pre, x, post) std::cout << pre << x << post << std::endl

#define VAR_NAME(x) #x

#define DISPLAY_VAR(x) std::cout << #x " = "<< x << ";" << std::endl

#define DISPLAY_VAR_SC(x) { \
	std::cout << #x " = " << std::scientific << x << ";" << std::endl; \
	std::cout.unsetf(ios::scientific);}

#define DISPLAY_BVAR(x) std::cout << #x " = "<< ((x) ? "true" : "false") << std::endl

#define DISPLAY_BVAL(x) std::cout << ((x) ? "true" : "false")

#define DISPLAY_POOMA(x) std::cout << #x << x << ";" << std::endl

#define DISPLAY_POOMAV(x) std::cout << #x " = " << std::endl << "{" << x(0) << ", " << x(1) << ", " << x(2) << "}" << std::endl

#define DISPLAY_POOMAV2D(x) std::cout << #x " = " << std::endl << "{" << x(0) << ", " << x(1) << "}" << std::endl

#define DISPLAY_POOMAV_WOH(x) std::cout << "{" << x(0) << ", " << x(1) << ", " << x(2) << "}"

#define DISPLAY_POOMA2_wo_HEAD(x){             \
	std::cout << "{"; \
	if (x.length(0) > 0) std::cout << "{" << x(0)(0) << ", " << x(0)(1) << ", " << x(0)(2) << "}"; \
	for (int di = 1; di < x.length(0); di++) std::cout << ", {" << x(di)(0) << ", " << x(di)(1) << ", " << x(di)(2) << "}"; \
	std::cout << "}" << std::endl;}

#define DISPLAY_POOMA2_FLAT(x){             \
	if (x.length(0) > 0) std::cout << x(0)(0) << " " << x(0)(1) << " " << x(0)(2); \
	for (int di = 1; di < x.length(0); di++) std::cout << " " << x(di)(0) << " " << x(di)(1) << " " << x(di)(2); \
	std::cout << std::endl;}

#define DISPLAY_POOMA2Dv(x){             \
	std::cout << #x " = "; \
	std::cout << "{"; \
	if (x.length(0) > 0) std::cout << "{" << x(0)(0) << ", " << x(0)(1) << "}"; \
	for (int di = 1; di < x.length(0); di++) std::cout << ", {" << x(di)(0) << ", " << x(di)(1) << "}"; \
	std::cout << "};" << std::endl;}

#define DISPLAY_POOMA3Dv(x){             \
	std::cout << #x " = "; \
	std::cout << "{"; \
	if (x.length(0) > 0) std::cout << "{" << x(0)(0) << ", " << x(0)(1) <<  ", " << x(0)(2) << "}"; \
	for (int di = 1; di < x.length(0); di++) std::cout << ", {" << x(di)(0) << ", " << x(di)(1) << ", " << x(di)(2) << "}"; \
	std::cout << "};" << std::endl;}

#define DISPLAY_POOMA2(x){             \
	std::cout << #x " = "; \
	std::cout << "{"; \
	if (x.length(0) > 0) std::cout << "{" << x(0)(0) << ", " << x(0)(1) << ", " << x(0)(2) << "}"; \
	for (int di = 1; di < x.length(0); di++) std::cout << ", {" << x(di)(0) << ", " << x(di)(1) << ", " << x(di)(2) << "}"; \
	std::cout << "};" << std::endl;}

#define PRINT_POOMA1DARRAY(x){ \
	std::cout << "{"; \
	if (x.length(0) > 0) std::cout << x(0); \
	for (int di = 1; di < x.length(0); di++) std::cout << ", " << x(di); \
	std::cout << "}";\
}

#define DISPLAY_POOMA2DARRAY(x){ \
	std::cout << #x " = " << std::endl; \
	std::cout << "{"; \
	if (x.length(0) > 0) { \
		std::cout << "{"; \
		std::cout << x(0, 0); \
		for (int di = 1; di < x.length(1); di++) std::cout << ", " << x(0, di); \
		std::cout << "}";\
	}\
	for (int di1 = 1; di1 < x.length(0); di1++) { \
		std::cout << "{"; \
		std::cout << x(di1, 0); \
		for (int di = 1; di < x.length(1); di++) std::cout << ", " << x(di1, di); \
		std::cout << "}";\
	} \
	std::cout << "};" << std::endl;}

#define DISPLAY_POOMA2DARRAY_WOH(x){ \
	std::cout << "{"; \
	if (x.length(0) > 0) { \
		std::cout << "{"; \
		std::cout << x(0, 0); \
		for (int di = 1; di < x.length(1); di++) std::cout << ", " << x(0, di); \
		std::cout << "}";\
	}\
	for (int di1 = 1; di1 < x.length(0); di1++) { \
		std::cout << "{"; \
		std::cout << x(di1, 0); \
		for (int di = 1; di < x.length(1); di++) std::cout << ", " << x(di1, di); \
		std::cout << "}";\
	} \
	std::cout << "};" << std::endl;}

#define DISPLAY_STL1(x) { \
	cout << #x " = {"; \
	if (x.size() == 0) std::cout << ""; \
	else {copy(x.begin(), --x.end(), ostream_iterator<int>(std::cout, ", ")); \
		std::cout << *(x.rbegin());} \
	std::cout << "};" << std::endl; }

#define DISPLAY_STL_D1(x) { \
	cout << #x " = {"; \
	if (x.size() == 0) std::cout << ""; \
	else {std::copy(x.begin(), --x.end(), ostream_iterator<double>(std::cout, ", ")); \
		std::cout << *(x.rbegin());} \
	std::cout << "};" << std::endl; }

#define DISPLAY_STL_D1_WOH(x) { \
	cout << "{"; \
	if (x.size() == 0) std::cout << ""; \
	else {std::copy(x.begin(), --x.end(), ostream_iterator<double>(std::cout, ", ")); \
		std::cout << *(x.rbegin());} \
	std::cout << "};" << std::endl; }

#define DISPLAY_STL_D1_WB(x) { \
	cout << "{"; \
	if (x.size() == 0) std::cout << ""; \
	else {std::copy(x.begin(), --x.end(), ostream_iterator<double>(std::cout, ", ")); \
		std::cout << *(x.rbegin());} \
	std::cout << "}"; }

#define DISPLAY_STL_D2(x) { \
	std::cout << #x " = {"; \
	if (x.size() > 0) DISPLAY_STL_D1_WB(x[0]); \
	for(unsigned int pk2 = 1; pk2 < x.size(); pk2++) \
	{ \
		std::cout << ", "; \
		DISPLAY_STL_D1_WB(x[pk2]); \
	} \
	std::cout << "};" << std::endl; }

#define DISPLAY_STL_D2_WB(x) { \
	std::cout << "{"; \
	if (x.size() > 0) DISPLAY_STL_D1_WB(x[0]); \
	for(unsigned int pk2 = 1; pk2 < x.size(); pk2++) \
	{ \
		std::cout << ", "; \
		DISPLAY_STL_D1_WB(x[pk2]); \
	} \
	std::cout << "}" << std::endl; }

#define DISPLAY_STL_D3(x) { \
	std::cout << #x " = {"; \
	if (x.size() > 0) DISPLAY_STL_D2_WB(x[0]); \
	for(unsigned int pk3 = 1; pk3 < x.size(); pk3++) \
	{ \
		std::cout << ", "; \
		DISPLAY_STL_D2_WB(x[pk3]); \
	} \
	std::cout << "};" << std::endl; }

#define DISPLAY_STL_I1_WB(x) { \
	cout << "{"; \
	if (x.size() == 0) std::cout << ""; \
	else {copy(x.begin(), --x.end(), ostream_iterator<int>(std::cout, ", ")); \
		std::cout << *(x.rbegin());} \
	std::cout << "}"; }

#define DISPLAY_STL_I2_WB(x) { \
	std::cout << "{"; \
	if (x.size() > 0) DISPLAY_STL_I1_WB(x[0]); \
	for(unsigned int pk = 1; pk < x.size(); pk++) \
	{ \
		std::cout << ", "; \
		DISPLAY_STL_I1_WB(x[pk]); \
	} \
	cout << "}"; \
}

#define DISPLAY_STL_I3_WB(x) { \
	std::cout << "{"; \
	if (x.size() > 0) DISPLAY_STL_I2_WB(x[0]); \
	for(unsigned int pk2 = 1; pk2 < x.size(); pk2++) \
	{ \
		std::cout << ", "; \
		DISPLAY_STL_I2_WB(x[pk2]); \
	} \
	cout << "}"; \
}

#define DISPLAY_STL_I4_WB(x) { \
	std::cout << "{"; \
	if (x.size() > 0) DISPLAY_STL_I3_WB(x[0]); \
	for(unsigned int pk3 = 1; pk3 < x.size(); pk3++) \
	{ \
		std::cout << ", "; \
		DISPLAY_STL_I3_WB(x[pk3]); \
	} \
	cout << "}"; \
}

#define DISPLAY_STL_I4(x) { \
	std::cout << #x " = "; \
	DISPLAY_STL_I4_WB(x); \
	std::cout << ";" << std::endl; \
}

#define DISPLAY_STL1_FLAT(x) { \
	if (x.size() == 0) std::cout << ""; \
	else {copy(x.begin(), --x.end(), ostream_iterator<int>(std::cout, " ")); \
		std::cout << *(x.rbegin());} \
	std::cout << std::endl; \
}

#define DISPLAY_STL2(x) { \
	std::cout << #x " = {"; \
	std::cout << "{"; \
	if (x[0].size() > 0) { \
	copy(x[0].begin(), --x[0].end(), ostream_iterator<int>(std::cout, ", ")); \
	std::cout << *(x[0].rbegin()); } \
	std::cout << "}"; \
	for(unsigned int pk = 1; pk < x.size(); pk++) \
	{ \
		std::cout << ", {"; \
		if (x[pk].size() > 0) { \
		copy(x[pk].begin(), --x[pk].end(), ostream_iterator<int>(cout, ", ")); \
		std::cout << *(x[pk].rbegin()); } \
		std::cout << "}"; \
	} \
	cout << "};" << std::endl; \
}

#define DISPLAY_ARRAY1(x) { \
	int SIZE = sizeof(x)/sizeof(x[0]); \
	cout << #x " = {"; \
	if (SIZE == 0) std::cout << ""; \
	else {copy(x, x + SIZE - 1, ostream_iterator<int>(std::cout, ", ")); \
		std::cout << x[SIZE - 1];} \
	std::cout << "};" << std::endl; \
}

#define DISPLAY_ARRAY1n(x, SIZE) { \
	std::cout << #x " = {"; \
	if (SIZE == 0) std::cout << ""; \
	else if (SIZE == 1) std::cout << x[0]; \
	else {copy(x, x + SIZE - 1, ostream_iterator<int>(std::cout, ", ")); \
		std::cout << x[SIZE - 1];} \
	std::cout << "};" << std::endl; \
}

#define DISPLAY_DARRAY1(x) { \
	int SIZE = sizeof(x)/sizeof(x[0]); \
	cout << #x " = {"; \
	if (SIZE == 0) std::cout << ""; \
	else {copy(x, x + SIZE - 1, ostream_iterator<double>(std::cout, ", ")); \
		std::cout << x[SIZE - 1];} \
	std::cout << "};" << std::endl; \
}

#define DISPLAY_DARRAY1n(x, SIZE) { \
	std::cout << #x " = {"; \
	if (SIZE == 0) std::cout << ""; \
	else if (SIZE == 1) std::cout << x[0]; \
	else {copy(x, x + SIZE - 1, ostream_iterator<double>(std::cout, ", ")); \
		std::cout << x[SIZE - 1];} \
	std::cout << "};" << std::endl; \
}

#define DISPLAY_DARRAY2n(xx, nROWS, nCOLS) { \
	std::cout << #xx " = {"; \
	int qq; \
	for (qq = 0; qq < nROWS - 1; qq++){ \
		cout << "{"; \
		if (nCOLS == 0) std::cout << ""; \
		else {copy(xx + qq*nCOLS, xx + qq*nCOLS + nCOLS - 1, ostream_iterator<double>(std::cout, ", ")); \
		std::cout << (xx + qq*nCOLS)[nCOLS - 1];} \
		std::cout << "}, "; \
	}\
	qq = nROWS - 1; \
	cout << "{"; \
	if (nCOLS == 0) std::cout << ""; \
	else {copy(xx + qq*nCOLS, xx + qq*nCOLS + nCOLS - 1, ostream_iterator<double>(std::cout, ", ")); \
	std::cout << (xx + qq*nCOLS)[nCOLS - 1];} \
	std::cout << "}"; \
	std::cout << "};" << std::endl; \
}

#define DISPLAY_DARRAY2n0(xx, nROWS, nCOLS) { \
	std::cout << #xx " = {"; \
	int qq; \
	for (qq = 0; qq < nROWS - 1; qq++){ \
		cout << "{"; \
		if (nCOLS == 0) std::cout << ""; \
		else {copy(&(xx[qq][0]), &(xx[qq][0]) + nCOLS - 1, ostream_iterator<double>(std::cout, ", ")); \
		std::cout << xx[qq][nCOLS - 1];} \
		std::cout << "}, "; \
	}\
	qq = nROWS - 1; \
	cout << "{"; \
	if (nCOLS == 0) std::cout << ""; \
	else {copy(&(xx[qq][0]), &(xx[qq][0]) + nCOLS - 1, ostream_iterator<double>(std::cout, ", ")); \
	std::cout << xx[qq][nCOLS - 1];} \
	std::cout << "}"; \
	std::cout << "};" << std::endl; \
}

#define NORMALIZE(x) (x)/norm(x)

#define pooma_CROSS(dest, v1, v2){                 \
          dest(0) = v1(1) * v2(2) - v1(2) * v2(1); \
          dest(1) = v1(2) * v2(0) - v1(0) * v2(2); \
          dest(2) = v1(0) * v2(1) - v1(1) * v2(0);}

#define pooma_DOT(v1, v2) (v1(0) * v2(0) + v1(1) * v2(1) + v1(2) * v2(2))

#define pooma_SUB(dest, v1, v2){       \
          dest(0) = v1(0) - v2(0); \
          dest(1) = v1(1) - v2(1); \
          dest(2) = v1(2) - v2(2);}

inline double pow_di(double &ap, int bp)
{
	double pow, x;
	int n;
	unsigned long u;

	pow = 1;
	x = ap;
	n = bp;

	if(n != 0)
	{
		if(n < 0)
		{
			n = -n;
			x = 1/x;
		}
		for(u = n; ; )
		{
			if(u & 01)
				pow *= x;
			if(u >>= 1)
				x *= x;
			else
				break;
		}
	}
	return(pow);
}




#ifdef dummy
struct overlapTest {
	bool result;
	double d;
	double theta;
	double det;
	Vector<3> tpivot1;
	Vector<3> tpivot2;
	Vector<3> tnormal11;
	Vector<3> tnormal12;
	Vector<3> tnormal21;
	Vector<3> tnormal22;
	Vector<3> tP11;
	Vector<3> tP12;
	Vector<3> tP21;
	Vector<3> tP22;
	//Vector<3> otP11;
	//Vector<3> otP12;
	//Vector<3> otP21;
	//Vector<3> otP22;
};

struct bond{
	int n1, n2;// ensure n1 < n2 yourself
	bool operator==(const bond & p) const;
	bool operator<(const bond &p) const ;
	void setValue(const int &, const int &);
};

struct bondPair{
	bond n1, n2;// ensure n1 < n1 yourself
	double overlap;
	overlapTest Q;
	Vector<3> P11;
	Vector<3> P12;
	Vector<3> P21;
	Vector<3> P22;
	bool operator==(const bondPair & p) const;
	bool operator<(const bondPair &p) const ;
	void setValue(const bond &, const bond &);
};

inline bool bond::operator<(const bond & p) const {
		return this->n1 < p.n1 || (this->n1 == p.n1 && this->n2 < p.n2);
}

inline void bond::setValue(const int & m1, const int & m2){
	this->n1 = m1;
	this->n2 = m2;
	if (m2 < m1) std::swap(this->n1, this->n2);
}

inline bool bond::operator==(const bond & p) const {
		return (this->n1 == p.n1) && (p.n2 == this->n2);
}

inline bool bondPair::operator<(const bondPair & p) const {
		return this->n1 < p.n1 || (this->n1 == p.n1 && this->n2 < p.n2);
}

inline bool bondPair::operator==(const bondPair & p) const {
		return (this->n1 == p.n1) && (p.n2 == this->n2);
}

inline void bondPair::setValue(const bond & m1, const bond & m2){
	this->n1 = m1;
	this->n2 = m2;
	if (m2 < m1) std::swap(this->n1, this->n2);
}

inline void printSTLvectorOfPOOMAVectors(const std::vector<Vector<3, double, Full> > & Coords, const std::string & Name)
{
	std::cout << Name << " = {";
	for (unsigned int jj = 0; jj < Coords.size() - 1; jj++){
		DISPLAY_POOMAV_WOH(Coords[jj]);
		std::cout << ", ";
	}
	DISPLAY_POOMAV_WOH(Coords[Coords.size() - 1]);
	std::cout << "};\n";
}

extern Vector<3> Ox, Oy, Oz, ZeroVector;
extern Tensor<3, double, Full> zeroTensor3;
extern Tensor<3, double, Diagonal> identity3;
extern Tensor<2, double, Diagonal> identity2;
extern std::vector<int> BINCmat;


inline double Chop(const double & x)
{
	return (fabs(x) < MY_EPSILON4CHOP) ? 0 : x;
}

inline double Chop(const double & x, const double & tol)
{
	return (fabs(x) < tol) ? 0 : x;
}

int prod (int x, int y);

template<class T>
struct index_cmp
{
	index_cmp(const T arr) : arr(arr) {}
	bool operator()(const size_t a, const size_t b) const
	{
		return arr[a] > arr[b];
	}
	const T arr;
};

template<class E1, class E2>
inline Vector<3> Cross(const Vector<3, double, E1> & v1, const Vector<3, double, E2> & v2)
{
	Vector<3> temp;
	temp(0) = v1(1) * v2(2) - v1(2) * v2(1);
	temp(1) = v1(2) * v2(0) - v1(0) * v2(2);
	temp(2) = v1(0) * v2(1) - v1(1) * v2(0);
	return temp;
}

template<class E1, class E2>
inline void axisRotation(const double & theta, const Vector<3, double, E1> & axis, Tensor<3, double, E2> & refl)
{
	double c = cos(theta);
	double s = sin(theta);
	double t = 1 - c;
	double x = axis(0);
	double y = axis(1);
	double z = axis(2);

	// build the rotation matrix
	refl(0, 0) = t*x*x + c;
	refl(0, 1) = t*x*y - s*z;
	refl(0, 2) = t*x*z + s*y;

	refl(1, 0) = t*x*y + s*z;
	refl(1, 1) = t*y*y + c;
	refl(1, 2) = t*y*z - s*x;

	refl(2, 0) = t*x*z - s*y;
	refl(2, 1) = t*y*z + s*x;
	refl(2, 2) = t*z*z + c;

}

template<class E>
inline Vector<3> Normalize(const Vector<3, double, E> & V)
{
	double d = norm(V);
	Vector<3> res = (d < ROT_EPSILON) ? V : V/d;
	return res;
}


template<class E1, class E2, class E3>
inline void fromToRotation2(const Vector<3, double, E1> & from1, const Vector<3, double, E2> & to1, Tensor<3, double, E3> & mtx)
{
	Vector<3> v;
	double e, h, f;
	Vector<3> from(Normalize(from1)), to(Normalize(to1));
	pooma_CROSS(v, from, to);
	e = pooma_DOT(from, to);
	h = 1.0/(1.0 + e);  	/* My addition */
	f = (e < 0)? -e:e;
	if (f > 1.0 - ROT_EPSILON || std::isinf(h) || std::isnan(h))     /* "from" and "to"-vector almost parallel */
	{
		Vector<3> u, v; /* temporary storage vectors */
		Vector<3> x;       /* vector most nearly orthogonal to "from" */
		double c1, c2, c3; /* coefficients for later use */
		int i, j;

		x(0) = (from(0) > 0.0)? from(0) : -from(0);
		x(1) = (from(1) > 0.0)? from(1) : -from(1);
		x(2) = (from(2) > 0.0)? from(2) : -from(2);

		if (x(0) < x(1))
		{
			if (x(0) < x(2))
			{
				x(0) = 1.0; x(1) = x(2) = 0.0;
			}
			else
			{
				x(2) = 1.0; x(0) = x(1) = 0.0;
			}
		}
		else
		{
			if (x(1) < x(2))
			{
				x(1) = 1.0; x(0) = x(2) = 0.0;
			}
			else
			{
				x(2) = 1.0; x(0) = x(1) = 0.0;
			}
		}

		u(0) = x(0) - from(0); u(1) = x(1) - from(1); u(2) = x(2) - from(2);
		v(0) = x(0) - to(0);   v(1) = x(1) - to(1);   v(2) = x(2) - to(2);

		c1 = 2.0 / pooma_DOT(u, u);
		c2 = 2.0 / pooma_DOT(v, v);
		c3 = c1 * c2  * pooma_DOT(u, v);

		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				mtx(i, j) =  - c1 * u(i) * u(j)
                    		 - c2 * v(i) * v(j)
                    		 + c3 * v(i) * u(j);
			}
			mtx(i, i) += 1.0;
		}
	}
	else  /* the most common case, unless "from"="to", or "from"=-"to" */
	{
#if 0
		/* unoptimized version - a good compiler will optimize this. */
		/* h = (1.0 - e)/pooma_DOT(v, v); old code */
		h = 1.0/(1.0 + e);      /* optimization by Gottfried Chen */
		mtx(0, 0) = e + h * v(0) * v(0);
		mtx(0, 1) = h * v(0) * v(1) - v(2);
		mtx(0, 2) = h * v(0) * v(2) + v(1);

		mtx(1, 0) = h * v(0) * v(1) + v(2);
		mtx(1, 1) = e + h * v(1) * v(1);
		mtx(1, 2) = h * v(1) * v(2) - v(0);

		mtx(2, 0) = h * v(0) * v(2) - v(1);
		mtx(2, 1) = h * v(1) * v(2) + v(0);
		mtx(2, 2) = e + h * v(2) * v(2);
#else
	/* ...otherwise use this hand optimized version (9 mults less) */
		double hvx, hvz, hvxy, hvxz, hvyz;
		/* h = (1.0 - e)/pooma_DOT(v, v); old code */
		//h = 1.0/(1.0 + e);      /* optimization by Gottfried Chen */
		hvx = h * v(0);
		hvz = h * v(2);
		hvxy = hvx * v(1);
		hvxz = hvx * v(2);
		hvyz = hvz * v(1);
		mtx(0, 0) = e + hvx * v(0);
		mtx(0, 1) = hvxy - v(2);
		mtx(0, 2) = hvxz + v(1);

		mtx(1, 0) = hvxy + v(2);
		mtx(1, 1) = e + h * v(1) * v(1);
		mtx(1, 2) = hvyz - v(0);

		mtx(2, 0) = hvxz - v(1);
		mtx(2, 1) = hvyz + v(0);
		mtx(2, 2) = e + hvz * v(2);
#endif
	}
}

inline int PQselector(const int & S, const int & P, const int & O)
{
	return (S == P) - (S == O);
}

template<class E1, class E2>
inline Tensor<3> cross(const Tensor<3, double, E1> & M, const Vector<3, double, E2> & v)
{
	Vector<3> m;
	Tensor<3> Mv;
	m(0) = M(0, 0);
	m(1) = M(1, 0);
	m(2) = M(2, 0);
	m = Cross(m, v);
	Mv(0, 0) = m(0);
	Mv(1, 0) = m(1);
	Mv(2, 0) = m(2);
	m(0) = M(0, 1);
	m(1) = M(1, 1);
	m(2) = M(2, 1);
	m = Cross(m, v);
	Mv(0, 1) = m(0);
	Mv(1, 1) = m(1);
	Mv(2, 1) = m(2);
	m(0) = M(0, 2);
	m(1) = M(1, 2);
	m(2) = M(2, 2);
	m = Cross(m, v);
	Mv(0, 2) = m(0);
	Mv(1, 2) = m(1);
	Mv(2, 2) = m(2);
	return Mv;
}

template<class E>
inline Tensor<3> cross(const Vector<3, double, E> & u)
{
	Tensor<3> Mx;
	Mx(0 , 0) = 0;
	Mx(0 , 1) = -u(2);
	Mx(0 , 2) = u(1);
	Mx(1 , 0) = u(2);
	Mx(1 , 1) = 0;
	Mx(1 , 2) = -u(0);
	Mx(2 , 0) = -u(1);
	Mx(2 , 1) = u(0);
	Mx(2 , 2) = 0;
	return Mx;
}

template<class E1, class E2, class E3, class E4>
inline Tensor<3> d_SCR_axisUnitVector__d_rS(const int & S, const int & P, const int & O, const Vector<3, double, E1> & cv, const Vector<3, double, E2> & bv, const Vector<3, double, E3> & SCR_X, const Vector<3, double, E4> & SCR_Y, const double & cphi, const double & sphi, const double & bc, const double & bb, const double & bxc, const double & bxbxc)
{
	int pqs = PQselector(S, P, O);
#ifdef CHECK_RG
	DISPLAY_VAR(pqs);
	if (pqs == 0) abort();
#endif
	Tensor<3> bvXcv = outerProduct(bv, cv);
	Tensor<3> cvXbv = transpose(bvXcv);
	if (pqs == 0) return zeroTensor3;
	else return pqs*(cphi/bxc*(cross(-cv) + outerProduct(SCR_X/bxc, bc*cv - bv)) + sphi/bxbxc*(bc*identity3 + bvXcv - 2*cvXbv - outerProduct(SCR_Y/bxbxc, (2*bb - bc*bc)*bv - bb*bc*cv)));
	return zeroTensor3;
}

template<class E1, class E2, class E3, class E4, class E5, class E6>
inline Tensor<3> forgottenSCRJacobian(const Vector<3, double, E1> & axisUnitVector, const double & c, const double & s, const double & t, Vector<3, double, E2> & victim, const int & S, const int & P, const int & O, const Vector<3, double, E3> & cv, const Vector<3, double, E4> & bv, const Vector<3, double, E5> & SCR_X, const Vector<3, double, E6> & SCR_Y, const double & cphi, const double & sphi, const double & bc, const double & bb, const double & bxc, const double & bxbxc)
{
	Tensor<3> M, Mt;
	M = d_SCR_axisUnitVector__d_rS(S, P, O, cv, bv, SCR_X, SCR_Y, cphi, sphi, bc, bb, bxc, bxbxc);
	Mt = transpose(M);
	return s*cross(M, victim) + t*(dot(axisUnitVector, victim)*M + outerProduct(axisUnitVector, dot(Mt, victim)));
}

template<class E1>
inline Tensor<3> d_CSR_axisUnitVector__d_rS(const int & S, const int & P, const int & O, const Vector<3, double, E1> & uv, const double & el)
{
	int pqs = PQselector(S, P, O);
	if (pqs == 0) return zeroTensor3;
	else {
		return pqs/el*(identity3 - outerProduct(uv, uv));
	}
	return zeroTensor3;
}

template<class E1, class E2>
inline Tensor<3> forgottenCSRJacobian(const Vector<3, double, E1> & axisUnitVector, const double & c, const double & s, const double & t, const Vector<3, double, E2> & victim, const int & S, const int & P, const int & O, const double & el)
{
	Tensor<3> M, Mt;
	Mt = transpose(M);
	M = d_CSR_axisUnitVector__d_rS(S, P, O, axisUnitVector, el);
	return s*cross(M, victim) + t*(dot(axisUnitVector, victim)*M + outerProduct(axisUnitVector, dot(Mt, victim)));
}

template<class E1, class E2, class E3, class E4, class E5>
int getr1(Vector<3, double, E1> & r1p, Tensor<3, double, E2> & jac, Array<1, Vector<3, double, E3>, E4> trimer, const Vector<3, double, E5> & r2p, const double & reflen)
{
	Tensor<3> rot;
	Vector<3> r02 = trimer(2) - trimer(0);
	Vector<3> r02p = r2p - trimer(0);
	Vector<3> r01 = trimer(1) - trimer(0);
	double d2 = norm2(r02p);
	if (Chop(d2) == 0 || Chop(norm2(r02)) == 0) {
		return 0;
	}
	else {
		fromToRotation2(r02, r02p, rot);
	}
	Vector<3> r0 = trimer(0);
	Vector<3> r1 = r0 + dot(rot, r01);
	Vector<3> r2 = r0 + dot(rot, r02);
	double d = norm(r02);
	double dp = norm(r2p - r0);
	r01 = r1 - r0;
	Vector<3> ur02 = Normalize(r02p);
	//double b1 = norm(r01);
	//double b2 = norm(r2 - r1);
	double b1 = reflen;
	double b2 = reflen;
	double s = (b1 + b2 + dp)/2.;
	double h = 2.*sqrt(s*(s - b1)*(s - b2)*(s - dp))/dp;
	if (Chop(h) == 0) {
		r1p = r0 + b1*ur02;
#ifdef CHECK_NEW_JAC
		if (std::isnan(r1p(0)) || std::isnan(r1p(1)) || std::isnan(r1p(2)) \
				|| std::isinf(r1p(0)) || std::isinf(r1p(1)) || std::isinf(r1p(2))) {
			DISPLAY_VAR(h);
			DISPLAY_VAR(s);
			DISPLAY_VAR(d);
			DISPLAY_VAR(dp);
			DISPLAY_VAR(norm(r01));
			DISPLAY_VAR(norm(r2 - r1));
			DISPLAY_POOMA2(trimer);
			DISPLAY_POOMAV(r1p);
			abort();
		}
#endif
		jac = cbrt(d/dp)*identity3;  // Not really correct, but determinant is.
		return 1;
	}
	Vector<3> rperp = r01 - dot(r01, ur02)*ur02;
	double h0 = norm(rperp);
	if (h0 == 0) return 0;
	Vector<3> urperp = Normalize(rperp);
	Tensor<3> ur02Xur02 = outerProduct(ur02, ur02);
	Tensor<3> urperpXurperp = outerProduct(urperp, urperp);
	if (Chop(b1 - b2) == 0) {
#ifdef CHECK_RG
		printf("Using NEW transformation for trimers. \n");
#endif
		r1p = r0 + ur02*dp/2. + h*urperp;
#ifdef CHECK_NEW_JAC
		if (std::isnan(r1p(0)) || std::isnan(r1p(1)) || std::isnan(r1p(2)) \
				|| std::isinf(r1p(0)) || std::isinf(r1p(1)) || std::isinf(r1p(2))) {
			DISPLAY_VAR(h);
			DISPLAY_VAR(s);
			DISPLAY_VAR(d);
			DISPLAY_VAR(dp);
			DISPLAY_VAR(norm(r01));
			DISPLAY_VAR(norm(r2 - r1));
			DISPLAY_POOMA2(trimer);
			DISPLAY_POOMAV(r1p);
			abort();
		}
#endif
		jac = (h/h0)*identity3 + ((d/dp) - (h/h0))*ur02Xur02 + ((h0/h) - (h/h0))*urperpXurperp;
#ifdef CHECK_RG
		printf("Finished NEW transformation for trimers. \n");
		DISPLAY_POOMAV(r1p);
#endif
		return 1;
	}
	return 0;
}


inline void SLA(double & area, const double &c, const double & r1, const double & r2)
{
	// DISJOINT CIRCLES
	if (Chop(c - r1 - r2) >= 0){
#ifdef CHECK_RG
		printf("Circles are disjoint.\n");
#endif
	area = 0;
		return;
	}
	// CIRCLE 1 IS WITHIN CIRCLE 2
	else if (Chop(c + r1 - r2) <= 0){
#ifdef CHECK_RG
		printf("Circle 1 is within Circle 2.\n");
#endif
		area = 1.0 - cos(r1);
		area *= 2.0*PI;
		return;
	}
	// CIRCLE 2 IS WITHIN CIRCLE 1
	else if (Chop(c + r2 - r1) <= 0){
#ifdef CHECK_RG
		printf("Circle 2 is within Circle 1.\n");
#endif
		area = 1.0 - cos(r2);
		area *= 2.0*PI;
		return;
	}
	// DISJOINT CIRCLE-COMPLEMENTS
	else if (c > PI - r2 + PI - r1){
#ifdef CHECK_RG
		printf("Complements of both circles are disjoint.\n");
#endif
		area = 4.0*PI;
		area -= 2.0*PI*(1.0 - cos(PI - r1));
		area -= 2.0*PI*(1.0 - cos(PI - r2));
		return;
	}
/*	// COMPLEMENT OF CIRCLE 2 IS WITHIN CIRCLE 1
	else if (Chop(PI - c + PI - r2 - r1) <= 0){
		area = 1.0 - cos(r1);
		area *= 2.0*PI;
		area -= 2.0*PI*(1.0 - cos(PI - r2));
		return;
	}
	// COMPLEMENT OF CIRCLE 1 IS WITHIN CIRCLE 2
	else if (Chop(PI - c + PI - r1 - r2) <= 0){
		area = 1.0 - cos(r2);
		area *= 2.0*PI;
		area -= 2.0*PI*(1.0 - cos(PI - r1));
		return;
	}
*/
	// semi-perimeter of spherical triangle:
	double s = 0.5*(c + r1 + r2);
	// spherical excess/solid angle/area of spherical triangle:
	double tr = 4.0*atan(sqrt(tan(0.5*s)*tan(0.5*(s - c))*tan(0.5*(s - r1))*tan(0.5*(s - r2))));
	tr += ((tr < 0)? 4.0*PI : 0);
	// area of sector of small circle 1
	double area1 = 1.0 - cos(r1);
	double theta1 = 2.0*acos(sqrt(sin(s)*sin(s - r2)/sin(r1)/sin(c)));
	area1 *= theta1;
	// area of sector of small circle 2
	double area2 = 1.0 - cos(r2);
	double theta2 = 2.0*acos(sqrt(sin(s)*sin(s - r1)/sin(r2)/sin(c)));
	area2 *= theta2;
	// area of intersection of small circles 1 and 2
	area = 2.0*(area1 + area2 - tr);
#ifdef CHECK_RG
	printf("Circle borders intersect.\n");
	DISPLAY_VAR(c);
	DISPLAY_VAR(r1);
	DISPLAY_VAR(r2);
	DISPLAY_VAR(s);
	DISPLAY_VAR(tr);
	DISPLAY_VAR(area1);
	DISPLAY_VAR(area2);
#endif
	return;
}


inline int lti(const int & i, const int & j)
{
	/*	Index adaptor for a lower triangular
	 *  matrix.
	 */
	return i*(i + 1)/2 + j;
}


inline int binC(const int & n, const int & r)
{
	return BINCmat[lti(n, r)];
}


inline double pgen(const int & pos, const double & area0, const double & area1, const int & nSamples, const int & nTrials)
{
	double pg;
	int maxFailures = nTrials - nSamples;
	int n = nTrials, k = nSamples, p = pos;
	double s = area1/area0;	// probability of success in one trial
	double f = 1. - s;	// probability of failure in one trial
#ifdef CHECK_RG
	DISPLAY_VAR(s);
	DISPLAY_VAR(f);
	DISPLAY_VAR(n);
	DISPLAY_VAR(k);
	DISPLAY_VAR(p);
	DISPLAY_VAR(area0);
	DISPLAY_VAR(area1);
#endif
	if (s - 1. == 0){
		pg = 1./area0;
	}
	else if (nSamples == 1){
		pg = (1. - pow(f, (double) n))/area1;
	}
	else if (pos == 0 || pos == (nSamples - 2)){// first and last but one
		double sum = 0.;
		for(int j = 0; j <= maxFailures; j++){
			sum += (pow(f, (double)j)*binC(k - 2 + j, k - 2));
#ifdef CHECK_RG
			DISPLAY_VAR(j);
			DISPLAY_VAR(sum);
			DISPLAY_VAR(k - 2 + j);
			DISPLAY_VAR(k - 2);
			DISPLAY_VAR(binC(k - 2 + j, k - 2));
#endif
		}
#ifdef CHECK_RG
		DISPLAY_VAR(n - 1);
		DISPLAY_VAR(k - 1);
		DISPLAY_VAR(binC(n - 1, k - 1));
		DISPLAY_VAR(sum);
#endif
		pg = pow(s, (double)(k - 2))*(sum - pow(f, (double)(n - k + 1))*binC(n - 1, k - 1))/area0;
	}
	else if (pos == (nSamples - 1)){// last sample
		double sum = 0.;
		for(int j = 0; j <= maxFailures; j++){
			sum += pow(f, (double)j)*binC(k - 1 + j, k - 1);
#ifdef CHECK_RG
			DISPLAY_VAR(j);
			DISPLAY_VAR(k - 1 + j);
			DISPLAY_VAR(k - 1);
			DISPLAY_VAR(binC(k - 1 + j, k - 1));
#endif
		}
		pg = pow(s, (double)(k - 1))*sum/area0;
	}
	else {
		double sum = 0.;
		for(int j = 0; j <= maxFailures; j++){
			for(int m = 0; m <= maxFailures - j; m++){
				sum += pow(f, (double)(j + m))*binC(p + j, p)*binC(m + k - p - 2, k - p - 2);
#ifdef CHECK_RG
				DISPLAY_VAR(j);
				DISPLAY_VAR(m);
				DISPLAY_VAR(m + k - p - 2);
				DISPLAY_VAR(k - p - 2);
				DISPLAY_VAR(binC(m + k - p - 2, k - p - 2));
#endif
			}
		}
		pg = pow(s, (double)(k - 1))*sum/area0;
	}
#ifdef CHECK_RG
	DISPLAY_VAR(pg);
#endif
	return pg;
}

template<class T, class E>
struct ltVector
{
  bool operator()(const Vector<3, T, E> & s1, const Vector<3, T, E> & s2) const
  {
    return (s1(0) <  s2(0)) \
    		|| (s1(0) == s2(0) && s1(1) < s2(1)) \
    		|| (s1(0) ==  s2(0) && s1(1) == s2(1) && s1(2) < s2(2));
  }
};


#endif /*0*/
#endif /* MYDEFS_HPP_ */
