#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#ifndef SPHEROCYLINDER_H
#define SPHEROCYLINDER_H

using namespace std;
// using namespace boost::numeric::ublas;
namespace ublas = boost::numeric::ublas;

class Crosslink {
public:
	//Crosslink(int id, ublas::vector<double> r);
	Crosslink(int id, ublas::vector<double> r, ublas::matrix<int> conTabSelf, ublas::matrix<int> conTabOther);

    void updatePosition(ublas::vector<double> r);

	ublas::vector<double> incrPosition(ublas::vector<double> r);

	void updateForce(ublas::vector<double> f);

	void incrForce(ublas::vector<double> f, int filId);

	ublas::vector<double> getPosition();


	ublas::vector<double> getForce();
	ublas::vector<double> getFilForce(int filId);
	int getId();

	ublas::matrix<int> getConTabSelf();
	ublas::matrix<int> getConTabOther();

	void updateFilForce(ublas::vector<double> f, int filId);

private:

	ublas::vector<double> r_;
	
	ublas::vector<double> f_;

	ublas::matrix<int> conTabSelf_; // gives indices for WBX,WBY,WBZ tables that stand for the connections of the crosslink and its neighbors
	ublas::matrix<int> conTabOther_;


	int id_;

	ublas::vector<double> f1_; //force from filament 1
	ublas::vector<double> f2_; //force from filament 2

};


#endif