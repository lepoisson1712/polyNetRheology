#include <typeinfo>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "crosslink.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;


Crosslink::Crosslink(int id, ublas::vector<double> r, ublas::matrix<int> conTabSelf, ublas::matrix<int> conTabOther):
	r_(3,0), f_(3,0),f1_(3,0),f2_(3,0)
	{
		id_ = id;
		r_ = r;

		conTabSelf_ = conTabSelf;
		conTabOther_ = conTabOther;
	}



void Crosslink::updatePosition(ublas::vector<double> r) {
	r_ = r;
}

ublas::vector<double> Crosslink::incrPosition(ublas::vector<double> r) {
	r_ += r;
	return r_;
}

void Crosslink::updateForce(ublas::vector<double> f) {
	f_ = f;
}

void Crosslink::updateFilForce(ublas::vector<double> f, int filId){
	if(filId==1){ f1_=f;}
	if(filId==2){ f2_=f;}
}

void Crosslink::incrForce(ublas::vector<double> f, int filId) {
	f_ += f;
	if(filId==1){f1_+=f;}
	else if(filId==2){f2_+=f;}
	else {cout << "FilId error"<< "\n";}

}

ublas::vector<double> Crosslink::getForce() {
	return f_;
}
	
ublas::vector<double> Crosslink::getFilForce(int filId){
	if(filId==1){ return f1_; }
	else if (filId==2){ return f2_; }
	else {cout << "FilId error2"<< "\n";}
}



ublas::vector<double> Crosslink::getPosition() {
	return r_;
}

int Crosslink::getId() {
	return id_;
}

ublas::matrix<int> Crosslink::getConTabSelf(){
	return conTabSelf_;
}

ublas::matrix<int> Crosslink::getConTabOther(){
	return conTabOther_;
}
