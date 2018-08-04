#include "crosslinkPair2.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

CrosslinkPair2::CrosslinkPair2(Crosslink &cl1, Crosslink &cl2, double &RRij, int &WBXij, int &WBYij, int &WBZij, int filId1, int filId2):
epotPair_(0), k2(0), Leq(0), r12Abs(0), f1(0)
{	
	cl1_=&cl1;
	cl2_=&cl2;
	contourLength_ = &RRij;
	WBX_ = &WBXij;
	WBY_ = &WBYij;
	WBZ_ = &WBZij;
	filIds_.push_back(filId1);
	filIds_.push_back(filId2);
	//eqLength_=ublas::norm_2(cl1_->getPosition()-cl2_->getPosition());
}

void CrosslinkPair2::updateEpotPair(double epotPair){
	epotPair_=epotPair;
}

/*void CrosslinkPair2::updateContourLength(double contourLength){
	contourLength_=contourLength;
}*/

double CrosslinkPair2::getContourLength(){
	return *contourLength_;
}

Crosslink* CrosslinkPair2::getCrosslink(int i){
 	if(i==1){
 		return cl1_;
 	}
 	else{
 		return cl2_;
 	}
}

/*void CrosslinkPair2::updateWB(int dWBX, int dWBY, int dWBZ){
	WBX_+=dWBX;
	WBY_+=dWBY;
	WBZ_+=dWBZ;
}*/

int CrosslinkPair2::getWBX(){
	return *WBX_;
}
int CrosslinkPair2::getWBY(){
	return *WBY_;
}
int CrosslinkPair2::getWBZ(){
	return *WBZ_;
}

int CrosslinkPair2::getFilId(int i){
	return filIds_[i-1];
}