#include "crosslinkTriple2.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

CrosslinkTriple2::CrosslinkTriple2(Crosslink &cl1, Crosslink &cl2, Crosslink &cl3, double &RR12a, double &RR23a, double &RR12b, double &RR23b, int &WBX12, int &WBY12, int &WBZ12,
					 int &WBX23, int &WBY23, int &WBZ23, int iDummy, int filId1,int filId2, int filId3):
angle_(0),epotTriple_(0), fPull_(0)
{
	cl1_ = &cl1;
	cl2_ = &cl2;
	cl3_ = &cl3;
	contourLength12a_ = &RR12a;
	contourLength23a_ = &RR23a;
	contourLength12b_ = &RR12b;
	contourLength23b_ = &RR23b;	
	WBX12_ = &WBX12;
	WBY12_ = &WBY12;
	WBZ12_ = &WBZ12;
	WBX23_ = &WBX23;
	WBY23_ = &WBY23;
	WBZ23_ = &WBZ23;

	iDummy_=iDummy;
	filIds_.push_back(filId1);
	filIds_.push_back(filId2);
	filIds_.push_back(filId3);
	//filId_ = filId;
}


void CrosslinkTriple2::updateEpotTriple(double epotTriple){
	epotTriple_=epotTriple;
}

void CrosslinkTriple2::updateAngle(double angle){
	angle_=angle;
}

double CrosslinkTriple2::getAngle(){
	return angle_;
}

Crosslink* CrosslinkTriple2::getCrosslink(int i){
 	if(i==1){
 		return cl1_;
 	}
 	else if (i==2){
 		return cl2_;
	}
	else{
		return cl3_;
	}
}

/*void CrosslinkTriple2::updateContourLength(double contourLength, int i){
	if(i==1){
		contourLength1_=contourLength;
	}
	else{
		contourLength2_=contourLength;
	}
}*/

double CrosslinkTriple2::getContourLength(int i){
	if(i==12){
		return *contourLength12a_;
	}
	else if(i==23){
		return *contourLength23a_;
	}
}

void CrosslinkTriple2::incrContourLength(double ds,int i){
	if(i==12){
		*contourLength12a_+=ds;
		if(iDummy_!=1)
		{*contourLength12b_+=ds;}
	}
	else if(i==23){
		*contourLength23a_+=ds;
		if(iDummy_!=3)
		{*contourLength23b_+=ds;}
	}
}

/*void CrosslinkTriple2::updateWB(ublas::vector<int> WB, int i){
	if(i==1){
		WB1_=WB;
	}
	else{
		WB2_=WB;
	}
}*/

int CrosslinkTriple2::getWBX(int i){
	if(i==12){
		return *WBX12_;
	}
	else if(i==23){
		return *WBX23_;
	}
}

int CrosslinkTriple2::getWBY(int i){
	if(i==12){
		return *WBY12_;
	}
	else if(i==23){
		return *WBY23_;
	}
}

int CrosslinkTriple2::getWBZ(int i){
	if(i==12){
		return *WBZ12_;
	}
	else if(i==23){
		return *WBZ23_;
	}
}

double CrosslinkTriple2::getfPull(){
	return fPull_;
}
void CrosslinkTriple2::updatefPull(double fPull){
	fPull_=fPull;
}

int CrosslinkTriple2::getDummy(){
	return iDummy_;
}

int CrosslinkTriple2::getFilId(int i){
	return filIds_[i-1];
}