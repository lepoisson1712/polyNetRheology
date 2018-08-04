#include <boost/numeric/ublas/vector.hpp>
#include "crosslink.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

class CrosslinkTriple2 {
public:
	CrosslinkTriple2(Crosslink &cl1, Crosslink &cl2, Crosslink &cl3, double &RR12a, double &RR23a, double &RR12b, double &RR23b, int &WBX12, int &WBY12, int &WBZ12,
					 int &WBX23, int &WBY23, int &WBZ23, int iDummy, int filId1, int filId2, int filId3);

	void updateEpotTriple(double epotTriple);
	void updateAngle(double angle);
	double getAngle();
	Crosslink* getCrosslink(int i);
	//void updateContourLength(double contourLength, int i);
	double getContourLength(int i);
	//void updateWB(ublas::vector<int> WB, int i);
	int getWBX(int i);
	int getWBY(int i);
	int getWBZ(int i);

	double getfPull();
	void updatefPull(double fPull);

	void incrContourLength(double ds,int i);
	int getFilId(int i);

	int getDummy();

private:
	Crosslink *cl1_;
	Crosslink *cl2_;
	Crosslink *cl3_;

	
	double angle_;
	double epotTriple_;
	double *contourLength12a_;
	double *contourLength23a_;
	double *contourLength12b_;
	double *contourLength23b_;	
	int *WBX12_, *WBY12_, *WBZ12_;
	int *WBX23_, *WBY23_, *WBZ23_;

	double fPull_;
	int iDummy_;

	std::vector<int> filIds_; 
};
