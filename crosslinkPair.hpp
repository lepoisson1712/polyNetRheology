#include <boost/numeric/ublas/vector.hpp>
#include "crosslink.hpp"
//#include "simulation.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

class CrosslinkPair2 {
public:	
	CrosslinkPair2(Crosslink &cl1, Crosslink &cl2, double &RRij, int &WBXij, int &WBYij, int &WBZij, int filId1, int filId2);

	void updateEpotPair(double epotPair);
	//void updateContourLength(double contourLength);
	double getContourLength();
	Crosslink* getCrosslink(int i);
	//void updateWB(ublas::vector<int> WB);
	int getWBX();
	int getWBY();
	int getWBZ();
	
	int getFilId(int i);	

	double k2;
	double Leq;
	double r12Abs;
	double f1;


private:
	Crosslink *cl1_;
	Crosslink *cl2_;
	std::vector<int> filIds_; //filId1_= 1 or 2
	double *contourLength_; // length of polymer between crosslinks
	double epotPair_;
	//ublas::vector<int*> WB_; 
	int *WBX_, *WBY_, *WBZ_; // tells you in which box cl2 is located (periodic boundary conditions), WBX=WBY=WBZ=0 -> cl2 in real simulation box
						  // WBX/WBY/WBZ=+-1 cl2 shifted one (virtual) box length in positive/negative x/y/z-direction

};