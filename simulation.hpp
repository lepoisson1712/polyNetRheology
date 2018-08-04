#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "crosslink.hpp"
#include "crosslinkPair.hpp"
#include "crosslinkTriple.hpp"

#include "crosslinkPair2.hpp"
#include "crosslinkTriple2.hpp"

//#include <random>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/random/uniform_real.hpp>

#include "rngmit.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

class Simulation {
public:
	Simulation(unsigned int, double, double, double, double, int, string, bool);

	void integrate(unsigned int);

	void sampleCoordinates();
	void sampleEnergy();
	void sampleTriples();
	void samplePairs();
	void sampleRR();
	void sampleNetwork();
	void sampleFilLength(vector<double> contLength);
	void sampleAngles(vector<double> theta);
	void sampleClDist(vector<double> clDist);

	double epot2;
	double epot3;
	double msd;

	//####### properties read in from Henry's network files ######//
	int NATOMS; // number of crosslinks
	int FEC,rp,rpp; // ????
	unsigned long int rIA[55];
	double FXX,FYY,FZZ,FXY,FXZ,FYX,FYZ,FZX,FZY; // matrix elements defining box shape
	double BENDING_PF,StretchEps; // ????
	double L, Lc, Lp; // L: length of (quadratic) box, Lc: contour length of polymers, Lp: persistence length
	ublas::vector<double> posX,posY,posZ; // particle positions

	ublas::matrix<int> NN,GB,WBX,WBY,WBZ;; // connectivity table: NN: indices of particles, GB: ghost bonds, WBX/Y/Z: gives (virtual) in which connected particle is
	ublas::matrix<double> RR; // contour length for connections in connectivity table
	//################################################//

	//double k2; // longitudinal spring constant of WLC
	//double Leq; // equilibrium length of WLC

	//double k3; // force constant for three particle force

	double stress; // total stress of the system
	double stress2,stress3; // stress from pair/triple interaction

	double gStorage, gLoss; // storage and loss modulus

	double eTotAveraged; // energy averaged over nAv time steps
	int nAv;

	int nCross;

	string outFolder;
	string outFolderNW;

	string inFile;

	bool CONTINUE;
	
private:

	void createFilament_(int filId, std::vector<int> fil);
	void initialize_();
	void initStreams_();

	void resetForce_();
	void force_();
	void force2_(CrosslinkPair2 &clPair);
	void force3_(CrosslinkTriple2 &clTriple);
	void euler_();
	void applyPeriodicBC_();

	double gaussian_random_number_();
	ublas::vector<double> gaussian_random_vector_();
	ublas::vector<double> gaussian_random_vector2_();

	void mean_square_displacement();

	void calcPullProperties();

	void loadNetwork(string file);

	///!!! NEW !!!
	std::vector<ublas::matrix<int>> createConTabSO(int i);
	void createCrosslinks();
	void createPairList();
	void createTripleList();
	///!!!

	void shear_();

	void correctRR_(int i);

	void compareForces_();



	std::vector<Crosslink> particles_;
	std::vector<CrosslinkPair> particlePairs_;
	std::vector<CrosslinkTriple> particleTriples_;

	///!!! NEW !!!
	std::vector<CrosslinkPair2> particlePairs2_;
	std::vector<CrosslinkTriple2> particleTriples2_;
	///!!!

	int nParticles_;
	double t_;
	double dt_;
	ofstream coordinatesStream_;
	ofstream energyStream_;
	ofstream energyAveragedStream_;
	ofstream pairStream_;
	ofstream tripleStream_;
	ofstream rrStream_;
	ofstream networkStream_;

	ofstream debugStream;

	ofstream filLengthStream_;
	ofstream angleStream_;
	ofstream clDistStream_;
	
	//typedef std::vector<Crosslink> particlevector;
	

	int seed_;
	
	boost::random::mt19937 rng_;
	boost::random::normal_distribution<> rndm_;
	boost::variate_generator< boost::random::mt19937, boost::random::normal_distribution<> >
                    dice_;
    //boost::random::uniform_real_distribution<> rndmU_;
    //boost::variate_generator< boost::random::mt19937, boost::random::uniform_real_distribution<> > diceU_;
	//mt19937 generator_;

	double T_;
	double gamma1_,gamma2_;

	double t_startShear_;

	double relGamma_; // gamma_crosslink/gamma_filament


	double gamma0_; // shear amplitude
	double omega_; // shear frequency

	double cut_; //maximum filament length between 2 crosslinks

};