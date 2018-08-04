#include <iostream>
#include <fstream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>

#include "simulation.hpp"

#include "string"

#include "rngmit.hpp"
// this is polymer

using namespace std;
namespace po = boost::program_options;

extern unsigned long int rng_ia[55];
extern int rng_p,rng_pp;
extern int rng_mod[55];

string inFile;

bool CONTINUE;

int main(int argc , char** argv) {
	// Get Simulation Parameters and options
	double optDt,optRelGamma,optGamma0,optOmega;
	unsigned int optSteps;

	if(argc!=13 && argc!=3){cout << "ERROR: Wrong number of arguments" << "\n"; return 0;}
	else if(argc==3)
	{	
		if (string(argv[1]) == "-con")
		{
			CONTINUE=1;
		}
		else{cout << "ERROR: only continue possible --> -con"<<"\n"; return 0;}
		
		inFile=argv[2];
	}
	else if(argc==13)
	{
		if (string(argv[1]) == "-con")
		{
			CONTINUE=1;
		}
		else if (string(argv[1]) == "-start")
		{
			CONTINUE=0;
		}
		else{cout << "ERROR: -start or -con!" << "\n"; return 0;}
		
		inFile=argv[2];

		for (int i = 3; i < argc; i++) 
		{
			if (string(argv[i]) == "-dt") {
				optDt=stod( argv[i + 1]); i++;
			}
			else if (string(argv[i]) == "-s"){
				optSteps=(int)(stod( argv[i + 1])); i++;
			}
			else if (string(argv[i]) == "-rg"){
				optRelGamma=stod( argv[i + 1]); i++;
			}
			else if (string(argv[i]) == "-g0"){
				optGamma0=stod( argv[i + 1]); i++;
			}
			else if (string(argv[i]) == "-w"){
				optOmega=stod( argv[i + 1]); i++;
			}
			else {cout << "ERROR: wrong order of arguments!" << "\n"; return 0;}
		}
}

	//double optDt = 1e-11; //0.001 gives artifacts, default: 0.00001
	//unsigned int optSteps = 1e8;//4100000

	unsigned int optN = 1;
	int seed = 1;

  	Simulation simulation(optN, optDt, optRelGamma, optGamma0, optOmega, seed, inFile,CONTINUE);

  	simulation.integrate(optSteps);

	cout << "simulation finished" << endl;

	return 0;
}
