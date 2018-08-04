#include <fstream>
#include <algorithm>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/assignment.hpp> //for easy initialization of vectors with <<=
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>

#include <ctime>
#include <math.h> /*acos*/
#include <sstream> // for concatenating strings
#include <iomanip>

#include "simulation.hpp"

#include "defs.hpp"

#include <boost/algorithm/string.hpp>


namespace ublas = boost::numeric::ublas;
using namespace std;

Simulation::Simulation(unsigned int nParticles, double dt, double relGamma, double gamma0, double omega, int seed, string file, bool cont): 
	 epot2(0),epot3(0), msd(0),
	 nParticles_(nParticles), t_(0), dt_(dt),
	 rng_(seed),rndm_(0,1),dice_(rng_,rndm_),
	 T_(0.001),gamma1_(1),gamma2_(0.5),
	 NN(1,4),GB(1,4),WBX(1,4),WBY(1,4),WBZ(1,4),RR(1,4),
	 posX(1,3),posY(1,3),posZ(1,3),
	 gamma0_(gamma0), omega_(omega),
	 stress(0), gStorage(0), gLoss(0),
	 t_startShear_(0),
	 nAv(40), eTotAveraged(0),
	 relGamma_(relGamma),
	 nCross(0),inFile(file), CONTINUE(cont),
	 cut_(0.05)
     {
     	string folder = "/usr/scratch1/fischer.andreas/results/runs/";
    	//string folder = "/home/andreas/polymer/results/";
     	outFolder = createTimeDir(folder);

	 	initStreams_();
	 	//initialize_();
	 	
	 	//inFile = "xonly_lc1.000000_L1.000000_xp6_nx1000.0_relGamma1_equilibrated.txt";

	 	//"xonly_lc1.000000_L1.000000_xp6_nx1000.0.txt";
	 	//"xonly_lc1.000000_L1.000000_xp6_nx1000.0_relD1_equilibrated.txt";
	 	//"xonly_lc1.000000_L1.000000_xp6_nx1000.0.txt";
	 	//"xonly_lc1.000000_L1.000000_xp6_nx1000.0_equilibrated.txt"
	 	//"xonly_lc1.000000_L1.000000_xp6_nx1000.0_relD0.001_equilibrated.txt";

	 	loadNetwork(file);

	 	double totFil=0;
	 	for(int i=0; i<NATOMS;i++)
	 	{
	 		for(int j=0; j<4;j++)
	 		{
	 			totFil+=RR(i,j);
	 		}
	 	}
	 	cout << "totFil1: " << totFil << "\n";


	 	if(CONTINUE==false)
	 	{
	 		correctRR_(0); // corrects RR=0 entries to a random entry between 0 and 0.1
	 	}

	 	totFil=0;
	 	for(int i=0; i<NATOMS;i++)
	 	{
	 		for(int j=0; j<4;j++)
	 		{
	 			totFil+=RR(i,j);
	 		}
	 	}
	 	cout << "totFil2: " << totFil << "\n";

	 	createCrosslinks();

	 	createPairList();

	 	createTripleList();

	 	totFil=0;

	 	for(unsigned int i=0; i<particlePairs2_.size(); i++)
	 	{
	 		totFil+=particlePairs2_.at(i).getContourLength();
	 	}
	 	cout << "totFil: " << totFil << "\t" << "nPairs:"<< particlePairs2_.size() << "\n";

	 	int NN2=0, NN3=0, NN4=0;
	 	int nn;
	 	for(int i=0; i<NATOMS;i++)
	 	{
	 		nn=0;
	 		for(int j=0; j<4;j++)
	 		{
	 			nn+=1-GB(i,j);
	 		}
	 		if(nn==2)
	 		{
	 			NN2+=1;
	 		}
	 		else if(nn==3)
	 		{
	 			NN3+=1;
	 		}
	 		else if(nn==4)
	 		{
	 			NN4+=1;
	 		}
	 	}
	 	cout << NN2 << "\t" << NN3 << "\t" << NN4 << "\n";

	 	//sampleCoordinates();
	 	sampleRR();

	 	cout << "continue: " << CONTINUE << "\n";

	 	cout << "gamma0: " << gamma0 << ", omega: " << omega << ", relGamma: " << relGamma_ << ", dt: " << dt << "\n";


	 	if (CONTINUE==false)
	 	{
	 		rngseed(seed);
	 	}
	 	else
	 	{
	 		rnginit(rIA, rpp, rp);
	 	}


	 	stringstream ss;
		ss << outFolder << "network" << "_str" << gamma0_ << "_w" << omega_ << "_dt" << dt_ << "_rg" << relGamma_ ;
		outFolderNW = ss.str();
		string sys_command = "mkdir " + outFolderNW;

		const char * comm = sys_command.c_str();

	 	system(comm);

	 	sampleNetwork();
	 	sampleTriples();
	 }

void Simulation::initStreams_() {
	stringstream ss;
	ss << "_str" << gamma0_ << "_w" << omega_ << "_dt" << dt_ << "_rg" << relGamma_ << "_cut" << cut_;
	string s = ss.str();

	string coordinateFile, energyFile, energyAvFile, rrFile, pairFile, tripleFile, debugFile, filLengthFile, angleFile, clDistFile;//, networkFile;
	coordinateFile=outFolder + "coord" + s +".dat";
	energyFile=outFolder + "energy" + s +".dat";
	energyAvFile= outFolder + "energyAv" + s +".dat";

	rrFile= outFolder + "rr" + s + ".dat";

	pairFile= outFolder + "pairs" + s + ".dat";
	tripleFile= outFolder + "triples" + s + ".dat";

	debugFile= outFolder + "debug" + s + ".dat";

	filLengthFile = outFolder + "filLength" + s + ".dat";

	angleFile = outFolder + "angles" + s + ".dat";

	clDistFile = outFolder + "clDist" + s + ".dat";

	//coordinatesStream_.open (coordinateFile);
  	//coordinatesStream_<< "# id \t r[0] \t r[1] \t r[2]" << endl;

	energyStream_.open(energyFile);	
  	energyStream_<< "# t \t epot2 \t epot3  \t etot \t FXY \t stress2 \t stress3 \t stress \t ghostFilLength \t pairFilLength \t totFilLength" << endl;
  	energyStream_.precision(9);

    debugStream.open(debugFile);

    filLengthStream_.open(filLengthFile);

    angleStream_.open(angleFile);

  	rrStream_.open(rrFile);
  	rrStream_<< "# id \t RR0 \t RR1 \t RR2 \t RR3" << endl;

  	clDistStream_.open(clDistFile);

  	//pairStream_.open(pairFile);
  	//pairStream_ << "# id1 \t id2 \t l12 \t k2 \t Leq \t r12Abs" << endl;
   	
   	tripleStream_.open(tripleFile);
  	tripleStream_ << "# id1 \t id2 \t id3 \t angle" << endl; 	
  	
  	//energyAveragedStream_.open (energyAvFile);
  	//energyAveragedStream_ << "# t  \t etot \t FXY"<< endl;
  	//energyAveragedStream_.precision(9);
}


void Simulation::integrate(unsigned int nSteps) {
	ublas::vector<double> r(3,0);


	for(size_t i=1;i<=nSteps;i++) {

		/*if(i==1){
			t_startShear_=t_;
		}*/

		//if(i>400000){
		shear_();
			//cout << FXY << "\n";
		//}
		stress=0;
		stress2=0;
		stress3=0;

        force_(); // force calculation
       
        if (relGamma_>0){
        	calcPullProperties();
        }
        //cout << RR(10,3) << "\n";
        //if(i>50000){
        stress=stress2+stress3;
        //}

        sampleEnergy();

        vector<double> elong,springC,contLength,fPull,theta,clDist; 
        ublas::vector<double> r1(3,0),r2(3,0);
        ublas::vector<int> WB(3,0);
		ublas::vector<double> FtimesWB(3,0);
        double r12Abs;
        if(i%1 == 0)
        {
        	for( unsigned int i=0; i<particlePairs2_.size(); i++)
        	{	
        		//cout << particlePairs2_.size();
        		CrosslinkPair2 clpair = particlePairs2_.at(i);
        		elong.push_back(clpair.r12Abs-clpair.Leq);
        		springC.push_back(clpair.k2);
        		contLength.push_back(clpair.getContourLength());

        		r1=clpair.getCrosslink(1)->getPosition();
        		r2=clpair.getCrosslink(2)->getPosition();
        		WB[0]=clpair.getWBX(); WB[1]=clpair.getWBY(); WB[2]=clpair.getWBZ();
				FtimesWB[0]=FXX*WB[0]+FXY*WB[1]; FtimesWB[1]=FYY*WB[1]; FtimesWB[2]=FZZ*WB[2]; // F*WB (F=matrix, WB=vector)
				r2 += FtimesWB;  
        		r12Abs=ublas::norm_2(r1-r2);
        		clDist.push_back(r12Abs);
	        }

	        for( unsigned int i=0;i<particleTriples2_.size(); i++)
	        {
	        	CrosslinkTriple2 cltriple = particleTriples2_.at(i);
	        	fPull.push_back(abs(cltriple.getfPull()));
	        	if(cltriple.getDummy()==0)
	        	{
	        		theta.push_back(cltriple.getAngle());
	        	}
	        }

      		debugStream << t_ << "\t" << *max_element(elong.begin(),elong.end()) << "\t" << *max_element(springC.begin(),springC.end()) << "\t" <<
      	    *min_element(contLength.begin(),contLength.end()) << "\t" << *max_element(contLength.begin(),contLength.end()) <<"\t" << *max_element(fPull.begin(),fPull.end()) << "\t" << epot2 << "\n";
    	}
        

        euler_(); // integrate equations of motion;

        t_ += dt_;
        //applyPeriodicBC_();
        
        if (i%1000 == 0)
        {
        	//sampleCoordinates();
        	sampleRR();
        	//sampleTriples();
        	//samplePairs();
        	sampleNetwork();
        	sampleFilLength(contLength);
        	sampleAngles(theta);
        	sampleClDist(clDist);
        	cout << t_ << "\n";	
        	//cout << nCross << "\n";
        	//cout << FXY << "\n";
        	//compareForces_();

        }

        /*for(int i=0;i<NATOMS;i++)
        {
        	for(int j=0;j<4;j++)
        		if(RR(i,j)<0){cout << "some negative RR" << "\n";}
        }*/
        // check if arclength is conserved
        /*if (i%100 == 0){
        double RRSUM=0;
	        for(int i=0; i<NATOMS;i++)
	        {
	        	for(int j=0; j<4;j++)
	        	{
	        		RRSUM+=RR(i,j);
	        	}
	        }
	        cout << RRSUM << "\n";
	    }*/
	}

}


void Simulation::resetForce_() {
	ublas::vector<double> f(3,0);
	for (std::vector<Crosslink>::iterator particlei = particles_.begin(); particlei != particles_.end(); ++particlei)
	{
		particlei->updateForce(f);
		particlei->updateFilForce(f,1);
		particlei->updateFilForce(f,2);
		//particlei->filament1_.updateForce(f);
		//particlei->filament2_.updateForce(f);
	}
}


void Simulation::force_() {
	resetForce_(); //sets force to zero;

	epot2=0;
	epot3=0;

	//for (std::vector<Crosslink>::iterator particlei = particles_.begin(); particlei != particles_.end(); ++particlei)
	for ( unsigned int i=0; i<particlePairs2_.size(); i++)
	{	

		force2_(particlePairs2_.at(i));

		
	}
	for ( unsigned int i=0; i<particleTriples2_.size(); i++)
	{	
		force3_(particleTriples2_.at(i));

	}
}


void Simulation::force2_(CrosslinkPair2 &clPair){
	
	ublas::vector<double> r1(3,0),r2(3,0);
	ublas::vector<double> f1(3,0),f2(3,0);
	ublas::vector<double> r12(3,0);
	double r12Abs=0;
	double epotPair=0;
	ublas::vector<int> WB(3,0);
	ublas::vector<double> FtimesWB(3,0);
	int filId1, filId2;

	double l12; // length of polymer between cl 1 and 2
	double k2; // spring constant
	double Leq; //equilibrium length

	Crosslink *cl1 = clPair.getCrosslink(1);
	Crosslink *cl2 = clPair.getCrosslink(2);
	//int filId = clPair.getFilId();

	r1=cl1->getPosition();
	r2=cl2->getPosition();
	//!!!!periodic boundary conditions!!!!
	WB[0]=clPair.getWBX(); WB[1]=clPair.getWBY(); WB[2]=clPair.getWBZ();
	
	FtimesWB[0]=FXX*WB[0]+FXY*WB[1]; FtimesWB[1]=FYY*WB[1]; FtimesWB[2]=FZZ*WB[2]; // F*WB (F=matrix, WB=vector)
	//cout << FtimesWB << "\n";
	r2 += FtimesWB;  
	//!!!!!//

	r12 = r2-r1;
	r12Abs=ublas::norm_2(r12);
	//calculate k2 and Leq
	l12=clPair.getContourLength();
	k2=90./pow(l12,4);
	Leq=l12*(1-l12/6.);

	clPair.k2=k2;
	clPair.r12Abs=r12Abs;
	clPair.Leq=Leq;

	// pair force: harmonic spring with spring constant k2
	f1=k2*(r12Abs-Leq)*r12/r12Abs; 	
	f2=-f1;

	clPair.f1=ublas::norm_2(f1);


	filId1=clPair.getFilId(1);
	filId2=clPair.getFilId(2);

	cl2->incrForce(f2,filId2);	

	//cl2->incrForceFil(f,filId);
	cl1->incrForce(f1,filId1);
	//cl1->incrForceFil(-1*f,filId);
	
	epotPair=0.5*k2*pow(r12Abs-Leq,2);
	clPair.updateEpotPair(epotPair);
	epot2+=epotPair;

	//cout << "k2: " << k2 << "  Leq: " << Leq << "    deltaRabs: " << deltaRabs << "    f: " << f << "\n";

	// compute stress sigma_xy=- (r2[0]-r1[0])*f1[1]	
	stress2+=r12[0]*f1[1];
}

void Simulation::force3_(CrosslinkTriple2 &clTriple){
	
	int iDummy;
	double angle;//,r12abs,r13abs,r23abs;
	ublas::vector<double> f1(3,0),f2(3,0),f3(3,0);
	ublas::vector<double> r1(3,0),r2(3,0),r3(3,0);
	ublas::vector<double> r12(3,0),r23(3,0),e12(3,0),e23(3,0);
	double r12abs, r23abs;
	double epotTriple=0;
	ublas::vector<int> WB12(3,0),WB23(3,0);
	ublas::vector<double> FtimesWB12(3,0);
	ublas::vector<double> FtimesWB23(3,0);

	double l12,l23; // contour lengths between 1,2 and 2,3
	double k3; // prefactor of force

	int filId1,filId2,filId3;

	double pref;
	double angTol=M_PI/100.;

	Crosslink *cl1 = clTriple.getCrosslink(1);
	Crosslink *cl2 = clTriple.getCrosslink(2);
	Crosslink *cl3 = clTriple.getCrosslink(3);

	//int filId = clTriple.getFilId();
	

	iDummy=clTriple.getDummy();
	if(iDummy==0) // 3-cl force only if no dummies in triple
	{
		r1=cl1->getPosition();
		r2=cl2->getPosition();
		r3=cl3->getPosition();

		///!!!! periodic boundary conditions !!!!
		WB12[0]=clTriple.getWBX(12); WB12[1]=clTriple.getWBY(12); WB12[2]=clTriple.getWBZ(12);
		WB23[0]=clTriple.getWBX(23); WB23[1]=clTriple.getWBY(23); WB23[2]=clTriple.getWBZ(23);

		FtimesWB12[0]=FXX*WB12[0]+FXY*WB12[1]; FtimesWB12[1]=FYY*WB12[1]; FtimesWB12[2]=FZZ*WB12[2]; // F*WB (F=matrix, WB=vector)
		FtimesWB23[0]=FXX*WB23[0]+FXY*WB23[1]; FtimesWB23[1]=FYY*WB23[1]; FtimesWB23[2]=FZZ*WB23[2]; // F*WB (F=matrix, WB=vector)


		r1+=FtimesWB12; //was r2!!
		r3+=FtimesWB23;
		///!!!

		r12=r2-r1;
		r23=r3-r2;
		//cout << ublas::inner_prod(r12,r23) << "\n";
		r12abs=ublas::norm_2(r12);
		r23abs=ublas::norm_2(r23);

		angle= acos(ublas::inner_prod(r12,r23)/(r12abs*r23abs));
		//cout << acos(1) << "\n";

		l12=clTriple.getContourLength(12);
		l23=clTriple.getContourLength(23);

		k3=1.5/(l12+l23);

		if(angle>angTol)
		{  // need this condition, because for angle=0, divison by 0 (sin(0)=0)
			pref=angle/sin(angle);
		}
		else //taylor expand sin(angle) for small angles
		{
			pref=1./(1-1./6.*angle*angle+1./120.*angle*angle*angle*angle-1./5040.*angle*angle*angle*angle*angle*angle);
		}


		e12=r12/r12abs;
		e23=r23/r23abs;

		f1=(-e23+cos(angle)*e12)/r12abs;
		f3=(e12-cos(angle)*e23)/r23abs;
		//f2=k3*2*pref/r12abs/r23abs*(r23-r12-cos(angle)*(r12/r12abs-r23/r23abs));
		f2=(e23-cos(angle)*e12)/r12abs+(e23*cos(angle)-e12)/r23abs; //10.03.2016

		f1=f1*k3*2*pref;
		f2=f2*k3*2*pref;
		f3=f3*k3*2*pref;
		//cout << f1 << "\t" << f2 << "\t" << f3 << "\n";

		// new forces (!!!WRONG!!!)
		//f1=-k3*(r12/r12abs/pow((r12abs+r23abs),2)*angle*angle-2*angle/sin(angle)/r12abs*(-r23/r23abs+cos(angle)*r12/r12abs));
		//f3=-k3*(r23/r23abs/pow((r12abs+r23abs),2)*angle*angle-2*angle/sin(angle)/r23abs*(r12/r12abs-cos(angle)*r23/r23abs));
		//f2=-k3*((r23/r23abs-r12/r12abs)/pow((r12abs+r23abs),2)*angle*angle-2*angle/sin(angle)/r12abs/r23abs*(r23-r12-cos(angle)*(r12/r12abs-r23/r23abs)));

		filId1=clTriple.getFilId(1);
		filId2=clTriple.getFilId(2);
		filId3=clTriple.getFilId(3);

		cl1->incrForce(f1,filId1);
		//cl1->incrForceFil(f1,filId);
		cl2->incrForce(f2,filId2);
		//cl2->incrForceFil(f2,filId);
		cl3->incrForce(f3,filId3);
		//cl3->incrForceFil(f3,filId);

		//epotTriple=k3/(r12abs+r23abs)*angle*angle;
		epotTriple=k3*angle*angle;

		clTriple.updateEpotTriple(epotTriple);
		clTriple.updateAngle(angle);
		epot3+=epotTriple;
		//cout << angle << "\n";

		// compute stress sigma_xy=-((r1[0]-r2[0])*f1[1]+(r3[0]-r2[0])*f3[1]), WRONG!! virial stress is zero for 3-body forces that depend on bond angle only
		stress3+= (r2[0]-r1[0])*f1[1]+(r2[0]-r3[0])*f3[1]; 
	}
}


void Simulation::euler_() {
	ublas::vector<double>  noise(3,0);
	double noise2;
	double ds;
	double sigma1;
	double sigma2=sqrt(2*dt_);//*T_/gamma1_);
	double l12, l23;
	ublas::vector<double>  r1(3,0),r2(3,0),r3(3,0);
	double r12,r23;
	int iDummy;
	//double cut=0.05;

	
	if(relGamma_!=-1)
	{	
		int d=0;
		for (std::vector<CrosslinkTriple2> ::iterator triplei = particleTriples2_.begin(); triplei != particleTriples2_.end(); triplei++)
		{	
			iDummy=triplei->getDummy();
			
			//if(triplei->getDummy() == 0){
			//noise2=dice_();
			noise2=gaussian_random_number_();
			ds=dt_*triplei->getfPull()+sigma2*noise2;
			//cout << ds << "\n";

			// FOR DEBUGGING ONLY
			r1=triplei->getCrosslink(1)->getPosition();
			r2=triplei->getCrosslink(2)->getPosition();
			r3=triplei->getCrosslink(3)->getPosition();
			r12=ublas::norm_2(r1-r2);
			r23=ublas::norm_2(r3-r2);
			// FOR DEBUGGING ONLY
			l12=triplei->getContourLength(12);
			l23=triplei->getContourLength(23);
			//cout << l12 << "\t" << l23 << "\n";
			//if( (l12+ds >= r12) && (l23-ds >= r23))
			//{
		
			if(iDummy==0)
			{
				//bool a = (l23-ds>cut_) && (l12+ds>=cut_);

				//cout << triplei->getCrosslink(1)->getId() << "\t"<< triplei->getCrosslink(2)->getId()<< "\t" << triplei->getCrosslink(3)->getId() << "\t" << a << "\t" << l12 << "\t" << l23 << "\t" << ds << "\n";
				if( (l12+ds>cut_) && (l23-ds>cut_) ) // contour length between two crosslinks must be > cut_
				{
					triplei->incrContourLength(ds,12);
					triplei->incrContourLength(-ds,23);	
				}
			}
			if(iDummy==1)
			{	
				//bool a = (l23-ds>cut_) && (l12+ds>=0);
				//cout << triplei->getCrosslink(1)->getId() << "\t"<< triplei->getCrosslink(2)->getId()<< "\t" << triplei->getCrosslink(3)->getId() << "\t" << a << "\t" << l12 << "\t" << l23 << "\t" << ds << "\n";
				if((l23-ds>cut_) && (l12+ds>=0) )
				{	
					//cout << l12+ds << "\n";
					//cout << l23-ds << "\n";
					triplei->incrContourLength(ds,12);
					triplei->incrContourLength(-ds,23);
					//cout << triplei->getContourLength(12) << "\n";
					//cout << triplei->getContourLength(23) << "\n";
				}
			}

			if(iDummy==3)
			{
				if( (l12+ds>cut_) && (l23-ds>=0))
				{
					triplei->incrContourLength(ds,12);
					triplei->incrContourLength(-ds,23);	
				}
			}
			

			l12=triplei->getContourLength(12);
			l23=triplei->getContourLength(23);

			if(l12 < 0 || l23 <0){cout << "some RR smaller cut" << "\t" << ds << "\n";}
			if( (l12 < r12) || (l23 < r23))
			{
				nCross++;
			}
		//}
			
		}
	}

	int id;
	ublas::vector<double> dr,r,rN;
	double drN;
	int check;

	double rG;

	for (std::vector<Crosslink>::iterator particlei = particles_.begin(); particlei != particles_.end(); particlei++)
	{	
		check=1;
		noise=gaussian_random_vector2_();

		if(relGamma_==-1){rG=1;}
		else{rG=relGamma_;}

		sigma1=sqrt(2*rG*dt_);

		dr=	rG*dt_*particlei->getForce()+sigma1*noise;	

		//cout << sigma1*noise << "\t" << noise << "\n";
		//cout << relGamma_*dt_*particlei->getForce() << "\n";
		// check whether step doesn't bring cls too close
		id=particlei->getId();
		r=particlei->getPosition();
		for(int i=0;i<4;i++)
		{
			if(GB(id,i)==0)
			{
				rN=particles_.at(NN(id,i)).getPosition();
				drN=ublas::norm_2(rN-(r+dr));
				//cout << drN << "\n";
				if(drN<0.001)
				{
					check=0;
					
				}
			}
		}
		if(check==1)
		{
			particlei->incrPosition(dr);
		}
		
		
		//particlei->incrPosition(sigma1*noise);
		//DeltaS1=1/gamma*Force1+noise1;
		//DeltaS2=1/gamma*Force2+noise2;
		//DeltaR=0.5*(DeltaS1*particlei->filament1_.getPullDirection()+DeltaS2*particlei->filament2_.getPullDirection());
		//particlei->incrPosition(deltaR);		
	}
}

void Simulation::applyPeriodicBC_(){
	ublas::vector<double> r(3,0);
	ublas::matrix<int> conTabSelf;
	ublas::matrix<int> conTabOther;
	int k,l;
	int deltaWBX, deltaWBY, deltaWBZ;
	ublas::vector<double> eX(3), eY(3), eZ(3);
	eX[0]=1;eX[1]=0;ex[2]=0; eY[0]=0;eY[1]=1;eY[2]=0; eZ[0]=0;eZ[1]=0;eZ[2]=1;  
	double dxPbc, dyPbc, dzPbc;
	ublas::vector<double> drPbc(3,0);

	for (std::vector<Crosslink>::iterator particlei = particles_.begin(); particlei != particles_.end(); particlei++)
	{	
		deltaWBX=0; deltaWBY=0; deltaWBZ=0;
		dxPbc=0; dyPbc=0; dzPbc=0;
		r=particlei->getPosition();
		if(r[0] > 0.5*L || r[0] < -0.5*L || r[1] > 0.5*L || r[1] < -0.5*L || r[2] > 0.5*L || r[2] < -0.5*L)
		{
			conTabSelf=particlei->getConTabSelf();
			conTabOther=particlei->getConTabOther();
			if(r[0]>0.5*L){
				deltaWBX=1;
				dxPbc= -L;
			}
			else if(r[0]<-0.5*L){
				deltaWBX=-1;
				dxPbc=L;
			}
			if(r[1]>0.5*L){
				deltaWBY=1;
				dyPbc= -L;
			}
			else if(r[1]<-0.5*L){
				deltaWBY=-1;
				dyPbc=L;
			}	
			if(r[2]>0.5*L){
				deltaWBZ=1;
				dzPbc = -L;
			}
			else if(r[2]<-0.5*L){
				deltaWBZ=-1;
				dzPbc = L;
			}	
			drPbc <<= dxPbc, dyPbc, dzPbc;
			particlei->incrPosition(drPbc);

			for(int i=0; i<conTabSelf.size1(); i++){
				k=conTabSelf(i,0); 
				l=conTabSelf(i,1);
				WBX(k,l) -= deltaWBX;
				WBY(k,l) -= deltaWBY;
				WBZ(k,l) -= deltaWBZ;
				if(abs(WBX(k,l))>1 || abs(WBY(k,l)) > 1 || abs(WBZ(k,l)) > 1){
					cout << "some abs(WB) > 2, probably something went wrong!!" <<"\n";
				}
			}

			for(int i=0; i<conTabOther.size1(); i++){
				k=conTabOther(i,0);
				l=conTabOther(i,1);
				WBX(k,l) += deltaWBX;
				WBY(k,l) += deltaWBY;
				WBZ(k,l) += deltaWBZ;
				if(abs(WBX(k,l))>1 || abs(WBY(k,l)) > 1 || abs(WBZ(k,l)) > 1){
					cout << "some abs(WB) > 2, probably something went wrong!!" <<"\n";
				}
			}
		}
	}
}

void Simulation::sampleCoordinates() {
	ublas::vector<double> r(3);//, forceFil1(3), forceFil2(3), force(3), pullDir1(3), pullDir2(3);
	int id;//,filId1,filId2,filClId11,filClId12,filClId21,filClId22;
	//double pullForce1,pullForce2;
	//Crosslink* particlei;
	coordinatesStream_ << "# t=" << t_ << "##################" << "\n"; 
	for (std::vector<Crosslink>::iterator particlei = particles_.begin(); particlei != particles_.end(); particlei++)
	{		
		id = particlei->getId();
		r = particlei->getPosition();
		//coordinatesStream_ << id << "\t" << r[0] << "\t" << r[1] << "\t" << r[2] << "\n";
		coordinatesStream_ << r[0] << "\t" << r[1] << "\t" << r[2] << "\n";
		/*if(abs(r[0])>0.5 || abs(r[1])>0.5 || abs(r[2])>0.5){
			cout << "some coordinate outside box --> problem with periodic boundaries" << "\n";
		}*/
	}
}


void Simulation::sampleEnergy() {
	double ghostFilLength=0, pairFilLength=0, totFilLength=0;
	for(int i=0;i<NATOMS;i++)
	{
		for(int j=0;j<4;j++)
		{
			if(GB(i,j)==1)
			{
				ghostFilLength += RR(i,j);
			}
			else
			{
				pairFilLength += RR(i,j);
			}
		}
	}
	pairFilLength/=2.;

	totFilLength=ghostFilLength+pairFilLength;

	energyStream_ << t_ << "\t" << epot2 << "\t" << epot3 << "\t" << epot2+epot3 << "\t" << FXY << "\t" << stress2 << "\t" << stress3 << "\t" << stress << "\t" << ghostFilLength << "\t" << pairFilLength << "\t" << totFilLength << endl;
	//if(i % nAv == 0){
	//energyAveragedStream_ << t_ << "\t" << epot2+epot3 << "\t" << FXY << endl;
	//}
}

void Simulation::sampleTriples() {
	int id1,id2,id3;
	double angle;
	int iDummy;

	for (std::vector<CrosslinkTriple2>::iterator triplei = particleTriples2_.begin(); triplei != particleTriples2_.end(); triplei++)
	{
		iDummy=triplei->getDummy();

		//if(iDummy==0)
		//{
			id1=triplei->getCrosslink(1)->getId();id2=triplei->getCrosslink(2)->getId();id3=triplei->getCrosslink(3)->getId();
			angle=triplei->getAngle();
			tripleStream_ << t_ << "\t" << id1 << "\t" << id2 << "\t" << id3 << "\t" << angle << "\t" << triplei->getContourLength(12) << "\t" << triplei->getContourLength(23) << endl;
		//}	
	}
}

void Simulation::samplePairs() {
	int id1,id2;
	pairStream_ << "# t=" << t_ << "##################" << "\n"; 
	for (std::vector<CrosslinkPair2>::iterator pairi = particlePairs2_.begin(); pairi != particlePairs2_.end(); pairi++)
	{
		id1=pairi->getCrosslink(1)->getId();id2=pairi->getCrosslink(2)->getId();
		pairStream_ << id1 << "\t" << id2 << "\t" << pairi->getContourLength() << "\t" << pairi->k2 << "\t" << pairi->r12Abs - pairi->Leq << "\t"
		<< pairi->f1 << "\t" << pairi->Leq << "\t" << pairi-> r12Abs << endl;
	}
}

void Simulation::sampleRR() {
	rrStream_ << "#t=" << t_ << "##################" << "\n"; 
	for(int i=0;i<NATOMS;i++)
	{
		rrStream_ << RR(i,0) << "\t" << RR(i,1) << "\t" << RR(i,2) << "\t" << RR(i,3) << endl;
	}
}

void Simulation::sampleFilLength(vector<double> contLength)
{	
	filLengthStream_ << "#t=" << t_ << "###############" << endl;
	double binSize=0.005;
	int nBins=(int) (1./binSize);
	//cout << nBins << "\n";
	double xHist [nBins];
	for(int i=0;i<nBins;i++)
	{
		xHist[i] = (i+1)*binSize;
	}
	int yHist [nBins];
	for(int i=0;i<nBins;i++)
	{
		yHist[i]=0;
	}
	int n;

	for (int i=0; i<contLength.size(); i++)
	{
		n=floor(contLength[i]/binSize);
		yHist[n]++;
	}

	for(int i=0;i<nBins;i++)
	{
		filLengthStream_ << xHist[i] << "\t" << yHist[i] << endl;
	}

}

void Simulation::sampleAngles(vector<double> theta)
{
	angleStream_ << "#t=" << t_ << "###############" << endl;


	double binSize = 0.02;
	int nBins = (int) (1.6/binSize);

	double xHist [nBins];
	for(int i=0;i<nBins;i++)
	{
		xHist[i] = (i+1)*binSize;
	}
	int yHist [nBins];
	for(int i=0;i<nBins;i++)
	{
		yHist[i]=0;
	}
	int n;

	for (int i=0; i<theta.size(); i++)
	{
		n=floor(theta[i]/binSize);
		yHist[n]++;
	}

	for(int i=0;i<nBins;i++)
	{
		angleStream_ << xHist[i] << "\t" << yHist[i] << endl;
	}
}

void Simulation::sampleClDist(vector<double> clDist)
{
	clDistStream_ << "#t=" << t_ << "###############" << endl;

	double binSize = 0.005;
	int nBins = (int) (1./binSize);

	double xHist [nBins];
	for(int i=0;i<nBins;i++)
	{
		xHist[i] = (i+1)*binSize;
	}
	int yHist [nBins];
	for(int i=0;i<nBins;i++)
	{
		yHist[i]=0;
	}
	int n;

	for (int i=0; i<clDist.size(); i++)
	{
		n=floor(clDist[i]/binSize);
		yHist[n]++;
	}

	for(int i=0;i<nBins;i++)
	{
		clDistStream_ << xHist[i] << "\t" << yHist[i] << endl;
	}
}



void Simulation::sampleNetwork() 
{
	string outFile;

	stringstream ss;

	//inFile = "xonly_lc1.000000_L1.000000_xp6_nx1000.0_relGamma1_equilibrated.txt";

	ss.precision(9);

	size_t txtPos = inFile.find(".txt");
	string nw_str = inFile.substr(0,txtPos);

	ss << nw_str << "_t" << fixed << t_;
	ss.precision(4);
	ss << defaultfloat << "_str" << gamma0_ << "_w" << defaultfloat << omega_ << "_dt" << dt_ << "_rg" << relGamma_ ;
	string s = ss.str();

	outFile= outFolderNW + "/" +  s + ".dat";
	networkStream_.open(outFile);

	//networkStream_ << "###############################" << endl; 
	//networkStream_ << defaultfloat << "#t=" << t_  << endl; 

	networkStream_ << "#NATOMS=" << dec << NATOMS << endl;
	networkStream_ << "#FEC=1" << endl;

	networkStream_<<  "#FXX=" << std::hexfloat << FXX << endl;
	networkStream_<< "#FYY=" <<  FYY << endl;
	networkStream_<< "#FZZ=" <<  FZZ << endl;
	networkStream_<< "#FXY=" <<  FXY << endl;
	networkStream_<< "#FXZ=" <<  FXZ << endl;
	networkStream_<< "#FYZ=" <<  FYZ << endl;
	networkStream_<< "#FYX=" <<  FYX << endl;
	networkStream_<< "#FZX=" <<  FZX << endl;
	networkStream_<< "#FZY=" <<  FZY << endl;
	networkStream_<< "Lp=" << dec << Lp << endl;
	networkStream_<< "#BENDING_PF=0x1.8p+0" << endl;
	networkStream_<< "#StretchEps=0x1.0c6f7a0b5ed8dp-20" << endl;
	networkStream_<< "#rp=" << rng_p << endl;
	networkStream_<< "#rpp=" << rng_pp << endl;

	networkStream_<< "#rIA=";
	for(int i=0; i<54;i++)
	{
		networkStream_<< rng_ia[i] << " ";
	}
	networkStream_<< rng_ia[54] << endl; 


	ublas::vector<double> r(3,0);

	for (std::vector<Crosslink>::iterator particlei = particles_.begin(); particlei != particles_.end(); particlei++)
	{		
		r = particlei->getPosition();
		networkStream_ << hexfloat << r[0] << "\t" << r[1] << "\t" << r[2] << "\n";
	}

	for(int i=0;i<NATOMS;i++)
	{
		networkStream_ << dec << NN(i,0) << " " << NN(i,1) << " " << NN(i,2) << " " << NN(i,3) << endl;
	}

	for(int i=0;i<NATOMS;i++)
	{
		networkStream_ << dec << GB(i,0) << " " << GB(i,1) << " " << GB(i,2) << " " << GB(i,3) << endl;
	}

	for(int i=0;i<NATOMS;i++)
	{
		networkStream_ << hexfloat << RR(i,0) << " " << RR(i,1) << " " << RR(i,2) << " " << RR(i,3) << endl;
	}

	for(int i=0;i<NATOMS;i++)
	{
		networkStream_ << dec << WBX(i,0) << " " << WBX(i,1) << " " << WBX(i,2) << " " << WBX(i,3) << endl;
	}

	for(int i=0;i<NATOMS;i++)
	{
		networkStream_ << dec << WBY(i,0) << " " << WBY(i,1) << " " << WBY(i,2) << " " << WBY(i,3) << endl;
	}

	for(int i=0;i<NATOMS;i++)
	{
		networkStream_ << dec << WBZ(i,0) << " " << WBZ(i,1) << " " << WBZ(i,2) << " " << WBZ(i,3) << endl;
	}

	networkStream_.close();


}




/* void Simulation::gaussian_random_number_(){
	ublas::vector<double> randVec(3,0);
    double f1,f2,p;

    //uniform_real_distribution<double> dis(0.0, 1.0);
	
    //f1=dis(generator_);
    //f2=dis(generator_);
    
    f1=dice_();
    f2=dice_();
    p=sqrt(-2.*log(f1));

    randVec[0]=p*cos(2.*M_PI*f2);
    randVec[1]=p*sin(2.*M_PI*f2);

    //rand_no1_=p*cos(2.*M_PI*f2);
    //rand_no2_=p*sin(2.*M_PI*f2);
    // randomStream << rand_no[1] << "\n";
    // randomStream << rand_no[2] << "\n";
}*/

ublas::vector<double> Simulation::gaussian_random_vector_(){
	ublas::vector<double> randVector(3,0);
	randVector[0]=dice_();
	randVector[1]=dice_();
	randVector[2]=dice_();
	return randVector;
}

 ublas::vector<double> Simulation::gaussian_random_vector2_(){

 	ublas::vector<double> randVec(3,0);

     double f1,f2,f3,f4,p1,p2;
     

     //f1=dice_();
     f1=rngmit;
     //f2=dice_();
     f2=rngmit;

     p1=sqrt(-2.*log(f1));

     f3=rngmit;
     //f2=dice_();
     f4=rngmit;

     p2=sqrt(-2.*log(f3));

     randVec[0]=p1*cos(2.*M_PI*f2);

     randVec[1]=p1*sin(2.*M_PI*f2);

     randVec[2]=p2*cos(2.*M_PI*f4);

     //cout << randVec << "\n";

     return randVec;

 

     //rand_no1_=p*cos(2.*M_PI*f2);

     //rand_no2_=p*sin(2.*M_PI*f2);

     // randomStream << rand_no[1] << "\n";

     // randomStream << rand_no[2] << "\n";

 }

 double Simulation::gaussian_random_number_()
 {
 	double f1,f2,p;
 	double randNo;

 	f1=rngmit;
 	f2=rngmit;

 	p=sqrt(-2.*log(f1));

 	randNo=p*cos(2.*M_PI*f2);

 	return randNo;
 }


void Simulation::mean_square_displacement(){
	ublas::vector<double> r(3),deltaR(3),r0(3,1);
	double msd1;
	int i=0;
	for (std::vector<Crosslink>::iterator particle = particles_.begin(); particle != particles_.end(); particle++)
	{
		i++;
		r = particle->getPosition();
		deltaR = r-i*r0;
		msd1=pow(ublas::norm_2(deltaR),2);
		msd+=msd1;
	}		
	msd=msd/3000;
}


void Simulation::calcPullProperties(){
	ublas::vector<double> r1(3,0),r2(3,0),r3(3,0);
	ublas::vector<double> r21(3,0),r23(3,0),r13(3,0);
	ublas::vector<double> e21(3,0),e23(3,0),e13(3,0);
	ublas::vector<int> WB12(3,0),WB23(3,0);
	ublas::vector<double> FtimesWB12(3,0);
	ublas::vector<double> FtimesWB23(3,0);
	ublas::vector<double> f2(3,0);	
	double f21,f23, f13;
	double pullforce;
	int iDummy;
	int filId2;
	ublas::vector<double> fOtherFil(3,0);
	
	for (std::vector<CrosslinkTriple2> ::iterator clTriple = particleTriples2_.begin(); clTriple != particleTriples2_.end(); clTriple++)
	{	
		Crosslink *cl1 = clTriple->getCrosslink(1);
		Crosslink *cl2 = clTriple->getCrosslink(2);
		Crosslink *cl3 = clTriple->getCrosslink(3);
	
		r1=cl1->getPosition();
		r2=cl2->getPosition();
		r3=cl3->getPosition();

		iDummy=clTriple->getDummy();

		//if(iDummy==0){
		///!!!! periodic boundary conditions !!!!
		WB12[0]=clTriple->getWBX(12); WB12[1]=clTriple->getWBY(12); WB12[2]=clTriple->getWBZ(12);
		WB23[0]=clTriple->getWBX(23); WB23[1]=clTriple->getWBY(23); WB23[2]=clTriple->getWBZ(23);

		FtimesWB12[0]=FXX*WB12[0]+FXY*WB12[1]; FtimesWB12[1]=FYY*WB12[1]; FtimesWB12[2]=FZZ*WB12[2]; // F*WB (F=matrix, WB=vector)
		FtimesWB23[0]=FXX*WB23[0]+FXY*WB23[1]; FtimesWB23[1]=FYY*WB23[1]; FtimesWB23[2]=FZZ*WB23[2]; // F*WB (F=matrix, WB=vector)


		r1+=FtimesWB12; // was r2!! 
		r3+=FtimesWB23;


		// !!! WHEN DUMMY TRIPLE, TAKE CARE!!
		
		r21=r1-r2;
		r23=r3-r2;

		if(iDummy==0)
		{ 	
			e21=r21/ublas::norm_2(r21);
			e23=r23/ublas::norm_2(r23);
		}
		// set projection direction for dummy crosslinks
		else if(iDummy==1) 
		{
			e23=r23/ublas::norm_2(r23);
			e21=-e23;
		}
		else if(iDummy==3)
		{
			e21=r21/ublas::norm_2(r21);
			e23=-e21;
		}
		// !!! CHANGE TO FORCE FROM OTHER FILAMENT !!!
		filId2=clTriple->getFilId(2);

		if(filId2==1)
		{
			fOtherFil=(cl2->getFilForce(2)-cl2->getFilForce(1))*0.5;
			//cout << fOtherFil << "\t" << cl2->getFilForce(2) << "\t" << cl2->getFilForce(1) << "\t" << cl2->getForce() << "\n";

		}
		else if(filId2==2)
		{
			fOtherFil=(cl2->getFilForce(1)-cl2->getFilForce(2))*0.5;
		}
		else
		{
			cout << "filId Error3" << "\n";
		}
		//fOtherFil=cl2->getForce();
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		r13=r3-r1;
		e13=r13/ublas::norm_2(r13);

		f13=ublas::inner_prod(fOtherFil,e13);

		f21=ublas::inner_prod(fOtherFil,e21);  // project pull force
		f23=ublas::inner_prod(fOtherFil,e23);


		pullforce=f13;

		/*if(f23>f21){
			//pullforce=f23*1;
			pullforce=ublas::norm_2(fOtherFil);
		}
		else{
			//pullforce=f21*(-1);
			pullforce=ublas::norm_2(fOtherFil)*(-1);
		}*/

		clTriple->updatefPull(pullforce);
		///!! case f23==f21 !! missing
		/*if(iDummy==0)
		{
			cout << pullforce << "\n";
		}*/
	//}
	}

}



void Simulation::loadNetwork(string file){	
	
	// read Lc, L from filename
	size_t lcPos = file.find("lc");
	string Lc_str = file.substr(lcPos+2,8);
	Lc = stod(Lc_str); // convert string to double
	cout << Lc << "\n";

	size_t LPos = file.find("L");
	string L_str = file.substr(LPos+1,8);
	L = stod(L_str);
	cout << L << "\n";	

	if(CONTINUE==1){
		size_t tPos = file.find("_t");
		size_t strPos = file.find("_str");
		string t_str = file.substr(tPos+2,strPos-tPos);
		t_=stod(t_str);
		cout << t_ << "\n";
	}

	// read parameters from file
	ifstream infile(file);
	string str;


	str=readNumber(infile);
	NATOMS = stoi(str);  // convert string to int
	cout << NATOMS << "\n";
	str=readNumber(infile);
	FEC=stoi(str);
	cout << FEC << "\n";
	str=readNumber(infile);
	::sscanf(str.c_str(), "%lA", &FXX); // convert hexa-string to double
	cout << FXX << "\n";
	str=readNumber(infile);
	::sscanf(str.c_str(), "%lA", &FYY);
	cout << FYY << "\n";
	str=readNumber(infile);
	::sscanf(str.c_str(), "%lA", &FZZ);
	cout << FZZ << "\n";
	str=readNumber(infile);
	::sscanf(str.c_str(), "%lA", &FXY);
	cout << FXY << "\n";	
	str=readNumber(infile);
	::sscanf(str.c_str(), "%lA", &FXZ); 
	cout << FXZ << "\n";
	str=readNumber(infile);
	::sscanf(str.c_str(), "%lA", &FYX);
	cout << FYX << "\n";
	str=readNumber(infile);
	::sscanf(str.c_str(), "%lA", &FYZ);
	cout << FYZ << "\n";
	str=readNumber(infile);
	::sscanf(str.c_str(), "%lA", &FZX);
	cout << FZX << "\n";
	str=readNumber(infile);
	::sscanf(str.c_str(), "%lA", &FZY);
	cout << FZY << "\n";	
	str=readNumber(infile);
	//Lp=stod(str);
	::sscanf(str.c_str(), "%lA", &Lp);
	cout << Lp << "\n";	
	str=readNumber(infile);
	::sscanf(str.c_str(), "%lA", &BENDING_PF);
	cout << BENDING_PF << "\n";	
	str=readNumber(infile);
	::sscanf(str.c_str(), "%lA", &StretchEps);
	cout << StretchEps << "\n";	
	str=readNumber(infile);
	rp=stoi(str);
	cout << rp << "\n";
	str=readNumber(infile);
	rpp=stoi(str);
	cout << rpp << "\n";
	
	str=readNumber(infile);
	cout << str << "\n";
	vector<string> strs;
  	boost::split(strs,str,boost::is_any_of(" "));
  	const char * c;
  	for (size_t i = 0; i < strs.size(); i++){
  		c = strs[i].c_str();
  		rIA[i]=strtoul(c,NULL,10);
  		//cout << rIA[i] << endl;
  	}
    //

	//size_t spacePos=str.find(" ");
	//strtoul

	//getline(infile,str); // not reading rIA yet!!!


	posX.resize(NATOMS,0); posY.resize(NATOMS,0); posZ.resize(NATOMS,0);

	string posX_hex,posY_hex,posZ_hex;

	NN.resize(NATOMS,4);
	GB.resize(NATOMS,4);
	RR.resize(NATOMS,4);
	WBX.resize(NATOMS,4);WBY.resize(NATOMS,4);WBZ.resize(NATOMS,4);


	string RR1_hex, RR2_hex, RR3_hex, RR4_hex;

	for(int i=0; i<NATOMS; i++){
		infile >> posX_hex >> posY_hex >> posZ_hex;
		//cout << posX_hex << "\t" <<  posY_hex << "\t" << posZ_hex << "\n";
		
		//posX(i)=stod(posX_hex);
		//posY(i)=stod(posY_hex);
		//posZ(i)=stod(posZ_hex);
		
		::sscanf(posX_hex.c_str(), "%lA", &posX(i));
		::sscanf(posY_hex.c_str(), "%lA", &posY(i));
		::sscanf(posZ_hex.c_str(), "%lA", &posZ(i));
		
		//cout << posX(i) << "\t" <<  posY(i) << "\t" << posZ(i) << "\n";
	}
	for(int i=0; i<NATOMS; i++){
		infile >> NN(i,0) >> NN(i,1) >> NN(i,2) >> NN(i,3);

		//cout << NN(i,0) << "\t" << NN(i,1) << "\t" << NN(i,2) << "\t" << NN(i,3) << "\n";
	}

	for(int i=0; i<NATOMS; i++){
		infile >>  GB(i,0) >> GB(i,1) >> GB(i,2) >> GB(i,3);

		//cout << GB(i,0) << "\t" << GB(i,1) << "\t" << GB(i,2) << "\t" << GB(i,3) << "\n";
	}

	for(int i=0; i<NATOMS; i++){
		infile >> RR1_hex >> RR2_hex >> RR3_hex >> RR4_hex;

		//cout << RR1_hex << "\t" << RR2_hex<< "\t" << RR3_hex << "\t" << RR4_hex << "\n";
		//RR(i,0)=stod(RR1_hex);
		//RR(i,1)=stod(RR2_hex);
		//RR(i,2)=stod(RR3_hex);
		//RR(i,3)=stod(RR4_hex);

		::sscanf(RR1_hex.c_str(), "%lA", &RR(i,0));
		::sscanf(RR2_hex.c_str(), "%lA", &RR(i,1));
		::sscanf(RR3_hex.c_str(), "%lA", &RR(i,2));
		::sscanf(RR4_hex.c_str(), "%lA", &RR(i,3));

		//cout << RR(i,0) << "\t" << RR(i,1) << "\t" << RR(i,2) << "\t" << RR(i,3) << "\n";
	}

	for(int i=0; i<NATOMS; i++){
		infile >> WBX(i,0) >> WBX(i,1) >> WBX(i,2) >> WBX(i,3);
		//cout << WBX(i,0) << "\t" << WBX(i,1) << "\t" << WBX(i,2) << "\t" << WBX(i,3) << "\n";
	}
	
	for(int i=0; i<NATOMS; i++){
		infile >> WBY(i,0) >> WBY(i,1) >> WBY(i,2) >> WBY(i,3);
		//cout << WBY(i,0) << "\t" << WBY(i,1) << "\t" << WBY(i,2) << "\t" << WBY(i,3) << "\n";
	}

	for(int i=0; i<NATOMS; i++){
		infile >> WBZ(i,0) >> WBZ(i,1) >> WBZ(i,2) >> WBZ(i,3);
		//cout << WBZ(i,0) << "\t" << WBZ(i,1) << "\t" << WBZ(i,2) << "\t" << WBZ(i,3) << "\n";
	}
}

std::vector<ublas::matrix<int>> Simulation::createConTabSO(int i){
	std::vector<ublas::matrix<int>> result;
	ublas::matrix<int> conTabSelf(1,2), conTabOther(1,2);
	int numbCon = 0;
	int line = 0;
	int neighLine = 0;
	for(int j=0; j<4; j++)
		{
			if(GB(i,j)==0)
			{
				numbCon++;
			}	
		}
		conTabSelf.resize(numbCon,2);
		conTabOther.resize(numbCon,2);
		for(int j=0;j<4;j++)
		{
			if(GB(i,j)==0)
			{			
				conTabSelf(line,0)=i;
				conTabSelf(line,1)=j;
				neighLine=NN(i,j);
				conTabOther(line,0)=neighLine;
				for(int k=0;k<4;k++)
				{
					if(NN(neighLine,k)==i)
					{
						conTabOther(line,1)=k;
					}
				}
				line++;
			}
		}
		result.push_back(conTabSelf);
		result.push_back(conTabOther);
		return result;
}

void Simulation::createCrosslinks(){
	ublas::vector<double> r(3,0);
	std::vector<ublas::matrix<int>> conTabs;
	ublas::matrix<int> conTabSelf, conTabOther;

	for(int i=0;i<NATOMS;i++){
		r[0]=posX[i]; r[1]=posY[i]; r[2]=posZ[i];

		conTabs=createConTabSO(i);

		conTabSelf = conTabs.at(0);
		conTabOther = conTabs.at(1);

		Crosslink cl = Crosslink(i,r,conTabSelf,conTabOther);
		particles_.push_back(cl);
		//cout << cl.getPosition() << "\n";
		//cout << cl.getForce() << "\n";
	}

	/*for(int i=0; i<4;i++){

		cout << particles_.at(i).getConTabSelf() <<  "\t" << particles_.at(i).getConTabOther() << "\n";
	}*/
}

void Simulation::createPairList(){
	int clId1, clId2;
	int nn, gb;
	ublas::vector<int> wb(3,0);
	int filId1,filId2;
	for(int i=0;i<NATOMS;i++)
	{
		for(int j=0;j<4;j++)
		{
			clId1=i;
			nn=NN(i,j);
			gb=GB(i,j);
			if(nn>clId1 && gb==0)
			{
				clId2=nn;
				if(j<2){
					filId1=1;
				}
				else{
					filId1=2;
				}
				// check which filament is cl1 for cl2 on
				for(int k=0;k<4;k++)
				{
					if (NN(nn,k)==i && GB(nn,k)==0) //find cl1 in cl2-line of NN table
					{
						if(k<2){filId2=1;}
						
						else{filId2=2;}
					}
				}

				CrosslinkPair2 clpair=CrosslinkPair2(particles_.at(clId1), particles_.at(clId2),RR(i,j),WBX(i,j),WBY(i,j),WBZ(i,j), filId1, filId2);
				particlePairs2_.push_back(clpair);		
			}
		}
	}
	//WBX(1,1)=5;
	//RR(1,1)=7;
	/*for(int i=0; i<particlePairs2_.size(); i++){
		cout << particlePairs2_.at(i).getCrosslink(1)->getId() << "\t" << particlePairs2_.at(i).getCrosslink(2)->getId() << "\t"
			 << particlePairs2_.at(i).getContourLength() << "\t" << particlePairs2_.at(i).getWBX() << "\t" << particlePairs2_.at(i).getWBY() <<
			 "\t" << particlePairs2_.at(i).getWBZ() << "\t" << particlePairs2_.at(i).getFilId(1) << "\t" << particlePairs2_.at(i).getFilId(2) <<  
			 "\t" << particlePairs2_.at(i).getContourLength() << "\n";
	}*/
}

void Simulation::createTripleList(){
	int clId1, clId2, clId3;
	ublas::vector<int> wb1(3,0),wb2(3,0);
	int iDummy; // index of dummy crosslink, 0=no dummy
	int filId1,filId2,filId3;

	for(int i=0; i<NATOMS; i++)
	{
		if(GB(i,0)==0 && GB(i,1)==0){
			clId1=NN(i,0);
			clId2=i;
			clId3=NN(i,1);

			iDummy=0;
		}
		// if cl is the end of a polyer, then one GB=1, that crosslink becomes dummy crosslink
		else if(GB(i,0)==1 && GB(i,1)==0){
			clId1=i; // dummy crosslink
			clId2=i;
			clId3=NN(i,1);

			iDummy=1;
		}
		else if(GB(i,0)==0 && GB(i,1)==1){
			clId1=NN(i,0); 
			clId2=i;
			clId3=i;// dummy crosslink

			iDummy=3;
		}

		filId2=1;

		filId1=-1;
		filId3=-1;


		for(int k=0;k<4;k++)
		{	if(iDummy!=1)
			{
				if (NN(clId1,k)==i && GB(clId1,k)==0) 
				{
					if(k<2){filId1=1;}
							
					else{filId1=2;}
				}
			}
			if(iDummy!=3)
			{
				if (NN(clId3,k)==i && GB(clId3,k)==0) 
				{
					if(k<2){filId3=1;}
							
					else{filId3=2;}
				}
			}

		}

		int n1a, n1b, n3a, n3b;

		if (iDummy==0)
		{
			n1a=clId1;
			n3a=clId3;
			for(int k=0;k<4;k++)
			{
				if (NN(clId1,k)==i && GB(clId1,k)==0){ n1b=k;}
				if (NN(clId3,k)==i && GB(clId3,k)==0){ n3b=k;}
			}
		}

		if (iDummy==1)
		{
			n1a=i;
			n1b=0;
			n3a=clId3;
			for(int k=0;k<4;k++)
			{
				if (NN(clId3,k)==i && GB(clId3,k)==0){ n3b=k;}
			}
		}

		if (iDummy==3)
		{
			n3a=i;
			n3b=1;
			n1a=clId1;
			for(int k=0;k<4;k++)
			{
				if (NN(clId1,k)==i && GB(clId1,k)==0){ n1b=k;}
			}
		}


		//CrosslinkTriple2 cltriple=CrosslinkTriple2(particles_.at(clId1), particles_.at(clId2), particles_.at(clId3));
		CrosslinkTriple2 cltriple=CrosslinkTriple2(particles_.at(clId1), particles_.at(clId2), particles_.at(clId3), RR(i,0), RR(i,1), RR(n1a,n1b), RR(n3a,n3b),
				WBX(i,0), WBY(i,0), WBZ(i,0), WBX(i,1), WBY(i,1), WBZ(i,1),iDummy,filId1,filId2,filId3);

		particleTriples2_.push_back(cltriple);	
		

		if(GB(i,2)==0 && GB(i,3)==0){
			clId1=NN(i,2);
			clId2=i;
			clId3=NN(i,3);

			iDummy=0;
		}
		// if cl is the end of a polyer, then one GB=1, that crosslink becomes dummy crosslink
		else if(GB(i,2)==1 && GB(i,3)==0){
			clId1=i; // dummy crosslink
			clId2=i;
			clId3=NN(i,3);

			iDummy=1;
		}
		else if(GB(i,2)==0 && GB(i,3)==1){
			clId1=NN(i,2); 
			clId2=i;
			clId3=i;// dummy crosslink

			iDummy=3;
		}
			//CrosslinkTriple2 cltriple=CrosslinkTriple2(particles_.at(clId1), particles_.at(clId2), particles_.at(clId3));

		filId2=2;

		filId1=-1;
		filId3=-1;

		for(int k=0;k<4;k++)
		{	if(iDummy!=1)
			{
				if (NN(clId1,k)==i && GB(clId1,k)==0) 
				{
					if(k<2){filId1=1;}
							
					else{filId1=2;}
				}
			}
			if(iDummy!=3)
			{
				if (NN(clId3,k)==i && GB(clId3,k)==0) 
				{
					if(k<2){filId3=1;}
							
					else{filId3=2;}
				}
			}

		}

		if (iDummy==0)
		{
			n1a=clId1;
			n3a=clId3;
			for(int k=0;k<4;k++)
			{
				if (NN(clId1,k)==i && GB(clId1,k)==0){ n1b=k;}
				if (NN(clId3,k)==i && GB(clId3,k)==0){ n3b=k;}
			}
		}

		if (iDummy==1)
		{
			n1a=i;
			n1b=2;
			n3a=clId3;
			for(int k=0;k<4;k++)
			{
				if (NN(clId3,k)==i && GB(clId3,k)==0){ n3b=k;}
			}
		}

		if (iDummy==3)
		{
			n3a=i;
			n3b=3;
			n1a=clId1;
			for(int k=0;k<4;k++)
			{
				if (NN(clId1,k)==i && GB(clId1,k)==0){ n1b=k;}
			}
		}

		CrosslinkTriple2 cltriple2=CrosslinkTriple2(particles_.at(clId1), particles_.at(clId2), particles_.at(clId3), RR(i,2), RR(i,3), RR(n1a,n1b), RR(n3a,n3b),
			WBX(i,2), WBY(i,2), WBZ(i,2), WBX(i,3), WBY(i,3), WBZ(i,3),iDummy,filId1,filId2,filId3);

		particleTriples2_.push_back(cltriple2);	
	}	
	
	
	
	/*for(int i=0; i<particleTriples2_.size(); i++){
		cout << particleTriples2_.at(i).getCrosslink(1)->getId() << "\t" << particleTriples2_.at(i).getCrosslink(2)->getId() 
		<< "\t" << particleTriples2_.at(i).getCrosslink(3)->getId()<< "\t" << particleTriples2_.at(i).getContourLength(12)
		<< "\t" << particleTriples2_.at(i).getContourLength(23) << "\t" << particleTriples2_.at(i).getWBX(12) << "\t" <<
		particleTriples2_.at(i).getWBY(12) <<  "\t" << particleTriples2_.at(i).getWBZ(12) << "\t" << particleTriples2_.at(i).getWBX(23) << "\t" <<
		particleTriples2_.at(i).getWBY(23) <<  "\t" << particleTriples2_.at(i).getWBZ(23) << "\t" << particleTriples2_.at(i).getDummy() << "\t" <<
		particleTriples2_.at(i).getFilId(1)<<  "\t" << particleTriples2_.at(i).getFilId(2)<< "\t" << particleTriples2_.at(i).getFilId(3) << 
		"\t" << particleTriples2_.at(i).getContourLength(12) << "\t" <<particleTriples2_.at(i).getContourLength(23) <<"\n";
	}*/
}


void Simulation::shear_(){
	FXY=gamma0_*sin(omega_*t_);
}

void Simulation::correctRR_(int opt)
{
	double randNo;
	for(int i=0;i<NATOMS;i++)
	{
		for(int j=0;j<4;j++)
		{
			if(GB(i,j)==1){
				if(opt==1)
				{
					randNo=0.001 + 0.09*((double) rand() / (RAND_MAX));
				}
				else{randNo=0;}
				RR(i,j)=randNo;
			}
		}
	}
}


void Simulation::compareForces_()
{
	ublas::vector<double> Ftot(3,0);
	ublas::vector<double> r_i, r_j;
	double Ftot_ij=0;
	double F=0;
	for(int i=0;i<NATOMS;i++)
	{
		//Ftot += ublas::norm_2(particles_.at(i).getForce());
		Ftot = particles_.at(i).getForce();
		r_i = particles_.at(i).getPosition();
		for(int j=0;j<4;j++)
		{
			if(GB(i,j)==0)
			{
				r_j=particles_.at(NN(i,j)).getPosition();
				Ftot_ij=abs(ublas::inner_prod(Ftot,r_i-r_j));
				F+=Ftot_ij;
			}
		}
	}
	F=F/NATOMS;

	double Fpull=0;
	int nTriples=0;
	for (std::vector<CrosslinkTriple2> ::iterator triplei = particleTriples2_.begin(); triplei != particleTriples2_.end(); triplei++)
	{
		Fpull+=abs(triplei->getfPull());
		nTriples++;
	}
	Fpull=Fpull/nTriples;

	cout << F << "\t" << Fpull << "\n";
}