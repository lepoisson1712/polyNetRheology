#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <iostream>
#include "filelock.hpp"
#include "rngmit.hpp"	// header file for the random number generator rngmit
#include "defs.hpp"

using namespace std;

// the following are external variables declared in rngmit source files
extern unsigned long int rng_ia[55];
extern int rng_p,rng_pp;
extern int rng_mod[55];

char inputFile[1024];
char outputDirectory[1024] = "";
bool NOT_READ;
bool rngSet = false;
long SEED;
FILE *OUT;

int NX=1;
double FXX=1, FXY=0, FXZ=0, FYX=0, FYY=1, FYZ=0, FZX=0, FZY=0, FZZ=1;
double Lp=1;

void spit_msg(void);
void init_rngmit_state_from_file(char inputFile[]);
void load_exact_dbl_numbers(char *filename);
void save_rngmit_state_and_exact_dbl_numbers(char *filename);

int main(int argc,char **argv)
{
/*
	Usage:
		forandy -rng <network file>
		OR
		forandy -s <network file>
	Options:
        forandy -rng <network file>
			(<network file> should have random number generator state
			and hexadecimal representation of floating-point numbers),
	    forandy -s <network file>
			To save some arbitrary data to <network file> in the
			required format.
*/

	// The following tests and reads the program arguments
	if (argc < 2) spit_msg();
	else {
		for (int i = 1; i < argc; i++) {
			if (string(argv[i]) == "-rng") {
				rngSet = true;
				strcpy(inputFile, argv[i + 1]); i++;
			}
			else if (string(argv[i]) == "-s") {
				strcpy(inputFile, argv[i + 1]); i++;
			}
			else spit_msg();
		}
	}

	if (rngSet){
		// The following loads the state of the random number generator (rngmit)
		// from the file 'inputFile' into the program and initializes rngmit.
		// If the file does not contain the state
		// it initializes rngmit with the global variableSEED
		init_rngmit_state_from_file(inputFile);

		// The following loads the values of double-precision global variables
		// from the file 'inputFile' into the program.
		load_exact_dbl_numbers(inputFile);
	}
	else rngseed(SEED);

	// to generate a few random numbers between 0 and 1:
	printf("Some random numbers from rngmit before saving state:\n");
	printf("%lf\n", rngmit);
	printf("%lf\n", rngmit);
	printf("%lf\n", rngmit);

	// The following saves the current state of the random number generator (rngmit)
	// and some double-precision numbers into the file 'inputFile'
	save_rngmit_state_and_exact_dbl_numbers(inputFile);

	printf("Some random numbers from rngmit after saving state:\n");
	printf("%lf\n", rngmit);
	printf("%lf\n", rngmit);
	printf("%lf\n", rngmit);
	return 0;
}

void spit_msg(void)
{
	std::cout << "Not enough or invalid arguments, please try again.\n";
	std::cout << "Usage:\n";
	std::cout << "   forandy -rng <network file> \n"
			"(<network file> should have random number generator state \n"
			"and hexadecimal representation of floating-point numbers),\n";
			std::cout << "OR   forandy -s <network file> \n"
			"to save some arbitrary data to <network file> in the required format.\n";
	std::cout.flush();
	std::cin.get();
	exit(0);
}

void init_rngmit_state_from_file(char inputFile[])
{
	int rp, rpp;
	unsigned long int rIA[55];
	SETNBVAL(int,rp,inputFile);
	SETNBVAL(int,rpp,inputFile);
	SETARRVAL(rIA,55,inputFile);
	if (NOT_READ) rngseed(SEED);
	else rnginit(rIA, rpp, rp);
}


void load_exact_dbl_numbers(char *filename)
{
	SETNBVAL(int,NX,inputFile);
	SETVAL(FXX,inputFile);
	SETVAL(FYY,inputFile);
	SETVAL(FZZ,inputFile);
	SETVAL(FXY,inputFile);
	SETVAL(FXZ,inputFile);
	SETVAL(FYZ,inputFile);
	SETVAL(FYX,inputFile);
	SETVAL(FZX,inputFile);
	SETVAL(FZY,inputFile);
	SETVAL(Lp,inputFile);
}


void save_rngmit_state_and_exact_dbl_numbers(char *filename)
	{
	FILE *fp; char fname[1024];
	sprintf(fname,"%s",filename); fp=savewrite(fname);

	fprintf(fp,"#NATOMS=%i\n",NX);

	// The following saves some double precision numbers in hexadecimal format
	fprintf(fp,"#FXX=%a\n",FXX);
	fprintf(fp,"#FYY=%a\n",FYY);
	fprintf(fp,"#FZZ=%a\n",FZZ);
	fprintf(fp,"#FXY=%a\n",FXY);
	fprintf(fp,"#FXZ=%a\n",FXZ);
	fprintf(fp,"#FYZ=%a\n",FYZ);
	fprintf(fp,"#FYX=%a\n",FYX);
	fprintf(fp,"#FZX=%a\n",FZX);
	fprintf(fp,"#FZY=%a\n",FZY);
	fprintf(fp,"#Lp=%a\n",Lp);

	// The following saves the rngmit state
	fprintf(fp,"#rp=%d\n",rng_p);
	fprintf(fp,"#rpp=%d\n",rng_pp);
	fprintf(fp,"#rIA=%lu",rng_ia[0]);
	for(int j=1;j<55;j++) fprintf(fp," %lu",rng_ia[j]);
	fprintf(fp,"\n");

	// The following closes the file
	fflush(fp);
	saveclose(fp);

}
