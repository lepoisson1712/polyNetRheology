#ifndef DEFS2_H_
#define DEFS2_H_


/* some name to describe the program */
#define IDSTRING        "WWW-X-link.2010"

/* where errors/warnings are written */
#define ERRORFILE	"JOBINFO"

/*************************
* do not edit below here *
*************************/

/* structure defs go here */

/* end-of-structure defs */

#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/vfs.h>
#include "rngmit.hpp"
//#include "head.hpp"
//#include "tools.hpp"

// FUNCTION HEADERS
double readvalue(const char *str,const char *file);
void readLintarray(const char *str,const int sz, const char *file, unsigned long int arr[]);
void printerror1(const char *msg, ...);
void init_rngmit_state_from_file(char inputFile[]);
void save_rngmit_state_and_exact_dbl_numbers(char *filename);


// CUSTOM MACROS
#define SETVAL(var,file) var=readvalue(#var,file)
#define SETARRVAL(var,sz,file) readLintarray(#var,sz,file,var)
#define SETNBVAL(type,var,file) var=(type)readvalue(#var,file)
#define printm(...) { printf("# %s: ",__FUNCTION__); printf(__VA_ARGS__); }
#define printerror(...) { printf("# %s: ",__FUNCTION__); printerror1(__VA_ARGS__); }

#define IFSETVAL(var,file) var=ifreadvalue(#var,var,file)

/* print meaningful messages to screen */


#endif // DEFS2_H_
