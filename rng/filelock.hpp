/*
 * filelock.hpp
 *
 *  Created on: Aug 10, 2012
 *      Author: amuasi
 */

#ifndef FILELOCK_HPP_
#define FILELOCK_HPP_
#include "defs.hpp"

int lockfile(char *fname);
void unlockfile(char *fname);
void getrunid(void);
void ungetrunid(void);
FILE *savewrite(char *fname);
FILE* appendwrite(char *fname);
void saveclose(FILE *fp);

#endif /* FILELOCK_HPP_ */
