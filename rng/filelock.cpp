#include "filelock.hpp"


/* 28 APRIL 2008 :: FILE LOCKING ROUTINES

lock a file by hard-linking to it : this should also work 
over NFS but, of course, here one never knows */

/* how long to wait in seconds before giving up */
#define TIMEOUT	120.00

/* some unlikely string */
#define LOCK	"lock"

/* upon obtaining a successfull lock, return 1 if 
   fname already existed; otherwise return 0 */
extern char outputDirectory[1024];
int RUNID;

int lockfile(char *fname)
{ int ok,dt,status; 
struct stat info;
char lname[1024];
FILE *fp;

sprintf(lname,"%s.%s",fname,LOCK); dt=0; status=1;

while(1) {

	ok=link(fname,lname);

	if(ok==0) break;

	else if(errno==ENOENT) {

		if(status==0) printerror("E: filelock: files are created only once\n");

		/* file fname does not exist yet so we safely create it */
		fp=fopen(fname,"a");
		if(fp==NULL) printerror("E: filelock: failed to create file\n");
		fclose(fp); status=0;

	} else if(errno==EEXIST) {

		/* file is locked by another process */
		printf("W: %s is locked by another process\n",fname); sleep(1); dt++;
		if(dt>TIMEOUT) printerror("E: lockfile: timeout\n");

	} else printerror("E: lockfile: errno=%i %s\n",errno,strerror(errno));

}

/* some final additional check, just to be safe */
if(stat(lname,&info)!=0) printerror("E: lockfile: magic_mystery_error\n");

return status;
}

/* release the lock */

void unlockfile(char *fname)
{ char lname[1024];

sprintf(lname,"%s.%s",fname,LOCK);
if(unlink(lname)!=0) printerror("E: unlockfile: failed to unlock %s\n",fname);
}

/* set runid safely using file locking and set the output file */


#undef LOCK
#undef TIMEOUT

/* 18 JULY 2006 :: Write data safely, that is, to a tempfile first. Only
after the tempfile has been written to, is it renamed to the requested
file name */

static char FNAME[1024],BNAME[1024],DNAME[1024],DNAME2[1024];

FILE* savewrite(char *fname)
{
	int fd; FILE *fp=NULL;

	sprintf(FNAME,"%s",fname);
	sprintf(BNAME,"%s.%i.bak.XXXXXX",fname,RUNID);

	int size = 1024;
	if (strcmp(outputDirectory, "") != 0){
		getcwd(DNAME, size);
		printm("pwd: %s\n", DNAME);
		if (chdir(outputDirectory) != 0) printerror("Unable to chdir to %s\n", outputDirectory);
	}
	printm("BNAME = %s\n", BNAME);
	char bn[1024];
	strcpy(bn,BNAME);
    fd=mkstemp(bn);
	if(fd==-1) {printerror("ERROR: NOMKSTEMP\n");}
	else close(fd);
	strcpy(BNAME,bn);

	fp=fopen(BNAME,"w");
	if (fp==NULL) printerror("Unable to create file %s\n", FNAME);
	getcwd(DNAME2, size);
	printm("pwd: %s\n", DNAME2);
	return fp;
}

FILE* appendwrite(char *fname)
{
	int fd; FILE *fp=NULL;

	sprintf(FNAME,"%s",fname);
	sprintf(BNAME,"%s.%i.bak.XXXXXX",fname,RUNID);

	int size = 1024;
	if (strcmp(outputDirectory, "") != 0){
		getcwd(DNAME, size);
		printm("pwd: %s\n", DNAME);
		if (chdir(outputDirectory) != 0) printerror("Unable to chdir to %s\n", outputDirectory);
	}
	printm("BNAME = %s\n", BNAME);
	char bn[1024];
	strcpy(bn,BNAME);
    fd=mkstemp(bn);
	if(fd==-1) {printerror("ERROR: NOMKSTEMP\n");}
	else close(fd);
	strcpy(BNAME,bn);

	fp=fopen(BNAME,"a");
	if (fp==NULL) printerror("Unable to create file %s\n", FNAME);
	getcwd(DNAME2, size);
	printm("pwd: %s\n", DNAME2);
	return fp;
}

void saveclose(FILE *fp)
{
	printm("pwd: %s\n", DNAME2);
	fflush(fp); fclose(fp);
	rename(BNAME,FNAME);
	if (strcmp(outputDirectory, "") != 0){
		if (chdir(DNAME) != 0) printerror("Unable to chdir to %s\n", DNAME);
		int size = 1024;
		getcwd(DNAME, size);
		printm("pwd: %s\n", DNAME);
	}

}
