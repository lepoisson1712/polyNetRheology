#include "defs.hpp"

extern bool NOT_READ;

int contain(char *s1,char *s2)
{
	int i,j,ok,l1,l2;

	l1=strlen(s1); l2=strlen(s2);

	/* detect blank lines */
	for(i=0;i<l1;i++) if(s1[i]!='\n' && isspace(s1[i])==0) goto next;
	return -1; next:

	for(i=0;i<l1-l2;i++) {
		for(j=0,ok=1;j<l2;j++) if( s1[i+j]!=s2[j] ) ok=0;
		if(ok) return i+l2;
	}
	return 0;
}




double readvalue(const char *str,const char *file)
{
	int n; FILE *fp;
	double tmp;
	char head[1024];
	char line[55*1024];

	sprintf(head,"#%s=",str); fp=fopen(file,"r");

	if(fp==NULL) {
		NOT_READ = true;
		printerror("ERROR: readvalue: file %s not found while trying to read %s\n",file,str);
	}

	while( fgets(line,1024,fp)!=NULL ) {
		n=contain(line,head); if(n>0) {
			// sscanf(&line[n],"%lf",&tmp);
			tmp = strtod(&line[n], NULL);
			fclose(fp);
			printm("%s: %s=%f\n",file,str,tmp);
			// if(OUT!=NULL) fprintf(OUT,"#%s=%f\n",str,tmp);
			NOT_READ = false;
			return tmp;
		}
	}
	NOT_READ = true;
	fclose(fp);
	printm("ERROR: could not read %s from %s\n",str,file);

	return -1.0;
}


void readLintarray(const char *str,const int sz, const char *file, unsigned long int arr[])
{
	int i,n; FILE *fp;
	char head[1024];
	char line[55*1024];

	sprintf(head,"#%s=",str); fp=fopen(file,"r");

	if(fp==NULL) {
		NOT_READ = true;
		printerror("ERROR: readvalue: file %s not found while trying to read %s\n",file,str);
	}

	while( fgets(line,55*1024,fp)!=NULL ) {
		n=contain(line,head); if(n>0) {
			char* ptr = strtok(&line[n], " ");
			i=0;
			while(ptr != NULL)
			{
				arr[i] = strtoul(ptr,NULL,0);i++;
				printf("%s[%d] = %lu \n", str, i-1, strtoul(ptr,NULL,0));
				// if(OUT!=NULL) fprintf(OUT,"%s[%d] = %lu \n", str, i, strtoul(ptr,NULL,0));
				ptr = strtok(NULL, " ");
			}
			fclose(fp);
			NOT_READ = false;
			return;
		}
	}
	NOT_READ = true;
	fclose(fp);
	printm("ERROR: could not read %s from %s\n",str,file);

}

void printerror1(const char *msg, ...)
{
	char dd[1024];
	struct tm *ptr;
	time_t tm;
	va_list argp;

	/* generate date label dd */

	tm=time(NULL); ptr=localtime(&tm);
	strftime(dd,100,"%b-%d-%Y-%H.%M",ptr);

	/* print error to screen... */

	va_start(argp,msg);
	vfprintf(stderr,msg,argp);
	va_end(argp);

	/* ...and to logfile... */

	exit(1);
}


