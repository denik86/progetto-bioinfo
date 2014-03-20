#include <stdlib.h>
#include "utility.h"
#include "bio.h"

void SamLine::readInfo(string l)
{
	char * line = stc(l);
	/*qname = */strtok(line, "\t");
	flag = atoi(strtok(NULL, "\t"));
	flagToBinaryVector();
	/*chromo = */strtok(NULL, "\t");
	pos = atoi(strtok(NULL, "\t"));
	mapq = atoi(strtok(NULL, "\t"));
	cigar = strtok(NULL, "\t");
	/*mrnm = */strtok(NULL, "\t");
	mpos = atoi(strtok(NULL, "\t"));
	size = atoi(strtok(NULL, "\t"));
	//seq = strtok(NULL, "\t");
	//qual = strtok(NULL, "\t");

}

void SamLine::flagToBinaryVector()
{
	for (int i = 0; i < 12; ++i) 
	{  // assuming a 7 bit int
    	vFlag[i] = flag & (1 << i) ? 1 : 0;
    }
}
