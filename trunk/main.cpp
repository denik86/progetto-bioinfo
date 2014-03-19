// main.cpp
// Il primo esempio in C++
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <stdlib.h>
#include "utility.h"
using namespace std;

class SamLine{
	public:
		// Info reads
		//char *qname;
		int flag;
		//char *chromo;
		int pos;
		int mapq;
		int cigar;
		//char *mrnm;
		int mpos;
		int size;
		//char *seq;
		//char *qual;

		int vFlag[7];

		void readInfo(string l);
		void flagToBinaryVector();
};


void SamLine::readInfo(string l)
{
	char * line = stc(l);
	/*qname = */strtok(line, "\t");
	flag = atoi(strtok(NULL, "\t"));
	/*chromo = */strtok(NULL, "\t");
	pos = atoi(strtok(NULL, "\t"));
	mapq = atoi(strtok(NULL, "\t"));
	cigar = atoi(strtok(NULL, "\t"));
	/*mrnm = */strtok(NULL, "\t");
	mpos = atoi(strtok(NULL, "\t"));
	size = atoi(strtok(NULL, "\t"));
	//seq = strtok(NULL, "\t");
	//qual = strtok(NULL, "\t");
}

void SamLine::flagToBinaryVector()
{
	for (int i = 0; i < 7; ++i) 
	{  // assuming a 7 bit int
    	vFlag[i] = flag & (1 << i) ? 1 : 0;
    }
}



int main(int argc, char *argv[])
{
	
	ifstream file ("../bioinfo/illumina.sam");
	if(file.is_open())
	{
		cout << "*** Sam file opened! ***" << endl;
		string line;
		getline(file, line);
		cout << "HEADER: " << line << endl;

		// number of reads
		int ln;
		strtok(stc(line), ":");
		strtok(NULL, ":");
		ln = atoi(strtok(NULL, ":"));
		
		// Load total lines into array sam
		
		vector<SamLine> sam(ln);
		getline(file, line);
		sam.at(0).readInfo(line);
		cout << sam.at(0).flag << " KByte" << endl;


		
		int i = 0;
		while(getline(file, line))
		{	/*
			sam[i].readInfo(line);
			i++;*/
		}

		/*
		// Statistic informations
		long int totalReads = 0;
		long int pairedReads = 0;
		long int notMappedReads = 0;
		
		while(getline(file,line))
		{
			sam[i-1].readInfo(line);

		
			//SamLine sl;
			//sl.readInfo(line);
			
			if(sl.size > 0)
				sumSize = sumSize + sl.size;
			else if (sl.size < 0)
				sumSize = sumSize - sl.size;
				
			else
			{
				continue;
			} 
			
			
			totalLines++;
			
		}*/
		file.close();

		//int media = sumSize / totalLines;
		//cout << sumSize << "/" << totalLines << " = " << media << endl;
	}
	return 0;
}