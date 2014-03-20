// main.cpp
// Il primo esempio in C++
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <stdlib.h>
#include "utility.h"
#include <math.h>
#include "bio.h"
using namespace std;

int main(int argc, char *argv[])
{
	ifstream file ("../bioinfo/illumina.sam");
	if(file.is_open())
	{
		cout << "*** Sam file opened! ***" << endl;
		string line;
		getline(file, line);
		cout << "HEADER: " << line << endl;

		// number of reads for size of array
		strtok(stc(line), ":");
		strtok(NULL, ":");
		int ln = atoi(strtok(NULL, ":"));
		
		// Statistic informations
		long int sumSizes = 0;
		int totalReadsLoads = 0;
		int mappedReads = 0;
		int notMappedReads = 0;
		float st_dev;
		float mean;
		
		// Load total lines into vector sam
		vector<SamLine> sam(ln);
		int i = 0;
		int count = 0;
		while(getline(file, line))
		{	
			sam.at(i).readInfo(line);

			// Conditions
			if(strcmp(sam.at(i).cigar, "*"))
			{ // Mapped read
				mappedReads++;
				
				if(i == 0 || abs(sam.at(i).size) < 2*mean)
				{
					count++;
					sumSizes += abs(sam.at(i).size);
					mean = (float) sumSizes / count;
				}
			}
			i++;
		}

		// Total reads load from file
		totalReadsLoads = i;

		//print info
		cout << "Total reads = " << totalReadsLoads << endl;
		cout << "Mapped reads = " << mappedReads << endl;
		cout << "Mean size = " << mean << endl;

		// Calculate st_dev

		int * bases = new int[ln]();
		long int sum_st_dev = 0;
		for(i = 0; i < totalReadsLoads; i++)
		{
			if(strcmp(sam.at(i).cigar, "*") &&
				abs(sam.at(i).size) < 2*mean )
				sum_st_dev += pow(abs(sam.at(i).size) - mean  ,2);
			
			

		}
		st_dev = sqrt(sum_st_dev / mappedReads);
		cout << "Standard deviation = " << st_dev << endl;

		file.close();
	}
	return 0;
}