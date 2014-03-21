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

// Statistic informations
int ln; // number of reads writed on file
int totalReadsLoads = 0; // number of reads loaded
int mappedReads = 0;
int notMappedReads = 0;
float st_dev;
float mean;

//Data Structures
vector<SamLine> sam; // Vector of lines from Sam file
int * bases; // quantitative info of each base

void loadData(const char* path, int valueOutlier)
{
	ifstream file (path);
	if(file.is_open())
	{
		cout << "*** Sam file opened! ***" << endl;
		string line;
		getline(file, line);
		cout << "HEADER: " << line << endl;

		// number of reads writed on file Sam
		strtok(stc(line), ":");
		strtok(NULL, ":");
		ln = atoi(strtok(NULL, ":"));
		
		// Load total lines into vector sam
		int count = 0;
		long int sumSizes = 0;
		while(getline(file, line))
		{	
			SamLine sl;
			sl.readInfo(line);
			sam.push_back(sl);

			// Conditions
			if(strcmp(sl.cigar, "*"))
			{ // Mapped read
				mappedReads++;
				if(count == 0 || abs(sl.size) < valueOutlier*mean)
				{
					count++;
					sumSizes += abs(sl.size);
					mean = (float) sumSizes / count;
				}
			}
			else
				notMappedReads++;
			
		}

		// Total reads load from file
		totalReadsLoads = sam.size();

		//print info
		cout << "Total reads = " << totalReadsLoads << endl;
		cout << "Mapped reads = " << mappedReads << endl;
		cout << "Not mapped reads = " << notMappedReads << endl;
		cout << "Mean size = " << mean << endl;

		file.close();
	}
	else
	{
		cout << "Errore file" << endl;
	}
}



int main(int argc, char *argv[])
{
	
		loadData("../bioinfo/illumina.sam", 2);
		
		
		// Calculate st_dev
		long int sum_st_dev = 0;
		int count = 0;
		for(int i = 0; i < totalReadsLoads; i++)
		{
			if(strcmp(sam.at(i).cigar, "*") &&
				abs(sam.at(i).size) < 1.8*mean )
				 {
				 	sum_st_dev += pow(abs(sam.at(i).size) - mean  ,2);
				 	count++;

				 }
		}
		st_dev = sqrt(sum_st_dev / count);
		cout << "Standard deviation = " << st_dev << endl;

	
	return 0; 
}