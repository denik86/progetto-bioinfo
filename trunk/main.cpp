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
int ln; // lenght of reference genome
int nReads = 0; // number of reads loaded
int l;
int mappedReads = 0;
int notMappedReads = 0;
float mean;
float sd;

//Data Structures
vector<SamLine> sam; // Vector of lines from Sam file

void loadData(const char* path, int valueOutlier, int lenghtRead)
{
	ifstream file (path);
	if(file.is_open())
	{
		cout << "********************************************" << endl;
		cout << "*** File '" << path << "' opened! ***" << endl;
		string line;
		bool header = true;
		while(header && getline(file, line))
		{
			if(line.substr(0,1).compare("@") == 0)
			{
				if(line.substr(1,2).compare("SQ") == 0)
				{// lenght of reference genome
					strtok(stc(line), ":");
					strtok(NULL, ":");
					ln = atoi(strtok(NULL, ":"));
				}
				cout << "HEADER: " << line << endl;
			}
			else
				header = false;
		}

		// Load total lines into vector sam
		int count = 0;
		long int sumSizes = 0;
		while(getline(file, line)) {	
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
		nReads = sam.size();
		l = lenghtRead;

		//print info
		cout << "Total reads = " << nReads << endl;
		cout << "Mapped reads = " << mappedReads << endl;
		cout << "Not mapped reads = " << notMappedReads << endl;
		cout << "Mean size = " << mean << endl;

		file.close();
	}
	else {
		cout << "Error file" << endl;
	}
}

void st_dev()
{
	long int sum_st_dev = 0;
	int count = 0;
	for(int i = 0; i < nReads; i++)
	{
		if(strcmp(sam.at(i).cigar, "*") &&
			abs(sam.at(i).size) < 1.8*mean )
			 {
			 	sum_st_dev += pow(abs(sam.at(i).size) - mean  ,2);
			 	count++;
			 }
	}
	sd = sqrt(sum_st_dev / count);
	cout << "Standard deviation = " << sd << endl;
}


void physicalCov(const char * path)
{
	// quantitative info of each base
	vector<int> bases(ln, 0);

	// iterate all the reads
	for(int i = 0; i < sam.size(); i++)
	{
		if(sam[i].valid())
		{
			// for each read iterate all position base
			int pos = sam.at(i).pos;
			for(int j = pos; j < pos+l; j++)
				bases[j]++;
		}
	}

	wig(bases, ln, path);
	bases.clear();
}

/* This function IMPLIED that the reads in the sam file
 * are sorted by name */
void multiCoverage(const char * path)
{
	cout << "Analisys of multiple reads coverage ... " << flush;
	vector<int> bases(ln, 0);
	for(int i = 0; i < sam.size(); i++)
	{
		if(i > 0 && strcmp(sam[i].name, sam[i-1].name) == 0 &&
			sam[i].num == sam[i-1].num)
		{ // two reads are the same
			int pos1 = sam[i-1].pos;
			int pos2 = sam[i].pos;
			for (int i = 0; i < l; i++)
			{
				bases[i+pos1]++;
				bases[i+pos2]++;
			}
		}
	}

	wig(bases, ln, path);
	bases.clear();
	cout << "DONE." << endl;
}

void orientation(const char * path, int orientFlag)
{
	if(orientFlag == 0 || orientFlag == 1)
	{
		cout << "Analisys of reads orientation " <<
		orientFlag << " ... " << flush;

		// quantitative info of each base
		vector<int> totBase(ln, 0);
		vector<int> orBase(totBase);

		// iterate all the reads
		for(int i = 0; i < sam.size(); i++)
			if(sam[i].valid())
			{  
				int pos = sam[i].pos;
				for(int j = pos; j < pos+l; j++)
				{
					orBase[j] += (sam[i].vFlag[4] == orientFlag);
					totBase[j]++;
				}
			}

		vector<float> percent(ln, 0);
		for(int j = 0; j < ln; j++)
		{
			if(totBase[j] == 0)
				continue;
			percent[j] = ((float) orBase[j] / totBase[j]) * 100;
		}

		wig(percent, ln, path);
		totBase.clear();
		orBase.clear();
		percent.clear();

		cout << "DONE." << endl;
	}
	else
		cout << "ERROR orientFlag"<< endl;

}

void wrongMate(const char *path)
{
	cout << "Analisys of wrong mate pair orientation ... " << flush;
	vector<int> totBase(ln, 0);
	vector<int> orBase(totBase);

	// iterate all the reads
	for(int i = 0; i < sam.size(); i++)
		if(strcmp(sam[i].cigar, "*"))
		{  


		}
}

int main(int argc, char *argv[])
{
	cout << "cls()" << endl;
	loadData("../bioinfo/allname.sam", 2, 50);
	multiCoverage("multiple.wig");
	//st_dev();
	physicalCov("physicalCoverage.wig");
	//orientation("right.wig", 0);
	//orientation("left.wig", 1);

	return 0; 
}