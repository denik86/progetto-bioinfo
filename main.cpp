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
int mappedReads = 0;
int notMappedReads = 0;
float mean;
float sd;

//Data Structures
vector<SamLine> sam; // Vector of lines from Sam file

void loadData(const char* path, int valueOutlier)
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
			if(sl.valid())
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
		if(sam[i].valid() && sam[i].size > 0 &&
			sam[i].size < 2*mean )
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
			for(int j = 0; j < sam[i].lenght; j++)
				bases[sam[i].pos + j] += sam[i].match[j];
		}
	}
	wig(bases, ln, path);
	bases.clear();
}

void sequenceCov(const char * pathCoverage, const char *pathAvLenght)
{
	cout << "Analisys of sequence coverage and average lenght inserts ... " << flush;
	// quantitative info of each base
	vector<int> count(ln, 0);
	vector<int> lenght(ln, 0);
	
	// iterate all the reads
	for(int i = 0; i < sam.size(); i++)
	{
		if(sam[i].valid() && sam[i].size > 0)
		{
			// for each read iterate all position base
			for(int j = 0; j < sam[i].size; j++)
			{
				count[sam[i].pos + j]++;
				lenght[sam[i].pos + j] += sam[i].size;
			}
		}
	}

	for(int i = 0; i < ln; i++)
	{
		if (count[i] == 0) continue;
		lenght[i] = lenght[i] /count[i];
	}
	wig(count, ln, pathCoverage);
	wig(lenght, ln, pathAvLenght);
	count.clear();
	lenght.clear();
	cout << "DONE." << endl;

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
			for (int j = 0; j < sam[i].lenght; j++)
			{
				bases[sam[i].pos + j] += sam[i].match[j];
			}
			for (int j = 0; j < sam[i-1].lenght; j++)
			{
				bases[sam[i-1].pos + j] += sam[i-1].match[j];
			}
		}
	}

	wig(bases, ln, path);
	bases.clear();
	cout << "DONE." << endl;
}

/* La funzione esegue correttamente il calcolo se:
- viene caricato un file sam che NON contiene /n dopo il qname
  e le read sono accoppiate una dopo l'altra
  */
void wrongMate(const char *path)
{
	cout << "Analisys of WRONG mate pair orientation ... " << flush;
	vector<int> base(ln, 0);

	// iterate all the reads
	int i = 0;
	while (i < sam.size()-1)
		if(sam[i].valid() && sam[i+1].valid()
			&& !strcmp(sam[i].name, sam[i+1].name))
			//&& sam[i].vFlag[2] == 0 
			//&& sam[i].vFlag[3] == 0 )
		{
			if(sam[i].vFlag[4] == sam[i+1].vFlag[4])
			{
				for(int j = 0; j < sam[i].lenght; j++)
				{
					base[sam[i].pos + j] += sam[i].match[j];
				}
				for(int j = 0; j < sam[i+1].lenght; j++)
				{
					base[sam[i+1].pos + j] += sam[i+1].match[j];
				}
			}
			i += 2;
		}
		else
			i++;
	wig(base, ln, path);
	base.clear();
	cout << "DONE." << endl;
}

void singleMate(const char *path)
{
	cout << "Analisys of SINGLE mate ... " << flush;
	vector<int> base(ln, 0);

	int i = 0;
	for(int i = 0; i < sam.size(); i++)
	{
		if(sam[i].valid() && sam[i].vFlag[3] == 1)
		{
			base[sam[i].pos]++;
		}
	}
	wig(base, ln, path);
	base.clear();
	cout << "DONE." << endl;
}


int main(int argc, char *argv[])
{
	loadData("../bioinfo/illumina.sam", 2);
	//multiCoverage("multiple.wig");
	//st_dev();
	physicalCov("physicalCoverage.wig");
	sequenceCov("sequenceCoverage.wig", "averageLenght.wig");
	wrongMate("wrong.wig");
	
	singleMate("single.wig");

	return 0; 
}