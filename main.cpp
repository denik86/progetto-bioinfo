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

		//print info
		cout << "** Total reads = " << sam.size() << endl;
		cout << "** Mapped reads = " << mappedReads << endl;
		cout << "** Not mapped reads = " << notMappedReads << endl;
		cout << "** Mean size = " << mean << "\n" << endl;

		file.close();
	}
	else {
		cout << "Error file" << endl;
	}
}

void st_dev()
{
	if(mean == 0)
	{
		cout <<"ERROR: La media è nulla" << endl;
	}
	else
	{
		long int sum_st_dev = 0;
		int count = 0;
		for(int i = 0; i < sam.size(); i++)
		{
			if(sam[i].valid() && sam[i].size > 0 &&
				sam[i].size < 2*mean )
				 {
				 	sum_st_dev += pow(abs(sam.at(i).size) - mean  ,2);
				 	count++;
				 }
		}
		sd = sqrt(sum_st_dev / count);
		cout << "** Standard deviation = " << sd << endl;
	}
}

void physicalCov(const char * path)
{
	cout << "Analisys of physical coverage ... " << flush;
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
	cout << "DONE.\n" << endl;
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
	cout << "DONE.\n" << endl;

}

/* This function IMPLIED that the reads in the sam file
 * are necessarily sorted by name */
void multiCoverage(const char * path)
{
	cout << "Analisys of MULTIPLE READ coverage ... " << flush;
	vector<int> bases(ln, 0);
	for(int i = 0; i < sam.size(); i++)
	{
		if(i > 0 && strcmp(sam[i].name, sam[i-1].name) == 0 &&
			sam[i].num == sam[i-1].num)
		{ // two reads are the same, also in the number /n
			for (int j = 0; j < sam[i].lenght; j++)
				bases[sam[i].pos + j] += sam[i].match[j];
			for (int j = 0; j < sam[i-1].lenght; j++)
				bases[sam[i-1].pos + j] += sam[i-1].match[j];
		}
	}

	wig(bases, ln, path);
	bases.clear();
	cout << "DONE.\n" << endl;
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
		{
			if(sam[i].vFlag[4] == sam[i+1].vFlag[4])
			{  // reads are in same directions
				for(int j = 0; j < sam[i].lenght; j++)
					base[sam[i].pos + j] += sam[i].match[j];
				for(int j = 0; j < sam[i+1].lenght; j++)
					base[sam[i+1].pos + j] += sam[i+1].match[j];
			}
			i += 2;
		}
		else
			i++;
	wig(base, ln, path);
	base.clear();
	cout << "DONE.\n" << endl;
}

void singleMate(const char *path)
{
	cout << "Analisys of SINGLE mate ... " << flush;
	vector<int> base(ln, 0);

	// check the 3 bit, that indicates if his
	// mate is mapped
	for(int i = 0; i < sam.size(); i++)
		if(sam[i].valid() && sam[i].vFlag[3] == 1)
			base[sam[i].pos]++;

	wig(base, ln, path);
	base.clear();
	cout << "DONE.\n" << endl;
}


int main(int argc, char *argv[])
{
	bool caricato = false;
	bool fine = false;
	cout << "\n*****************************************" << endl;
	cout << "******** Resequencing Project ***********" << endl;
	cout << "*****************************************" << endl;
	while(!fine)
	{
		cout << "----- SEGLI UNA FUNZIONE: (qualunque lettera per terminare) ----" << endl;
		if(!caricato) cout << "ATTENZIONE: CARICARE PRIMA UN FILE SAM (1)!!"<< endl;
		int in;
		char *path = new char[100];
		char *path2 = new char[100];
		cout << "1. Carica File SAM" << endl;
		cout << "2. Deviazione Standard" << endl;
		cout << "3. Physical Coverage (file .wig)" << endl;
		cout << "4. Sequence Coverage / Average Lenght (2 file .wig)" << endl;
		cout << "5. Read multipli (file.wig)" << endl;
		cout << "6. Stesso orientamento mate (file. wig)" << endl;
		cout << "7. Mate mappati singolarmente (file .wig)" << endl;
		cout << "COMANDO: ";
		cin >> in;

		switch(in)
		{
			case 1:
				cout << "Scegli il percorso del file: ";
				cin >> path;
				loadData(path, 2);
				caricato = true;
				break;

			case 2:
				st_dev();
				break;
		
			case 3:
				cout << "Scegli il percorso del file .wig: ";
				cin >> path;
				physicalCov(path);
				break;
			
			case 4:
				cout << "Scegli il percorso del file (Sequence coverage) .wig: ";
				cin >> path;
				cout << "Scegli il percorso del file (Average Lenght) .wig: ";
				cin >> path2;
				sequenceCov(path, path2);
				break;
				
			case 5:
				cout << "Scegli il percorso del file .wig: ";
				cin >> path;
				multiCoverage(path);
				break;

			case 6:
				cout << "Scegli il percorso del file .wig: ";
				cin >> path;
				wrongMate(path);
				break;

			case 7:
				cout << "Scegli il percorso del file .wig: ";
				cin >> path;
				singleMate(path);
				break;
	
			default:
				cout << "Programma Terminato\n" << endl;
				fine = true;
				break;
		}
	}
	return 0; 
}