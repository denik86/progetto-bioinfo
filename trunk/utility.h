#include <fstream>
#include <string.h>
using namespace std;

// Trasforma un "string" in un "char[]"
char * stc(string s)
{
	char * str = new char [s.length()+1];
 	strcpy (str, s.c_str());
 	return str;
}

void wig(vector<int> array, int lenght, const char *path)
{
	ofstream file;
	file.open(path);
	file << "fixedStep chrom=Chromosome start=1 step=1" << endl;
	
	for(int i = 0; i < lenght; i++)
	{
		file << array[i] << endl; 
	}
	file.close();
}

void wig(vector<float> array, int lenght, const char *path)
{
	ofstream file;
	file.open(path);
	file << "fixedStep chrom=Chromosome start=1 step=1" << endl;
	
	for(int i = 0; i < lenght; i++)
	{
		file << array[i] << endl; 
	}
	file.close();
}


// Restituisce una stringa con il binario del numero number
string decToBin(int number)
{
	if ( number == 0 ) return "0";
    if ( number == 1 ) return "1";

    if ( number % 2 == 0 )
        return decToBin(number / 2) + "0";
    else
        return decToBin(number / 2) + "1";
}