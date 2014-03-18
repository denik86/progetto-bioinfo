// main.cpp
// Il primo esempio in C++
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
using namespace std;

int main(int argc, char *argv[])
{
	ifstream file ("../bioinfo/illumina.sam");
	char delim[] = "\t";
	string line;
	if(file.is_open())
	{
		int i = 0;
		cout << "file letto e aperto" << endl;
		while(getline (file,line) && i < 2)
		{
			cout << "\nLinea " << i << ": ";
			cout << line << endl;
			char * str = new char [line.length()+1];
 			strcpy (str, line.c_str());

 			str = strtok (str, delim);
 			
 			printf("Lista degli elementi:\n");
		 	while(str != NULL)
		 	{

		 		printf("%s\n", str);
		 		str = strtok (NULL, delim); 
		 	}
			i++;
		}
		file.close();
	}
	return 0;
}