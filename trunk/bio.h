// Classes of structure used in the project

class SamLine{
public:
	// Info reads
	char *name;
	int num;
	int flag;
	char *chromo;
	int pos;
	int mapq;
	string cigar;
	char *mrnm;
	int mpos;
	int size;
	char *seq;
	char *qual;

	int lenght;
	vector<int> match;

	int vFlag[12];

	void readInfo(string l);
	void flagToBinaryVector();
	bool valid();
	int cigarValues();
};

// Verify that read is valid (cigar != *)
bool SamLine::valid()
{
	return strcmp(chromo, "*");
}

void SamLine::readInfo(string l)
{
	char * line = stc(l);
	name = strtok(line, "\t");
	flag = atoi(strtok(NULL, "\t"));
	flagToBinaryVector();
	chromo = strtok(NULL, "\t");
	pos = atoi(strtok(NULL, "\t"));
	mapq = atoi(strtok(NULL, "\t"));
	cigar = strtok(NULL, "\t");
	mrnm = strtok(NULL, "\t");
	mpos = atoi(strtok(NULL, "\t"));
	size = atoi(strtok(NULL, "\t"));
	seq = strtok(NULL, "\t");
	qual = strtok(NULL, "\t");
	
	name = strtok(name, "/");
	num = atoi(strtok(name, "/"));

	lenght = cigarValues();

}

void SamLine::flagToBinaryVector()
{
	for (int i = 0; i < 12; ++i) 
	{  // assuming a 7 bit int
    	vFlag[i] = flag & (1 << i) ? 1 : 0;
    }
}

int SamLine::cigarValues()
{ 
	int count = 0;
	for(int i = 0; i < cigar.size(); i++)
	{
		if (isalpha(cigar[i]))
		{
			char letter = cigar[i];
		  	int number = atoi(cigar.substr(i-count, i).c_str());
		  	
		  	if(letter == 'M')
		  	{ 
		  		for(int j = 0; j < number; j++)
		  		{
		  			match.push_back(1);
		  		}
		  	}
		  	else if(letter == 'D')
		  	{
		  		for(int j = 0; j < number; j++)
		  		{
		  			match.push_back(0);
		  		}
		  	}
		  	else {  /* non fai nulla */ }
		  	count = 0;
		}
		else 
			count++;

	}
	return match.size();
	
}
