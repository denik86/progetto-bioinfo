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
	char *cigar;
	char *mrnm;
	int mpos;
	int size;
	char *seq;
	char *qual;

	int vFlag[12];

	void readInfo(string l);
	void flagToBinaryVector();
	bool valid();
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


}

void SamLine::flagToBinaryVector()
{
	for (int i = 0; i < 12; ++i) 
	{  // assuming a 7 bit int
    	vFlag[i] = flag & (1 << i) ? 1 : 0;
    }
}
