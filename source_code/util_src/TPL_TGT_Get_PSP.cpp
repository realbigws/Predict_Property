#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

//----------- get PSP ----------//
int TPL_TGT_Get_PSP(string &infile,FILE *fp)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"%s not found!\n",infile.c_str());
		exit(-1);
	}
	//skip1
	int len;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"%s format bad!\n",infile.c_str());
			exit(-1);
		}
		len=(int)buf.length();
		if(len<30)continue;
		temp=buf.substr(0,30);
		if(temp=="//////////// Original PSP file")break;
	}
	if(!getline(fin,buf,'\n'))
	{
		fprintf(stderr,"%s format bad!\n",infile.c_str());
		exit(-1);
	}
	//init
	fprintf(fp,"\n");
	fprintf(fp,"Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts\n");
	fprintf(fp,"           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V\n");
	//PSP
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"%s format bad!\n",infile.c_str());
			exit(-1);
		}
		len=(int)buf.length();
		if(len<10)break;
		fprintf(fp,"%s\n",buf.c_str());
		count++;
	}
	fprintf(fp,"\n");
	return count;
}

//------ main ------//
int main(int argc,char **argv)
{
	//---- TPL_TGT_Get_PSP ----//
	{
		if(argc<3)
		{
			fprintf(stderr,"TPL_TGT_Get_PSP <tpl/tgt_file> <outfile>\n");
			exit(-1);
		}
		string file=argv[1];
		string outfile=argv[2];
		FILE *fp=fopen(outfile.c_str(),"wb");
		TPL_TGT_Get_PSP(file,fp);
		fclose(fp);
		exit(0);
	}
}
