#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

//----------- get PSM ----------//
int TPL_TGT_Get_PSM(string &infile,FILE *fp)
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
		if(temp=="//////////// Original PSM file")break;
	}
	if(!getline(fin,buf,'\n'))
	{
		fprintf(stderr,"%s format bad!\n",infile.c_str());
		exit(-1);
	}
	//PSM
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))
		{
			printf("%s format bad!\n",infile.c_str());
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
	//---- TPL_TGT_Get_PSM ----//
	{
		if(argc<3)
		{
			fprintf(stderr,"TPL_TGT_Get_PSM <tpl/tgt_file> <outfile>\n");
			exit(-1);
		}
		string file=argv[1];
		string outfile=argv[2];
		FILE *fp=fopen(outfile.c_str(),"wb");
		TPL_TGT_Get_PSM(file,fp);
		fclose(fp);
		exit(0);
	}
}

