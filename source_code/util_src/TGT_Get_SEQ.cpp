#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}

//----------- get SEQ ----------//
int TGT_Get_SEQ(string &infile,FILE *fp)
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
	string name;
	getBaseName(infile,name,'/','.');
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
		if(len<8)continue;
		temp=buf.substr(0,10);
		if(temp=="Sequence =")
		{
			istringstream www(buf);
			www>>temp>>temp>>temp;
			fprintf(fp,">%s\n",name.c_str());
			fprintf(fp,"%s\n",temp.c_str());
			exit(0);
		}
	}
}

//------ main ------//
int main(int argc,char **argv)
{
	//---- TGT_Get_SEQ ----//
	{
		if(argc<3)
		{
			fprintf(stderr,"TGT_Get_SEQ <tgt_file> <outfile>\n");
			exit(-1);
		}
		string file=argv[1];
		string outfile=argv[2];
		FILE *fp=fopen(outfile.c_str(),"wb");
		TGT_Get_SEQ(file,fp);
		fclose(fp);
		exit(0);
	}
}
