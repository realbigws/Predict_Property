#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <vector> 
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <cmath>
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

//-------- load FASTA file -------//
int Read_FASTA_SEQRES(string &seqfile,string &seqres,int skip=1) //->from .fasta file
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(seqfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such file! %s \n",seqfile.c_str());
		exit(-1);
	}
	//skip
	int i;
	for(i=0;i<skip;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file bad! %s \n",seqfile.c_str());
			exit(-1);
		}
	}
	//process
	temp="";
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		temp+=buf;
	}
	seqres=temp;
	//return
	return (int)seqres.length();
}


//-------- load predicted results -------------//
/*
#-> 123
-1 -> 0.909218 0.090782 -> 0.909218  0
-1 -> 0.958281 0.041719 -> 0.958281  0
-1 -> 0.970875 0.029125 -> 0.970875  0
-1 -> 0.974167 0.025833 -> 0.974167  0
...
*/
int Parse_String(string &in,vector <double> &out)
{
	string buf;
	istringstream www(in);
	out.clear();
	for(;;)
	{
		if(!getline(www,buf,' '))break;
		if(buf=="")continue;
		if(buf=="-")break;
		double value=atof(buf.c_str());
		out.push_back(value);
	}
	return (int)out.size();
}
int Load_Predicted_Results(string &input_file,vector <vector <double> > &out_reso)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such file! %s \n",input_file.c_str());
		exit(-1);
	}
	//process
	out_reso.clear();
	int first=1;
	int first_num;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='#')continue;
		istringstream www(buf);
		getline(www,temp,'>');
		getline(www,temp,'>');
		//--- parse ----//
		vector <double> out;
		int retv=Parse_String(temp,out);
		if(first==1)
		{
			first_num=retv;
			first=0;
		}
		else
		{
			if(retv!=first_num)
			{
				fprintf(stderr,"file %s format bad at %s \n",input_file.c_str(),buf.c_str());
				exit(-1);
			}
		}
		//--- push_back ---//
		out_reso.push_back(out);
	}
	//return
	return (int)out_reso.size();
}


//-------- output to DISOPRED format ---------//
/*
#AUCpreD: order/disorder state prediction results
#Disorder residues are marked with asterisks (*) above threshold 0.200
# Ordered residues are marked with dots (.) below threshold 0.200

   1 M * 0.988
   2 S * 0.977
   3 Q * 0.970
...
*/
void Output_DisoPred_Format(string &seq_file,string &reso_file,double thres,int mode)
{
	//-> load seq_file
	string sequence;
	int seq_len=Read_FASTA_SEQRES(seq_file,sequence);
	//-> load reso_file
	vector <vector <double> > out_reso;
	int reso_len=Load_Predicted_Results(reso_file,out_reso);
	//check length
	if(seq_len!=reso_len)
	{
		fprintf(stderr,"seq_len %d not equal to reso_len %d \n",seq_len,reso_len);
		exit(-1);
	}
	//mode
	string mode_str="";
	if(mode==0) //profile mode
	{
		mode_str="by profile mode";
	}
	else
	{
		mode_str="by no_profile mode";
	}
	//output
	printf("#AUCpreD: order/disorder state prediction results %s \n",mode_str.c_str());
	printf("#Disorder residues are marked with asterisks (*) above threshold %5.3f \n",thres);
	printf("# Ordered residues are marked with dots (.) below threshold %5.3f \n",thres);
	for(int i=0;i<seq_len;i++)
	{
		char dot;
		double value=1-out_reso[i][0];
		//-> sigmoid transform
		double value_= 1.0 - 1.0/(1.0 + exp( 1.0*(value-0.23)/0.05 ) );
		if(value_>1)value_=1;
		if(value_<0)value_=0;
		if(value_>thres)dot='*';
		else dot='.';
		printf("%4d %c %c %5.3f \n",i+1,sequence[i],dot,value_);
	}
}


//------------ main -------------//
int main(int argc, char** argv)
{
	//---- DisoPred_Trans ----//
	{
		if(argc<5)
		{
			fprintf(stderr,"DisoPred_Trans <seq_file> <reso_file> <threshold> <mode> \n");
			exit(-1);
		}
		string seq_file=argv[1];
		string reso_file=argv[2];
		double thres=atof(argv[3]);
		int mode=atoi(argv[4]);
		//process
		Output_DisoPred_Format(seq_file,reso_file,thres,mode);
		//exit
		exit(0);
	}
}
