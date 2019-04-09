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


//---------- Membrane Protein Labels ----------//
/*
1  -> Side1
2  -> Side2
B  -> Beta-strand
H  -> alpha-helix
C  -> coil
I  -> membrane-inside
L  -> membrane-loop
F  -> interfacial helix
U  -> unknown localizations
*/
char TM8_To_Char(int in)
{
	switch(in)
	{
		case 0: return 'H'; // H
		case 1: return 'E'; // E
		case 2: return 'C'; // C
		case 3: return 'I'; // I
		case 4: return 'L'; // L
		case 5: return 'F'; // F
		case 6: return 'X'; // X
		case 7: return '_'; // _ (=1+2)
		default: return '_';
	}
}
void TM9_To_TM8_Prob(vector <double> &in, vector <double> &out)
{
	out.resize(8);
	out[0]=in[3];
	out[1]=in[2];
	for(int i=2;i<=6;i++)out[i]=in[i+2];
	out[7]=in[0]+in[1];
}
int Return_Label(vector <double> &in,double &maxval)
{
	int i;
	int size=(int)in.size();
	int label=0;
	double wsmax=0;
	for(i=0;i<size;i++)
	{
		if(in[i]>wsmax)
		{
			wsmax=in[i];
			label=i;
		}
	}
	maxval=wsmax;
	return label;
}

//-------- determine transmembrane type ------//
int Determine_Transmembrane_Type(vector <vector <double> > &in,double thres)
{
	int alpha_type=0;
	int beta_type=0;
	for(int i=0;i<(int)in.size();i++)
	{
		if(in[i][3]>thres)alpha_type++;
		if(in[i][2]>thres)beta_type++;
	}
	if(alpha_type==0 && beta_type==0)return 0; //-> NON-TM
	if(alpha_type>beta_type)return 1; //-> Alpha Type
	else return 2; //-> Beat Type
}


//-------- output to TopoPred format ---------//
/*
#TopoPred: transmembrane topology prediction results
#Transmembrane residues are marked with asterisks (*) above threshold 0.500
# Non-transmembrane residues are marked with dots (.) below threshold 0.500

   1 M * 0.988
   2 S * 0.977
   3 Q * 0.970
...
*/
void Output_TopoPred_Format(string &seq_file,string &reso_file,
	double thres,int mode, int type)
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
	if(type==0) //-> 2-state
	{
		//--- determine type ---//
		char type_char;
		int TM_Type=Determine_Transmembrane_Type(out_reso,thres);
		if(TM_Type==0)type_char='_';
		else if(TM_Type==1)type_char='H';
		else type_char='B';
		//--- output to file ---//
		printf("#TopoPred_TM2: 2-state transmembrane topology prediction results %s \n",mode_str.c_str());
		printf("#Transmembrane residues are marked with 'H' or 'E' above threshold %5.3f \n",thres);
		printf("# Non-transmembrane residues are marked with '_' below threshold %5.3f \n",thres);
		for(int i=0;i<seq_len;i++)
		{
			char dot;
			double value;
			if(TM_Type==0)value=out_reso[i][3]>out_reso[i][2]?out_reso[i][3]:out_reso[i][2];
			else if(TM_Type==1)value=out_reso[i][3];
			else value=out_reso[i][2];
			if(value>thres)dot=type_char;
			else dot='_';
			printf("%4d %c %c %5.3f \n",i+1,sequence[i],dot,value);
		}
	}
	else        //-> 8-state
	{
		printf("#TopoPred_TM8: 8-state transmembrane topology prediction results %s \n",mode_str.c_str());
		printf("#probabilities are in the order of H E C I L F X _, the 8 transmembrane topology types used in PDBTM \n");
		for(int i=0;i<seq_len;i++)
		{
			int state_num=8;
			//transfer to 8-TM
			vector <double> output_tm;
			TM9_To_TM8_Prob(out_reso[i],output_tm);
			//return maximal value
			double maxval;
			int label=Return_Label(output_tm,maxval);
			char lab_c=TM8_To_Char(label);
			//output
			printf("%4d %c %c   ",i+1,sequence[i],lab_c);
			for(int j=0;j<state_num;j++)printf("%5.3f ",output_tm[j]);
			printf("\n");
		}
	}
}


//------------ main -------------//
int main(int argc, char** argv)
{
	//---- TopoPred_Trans ----//
	{
		if(argc<6)
		{
			fprintf(stderr,"TopoPred_Trans <seq_file> <reso_file> <threshold> <mode> <type> \n");
			fprintf(stderr,"[note]: type=0 for 2-state; type=1 for 8-state \n");
			exit(-1);
		}
		string seq_file=argv[1];
		string reso_file=argv[2];
		double thres=atof(argv[3]);
		int mode=atoi(argv[4]);
		int type=atoi(argv[5]);
		//process
		Output_TopoPred_Format(seq_file,reso_file,thres,mode,type);
		//exit
		exit(0);
	}
}

