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

//------- return label ---------//
/*
Secondary structure labels, with the sequence of 'L', 'B', 'E', 'G', 'I', 'H', 'S', 'T'
*/
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
int SS8_To_SS3(int in)
{
	switch(in) 
	{
		case 0: return 0; // H->H
		case 1: return 2; // G->C
		case 2: return 2; // I->C
		case 3: return 1; // E->E
		case 4: return 2; // B->C
		case 5: return 2; // T->C
		case 6: return 2; // S->C
		case 7: return 2; // L->C
		default: return 2;
	}
}
void SS8_To_SS3_Prob(vector <double> &in, vector <double> &out)
{
	out.resize(3);
	out[0]=in[0];
	out[1]=in[3];
	out[2]=in[5]+in[6]+in[7]+in[1]+in[2]+in[4];
}

//--- SS to Char ---//
char SS8_To_Char(int in)
{
	switch(in) 
	{
		case 0: return 'H'; // H
		case 1: return 'G'; // G
		case 2: return 'I'; // I
		case 3: return 'E'; // E
		case 4: return 'B'; // B
		case 5: return 'T'; // T
		case 6: return 'S'; // S
		case 7: return 'L'; // L
		default: return 'L';
	}
}
char SS3_To_Char(int in)
{
	switch(in) 
	{
		case 0: return 'H'; // H
		case 1: return 'E'; // E
		case 2: return 'C'; // C
		default: return 'C';
	}
}
char ACC_To_Char(int in)
{
	switch(in) 
	{
		case 0: return 'B'; // B
		case 1: return 'M'; // M
		case 2: return 'E'; // E
		default: return 'M';
	}
}

//=========== Label Parser ===============//
//-> lab_type
//   [0=SS8], [1=SS3], [2=ACC], [3=CN]
void Label_Parser(string &seq_file,string &pred_file,int SS_type)
{
	//-> load seq
	int seq_len;
	string ami_seq;
	seq_len=Read_FASTA_SEQRES(seq_file,ami_seq);
	//-> load pred
	int pred_len;
	vector <vector <double> > output;
	pred_len=Load_Predicted_Results(pred_file,output);
	//-> check
	if(seq_len!=pred_len)
	{
		fprintf(stderr,"seq_len %d not equal to pred_len %d \n",
			seq_len,pred_len);
		exit(-1);
	}

	//-> output header
	if(SS_type==0)      //--> output SS8
	{
		printf("#DeepConCNF_SS8: eight-class secondary structure prediction results \n");
		printf("#probabilities are in the order of H G I E B T S L(loops), the 8 secondary structure types used in DSSP \n\n");
	}
	else if(SS_type==1) //--> output SS3
	{
		printf("#DeepConCNF_SS3: three-class secondary structure prediction results \n");
		printf("#probabilities are in the order of H E C, the 3 secondary structure types used in DSSP \n\n");
	}
	else if(SS_type==2) //--> output ACC
	{
		printf("#DeepConCNF_SAS: three-state solvent accessibility prediction results \n");
		printf("#probabilities are in the order of B (Bury, pACC: 0-10), M (Medium, pACC: 11-40) and E (Exposed, pACC: 41-100), \n");
		printf("#where pACC is the relative solvent accessibility value calculated by DSSP \n");
	}
	else                //--> output CN
	{
		printf("# AcconPred: 15-state Cb contact number prediction results \n");
		printf("# probabilities are in the order of 0 to 14 of the Cb contact number. \n");
		printf("# The original training data are made from Cb with separation=6 (means +-5 residues are excluded), and distance=7.5 (means Cb-Cb distance should within 7.5A) \n");
		printf("# output format: <position> <one letter amino acid> <predicted label> ... <probabilities> ... -> <maximal probability> \n");
	}

	//-> output content
	//   1 E L   0.000 0.000 0.000 0.000 0.000 0.000 0.000 1.000
	//   2 N L   0.000 0.000 0.000 0.379 0.004 0.000 0.000 0.617
	for(int i=0;i<(int)output.size();i++)
	{
		if(SS_type==0)      //--> output SS8
		{
			int state_num=8;
			//return maximal value
			double maxval;
			int label=Return_Label(output[i],maxval);
			char lab_c=SS8_To_Char(label);
			//output
			printf("%4d %c %c   ",i+1,ami_seq[i],lab_c);
			for(int j=0;j<state_num;j++)printf("%5.3f ",output[i][j]);
			printf("\n");
		}
		else if(SS_type==1) //--> output SS3
		{
			int state_num=3;
			//transfer to 3-SSE
			vector <double> output_sse;
			SS8_To_SS3_Prob(output[i],output_sse);
			//return maximal value
			double maxval;
			int label=Return_Label(output_sse,maxval);
			char lab_c=SS3_To_Char(label);
			//output
			printf("%4d %c %c   ",i+1,ami_seq[i],lab_c);
			for(int j=0;j<state_num;j++)printf("%5.3f ",output_sse[j]);
			printf("\n");
		}
		else if(SS_type==2) //--> output ACC
		{
			int state_num=3;
			//return maximal value
			double maxval;
			int label=Return_Label(output[i],maxval);
			char lab_c=ACC_To_Char(label);
			//output
			printf("%4d %c %c   ",i+1,ami_seq[i],lab_c);
			for(int j=0;j<state_num;j++)printf("%5.3f ",output[i][j]);
			printf("\n");
		}
		else               //--> output CN
		{
			int state_num=15;
			//return maximal value
			double maxval;
			int label=Return_Label(output[i],maxval);
			//output
			printf(" %d %c %d ",i+1,ami_seq[i],label);
			for(int j=0;j<state_num;j++)printf("%lf ",output[i][j]);
			printf(" -> %lf \n",maxval);
		}
	}
}


//------------ main -------------//
int main(int argc, char** argv)
{
	//---- Label_Parser ----//
	{
		if(argc<4)
		{
			fprintf(stderr,"Label_Parser <seq_file> <pred_file> <mode> \n");
			fprintf(stderr,"[note]: for <mode>: [0=SS8], [1=SS3], [2=ACC], [3=CN] \n");
			exit(-1);
		}
		string seq_file=argv[1];
		string pred_file=argv[2];
		int mode=atoi(argv[3]);
		//process
		Label_Parser(seq_file,pred_file,mode);
		//exit
		exit(0);
	}
}
