#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <vector>
#include <string>
#include "seq.h"
#include "template.h"
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

//-------- load SSE file ----------//
//file format
/*
#DeepConCNF_SS8: eight-class secondary structure prediction results
#probabilities are in the order of H G I E B T S L(loops), the 8 secondary structure types used in DSSP

   1 E L   0.000 0.000 0.000 0.003 0.000 0.000 0.001 0.995
   2 N L   0.000 0.001 0.000 0.386 0.009 0.004 0.001 0.598
   3 I E   0.000 0.001 0.000 0.901 0.005 0.002 0.012 0.078
   4 E E   0.000 0.001 0.000 0.971 0.001 0.000 0.002 0.024
   5 V E   0.000 0.001 0.000 0.983 0.001 0.001 0.001 0.013
   6 H E   0.001 0.001 0.000 0.946 0.002 0.001 0.004 0.046
   7 M E   0.001 0.002 0.000 0.819 0.001 0.003 0.006 0.167
....
*/

int Load_Prob_File(string &prob_file,int lab_num,vector < vector <double> > &prob_out)
{
	//start
	ifstream fin;
	string buf,temp;
	//read
	fin.open(prob_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",prob_file.c_str());
		exit(-1);
	}
	//skip
	for(;;)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file %s format bad!\n",prob_file.c_str());
			exit(-1);
		}
		if(buf=="")break;
		if(buf[0]=='#')continue;
	}
	//load
	int count=0;
	prob_out.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		for(int i=0;i<3;i++)www>>temp;  //skip
		vector <double> tmp_rec;
		for(int i=0;i<lab_num;i++)
		{
			double tmp_prob;
			www>>tmp_prob;
			tmp_rec.push_back(tmp_prob);
		}
		prob_out.push_back(tmp_rec);
		count++;
	}
	//return
	return count;
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
X  -> unassigned
*/

//-------- Memb<->Int --------//
int Memb_To_Int(char c)
{
	switch(c)
	{
		case '1': return 0;
		case '2': return 1;
		case 'B': return 2;
		case 'H': return 3;
		case 'C': return 4;
		case 'I': return 5;
		case 'L': return 6;
		case 'F': return 7;
		default: return 8;
	}
}
char Int_To_Memb(int c)
{
	switch(c)
	{
		case 0: return '1';
		case 1: return '2';
		case 2: return 'B';
		case 3: return 'H';
		case 4: return 'C';
		case 5: return 'I';
		case 6: return 'L';
		case 7: return 'F';
		default: return 'X';
	}
}

int AA4SUB[26] = { 1, -1,  5,  4,  7, 14,  8, 9,  10,-1, 12,11,  13, 3, -1, 15,  6,  2,  16,  17,-1, 20, 18, -1, 19, -1};

//------ load label file -------//
//-> example
/*
>1a0sP
SGFEFHGYARSGVIMNDSGASTKSGAYITPAGETGGAIGRLGNQADTYVEMNLEHKQTLDNGATTRFKVMVADGQTSYNDWTASTSDLNVRQAFVELGNLPTFAGPFKGSTLWAGKRFDRDNFDIHWIDSDVVFLAGTGGGIYDVKWNDGLRSNFSLYGRNFGDIDDSSNSVQNYILTMNHFAGPLQMMVSGLRAKDNDERKDSNGNLAKGDAANTGVHALLGLHNDSFYGLRDGSSKTALLYGHGLGAEVKGIGSDGALRPGADTWRIASYGTTPLSENWSVAPAMLAQRSKDRYADGDSYQWATFNLRLIQAINQNFALAYEGSYQYMDLKPEGYNDRQAVNGSFYKLTFAPTFKVGSIGDFFSRPEIRFYTSWMDWSKKLNNYASDDALGSDGFNSGGEWSFGVQMETWF
222222BBBBBBBBB1111111111111111111IIIIIIIIIIIIBBBBBBB22222222222222BBBBBBBBBBBBBBBB11BBBBBBBBBB22222222222222222BBBBBBIIIIIIIIIIIIIIIBBBBBBBBB22222222222BBBBBBB11111111111111BBBBBBB22222BBBBBBB111111111111111111111111BBBBBBB2222222222222BBBBBBB1111111111111111111111BBBBBBB222222222BBBBBBBB11111111111111BBBBBBB222222222BBBBBBBB1111111111111111111BBBBBBBB22222222222222BBBBBBB11111111111111111111111111BBBBBBBBBB2
*/
int Load_LAB_File(string &cn_file,vector <int> &lab_number)
{
	//start
	ifstream fin;
	string buf,temp;
	//read
	fin.open(cn_file.c_str(), ios::in);
	if(fin.fail()!=0)return -1;
//	{
//		fprintf(stderr,"file %s not found!\n",cn_file.c_str());
//		exit(-1);
//	}
	//skip
	for(int i=0;i<3;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file %s format bad!\n",cn_file.c_str());
			exit(-1);
		}
	}
	//load
	lab_number.clear();
	int count=0;
	for(int i=0;i<(int)buf.length();i++)
	{
		char c=buf[i];
		int lab=Memb_To_Int(c);
		lab_number.push_back(lab);
	}
	//return
	return (int)buf.length();
}


//------------ Load TGT, SS3, Label, output MemProt feature -------//
void MemProt_Feat(string &tgt_file, string &sse_file, string &label_file)
{
	//-- read tgt ---//
	string tgt_name,tgt_root;
	getBaseName(tgt_file,tgt_name,'/','.');
	getRootName(tgt_file,tgt_root,'/');
	SEQUENCE *s=new SEQUENCE(tgt_name,tgt_root,1,1);
	//-- load sse files --//
	int sse_label=8;
	vector < vector <double> > sse_prob_out;
	int sse_len=Load_Prob_File(sse_file,sse_label,sse_prob_out);
	//-- load label string --//
	vector <int> lab_number;
	int label_len=Load_LAB_File(label_file,lab_number);
	if(label_len==-1)
	{
		label_len=sse_len;
		lab_number.resize(label_len);
		for(int k=0;k<label_len;k++)lab_number[k]=-1;
	}
	//-- check --//
	if(sse_len!=s->length || label_len!=s->length )
	{
		fprintf(stderr,"s->length %d not equal to sse_len %d or label_len %d \n",s->length,sse_len,label_len);
		delete s;
		exit(-1);
	}

	//======== process feature =======//
	vector <string> output;
	output.clear();
	for(int k=0;k<s->length;k++)
	{
		//----- feature -----//
		vector <double> features;
		features.clear();
		//-- profile realted --//
		//emission score
		for(int i=0;i<20;i++)
		{
			float template_amino_Score = 1.0*s->EmissionProb[k][i];
			features.push_back(template_amino_Score);
		}
		//PSM score
		for(int i=0;i<20;i++)
		{
			double template_amino_Prob = 1.0*s->PSM[k][i];
			template_amino_Prob=1.0/(1.0+exp(-1.0*template_amino_Prob));
			features.push_back(template_amino_Prob);
		}
		//-- amino acid realted --//
		//zy ami model
		int pos=AA4SUB[s->sequence[k]-'A']-1;
		vector <int> aaind(20,0);
		if(pos<0 || pos>=20)
		{
			for(int i=0;i<20;i++)
			{
				features.push_back(aaind[i]);
			}
		}
		else
		{
			aaind[pos]=1;
			for(int i=0;i<20;i++)
			{
				features.push_back(aaind[i]);
			}
		}
		//-- secondary structure ---//
		for(int i=0;i<sse_label;i++)
		{
			double psipred_reso = 1.0*sse_prob_out[k][i];
			features.push_back(psipred_reso);
		}

		//------ output ------//
		int featdim=(int)features.size();
		stringstream oss;
		for(int i=0;i<featdim;i++)
		{
			int wsiii=(int)features[i];
			if(wsiii!=features[i])oss << features[i] << " ";
			else oss << wsiii << " ";
		}
		string wsbuf=oss.str();
		output.push_back(wsbuf);
	}
	delete s;

	//====== output =======//
	printf("%d\n",label_len);
	for(int k=0;k<label_len;k++)printf("%s\n",output[k].c_str());
	for(int k=0;k<label_len;k++)printf("%d\n",lab_number[k]);
}

//----------- main -------------//
int main(int argc,char **argv)
{
	//------- Process Membrane Proteins Features --------// 
	{
		if(argc<4)
		{
			printf("MemProt_Feat <tgt_file> <ss8_file> <label_file> \n");
			printf("[note]: the label_file should be three lines !!! \n");
			exit(-1);
		}
		string tgt_file=argv[1];
		string ss8_file=argv[2];
		string label_file=argv[3];
		//process
		MemProt_Feat(tgt_file,ss8_file,label_file);
		exit(0);
	}
}
