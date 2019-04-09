#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <iomanip>
#include <time.h>
#include "seq.h"
using namespace std;


//======================= I/O related ==========================//
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

//===================== ACC_Predict ==================//
double ACC_ZY_AmiModel[20][58]={
{0.42,0.8,0.9,0.85,1.05,1.06,1.18,0.73,0.93,0.77,0.87,1.28,0.05,1,0.31,6.11,0.42,0.23,1,0.58,0.514,0.266,0.822,0.709,0.463,0.736,0.723,0.962,0.956,0.513,0.971,0.976,0.861,0.71,0.88,0.949,0.959,0.968,10,6.6,0,0,0,6.6,9,0,0,0,5,6.6,5,0,0,0,0,0,0,0},
{0.61,0.87,0.78,0.65,1.05,1.08,1.06,1.09,1.27,0.86,0.8,2.34,0.29,6.13,-1.01,10.74,0.36,0.25,0.58,1,0.814,0.566,0.605,0.825,0.601,0.82,0.859,0.405,0.389,0.959,0.579,0.465,0.798,0.879,0.809,0.658,0.687,0.417,6.6,10,0,0,0,9,6.6,0,5,0,6.6,9,9,9,0,5,5,6.6,5,0},
{2.24,0.78,1.16,0.94,0.63,1.52,1.68,1.58,0.98,1.55,1.71,1.6,0.13,2.95,-0.6,6.52,0.21,0.22,0.514,0.814,1,0.844,0.65,0.93,0.808,0.926,0.883,0.292,0.269,0.821,0.568,0.393,0.822,0.92,0.806,0.641,0.63,0.31,0,0,10,9,0,0,0,0,5,0,0,0,0,0,0,6.6,0,0,0,0},
{2.56,1.07,2.06,1.95,0.47,0.8,0.78,1.65,0.73,2.06,1.59,1.6,0.11,2.78,-0.77,2.95,0.25,0.2,0.266,0.566,0.844,1,0.431,0.766,0.932,0.724,0.689,0.068,0.053,0.571,0.311,0.152,0.604,0.754,0.6,0.404,0.403,0.092,0,0,9,10,0,0,0,0,6.6,0,0,5,0,0,0,6.6,0,0,0,0},
{0.95,0.68,0.6,0.74,1.21,1.06,1.11,0.82,0.8,0.92,1.2,1.77,0.13,2.43,1.54,6.35,0.17,0.41,0.822,0.605,0.65,0.431,1,0.777,0.511,0.829,0.822,0.73,0.711,0.549,0.856,0.802,0.852,0.75,0.825,0.892,0.845,0.737,0,0,0,0,10,0,0,0,9,0,0,0,0,5,0,9,6.6,5,5,0},
{0.6,0.78,1.04,1.38,0.83,1,1.21,0.95,1.27,0.93,0.67,1.56,0.18,3.95,-0.22,5.65,0.36,0.25,0.709,0.825,0.93,0.766,0.777,1,0.821,0.954,0.922,0.54,0.52,0.819,0.752,0.619,0.93,0.951,0.923,0.814,0.809,0.556,6.6,9,0,0,0,10,6.6,0,5,0,9,9,5,6.6,0,0,5,5,5,0},
{0.58,0.93,1.94,1.63,0.4,0.93,0.8,0.83,1.12,0.93,0.74,1.56,0.15,3.78,-0.64,3.09,0.42,0.21,0.463,0.601,0.808,0.932,0.511,0.821,1,0.739,0.724,0.314,0.303,0.586,0.49,0.369,0.727,0.792,0.719,0.576,0.594,0.335,9,6.6,0,0,0,6.6,10,0,0,0,5,9,0,0,0,0,0,0,0,0},
{1.41,1.66,1.76,1.31,0.77,0.72,2.36,1.77,1.23,0.94,1.72,0,0,0,0,6.07,0.13,0.15,0.736,0.82,0.926,0.724,0.829,0.954,0.739,1,0.935,0.565,0.546,0.789,0.793,0.665,0.929,0.939,0.91,0.845,0.825,0.576,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,0},
{1.18,0.88,1.05,1.17,0.85,1.28,1.32,1.06,1.1,1.07,1.08,2.99,0.23,4.66,0.13,7.69,0.27,0.3,0.723,0.859,0.883,0.689,0.822,0.922,0.724,0.935,1,0.552,0.53,0.832,0.773,0.634,0.909,0.96,0.907,0.823,0.817,0.569,0,5,5,6.6,9,5,0,0,10,0,0,0,0,5,0,9,5,5,6.6,0},
{0.29,0.89,0.53,0.71,1.78,0.73,0.48,0.5,0.85,0.62,0.6,4.19,0.19,4,1.8,6.04,0.3,0.45,0.962,0.405,0.292,0.068,0.73,0.54,0.314,0.565,0.552,1,0.998,0.336,0.93,0.983,0.728,0.526,0.76,0.886,0.899,0.997,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,9},
{0.36,0.88,0.51,0.68,1.45,1.17,1.05,0.51,0.79,0.68,0.78,2.59,0.19,4,1.7,6.04,0.39,0.31,0.956,0.389,0.269,0.053,0.711,0.52,0.303,0.546,0.53,0.998,1,0.313,0.922,0.982,0.716,0.503,0.742,0.875,0.893,0.997,5,6.6,0,0,0,9,5,0,0,0,10,9,6.6,5,0,0,0,5,0,0},
{0.56,0.93,0.95,0.73,0.69,1.23,1.14,1.15,1.52,0.86,0.88,1.89,0.22,4.77,-0.99,9.99,0.32,0.27,0.513,0.959,0.821,0.571,0.549,0.819,0.586,0.789,0.832,0.336,0.313,1,0.523,0.388,0.727,0.873,0.785,0.584,0.614,0.347,6.6,9,0,5,0,9,9,0,0,0,9,10,5,5,0,5,5,0,0,0},
{0.42,0.82,0.52,0.72,1.35,1.11,1.07,0.65,0.95,0.67,0.77,2.35,0.22,4.43,1.23,5.71,0.38,0.32,0.971,0.579,0.568,0.311,0.856,0.752,0.49,0.793,0.773,0.93,0.922,0.523,1,0.97,0.891,0.733,0.877,0.981,0.974,0.932,5,9,0,0,0,5,0,0,0,0,6.6,5,10,6.6,0,0,0,6.6,0,0},
{0.52,1.07,0.67,0.91,1.52,1.07,0.9,0.66,0.9,0.74,0.72,2.94,0.29,5.89,1.79,5.67,0.3,0.38,0.976,0.465,0.393,0.152,0.802,0.619,0.369,0.665,0.634,0.983,0.982,0.388,0.97,1,0.8,0.6,0.799,0.937,0.938,0.984,0,9,0,0,5,6.6,0,0,5,0,5,5,6.6,10,0,5,5,9,9,0},
{1,1,1,1,1,1,1,1.84,0.63,2.15,1.05,2.67,0,2.72,0.72,6.8,0.13,0.34,0.861,0.798,0.822,0.604,0.852,0.93,0.727,0.929,0.909,0.728,0.716,0.727,0.891,0.8,1,0.898,0.927,0.941,0.94,0.74,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0},
{2.56,1.13,1.49,1.1,0.65,1.43,1.19,0.98,1.03,1.23,1.28,1.31,0.06,1.6,-0.04,5.7,0.2,0.28,0.71,0.879,0.92,0.754,0.75,0.951,0.792,0.939,0.96,0.526,0.503,0.873,0.733,0.6,0.898,1,0.923,0.792,0.789,0.542,0,5,6.6,6.6,9,0,0,0,9,0,0,5,0,5,0,10,6.6,5,6.6,0},
{2.18,1.11,1.14,1.45,0.93,1.2,0.7,1.02,1.14,0.91,1.29,3.03,0.11,2.6,0.26,5.6,0.21,0.36,0.88,0.809,0.806,0.6,0.825,0.923,0.719,0.91,0.907,0.76,0.742,0.785,0.877,0.799,0.927,0.923,1,0.897,0.908,0.773,0,5,0,0,6.6,5,0,0,5,0,0,5,0,5,0,6.6,10,5,9,0},
{0.42,1.34,0.9,0.93,1.38,0.73,0.55,0.87,0.9,0.78,0.7,3.21,0.41,8.08,2.25,5.94,0.32,0.42,0.949,0.658,0.641,0.404,0.892,0.814,0.576,0.845,0.823,0.886,0.875,0.584,0.981,0.937,0.941,0.792,0.897,1,0.983,0.888,0,6.6,0,0,5,5,0,0,5,0,5,0,6.6,9,0,5,5,10,5,0},
{0.63,1.09,0.82,0.89,1.3,1.2,0.98,0.73,1,0.73,0.77,2.94,0.3,6.47,0.96,5.66,0.25,0.41,0.959,0.687,0.63,0.403,0.845,0.809,0.594,0.825,0.817,0.899,0.893,0.614,0.974,0.938,0.94,0.789,0.908,0.983,1,0.904,0,5,0,0,5,5,0,0,6.6,0,0,0,0,9,0,6.6,9,5,10,0},
{0.26,1.06,0.66,0.99,1.7,0.68,0.46,0.61,0.87,0.61,0.77,3.67,0.14,3,1.22,6.02,0.27,0.49,0.968,0.417,0.31,0.092,0.737,0.556,0.335,0.576,0.569,0.997,0.997,0.347,0.932,0.984,0.74,0.542,0.773,0.888,0.904,1,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,10}};

//------ from the following PNAS paper related to AAindex ------//
/*
Solving the protein sequence metric problem
*/
double PNAS_AmiModel[20][5]={
{-0.591,-1.302,-0.733,1.570,-0.146},  //A
{1.538,-0.055,1.502,0.440,2.897},     //R
{0.945,0.828,1.299,-0.169,0.933},     //N
{1.050,0.302,-3.656,-0.259,-3.242},   //D
{-1.343,0.465,-0.862,-1.020,-0.255},  //C
{0.931,-0.179,-3.005,-0.503,-1.853},  //Q
{1.357,-1.453,1.477,0.113,-0.837},    //E
{-0.384,1.652,1.330,1.045,2.064},     //G
{0.336,-0.417,-1.673,-1.474,-0.078},  //H
{-1.239,-0.547,2.131,0.393,0.816},    //I
{-1.019,-0.987,-1.505,1.266,-0.912},  //L
{1.831,-0.561,0.533,-0.277,1.648},    //K
{-0.663,-1.524,2.219,-1.005,1.212},   //M
{-1.006,-0.590,1.891,-0.397,0.412},   //F
{0.189,2.081,-1.628,0.421,-1.392},    //P
{-0.228,1.399,-4.760,0.670,-2.647},   //S
{-0.032,0.326,2.213,0.908,1.313},     //T
{-0.595,0.009,0.672,-2.128,-0.184},   //W
{0.260,0.830,3.097,-0.838,1.512},     //Y
{-1.337,-0.279,-0.544,1.242,-1.262},  //V
};

//-------- process_part --------//
//                 A   B   C   D   E   F   G  H    I  J  K  L    M   N   O   P   Q   R    S    T  U   V   W   X   Y   Z
int AA4SUB[26] = { 1, -1,  5,  4,  7, 14,  8, 9,  10,-1, 12,11,  13, 3, -1, 15,  6,  2,  16,  17,-1, 20, 18, -1, 19, -1};


//------ load SS8 prob file -------//
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

//------- load PSIPRED and SolvPRED ---------//
//-> PSIPRED
/*
# PSIPRED VFORMAT (PSIPRED V3.5)

   1 S C   0.998  0.000  0.003
   2 K C   0.869  0.021  0.056
   3 E C   0.632  0.020  0.232
   4 L E   0.338  0.013  0.569
   5 K E   0.063  0.010  0.915
....
*/
//-> SolvPRED
/*
   1 S  0.925
   2 K  0.619
   3 E  0.512
   4 L  0.123
   5 K  0.341
   6 V  0.092
   7 L  0.065
   8 V  0.032
....
*/
int Load_SolvPRED(string &input_file, vector <double> &prob_out)
{
	//start
	ifstream fin;
	string buf,temp;
	//read
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",input_file.c_str());
		exit(-1);
	}
	//load
	int count=0;
	prob_out.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		for(int i=0;i<2;i++)www>>temp;  //skip
		double tmp_prob;
		www>>tmp_prob;
		prob_out.push_back(tmp_prob);
		count++;
	}
	//return
	return count;
}


//------ load LABEL_file for predicted contact number ------//
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
		int lab=c-'0';
		lab_number.push_back(lab);
	}
	//return
	return (int)buf.length();
}

//----- proc label -----//
void Proc_Label(vector <int> &lab_number)
{
	int i;
	int size=(int)lab_number.size();
	for(i=0;i<size;i++)
	{
		if(lab_number[i]==1)lab_number[i]=2;
		else break;
	}
	for(i=size-1;i>=0;i--)
	{
		if(lab_number[i]==1)lab_number[i]=2;
		else break;
	}
}

//--------- read FASTA sequence --------//
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

//-------------- for prediction -------------//
//given one TGT file, output it's feature as well as the label 
void Feature_Make(string &seq_file,string &ss3_file,string &acc_file,string &label_file)
{
	//-- read tpl ---//
//	string tgt_name,tgt_root;
//	getBaseName(tgt_file,tgt_name,'/','.');
//	getRootName(tgt_file,tgt_root,'/');
//	SEQUENCE *s=new SEQUENCE(tgt_name,tgt_root,0,1);
	string sequence;
	int length=Read_FASTA_SEQRES(seq_file,sequence);
	//-- load cn_file --//
	vector <int> cn_number;
	int cn_len=Load_LAB_File(label_file,cn_number);
	//-- check --//
	if(cn_len!=length)
	{
//		fprintf(stderr,"CN_len %d not equal to tgt_len %d \n",cn_len,s->length);
//		delete s;
//		exit(-1);
		for(int i=0;i<length;i++)cn_number.push_back(-1);
		cn_len=length;
	}
	//-- load other files --//
	vector < vector <double> > ss3_prob_out;
	int ss3_len=Load_Prob_File(ss3_file,3,ss3_prob_out);
	vector <double> acc_prob_out;
	int acc_len=Load_SolvPRED(acc_file,acc_prob_out);
	//-- check --//
	if(cn_len!=ss3_len || cn_len!=acc_len)
	{
		fprintf(stderr,"CN_len %d not equal to ss3_len %d or acc_len %d \n",cn_len,ss3_len,acc_len);
//		delete s;
		exit(-1);
	}

	//==== generate feature =====//
	vector <string> output;
	output.clear();
	for(int k=0;k<length;k++)
	{
		//----- feature -----//
		vector <double> features;
		features.clear();
		//-- profile realted --//
		//emission score
/*
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
*/
		//-- amino acid realted --//
		//zy ami model
		int pos=AA4SUB[sequence[k]-'A']-1;
		vector <int> aaind(20,0);
		if(pos<0 || pos>=20)
		{
			for(int i=0;i<58;i++)
			{
				features.push_back(0);
			}
			for(int i=0;i<5;i++)
			{
				features.push_back(0);
			}
			for(int i=0;i<20;i++)
			{
				features.push_back(aaind[i]);
			}
		}
		else
		{
			for(int i=0;i<58;i++)
			{
				double zy_model=ACC_ZY_AmiModel[pos][i];
				features.push_back(zy_model);
			}
			for(int i=0;i<5;i++)
			{
				double pnas_value=PNAS_AmiModel[pos][i];
				features.push_back(pnas_value);
			}
			aaind[pos]=1;
			for(int i=0;i<20;i++)
			{
				features.push_back(aaind[i]);
			}
		}

		//-- secondary structure ---//
		for(int i=0;i<3;i++)
		{
			double prob_reso = 1.0*ss3_prob_out[k][i];
			features.push_back(prob_reso);
		}
		//-- solvent accessibility ---//
//		for(int i=0;i<3;i++)
		{
			double prob_reso = 1.0*acc_prob_out[k];
			features.push_back(prob_reso);
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

	//-> printf
	printf("%d\n",length);
	for(int k=0;k<length;k++)printf("%s\n",output[k].c_str());
	Proc_Label(cn_number);
	for(int k=0;k<length;k++)printf("%d\n",cn_number[k]);
	
	//delete
//	delete s;
}



//------------ main -------------//
int main(int argc, char** argv)
{
	//------- Disorder Feature Process -----//
	{
		if(argc<5)
		{
			fprintf(stderr,"Diso_Feature_Make_amionly <seq_file> <ss3_file> <acc_file> <lab_file> \n");
			exit(-1);
		}
		string seq_file=argv[1];
		string ss8_file=argv[2];
		string acc_file=argv[3];
		string lab_file=argv[4];
		//process
		Feature_Make(seq_file,ss8_file,acc_file,lab_file);
		//exit
		exit(0);
	}
}


