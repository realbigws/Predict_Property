#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <vector>
#include <string>
#include "getopt.h"
#include "seq.h"
#include "DeepCNF.h"
using namespace std;

//--- model data related ----# start
#include "DeepCNF_model_1.h"
#include "DeepCNF_model_2.h"
#include "DeepCNF_model_3.h"
#include "DeepCNF_model_4.h"
#include "DeepCNF_model_5.h"
#include "DeepCNF_model_6.h"
#include "DeepCNF_model_7.h"
#include "DeepCNF_model_8.h"
#include "DeepCNF_model_con.h"
//--- model data related ----# end

//====== determine model ========//
void Determine_DeepCNF_Model(double * &model_weight, int &feat_num,int feat_lab)
{
	switch (feat_lab) 
	{
		case 1:
			feat_num = Feature_Model_1_size;
			model_weight = (double *)Feature_Model_1;
			break;
		case 2:
			feat_num = Feature_Model_2_size;
			model_weight = (double *)Feature_Model_2;
			break;
		case 3:
			feat_num = Feature_Model_3_size;
			model_weight = (double *)Feature_Model_3;
			break;
		case 4:
			feat_num = Feature_Model_4_size;
			model_weight = (double *)Feature_Model_4;
			break;
		case 5:
			feat_num = Feature_Model_5_size;
			model_weight = (double *)Feature_Model_5;
			break;
		case 6:
			feat_num = Feature_Model_6_size;
			model_weight = (double *)Feature_Model_6;
			break;
		case 7:
			feat_num = Feature_Model_7_size;
			model_weight = (double *)Feature_Model_7;
			break;
		case 8:
			feat_num = Feature_Model_8_size;
			model_weight = (double *)Feature_Model_8;
			break;
		case 0:
			feat_num = Feature_Model_Con_size;
			model_weight = (double *)Feature_Model_Con;
			break;
		default:
			exit(-1);
	}
};


//=========================== TGT file part =========================//

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
int AA4SUB[26] = { 1, -1,  5,  4,  7, 14,  8, 9,  10,-1, 12,11,  13, 3, -1, 15,  6,  2,  16,  17,-1, 20, 18, -1, 19, -1};

//-------------- for probability make -------------//
//given 1st layer feature, output 8 label probability
void Prob_Make(DeepCNF_Model *m_pModel, vector <string> &input,vector <string> &out_str)
{
	//-> construct a new sequence
	DeepCNF_Seq *seq = new DeepCNF_Seq(input.size(), m_pModel);
	seq->Read_Init_Features(input);
	//-> predict ss8 label
	vector <vector <double> > output;
	seq->MAP_Probability(output);
	//-> output
	out_str.clear();
	for(int k=0;k<(int)output.size();k++)
	{
		//----- feature -----//
		vector <double> features=output[k];
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
		out_str.push_back(wsbuf);
	}
	//-> delete
	delete seq;
}

//-------------- for primary feature make -------------//
//given one TGT file, output it's primary feature 
int Feature_Make(string &tgt_file,vector <string> &output,string &ami_seq)
{
	//-- read tgt ---//
	string tgt_name,tgt_root;
	getBaseName(tgt_file,tgt_name,'/','.');
	getRootName(tgt_file,tgt_root,'/');
	SEQUENCE *s=new SEQUENCE(tgt_name,tgt_root,0,1);

	//==== generate feature =====//
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
		for(int i=0;i<3;i++)
		{
			double psipred_reso = 1.0*s->SS2[k][i];
			features.push_back(psipred_reso);
		}
		//-- solvant accessibility --//-> for SAS only
		for(int i=0;i<3;i++)
		{
			double acc_reso = 1.0*s->acc_our_10_42[k][i];
			features.push_back(acc_reso);
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

	//delete
	ami_seq=s->sequence;
	int length=s->length;
	delete s;
	//return
	return length;
}

//------- for secondary feature make ------//
//given one TGT file, output it's secondary feature 
void Feature_Make_II(string &tgt_file,vector <string> &feat_out_ii,string &ami_seq)
{
	//--- generate primary feature ----//
	vector <string> feat_out;
	int length=Feature_Make(tgt_file,feat_out,ami_seq);
	//--- generate secondary feature ---//
	int model_layer=5;
	string window_str = "5,5,5,5,5";
	string node_str = "100,100,100,100,100";
	int state_num = 3;
	int local_num = 66;
	DeepCNF_Model *m_pModel = new DeepCNF_Model(model_layer,window_str,node_str,state_num,local_num,0);
	vector <vector <string> > prob_out_total;
	prob_out_total.clear();
	//for each DeepCNF model
	for(int s=1;s<=8;s++)
	{
		//-> load model
		double * model_weight;
		int feat_num;
		Determine_DeepCNF_Model(model_weight, feat_num,s);
		m_pModel->cnf_model_load(model_weight,feat_num);
		//-> calculate prob
		vector <string> prob_out;
		Prob_Make(m_pModel,feat_out,prob_out);
		prob_out_total.push_back(prob_out);
	}
	delete m_pModel;
	//output
	feat_out_ii.clear();
	for(int k=0;k<length;k++)
	{
		string tmp_str="";
		for(int s=0;s<(int)prob_out_total.size();s++)tmp_str=tmp_str+prob_out_total[s][k]+" ";
		feat_out_ii.push_back(tmp_str);
	}
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
//--- SAS to Char ---//
char SAS_To_Char(int c)
{
	switch(c)
	{
		case 0: return 'B';
		case 1: return 'M';
		case 2: return 'E';
		default: return 'M';
	}
}


//------------ usage -------------//
void Usage() {
	cerr << "DeepCNF_SAS v1.00 [Apr-08-2016] \n\n";
	cerr << "./DeepCNF_SAS -t tgt_file [-m model] \n";
	cerr << "Options:\n\n";
	cerr << "-t tgt_file : input tgt file. \n";
	cerr << "-m model :    1-8 for primary model, 0 for consensus. \n\n";
}


//------------ main -------------//
int main(int argc, char** argv)
{
	//-- help --//
	if (argc < 2)
	{
		Usage();
		exit(0);
	}

	//---- init parameter ----//
	string input_tgt_file = "";
	int model_num = 0;  //-> consensus model

	//command-line arguments process
	extern char* optarg;
	char c = 0;
	while ((c = getopt(argc, argv, "t:m:")) != EOF) 
	{
		switch (c) 
		{
		//-> input tgt file
		case 't':
			input_tgt_file = optarg;
			break;
		//-> input model
		case 'm':
			model_num = atoi(optarg);
			break;
		//-> default
		default:
			Usage();
			exit(-1);
		}
	}

	//----- check parameter -----//
	//-> check input file
	if(input_tgt_file=="")
	{
		fprintf(stderr,"input tgt %s is NULL\n",input_tgt_file.c_str());
		exit(-1);
	}
	//-> check input model
	if(model_num<0 || model_num>8)
	{
		fprintf(stderr,"model_num %d should be 0 to 8 \n",model_num);
		exit(-1);
	}

	//===================== initilize weights ================//start
	//----------- init model to get the dimension of weights ------------//
	DeepCNF_Model *m_pModel;
	vector <string> feat_out;
	string ami_seq;
	if(model_num==0)
	{
		//-> consensus model
		int model_layer=5;
		string window_str = "5,5,5,5,5";
		string node_str = "20,20,20,20,20";
		int state_num = 3;
		int local_num = 24; //-> consensus mode, 8*3=24 input features
		//-> load model
		double * model_weight;
		int feat_num;
		Determine_DeepCNF_Model(model_weight, feat_num, model_num);
		//-> generate feature
		m_pModel = new DeepCNF_Model(model_layer,window_str,node_str,state_num,local_num,0);
		m_pModel->cnf_model_load(model_weight,feat_num);
		Feature_Make_II(input_tgt_file,feat_out,ami_seq);
	}
	else
	{
		//-> primary model
		int model_layer = 5;
		string window_str = "5,5,5,5,5";
		string node_str = "100,100,100,100,100";
		int state_num = 3;
		int local_num = 66;
		//-> load model
		double * model_weight;
		int feat_num;
		Determine_DeepCNF_Model(model_weight, feat_num, model_num);
		//-> generate feature
		m_pModel = new DeepCNF_Model(model_layer,window_str,node_str,state_num,local_num,0);
		m_pModel->cnf_model_load(model_weight,feat_num);
		Feature_Make(input_tgt_file,feat_out,ami_seq);
	}

	//-> construct a new sequence
	DeepCNF_Seq *seq = new DeepCNF_Seq(feat_out.size(), m_pModel);
	seq->Read_Init_Features(feat_out);
	//-> predict SAS label
	vector <vector <double> > output;
	seq->MAP_Probability(output);
	//-> output header
	printf("#DeepConCNF_SAS: three-state solvent accessibility prediction results \n");
	printf("#probabilities are in the order of B (Bury, pACC: 0-10), M (Medium, pACC: 11-40) and E (Exposed, pACC: 41-100), \n");
	printf("#where pACC is the relative solvent accessibility value calculated by DSSP \n\n");
	//-> output content
//   1 E L   0.000 0.000 0.000 0.000 0.000 0.000 0.000 1.000
//   2 N L   0.000 0.000 0.000 0.379 0.004 0.000 0.000 0.617
	for(int i=0;i<(int)output.size();i++)
	{
		//return maximal value
		double maxval;
		int label=Return_Label(output[i],maxval);
		char lab_c=SAS_To_Char(label);
		//output
		printf("%4d %c %c   ",i+1,ami_seq[i],lab_c);
		int state_num = 3;
		for(int j=0;j<state_num;j++)printf("%5.3f ",output[i][j]);
		printf("\n");
	}

	//delete
	delete seq;
	delete m_pModel;

	//exit
	exit(0);
}
