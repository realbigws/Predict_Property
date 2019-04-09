#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;


//========= generate simp Summary File =========//
//-> required the following files,
//[1] 28079.diso     -> for disorder
//[2] 28079.ss3      -> for SSE percentage
//[3] 28079.acc      -> for ACC percentage


//-> binary classification file
int Load_BinClass_File(string &infile,string &amiseq,int &class_len,double thres)
{
	ifstream fin;
	string buf,temp,ami;
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found !!\n",infile.c_str());
		exit(-1);
	}
	//skip
	amiseq="";
	class_len=0;
	int totlen=0;
	double class_val;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file %s not found !!\n",infile.c_str());
			exit(-1);
		}
		if(buf=="")continue;
		if(buf[0]!='#')goto start;
	}
	//process
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
start:
		istringstream www(buf);
		www>>temp>>ami>>temp>>class_val;
		if(class_val>thres)class_len++;
		totlen++;
		amiseq+=ami;
	}
	return totlen;
}

//-> SS3 file
int Load_SS3_File(string &infile,string &amiseq,int &h_len,int &e_len,int &c_len )
{
	ifstream fin;
	string buf,temp,ami;
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found !!\n",infile.c_str());
		exit(-1);
	}
	//getline
	amiseq="";
	h_len=0;
	e_len=0;
	c_len=0;
	int totlen=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file %s not found !!\n",infile.c_str());
			exit(-1);
		}
		if(buf=="")continue;
		if(buf[0]!='#')goto start;
	}
	//process
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
start:
		istringstream www(buf);
		www>>temp>>ami>>temp;
		if(temp=="H")h_len++;
		if(temp=="E")e_len++;
		if(temp=="C")c_len++;
		totlen++;
		amiseq+=ami;
	}
	return totlen;
}

//-> ACC file
int Load_ACC_File(string &infile,string &amiseq,int &b_len,int &m_len,int &e_len )
{
	ifstream fin;
	string buf,temp,ami;
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found !!\n",infile.c_str());
		exit(-1);
	}
	//getline
	amiseq="";
	b_len=0;
	m_len=0;
	e_len=0;
	int totlen=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file %s not found !!\n",infile.c_str());
			exit(-1);
		}
		if(buf=="")continue;
		if(buf[0]!='#')goto start;
	}
	//process
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
start:
		istringstream www(buf);
		www>>temp>>ami>>temp;
		if(temp=="B")b_len++;
		if(temp=="M")m_len++;
		if(temp=="E")e_len++;
		totlen++;
		amiseq+=ami;
	}
	return totlen;
}


//--------- main --------//
int main(int argc,char **argv)
{
	//---- Generate Summary ----//
	{
		if(argc<8)
		{
			fprintf(stderr,"Generate_Summary_Simp <diso_file> <topo_file> <ss3_file> <acc_file> \n");
			fprintf(stderr,"                      <diso_thres> <topo_thres> <outfile> \n");
			exit(-1);
		}
		string diso=argv[1];
		string topo=argv[2];
		string ss3=argv[3];
		string acc=argv[4];
		double diso_thres=atof(argv[5]);   //default 0.5
		double topo_thres=atof(argv[6]);   //default 0.5
		string out_file=argv[7];
		//process
		int totlen2,totlen3,totlen4,totlen5;
		string amiseq2,amiseq3,amiseq4,amiseq5;
		int diso_len,topo_len;
		int h_len,e_len,c_len,B_len,M_len,E_len;
		totlen2=Load_BinClass_File(diso,amiseq2,diso_len,diso_thres);
		totlen3=Load_SS3_File(ss3,amiseq3,h_len,e_len,c_len);
		totlen4=Load_ACC_File(acc,amiseq4,B_len,M_len,E_len);
		totlen5=Load_BinClass_File(topo,amiseq5,topo_len,topo_thres);
		//check
		if(totlen2!=totlen3 || totlen2!=totlen4 || totlen2!=totlen5)
		{
			fprintf(stderr,"totlen not equal !! [%d %d %d %d] \n",
				totlen2,totlen3,totlen4,totlen5);
			exit(-1);
		}
		if(amiseq2!=amiseq3 || amiseq2!=amiseq4 || amiseq2!=amiseq5)
		{
			fprintf(stderr,"amiseq not equal !!\n%s\n%s\n%s\n%s\n",
				amiseq2.c_str(),amiseq3.c_str(),amiseq4.c_str(),amiseq5.c_str());
			exit(-1);
		}
		//output
		int totlen1=totlen2;
		FILE *fp=fopen(out_file.c_str(),"wb");
		fprintf(fp,"%d \n%d \n%d \n%d \n%d \n%d \n%d \n%d \n%d \n%d \n%d \n",
			totlen2,
			diso_len,(int)(100.0*diso_len/totlen1),
			(int)(100.0*h_len/totlen1),(int)(100.0*e_len/totlen1),(int)(100.0*c_len/totlen1),
			(int)(100.0*E_len/totlen1),(int)(100.0*M_len/totlen1),(int)(100.0*B_len/totlen1),
			topo_len,(int)(100.0*topo_len/totlen1));
		fclose(fp);
		exit(0);
	}
}

