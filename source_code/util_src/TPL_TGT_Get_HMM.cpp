#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

//----------- get HMM ----------//
//-> TPL format
/*
Version 1.0
Template Name  = 1pazA
Chain ID = A.1-123
Length   = 123
SEQRES sequence = ENIEVHMLNKGAEGAMVFEPAYIKANPGDTVTFIPVDKGHNVESIKDMIPEGAEKFKSKINENYVLTVTQPGAYLVKCTPHYAMGMIALIAVGDSPANLDQIVSAKKPKIVQERLEKVIASAK
DSSP   sequence = ENIEVHMLNKGAEGAMVFEPAYIKANPGDTVTFIPVDKGHNVESIKDMIPEGAEKFKSKINENYVLTVTQPGAYLVKCTPHYAMGMIALIAVGDSPANLDQIVSAKKPKIVQERLEKVIA---
NEFF = 6.4
Date = 2011-10-29  0:0:30
...
*/

//-> TGT format
/*
Version 1.1 -> user command = ./TPL_To_TGT2 -i ../TEMPLATE/1pazA.tpl
Sequence Name  = 1pazA
Length   = 123
Sequence = ENIEVHMLNKGAEGAMVFEPAYIKANPGDTVTFIPVDKGHNVESIKDMIPEGAEKFKSKINENYVLTVTQPGAYLVKCTPHYAMGMIALIAVGDSPANLDQIVSAKKPKIVQERLEKVIASAK
SSEseq   = LLEEEEEEELLLLLLEEEELLEEEELLLLEEEEEELLLLLLEEELLLLLLLLLLLLLLLLLLEEEEEELLLLEEEEELLLLLLLLLEEEEEELLLLLLHHHHLLLLLLHHHHHHHHHHHHHLL
SSEconf  = 979999777599986789869689979997999959999759996999988865556789965999976895699977555678869999987999996665578999578999999999859
ACCseq   = EEBEBMBMEEEEEMMBBBMMEEBEBMEEMEBMBMBMEEEBBBMBMEEEBEEEMEEMEMEEEEEMEBMBEEEEEBMBMBMBMMEMEBMBBBBBEEEEEEMEEBMEEEMEEMBMEMBEEBMEEME
ACCconf  = 974486865756844858446576658844758465695556454885446558757585866564478785576868435555574858464767665866565956966486679758849
NEFF = 6.9
EVD =
Date = 2014-4-30  8:8:35
...
*/
int TPL_TGT_Get_HMM(string &infile,FILE *fp)
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
	int length;
	string neff;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"%s format bad!\n",infile.c_str());
			exit(-1);
		}
		len=(int)buf.length();
		if(len<6)continue;
		temp=buf.substr(0,6);
		if(temp=="Length")
		{
			istringstream www(buf);
			www>>temp>>temp>>temp;
			length=atoi(temp.c_str());
		}
		if(temp=="NEFF =")
		{
			istringstream www(buf);
			www>>temp>>temp>>neff;
		}
		if(len<30)continue;
		temp=buf.substr(0,30);
		if(temp=="//////////// Original HHM file")break;
	}
	if(!getline(fin,buf,'\n'))
	{
		fprintf(stderr,"%s format bad!\n",infile.c_str());
		exit(-1);
	}
	//init
	fprintf(fp,"NEFF  %s \n",neff.c_str());
	fprintf(fp,"LENG  %d \n",length);
	//HMM
	int i;
	int count=0;
	fprintf(fp,"#\n");
	for(i=0;i<4;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"%s format bad!\n",infile.c_str());
			exit(-1);
		}
		fprintf(fp,"%s\n",buf.c_str());
	}
	for(i=0;i<3*length;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"%s format bad!\n",infile.c_str());
			exit(-1);
		}
		fprintf(fp,"%s\n",buf.c_str());
	}
	fprintf(fp,"//\n");
	return count;
}

//------ main ------//
int main(int argc,char **argv)
{
	//---- TPL_TGT_Get_HMM ----//
	{
		if(argc<3)
		{
			fprintf(stderr,"TPL_TGT_Get_HMM <tpl/tgt_file> <outfile>\n");
			exit(-1);
		}
		string file=argv[1];
		string outfile=argv[2];
		FILE *fp=fopen(outfile.c_str(),"wb");
		TPL_TGT_Get_HMM(file,fp);
		fclose(fp);
		exit(0);
	}
}
