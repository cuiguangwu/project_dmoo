#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <vector>
#include<map>
#include <cassert>
#include <iomanip>
#include <numeric>
#include<algorithm>

#include "config.h"
using namespace std;

class evaluate
{
public:
	vector < vector<double>>  pf;
	vector < vector<double>>  p;
	vector <double>  in_avg; //每个指标的平均值
	vector <double>  igd_avg;
	
	double dist_vector(vector <double> &vec1, vector <double> &vec2);
	void loadpfront(char *filename, vector< vector<double>> &pf,int nobj);
	void newPf();
	
	double STD(vector<double> v);
	double Avg(vector<double> v);	
	double indicator_GD();
	double indicator_IGD();
	double indicator_SP(vector < vector<double>>  objectf);
	double indicator_MS(vector<vector<double>> v1, vector<vector<double>> v2,int k);
	double indicator_hvd(vector<vector<double>> vd, int gen,string instances,int nt);
	void   indicator_AVG(vector < vector<double>>  f, vector<double> &avg);
	double main_hv(vector<vector<double>> vd);
	void frontcheck(vector<vector<double>>&vd);
	double pfhv(string instances,int gen,int nt);
	int dominates(int index,vector<double> v1,vector<double> v2);
	double hv(vector<vector<double>>&vd);
	double inclhv(vector<vector<double>> &vd);
	double inclhv2(vector<vector<double>> &vd);
	double inclhv3(vector<vector<double>> &vd);
	double distance;                   // generational distance from PF to solutions found
	void getPOF(char * problem, int nt,int gen,vector< vector<double>> &pf);//得到真实的POF
	void loadIGD(char *filename,vector<double> &pf, int nobj);
	void nodominate(vector< vector<double>> &pf); //remove donamited individuals
	int Dominate(vector<double>index1,vector<double>index2);
	double area_sum(vector<vector<double>> vd,int i); // hv切面的面积和


	void GetAvgIGD(char * problem,int run,int gener);
	void Statistics(char * problem,int metrics,int run);//统计数据
	
	void loads(char *filename, vector< vector<double>> &pf, int nobj);
	void readToScience(vector< vector<double>> &pf,char * name);
	void JY1();
	void JY2();//取点
	void JY5();
	void dmop1();
	double JY_area(double x,double A,double W);
};

void evaluate ::indicator_AVG(vector < vector<double>>  f,vector<double> &avg){
	int size =f.size();
	int obj=f[0].size();
	int i,j;
	double sum;
	avg.clear();
	for(i=0;i<obj;i++){
		sum=0;
		for(j=0;j<size;j++){
			sum+=f[j][i];			
		}	
		sum/=size;
		avg.push_back(sum);
	}
}


double evaluate ::Avg(vector<double> v)
{
	double sum = accumulate(v.begin(), v.end(), 0.0);
	if (v.size()) return sum / v.size();
	else return 0;
}

double evaluate ::STD(vector<double> v)
{
	double avgv = Avg(v);
	int size = v.size();
	double stdv = 0;
	for (int i = 0; i < size; i++)
		stdv += pow(v[i] - avgv, 2.0);

	if (size>1) return sqrt(stdv / (size - 1));
	else return 0;
}

double evaluate ::indicator_SP(vector < vector<double>>  objectf){
	int size = objectf.size();
	vector<double> dist(size, 1.0e+30);
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
		{
			if (j != i)
			{
				double tmp = dist_vector(objectf[i], objectf[j]);
				if (tmp < dist[i]) dist[i] = tmp;
			}
		}
	return STD(dist);
}


double evaluate::indicator_MS(vector<vector<double>> v1, vector<vector<double>> v2,int nobj){
	vector<double> v1_max(nobj, -1.0e+30), v1_min(nobj, 1.0e+30);
	vector<double> v2_max(nobj, -1.0e+30), v2_min(nobj, 1.0e+30);
	for (int n = 0; n < nobj; n++)
	{
		for (int i = 0; i < v1.size(); i++)
		{
			if (v1_max[n] < v1[i][n]) v1_max[n] = v1[i][n];
			if (v1_min[n] > v1[i][n]) v1_min[n] = v1[i][n];
		}
		for (int i = 0; i < v2.size(); i++)
		{
			if (v2_max[n] < v2[i][n]) v2_max[n] = v2[i][n];
			if (v2_min[n] > v2[i][n]) v2_min[n] = v2[i][n];
		}

	}
	double ms = 0;
	for (int n = 0; n < nobj; n++)
	{
		double lebesgue = 0;
		if (v2_min[n] < v1_max[n])
		{
			double minv = min(v1_max[n], v2_max[n]);
			double maxv = max(v1_min[n], v2_min[n]);
			lebesgue = minv - maxv;
		}
		ms += pow(lebesgue / (v1_max[n] - v1_min[n]), 2.0);
	}
	return sqrt(ms / nobj);
}

double evaluate::indicator_hvd(vector<vector<double>>vd, int gen,string instances,int nt){
	int size1 = vd.size();
	int size2 = vd[0].size();
	if(instances=="FDA4" || instances=="FDA5" ||instances=="DMOPD"){
		probreference.resize(3);
		probreference[0]=-1;probreference[1]=-1;probreference[2]=-1;
	}else{
		probreference.resize(2,-1);probreference[0]=-1;probreference[1]=-1;
	}
	
	for (int i = 0; i < size1; i++){
		for (int j = 0; j < size2; j++){		
			if((probreference[j]-0.5)<vd[i][j]){
				probreference[j]=vd[i][j]+0.5;
			}
		}
	}
	double pfv=pfhv(instances,gen,nt);

	double hvd = fabs(main_hv(vd) - pfv);

	return hvd;
}

double evaluate::main_hv(vector<vector<double>> vd)
{
	double hvs=0;
	double f1_d=0, f2_d=0, f3_d=0;
	double sum_s=0;
	double sum_hv=0;
	map<int,vector<double>> f1;
	map<int,vector<double>> f2;
	frontcheck(vd);
	
	int maxm=vd.size();
	int maxn=vd[0].size();
	for(int i=0;i<maxm;i++)  //第一维降序排列
	{
		for(int j=i+1;j<maxm;j++) 
		{
			if(vd[i][0]<vd[j][0]){
				std::swap(vd[i],vd[j]);					
			}
		}
	}	
	
	if(maxn==2)
	{
		hvs+=fabs((probreference[0]-vd[0][0])*(probreference[1]-vd[0][1]));
		for(int i=1;i<maxm;i++)
		{
			hvs+= fabs((probreference[1]-vd[i][1])*(vd[i-1][0]-vd[i][0]));
		}
	}
	else
	{
		double slice=fabs(probreference[0]-vd[0][0]);
		hvs += slice*area_sum(vd,0);
		for(int i=1;i<vd.size();i++){
			hvs += area_sum(vd,i)*fabs(vd[i][0]-vd[i-1][0]);
		}
	}
		
	return hvs;
}

double evaluate::area_sum(vector<vector<double>> vd,int number){
	int size=vd.size();
	vector<vector<double>> efficient;
	for(int i=number;i<size;i++){
		int flg=1;
		for(int j=i+1;j<size;j++){
			int dominate=dominates(1,vd[i],vd[j]);
			if(dominate==-1)
			{
				flg=0;			
				break;
			}
		}
		if(flg==1)efficient.push_back(vd[i]);
	}
	int maxm=efficient.size();
	for(int i=0;i<maxm;i++)  //第二维降序排列
	{
		for(int j=i+1;j<maxm;j++) 
		{
			if(efficient[i][1]<efficient[j][1]){
				std::swap(efficient[i],efficient[j]);					
			}
		}
	}
	double area=0;
	//std::cout<<probreference[1]<<efficient[0][1];
	area +=fabs((probreference[1]-efficient[0][1])*(probreference[1]-efficient[0][2]));
	for(int i=1;i<maxm;i++)
	{
		area+= fabs((probreference[2]-vd[i][2])*(efficient[i-1][1]-efficient[i][1]));
	}
	efficient.clear();
	return area;
}

double evaluate::hv(vector<vector<double>> &vd){
	int size1=vd.size();
	switch (size1)
	{
		case 0 :return 0;
		case 1: return inclhv(vd);
		case 2: return inclhv2(vd);
		case 3: return inclhv3(vd);
	}
	return 0;
}

double evaluate::inclhv3(vector<vector<double>> &vd)
{
	double vp = 1; double vq = 1; double vr = 1;
	double v01 = 1; double v02 = 1; double v12 = 1;
	double v012 = 1;
	int size1=vd[0].size();
	for(int i=0;i<size1;i++)
	{
		vp *= vd[0][i];
		vq *= vd[1][i];
		vr *= vd[2][i];
		if(vd[0][i]<vd[1][i])
		{
			if(vd[1][i]<vd[2][i]){
				v01 *= vd[0][i];
				v02 *= vd[0][i];
				v12 *= vd[1][i];
				v012 *=vd[0][i];
			}
			else
			{
				if(vd[0][i]<vd[2][i])
				{
					v01 *= vd[0][i];
					v02 *= vd[0][i];
					v12 *= vd[2][i];
					v012 *=vd[0][i];
				}else
				{
					v01 *= vd[0][i];
					v02 *= vd[2][i];
					v12 *= vd[2][i];
					v012 *= vd[2][i];
				}
			}
		}
		else
		{
			if(vd[0][i]<vd[2][i])
			{
				v01 *= vd[1][i];
				v02 *= vd[0][i];
				v12 *= vd[1][i];
				v012 *= vd[1][i];
			}else
			{
				if(vd[1][i]<vd[2][i])
				{
					v01 *= vd[1][i];
					v02 *= vd[2][i];
					v12 *= vd[1][i];
					v012 *= vd[1][i];
				}else
				{
					v01 *= vd[1][i];
					v02 *= vd[2][i];
					v12 *= vd[2][i];
					v012 *= vd[2][i];
				}
			}
		}
	}
	return vp + vq + vr - v01 - v02 - v12 + v012;
}

double evaluate::inclhv2(vector<vector<double>> &vd)
{
	
	int maxn=vd[0].size();
	double vp = 1; double vq = 1;
	double vpq = 1;
	for(int i=0;i<maxn;i++)
	{
		vp *= vd[0][i];
		vq *= vd[1][i];
		if(vd[0][i]<vd[1][i])
		{
			vpq *=vd[0][i];
		}
		else
		{
			vpq *=vd[1][i];
		}
	}
	return vp + vq - vpq;
}

double evaluate::inclhv(vector<vector<double>> &vd){
	double volume = 1;
	int size1=vd[0].size();
	for(int i=0;i<size1;i++){
		volume*=vd[0][i];
	}
	return  volume;
}

void evaluate::frontcheck(vector<vector<double>>&vd)
{
	vector<vector<double>> temp;
	int flg;
	for(int i=0;i<vd.size();i++){
		flg=1;
		for(int j=0;j<vd.size();j++){
			int dominate=dominates(0,vd[i],vd[j]);
			if(dominate==-1)
			{
				flg=0;			
				break;
			}
		}
		if(flg==1)temp.push_back(vd[i]);
	}
	vd.clear();
	vd=temp;
}

int evaluate::dominates(int index,vector<double> v1,vector<double> v2){
	int size1=v1.size();
	int better=0,worse=0;
	for(int i=index;i<size1;i++)
	{
		if(v1[i]<=v2[i]) better++;
		if(v1[i]>=v2[i]) worse++;
	}
	if(better==size1) return 1;
	else if(worse==size1) return -1;
	return 0;
}

double evaluate::pfhv(string instances,int gen,int nt)
{
	double t, chv = 1, phv = 0, thv=0, G, H;
	//int nt=10;	
	t = 1.0 / nt*(gen-1);

	if(instances == "FDA1")
	{		
		phv = 1.0 / 3;
		chv =probreference[0]*probreference[1];
		thv= chv - phv;		
	}else if(instances == "FDA2")
	{		
		H=0.75+0.7*sin(0.5*PI*t);
	
		G=H+10*(H+1)*(H+1);
		G=pow(G,-1);
		phv=1- 1.0 /(G+1);
		chv = probreference[0]*probreference[1];
		thv= chv - phv;		
	}else if(instances == "FDA3")
	{
		G = fabs(sin(0.5*PI*t));
		double g = G+1;
		double ft=pow(10,2*sin(0.5*PI*t));
		chv=probreference[0]*probreference[1];
		phv = (1+G)*(1-2.0/(sqrt(g)*(ft+2)));  //???表示怀疑
		thv= chv - phv;
		double gts=fabs(sin(0.5*PI*t))+1+0.5;		

	}else if(instances == "FDA4" || instances == "FDA4"|| instances == "DMOPD")
	{
		chv=probreference[0]*probreference[1]*probreference[2];
		phv = (1.0/8)*4.0*PI / 3;
		thv=chv-phv;		
	}else if(instances == "DMOP1")
	{
		H=1.25+0.75*sin(0.5*PI*t);
		phv= H / (H+1);
		chv = probreference[0]*probreference[1];
		thv = chv-phv;
		probreference.resize(2,1.5);
	}else if(instances == "DMOP2"){
		H=1.25+0.75*sin(0.5*PI*t);
		phv= H / (H+1);
		chv = probreference[0]*probreference[1];
		thv = chv-phv;
	}else if(instances == "DMOP3"){
		phv = 1.0 / 3;
		chv = probreference[0]*probreference[1];
		thv= chv - phv;
	}
	else if(instances == "DMOPA" ||instances == "DMOPB" ||instances == "DMOPC" ||instances == "DMOPE" ||instances == "DMOPF"){
		H=1.25+0.75*sin(PI*t);
		double f1=0;
		double f2=0;
		double s=0.001; 
		double dx=0;
		double area=0;
		double fx=1;
		while(s<=1){
			f1=pow(s,H);
			f2=pow((1-s),H);
			area += (f1-dx)*(fx+f2)/2;
			dx=f1;
			fx=f2;
			s+=0.001;
		}
		phv=area;
		chv = probreference[0]*probreference[1];
		thv = chv-phv;
	}else if(instances == "JY1"){
		phv=JY_area(1,0.05,6)-JY_area(0,0.05,6);
		chv = probreference[0]*probreference[1];
		thv = chv-phv;
	}else if(instances == "JY2"){
		double A=0.05;
		double W=floor(6*sin(0.5*PI*(t-1)));
		phv=JY_area(1,A,W)-JY_area(0,A,W);
		chv = probreference[0]*probreference[1];
		thv = chv-phv;
	}else if(instances == "JY3"){
		double A=0.05;
		double W=floor(6*sin(0.5*PI*(t-1)));
		phv=JY_area(1,A,W)-JY_area(0,A,W);
		chv = probreference[0]*probreference[1];
		thv = chv-phv;
	}
	else if(instances == "JY4"){  //非连续的图形
		phv=main_hv(pf); //去点求面积
		//chv = probreference[0]*probreference[1];
		thv = phv;
	}else if(instances == "JY5"){
		double A=0.3*sin(0.5*PI*(t-1));
		double W=1;
		phv=JY_area(1,A,W)-JY_area(0,A,W);
		chv = probreference[0]*probreference[1];
		thv = chv-phv;
	}else if(instances == "JY6"){
		double A=0.1;
		double W=3;
		phv=JY_area(1,A,W)-JY_area(0,A,W);
		chv = probreference[0]*probreference[1];
		thv = chv-phv;
	}else if(instances == "JY7"){
		double A=0.1;
		double W=3;
		phv=JY_area(1,A,W)-JY_area(0,A,W);
		chv = probreference[0]*probreference[1];
		thv = chv-phv;
	}
	else if(instances == "JY8"){
		double A=0.05;
		double W=6;
		phv=JY_area(1,A,W)-JY_area(0,A,W);
		chv = probreference[0]*probreference[1];
		thv = chv-phv;
	}else if(instances == "JY9"){
		double A=0.05;
		int sigma = (int) floor(t/5)%3; 
		double W=floor(6*pow(sin(0.5*PI*(t-1)),sigma));
		phv=JY_area(1,A,W)-JY_area(0,A,W);
		chv = probreference[0]*probreference[1];
		thv = chv-phv;
	}
	return thv;
}

double evaluate::JY_area(double x,double A,double W){
	if(W!=0){
		return x-0.5*x*x -A*cos(W*PI*x)/(W*PI)+A*sin(W*PI*x)-A*x*sin(W*PI*x)-A*cos(W*PI*x)/(W*PI)+A*A*pow(sin(W*PI*x),2)/2;
	}else
	{
		return x-0.5*x*x;
	}
}

double evaluate::indicator_GD(){
	distance=0.0;

	for(int i=0; i<p.size(); i++)
	{
		double min_d = 1.0e+10;
		for(int j=0; j<pf.size(); j++)
		{
			double d = dist_vector(p[i], pf[j]);
			if(d<min_d)  min_d = d;
		}
		distance+= min_d;
	}

	return distance/=p.size();
}

//IGD
double evaluate::indicator_IGD()
{
	distance=0.0;
	for(int i=0; i<pf.size(); i++)
	{
		double min_d = 1.0e+10;
		for(int j=0; j<p.size(); j++)
		{
			double d = dist_vector(pf[i], p[j]);
			if(d<min_d)  min_d = d;
		}
		distance+= min_d;
	}
	distance/=pf.size();
	return distance;
}

void evaluate::loadpfront(char *filename, vector< vector<double>> &pf, int nobj)
{
	std::fstream fin;
	int line=0;
	char str[100]=" ";
	fin.open(filename,std::ios::in);
	pf.clear();
	if(fin.is_open())
	{
		const char* pch2;
		char  a[20],b[20],c[20],d[20],e[20];
		std::string str;

		while(!fin.eof())
		{
			vector<double> data;
			std::getline(fin,str,'\n');
			pch2 = str.c_str();
			sscanf(pch2,"%s %s %s %s %s",a,b,c,d,e);

			data.push_back(atof(a));
			data.push_back(atof(b));
			if (nobj==3)
			   data.push_back(atof(c));
			if(nobj==4) 
			{
				data.push_back(atof(c));
				data.push_back(atof(d));
			}
			if(nobj==5){
				data.push_back(atof(c));
				data.push_back(atof(d));
				data.push_back(atof(e));
			}
			line++;
			pf.push_back(data);
		}
	} //end if

	else
		std::cout<<"failed to open "<<filename<<endl;
	fin.close();
}

void evaluate::loads(char *filename, vector< vector<double>> &pf, int nobj){
	std::fstream fin;
	int line=0;
	char str[100]=" ";
	fin.open(filename,std::ios::in);
	pf.clear();
	if(fin.is_open())
	{
		const char* pch2;
		char  a1[20],a2[20],a3[20],a4[20],a5[20],a6[20],a7[20],a8[20],a9[20],a10[20],a11[20],a12[20];
		std::string str;

		while(!fin.eof())
		{
			vector<double> data;
			std::getline(fin,str,'\n');
			pch2 = str.c_str();
			sscanf(pch2,"%s %s %s %s %s %s %s %s %s %s %s %s",a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12);

			data.push_back(atof(a1));
			data.push_back(atof(a2));
			data.push_back(atof(a3));
			data.push_back(atof(a4));
			data.push_back(atof(a5));
			data.push_back(atof(a6));
			data.push_back(atof(a7));
			data.push_back(atof(a8));
			data.push_back(atof(a9));
			data.push_back(atof(a10));
			data.push_back(atof(a11));
			data.push_back(atof(a12));
			line++;
			pf.push_back(data);
		}
	} //end if

	else
		std::cout<<"failed to open "<<filename<<endl;
	fin.close();
}
void evaluate::readToScience(vector< vector<double>> &pf,char * name){
	std::ofstream files;
	char fis[1023];
	sprintf(fis,"%s",name);
	int length = pf.size();
	length/=2;
	files.open(fis,ios::trunc);
	files<<setprecision(4)<<scientific;
	for(int i=0;i<length;i++){
		files<<"&"<<pf[2*i][0]<<"("<<pf[2*i+1][0]<<")$\\ddagger$  "<<"&"<<pf[2*i][1]<<"("<<pf[2*i+1][1]<<")$\\ddagger$  "<<"&"<<pf[2*i][2]<<"("<<pf[2*i+1][2]<<")$\\ddagger$  "<<"&\\textbf{"<<pf[2*i][3]<<"("<<pf[2*i+1][3]<<")}\\"<<"\\"<<"\n";
		files<<"&"<<pf[2*i][4]<<"("<<pf[2*i+1][4]<<")$\\ddagger$  "<<"&"<<pf[2*i][5]<<"("<<pf[2*i+1][5]<<")$\\ddagger$  "<<"&"<<pf[2*i][6]<<"("<<pf[2*i+1][6]<<")$\\ddagger$  "<<"&\\textbf{"<<pf[2*i][7]<<"("<<pf[2*i+1][7]<<")}\\"<<"\\"<<"\n";
		files<<"&"<<pf[2*i][8]<<"("<<pf[2*i+1][8]<<")$\\ddagger$  "<<"&"<<pf[2*i][9]<<"("<<pf[2*i+1][9]<<")$\\ddagger$  "<<"&"<<pf[2*i][10]<<"("<<pf[2*i+1][10]<<")$\\ddagger$  "<<"&\\textbf{"<<pf[2*i][11]<<"("<<pf[2*i+1][11]<<")}\\"<<"\\"<<"\n";
		
	}
}

void evaluate::loadIGD(char *filename,vector<double> &pf, int nobj)
{
	std::fstream fin;
	int line=0;
	char str[100]=" ";
	fin.open(filename,std::ios::in);
	pf.clear();
	if(fin.is_open())
	{
		const char* pch2;
		char  a[20],b[20],c[20],d[20],e[20];
		std::string str;

		while(!fin.eof())
		{	
			std::getline(fin,str,'\n');
			pch2 = str.c_str();
			sscanf(pch2,"%s %s %s %s %s",a,b,c,d,e);
			
			pf.push_back(atof(b));	
			
			line++;
			
		}
	} //end if

	else
		std::cout<<"failed to open "<<filename<<endl;
	fin.close();
}

//求平均IGD
void evaluate::GetAvgIGD(char * problem,int run,int predict){
	std::ofstream pf;
	char avgIGD[1024];
	char strTest[1024];
	vector<double> data(predict,0);

	sprintf(avgIGD,"evaluate/avgIGD/%s.dat",problem);  //平均IGD路径
	pf.open(avgIGD,ios::trunc);


	for(int ccc=0;ccc<run;ccc++){//样本均值	
		sprintf(strTest,"evaluate/data/%s_%d.dat",problem,ccc);
		vector <double>  igds;
		loadIGD(strTest,igds,4);
				
		for(int len=0;len<predict;len++){
			data[len]+=igds[len];
		}
	}
	for(int len=0;len<predict;len++){
		data[len] =data[len]/run;
		pf<<data[len]<<"      "<<"\n";
	}		
	pf.close();
}

void evaluate::Statistics(char * problem,int metrics,int run){
	std::ofstream pf;
	char files[1024];
	char readProblem[1024];
	sprintf(readProblem,"evaluate/avg/%s.dat",problem);
	vector<vector<double>> datas;
	loadpfront(readProblem,datas,metrics);
	datas.erase(datas.begin(),datas.begin()+run);
	pf<<setprecision(3)<<setiosflags(ios::scientific);
	for(int i=0;i<metrics;i++){
		if(i==0)
		sprintf(files,"evaluate/statistics/%s.dat","GD");
		else if(i==1)
			sprintf(files,"evaluate/statistics/%s.dat","IGD");
		else if(i==2)
			sprintf(files,"evaluate/statistics/%s.dat","SP");
		else if(i==3)
			sprintf(files,"evaluate/statistics/%s.dat","MS");
		else if(i==4)
			sprintf(files,"evaluate/statistics/%s.dat","HV");
		pf.open(files,ios::app);
		if(problem =="FDA1" || problem =="FDA2" || problem =="FDA3" || problem =="FDA4" || problem =="FDA5")
			pf<<problem<<"    "<<datas[0][i]<<"     "<<datas[1][i]<<" "<<"\n";
		else if(problem =="DMOP1"||problem =="DMOP2" ||problem =="DMOP3"){
			pf<<problem<<"   "<<datas[0][i]<<"     "<<datas[1][i]<<" "<<"\n";
		}else{
			pf<<problem<<"     "<<datas[0][i]<<"     "<<datas[1][i]<<" "<<"\n";
		}
		pf.close();
	}	
}

//取点
void evaluate::getPOF(char * problem, int nt,int gen,vector< vector<double>> &pf){
	char  filename[1024];
	double x = 0 ; // 精度
	double t = 0;
	t = (double)(gen-1)/nt;
	while(x<=1){
		vector<double> data;
		if(problem == "FDA1")
		{
			data.push_back(x);
			data.push_back(1-sqrt(x));
		}else if(problem == "FDA2"){
			double ht=0.75+0.7*sin(0.5*PI*t);
			data.push_back(x);
			data.push_back(1-pow(x,pow(ht+10*(1+ht)*(1+ht),-1)));
		}else if(problem == "FDA3"){
			double ft = pow(10,2*sin(0.5*PI*t));
			double gt=fabs(sin(0.5*PI*t));
			data.push_back(x);
			data.push_back((1-sqrt(x/(1+gt)))*(1+gt));
		}else if(problem == "FDA4" || problem == "FDA5" || problem =="DMOPF"){ //取 1000个点
			int i,j;
			for(i=0;i<100;i++){
				for(j=0;j<100;j++){
					double theta=(double)i/99;
					double varphi=(double)j/99;
					data.push_back(sin(theta*PI/2)*cos(varphi*PI/2));
					data.push_back(sin(theta*PI/2)*sin(varphi*PI/2));
					data.push_back(cos(theta*PI/2));
					pf.push_back(data);
					data.clear();
					if(i==0)break;
				}
			}
			return;
		}else if(problem == "DMOP1"){
			double h=1.25+0.75*sin(0.5*PI*t);
			data.push_back(x);
			data.push_back(1-pow(x,h));
		}else if(problem == "DMOP2")
		{
			double h=1.25+0.75*sin(0.5*PI*t);
			data.push_back(x);
			data.push_back(1-pow(x,h));
		}else if(problem == "DMOP3"){
			data.push_back(x);
			data.push_back(1-sqrt(x));
		}
		
		else if(problem == "DMOPA"){
			double H = 1.25+0.75*sin(PI*t);
			data.push_back(pow(x,H));
			data.push_back(pow((1-x),H));
			//pfile<<pow(x,H)<<"    "<<setprecision(5)<<pow((1-x),H)<<"\n";
		}else if(problem == "DMOPB")
		{
			double H = 1.25+0.75*sin(PI*t);
			data.push_back(pow(x,H));
			data.push_back(pow((1-x),H));
		}else if(problem == "DMOPC")
		{
			double H = 1.25+0.75*sin(PI*t);
			data.push_back(pow(x,H));
			data.push_back(pow((1-x),H));
		}else if(problem == "DMOPD")
		{
			double H = 1.25+0.75*sin(PI*t);
			data.push_back(pow(x,H));
			data.push_back(pow((1-x),H));
		}else if(problem == "DMOPE")
		{
			double H = 1.25+0.75*sin(PI*t);
			data.push_back(pow(x,H));
			data.push_back(pow((1-x),H));
		}else if(problem == "DMOPF")
		{
			double H = 1.25+0.75*sin(PI*t);
			data.push_back(pow(x,H));
			data.push_back(pow((1-x),H));
		}else if(problem == "JY1"){
			data.push_back(x+0.05*sin(6*PI*x));
			data.push_back(1-x+0.05*sin(6*PI*x));
		}
		else if(problem == "JY2"){
			double wt;
			wt=floor(6*sin(0.5*(t-1)));
			data.push_back(x+0.05*sin(wt*PI*x));
			data.push_back(1-x+0.05*sin(wt*PI*x));
		}else if(problem == "JY5"){
			double at;
			at=0.3*sin(0.5*PI*(t-1));
			data.push_back(x+at*sin(PI*x));
			data.push_back(1-x+at*sin(PI*x));
		}else if(problem == "JY3"){
			double alphf,wt,at,y1;
			alphf = floor(100*pow(sin(0.5*PI*t),2));
			wt=floor(6*sin((2*alphf+0.5)*PI*x));
			at=0.05;
			y1=fabs(x*sin(2*alphf+0.5)*PI*x);
			data.push_back(y1+at*sin(wt*PI*y1));
			data.push_back(1-y1+at*sin(wt*PI*y1));
		}else if(problem == "JY4"){
			double at,wt;
			at=0.05;
			wt=pow(10,(1+fabs(sin(0.5*PI*t))));
			data.push_back(x+at*sin(wt*PI*x));
			data.push_back(1-x+at*sin(wt*PI*x));
		}else if(problem == "JY6"){
			double at,wt;
			at=0.1;
			wt=3;
			data.push_back(x+at*sin(wt*PI*x));
			data.push_back(1-x+at*sin(wt*PI*x));
		}else if(problem == "JY7"){
			double at,wt,alphf,beta;
			at=0.1;
			wt=3;
			alphf=0.2+2.8*fabs(sin(0.5*PI*t));
			beta=0.2+2.8*fabs(sin(0.5*PI*t));
			data.push_back(pow((x+at*sin(wt*PI*x)),alphf));
			data.push_back(pow((1-x+at*sin(wt*PI*x)),beta));
		}
		else if(problem == "JY8"){
			double at,wt,alphf,beta;
			at=0.05;
			wt=6;
			beta=10-9.8*fabs(sin(0.5*PI*t));
			alphf=2.0/beta;
			data.push_back(pow((x+at*sin(wt*PI*x)),alphf));
			data.push_back(pow((1-x+at*sin(wt*PI*x)),beta));
		}else if(problem == "JY9"){
			double at,wt;
			int sigma;
			at=0.05;
			sigma = (int) floor(t/5)%3; 
			wt=floor(6*pow(sin(0.5*PI*(t-1)),sigma));
			data.push_back(x+at*sin(wt*PI*x));
			data.push_back(1-x+at*sin(wt*PI*x));

		}else if(problem == "JY10"){
		
		}
		pf.push_back(data);
		x+=0.002;
	}
	if(problem=="JY4" ){
		nodominate(pf);
	}
	//pfile.close();
}

void evaluate::nodominate(vector< vector<double>> &pf){
	int size=pf.size();
	int object=pf[0].size();
	vector<vector<double>> datas;
	int flg;
	for(int i=0;i<size;i++){
		flg=0;
		for(int j=0;j<size && i!=j;j++){
			flg=Dominate(pf[i],pf[j]);
			if(flg==-1){
				break;
			}
		}
		if(flg!=-1){
			datas.push_back(pf[i]);
		}
	}
	pf.clear();
	pf=datas;
}

int evaluate::Dominate(vector<double>index1,vector<double>index2)
{
	int better = 0, wbetter = 0, worse = 0, wworse = 0;
	for(int i=0; i<index1.size(); i++)
	{
		if(	index1[i] <= index2[i])
		{
			better++;
		}
		
		if(index2[i] <= index1[i])
		{
			worse++;
		}
	}

	if(better == index1.size())	return  1;
	else if(worse == index1.size())	return -1;
	else return 0;
}

double evaluate::dist_vector(vector <double> &vec1, vector <double> &vec2)
{
	int dim = vec1.size();
	double sum = 0;
	for(int n=0; n<dim; n++)
		sum+=(vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

void evaluate::newPf(){
	dmop1();
}

void evaluate::JY1(){
	char  filename[1024];
	double rate=0.02;
	double i=0;
	std::ofstream pfile;
	sprintf(filename,"news/pf_%s.dat","JY1");
	//sprintf(filename,"news/FDA2.dat","FDA1");
	pfile.open(filename,ios::out);
	/*while(i<=1.002){
		pfile<<i+0.05*sin(-3*PI*i)<<"    "<<setprecision(5)<<setiosflags(ios::scientific)<<1-i+0.05*sin(-3*PI*i)<<"\n";
		i+=0.002;
	}*/
	int j;
	int P=10000;
	for(i=0;i<=100;i++){
		for(j=0;j<=100;j++){
			double theta=(double)i/99;
			double varphi=(double)j/99;
			pfile<<sin(theta*PI/2)*cos(varphi*PI/2)<<"    "<<sin(theta*PI/2)*sin(varphi*PI/2)<<"   "<<cos(theta*PI/2)<<"\n";
		}
		
	}
	pfile.close();
}
/*


int i,j;
			for(i=0;i<100;i++){
				for(j=0;j<100;j++){
					double theta=(double)i/99;
					double varphi=(double)j/99;
					data.push_back(sin(theta*PI/4)*cos(varphi*PI/4));
					data.push_back(sin(theta*PI/4)*sin(varphi*PI/4));
					data.push_back(cos(theta*PI/4));
					pf.push_back(data);
					data.clear();
				}
			}
*/

void evaluate::JY2(){
	char  filename[1024];
	double rate=0.02;
	double i=0,wt,nt;
	int j;
	nt=0.1;
	for(j=1;j<=120;j++){
		std::ofstream pfile;
		sprintf(filename,"news/pf_%s_%d.dat","JY2",j);
		//sprintf(filename,"news/FDA2.dat","FDA1");
		pfile.open(filename,ios::out);		
		wt=floor(6*sin(0.5*(nt-1)));
		//cout<<wt<<",";
		i=0;
		while(i<=1.002){
			pfile<<i+0.05*sin(wt*PI*i)<<"    "<<setprecision(8)<<1-i+0.05*sin(wt*PI*i)<<"\n";
			i+=0.002;
		}
		nt+=0.1;
		pfile.close();
	}
		
}

void evaluate::JY5(){
	char  filename[1024];
	double i=0,at,nt;
	int j;
	nt=0.1;
	for(j=1;j<=120;j++){
		std::ofstream pfile;
		sprintf(filename,"news/pf_%s_%d.dat","JY5",j);
		pfile.open(filename,ios::out);		
		at=0.3*sin(0.5*PI*(nt-1));
		cout<<at<<",";
		i=0;
		while(i<=1.002){
			pfile<<i+at*sin(PI*i)<<"    "<<setprecision(8)<<1-i+at*sin(PI*i)<<"\n";
			i+=0.002;
		}
		nt+=0.1;
		pfile.close();
	}
		
}

void evaluate::dmop1(){
	char  filename[1024];
	std::ofstream pfile;
	sprintf(filename,"news/pf_%s_%d.dat","fda3",2);
	pfile.open(filename,ios::out);
	int gen=2;
	int nt=10;
	char * problem="fda3";
	double x = 0 ; // 精度
	double t = 0;
	t = (double)(gen-1)/nt;
	while(x<=1){
		vector<double> data;
		if(problem == "fda3"){
			double ft = pow(10,2*sin(0.5*PI*t));
			double gt=fabs(sin(0.5*PI*t));
			pfile<<pow(x,ft)<<"    "<<(1-sqrt(pow(x,ft)/(1+gt)))*(1+gt)<<"\n";
		}
		//pf.push_back(data);
		x+=0.002;
	}
	pfile.close();
}

double l2_norm(vector <double> v)
{
	int size = v.size();
	double sum = 0;
	for (int i = 0; i < size; i++)
		sum += pow(v[i], 2.0);
	return sqrt(sum);
}
