//	AlgD.cpp
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include <random>
#include "AlgD.h"
#include "emo/Sel.h"
#include "alg/AR.h"
#include "alg/Matrix.h"
#include "emo/GenMod.h"
#include "alg/normal.h"
#include "assert.h"
#include "../gd_generate.h"
#include <map>


#include <algorithm>
#include <yvals.h>
#if (_MSC_VER >= 1600)
#define __STDC_UTF_16__
#endif
//#include"engine.h"
#include"memory.h"
//#pragma comment(lib,"libmat.lib") 
//#pragma comment(lib,"libmx.lib") 
//#pragma comment(lib,"libeng.lib") 
//#define SAVE_CEN 1

//!\brief	az namespace, the top namespace
namespace az
{
//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
//!\brief namespace of dynamic evolutionary algoirhtm
namespace dea
{
//const double PI = 3.141592653589793;
double		 T	= 0.0, T0;

// FDAs are from M. Farina, K. Deb and P. Amato's paper
//FDA1
void FDA1(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0, G = sin(0.5*PI*T);
	for(unsigned int i=1; i<X.size(); i++)
		gx += (X[i]-G)*(X[i]-G);
	F[0] = X[0];
	F[1] = gx*(1-sqrt(F[0]/gx));
}
//FDA2
void FDA2(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0, hx = 0.0, H = 0.75+0.7*sin(0.5*PI*T);
	unsigned int i,xii = (unsigned int)(X.size()/2);
	for(i=1; i<xii; i++) gx += X[i]*X[i];
	for(i=xii; i<X.size(); i++) hx += (X[i]-H)*(X[i]-H); hx += H;
	F[0] = X[0];
	F[1] = gx*(1-pow(F[0]/gx, 1.0/hx));
}
//FDA3
void FDA3(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0, Ft = pow(10.0,2.0*sin(0.5*PI*T)), Gt = fabs(sin(0.5*PI*T));
	unsigned int i;
	for(i=1; i<X.size(); i++) gx += (X[i]-Gt)*(X[i]-Gt); gx += Gt;
	F[0] = pow(X[0],Ft);
	F[1] = gx*(1-pow(F[0]/gx, 0.5));
}
//FDA4
void FDA4(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 0.0, G = fabs(sin(0.5*PI*T));
	unsigned int i;
	for(i=2; i<X.size(); i++) gx += (X[i]-G)*(X[i]-G);
	F[0] = (1.0+gx)*cos(0.5*PI*X[0])*cos(0.5*PI*X[1]);
	F[1] = (1.0+gx)*cos(0.5*PI*X[0])*sin(0.5*PI*X[1]);
	F[2] = (1.0+gx)*sin(0.5*PI*X[0]);
}
// == DMOPs are from C.K Goh and K.C Tan 's paper
//dMOP1
void DMOP1(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 0.0, H = 0.75*sin(0.5*PI*T)+1.25;
	unsigned int i;
	for(i=1; i<X.size(); i++) gx += X[i]*X[i]; gx = 1.0 + 9.0*gx;
	F[0] = X[0];
	F[1] = gx*(1-pow(F[0]/gx, H));
}
//dMOP2
void DMOP2(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0, H = 0.75*sin(0.5*PI*T)+1.25, G = sin(0.5*PI*T);
	unsigned int i;
	for(i=1; i<X.size(); i++) gx += (X[i]-G)*(X[i]-G);
	F[0] = X[0];
	F[1] = gx*(1-pow(F[0]/gx, H));
}
//dMOP3
void DMOP3(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0, H = 0.75*sin(0.5*PI*T)+1.25, G = sin(0.5*PI*T);
	unsigned int i;
	for(i=1; i<X.size(); i++) gx += (X[i]-G)*(X[i]-G);
	F[0] = X[0];
	F[1] = gx*(1-pow(F[0]/gx, 0.5));
}

// newly designed DMOPs problems
//centre moving functions
void LT(double t, unsigned int index, double& x1, double& x2)   
{
	t *= 0.5;
	switch(index)
	{	
	case 0:	// Lissajous curve
		x1 = (cos(t*PI)+1.0)*2.0;
		x2 = (sin(2.0*t*PI)+1.0)*2.0;
		break;
	case 1:	// Rose curve
		x1 = (cos(1.5*t*PI)*sin(0.5*t*PI)+1.0)*2.0;
		x2 = (cos(1.5*t*PI)*cos(0.5*t*PI)+1.0)*2.0;
		break;
	case 2: // Heart curve
		x1 = ((1.0-sin(t*PI))*sin(t*PI)+2.0)*1.7;
		x2 = ((1.0-sin(t*PI))*cos(t*PI)+1.5)*1.4;
		break;
	case 3: // discontinus Lissajous curve
		t  = t-floor(t);
		x1 = (cos(t*PI)+1.0)*2.0;
		x2 = (sin(2.0*t*PI)+1.0)*2.0;
		break;
	default:
		x1 = x2 = 0.0;
		break;
	}
}
//shape moving functions
double HT(double t)
{
	return 1.25+0.75*sin(t*PI);
}
// problem framework
void DMOP(std::vector< double >& F, std::vector< double >& X, unsigned int index, bool turn)
{
	double a, b, Gi, gx1=0.0, gx2=0.0, ht=HT(T);
	
	bool old = (((unsigned int)(T*10.0+0.001)) % 2 == 1);
	
	LT(T, index, a, b);
	
	for(unsigned int i=1; i<X.size(); i++)
	{
		//if(turn && old)
		//	Gi = pow(fabs(2.0*X[0]-2.0*a-1.0), ht+double(i-1.0)/double(X.size()-2.0));
		//else	
		//	Gi = 1.0 - pow(fabs(2.0*X[0]-2.0*a-1.0), ht+double(i-1.0)/double(X.size()-2.0));

		if(turn && old)
			Gi = pow(fabs(X[0]-a), ht+double(i+1.0)/double(X.size()));
		else	
			Gi = 1.0 - pow(fabs(X[0]-a), ht+double(i+1.0)/double(X.size()));


		if(i % 2 == 1)
			gx1 += pow(X[i] - b - Gi ,2.0);
		else 
			gx2 += pow(X[i] - b - Gi ,2.0);
	}
	F[0] = pow(fabs(X[0]-a), ht)     + 0.5*gx1;
	F[1] = pow(fabs(X[0]-a-1.0), ht) + 0.5*gx2;
}
void DMOPA(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X) //F5
{
	DMOP(F,X,0,false);
}
void DMOPB(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)//F6
{
	DMOP(F,X,1,false);
}
void DMOPC(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)//F7
{
	DMOP(F,X,2,false);
}
void DMOPD(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X) //F9
{
	DMOP(F,X,3,false);
}
void DMOPE(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)//F10
{
	DMOP(F,X,0,true);
}
void DMOPF(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)//F8
{	
	double gx = 0.0, G = sin(0.5*PI*T);
	for(unsigned int i=2; i<X.size(); i++) gx += pow(X[i]-pow(0.5*(X[0]+X[1]), HT(T)+double(i+1.0)/double(2.0*X.size()))-G,2.0);
	F[0] = (1.0+gx)*cos(0.5*PI*X[0])*cos(0.5*PI*X[1]);
	F[1] = (1.0+gx)*cos(0.5*PI*X[0])*sin(0.5*PI*X[1]);
	F[2] = (1.0+gx)*sin(0.5*PI*X[0]);
}
//多项式偏转：P(t)
double Pt()
{
	return 1.0 + 16.0 * pow(cos(0.5 * PI * T) , 4);
}

//欺骗函数中的A(t)与多模函数中的C(t)
double tempT()
{
	return 0.5+0.35*sin(PI*T);
}

//欺骗函数  A(t),B=0.001,C=0.1
double S_decept(double x)
{
	double A = tempT(),B = 0.001,C = 0.1;
	assert(!(A-B));assert(!(1-A-B));
	double temp1 = 0,temp2 = 0;
	temp1 = floor(x-A+B)*(1-C+(A-B)/B)/(A-B);
	temp2 =  floor(A+B-x)*(1-C+(1-A-B)/B)/(1-A-B);
	return 1+(fabs(x-A)-B)*(temp1+temp2+1/B);
}

//多模函数  A=20,B=10,C(t)
double S_multi(double x)
{
	double A = 20,B = 10,C = tempT();
	double temp = fabs(x-C)/(2*(floor(C-x)+C));
	return (1+cos((4*A+2)*PI*(0.5-temp))+4*B*pow(temp,2))/(B+2);
}

//SDMOP1
void SDMOP1(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0,Gt = cos(0.5*PI*T);
	for(unsigned int i = 1;i<X.size();i++)
		gx += pow(X[i]-Gt,2);
	F[0] = X[0];
	F[1] = gx*(1-pow(F[0]/gx,0.5));
}

//SDMOP2
void SDMOP2(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0,Ht = 1.2+0.8*cos(0.5*PI*T);
	for(unsigned int i=1;i<X.size();i++)
		gx += 9*X[i]*X[i];
	F[0] = pow(X[0],Pt());
	//F[0]=X[0];
	F[1] = gx*(1-pow(F[0]/gx,Ht));
}

//SDMOP3
void SDMOP3(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0,Gt = cos(0.5*PI*T),Ht = 1.2+0.8*cos(0.5*PI*T);
	for(unsigned int i=1;i<X.size();i++)
		gx += (X[i]-Gt)*(X[i]-Gt);
	F[0] = S_decept(X[0]);
	F[1] = gx*(1-pow(F[0]/gx,Ht));
}

//SDMOP4
void SDMOP4(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	double Gt = cos(0.5*PI*T),gx = 1.0+fabs(Gt),sum = 0,mulit = 0;
	unsigned int i;
	for(i=0;i<5;i++)
	{sum += i;mulit += X[i]*i;}
	for(i=5;i<X.size();i++)
		gx += (X[i]-Gt)*(X[i]-Gt);
	F[0] = mulit/sum;
	//F[0] = X[0];
	F[1] = gx*(1-F[0]-cos(10*PI*F[0]+0.5*PI)/(10*PI));
}

//SDMOP5
void SDMOP5(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	assert(1);
	unsigned int M = F.size(),N = X.size(),i,j;
	double gx = 0.0,temp,Gt = cos(0.5*PI*T);
	for(i=M-1;i<N;i++)
		//gx += (X[i]-Gt)*(X[i]-Gt);
		gx += S_multi(X[i])*S_multi(X[i]);
	for(i=0;i<M;i++)
	{
		temp = 1.0;
		if(i==0)
			for(j=0;j<M-1;j++)
				temp *= cos(0.5*PI*X[j]);
		else
			if(i==M-1)
				temp *= sin(0.5*PI*X[0]);
			else
				{
					temp = sin(0.5*PI*X[M-i-1]);
					for(j=0;j<M-i-1;j++)    temp *= cos(0.5*PI*X[j]); 
				}
		F[i] = (1+gx)*temp;
	}
	/*F[0] = (1.0+gx)*cos(0.5*PI*X[0])*cos(0.5*PI*X[1]);
	F[1] = (1.0+gx)*cos(0.5*PI*X[0])*sin(0.5*PI*X[1]);
	F[2] = (1.0+gx)*sin(0.5*PI*X[0]);*/
}

//SDMOP6
void SDMOP6(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	assert(1);
	unsigned int M = F.size(),N = X.size(),i,j;
	double temp,Gt = cos(0.5*PI*T),gx = 0.0;
	std::vector< double > Z(M);
	for(i=0;i<M-1;i++)
		Z[i] = pow(X[i],Pt());
	for(i=M-1;i<N;i++)
		gx += (X[i]-Gt)*(X[i]-Gt);
	for(i=0;i<M;i++)
	{
		temp = 1.0;
		if(i==0)
			for(j=0;j<M-1;j++)
				temp *= cos(0.5*PI*Z[j]);
		else
			if(i==M-1)
				temp *= sin(0.5*PI*Z[0]);
			else
				{
					temp = sin(0.5*PI*Z[M-i-1]);
					for(j=0;j<M-i-1;j++)    temp *= cos(0.5*PI*Z[j]); 
				}
		F[i] = (1+gx)*temp;
	}
}

//SDMOP7
void SDMOP7(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	assert(1);
	unsigned int M = F.size(),N = X.size(),i,j;
	double temp,Gt = cos(0.5*PI*T),gx = fabs(Gt);
	std::vector< double > Z(M);
	for(i=0;i<M-1;i++)
		Z[i] = S_decept(X[i]);
	for(i=M-1;i<N;i++)
		gx += (X[i]-Gt)*(X[i]-Gt);
	for(i=0;i<M;i++)
	{
		temp = 1.0;
		if(i==0)
			for(j=0;j<M-1;j++)
				temp *= cos(0.5*PI*Z[j]);
		else
			if(i==M-1)
				temp *= sin(0.5*PI*Z[0]);
			else
				{
					temp = sin(0.5*PI*Z[M-i-1]);
					for(j=0;j<M-i-1;j++)    temp *= cos(0.5*PI*Z[j]); 
				}
		F[i] = (1+gx)*temp;
	}
}

//SDMOP8
void SDMOP8(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	assert(1);
	unsigned int M = F.size(),N = X.size(),i,j;
	double temp = 0,Ht = fabs(cos(0.5*PI*T)),gx = 0;
	for(i=2;i<N;i++)
		gx += 9*X[i]*X[i];
	F[0] = X[0];F[1] = X[1];

	for(i=2;i<M;i++)
	{
		temp = (X[0]/(1+gx))*(1+Ht*cos(PI*X[0]*X[0]*X[0])*cos(PI*X[0]*X[0]*X[0])) + (X[1]/(1+gx))*(1+Ht*cos(PI*X[1]*X[1]*X[1])*cos(PI*X[1]*X[1]*X[1]));
		F[i] = (1+gx)*(i-1)*(2-0.5*temp);
	}
}

void JY1(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	unsigned int j, count1, count2;
		double a, w, g, sum1;

		a = 0.05, w = 6;
		g = sin(0.5*PI*T);
		sum1 = 0;
		for (j = 1; j < X.size(); j++)
		{
			double x = X[j];//2 * X[j] - 1.0;
			double yj = x - g;
			yj = yj*yj;
			sum1 += yj;
		}
		F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
		F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));

		return;
}

void JY2(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	unsigned int j, count1, count2, nvar = X.size();
		double a, w, g, sum1;

		a = 0.05, w = floor(6 * sin(0.5*PI*(T - 1)));
		g = sin(0.5*PI*T);
		sum1 = 0;

		for (j = 1; j < nvar; j++)
		{
			double x = 2 * X[j];//X[j] - 1.0;
			double yj = x - g;
			yj = yj*yj;
			sum1 += yj;
		}
		F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
		F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));

		return;
}

void JY3(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	unsigned int j, count1, count2, nvar = X.size();
		double a, aa, y1, w, g, sum1;

		aa = floor(100 * sin(0.5*PI*T)*sin(0.5*PI*T));
		y1 = fabs(X[0] * sin((2 * aa + 1)*PI*X[0]));

		a = 0.05, w = floor(6 * sin(0.5*PI*(T - 1)));
		//g=sin(0.5*PI*Tstep);

		sum1 = 0;

		for (j = 1; j < nvar; j++)
		{
			double x = 2 * X[j];//X[j] - 1.0;
			double x0 = 2 * X[j - 1] - 1.0;
			double yj;
			if (j == 1){ yj = x*x - y1; }
			else{ yj = x*x - x0; }

			yj = yj*yj;
			sum1 += yj;
		}
		F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
		F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));

		return;
}

void JY4(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	unsigned int j, count1, count2, nvar = X.size();
		double a, w, g, sum1;

		g = sin(0.5*PI*T);
		a = 0.05, w = pow(10.0, 1.0 + fabs(g));

		sum1 = 0;

		for (j = 1; j < nvar; j++)
		{
			double x = 2 * X[j];//X[j] - 1.0;
			double yj = x - g;
			yj = yj*yj;
			sum1 += yj;
		}
		F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
		F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));

		return;
}

void JY5(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	unsigned int j, count1, count2, nvar = X.size();
		double a, w, g, sum1;

		g = sin(0.5*PI*T);
		a = 0.3*sin(0.5*PI*(T - 1)), w = 1.0;

		sum1 = 0;

		for (j = 1; j < nvar; j++)
		{
			double x = 2 * X[j];//X[j] - 1.0;
			double yj = x;
			yj = yj*yj;
			sum1 += yj;
		}
		F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
		F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));

		return;
}

void JY6(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	unsigned int j, count1, count2, nvar = X.size();
		double a, w, g, k, sum1;

		g = sin(0.5*PI*T);
		a = 0.1, w = 3;
		k = 2 * floor(10 * fabs(g));

		sum1 = 0;

		for (j = 1; j < nvar; j++)
		{
			double x = 2 * X[j];//X[j] - 1.0;
			double yj = x - g;
			yj = 4 * yj*yj - cos(k*PI*yj) + 1;
			sum1 += yj;
		}
		F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
		F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));

		return;
}

void JY7(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	unsigned int j, count1, count2, nvar = X.size();
		double a, w, g, s, t, sum1;

		g = sin(0.5*PI*T);
		a = 0.1, w = 3;
		s = t = 0.2 + 2.8*fabs(g);

		sum1 = 0;

		for (j = 1; j < nvar; j++)
		{
			double x = 2 * X[j];//X[j] - 1.0;
			double yj = x - g;
			yj = yj*yj - 10 * cos(2 * PI*yj) + 10;
			sum1 += yj;
		}
		F[0] = (1 + sum1)*pow(X[0] + a*sin(w*PI*X[0]), s);
		F[1] = (1 + sum1)*pow(1 - X[0] + a*sin(w*PI*X[0]), t);

		return;
}

void JY8(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	unsigned int j, count1, count2, nvar = X.size();
		double a, w, g, s, t, sum1;

		g = sin(0.5*PI*T);
		a = 0.05, w = 6;
		t = 10.0 - 9.8*fabs(g);
		s = 2.0 / t;

		sum1 = 0;

		for (j = 1; j < nvar; j++)
		{
			double x = 2 * X[j];//X[j] - 1.0;
			double yj = x;
			yj = yj*yj;
			sum1 += yj;
		}
		F[0] = (1 + sum1)*pow(X[0] + a*sin(w*PI*X[0]), s);
		F[1] = (1 + sum1)*pow(1 - X[0] + a*sin(w*PI*X[0]), t);

		return;
}

void JY9(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
{
	unsigned int j, d, count1, count2, nvar = X.size();
		double a, w, g, sum1;

		d = (int)floor(T/5) % 3;
		g = floor(sin(0.5*PI*T));
		a = 0.05, w = floor(6 * pow(sin(0.5*PI*(T - 1)), d));

		sum1 = 0;

		for (j = 1; j < nvar; j++)
		{
			double x = X[j];//2 * X[j] - 1.0;
			double yj = x + d - g;
			yj = yj*yj;
			sum1 += yj;
		}
		F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
		F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));

		return;
}

void DFRange(std::vector<double>& low, std::vector<double>& upp, std::string& name ,unsigned int mObj)
{
	if(	name == std::string("FDA1")   || name == std::string("FDA2")   || name == std::string("FDA3")||
		name == std::string("DMOP1")  || name == std::string("DMOP2")  || name == std::string("DMOP3")||
		name == std::string("SDMOP1") || name == std::string("SDMOP2") || name == std::string("SDMOP3") ||
		name == std::string("JY1")    || name == std::string("JY2")    || name == std::string("JY3")    ||
		name == std::string("JY4")    || name == std::string("JY5")    || name == std::string("JY6")   
		|| name == std::string("JY7")  || name == std::string("JY8")   || name == std::string("JY9")
		)
	{
		low[0] = 0;upp[0] =  1;
		for(unsigned int i=1; i<(unsigned int)(low.size()); i++){low[i] = -1; upp[i] = 1;}
	}
	else if(name == std::string("FDA4") || name == std::string("SDMOP5") || name == std::string("SDMOP8"))
	{
		low[0] = 0;upp[0] =  1;
		for(unsigned int i=1; i<(unsigned int)(low.size()); i++){low[i] =  0; upp[i] = 1;}
	}
	else if(name == std::string("DMOPF"))
	{
		low[0] = 0; upp[0] =  1;
		low[1] = 0; upp[1] =  1;
		for(unsigned int i=2; i<(unsigned int)(low.size()); i++){low[i] = -1.0; upp[i] = 2.0;}
	}
	else if(name == std::string("SDMOP4"))
	{
		for(unsigned int i=0;i<5;i++){low[i] = 0;upp[i] = 1.0;}
		for(unsigned int i=5;i<(unsigned int)(low.size());i++){low[i] = -1.0;upp[i] = 1.0;}
	}
	else if(name == std::string("SDMOP6") || name == std::string("SDMOP7"))
	{
		for(unsigned int i=0;i<mObj-1;i++){low[i] = 0;upp[i] = 1.0;}
		for(unsigned int i=mObj-1;i<(unsigned int)(low.size());i++){low[i] = -1.0;upp[i] = 1.0;}
	}
	else
	{	//DMOPA - DMOPE
		for(unsigned int i=0; i<(unsigned int)(low.size()); i++){low[i] = 0.0; upp[i] = 5.0;}
	}
}

DMOO::DMOO(
	unsigned int	strategy,
	std::string&	optimizer,
	unsigned int	popsize	,
	unsigned int	stepmax	,
	unsigned int	taot	,
	unsigned int	nt		,
	unsigned int	torder	,
	double			t0		,
	double			alpha   ,
	CParameter&		par		,
	std::vector<int> &t_order)
	:mPop(par),mBest0(par),mBest1(par),observe(par),memory(par)
{
	switch(strategy)
	{
	case 1:
		mStrategy = INIRIS;
		break;
	case 2:
		mStrategy = INIFPS;
		break;
	case 3:
		mStrategy = INIPPS;
		break;
	case 4:
		mStrategy = INIPZ;
		break;
	case 5:
		mStrategy = INIEGS;
		break;
	case 6:
		mStrategy = INIRG;
		break;
	case 7:
		mStrategy = INIDSP;
		break;
	default:
		mStrategy = INIRIS;
		break;
	}
	mOptimizer = optimizer;
	numMemory = 0;
	mPopSize = popsize;
	mMaxStep = stepmax;
	mTaoT	 = taot;
	mDelT	 = 1.0/nt;
	T0		 = mT0 = t0;
	mMaxOrder= torder;
	mAlpha   = alpha;
	pPar	 = &par;
	for(unsigned int i=0;i<t_order.size();i++)
		T_order.push_back(t_order[i]);

	if(pPar->Problem() == std::string("FDA1"))
	{
		P().Evaluator( FDA1 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("FDA2"))
	{
		P().Evaluator( FDA2 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("FDA3"))
	{
		P().Evaluator( FDA3 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("FDA4"))
	{
		P().Evaluator( FDA4 );
		P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOP1"))
	{
		P().Evaluator( DMOP1 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOP2"))
	{
		P().Evaluator( DMOP2 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOP3"))
	{
		P().Evaluator( DMOP3 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOPA"))
	{
		P().Evaluator( DMOPA );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOPB"))
	{
		P().Evaluator( DMOPB );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOPC"))
	{
		P().Evaluator( DMOPC );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOPD"))
	{
		P().Evaluator( DMOPD );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOPE"))
	{
		P().Evaluator( DMOPE );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOPF"))
	{
		P().Evaluator( DMOPF );
		P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("SDMOP1"))
	{
		P().Evaluator( SDMOP1 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("SDMOP2"))
	{
		P().Evaluator( SDMOP2 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("SDMOP3"))
	{
		P().Evaluator( SDMOP3 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("SDMOP4"))
	{
		P().Evaluator( SDMOP4 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("SDMOP5"))
	{
		P().Evaluator( SDMOP5 );
		P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("SDMOP6"))
	{
		P().Evaluator( SDMOP6 );
		P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("SDMOP7"))
	{
		P().Evaluator( SDMOP7 );
		P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("SDMOP8"))
	{
		P().Evaluator( SDMOP8 );
		P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("JY1"))
	{
		P().Evaluator( JY1 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("JY2"))
	{
		P().Evaluator( JY2 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("JY3"))
	{
		P().Evaluator( JY3 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("JY4"))
	{
		P().Evaluator( JY4 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("JY5"))
	{
		P().Evaluator( JY5 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("JY6"))
	{
		P().Evaluator( JY6 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("JY7"))
	{
		P().Evaluator( JY7 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("JY8"))
	{
		P().Evaluator( JY8 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("JY9"))
	{
		P().Evaluator( JY9 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	Reset();
}

void DMOO::Reset()
{
	mStep		= 0;
	mEvas		= 0;
	T			= mT0;
	mbToChange	= false;
	hC.clear();
	DFRange(P().XLow(), P().XUpp(), P().Problem(),P().FSize());	// set the original search space

	g.resize(P().XSize(),0.0);
	az::rnd::seed((long) time(NULL));
}

// main evolution step
unsigned int DMOO::Step()
{
	// environment change
	if(mStep%mTaoT == 0)
	//if(EnvironmentChange())
	{
		//std::cout<<"环境变化程度："<<EnvironmentChangeDegree()<<std::endl;
		//EnvironmentChangeDegree();
		//std::cout<<"环境变化"<<std::endl;
		//if(mStep == 0) {InitRIS(true);}
		switch(mStrategy)
		{
			case INIRIS:
				InitRIS(true);
				break;
			case INIFPS:
				InitFPS();
				break;
			case INIPPS:
				InitPPS();
				break;
			case ININN:
				InitNNPS();
			case INIPZ:
				{
					InitGIPS();
					/*InitPZ();
					CPopulationMO popnew(P());
					Generate(popnew, mPopSize);
					Check(popnew);
					popnew.Evaluate(); mEvas += popnew.Size();
					Population()=popnew;*/
					if (mStep!=0)
					{
						//ObservePZ();
						mStep = mStep + 1;
					}
					break;
				}
			case INIRG:
				{
					InitRG();
					break;
				}
			case INIEGS:
				InitEGS();
				mStep = mStep+2;
				break;
			case INIDSP:
				{
					InitDSP();
					break;
				}
		}
		Population().Evaluate(); 
		mEvas += Population().Size();
	}
	// no environment change
	else
	{
		//std::cout<<"环境不变化"<<std::endl;
		CPopulationMO popnew(P());

		// generate new solutions
		Generate(popnew, mPopSize);

		//check the range of new solutions
		Check(popnew);
		//remove these repeated solutions
		Check(popnew, Population());

		//evaluate new solutions
		popnew.Evaluate(); mEvas += popnew.Size();

		/*if(mStrategy == INIPPS || mStrategy == INIEGS)
			LocalSearchEGS(popnew);*/
		//environmental select
		Population().Combine(popnew);
	
		az::mea::sel::SCrowd2 sel;
		sel.Select(Population(), mPopSize);
	}
	mStep++;

	//!!!! CHANGE STATE IF POSSIBLE
	// the change happens in (mTaoT*k + 1)th generation, k=1,2,...
	if(mStep > 1 && (mStep % mTaoT == 0))
	{
		if(P().Problem().substr(0,5) == "SDMOP")
			T=mDelT*T_order[(mStep/mTaoT)%21];								//update time T in SDMOP series instances 
		else
			T += mDelT;														// update time T

		mbToChange = true;
	}
	else
	{
		mbToChange = false;
	}

	return mStep;
}

// check the range of each individual
void DMOO::Check(CPopulationMO& pop)
{
	for(unsigned int i=0; i<pop.Size(); i++) pop[i].Check();
}

// delete duplicate individuals
void DMOO::Check(CPopulationMO& popnew, CPopulationMO& pop)
{
	CPopulationMO tmp(P());
	for(unsigned int i=0; i<popnew.Size(); i++) if(!pop.IsContain(popnew[i])) tmp.Combine(popnew.At(i));
	popnew.Clear();
	popnew.Combine(tmp);
}

bool DMOO::EnvironmentChange()
{
	if(mStep == 0) return true;

	CIndividualMO ind(P());
	unsigned int i, j, in, cnum = (unsigned int)(Population().Size()*0.05); if(cnum<1) cnum = 1;
	for(i=0; i<cnum; i++)
	{
		in  = rnd::rand((unsigned int)0, Population().Size());
		ind = Population()[in];
		ind.Evaluate(); mEvas++;
		for(j=0; j<Population().P().FSize(); j++)
			if(fabs(ind.F(j)-Population()[in].F(j))>1.0E-10) 
				return true;
	}
	return false;
}

//环境变化程度
void DMOO::EnvironmentChangeDegree(CPopulationMO& pop)
{
	unsigned int i,j,k;
	double tempVar = 0.0, tempInd = 0.0;
	unsigned int N = pop.Size(),M = pop.P().FSize();
	CPopulationMO temp_pop(pop);
	pop.degreeAssign(0.0);
	std::vector<std::vector <double>> mObj;
	std::vector<double> temp;

	std::vector<double> tempVec(N);

	for(i=0;i<N;i++)
	{
		for(j=0;j<M;j++)
		{
			temp.push_back(temp_pop[i].F(j));
		}
		mObj.push_back(temp);
		temp.clear();
	}

	temp_pop.Evaluate();

	for(i=0;i<N;i++)
	{
		for(j=0;j<M;j++)
		{
			tempInd += fabs(temp_pop[i].F(j)-mObj[i][j]);
		}	
		pop[i].degreeAssign(tempInd);
		tempVar += pop[i].indDegree;
	}
	tempVar /= N;
	pop.degreeAssign(tempVar);
}

// predict the location of the next centre point
void DMOO::Predict()
{
	unsigned int i, k, d, order, dim=(unsigned int)(mC.size());
	std::list< std::vector<double> >::iterator it, it0;

	while(hC.size()>20+mMaxOrder) hC.pop_back();

	order = std::min((unsigned int)hC.size()-1, (unsigned int)mMaxOrder);
	double **px, **pa, *pv;
	px = new double*[dim]; for(i=0; i<dim; i++) px[i] = new double[hC.size()];
	pa = new double*[dim]; for(i=0; i<dim; i++) pa[i] = new double[order];
	pv = new double[dim];
	k  = (unsigned int)hC.size()-1;
	it = hC.begin();
	d  = 0;
	while(it!=hC.end())
	{
		for(i=0; i<dim; i++) px[i][k] = (*it)[i];
		k--;
		it++;
		d++;
	}
	alg::aruv(px, dim, (unsigned int)hC.size(), order, pa, pv);
	mStdC.resize(dim); for(i=0; i<dim; i++) mStdC[i] = pv[i];

	pC.resize(mC.size()); for(i=0; i<dim; i++) pC[i] = 0.0;
	it = hC.begin();
	for(k=0; k<order; k++)
	{
		for(i=0; i<dim; i++) pC[i] += (*it)[i]*pa[i][k];
		it++;
	}
	it = hC.begin();
	for(i=0; i<dim; i++)
	{
		// keep it in the boundary
		if(     pC[i]>P().XUpp(i%P().XSize())) pC[i] = (*it)[i]; //P().XUpp(i%P().XSize());
		else if(pC[i]<P().XLow(i%P().XSize())) pC[i] = (*it)[i]; //P().XLow(i%P().XSize());
		if(     pC[i]>P().XUpp(i%P().XSize())) pC[i] = P().XUpp(i%P().XSize());
		else if(pC[i]<P().XLow(i%P().XSize())) pC[i] = P().XLow(i%P().XSize());
	}

	for(i=0; i<dim; i++) delete []px[i]; delete []px;
	for(i=0; i<dim; i++) delete []pa[i]; delete []pa;
	delete []pv;
}

void DMOO::InitEGS()
{
	if(mStep == 0) {InitRIS(true); return;}

	unsigned int i, k, xdim=P().XSize(), psize = Population().Size(); 
	int seed_init;
	int num1 = 0;
	std::vector<std::vector<double>> infor;
	mC.resize(xdim);
	g.resize(xdim,0.0);
	variance.resize(xdim);
	seed_init = 123456789 + rnd::rand(0,1000);

	CPopulationMO	old_pop(P());
	old_pop = Population();
	for(k=0; k<xdim; k++) 
	{
		mC[k]  = 0.0;
		num1 = 0;
		for(i=0; i<psize; i++)
			if (old_pop[i].Rank()==1)
			{
				mC[k] += old_pop[i][k];
				num1++;
			}
			mC[k] /= double(num1);
	}
	for(k=0; k<xdim; k++) 
	{
		variance[k]  = 0.0;
		for(i=0; i<psize; i++)
			if (old_pop[i].Rank()==1)
				variance[k] += (old_pop[i][k]-mC[k])*(old_pop[i][k]-mC[k]);
		variance[k] /= double(num1);
	}
	if (mStep < 3*mTaoT)
	{
		oldC2.push_back(mC);
		infor.push_back(mC);
		infor.push_back(variance);
		memoryC.push_back(infor);
		if (memoryC.size()>psize)
		{
			memoryC.erase(memoryC.begin(),memoryC.begin()+1);
		}
		infor.clear();
		return;
	}

	//exploit(开采)
	az::mea::sel::SCrowd2 sel;
	sel.Select(Population(), 40);
	for (k = 0; k < 40; k++)
	{
		for (i = 0;i < xdim; i++)
		{
			Population()[k][i] += normal::r4_normal(0,0.01,&seed_init);   //N（0,0.01）
			if (Population()[k][i]>P().XUpp(i))
			{
				Population()[k][i] = P().XUpp(i) - (Population()[k][i]-P().XUpp(i))/10.0;
			}
			if (Population()[k][i]<P().XLow(i))
			{
				Population()[k][i] = P().XLow(i) + (P().XLow(i)-Population()[k][i])/10.0;
			}

		}
	}
	Population().Evaluate();

	//predict

	CPopulationMO	archive(P());
	std::vector<size_t> vnIndex(psize);
	for (size_t j = 0; j < vnIndex.size(); ++j)
		vnIndex[j] = j;
	std::random_shuffle(vnIndex.begin(), vnIndex.end());
	for (size_t j = 0; j < 30; ++j)
		archive.Copy(old_pop[vnIndex[j]]);
	int u =1,negative_times = 0;
	//有问题：k值此时为40，但是决策向量维数为20
	//std::cout<<"k值："<<k<<"\t"<<"决策维数"<<xdim<<std::endl;
	for(k=0;k<xdim;k++)
		g[k] = g[k]/2+(mC[k]-oldC2[1][k])/3+(oldC2[1][k]-oldC2[0][k])/6;
	for (i = 0; i < 30; i++)
	{
		for (k = 0; k < xdim; k++)
		{
			archive[i][k] +=  u*g[k];
			if (archive[i][k]>P().XUpp(k))
			{
				archive[i][k] = P().XUpp(k) - (archive[i][k]-P().XUpp(k))/10.0;
			}
			if (archive[i][k]<P().XLow(k))
			{
				archive[i][k] = P().XLow(k) + (P().XLow(k)-archive[i][k])/10.0;
			}
		}
		archive[i].Evaluate();
		if (old_pop[rnd::rand(0,(int)psize)].Dominate(archive[i]))
		{
			u=-u;
		}
		if (u<0)
		{
			negative_times++;
		}
	}
	if (negative_times > 15)
	{
		for(k=0;k<xdim;k++)
			g[k] = -g[k];
	}

	//memory
	std::vector<std::vector<std::vector<double>>>::iterator p,h,ip1,ip2;	//h is head poiter, ip1和ip2指向具有最短路径的两个点
	std::vector<std::vector<std::vector<double>>> duplicate;
	double min,min2,minb,dist;
	duplicate = memoryC;

	for(i=duplicate.size();i>10;i--)
	{
		min = 1.0e14;
		//Delete the crowdest item from Duplicate 
		for(h=duplicate.begin();h!=duplicate.end();++h)
		{
			p = h;
			for(++p;p!=duplicate.end();++p)
			{
				dist = 0.0;
				for(int j=0;j<xdim;j++)
				{
					dist += ((*h)[1][j]-(*p)[1][j]) * ((*h)[1][j]-(*p)[1][j]);
				}
				if (dist<min)
				{
					min = dist;
					ip1 = h;
					ip2 = p;
				}
			}
		}
		//计算ip1到其他点的第二极短距离min2
		min2 = 1.0e14;
		for(p=duplicate.begin(); p!=duplicate.end(); ++p)
		{
			if(p!=ip1 && p!=ip2)
			{
				dist = 0.0;
				for(int j=0;j<xdim;j++)
				{
					dist += ((*ip1)[1][j]-(*p)[1][j]) * ((*ip1)[1][j]-(*p)[1][j]);
				}
				if (dist<min2)
				{
					min2 = dist;
				}
			}
			else
				continue;
		}

		//计算ip2到其他点的第二极短距离minb
		minb = 1.0e14;
		for(p=duplicate.begin(); p!=duplicate.end(); ++p)
		{
			if(p!=ip1 && p!=ip2)
			{
				dist = 0.0;
				for(int j=0;j<xdim;j++)
				{
					dist += ((*ip2)[1][j]-(*p)[1][j]) * ((*ip2)[1][j]-(*p)[1][j]);
				}
				if (dist<minb)
				{
					minb = dist;
				}
			}
			else
				continue;
		}
		//比较min2和minb，找更短次短距离，将其删除
		if(min2<minb)
			duplicate.erase(ip1);
		else if(min2>minb)
			duplicate.erase(ip2);
		else
		{
			if(rnd::rand(0,2)==0 )
				duplicate.erase(ip1);
			else
				duplicate.erase(ip2);
		}	
	}

	CPopulationMO	old_pop2(P());
	old_pop2 = Population();
	seed_init = 123456789 + rnd::rand(0,1000);
	for(p=duplicate.begin();p!=duplicate.end();++p)
	{
		for (i=0;i<1;i++)	//Since Roffspring==1, only one loop in fact!!!
		{
			for(int j=0;j<xdim;j++)
			{
				do {
					Population()[0][j] = normal::r4_normal((*p)[0][j],sqrt((*p)[1][j]),&seed_init);
				} while (Population()[0][j]<P().XLow(j)|| Population()[0][j]>P().XUpp(j));	
			}
		}
		archive.Copy(Population()[0]);
	}
	archive.Evaluate();
	Population()=old_pop2;

	oldC2.push_back(mC);
	if (oldC2.size()>2)
	{
		oldC2.erase(oldC2.begin());
	}
	infor.push_back(mC);
	infor.push_back(variance);
	memoryC.push_back(infor);
	if (memoryC.size()>psize)
	{
		memoryC.erase(memoryC.begin(),memoryC.begin()+1);
	}
	infor.clear();

	Population().Combine(archive);
	Population().Combine(old_pop);
	Population().Evaluate();
	az::mea::sel::SCrowd2 sel2;
	sel2.Select(Population(), mPopSize);
}

void DMOO::InitPZ()
{	
	//Step 0: initialize the population

	if(mStep == 0) {InitRIS(true); return;}

	CPopulationMO	old_pop(P());
	old_pop = Population();
	unsigned int  xdim=P().XSize(), psize = Population().Size();
	int der = 40;
	//Step 1: save center point and shape points
	mC.resize(xdim);
	oldC.resize(xdim,0);
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	int num1=0,num_nodominance=0;
	for (int i=0; i<psize; i++)
	{
		if (Population()[i].Rank()==1)
			num1++;
	}
	num_nodominance=num1;
	if (num1<psize/20)
	{
		for (int i=0;i<psize/20;i++)
			memory.Copy(Population()[i]);
	}	
	else
	{
		std::vector<size_t> vnIndex(num1);
		for (size_t j = 0; j < vnIndex.size(); ++j)
			vnIndex[j] = j;
		std::random_shuffle(vnIndex.begin(), vnIndex.end());
		for (size_t j = 0; j < psize/20; ++j)
			memory.Copy(Population()[vnIndex[j]]);
	}
	if(memory.Size()>2*psize) 
	{
		memory.Pop(memory.Size()-2*psize);
	}
	for(int k=0; k<xdim; k++) 
	{
		mC[k]  = 0.0;
		for(int i=0; i<num1; i++)
			if (Population()[i].Rank()==1)
				mC[k] += Population()[i][k];
		mC[k] /= double(num1);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//DEE-PDMS  基于动态环境进化模型的种群多样性保持策略
	double temp_minF[2];
	double temp_maxF[2];
	for(int i=0;i<2;i++)
	{
		temp_minF[i]=Population()[0].F(i);
		temp_maxF[i]=Population()[0].F(i);
		for (int n = 1 ; n < psize ; n++)
		{
			if (temp_minF[i]>Population()[n].F(i))
				temp_minF[i] = Population()[n].F(i);
			if (temp_maxF[i]<Population()[n].F(i))
				temp_maxF[i] = Population()[n].F(i);
		}
	}

	/*std::cout<<temp_minF[0]<<"\t"<<temp_minF[1]<<std::endl
			 <<temp_maxF[0]<<"\t"<<temp_maxF[1]<<std::endl;*/

	old_pop.Evaluate();
	for(int i=0;i<2;i++)
	{
		for (int n = 1 ; n < psize ; n++)
		{
			if (temp_minF[i]>old_pop[n].F(i))
				temp_minF[i] = old_pop[n].F(i);
			if (temp_maxF[i]<old_pop[n].F(i))
				temp_maxF[i] =old_pop[n].F(i);
		}
	}

	/*std::cout<<std::endl<<temp_minF[0]<<"\t"<<temp_minF[1]<<std::endl
		 <<temp_maxF[0]<<"\t"<<temp_maxF[1]<<std::endl;*/
	for (int i = 0; i < 2; i++)
		for(int j = 0; j < psize; j++)
			Population()[j].old_GD(i) = (int)((Population()[j].F(i)-temp_minF[i])/((temp_maxF[i]-temp_minF[i])/der));

	Population().Evaluate();
	/*Population().IsSort(false);
	Population().RankSort(false);*/
	for (int i = 0; i < 2; i++)
		for(int j = 0; j < psize; j++)
			Population()[j].new_GD(i) = (int)((Population()[j].F(i)-temp_minF[i])/((temp_maxF[i]-temp_minF[i])/der));
	
	/*Engine *ep;
	static bool flag2 =false ;
	double f3[100];
	double f4[100];

	mxArray *T3 = mxCreateDoubleMatrix( 1 , 100 , mxREAL  ) ; 
	mxArray *T4 = mxCreateDoubleMatrix( 1 ,  100 , mxREAL  ) ;
	mxArray *M1 = mxCreateDoubleMatrix( 1 ,  1 , mxREAL  ) ;
	for(int n=0; n<100; n++)
	{ 
	f3[n]=Population()[n].F(0);
	f4[n]=Population()[n].F(1);
	}
	if(  (ep=engOpen(NULL)) )
	{
	memcpy( (char*)mxGetPr(T3),(char*)f3 , 100*sizeof(double));
	memcpy( (char*)mxGetPr(T4),(char*)f4 ,100*sizeof(double));
	engPutVariable(ep , "TX1", T3) ;
	engPutVariable(ep , "TY1", T4) ;
	engPutVariable(ep , "NT1", M1) ;
	if(flag2==false)
	{
	flag2 =true ;
	engEvalString( ep ,"h=plot(TX1,TY1,'ro');grid on");

	engEvalString(ep,"m=0:0.01:1,frontx=m.^1.25 ;fronty=(1-m).^1.25 ;"); 
	engEvalString(ep,"hold on");
	engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");		
	engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi/10)) ;fronty=(1-m).^(1.25+0.75*sin(pi/10)) ;"); 
	engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");		
	engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*2/10)) ;fronty=(1-m).^(1.25+0.75*sin(pi*2/10)) ;");  
	engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
	engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*3/10)) ;fronty=(1-m).^(1.25+0.75*sin(pi*3/10)) ;"); 
	engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
	engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*4/10)) ;fronty=(1-m).^(1.25+0.75*sin(pi*4/10)) ;"); 
	engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
	engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*5/10)) ;fronty=(1-m).^(1.25+0.75*sin(pi*5/10)) ;");
	engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
	engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*11/10)) ;fronty=(1-m).^(1.25+0.75*sin(pi*11/10)) ;");
	engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
	engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*12/10)) ;fronty=(1-m).^(1.25+0.75*sin(pi*12/10)) ;");
	engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
	engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*13/10)) ;fronty=(1-m).^(1.25+0.75*sin(pi*13/10)) ;");
	engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
	engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*14/10)) ;fronty=(1-m).^(1.25+0.75*sin(pi*14/10)) ;");
	engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
	engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*15/10)) ;fronty=(1-m).^(1.25+0.75*sin(pi*15/10)) ;");
	engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");						   
	}
	else
	{
	engEvalString( ep ,"set(h, 'XData',TX1     ,'YData',TY1   )");
	}
	}
	engEvalString(ep,"hold off");
	engClose(ep);*/

	double temp_minX[20];
	double temp_maxX[20];
	for(int i=0;i<xdim;i++)
	{
		temp_minX[i]=Population()[0][i];
		temp_maxX[i]=Population()[0][i];
		for (int n = 1 ; n < psize ; n++)
		{
			if (temp_minX[i]>Population()[n][i])
				temp_minX[i] = Population()[n][i];
			if (temp_maxX[i]<Population()[n][i])
				temp_maxX[i]=Population()[n][i];
		}
	}
	CPopulationMO	init_pop(P());
	init_pop = Population();

	if (mStep == mTaoT)
	{
		for(int i=0; i<psize; i++) 
		{
			for(int j=0; j<P().XSize(); j++) 
			{
				if ((temp_minX[j]-0.05>=P().XLow(j))&&(temp_maxX[j]+0.05<=P().XUpp(j)))
					init_pop[i][j] = rnd::rand(temp_minX[j]-0.05,temp_maxX[j]+0.05);
				if ((temp_minX[j]-0.05>=P().XLow(j))&&(temp_maxX[j]+0.05 > P().XUpp(j)))
					init_pop[i][j] = rnd::rand(temp_minX[j]-0.05,P().XUpp(j));
				if ((temp_minX[j]-0.05 <P().XLow(j))&&(temp_maxX[j]+0.05<=P().XUpp(j)))
					init_pop[i][j] = rnd::rand(P().XLow(j),temp_maxX[j]+0.05);
				if ((temp_minX[j]-0.05<P().XLow(j))&&(temp_maxX[j]+0.05>P().XUpp(j)))
					init_pop[i][j] = rnd::rand(P().XLow(j), P().XUpp(j));
			}
		}
		oldC=mC;
	}
	else
	{
		for (int n = 0 ; n < 100 ; n++)
		{
			for (int i=0;i<xdim;i++)
			{
				double x_random = abs((mC[i]-oldC[i])*rnd::gaussian());
				if (mC[i]-oldC[i]>0)
					init_pop[n][i] += x_random;
				if (mC[i]-oldC[i]<0)
					init_pop[n][i] -= x_random;
				if (init_pop[n][i]>P().XUpp(i))
				{
					init_pop[n][i] = P().XUpp(i) - (init_pop[n][i] - P().XUpp(i))/2;
				}
				if (init_pop[n][i] < P().XLow(i))
				{
					init_pop[n][i] = P().XLow(i) + (P().XLow(i) - init_pop[n][i])/2;
				}
			}
		}
		oldC = mC;
	}

	init_pop.Evaluate();
	

	//static bool flag =false ;
	//double f1[100];
	//double f2[100];
	//mxArray *T = mxCreateDoubleMatrix( 1 , 100 , mxREAL  ) ; 
	//mxArray *T2 = mxCreateDoubleMatrix( 1 ,  100 , mxREAL  ) ;
	//mxArray *M = mxCreateDoubleMatrix( 1 ,  1 , mxREAL  ) ;
	//for(int n=0; n<100; n++)
	//{ 
	//	f1[n]=init_pop[n].X(0);
	//	f2[n]=init_pop[n].X(1);
	//}
	//if(  (ep=engOpen(NULL)) )
	//{
	//	memcpy( (char*)mxGetPr(T),(char*)f1 , 100*sizeof(double));
	//	memcpy( (char*)mxGetPr(T2),(char*)f2 ,100*sizeof(double));
	//	engPutVariable(ep , "TX", T) ;
	//	engPutVariable(ep , "TY", T2) ;
	//	engPutVariable(ep , "NT", M) ;
	//	if(flag==false)
	//	{
	//		flag =true ;
	//		engEvalString( ep ,"h=plot(TX,TY,'ro');grid on");

	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*0/10);"); 
	//		engEvalString(ep,"hold on");
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");	
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*1/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");	
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*2/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");	
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*3/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*4/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*5/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");

	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*11/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*12/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*13/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*14/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*15/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");

	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*21/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*22/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*23/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*24/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*25/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");

	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*31/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*32/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*33/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*34/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");
	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*35/10);"); 
	//		engEvalString(ep,"plot(frontx,fronty,'b.')");	

	//	}
	//	else
	//	{
	//		engEvalString( ep ,"set(h, 'XData',TX     ,'YData',TY   )");
	//	}
	//}
	//engEvalString(ep,"hold off");
	//engClose(ep);

	for (int i = 0; i < 2; i++)
		for(int j = 0; j < psize; j++)
			init_pop[j].new_GD(i) = (int)((init_pop[j].F(i)-temp_minF[i])/((temp_maxF[i]-temp_minF[i])/der));
	int sum=0;bool fl=false;
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	/*for (int j = 50;j < 100; j++)
	{
		int maxfit=-1,m=0;
		fl=false;
		for (int k = 0;k < psize; k++)
		{
			if (init_pop[k].new_GD(0)<=Population()[j].new_GD(0)&&init_pop[k].new_GD(1)<=Population()[j].new_GD(1))
			{
				if (maxfit<abs(init_pop[k].new_GD(0)-Population()[j].new_GD(0))+abs(init_pop[k].new_GD(1)-Population()[j].new_GD(1)))
				{
					maxfit = abs(init_pop[k].new_GD(0)-Population()[j].new_GD(0))+abs(init_pop[k].new_GD(1)-Population()[j].new_GD(1));
					m=k;
				}
				if (maxfit==abs(init_pop[k].new_GD(0)-Population()[j].new_GD(0))+abs(init_pop[k].new_GD(1)-Population()[j].new_GD(1)))
				{
					int q = (init_pop[m].F(0)-temp_minF[0])/((temp_maxF[0]-temp_minF[0])/40)-init_pop[m].new_GD(0)+(init_pop[m].F(1)-temp_minF[1])/((temp_maxF[1]-temp_minF[1])/40)-init_pop[m].new_GD(1);
					int p = (init_pop[k].F(0)-temp_minF[0])/((temp_maxF[0]-temp_minF[0])/40)-init_pop[k].new_GD(0)+(init_pop[k].F(1)-temp_minF[1])/((temp_maxF[1]-temp_minF[1])/40)-init_pop[k].new_GD(1);
					if (q>p)
						m=k;
				}
				fl = true;
			}
		}
		if ( fl == true)
		sum++;
		for (int i = 0;i<xdim;i++)
			Population()[j][i] +=rnd::rand(0,1)*(init_pop[m].X(i)-Population()[j][i]);
	}
	printf("%d\n",sum);*/
	
	//Engine *ep1;
	//static bool flag3 =false ;
	//double f5[30];
	//double f6[30];
	//mxArray *T5 = mxCreateDoubleMatrix( 1 , 30 , mxREAL  ) ; 
	//mxArray *T6 = mxCreateDoubleMatrix( 1 ,  30 , mxREAL  ) ;
	//for(int n=40; n<70; n++)
	//{ 
	//	f5[n-40]=Population()[n].F(0);
	//	f6[n-40]=Population()[n].F(1);
	//}
	//if(  (ep1=engOpen(NULL)) )
	//{
	//	memcpy( (char*)mxGetPr(T5),(char*)f5 , 30*sizeof(double));
	//	memcpy( (char*)mxGetPr(T6),(char*)f6 ,30*sizeof(double));
	//	engPutVariable(ep1 , "TX2", T5) ;
	//	engPutVariable(ep1 , "TY2", T6) ;
	//	//engPutVariable(ep , "NT1", M1) ;
	//	if(flag3==false)
	//	{
	//		flag3 =true ;
	//		engEvalString( ep1 ,"h=plot(TX2,TY2,'ro');grid on");
	//		engEvalString(ep1,"frontx=0:0.01:1 ;fronty=1-frontx.^0.5 ;"); 
	//		engEvalString(ep1,"hold on");
	//		engEvalString(ep1,"plot(frontx,fronty,'Linewidth',2)");						   
	//	}
	//	else
	//	{
	//		engEvalString( ep1 ,"set(h, 'XData',TX2     ,'YData',TY1   )");
	//	}
	//}
	//engEvalString(ep1,"hold off");
	//engClose(ep1);

	int mincha[2]={der*der,der*der};
	for (int i = numsub1; i<psize; i++)
	{
		for (int j =i+1; j<psize; j++)
		{
			if (mincha[0]>abs(Population()[j].new_GD(0)-Population()[i].new_GD(0))&&Population()[j].new_GD(0)!=Population()[i].new_GD(0))
			{
				mincha[0] = abs(Population()[j].new_GD(0)-Population()[i].new_GD(0));
			}
			if (mincha[1]>abs(Population()[j].new_GD(1)-Population()[i].new_GD(1))&&Population()[j].new_GD(1)!=Population()[i].new_GD(1))
			{
				mincha[1] = abs(Population()[j].new_GD(1)-Population()[i].new_GD(1));
			}
		}
	}

	for (int i = 99; i >= numsub1+numsub2; i--)
	{
		for (int j = numsub1; j< i; j++)
		{
			if ((abs(Population()[i].new_GD(0)-Population()[j].new_GD(0))>=(der/30+1)*mincha[0])||(abs(Population()[i].new_GD(1)-Population()[j].new_GD(1))>=(der/30+1)*mincha[1]))
			{
				if (i-1!=j) sum++;
				CPopulationMO	swapindiv(P());
				swapindiv = Population();
				swapindiv[0]=Population()[i-1];
				Population()[i-1]=Population()[j];
				Population()[j]=swapindiv[0];
				break;
			}
			else
			{
				if ((abs(Population()[i].new_GD(0)-Population()[j].new_GD(0))<=mincha[0])&&(abs(Population()[i].new_GD(1)-Population()[j].new_GD(1))<=mincha[1]))
				{
					for (int n = j+1; n < i;n++)
					{
						if (Population()[n].new_GD(0)==Population()[j].new_GD(0)&&Population()[n].new_GD(1)==Population()[j].new_GD(1))
						{
							if (i-1!=j) 
								sum++;
							CPopulationMO	swapindiv(P());
							swapindiv = Population();
							swapindiv[0]=Population()[i-1];
							Population()[i-1]=Population()[j];
							Population()[j]=swapindiv[0];
							j = psize;
							break;
						}
					}
				}
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	/*for (int i = 40; i < 70; i++)
	{
		int flag = false;
		for (int j = i+1; j< psize; j++)
		{
			if ((abs(Population()[i].new_GD(0)-Population()[j].new_GD(0))>mincha[0]+der/30)||(abs(Population()[i].new_GD(1)-Population()[j].new_GD(1))>mincha[1]+der/30))
			{
				if (i+1!=j) sum++;
				flag = true;
				CPopulationMO	swapindiv(P());
				swapindiv = Population();
				swapindiv[0]=Population()[i+1];
				Population()[i+1]=Population()[j];
				Population()[j]=swapindiv[0];
				break;
			}
		}
		if (flag == false)
		{
			for (int j = i+1; j< psize; j++)
			{
				if ((abs(Population()[i].new_GD(0)-Population()[j].new_GD(0))<=mincha[0]+der/30)&&(abs(Population()[i].new_GD(1)-Population()[j].new_GD(1))<=mincha[1]+der/30))
				{
					for (int n = j+1; n < psize;n++)
					{
						if (abs(Population()[n].new_GD(0)-Population()[j].new_GD(0))<=mincha[0]&&abs(Population()[n].new_GD(1)-Population()[j].new_GD(1))<=mincha[1])
						{
								sum++;
							CPopulationMO	swapindiv(P());
							swapindiv = Population();
							swapindiv[0]=Population()[i+1];
							Population()[i+1]=Population()[j];
							Population()[j]=swapindiv[0];
							j = psize;
							break;
						}
					}
				}
			}
		}
	}*/
	//printf("%d\n",sum);
	/*Engine *ep;
	static bool flag2 =false ;
	double f3[30];
	double f4[30];
	mxArray *T3 = mxCreateDoubleMatrix( 1 , 30 , mxREAL  ) ; 
	mxArray *T4 = mxCreateDoubleMatrix( 1 ,  30 , mxREAL  ) ;
	for(int n=40; n<70; n++)
	{ 
	f3[n-40]=Population()[n].F(0);
	f4[n-40]=Population()[n].F(1);
	}
	if(  (ep=engOpen(NULL)) )
	{
	memcpy( (char*)mxGetPr(T3),(char*)f3 , 30*sizeof(double));
	memcpy( (char*)mxGetPr(T4),(char*)f4 ,30*sizeof(double));
	engPutVariable(ep , "TX1", T3) ;
	engPutVariable(ep , "TY1", T4) ;
	engPutVariable(ep , "NT1", M1) ;
	if(flag2==false)
	{
	flag2 =true ;
	engEvalString( ep ,"h=plot(TX1,TY1,'ro');grid on");
	engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^0.5 ;"); 
	engEvalString(ep,"hold on");
	engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");						   
	}
	else
	{
	engEvalString( ep ,"set(h, 'XData',TX1     ,'YData',TY1   )");
	}
	}
	engEvalString(ep,"hold off");
	engClose(ep);*/
	//for(int n=70; n<100; n++)
	//{ 
	//	f3[n-70]=Population()[n].F(0);
	//	f4[n-70]=Population()[n].F(1);
	//}
	//if(  (ep=engOpen(NULL)) )
	//{
	//	memcpy( (char*)mxGetPr(T3),(char*)f3 , 30*sizeof(double));
	//	memcpy( (char*)mxGetPr(T4),(char*)f4 ,30*sizeof(double));
	//	engPutVariable(ep , "TX1", T3) ;
	//	engPutVariable(ep , "TY1", T4) ;
	//	//engPutVariable(ep , "NT1", M1) ;
	//	if(flag2==false)
	//	{
	//		flag2 =true ;
	//		engEvalString( ep ,"h=plot(TX1,TY1,'ro');grid on");

	//		engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^0.5 ;"); 
	//		engEvalString(ep,"hold on");
	//		engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");						   
	//	}
	//	else
	//	{
	//		engEvalString( ep ,"set(h, 'XData',TX1     ,'YData',TY1   )");
	//	}
	//}
	//engEvalString(ep,"hold off");
	//engClose(ep);
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	for (int j = numsub1;j < numsub1+numsub2; j++)
	{
		int maxfit=65535,m=0;
		for (int k = 0;k < psize; k++)
		{
			if (((init_pop[k].new_GD(0)<=Population()[j].new_GD(0)&&init_pop[k].new_GD(0)>=Population()[j].old_GD(0))||
				(init_pop[k].new_GD(0)>=Population()[j].new_GD(0)&&init_pop[k].new_GD(0)<=Population()[j].old_GD(0)))&&
				((init_pop[k].new_GD(1)<=Population()[j].new_GD(1)&&init_pop[k].new_GD(1)>=Population()[j].old_GD(1))||
				(init_pop[k].new_GD(1)>=Population()[j].new_GD(1)&&init_pop[k].new_GD(1)<=Population()[j].old_GD(1))))
			{
				if (maxfit>abs(init_pop[k].new_GD(0)-Population()[j].old_GD(0))+abs(init_pop[k].new_GD(1)-Population()[j].old_GD(1)))
				{
					maxfit = abs(init_pop[k].new_GD(0)-Population()[j].old_GD(0))+abs(init_pop[k].new_GD(1)-Population()[j].old_GD(1));
					m=k;
				}
				if (maxfit==abs(init_pop[k].new_GD(0)-Population()[j].old_GD(0))+abs(init_pop[k].new_GD(1)-Population()[j].old_GD(1)))
				{
					int q = (init_pop[m].F(0)-temp_minF[0])/((temp_maxF[0]-temp_minF[0])/der)-init_pop[m].new_GD(0)+(init_pop[m].F(1)-temp_minF[1])/((temp_maxF[1]-temp_minF[1])/der)-init_pop[m].new_GD(1);
					int p = (init_pop[k].F(0)-temp_minF[0])/((temp_maxF[0]-temp_minF[0])/der)-init_pop[k].new_GD(0)+(init_pop[k].F(1)-temp_minF[1])/((temp_maxF[1]-temp_minF[1])/der)-init_pop[k].new_GD(1);
					if (q>p)
						m=k;
				}
			}
		}
		if (maxfit != 65535)
		{
			for (int i = 0;i<xdim;i++)
				Population()[j][i] +=rnd::rand(0.9,1.0)*(init_pop[m].X(i)-Population()[j][i]);
		}
		//else
		//	for (int i = 0;i<xdim;i++)
		//		Population()[j][i] +=rnd::rand(0,1)*(Population()[rnd::rand(0,99)][i]-Population()[j][i]);
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	for (int i = numsub1+numsub2;i < 100; i++)
	{
		//int a[2];
		/*a[0]=Population()[i].new_GD(0);
		a[1]=Population()[i].new_GD(1);
		if (Population()[i].new_GD(0)-Population()[i].old_GD(0)>0)
		a[0]+=1;
		if (Population()[i].new_GD(0)-Population()[i].old_GD(0)<0)
		a[0]-=1;
		if (Population()[i].new_GD(1)-Population()[i].old_GD(1)>0)
		a[1]+=1;
		if (Population()[i].new_GD(1)-Population()[i].old_GD(1)<0)
		a[1]-=1;
		int m=-1;
		for (int j =0;j<psize;j++)
		{
		if (init_pop[j].new_GD(0)==a[0]&&init_pop[j].new_GD(1)==a[1])
		{
		m=j;
		break;
		}
		}
		if (m != -1)
		for (int j = 0;j<xdim;j++)
		Population()[i][j] +=1*(init_pop[m].X(j)-Population()[i][j]);*/
		int minfit=-1,m=0;
		if (Population()[i].new_GD(0)>=Population()[i].old_GD(0)&&Population()[i].new_GD(1)>=Population()[i].old_GD(1))
		{
			for (int j =0;j<psize;j++)
			{
				if (init_pop[j].new_GD(0)>=Population()[i].old_GD(0)&&init_pop[j].new_GD(1)>=Population()[i].old_GD(1))
				{
					if (minfit<abs(init_pop[j].new_GD(0)-Population()[i].old_GD(0))+abs(init_pop[j].new_GD(1)-Population()[i].old_GD(1)))
					{
						minfit = abs(init_pop[j].new_GD(0)-Population()[i].old_GD(0))+abs(init_pop[j].new_GD(1)-Population()[i].old_GD(1));
						m=j;
					}
					if (minfit==abs(init_pop[j].new_GD(0)-Population()[i].old_GD(0))+abs(init_pop[j].new_GD(1)-Population()[i].old_GD(1)))
					{
						int q = (init_pop[m].F(0)-temp_minF[0])/((temp_maxF[0]-temp_minF[0])/der)-init_pop[m].new_GD(0)+(init_pop[m].F(1)-temp_minF[1])/((temp_maxF[1]-temp_minF[1])/der)-init_pop[m].new_GD(1);
						int p = (init_pop[j].F(0)-temp_minF[0])/((temp_maxF[0]-temp_minF[0])/der)-init_pop[j].new_GD(0)+(init_pop[j].F(1)-temp_minF[1])/((temp_maxF[1]-temp_minF[1])/der)-init_pop[j].new_GD(1);
						if (q>p)
							m=j;
					}
				}
			}
		}

		if (Population()[i].new_GD(0)<=Population()[i].old_GD(0)&&Population()[i].new_GD(1)<=Population()[i].old_GD(1))
		{
			for (int j =0;j<psize;j++)
			{
				if (init_pop[j].new_GD(0)<=Population()[i].old_GD(0)&&init_pop[j].new_GD(1)<=Population()[i].old_GD(1))
				{
					if (minfit<abs(init_pop[j].new_GD(0)-Population()[i].old_GD(0))+abs(init_pop[j].new_GD(1)-Population()[i].old_GD(1)))
					{
						minfit = abs(init_pop[j].new_GD(0)-Population()[i].old_GD(0))+abs(init_pop[j].new_GD(1)-Population()[i].old_GD(1));
						m=j;
					}
					if (minfit==abs(init_pop[j].new_GD(0)-Population()[i].old_GD(0))+abs(init_pop[j].new_GD(1)-Population()[i].old_GD(1)))
					{
						int q = (init_pop[m].F(0)-temp_minF[0])/((temp_maxF[0]-temp_minF[0])/der)-init_pop[m].new_GD(0)+(init_pop[m].F(1)-temp_minF[1])/((temp_maxF[1]-temp_minF[1])/der)-init_pop[m].new_GD(1);
						int p = (init_pop[j].F(0)-temp_minF[0])/((temp_maxF[0]-temp_minF[0])/der)-init_pop[j].new_GD(0)+(init_pop[j].F(1)-temp_minF[1])/((temp_maxF[1]-temp_minF[1])/der)-init_pop[j].new_GD(1);
						if (q>p)
							m=j;
					}
				}
			}
		}

		if (Population()[i].new_GD(0)>=Population()[i].old_GD(0)&&Population()[i].new_GD(1)<=Population()[i].old_GD(1))
		{
			for (int j =0;j<psize;j++)
			{
				if (init_pop[j].new_GD(0)>=Population()[i].old_GD(0)&&init_pop[j].new_GD(1)<=Population()[i].old_GD(1))
				{
					if (minfit<abs(init_pop[j].new_GD(0)-Population()[i].old_GD(0))+abs(init_pop[j].new_GD(1)-Population()[i].old_GD(1)))
					{
						minfit = abs(init_pop[j].new_GD(0)-Population()[i].old_GD(0))+abs(init_pop[j].new_GD(1)-Population()[i].old_GD(1));
						m=j;
					}
					if (minfit==abs(init_pop[j].new_GD(0)-Population()[i].old_GD(0))+abs(init_pop[j].new_GD(1)-Population()[i].old_GD(1)))
					{
						int q = (init_pop[m].F(0)-temp_minF[0])/((temp_maxF[0]-temp_minF[0])/der)-init_pop[m].new_GD(0)+(init_pop[m].F(1)-temp_minF[1])/((temp_maxF[1]-temp_minF[1])/der)-init_pop[m].new_GD(1);
						int p = (init_pop[j].F(0)-temp_minF[0])/((temp_maxF[0]-temp_minF[0])/der)-init_pop[j].new_GD(0)+(init_pop[j].F(1)-temp_minF[1])/((temp_maxF[1]-temp_minF[1])/der)-init_pop[j].new_GD(1);
						if (q>p)
							m=j;
					}
				}
			}
		}
		if (Population()[i].new_GD(0)<=Population()[i].old_GD(0)&&Population()[i].new_GD(1)>=Population()[i].old_GD(1))
		{
			for (int j =0;j<psize;j++)
			{
				if (init_pop[j].new_GD(0)>=Population()[i].old_GD(0)&&init_pop[j].new_GD(1)>=Population()[i].old_GD(1))
				{
					if (minfit<abs(init_pop[j].new_GD(0)-Population()[i].old_GD(0))+abs(init_pop[j].new_GD(1)-Population()[i].old_GD(1)))
					{
						minfit = abs(init_pop[j].new_GD(0)-Population()[i].old_GD(0))+abs(init_pop[j].new_GD(1)-Population()[i].old_GD(1));
						m=j;
					}
					if (minfit==abs(init_pop[j].new_GD(0)-Population()[i].old_GD(0))+abs(init_pop[j].new_GD(1)-Population()[i].old_GD(1)))
					{
						int q = (init_pop[m].F(0)-temp_minF[0])/((temp_maxF[0]-temp_minF[0])/der)-init_pop[m].new_GD(0)+(init_pop[m].F(1)-temp_minF[1])/((temp_maxF[1]-temp_minF[1])/der)-init_pop[m].new_GD(1);
						int p = (init_pop[j].F(0)-temp_minF[0])/((temp_maxF[0]-temp_minF[0])/der)-init_pop[j].new_GD(0)+(init_pop[j].F(1)-temp_minF[1])/((temp_maxF[1]-temp_minF[1])/der)-init_pop[j].new_GD(1);
						if (q>p)
							m=j;
					}
				}
			}
		}
		if (minfit != -1)
			for (int j = 0;j<xdim;j++)
				Population()[i][j] +=rnd::rand(0.9,1.0)*(init_pop[m].X(j)-Population()[i][j]);

		//else
		//	for (int j = 0;j<xdim;j++)
		//		Population()[i][j] +=rnd::rand(0,1)*(Population()[rnd::rand(0,99)][j]-Population()[i][j]);
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	Population().Evaluate();
	int num_do = 0; float sub1=0.0 ,sub2=0.0, sub3=0.0;
	CPopulationMO	fuben(P());
	fuben = Population();
	Population().IsSort(false);
	Population().RankSort(false);
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	for (int i=0; i<psize; i++)
	{
		if (Population()[i].Rank()==1)
			num_do++;
	}
	for (int j = 0 ; j < numsub1 ; j++)
	{
		for (int i = 0; i < num_do; i++)
		{
			if (Population()[i]==fuben[j])
			{
				sub1 = sub1 + 1.0;
				break;
			}
		}
	}

	for (int j = numsub1 ; j < numsub1+numsub2 ; j++)
	{
		for (int i = 0; i < num_do; i++)
		{
			if (Population()[i]==fuben[j])
			{
				sub2 = sub2 + 1.0;
				break;
			}
		}
	}

	for (int j = numsub1+numsub2 ; j < psize ; j++)
	{
		for (int i = 0; i < num_do; i++)
		{
			if (Population()[i]==fuben[j])
			{
				sub3 = sub3 + 1.0;
				break;
			}
		}
	}

	sub1 = sub1/numsub1;
	sub2 = sub2/numsub2;
	sub3 = sub3/numsub3;

	float a = sub1>sub2?sub1:sub2;
	float b = a>sub3?a:sub3;
	if (b == sub1)
	{
		if (numsub2 > 10 && numsub3 > 10 )
		{
		numsub1 += numsub2/5 + numsub3/5;
		numsub2 = numsub2 - numsub2/5;
		numsub3 = numsub3 - numsub3/5;
		}
		else if (numsub2 > 10)
		{
			numsub1 += numsub2/5;
			numsub2 = numsub2 - numsub2/5;
		}
		else if (numsub3 > 10)
		{
			numsub1 += numsub3/5;
			numsub3 = numsub3 - numsub3/5;
		}
	}
	if (b == sub2)
	{
		if (numsub3 > 10 && numsub1 > 10)
		{
		numsub2 += numsub1/5 + numsub3/5;
		numsub1 = numsub1 - numsub1/5;
		numsub3 = numsub3 - numsub3/5;
		}
		else if (numsub1 > 10)
		{
			numsub2 += numsub1/5;
			numsub1 = numsub1 - numsub1/5;
		}
		else if (numsub3 > 10)
		{
			numsub2 += numsub3/5;
			numsub3 = numsub3 - numsub3/5;
		}
	}
	if (b == sub3)
	{
		if (numsub1 > 10 && numsub2 > 10)
		{
		numsub3 += numsub1/5 + numsub2/5;
		numsub1 = numsub1 - numsub1/5;
		numsub2 = numsub2 - numsub2/5;
		}
		else if (numsub1 > 10)
		{
			numsub3 += numsub1/5;
			numsub1 = numsub1 - numsub1/5;
		}
		else if (numsub2 > 10)
		{
			numsub3 += numsub2/5;
			numsub2 = numsub2 - numsub2/5;
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	memory.Evaluate();
	CPopulationMO  copy_memory(memory);
	copy_memory = memory;
	memory.IsSort(false);
	memory.RankSort(false);
	for (int i = 0; i < psize/20; i++)
	{
		if (memory[i].Dominate(Population()[psize-i-1]))
		{
			for (int k = 0 ;k < xdim ;k++)
				Population()[psize-i-1][k] = memory[i][k];
		}
	}
	memory = copy_memory;
	/***
		基于预测和记忆相结合的双重策略  （+Observe()）
	**/
//	CPopulationMO	old_pop(P());
//	old_pop = Population();
////	unsigned int i, k, xdim=P().XSize(), psize = Population().Size();
//	//Step 1: save center point and shape points
//	mC.resize(xdim);
//	int num1=0,num_nodominance=0;
//	for (i=0; i<psize; i++)
//	{
//		if (Population()[i].Rank()==1)
//			num1++;
//	}
//	num_nodominance=num1;
//	if (num1<psize/20)
//	{
//		for (i=0;i<psize/20;i++)
//			memory.Copy(Population()[i]);
//	}	
//	else
//	{
//		std::vector<size_t> vnIndex(num1);
//		for (size_t j = 0; j < vnIndex.size(); ++j)
//			vnIndex[j] = j;
//		std::random_shuffle(vnIndex.begin(), vnIndex.end());
//		for (size_t j = 0; j < psize/20; ++j)
//			memory.Copy(Population()[vnIndex[j]]);
//	}
//	if(memory.Size()>2*psize) 
//	{
//		memory.Pop(memory.Size()-2*psize);
//	}
//	for(k=0; k<xdim; k++) 
//	{
//		mC[k]  = 0.0;
//		for(i=0; i<num1; i++)
//			if (Population()[i].Rank()==1)
//			    mC[k] += Population()[i][k];
//		mC[k] /= double(num1);
//	}

	Population().Evaluate(); mEvas += Population().Size();
}

void DMOO::InitGIPS()
{
		//Step 0: initialize the population
	if(mStep == 0) 
	{
		InitRIS(true); 
		return;
	}
	CPopulationMO	search_pop(P()),old_pop(P());
	search_pop = Population(),old_pop = Population();

	unsigned int i, k, xdim=P().XSize(), psize = Population().Size();

	//Step 1: save center point and shape points
	oldC.resize(xdim);
	mC.resize(xdim);
	newC.resize(xdim);
	D.resize(xdim);

	int num_nodominance=0;
	for (i=0; i<psize; i++)
	{
		if (Population()[i].Rank()==1)
			num_nodominance++;
	}
	if(mStep <= mTaoT)
	{
		for(k=0; k<xdim; k++) 
		{
			oldC[k]  = 0.0;
			mC[k] = 0.0;
			for(i=0; i<psize; i++)
				if (Population()[i].Rank()==1)
					mC[k] += Population()[i][k];
			mC[k] /= double(num_nodominance);
		}
	}
	else
	{
		for (k = 0; k <xdim; k++)
		{
			oldC[k] = mC[k];
			mC[k] = 0.0;
			for(i=0; i<psize; i++)
				if (Population()[i].Rank()==1)
					mC[k] += Population()[i][k];
			mC[k] /= double(num_nodominance);
		}
	}

	//判断种群进化方向的判断算子
	for (int n=0;n<2;n++)
	{
		CPopulationMO popnew(P());

		// generate new solutions using RM2
		az::mea::gen::mod::RM2 gen;
		gen.Set(1.0, 1.0, 0.8, 1,30);
		gen.Generate(mPopSize, popnew, Population());

		//check the range of new solutions
		Check(popnew);
		//remove these repeated solutions
		Check(popnew, Population());

		//evaluate new solutions
		popnew.Evaluate(); mEvas += popnew.Size();

		//environmental select
		Population().Combine(popnew);

		az::mea::sel::SCrowd2 sel;
		sel.Select(Population(), mPopSize);
	}

	int num2=0;
	for(int k=0; k<xdim; k++) 
	{
		num2=0;
		newC[k]  = 0.0;
		for(int i=0; i<psize; i++)
			if (Population()[i].Rank()==1)
			{
				newC[k] += Population()[i][k];
				num2++;
			}
			newC[k] /= double(num2);
	}

	Population() = old_pop;
	observe = Population();

	int num1=0;
	for (int i=0; i<psize; i++)
	{
		if (Population()[i].Rank()==1)
			num1++;
	}

	az::mea::sel::XCrowd selob;
	selob.Select(observe, 10);
	for (int n = 0 ; n < 10 ; n++)
	{
		for(int m = 0; m < 8; m++)
		{
			for(int i=0;i<xdim;i++)
			{
				if(newC[i]-mC[i]>0)
				{
					Population()[n*10+m+1][i] =observe[n][i] + (m+1)*(P().XUpp(i)-mC[i])/10.0;
				}
				if(newC[i]-mC[i]<0)
				{
					Population()[n*10+m+1][i] =observe[n][i] - (m+1)*(mC[i]-P().XLow(i))/10.0;
				}
				if (Population()[n*10+m+1][i]>P().XUpp(i))
				{
					Population()[n*10+m+1][i] = P().XUpp(i) - (P().XUpp(i)-mC[i])/50.0;
				}
				if (Population()[n*10+m+1][i]<P().XLow(i))
				{
					Population()[n*10+m+1][i] = P().XLow(i) + (mC[i]-P().XLow(i))/50.0;
				}
			}
			observe.Copy(Population()[n*10+m+1]);
		}
	}
	if (mStep == mTaoT)
	{
		oldC=mC;
	}
	else
	{
		for (int n = 0 ; n < 10 ; n++)
		{
			for (int i=0;i<xdim;i++)
			{
				Population()[0][i] = observe[n][i]+mC[i]-oldC[i];
				if (Population()[0][i]>P().XUpp(i))
				{
					Population()[0][i] = P().XUpp(i) - (P().XUpp(i)-mC[i])/10.0;
				}
				if (Population()[0][i]<P().XLow(i))
				{
					Population()[0][i] = P().XLow(i) + (mC[i]-P().XLow(i))/10.0;
				}
			}
			observe.Copy(Population()[0]);
		}
		oldC = mC;
	}
	Population() = old_pop;
	Check(observe);
	observe.Evaluate();mEvas += observe.Size();
	observe.IsSort(false);
	observe.RankSort(false);
	az::mea::sel::SCrowd2 sel2;
	sel2.Select(observe, 10);
	CPopulationMO newpop(Population());

	for (int k = 0 ;k < xdim ;k++)
	{
		for (int i = 0; i < observe.Size(); i++)
		{
			for(int j=0;j<5;j++)
			{
				if (j==0)
			       newpop[i*psize/observe.Size()][k] = observe[i][k];
				else
				{
					newpop[i*psize/observe.Size()+j*2][k] = observe[i][k]+rnd::rand(-0.05,0.05);
					if(newpop[i*psize/observe.Size()+j*2][k]>=P().XUpp(k)) 
						newpop[i*psize/observe.Size()+j*2][k] = P().XUpp(k)-rnd::rand(0.0,0.02);
					else if(newpop[i*psize/observe.Size()+j*2][k]<=P().XLow(k)) 
						newpop[i*psize/observe.Size()+j*2][k] = P().XLow(k)+rnd::rand(0.0,0.02);
				}
			}
		}
	}
	Check(newpop);
	newpop.Evaluate();mEvas += newpop.Size();
	newpop.IsSort(false);
	newpop.RankSort(false);
	az::mea::sel::SCrowd2 sel3;
	sel3.Select(newpop, 20);

	for (int i = 0; i < newpop.Size(); i++)
	{
		for (int k = 0 ;k < xdim ;k++)
		{
			Population()[psize-i-1][k] = newpop[i][k];
		}
	}
}

void DMOO::InitPMS()
{	//Population().Evaluate(); mEvas += Population().Size();
	memory.Evaluate();
    CPopulationMO	old_pop(P());
	CPopulationMO  copy_memory(memory);
	copy_memory = memory;
	memory.IsSort(false);
	memory.RankSort(false);
    old_pop = Population();
	unsigned int xdim=P().XSize(), psize = Population().Size(); 
	for (int n=0;n<2;n++)
	{
	//Population().Evaluate(); mEvas += Population().Size();
	//CPopulationMO	old_pop(P());
	//old_pop = Population();

	CPopulationMO popnew(P());

	// generate new solutions
	Generate(popnew, mPopSize);

	//check the range of new solutions
	Check(popnew);
	//remove these repeated solutions
	Check(popnew, Population());

	//evaluate new solutions
	popnew.Evaluate(); mEvas += popnew.Size();

	//environmental select
	Population().Combine(popnew);

	az::mea::sel::SCrowd2 sel;
	sel.Select(Population(), mPopSize);
	}

	int num2=0;
	newC.resize(xdim);
	for(int k=0; k<xdim; k++) 
	{
		num2=0;
		newC[k]  = 0.0;
		for(int i=0; i<psize; i++)
			if (Population()[i].Rank()==1)
			{
				newC[k] += Population()[i][k];
				num2++;
			}
			newC[k] /= double(num2);
	}

	Population() = old_pop;
	observe = Population();




	int num1=0;
	for (int i=0; i<psize; i++)
	{
		if (Population()[i].Rank()==1)
			num1++;
	}

	az::mea::sel::XCrowd selob;
	selob.Select(observe, 10);
	for (int n = 0 ; n < 10 ; n++)
	{
		for(int m = 0; m < 8; m++)
		{
			for(int i=0;i<xdim;i++)
			{
				if(newC[i]-mC[i]>0)
				{
					Population()[n*10+m+1][i] =observe[n][i] + (m+1)*(P().XUpp(i)-mC[i])/10.0;
				}
				if(newC[i]-mC[i]<0)
				{
					Population()[n*10+m+1][i] =observe[n][i] - (m+1)*(mC[i]-P().XLow(i))/10.0;
				}
				if (Population()[n*10+m+1][i]>P().XUpp(i))
				{
					Population()[n*10+m+1][i] = P().XUpp(i) - (P().XUpp(i)-mC[i])/50.0;
				}
				if (Population()[n*10+m+1][i]<P().XLow(i))
				{
					Population()[n*10+m+1][i] = P().XLow(i) + (mC[i]-P().XLow(i))/50.0;
				}
			}
			observe.Copy(Population()[n*10+m+1]);
		}
	}
	if (mStep == mTaoT)
	{
		oldC=mC;
	}
	else
	{
		for (int n = 0 ; n < 10 ; n++)
		{
			for (int i=0;i<xdim;i++)
			{
				Population()[0][i] = observe[n][i]+mC[i]-oldC[i];
				if (Population()[0][i]>P().XUpp(i))
				{
					Population()[0][i] = P().XUpp(i) - (P().XUpp(i)-mC[i])/10.0;
				}
				if (Population()[0][i]<P().XLow(i))
				{
					Population()[0][i] = P().XLow(i) + (mC[i]-P().XLow(i))/10.0;
				}
			}
			observe.Copy(Population()[0]);
		}
		oldC = mC;
	}
	Population() = old_pop;

	
	Check(observe);
	observe.Evaluate();mEvas += observe.Size();
	observe.IsSort(false);
	observe.RankSort(false);
	az::mea::sel::SCrowd2 sel2;
	sel2.Select(observe, 10);
	CPopulationMO newpop(Population());


	for (int k = 0 ;k < xdim ;k++)
	{
		for (int i = 0; i < observe.Size(); i++)
		{
			for(int j=0;j<5;j++)
			{
				if (j==0)
			       newpop[i*psize/observe.Size()][k] = observe[i][k];
				else
				{
					newpop[i*psize/observe.Size()+j*2][k] = observe[i][k]+rnd::rand(-0.05,0.05);
					if(newpop[i*psize/observe.Size()+j*2][k]>=P().XUpp(k)) newpop[i*psize/observe.Size()+j*2][k] = P().XUpp(k)-rnd::rand(0.0,0.02);
					else if(newpop[i*psize/observe.Size()+j*2][k]<=P().XLow(k)) newpop[i*psize/observe.Size()+j*2][k] = P().XLow(k)+rnd::rand(0.0,0.02);
				}
			}
		}
	}
	Check(newpop);
	newpop.Evaluate();mEvas += newpop.Size();
	newpop.IsSort(false);
	newpop.RankSort(false);
	az::mea::sel::SCrowd2 sel3;
	sel3.Select(newpop, 20);

	for (int i = 0; i < newpop.Size(); i++)
	{
		for (int k = 0 ;k < xdim ;k++)
		{
			Population()[psize-i-1][k] = newpop[i][k];
		}
	}

	for (int i = 0; i < psize/20; i++)
	{
		/*std::cout<<memory[i].Dominate(Population()[psize-i-21])<<"\t";*/
		if (memory[i].Dominate(Population()[psize-i-21]))
		{
			for (int k = 0 ;k < xdim ;k++)
				Population()[psize-i-21][k] = memory[i][k];
		}
	}
	//std::cout<<std::endl;
	memory = copy_memory;

}

//DMS
void DMOO::InitRG()
{
	//Step 0: initialize the population
	if(mStep == 0) 
	{
		InitRIS(true); 
		return;
	}
	CPopulationMO	search_pop(P()),old_pop(P()),temp_memory(P());
	search_pop = Population(),old_pop = Population(),observe = Population();

	unsigned int i, k, j, xdim=P().XSize(), psize = Population().Size();

	oldC.resize(xdim);
	mC.resize(xdim);
	D.resize(xdim);

	oldmin.resize(xdim);
	oldmax.resize(xdim);
	newmin.resize(xdim);
	newmax.resize(xdim);

	int num_nodominance=0;
	for (i=0; i<psize; i++)
	{
		if (Population()[i].Rank()==1)
			num_nodominance++;
	}
	if(mStep <= mTaoT)
	{
		for(k=0; k<xdim; k++) 
		{
			oldC[k]  = 0.0;
			mC[k] = 0.0;
			oldmin[k] = 1.0E+10;
			oldmax[k] = 1.0E-10;
			newmin[k] = Population()[0][k];
			newmax[k] = Population()[0][k];
			for(i=0; i<psize; i++)
			{
				if (Population()[i].Rank()==1)
				{
					mC[k] += Population()[i][k];
					if(Population()[i][k] < newmin[k])
						newmin[k] = Population()[i][k];
					if(Population()[i][k] > newmax[k])
						newmax[k] = Population()[i][k];
				}
			}
			mC[k] /= double(num_nodominance);
		}
	}
	else
	{
		for (k = 0; k <xdim; k++)
		{
			oldC[k] = mC[k];
			mC[k] = 0.0;
			oldmin[k] = newmin[k];
			oldmax[k] = newmax[k];
			newmin[k] = Population()[0][k];
			newmax[k] = Population()[0][k];
			for(i=0; i<psize; i++)
			{
				if (Population()[i].Rank()==1)
				{
					mC[k] += Population()[i][k];
					if(Population()[i][k] < newmin[k])
						newmin[k] = Population()[i][k];
					if(Population()[i][k] > newmax[k])
						newmax[k] = Population()[i][k];
				}
			}
			mC[k] /= double(num_nodominance);
		}
	}
	for (int i = 0; i < xdim; i++)
	{
		D[i] = mC[i] - oldC[i];
	}
	int a = 0,b = 0,c = 0;

	int Switch = 123;
	if(Switch == 123)
	{
		az::mea::sel::XCrowd sel1;
		sel1.Select(search_pop, psize/10);

		for(i = 0;i < observe.Size(); i++)
		{
			for (k = 0; k < xdim; k++)
				observe[i][k] = search_pop[i%10][k] + (int)( i / 10 ) * ( newmax[k] - oldmax[k] )/10 +rnd::rand(-0.01,0.01);
			observe[i].indexAssign(3);
		}
		Check(observe);

		az::mea::sel::SCrowd2 sel2;
		sel2.Select(old_pop,psize/2);

		for(i = 0;i < Population().Size();i++)
		{
			if(i < psize/2)
			{
				for (k = 0; k < xdim; k++)
					Population()[i][k] = old_pop[i][k] + (mC[k] - oldC[k]) ;
				Population()[i].indexAssign(1);
			}
			else
			{
				for (k = 0; k < xdim; k++)
					Population()[i][k] = rnd::rand(2 * newmin[k] - oldmin[k],2 * newmax[k] - oldmax[k]);
				Population()[i].indexAssign(2);
			}
		}

		Check(Population());

		Population().Combine(observe);

		az::mea::sel::SCrowd2 sel;
		sel.Select(Population(), mPopSize);
		for(int i=0;i<mPopSize;i++)
		{
			if(Population()[i].index==1)a++;
			if(Population()[i].index==2)b++;
			if(Population()[i].index==3)c++;
		}
		/*std::cout<<"第一部分个体数："<<a<<std::endl
			<<"第二部分个体数："<<b<<std::endl
			<<"第三部分个体数："<<c<<std::endl;*/
	}
	if(Switch ==1)
	{
		for(i = 0;i < Population().Size();i++)
		{
			if(i < psize)
				for (k = 0; k < xdim; k++)
					Population()[i][k] = Population()[i][k] + (mC[k] - oldC[k]) ;

		}
	}
	//random search of half the pop
	if(Switch ==11)
	{
		az::mea::sel::SCrowd2 sel2;
		sel2.Select(old_pop,psize/2);
		for(i = 0;i < Population().Size();i++)
		{
			if(i < psize/2)
				for (k = 0; k < xdim; k++)
					Population()[i][k] = old_pop[i][k] + (mC[k] - oldC[k]) ;
			else
				for (k = 0; k < xdim; k++)
					Population()[i][k] = rnd::rand(P().XLow(k),P().XUpp(k));
		}
	}
	if(Switch ==2)
	{
		for(i = 0;i < Population().Size();i++)
		{
			if(i < psize)
				for (k = 0; k < xdim; k++)
					Population()[i][k] = rnd::rand(2 * newmin[k] - oldmin[k],2 * newmax[k] - oldmax[k]) ;
		}
	}
	if(Switch ==3)
	{
		for(i = 0;i < Population().Size();i++)
		{
			if(i < psize)
				for (k = 0; k < xdim; k++)
					Population()[i][k] = Population()[i%10][k] + (int)( i / 10 ) * ( newmax[k] - oldmax[k] )/10 +rnd::rand(-0.01,0.01);
		}
	}
	if(Switch == 23)
	{
		for(i = 0;i < Population().Size();i++)
		{
			if(i < psize)
				for (k = 0; k < xdim; k++)
					Population()[i][k] = rnd::rand(2 * newmin[k] - oldmin[k],2 * newmax[k] - oldmax[k]) ;
		}
		for(i = 0;i < observe.Size();i++)
		{
			if(i < psize)
				for (k = 0; k < xdim; k++)
					observe[i][k] = observe[i%10][k] + (int)( i / 10 ) * ( newmax[k] - oldmax[k] )/10 +rnd::rand(-0.01,0.01);
		}
		Population().Combine(observe);
		az::mea::sel::SCrowd2 sel;
		sel.Select(Population(), mPopSize);
	}
	if(Switch ==12)
	{
		for(i = 0;i < Population().Size();i++)
		{
			if(i < psize/2)
				for (k = 0; k < xdim; k++)
					Population()[i][k] = Population()[i][k] + (mC[k] - oldC[k]) ;
		}
		for(i = 0;i < observe.Size();i++)
		{
			if(i < psize/2)
				for (k = 0; k < xdim; k++)
					observe[i][k] = rnd::rand(2 * newmin[k] - oldmin[k],2 * newmax[k] - oldmax[k]);
		}
		Population().Combine(observe);
		az::mea::sel::SCrowd2 sel;
		sel.Select(Population(), mPopSize);
	}
	if(Switch ==13)
	{
		for(i = 0;i < Population().Size();i++)
		{
			if(i < psize/2)
			{
				for (k = 0; k < xdim; k++)
					Population()[i][k] = Population()[i][k] + (mC[k] - oldC[k]);
				Population()[i].indexAssign(1);
			}
		}
		for(i = 0;i < observe.Size(); i++)
		{
			for (k = 0; k < xdim; k++)
				observe[i][k] = search_pop[i%10][k] + (int)( i / 10 ) * ( newmax[k] - oldmax[k] )/10 +rnd::rand(-0.01,0.01);
			observe[i].indexAssign(3);
		}
		Population().Combine(observe);
		az::mea::sel::SCrowd2 sel;
		sel.Select(Population(), mPopSize);
		for(int i=0;i<mPopSize;i++)
		{
			if(Population()[i].index==1)a++;
			//if(Population()[i].index==2)b++;
			if(Population()[i].index==3)c++;
		}
	}
	Check(Population());
}

//void DMOO::InitRG()
//{
//	//Step 0: initialize the population
//	if(mStep == 0) 
//	{
//		InitRIS(true); 
//		return;
//	}
//	CPopulationMO	search_pop(P()),old_pop(P()),temp_memory(P());
//	search_pop = Population(),old_pop = Population(),observe = Population();
//
//	unsigned int i, k, j, xdim=P().XSize(), psize = Population().Size();
//	unsigned int a = 0;
//	double maxDegree = 1.0E-6;
//
//	oldC.resize(xdim);
//	mC.resize(xdim);
//	D.resize(xdim);
//
//	oldmin.resize(xdim);
//	oldmax.resize(xdim);
//	newmin.resize(xdim);
//	newmax.resize(xdim);
//
//	int num_nodominance=0;
//	for (i=0; i<psize; i++)
//	{
//		if (Population()[i].Rank()==1)
//			num_nodominance++;
//	}
//	if(mStep <= mTaoT)
//	{
//		for(k=0; k<xdim; k++) 
//		{
//			oldC[k]  = 0.0;
//			mC[k] = 0.0;
//			oldmin[k] = 1.0E+10;
//			oldmax[k] = 1.0E-10;
//			newmin[k] = Population()[0][k];
//			newmax[k] = Population()[0][k];
//			for(i=0; i<psize; i++)
//			{
//				if (Population()[i].Rank()==1)
//				{
//					mC[k] += Population()[i][k];
//					if(Population()[i][k] < newmin[k])
//						newmin[k] = Population()[i][k];
//					if(Population()[i][k] > newmax[k])
//						newmax[k] = Population()[i][k];
//				}
//			}
//			mC[k] /= double(num_nodominance);
//		}
//	}
//	else
//	{
//		for (k = 0; k <xdim; k++)
//		{
//			oldC[k] = mC[k];
//			mC[k] = 0.0;
//			oldmin[k] = newmin[k];
//			oldmax[k] = newmax[k];
//			newmin[k] = Population()[0][k];
//			newmax[k] = Population()[0][k];
//			for(i=0; i<psize; i++)
//			{
//				if (Population()[i].Rank()==1)
//				{
//					mC[k] += Population()[i][k];
//					if(Population()[i][k] < newmin[k])
//						newmin[k] = Population()[i][k];
//					if(Population()[i][k] > newmax[k])
//						newmax[k] = Population()[i][k];
//				}
//			}
//			mC[k] /= double(num_nodominance);
//		}
//	}
//	for (int i = 0; i < xdim; i++)
//	{
//		D[i] = mC[i] - oldC[i];
//	}
//	EnvironmentChangeDegree(Population());
//	std::cout<<"环境变化程度"<<Population().degreeReturn()<<std::endl;
//
//	for(i = 0; i < Population().Size(); i++)
//		if(Population()[i].indDegree > maxDegree)
//			maxDegree = Population()[i].indDegree;
//
//	/*if(Population().degreeReturn() < 0.03)
//		Population() = observe;*/
//	//else
//	//{
//		for(i = 0;i < Population().Size();i++)
//		{
//			double ratio = 1 - (maxDegree - Population()[i].indDegree)/(maxDegree - 0.03);
//			//std::cout<<"gfhgggggggggggggggg"<<ratio<<std::endl;
//			if(Population()[i].indDegree < Population().degreeReturn())
//				a++;
//			else
//			{
//				for (k = 0; k < xdim; k++)
//					Population()[i][k] = Population()[i][k] + (mC[k] - oldC[k]) ;//+ rnd::gaussian();
//
//				if(Population()[i][k]>P().XUpp(k)) 
//					Population()[i][k] = Population()[i][k] = rnd::rand(P().XLow(k),P().XUpp(k));
//				else
//					if(Population()[i][k]<P().XLow(k)) Population()[i][k] = Population()[i][k] = rnd::rand(P().XLow(k),P().XUpp(k));
//			}
//		}
//		
//	//}
//	std::cout<<"环境变化小的个体数"<<a<<std::endl;
//}

#include "../neuralNetworkTrainer.h"
#include "../TrainingDataReader.h"
#include <iostream>

int DMOO::InitNNPS(std::string &trainingDataPath, std::string const &outDataPath, uint32_t const &numInputs,
	uint32_t const &numHidden, uint32_t const &numOutputs){

	if (mStep == 0) { InitRIS(true); return; }
	BPN::TrainingDataReader dataReader(trainingDataPath, outDataPath, numInputs, numOutputs);
	if (!dataReader.ReadData())
	{
		return 1;
	}

	// Create neural network
	BPN::Network::Settings networkSettings{ numInputs, numHidden, numOutputs };
	BPN::Network nn(networkSettings);

	// Create neural network trainer
	BPN::NetworkTrainer::Settings trainerSettings;
	trainerSettings.m_learningRate = 0.001;
	trainerSettings.m_momentum = 0.9;
	trainerSettings.m_useBatchLearning = false;
	trainerSettings.m_maxEpochs = 500;
	trainerSettings.m_desiredAccuracy = 90;

	BPN::NetworkTrainer trainer(trainerSettings, &nn);
	trainer.Train(dataReader.GetTrainingData());

	return 0;

	
}



void DMOO::InitPPS()
{	
	//Step 0: initialize the population
	if(mStep == 0) {InitRIS(true); return;}

	unsigned int i, j, k, xdim=P().XSize(), psize = Population().Size();
			
	//Step 1: save center point and shape points
	mC.resize(xdim); 

	for(k=0; k<xdim; k++) 
	{
		mC[k]  = 0.0;
		for(i=0; i<psize; i++) mC[k] += Population()[i][k]; 
		mC[k] /= double(psize);
	}
	
	//把当前中心点mC插入到历史中心点hC的最开始位置
	hC.insert(hC.begin(), mC);
	
	mBest1 = mBest0;
	mBest0 = Population();
	for(i=0; i<psize; i++) for(j=0; j<xdim; j++) mBest0[i][j] -= mC[j];

	//Step 2: if the history info. is not enough for prediction, then use random initialization
	if(hC.size()<2*mMaxOrder) {InitRIS(false); return;}

	//Step 3: find the parent point for each point in mBest0
	std::vector<unsigned int> pindex(psize);
	double dismin, dis;
	for(i=0; i<psize; i++)
	{
		pindex[i] = 0; dismin = 1.0E100;
		for(j=0; j<psize; j++)
		{
			dis = 0.0; for(k=0; k<xdim; k++) dis += (mBest0[i][k] - mBest1[j][k])*(mBest0[i][k] - mBest1[j][k]);
			if(dis<dismin) {dismin = dis; pindex[i] = j;}
		}
	}
	double mstd = 0.0;
	for(k=0; k<xdim; k++)
	{
		//for(i=0; i<psize; i++) mstd += (mBest0[i][k] - mBest1[pindex[i]][k] - mu[k])*(mBest0[i][k] - mBest1[pindex[i]][k] - mu[k]);
		for(i=0; i<psize; i++) mstd += (mBest0[i][k] - mBest1[pindex[i]][k])*(mBest0[i][k] - mBest1[pindex[i]][k]);
	}
	//mstd = sqrt(mstd/double(psize*xdim));
	mstd = mstd/double(psize*xdim);

	//Step 4: predict the center point
	Predict();	
	
	CPopulationMO pop(Population());
	//Step 5: sample new trial solutions
	for(k=0; k<xdim; k++)
	{
		double std = sqrt(mStdC[k]*mStdC[k] + mstd);
		for(i=0; i<psize; i++) 
		{
			Population()[i][k] = pC[k] + mBest0[i][k] + rnd::gaussian() * std;
			//Population()[i][k] = pC[k] + mBest0[i][k] + rnd::gaussian() * mAlpha*mStdC[k];// mAlpha*(mStdC[k] + mstd);
			//Population()[i][k] = pC[k] + 2*mBest0[i][k] - mBest1[pindex[i]][k] + rnd::gaussian() * mAlpha*mStdC[k];// mAlpha*(mStdC[k] + mstd);
			if(     Population()[i][k]>P().XUpp(k)) Population()[i][k] = 0.5*(pop[pindex[i]][k] + P().XUpp(k));
			else if(Population()[i][k]<P().XLow(k)) Population()[i][k] = 0.5*(pop[pindex[i]][k] + P().XLow(k));
		}
	}
}

void DMOO::InitFPS()
{
	//Step 0: initialize the population
	if(mStep == 0) {InitRIS(true); return;}

	unsigned int i, j, k, xdim=P().XSize(), cdim = (P().FSize()+1)*P().XSize(), psize = Population().Size();

	//Step 1: save extreme points
	mC.resize(cdim);
	double dismin, dis;
	std::vector<unsigned int> pp(P().FSize()+1); for(i=0; i<pp.size(); i++) pp[i] = 0;
	for(i=1; i<Population().Size(); i++) 
	{
		for(j=0; j<pp.size()-1; j++) 
			if(Population()[i].F(j) < Population()[pp[j]].F(j)) 
				pp[j] = i;
	}
	dismin = 1.0E100;
	for(i=0; i<Population().Size(); i++) 
	{
		dis = 0;
		for(j=0; j<P().FSize(); j++) dis += (Population()[i].F(j) - Population()[pp[j]].F(j))*(Population()[i].F(j) - Population()[pp[j]].F(j));
		if(dis<dismin) {pp[pp.size()-1] = i; dismin = dis;} // CTI (close-to-idea) point
	}
	for(i=0; i<pp.size(); i++) for(j=0; j<xdim; j++) mC[i*xdim+j] = Population()[pp[i]][j];
	hC.insert(hC.begin(), mC);

	//Step 2: if the history info. is not enough for prediction, then use random initialization
	if(hC.size()<2*mMaxOrder) {InitRIS(false); return;}

	//Step 3: predict the extreme points
	Predict();	

	//Step 4: generate new points
	// permutate the population
	Population().Shuffle();
	CPopulationMO pop(Population());
	// predict
	for(i=0; i<pp.size(); i++)
	{
		for(j=0; j<xdim; j++) Population()[3*i+0][j] = pC[i*xdim+j] + 1.6449*mStdC[j]*(rnd::rand()>0.5?1.0:-1.0);
		for(j=0; j<xdim; j++) Population()[3*i+1][j] = pC[i*xdim+j] + 1.6449*mStdC[j]*(rnd::rand()>0.5?1.0:-1.0);
		for(j=0; j<xdim; j++) Population()[3*i+2][j] = pC[i*xdim+j];
		for(j=0; j<xdim; j++) for(k=0; k<3; k++)
		{
			if(     Population()[3*i+k][j]>P().XUpp(j)) Population()[3*i+k][j] = 0.5*(pop[3*i+k][j] + P().XUpp(j));
			else if(Population()[3*i+k][j]<P().XLow(j)) Population()[3*i+k][j] = 0.5*(pop[3*i+k][j] + P().XLow(j));
		}
	}
	// 70% current + 30% random
	for(i=0; i<(unsigned int)(0.3*(Population().Size()-pp.size()*3)); i++)
	{
		for(j=0; j<xdim; j++) Population()[(unsigned int)(pp.size()*3+1+i)][j] = rnd::rand(P().XLow(j),P().XUpp(j));
	}
}

void DMOO::InitRIS(bool all)
{
	unsigned int i,j;
	if(all)
	{
		Population().Resize(mPopSize);
		for(i=0; i<mPopSize; i++) for(j=0; j<P().XSize(); j++) Population()[i][j] = rnd::rand(P().XLow(j),P().XUpp(j));
	}
	else
	{
		for(i=0; i<mPopSize/2; i++) for(j=0; j<P().XSize(); j++) Population()[i][j] = rnd::rand(P().XLow(j),P().XUpp(j));
	}
}

void DMOO::InitDSP()
{
	//Step 0: initialize the population
	if(mStep == 0) 
	{
		InitRIS(true); 
		return;
	}
	CPopulationMO	search_pop(P()),new_pop(P()),temp_pop(P()),test_pop;
	search_pop = Population(), observe = Population();
	
	unsigned int i, k, j, l, count, xdim=P().XSize(), psize = Population().Size();

	//Define external and inner space partition angles.
	double angle1 = acos(1/sqrt(xdim)), angle2 = angle1/(psize/xdim), beta;
	double minAngle, radius = sin(angle2/2);
	//if there exists individuals in Xlayer j of Xblock i
	std::vector<std::vector<bool>> ifExist;
	std::vector<std::vector<CPopulationMO>> partition_pop;

	ifExist.resize (xdim);
	partition_pop.resize (xdim);
	temp_pop.Resize (psize);

	for(j = 0; j < xdim; j++)
	{
		ifExist[j].resize (psize/xdim);
		partition_pop[j].resize (psize/xdim);
	}

	/*for(i = 0 ; i < psize; i++)
		test_pop.Copy(search_pop[i]);*/

	//std::cout<<test_pop.Size()<<std::endl;

	//call the constructor function of CPopulationMO through below way
	/*for(i = 0; i < xdim; i++)
		for(j = 0; j < psize/xdim; j++)
		{
			partition_pop[i][j].IsSort(true);
			partition_pop[i][j].P(P());
			partition_pop[i][j].IsGridSort(true);
		}*/

	for(i = 0; i< xdim; i++)
		for(j = 0; j < psize/xdim; j++)
			partition_pop[i][j] = temp_pop;

	//std::cout<<temp_pop.Size()<<std::endl;
	// There are xdim Xblocks and psize/xdim Xlayers:
	for(i = 0 ; i < psize; i++)
	{
		count = 0; minAngle = PI*2;
		//search_pop[i].computeCooAngle();
		for(j = 0; j < xdim; j++)
		{
			//std::cout<<i<<" "<<search_pop[i].X(j)<<" "<<search_pop[i].indXAngle[j]<<std::endl;
			//std::cout<<search_pop[i].indXAngle[j]<<" ";
			if(search_pop[i].indXAngle[j] < minAngle)
			{
				minAngle = search_pop[i].indXAngle[j];
				count = j;
			}
		}

		//std::cout<<minAngle<<std::endl;//<<count<<std::endl;

		search_pop[i].Xblock = count;
		k = std::ceil(minAngle/angle2);
		search_pop[i].Xlayer = k - 1;
		ifExist[count][k-1] = 1;
		
		//std::cout<<std::endl<<count<<" "<<k-1<<" "<<std::endl;//<<partition_pop[count][k].Size();
		partition_pop[count][k-1].Copy(search_pop[i]);
		//std::cout<<std::endl<<count<<" "<<partition_pop[count][k].Size()<<std::endl;
	}

	/*for(i = 0 ; i < psize; i++)
	{
		std::cout<<search_pop[i].Xblock<<"\t"<<search_pop[i].Xlayer<<std::endl;
	}*/

	//std::cout<<std::endl<<"abcd  "<<std::endl;
	az::mea::sel::SCrowd2 sel;
	for(i = 0; i < xdim; i++)
		for(j = 0; j < psize/xdim; j++)
		{
			if(ifExist[i][j] == 1)
			{
				partition_pop[i][j].Evaluate();
				sel.Select(partition_pop[i][j],1);
				new_pop.Copy(partition_pop[i][j]);
			}
			else
			{

				//Compute xdim-1 individuals  which every Xlayer intersects with the plane, for example with the xOy plane
				beta = (j + 0.5) * angle2;
				for(k = 0; k < xdim - 1; k++)
				{
					for(l = 0; l < xdim; l++)
					{
						temp_pop[k][l] = P().XLow(l);
						//partition_pop[i][j].Copy(temp_pop[k]);
					}
				}
						

				for(k = 0; k < xdim - 1; k++)
				{
					temp_pop[k][i] = P().XLow(i) + cos(beta);
				}

				for(k = 0; k < xdim - 1; k++)
				{
					if(k < i)
						temp_pop[k][k] = P().XLow(k) + sin(beta);
					else
						temp_pop[k][k+1] = P().XLow(k+1) + sin(beta);
				}
				//the xdim-th individual is the center of the 椎体
				for(k = 0; k < xdim; k++)
				{
					if(k == i)
						temp_pop[xdim - 1][k] = P().XLow(k) + cos(beta);
					else
						temp_pop[xdim - 1][k] = P().XLow(k) + sqrt(pow(sin(beta),2)/(xdim-1));
				}
				//Randomly produce xdim individuals according to the individual above within the radius
				for(k = xdim; k < 2 * xdim; k++)
					for(l = 0; l < xdim; l++)
						temp_pop[k][l] = P().XLow(l) + temp_pop[xdim - 1][l] + radius * rnd::gaussian();

				//Copy 2 * xdim individuals from temp_pop to population partition_pop[i][j]
				for(k = 0; k < 2 * xdim; k++)
					partition_pop[i][j].Copy(temp_pop[k]);

				//std::cout<<temp_pop.Size()<<std::endl;

				partition_pop[i][j].Evaluate();

				sel.Select(partition_pop[i][j],1);
				new_pop.Copy(partition_pop[i][j]);
			}
		}
		//std::cout<<new_pop.Size()<<std::endl;

	//	for(i = 0 ; i < psize; i++)
	//	{
	//	count = 0; minAngle = PI*2;
	//	new_pop[i].computeCooAngle();
	//	for(j = 0; j < xdim; j++)
	//	{
	//		//std::cout<<i<<" "<<search_pop[i].X(j)<<" "<<search_pop[i].indXAngle[j]<<std::endl;
	//		//std::cout<<search_pop[i].indXAngle[j]<<" ";
	//		if(new_pop[i].indXAngle[j] < minAngle)
	//		{
	//			minAngle = new_pop[i].indXAngle[j];
	//			count = j;
	//		}
	//	}

	//	//std::cout<<minAngle<<std::endl;//<<count<<std::endl;

	//	new_pop[i].Xblock = count;
	//	k = std::ceil(minAngle/angle2);
	//	new_pop[i].Xlayer = k - 1;
	//	
	//	std::cout<<std::endl<<new_pop[i].Xblock<<" "<<(minAngle/angle2)<<" ";
	//}

		Check(new_pop);

		Population() = new_pop;

		

}

void DMOO::LocalSearchEGS(CPopulationMO &pop)
{
	unsigned int i, j, k, N =4 ,xdim=P().XSize(), psize = pop.Size(),mObj = P().FSize();

	std::vector<double> F_max(mObj),F_min(mObj);

	double step_size = 1.0, sum = 0, fitness = 0;
	std::vector<double> v(xdim , 0.0);
	CIndividualMO g(P());
	CPopulationMO Q(P());
	

	Q.Resize(N);

	//pop.Evaluate();
	for (i = 0; i < psize; i++)
	{
		std::default_random_engine random(time(NULL));
		std::normal_distribution<double> dis1(0, step_size);

		for (j = 0; j < N; j++)
			for (k = 0; k < xdim; k++)
				Q[j][k] = dis1(random);
		Check(Q);
		Q.Evaluate();
		for (j = 0; j < N; j++)
		{
			fitness = EGS_GP( Q[j].F() , pop[i].F());
			for (k = 0; k < xdim; k++)
			{
				v[k] += fitness * ( Q[j].X(k) - pop[i].X(k) );
			}
		}

		//std::cout<<"第一次EGS_RW："<<std::endl;
		for (k = 0; k < xdim; k++)
			sum += v[k] * v[k];
		sum = sqrt(sum);
		for (k = 0; k < xdim; k++)
		{
			v[k] /= sum;
			g[k] = pop[i].X(k) - step_size * sqrt(psize) * v[k];
		}

		g.Check();
		g.Evaluate();
		step_size = EGS_GP( g.F(), pop[i].F()) < 0 ? step_size * 1.8 : step_size / 1.8;

		//std::cout<<"第二次EGS_RW："<<std::endl;
		if(g.Dominate(pop[i]) == 1)
			for (k = 0; k < xdim; k++)
			{
				pop[i].X(k) = g.X(k);
			}
		else
		if(g.Dominate(pop[i]) == -1)
		{
			if(rnd::rand() < 0.5)
			for (k = 0; k < xdim; k++)
			{
				pop[i].X(k) = g.X(k);
			}
		}
	}

}

double DMOO::EGS_RW(std::vector<double> &q, std::vector<double> &p)
{
	unsigned int i, mObj=P().FSize();
	double sum = 0, fitness = 0;
	std::vector<double> weight(mObj);
	//produce weight randomly
	do
	{
		for (i = 0; i < mObj-1; i++)
			weight[i] = az::rnd::rand();
		for (i = 0; i < mObj-1; i++)
			sum += weight[i];
		weight[mObj-1] = 1-sum;
	} while (weight[mObj-1] < 0);

	/*std::cout<<weight[0]<<"\t"<<weight[1];
	std::cout<<std::endl;*/
	//std::cout<<"目标数："<<mObj<<std::endl;
	for (i = 0; i < mObj; i++)
		fitness += weight[i] * q[i] - weight[i] * p[i];
	return fitness;
}

double DMOO::EGS_GP(std::vector<double> &q , std::vector<double> &p)
{
	unsigned int i, mObj=P().FSize();
	double sum = 0;
	std::vector<double> z(mObj, az::rnd::rand());
	for (i = 0; i < mObj; i++)
		sum += pow(q[i] - (p[i] - z[i]),2);
	return sqrt(sum);
}

double DMOO::EGS_RHV(std::vector<double> &q , std::vector<double> &p)
{
	unsigned int i, mObj=P().FSize();
	double temp = 1;
	for (i = 0; i < mObj; i++)
	{
		temp *= (p[i] + 1 - q[i]);
	}
	return temp;
}

// generate new trial solutions by any offspring generator
CPopulationMO& DMOO::Generate(CPopulationMO& popnew, unsigned int size)
{
	if(mOptimizer == std::string("GTM"))
	{
		az::mea::gen::mod::ModelGTM2 gen;
		gen.Set(25,2,popnew.P().FSize()-1,30,0.25);
		gen.Generate(size, popnew, Population());
	}
	else if(mOptimizer == std::string("PCA"))
	{
		az::mea::gen::mod::RM gen;
		gen.Set(popnew.P().FSize()-1, 5, 30, 0.25);
		gen.Generate(size, popnew, Population());
	}
	else if(mOptimizer == std::string("NSDE"))
	{
		popnew.P().ETA_PM() = 20.0;
		popnew.P().ETA_SBX()= 20.0;
		popnew.P().Pc()		= 0.8;
		popnew.P().Pm()		= 0.05;
		az::mea::gen::XNSDE	gen;
		gen.Set(0.5,1.0);
		gen.Generate(size, popnew, Population());
	}
	else if(mOptimizer == std::string("NSGA"))
	{
		popnew.P().ETA_PM() = 20.0;
		popnew.P().ETA_SBX()= 20.0;
		popnew.P().Pc()		= 0.8;
		popnew.P().Pm()		= 0.05;
		az::mea::gen::XSBX gen;
		gen.Generate(size, popnew, Population());
	}
	else if(mOptimizer == std::string("RM2"))
	{
		az::mea::gen::mod::RM2 gen;
		gen.Set(1.0, 1.0, 0.8, 1,30);
		gen.Generate(size, popnew, Population());
	}
	else if (mOptimizer == std::string("DEE"))
	{
		gd gen;
		gen.generater(size, popnew, Population());
		mStep = mStep + (mTaoT - 2);
	}

	return popnew;
}

} //namespace dea
} //namespace mea
} //namespace az
