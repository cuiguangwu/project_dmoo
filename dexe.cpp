#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "AlgD.h"
#include "emo/Parameter.h"
#include "emo/Config.h"
#include <algorithm>
#include "stdlib.h"
#include "../evaluate.h"
//#include "../Example/"


#include <yvals.h>
#if (_MSC_VER >= 1600)
#define __STDC_UTF_16__
#endif
#include"engine.h"
#include"memory.h"
#pragma comment(lib,"libmat.lib") 
#pragma comment(lib,"libmx.lib") 
#pragma comment(lib,"libeng.lib") 

#include "../neuralNetworkTrainer.h"
#include "../TrainingDataReader.h"
#include <iostream>

#if _MSC_VER
#pragma warning(push, 0)
#pragma warning(disable: 4702)
#endif

#include "../cmdParser.h"

#if _MSC_VER
#pragma warning(pop)
#endif


void average();
int runNN(cli::Parser &cmdParser);

int main(int argc, char* argv[])
{
	int ccc = 0;
	cli::Parser cmdParser(argc, argv);
	
	while(ccc<1)
	{
	std::string method, problem, fname, path;
	unsigned int strategy, runs, generation, popsize, dimension, gen, ir, i, j;
	unsigned int torder, taot, nt; 
	double		 t0, alpha;
	az::mea::CParameter mPar;
	char savefilename[1024];
	//Parameter setting
	//DMOPA-DMOPF : F5 F6 F7 F9 F10 F8  (7 8 9 10 11 12 )
	std::string instances[]  = {"FDA1","FDA2","FDA3","FDA4","DMOP1","DMOP2","DMOP3","DMOPA","DMOPB","DMOPC","DMOPD","DMOPE","DMOPF"}; // names of test instances
	std::string instances2[]  = {"SDMOP1","SDMOP2","SDMOP3","SDMOP4","SDMOP5","SDMOP6","SDMOP7","SDMOP8"};
	std::string instances3[]  = {"JY1","JY2","JY3","JY4","JY5","JY6","JY7","JY8","JY9"};
	std::string arraymethod[]  = {"PCA","GTM","NSDE","NSGA","RM2","DEE"}; // names of method
	int arraystrategy[] = {1 , 2 , 3, 4, 5 , 6, 7};
	//1, 2 , 3, 4 , 5 , 6 , and 7 respectively, corresponding to RIS--0, FPS--1 , PPS--2 , PZ--3 ,EGS--4 and DMS--5, DSP--6(Dimension Space Partition)
	problem = instances[3];
	method = arraymethod[4];
	strategy = arraystrategy[5];
	runs = 1; //Number of runs
	popsize = 100;

	alpha =1;
	taot = 30;
	t0 = 0;
	nt =10;
	torder = 3;
	dimension = 20;

	//int count_rank_1 = 0;

	//

	int a[21] = {10,14,0,19,11,16,13,5,2,8,17,12,15,4,20,9,6,1,3,7,18};
	std::vector<int> t_order;
	for(i=0;i<2*nt+1;i++) 
		t_order.push_back(a[i]);

	sprintf(savefilename,"data/%s_%s_%d_%d_%d",method.c_str(),problem.c_str(),strategy,taot,nt);
	mPar.TolF() = 1.0E-5;//TOLERANCEF
	mPar.TolX() = 1.0E-5;//TOLERANCEX
	mPar.TolC() = 1.0E-5;//TOLERANCEC

	mPar.XSize(dimension);
	mPar.Problem(problem);
	mPar.XCoding() = false;

	std::vector< std::vector<double> > PF0, PF1, PS_1,PF_1,PS1,PS2;//PS1和PS1存的是最后一代种群
	std::vector< std::vector<double> > POF;
	unsigned int ic,ipf0,ipf1, ps_1;

	bool justinit = true;

	std::cout<<savefilename<<std::endl;

	std::ofstream f0,f1,pf,f2,f3,f4;

	generation = 3600;

	//对t_order重新洗牌
	srand(123456);
	std::random_shuffle(t_order.begin(),t_order.end());

	/*for(std::vector<int>::size_type i=0 ; i != t_order.size(); i++)
		std::cout<<t_order[i]<<"\t";
	std::cout<<std::endl;*/

	az::mea::dea::DMOO* pEA = new az::mea::dea::DMOO(strategy, method, popsize, generation, taot, nt, torder, t0, alpha, mPar,t_order);
	//unsigned int xdim = (dimension < 3) ? dimension:3;
	unsigned int xdim = dimension;
	static bool flag =false ;
	int *count_rank_1 = new int[1000]();

	for(ir=0; ir<runs; ir++)
	{
		ps_1 =ic = ipf1 = ipf0 = 0;

		PF0.resize(popsize*20*nt);
		PF_1.resize(popsize*20*nt);
		PS_1.resize(popsize * 20 * nt);
		PF1.resize(popsize * 20 * nt);
		PS1.resize(popsize * 20 * nt);

		PS2.resize(popsize * 20 * nt);
		pEA->Reset();
		
		unsigned int cou = 0;
		
		int k = 0;
		while(!pEA->IsTerminate())
		{
			Engine *ep;
			gen = pEA->Step();
			if(justinit)
			{
				//if(cou != 0) cou = cou+2;
				justinit = false;
				pEA->Population().RankSort();
				
				int cr1 = 0; //count_rank_1
				for (i = 0; i < pEA->Population().Size(); i++)
				{
					PF1[ipf1 + i].resize(pEA->P().FSize());
					PS1[ipf1 + i].resize(xdim);

					PS2[i].resize(xdim);

					for (j = 0; j < pEA->P().FSize(); j++)
					{
						PF1[ipf1 + i][j] = pEA->Population()[i].F(j);
					}
						

					for (j = 0; j < xdim; j++) 
					{
						PS1[ipf1 + i][j] = pEA->Population()[i][j];
						PS2[i][j] = pEA->Population()[i][j];
					}
						

					if (pEA->Population()[i].Rank() == 1)
					{
						PS_1[ps_1 + cr1].resize(xdim);
						PF_1[ps_1 + cr1].resize(pEA->P().FSize());
						for (j = 0; j < xdim; j++) 
							PS_1[ps_1 + cr1][j] = pEA->Population()[i][j];
		
						for (j=0;j<pEA->P().FSize();j++)
							PF_1[ps_1 + cr1][j] = pEA->Population()[i].F(j);

						cr1++;
					}
					
					count_rank_1[k] = cr1;
					//for(j=0; j<pea->p().fsize(); j++) pf1[ipf1+i][j]		= pea->population()[i].f(j);
					
				}
	
				//跑PF需要注释的最上端
				/*double f1[100];
				double f2[100];
				double mNt = nt;
				mxArray *T = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ; 
				mxArray *T2 = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ;
				mxArray *M = mxCreateDoubleMatrix( 1 ,  1 , mxREAL  ) ;
				for(int n=0; n<pEA->Population().Size(); n++)
				{ 
				f1[n]=pEA->Population()[n].X(0);
				f2[n]=pEA->Population()[n].X(1);
				}
				if(  (ep=engOpen(NULL)) )
				{
				memcpy( (char*)mxGetPr(T),(char*)f1 , pEA->Population().Size()*sizeof(double));
				memcpy( (char*)mxGetPr(T2),(char*)f2 , pEA->Population().Size()*sizeof(double));
				memcpy( (char*)mxGetPr(M),(char*)&mNt ,1*sizeof(double));
				engPutVariable(ep , "TX", T) ;
				engPutVariable(ep , "TY", T2) ;
				engPutVariable(ep , "NT", M) ;
				if(flag==false)
				{
				flag =true ;
				engEvalString( ep ,"h=plot(TX,TY,'ro')");
				if(problem=="FDA1"||problem=="DMOP2"||problem=="DMOP3")
				{
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*0/NT);"); 
				engEvalString(ep,"hold on");
				engEvalString(ep,"plot(frontx,fronty,'b.')");	
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*1/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");	
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*2/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");	
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*3/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*4/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*5/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");

				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*11/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*12/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*13/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*14/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*15/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");

				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*21/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*22/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*23/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*24/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*25/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");

				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*31/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*32/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*33/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*34/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=sin(0.5*pi*35/NT);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				}
				if(problem=="FDA2"||problem=="DMOP1")
				{
				engEvalString(ep,"frontx=0:0.01:1 ;fronty=0;"); 
				engEvalString(ep,"hold on");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				}
				if(problem=="DMOPA")
				{
				engEvalString(ep,"frontx=2*cos(pi*0/NT)+2:0.05:2*cos(pi*0/NT)+3 ;fronty=2*sin(2*pi*0/NT)+3-(abs(frontx-2*cos(pi*0/NT)-2)).^(1.25+0.75*sin(pi*0/NT)+2/20);"); 
				engEvalString(ep,"hold on");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*1/NT)+2:0.05:2*cos(pi*1/NT)+3 ;fronty=2*sin(2*pi*1/NT)+3-(abs(frontx-2*cos(pi*1/NT)-2)).^(1.25+0.75*sin(pi*1/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*2/NT)+2:0.05:2*cos(pi*2/NT)+3 ;fronty=2*sin(2*pi*2/NT)+3-(abs(frontx-2*cos(pi*2/NT)-2)).^(1.25+0.75*sin(pi*2/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*3/NT)+2:0.05:2*cos(pi*3/NT)+3 ;fronty=2*sin(2*pi*3/NT)+3-(abs(frontx-2*cos(pi*3/NT)-2)).^(1.25+0.75*sin(pi*3/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*4/NT)+2:0.05:2*cos(pi*4/NT)+3 ;fronty=2*sin(2*pi*4/NT)+3-(abs(frontx-2*cos(pi*4/NT)-2)).^(1.25+0.75*sin(pi*4/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*5/NT)+2:0.05:2*cos(pi*5/NT)+3 ;fronty=2*sin(2*pi*5/NT)+3-(abs(frontx-2*cos(pi*5/NT)-2)).^(1.25+0.75*sin(pi*5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*6/NT)+2:0.05:2*cos(pi*6/NT)+3 ;fronty=2*sin(2*pi*6/NT)+3-(abs(frontx-2*cos(pi*6/NT)-2)).^(1.25+0.75*sin(pi*6/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*7/NT)+2:0.05:2*cos(pi*7/NT)+3 ;fronty=2*sin(2*pi*7/NT)+3-(abs(frontx-2*cos(pi*7/NT)-2)).^(1.25+0.75*sin(pi*7/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*8/NT)+2:0.05:2*cos(pi*8/NT)+3 ;fronty=2*sin(2*pi*8/NT)+3-(abs(frontx-2*cos(pi*8/NT)-2)).^(1.25+0.75*sin(pi*8/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*9/NT)+2:0.05:2*cos(pi*9/NT)+3 ;fronty=2*sin(2*pi*9/NT)+3-(abs(frontx-2*cos(pi*9/NT)-2)).^(1.25+0.75*sin(pi*9/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*10/NT)+2:0.05:2*cos(pi*10/NT)+3 ;fronty=2*sin(2*pi*10/NT)+3-(abs(frontx-2*cos(pi*10/NT)-2)).^(1.25+0.75*sin(pi*10/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");

				engEvalString(ep,"frontx=2*cos(pi*11/NT)+2:0.05:2*cos(pi*11/NT)+3 ;fronty=2*sin(2*pi*11/NT)+3-(abs(frontx-2*cos(pi*11/NT)-2)).^(1.25+0.75*sin(pi*11/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*12/NT)+2:0.05:2*cos(pi*12/NT)+3 ;fronty=2*sin(2*pi*12/NT)+3-(abs(frontx-2*cos(pi*12/NT)-2)).^(1.25+0.75*sin(pi*12/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*13/NT)+2:0.05:2*cos(pi*13/NT)+3 ;fronty=2*sin(2*pi*13/NT)+3-(abs(frontx-2*cos(pi*13/NT)-2)).^(1.25+0.75*sin(pi*13/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*14/NT)+2:0.05:2*cos(pi*14/NT)+3 ;fronty=2*sin(2*pi*14/NT)+3-(abs(frontx-2*cos(pi*14/NT)-2)).^(1.25+0.75*sin(pi*14/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*15/NT)+2:0.05:2*cos(pi*15/NT)+3 ;fronty=2*sin(2*pi*15/NT)+3-(abs(frontx-2*cos(pi*15/NT)-2)).^(1.25+0.75*sin(pi*15/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*16/NT)+2:0.05:2*cos(pi*16/NT)+3 ;fronty=2*sin(2*pi*16/NT)+3-(abs(frontx-2*cos(pi*16/NT)-2)).^(1.25+0.75*sin(pi*16/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*17/NT)+2:0.05:2*cos(pi*17/NT)+3 ;fronty=2*sin(2*pi*17/NT)+3-(abs(frontx-2*cos(pi*17/NT)-2)).^(1.25+0.75*sin(pi*17/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*18/NT)+2:0.05:2*cos(pi*18/NT)+3 ;fronty=2*sin(2*pi*18/NT)+3-(abs(frontx-2*cos(pi*18/NT)-2)).^(1.25+0.75*sin(pi*18/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*19/NT)+2:0.05:2*cos(pi*19/NT)+3 ;fronty=2*sin(2*pi*19/NT)+3-(abs(frontx-2*cos(pi*19/NT)-2)).^(1.25+0.75*sin(pi*19/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");

				engEvalString(ep,"frontx=2*cos(pi*1.5/NT)+2:0.05:2*cos(pi*1.5/NT)+3 ;fronty=2*sin(2*pi*1.5/NT)+3-(abs(frontx-2*cos(pi*1.5/NT)-2)).^(1.25+0.75*sin(pi*1.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*2.5/NT)+2:0.05:2*cos(pi*2.5/NT)+3 ;fronty=2*sin(2*pi*2.5/NT)+3-(abs(frontx-2*cos(pi*2.5/NT)-2)).^(1.25+0.75*sin(pi*2.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*3.5/NT)+2:0.05:2*cos(pi*3.5/NT)+3 ;fronty=2*sin(2*pi*3.5/NT)+3-(abs(frontx-2*cos(pi*3.5/NT)-2)).^(1.25+0.75*sin(pi*3.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*4.5/NT)+2:0.05:2*cos(pi*4.5/NT)+3 ;fronty=2*sin(2*pi*4.5/NT)+3-(abs(frontx-2*cos(pi*4.5/NT)-2)).^(1.25+0.75*sin(pi*4.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*5.5/NT)+2:0.05:2*cos(pi*5.5/NT)+3 ;fronty=2*sin(2*pi*5.5/NT)+3-(abs(frontx-2*cos(pi*5.5/NT)-2)).^(1.25+0.75*sin(pi*5.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*6.5/NT)+2:0.05:2*cos(pi*6.5/NT)+3 ;fronty=2*sin(2*pi*6.5/NT)+3-(abs(frontx-2*cos(pi*6.5/NT)-2)).^(1.25+0.75*sin(pi*6.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*7.5/NT)+2:0.05:2*cos(pi*7.5/NT)+3 ;fronty=2*sin(2*pi*7.5/NT)+3-(abs(frontx-2*cos(pi*7.5/NT)-2)).^(1.25+0.75*sin(pi*7.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*8.5/NT)+2:0.05:2*cos(pi*8.5/NT)+3 ;fronty=2*sin(2*pi*8.5/NT)+3-(abs(frontx-2*cos(pi*8.5/NT)-2)).^(1.25+0.75*sin(pi*8.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*9.5/NT)+2:0.05:2*cos(pi*9.5/NT)+3 ;fronty=2*sin(2*pi*9.5/NT)+3-(abs(frontx-2*cos(pi*9.5/NT)-2)).^(1.25+0.75*sin(pi*9.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*10.5/NT)+2:0.05:2*cos(pi*10.5/NT)+3 ;fronty=2*sin(2*pi*10.5/NT)+3-(abs(frontx-2*cos(pi*10.5/NT)-2)).^(1.25+0.75*sin(pi*10.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");

				engEvalString(ep,"frontx=2*cos(pi*11.5/NT)+2:0.05:2*cos(pi*11.5/NT)+3 ;fronty=2*sin(2*pi*11.5/NT)+3-(abs(frontx-2*cos(pi*11.5/NT)-2)).^(1.25+0.75*sin(pi*11.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*12.5/NT)+2:0.05:2*cos(pi*12.5/NT)+3 ;fronty=2*sin(2*pi*12.5/NT)+3-(abs(frontx-2*cos(pi*12.5/NT)-2)).^(1.25+0.75*sin(pi*12.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*13.5/NT)+2:0.05:2*cos(pi*13.5/NT)+3 ;fronty=2*sin(2*pi*13.5/NT)+3-(abs(frontx-2*cos(pi*13.5/NT)-2)).^(1.25+0.75*sin(pi*13.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*14.5/NT)+2:0.05:2*cos(pi*14.5/NT)+3 ;fronty=2*sin(2*pi*14.5/NT)+3-(abs(frontx-2*cos(pi*14.5/NT)-2)).^(1.25+0.75*sin(pi*14.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*15.5/NT)+2:0.05:2*cos(pi*15.5/NT)+3 ;fronty=2*sin(2*pi*15.5/NT)+3-(abs(frontx-2*cos(pi*15.5/NT)-2)).^(1.25+0.75*sin(pi*15.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*16.5/NT)+2:0.05:2*cos(pi*16.5/NT)+3 ;fronty=2*sin(2*pi*16.5/NT)+3-(abs(frontx-2*cos(pi*16.5/NT)-2)).^(1.25+0.75*sin(pi*16.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*17.5/NT)+2:0.05:2*cos(pi*17.5/NT)+3 ;fronty=2*sin(2*pi*17.5/NT)+3-(abs(frontx-2*cos(pi*17.5/NT)-2)).^(1.25+0.75*sin(pi*17.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*18.5/NT)+2:0.05:2*cos(pi*18.5/NT)+3 ;fronty=2*sin(2*pi*18.5/NT)+3-(abs(frontx-2*cos(pi*18.5/NT)-2)).^(1.25+0.75*sin(pi*18.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*19.5/NT)+2:0.05:2*cos(pi*19.5/NT)+3 ;fronty=2*sin(2*pi*19.5/NT)+3-(abs(frontx-2*cos(pi*19.5/NT)-2)).^(1.25+0.75*sin(pi*19.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(pi*0.5/NT)+2:0.05:2*cos(pi*0.5/NT)+3 ;fronty=2*sin(2*pi*0.5/NT)+3-(abs(frontx-2*cos(pi*0.5/NT)-2)).^(1.25+0.75*sin(pi*0.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				}
				if (problem == "DMOPB")
				{
				engEvalString(ep,"frontx=2*cos(1.5*pi*0/NT)*sin(0.5*pi*0/NT)+2:0.05:2*cos(1.5*pi*0/NT)*sin(0.5*pi*0/NT)+3 ;fronty=2*cos(1.5*pi*0/NT)*cos(0.5*pi*0/NT)+3-(abs(frontx-2*cos(1.5*pi*0/NT)*sin(0.5*pi*0/NT)-2)).^(1.25+0.75*sin(pi*0/NT)+2/20);"); 
				engEvalString(ep,"hold on");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*1/NT)*sin(0.5*pi*1/NT)+2:0.05:2*cos(1.5*pi*1/NT)*sin(0.5*pi*1/NT)+3 ;fronty=2*cos(1.5*pi*1/NT)*cos(0.5*pi*1/NT)+3-(abs(frontx-2*cos(1.5*pi*1/NT)*sin(0.5*pi*1/NT)-2)).^(1.25+0.75*sin(pi*1/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*2/NT)*sin(0.5*pi*2/NT)+2:0.05:2*cos(1.5*pi*2/NT)*sin(0.5*pi*2/NT)+3 ;fronty=2*cos(1.5*pi*2/NT)*cos(0.5*pi*2/NT)+3-(abs(frontx-2*cos(1.5*pi*2/NT)*sin(0.5*pi*2/NT)-2)).^(1.25+0.75*sin(pi*2/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*3/NT)*sin(0.5*pi*3/NT)+2:0.05:2*cos(1.5*pi*3/NT)*sin(0.5*pi*3/NT)+3 ;fronty=2*cos(1.5*pi*3/NT)*cos(0.5*pi*3/NT)+3-(abs(frontx-2*cos(1.5*pi*3/NT)*sin(0.5*pi*3/NT)-2)).^(1.25+0.75*sin(pi*3/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*4/NT)*sin(0.5*pi*4/NT)+2:0.05:2*cos(1.5*pi*4/NT)*sin(0.5*pi*4/NT)+3 ;fronty=2*cos(1.5*pi*4/NT)*cos(0.5*pi*4/NT)+3-(abs(frontx-2*cos(1.5*pi*4/NT)*sin(0.5*pi*4/NT)-2)).^(1.25+0.75*sin(pi*4/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*5/NT)*sin(0.5*pi*5/NT)+2:0.05:2*cos(1.5*pi*5/NT)*sin(0.5*pi*5/NT)+3 ;fronty=2*cos(1.5*pi*5/NT)*cos(0.5*pi*5/NT)+3-(abs(frontx-2*cos(1.5*pi*5/NT)*sin(0.5*pi*5/NT)-2)).^(1.25+0.75*sin(pi*5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*6/NT)*sin(0.5*pi*6/NT)+2:0.05:2*cos(1.5*pi*6/NT)*sin(0.5*pi*6/NT)+3 ;fronty=2*cos(1.5*pi*6/NT)*cos(0.5*pi*6/NT)+3-(abs(frontx-2*cos(1.5*pi*6/NT)*sin(0.5*pi*6/NT)-2)).^(1.25+0.75*sin(pi*6/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*7/NT)*sin(0.5*pi*7/NT)+2:0.05:2*cos(1.5*pi*7/NT)*sin(0.5*pi*7/NT)+3 ;fronty=2*cos(1.5*pi*7/NT)*cos(0.5*pi*7/NT)+3-(abs(frontx-2*cos(1.5*pi*7/NT)*sin(0.5*pi*7/NT)-2)).^(1.25+0.75*sin(pi*7/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*8/NT)*sin(0.5*pi*8/NT)+2:0.05:2*cos(1.5*pi*8/NT)*sin(0.5*pi*8/NT)+3 ;fronty=2*cos(1.5*pi*8/NT)*cos(0.5*pi*8/NT)+3-(abs(frontx-2*cos(1.5*pi*8/NT)*sin(0.5*pi*8/NT)-2)).^(1.25+0.75*sin(pi*8/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*9/NT)*sin(0.5*pi*9/NT)+2:0.05:2*cos(1.5*pi*9/NT)*sin(0.5*pi*9/NT)+3 ;fronty=2*cos(1.5*pi*9/NT)*cos(0.5*pi*9/NT)+3-(abs(frontx-2*cos(1.5*pi*9/NT)*sin(0.5*pi*9/NT)-2)).^(1.25+0.75*sin(pi*9/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*10/NT)*sin(0.5*pi*10/NT)+2:0.05:2*cos(1.5*pi*10/NT)*sin(0.5*pi*10/NT)+3 ;fronty=2*cos(1.5*pi*10/NT)*cos(0.5*pi*10/NT)+3-(abs(frontx-2*cos(1.5*pi*10/NT)*sin(0.5*pi*10/NT)-2)).^(1.25+0.75*sin(pi*10/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");

				engEvalString(ep,"frontx=2*cos(1.5*pi*11/NT)*sin(0.5*pi*11/NT)+2:0.05:2*cos(1.5*pi*11/NT)*sin(0.5*pi*11/NT)+3 ;fronty=2*cos(1.5*pi*11/NT)*cos(0.5*pi*11/NT)+3-(abs(frontx-2*cos(1.5*pi*11/NT)*sin(0.5*pi*11/NT)-2)).^(1.25+0.75*sin(pi*11/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*12/NT)*sin(0.5*pi*12/NT)+2:0.05:2*cos(1.5*pi*12/NT)*sin(0.5*pi*12/NT)+3 ;fronty=2*cos(1.5*pi*12/NT)*cos(0.5*pi*12/NT)+3-(abs(frontx-2*cos(1.5*pi*12/NT)*sin(0.5*pi*12/NT)-2)).^(1.25+0.75*sin(pi*12/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*13/NT)*sin(0.5*pi*13/NT)+2:0.05:2*cos(1.5*pi*13/NT)*sin(0.5*pi*13/NT)+3 ;fronty=2*cos(1.5*pi*13/NT)*cos(0.5*pi*13/NT)+3-(abs(frontx-2*cos(1.5*pi*13/NT)*sin(0.5*pi*13/NT)-2)).^(1.25+0.75*sin(pi*13/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*14/NT)*sin(0.5*pi*14/NT)+2:0.05:2*cos(1.5*pi*14/NT)*sin(0.5*pi*14/NT)+3 ;fronty=2*cos(1.5*pi*14/NT)*cos(0.5*pi*14/NT)+3-(abs(frontx-2*cos(1.5*pi*14/NT)*sin(0.5*pi*14/NT)-2)).^(1.25+0.75*sin(pi*14/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*15/NT)*sin(0.5*pi*15/NT)+2:0.05:2*cos(1.5*pi*15/NT)*sin(0.5*pi*15/NT)+3 ;fronty=2*cos(1.5*pi*15/NT)*cos(0.5*pi*15/NT)+3-(abs(frontx-2*cos(1.5*pi*15/NT)*sin(0.5*pi*15/NT)-2)).^(1.25+0.75*sin(pi*15/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*16/NT)*sin(0.5*pi*16/NT)+2:0.05:2*cos(1.5*pi*16/NT)*sin(0.5*pi*16/NT)+3 ;fronty=2*cos(1.5*pi*16/NT)*cos(0.5*pi*16/NT)+3-(abs(frontx-2*cos(1.5*pi*16/NT)*sin(0.5*pi*16/NT)-2)).^(1.25+0.75*sin(pi*16/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*17/NT)*sin(0.5*pi*17/NT)+2:0.05:2*cos(1.5*pi*17/NT)*sin(0.5*pi*17/NT)+3 ;fronty=2*cos(1.5*pi*17/NT)*cos(0.5*pi*17/NT)+3-(abs(frontx-2*cos(1.5*pi*17/NT)*sin(0.5*pi*17/NT)-2)).^(1.25+0.75*sin(pi*17/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*18/NT)*sin(0.5*pi*18/NT)+2:0.05:2*cos(1.5*pi*18/NT)*sin(0.5*pi*18/NT)+3 ;fronty=2*cos(1.5*pi*18/NT)*cos(0.5*pi*18/NT)+3-(abs(frontx-2*cos(1.5*pi*18/NT)*sin(0.5*pi*18/NT)-2)).^(1.25+0.75*sin(pi*18/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*19/NT)*sin(0.5*pi*19/NT)+2:0.05:2*cos(1.5*pi*19/NT)*sin(0.5*pi*19/NT)+3 ;fronty=2*cos(1.5*pi*19/NT)*cos(0.5*pi*19/NT)+3-(abs(frontx-2*cos(1.5*pi*19/NT)*sin(0.5*pi*19/NT)-2)).^(1.25+0.75*sin(pi*19/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");

				engEvalString(ep,"frontx=2*cos(1.5*pi*0.5/NT)*sin(0.5*pi*0.5/NT)+2:0.05:2*cos(1.5*pi*0.5/NT)*sin(0.5*pi*0.5/NT)+3 ;fronty=2*cos(1.5*pi*1.5/NT)*cos(0.5*pi*0.5/NT)+3-(abs(frontx-2*cos(1.5*pi*0.5/NT)*sin(0.5*pi*0.5/NT)-2)).^(1.25+0.75*sin(pi*0.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*1.5/NT)*sin(0.5*pi*1.5/NT)+2:0.05:2*cos(1.5*pi*1.5/NT)*sin(0.5*pi*1.5/NT)+3 ;fronty=2*cos(1.5*pi*1.5/NT)*cos(0.5*pi*1.5/NT)+3-(abs(frontx-2*cos(1.5*pi*1.5/NT)*sin(0.5*pi*1.5/NT)-2)).^(1.25+0.75*sin(pi*1.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*2.5/NT)*sin(0.5*pi*2.5/NT)+2:0.05:2*cos(1.5*pi*2.5/NT)*sin(0.5*pi*2.5/NT)+3 ;fronty=2*cos(1.5*pi*2.5/NT)*cos(0.5*pi*2.5/NT)+3-(abs(frontx-2*cos(1.5*pi*2.5/NT)*sin(0.5*pi*2.5/NT)-2)).^(1.25+0.75*sin(pi*2.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*3.5/NT)*sin(0.5*pi*3.5/NT)+2:0.05:2*cos(1.5*pi*3.5/NT)*sin(0.5*pi*3.5/NT)+3 ;fronty=2*cos(1.5*pi*3.5/NT)*cos(0.5*pi*3.5/NT)+3-(abs(frontx-2*cos(1.5*pi*3.5/NT)*sin(0.5*pi*3.5/NT)-2)).^(1.25+0.75*sin(pi*3.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*4.5/NT)*sin(0.5*pi*4.5/NT)+2:0.05:2*cos(1.5*pi*4.5/NT)*sin(0.5*pi*4.5/NT)+3 ;fronty=2*cos(1.5*pi*4.5/NT)*cos(0.5*pi*4.5/NT)+3-(abs(frontx-2*cos(1.5*pi*4.5/NT)*sin(0.5*pi*4.5/NT)-2)).^(1.25+0.75*sin(pi*4.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*5.5/NT)*sin(0.5*pi*5.5/NT)+2:0.05:2*cos(1.5*pi*5.5/NT)*sin(0.5*pi*5.5/NT)+3 ;fronty=2*cos(1.5*pi*5.5/NT)*cos(0.5*pi*5.5/NT)+3-(abs(frontx-2*cos(1.5*pi*5.5/NT)*sin(0.5*pi*5.5/NT)-2)).^(1.25+0.75*sin(pi*5.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*6.5/NT)*sin(0.5*pi*6.5/NT)+2:0.05:2*cos(1.5*pi*6.5/NT)*sin(0.5*pi*6.5/NT)+3 ;fronty=2*cos(1.5*pi*6.5/NT)*cos(0.5*pi*6.5/NT)+3-(abs(frontx-2*cos(1.5*pi*6.5/NT)*sin(0.5*pi*6.5/NT)-2)).^(1.25+0.75*sin(pi*6.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*7.5/NT)*sin(0.5*pi*7.5/NT)+2:0.05:2*cos(1.5*pi*7.5/NT)*sin(0.5*pi*7.5/NT)+3 ;fronty=2*cos(1.5*pi*7.5/NT)*cos(0.5*pi*7.5/NT)+3-(abs(frontx-2*cos(1.5*pi*7.5/NT)*sin(0.5*pi*7.5/NT)-2)).^(1.25+0.75*sin(pi*7.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*8.5/NT)*sin(0.5*pi*8.5/NT)+2:0.05:2*cos(1.5*pi*8.5/NT)*sin(0.5*pi*8.5/NT)+3 ;fronty=2*cos(1.5*pi*8.5/NT)*cos(0.5*pi*8.5/NT)+3-(abs(frontx-2*cos(1.5*pi*8.5/NT)*sin(0.5*pi*8.5/NT)-2)).^(1.25+0.75*sin(pi*8.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*9.5/NT)*sin(0.5*pi*9.5/NT)+2:0.05:2*cos(1.5*pi*9.5/NT)*sin(0.5*pi*9.5/NT)+3 ;fronty=2*cos(1.5*pi*9.5/NT)*cos(0.5*pi*9.5/NT)+3-(abs(frontx-2*cos(1.5*pi*9.5/NT)*sin(0.5*pi*9.5/NT)-2)).^(1.25+0.75*sin(pi*9.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*10.5/NT)*sin(0.5*pi*10.5/NT)+2:0.05:2*cos(1.5*pi*10.5/NT)*sin(0.5*pi*10.5/NT)+3 ;fronty=2*cos(1.5*pi*10.5/NT)*cos(0.5*pi*10.5/NT)+3-(abs(frontx-2*cos(1.5*pi*10.5/NT)*sin(0.5*pi*10.5/NT)-2)).^(1.25+0.75*sin(pi*10.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");

				engEvalString(ep,"frontx=2*cos(1.5*pi*11.5/NT)*sin(0.5*pi*11.5/NT)+2:0.05:2*cos(1.5*pi*11.5/NT)*sin(0.5*pi*11.5/NT)+3 ;fronty=2*cos(1.5*pi*11.5/NT)*cos(0.5*pi*11.5/NT)+3-(abs(frontx-2*cos(1.5*pi*11.5/NT)*sin(0.5*pi*11.5/NT)-2)).^(1.25+0.75*sin(pi*11.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*12.5/NT)*sin(0.5*pi*12.5/NT)+2:0.05:2*cos(1.5*pi*12.5/NT)*sin(0.5*pi*12.5/NT)+3 ;fronty=2*cos(1.5*pi*12.5/NT)*cos(0.5*pi*12.5/NT)+3-(abs(frontx-2*cos(1.5*pi*12.5/NT)*sin(0.5*pi*12.5/NT)-2)).^(1.25+0.75*sin(pi*12.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*13.5/NT)*sin(0.5*pi*13.5/NT)+2:0.05:2*cos(1.5*pi*13.5/NT)*sin(0.5*pi*13.5/NT)+3 ;fronty=2*cos(1.5*pi*13.5/NT)*cos(0.5*pi*13.5/NT)+3-(abs(frontx-2*cos(1.5*pi*13.5/NT)*sin(0.5*pi*13.5/NT)-2)).^(1.25+0.75*sin(pi*13.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*14.5/NT)*sin(0.5*pi*14.5/NT)+2:0.05:2*cos(1.5*pi*14.5/NT)*sin(0.5*pi*14.5/NT)+3 ;fronty=2*cos(1.5*pi*14.5/NT)*cos(0.5*pi*14.5/NT)+3-(abs(frontx-2*cos(1.5*pi*14.5/NT)*sin(0.5*pi*14.5/NT)-2)).^(1.25+0.75*sin(pi*14.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*15.5/NT)*sin(0.5*pi*15.5/NT)+2:0.05:2*cos(1.5*pi*15.5/NT)*sin(0.5*pi*15.5/NT)+3 ;fronty=2*cos(1.5*pi*15.5/NT)*cos(0.5*pi*15.5/NT)+3-(abs(frontx-2*cos(1.5*pi*15.5/NT)*sin(0.5*pi*15.5/NT)-2)).^(1.25+0.75*sin(pi*15.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*16.5/NT)*sin(0.5*pi*16.5/NT)+2:0.05:2*cos(1.5*pi*16.5/NT)*sin(0.5*pi*16.5/NT)+3 ;fronty=2*cos(1.5*pi*16.5/NT)*cos(0.5*pi*16.5/NT)+3-(abs(frontx-2*cos(1.5*pi*16.5/NT)*sin(0.5*pi*16.5/NT)-2)).^(1.25+0.75*sin(pi*16.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*17.5/NT)*sin(0.5*pi*17.5/NT)+2:0.05:2*cos(1.5*pi*17.5/NT)*sin(0.5*pi*17.5/NT)+3 ;fronty=2*cos(1.5*pi*17.5/NT)*cos(0.5*pi*17.5/NT)+3-(abs(frontx-2*cos(1.5*pi*17.5/NT)*sin(0.5*pi*17.5/NT)-2)).^(1.25+0.75*sin(pi*17.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*18.5/NT)*sin(0.5*pi*18.5/NT)+2:0.05:2*cos(1.5*pi*18.5/NT)*sin(0.5*pi*18.5/NT)+3 ;fronty=2*cos(1.5*pi*18.5/NT)*cos(0.5*pi*18.5/NT)+3-(abs(frontx-2*cos(1.5*pi*18.5/NT)*sin(0.5*pi*18.5/NT)-2)).^(1.25+0.75*sin(pi*18.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=2*cos(1.5*pi*19.5/NT)*sin(0.5*pi*19.5/NT)+2:0.05:2*cos(1.5*pi*19.5/NT)*sin(0.5*pi*19.5/NT)+3 ;fronty=2*cos(1.5*pi*19.5/NT)*cos(0.5*pi*19.5/NT)+3-(abs(frontx-2*cos(1.5*pi*19.5/NT)*sin(0.5*pi*19.5/NT)-2)).^(1.25+0.75*sin(pi*19.5/NT)+2/20);"); 
				engEvalString(ep,"plot(frontx,fronty,'b.')");


				}
				if(problem == "DMOPC")
				{
				engEvalString(ep,"frontx=1.7*(1-sin(pi*0/NT))*sin(pi*0/NT)+3.4:0.05:1.7*(1-sin(pi*0/NT))*sin(pi*0/NT)+4.4 ;fronty=1.4*(1-sin(pi*0/NT))*cos(pi*0/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*0/NT))*sin(pi*0/NT)-3.4)).^(1.25+0.75*sin(pi*0/NT)+2/20);"); 
				engEvalString(ep,"hold on");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*1/NT))*sin(pi*1/NT)+3.4:0.05:1.7*(1-sin(pi*1/NT))*sin(pi*1/NT)+4.4 ;fronty=1.4*(1-sin(pi*1/NT))*cos(pi*1/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*1/NT))*sin(pi*1/NT)-3.4)).^(1.25+0.75*sin(pi*1/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*2/NT))*sin(pi*2/NT)+3.4:0.05:1.7*(1-sin(pi*2/NT))*sin(pi*2/NT)+4.4 ;fronty=1.4*(1-sin(pi*2/NT))*cos(pi*2/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*2/NT))*sin(pi*2/NT)-3.4)).^(1.25+0.75*sin(pi*2/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*3/NT))*sin(pi*3/NT)+3.4:0.05:1.7*(1-sin(pi*3/NT))*sin(pi*3/NT)+4.4 ;fronty=1.4*(1-sin(pi*3/NT))*cos(pi*3/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*3/NT))*sin(pi*3/NT)-3.4)).^(1.25+0.75*sin(pi*3/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*4/NT))*sin(pi*4/NT)+3.4:0.05:1.7*(1-sin(pi*4/NT))*sin(pi*4/NT)+4.4 ;fronty=1.4*(1-sin(pi*4/NT))*cos(pi*4/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*4/NT))*sin(pi*4/NT)-3.4)).^(1.25+0.75*sin(pi*4/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*5/NT))*sin(pi*5/NT)+3.4:0.05:1.7*(1-sin(pi*5/NT))*sin(pi*5/NT)+4.4 ;fronty=1.4*(1-sin(pi*5/NT))*cos(pi*5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*5/NT))*sin(pi*5/NT)-3.4)).^(1.25+0.75*sin(pi*5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*6/NT))*sin(pi*6/NT)+3.4:0.05:1.7*(1-sin(pi*6/NT))*sin(pi*6/NT)+4.4 ;fronty=1.4*(1-sin(pi*6/NT))*cos(pi*6/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*6/NT))*sin(pi*6/NT)-3.4)).^(1.25+0.75*sin(pi*6/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*7/NT))*sin(pi*7/NT)+3.4:0.05:1.7*(1-sin(pi*7/NT))*sin(pi*7/NT)+4.4 ;fronty=1.4*(1-sin(pi*7/NT))*cos(pi*7/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*7/NT))*sin(pi*7/NT)-3.4)).^(1.25+0.75*sin(pi*7/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*8/NT))*sin(pi*8/NT)+3.4:0.05:1.7*(1-sin(pi*8/NT))*sin(pi*8/NT)+4.4 ;fronty=1.4*(1-sin(pi*8/NT))*cos(pi*8/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*8/NT))*sin(pi*8/NT)-3.4)).^(1.25+0.75*sin(pi*8/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*9/NT))*sin(pi*9/NT)+3.4:0.05:1.7*(1-sin(pi*9/NT))*sin(pi*9/NT)+4.4 ;fronty=1.4*(1-sin(pi*9/NT))*cos(pi*9/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*9/NT))*sin(pi*9/NT)-3.4)).^(1.25+0.75*sin(pi*9/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*10/NT))*sin(pi*10/NT)+3.4:0.05:1.7*(1-sin(pi*10/NT))*sin(pi*10/NT)+4.4 ;fronty=1.4*(1-sin(pi*10/NT))*cos(pi*10/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*10/NT))*sin(pi*10/NT)-3.4)).^(1.25+0.75*sin(pi*10/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");

				engEvalString(ep,"frontx=1.7*(1-sin(pi*0.5/NT))*sin(pi*0.5/NT)+3.4:0.05:1.7*(1-sin(pi*0.5/NT))*sin(pi*0.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*0.5/NT))*cos(pi*0.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*0.5/NT))*sin(pi*0.5/NT)-3.4)).^(1.25+0.75*sin(pi*0.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*1.5/NT))*sin(pi*1.5/NT)+3.4:0.05:1.7*(1-sin(pi*1.5/NT))*sin(pi*1.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*1.5/NT))*cos(pi*1.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*1.5/NT))*sin(pi*1.5/NT)-3.4)).^(1.25+0.75*sin(pi*1.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*2.5/NT))*sin(pi*2.5/NT)+3.4:0.05:1.7*(1-sin(pi*2.5/NT))*sin(pi*2.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*2.5/NT))*cos(pi*2.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*2.5/NT))*sin(pi*2.5/NT)-3.4)).^(1.25+0.75*sin(pi*2.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*3.5/NT))*sin(pi*3.5/NT)+3.4:0.05:1.7*(1-sin(pi*3.5/NT))*sin(pi*3.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*3.5/NT))*cos(pi*3.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*3.5/NT))*sin(pi*3.5/NT)-3.4)).^(1.25+0.75*sin(pi*3.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*4.5/NT))*sin(pi*4.5/NT)+3.4:0.05:1.7*(1-sin(pi*4.5/NT))*sin(pi*4.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*4.5/NT))*cos(pi*4.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*4.5/NT))*sin(pi*4.5/NT)-3.4)).^(1.25+0.75*sin(pi*4.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*5.5/NT))*sin(pi*5.5/NT)+3.4:0.05:1.7*(1-sin(pi*5.5/NT))*sin(pi*5.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*5.5/NT))*cos(pi*5.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*5.5/NT))*sin(pi*5.5/NT)-3.4)).^(1.25+0.75*sin(pi*5.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*6.5/NT))*sin(pi*6.5/NT)+3.4:0.05:1.7*(1-sin(pi*6.5/NT))*sin(pi*6.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*6.5/NT))*cos(pi*6.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*6.5/NT))*sin(pi*6.5/NT)-3.4)).^(1.25+0.75*sin(pi*6.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*7.5/NT))*sin(pi*7.5/NT)+3.4:0.05:1.7*(1-sin(pi*7.5/NT))*sin(pi*7.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*7.5/NT))*cos(pi*7.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*7.5/NT))*sin(pi*7.5/NT)-3.4)).^(1.25+0.75*sin(pi*7.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*8.5/NT))*sin(pi*8.5/NT)+3.4:0.05:1.7*(1-sin(pi*8.5/NT))*sin(pi*8.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*8.5/NT))*cos(pi*8.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*8.5/NT))*sin(pi*8.5/NT)-3.4)).^(1.25+0.75*sin(pi*8.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*9.5/NT))*sin(pi*9.5/NT)+3.4:0.05:1.7*(1-sin(pi*9.5/NT))*sin(pi*9.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*9.5/NT))*cos(pi*9.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*9.5/NT))*sin(pi*9.5/NT)-3.4)).^(1.25+0.75*sin(pi*9.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*10.5/NT))*sin(pi*10.5/NT)+3.4:0.05:1.7*(1-sin(pi*10.5/NT))*sin(pi*10.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*10.5/NT))*cos(pi*10.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*10.5/NT))*sin(pi*10.5/NT)-3.4)).^(1.25+0.75*sin(pi*10.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");

				engEvalString(ep,"frontx=1.7*(1-sin(pi*11/NT))*sin(pi*11/NT)+3.4:0.05:1.7*(1-sin(pi*11/NT))*sin(pi*11/NT)+4.4 ;fronty=1.4*(1-sin(pi*11/NT))*cos(pi*11/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*11/NT))*sin(pi*11/NT)-3.4)).^(1.25+0.75*sin(pi*11/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*12/NT))*sin(pi*12/NT)+3.4:0.05:1.7*(1-sin(pi*12/NT))*sin(pi*12/NT)+4.4 ;fronty=1.4*(1-sin(pi*12/NT))*cos(pi*12/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*12/NT))*sin(pi*12/NT)-3.4)).^(1.25+0.75*sin(pi*12/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*13/NT))*sin(pi*13/NT)+3.4:0.05:1.7*(1-sin(pi*13/NT))*sin(pi*13/NT)+4.4 ;fronty=1.4*(1-sin(pi*13/NT))*cos(pi*13/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*13/NT))*sin(pi*13/NT)-3.4)).^(1.25+0.75*sin(pi*13/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*14/NT))*sin(pi*14/NT)+3.4:0.05:1.7*(1-sin(pi*14/NT))*sin(pi*14/NT)+4.4 ;fronty=1.4*(1-sin(pi*14/NT))*cos(pi*14/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*14/NT))*sin(pi*14/NT)-3.4)).^(1.25+0.75*sin(pi*14/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*15/NT))*sin(pi*15/NT)+3.4:0.05:1.7*(1-sin(pi*15/NT))*sin(pi*15/NT)+4.4 ;fronty=1.4*(1-sin(pi*15/NT))*cos(pi*15/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*15/NT))*sin(pi*15/NT)-3.4)).^(1.25+0.75*sin(pi*15/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*16/NT))*sin(pi*16/NT)+3.4:0.05:1.7*(1-sin(pi*16/NT))*sin(pi*16/NT)+4.4 ;fronty=1.4*(1-sin(pi*16/NT))*cos(pi*16/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*16/NT))*sin(pi*16/NT)-3.4)).^(1.25+0.75*sin(pi*16/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*17/NT))*sin(pi*17/NT)+3.4:0.05:1.7*(1-sin(pi*17/NT))*sin(pi*17/NT)+4.4 ;fronty=1.4*(1-sin(pi*17/NT))*cos(pi*17/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*17/NT))*sin(pi*17/NT)-3.4)).^(1.25+0.75*sin(pi*17/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*18/NT))*sin(pi*18/NT)+3.4:0.05:1.7*(1-sin(pi*18/NT))*sin(pi*18/NT)+4.4 ;fronty=1.4*(1-sin(pi*18/NT))*cos(pi*18/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*18/NT))*sin(pi*18/NT)-3.4)).^(1.25+0.75*sin(pi*18/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*19/NT))*sin(pi*19/NT)+3.4:0.05:1.7*(1-sin(pi*19/NT))*sin(pi*19/NT)+4.4 ;fronty=1.4*(1-sin(pi*19/NT))*cos(pi*19/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*19/NT))*sin(pi*19/NT)-3.4)).^(1.25+0.75*sin(pi*19/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");

				engEvalString(ep,"frontx=1.7*(1-sin(pi*11.5/NT))*sin(pi*11.5/NT)+3.4:0.05:1.7*(1-sin(pi*11.5/NT))*sin(pi*11.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*11.5/NT))*cos(pi*11.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*11.5/NT))*sin(pi*11.5/NT)-3.4)).^(1.25+0.75*sin(pi*11.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*12.5/NT))*sin(pi*12.5/NT)+3.4:0.05:1.7*(1-sin(pi*12.5/NT))*sin(pi*12.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*12.5/NT))*cos(pi*12.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*12.5/NT))*sin(pi*12.5/NT)-3.4)).^(1.25+0.75*sin(pi*12.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*13.5/NT))*sin(pi*13.5/NT)+3.4:0.05:1.7*(1-sin(pi*13.5/NT))*sin(pi*13.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*13.5/NT))*cos(pi*13.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*13.5/NT))*sin(pi*13.5/NT)-3.4)).^(1.25+0.75*sin(pi*13.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*14.5/NT))*sin(pi*14.5/NT)+3.4:0.05:1.7*(1-sin(pi*14.5/NT))*sin(pi*14.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*14.5/NT))*cos(pi*14.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*14.5/NT))*sin(pi*14.5/NT)-3.4)).^(1.25+0.75*sin(pi*14.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*15.5/NT))*sin(pi*15.5/NT)+3.4:0.05:1.7*(1-sin(pi*15.5/NT))*sin(pi*15.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*15.5/NT))*cos(pi*15.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*15.5/NT))*sin(pi*15.5/NT)-3.4)).^(1.25+0.75*sin(pi*15.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*16.5/NT))*sin(pi*16.5/NT)+3.4:0.05:1.7*(1-sin(pi*16.5/NT))*sin(pi*16.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*16.5/NT))*cos(pi*16.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*16.5/NT))*sin(pi*16.5/NT)-3.4)).^(1.25+0.75*sin(pi*16.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*17.5/NT))*sin(pi*17.5/NT)+3.4:0.05:1.7*(1-sin(pi*17.5/NT))*sin(pi*17.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*17.5/NT))*cos(pi*17.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*17.5/NT))*sin(pi*17.5/NT)-3.4)).^(1.25+0.75*sin(pi*17.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*18.5/NT))*sin(pi*18.5/NT)+3.4:0.05:1.7*(1-sin(pi*18.5/NT))*sin(pi*18.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*18.5/NT))*cos(pi*18.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*18.5/NT))*sin(pi*18.5/NT)-3.4)).^(1.25+0.75*sin(pi*18.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				engEvalString(ep,"frontx=1.7*(1-sin(pi*19.5/NT))*sin(pi*19.5/NT)+3.4:0.05:1.7*(1-sin(pi*19.5/NT))*sin(pi*19.5/NT)+4.4 ;fronty=1.4*(1-sin(pi*19.5/NT))*cos(pi*19.5/NT)+3.1-(abs(frontx-1.7*(1-sin(pi*19.5/NT))*sin(pi*19.5/NT)-3.4)).^(1.25+0.75*sin(pi*19.5/NT)+2/20);");
				engEvalString(ep,"plot(frontx,fronty,'b.')");
				}
				}
				else
				{
				engEvalString( ep ,"set(h, 'XData',TX     ,'YData',TY    )");
				}
				}*///跑PF注释的最下端
				//engEvalString(ep,"hold off");
				//engClose(ep);
				//ipf1 += pEA->Population().Size(); //2019.3.2.15:21
				ps_1 += count_rank_1[k];
				k++;
				ipf1 += pEA->Population().Size();

			}

			if(pEA->IsToChange())
			{
				char pfname[1024];
				POF.resize(pEA->Population().Size());
				for(i=0; i<pEA->Population().Size(); i++) 
				{
					PF0[ipf0+i].resize(pEA->P().FSize()+xdim);
					POF[i].resize(pEA->P().FSize());
					for(j=0; j<pEA->P().FSize(); j++)
					{
						PF0[ipf0+i][j]		= pEA->Population()[i].F(j);
						POF[i][j]      = pEA->Population()[i].F(j);
					}
					for(j=0; j<xdim; j++) PF0[ipf0+i][j+pEA->P().FSize()]	= pEA->Population()[i][j];
				}
				//跑PS注释的最上端
				if (pEA->Population().P().FSize()==2)
				{
					double f1[100];
					double f2[100];
					double mNt = nt;
					mxArray *T = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ; 
					mxArray *T2 = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ;
					mxArray *M = mxCreateDoubleMatrix( 1 ,  1 , mxREAL  ) ;
					for(int n=0; n<pEA->Population().Size(); n++)
					{ 
						f1[n]=pEA->Population()[n].F(0);
						f2[n]=pEA->Population()[n].F(1);
					}
					if(  (ep=engOpen(NULL)) )
					{
						memcpy( (char*)mxGetPr(T),(char*)f1 , pEA->Population().Size()*sizeof(double));
						memcpy( (char*)mxGetPr(T2),(char*)f2 , pEA->Population().Size()*sizeof(double));
						memcpy( (char*)mxGetPr(M),(char*)&mNt ,1*sizeof(double));
						engPutVariable(ep , "TX", T) ;
						engPutVariable(ep , "TY", T2) ;
						engPutVariable(ep , "NT", M) ;
						if(flag==false)
						{
							flag =true ;
							engEvalString( ep ,"h=plot(TX,TY,'ro');grid on");
							if(problem=="FDA1"||problem=="DMOP3"||problem=="SDMOP1")
							{
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^0.5 ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");						   
							}
							if(problem=="FDA2")
							{
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1.75*1.75)*10+0.75)) ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");	
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*1/NT)))*(1+(0.75+0.7*sin(0.5*pi*1/NT)))*10+(0.75+0.7*sin(0.5*pi*1/NT)) )) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");	
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*2/NT)))*(1+(0.75+0.7*sin(0.5*pi*2/NT)))*10+(0.75+0.7*sin(0.5*pi*2/NT)) )) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");	
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*3/NT)))*(1+(0.75+0.7*sin(0.5*pi*3/NT)))*10+(0.75+0.7*sin(0.5*pi*3/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");	
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*4/NT)))*(1+(0.75+0.7*sin(0.5*pi*4/NT)))*10+(0.75+0.7*sin(0.5*pi*4/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");	
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*5/NT)))*(1+(0.75+0.7*sin(0.5*pi*5/NT)))*10+(0.75+0.7*sin(0.5*pi*5/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");	
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*11/NT)))*(1+(0.75+0.7*sin(0.5*pi*11/NT)))*10+(0.75+0.7*sin(0.5*pi*11/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*12/NT)))*(1+(0.75+0.7*sin(0.5*pi*12/NT)))*10+(0.75+0.7*sin(0.5*pi*12/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*13/NT)))*(1+(0.75+0.7*sin(0.5*pi*13/NT)))*10+(0.75+0.7*sin(0.5*pi*13/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*14/NT)))*(1+(0.75+0.7*sin(0.5*pi*14/NT)))*10+(0.75+0.7*sin(0.5*pi*14/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*15/NT)))*(1+(0.75+0.7*sin(0.5*pi*15/NT)))*10+(0.75+0.7*sin(0.5*pi*15/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*21/NT)))*(1+(0.75+0.7*sin(0.5*pi*21/NT)))*10+(0.75+0.7*sin(0.5*pi*21/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*22/NT)))*(1+(0.75+0.7*sin(0.5*pi*22/NT)))*10+(0.75+0.7*sin(0.5*pi*22/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*23/NT)))*(1+(0.75+0.7*sin(0.5*pi*23/NT)))*10+(0.75+0.7*sin(0.5*pi*23/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*24/NT)))*(1+(0.75+0.7*sin(0.5*pi*24/NT)))*10+(0.75+0.7*sin(0.5*pi*24/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*25/NT)))*(1+(0.75+0.7*sin(0.5*pi*25/NT)))*10+(0.75+0.7*sin(0.5*pi*25/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*31/NT)))*(1+(0.75+0.7*sin(0.5*pi*31/NT)))*10+(0.75+0.7*sin(0.5*pi*31/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*32/NT)))*(1+(0.75+0.7*sin(0.5*pi*32/NT)))*10+(0.75+0.7*sin(0.5*pi*32/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*33/NT)))*(1+(0.75+0.7*sin(0.5*pi*33/NT)))*10+(0.75+0.7*sin(0.5*pi*33/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*34/NT)))*(1+(0.75+0.7*sin(0.5*pi*34/NT)))*10+(0.75+0.7*sin(0.5*pi*34/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*35/NT)))*(1+(0.75+0.7*sin(0.5*pi*35/NT)))*10+(0.75+0.7*sin(0.5*pi*35/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
							}
							if(problem=="DMOP1"||problem=="DMOP2")
							{
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^1.25 ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");		
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*1/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*2/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*3/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*4/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*5/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*6/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*7/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*8/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*9/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*10/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");

								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*21/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*22/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*23/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*24/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*25/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*26/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*27/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*28/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*29/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*30/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
							}
							if(problem=="DMOPA"||problem=="DMOPB"||problem=="DMOPC"||problem=="DMOPD"||problem=="DMOPE")
							{
								engEvalString(ep,"m=0:0.01:1,frontx=m.^1.25 ;fronty=(1-m).^1.25 ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");		
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");		
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*2/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*2/NT)) ;");  
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*3/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*3/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*4/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*4/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*5/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*5/NT)) ;");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*11/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*11/NT)) ;");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*12/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*12/NT)) ;");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*13/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*13/NT)) ;");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*14/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*14/NT)) ;");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*15/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*15/NT)) ;");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
							}

							//SDMOP2与SDMOP3的动态变化的PF
							if(problem=="SDMOP2"||problem=="SDMOP3")
							{
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^2 ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");		
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*1/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*2/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*3/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*4/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*5/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*6/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*7/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*8/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*9/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*10/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");

								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*11/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*12/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*13/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*14/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*15/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*16/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*17/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*18/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*19/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*20/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
							}
							if(problem == "SDMOP4")
							{
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+1)*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*1/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*2/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*3/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*4/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*5/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*6/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*7/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*8/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*9/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*10/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
							}

							if(problem == "FDA3")
							{
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*0/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*0/NT)))))*(1+abs(sin(0.5*pi*0/NT))) ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*1/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*1/NT)))))*(1+abs(sin(0.5*pi*1/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*2/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*2/NT)))))*(1+abs(sin(0.5*pi*2/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*3/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*3/NT)))))*(1+abs(sin(0.5*pi*3/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*4/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*4/NT)))))*(1+abs(sin(0.5*pi*4/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*5/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*5/NT)))))*(1+abs(sin(0.5*pi*5/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*11/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*11/NT)))))*(1+abs(sin(0.5*pi*11/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*12/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*12/NT)))))*(1+abs(sin(0.5*pi*12/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*13/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*13/NT)))))*(1+abs(sin(0.5*pi*13/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*14/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*14/NT)))))*(1+abs(sin(0.5*pi*14/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*15/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*15/NT)))))*(1+abs(sin(0.5*pi*15/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");

								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*21/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*21/NT)))))*(1+abs(sin(0.5*pi*21/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*22/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*22/NT)))))*(1+abs(sin(0.5*pi*22/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*23/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*23/NT)))))*(1+abs(sin(0.5*pi*23/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*24/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*24/NT)))))*(1+abs(sin(0.5*pi*24/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*25/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*25/NT)))))*(1+abs(sin(0.5*pi*25/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*31/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*31/NT)))))*(1+abs(sin(0.5*pi*31/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*32/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*32/NT)))))*(1+abs(sin(0.5*pi*32/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*33/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*33/NT)))))*(1+abs(sin(0.5*pi*33/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*34/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*34/NT)))))*(1+abs(sin(0.5*pi*34/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*35/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*35/NT)))))*(1+abs(sin(0.5*pi*35/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
							}
						}
						else
						{
							engEvalString( ep ,"set(h, 'XData',TX     ,'YData',TY   )");
						}
					}
				}

				if (pEA->Population().P().FSize()==3)
				{
					double f1[200];
					double f2[200];
					double f3[200];
					double mNt = nt;
					mxArray *T = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ; 
					mxArray *T2 = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ;
					mxArray *T3 = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ;
					mxArray *M = mxCreateDoubleMatrix( 1 ,  1 , mxREAL  ) ;
					for(int n=0; n<pEA->Population().Size(); n++)
					{ 
						f1[n]=pEA->Population()[n].F(0);
						f2[n]=pEA->Population()[n].F(1);
						f3[n]=pEA->Population()[n].F(2);
					}
					if(  (ep=engOpen(NULL)) )
					{
						memcpy( (char*)mxGetPr(T),(char*)f1 , pEA->Population().Size()*sizeof(double));
						memcpy( (char*)mxGetPr(T2),(char*)f2 , pEA->Population().Size()*sizeof(double));
						memcpy( (char*)mxGetPr(T3),(char*)f3 , pEA->Population().Size()*sizeof(double));
						memcpy( (char*)mxGetPr(M),(char*)&mNt ,1*sizeof(double));
						engPutVariable(ep , "TX", T) ;
						engPutVariable(ep , "TY", T2) ;
						engPutVariable(ep , "TZ", T3) ;
						engPutVariable(ep , "NT", M) ;
						if(flag==false)
						{
							flag =true ;
							engEvalString( ep ,"h=plot3(TX,TY,TZ,'r.')");
							if (problem == "FDA4"|| problem == "DMOPF")
							{
								engEvalString(ep,"t=linspace(0,pi/2,25);p=linspace(0,pi/2,25);[theta,phi]=meshgrid(t,p);x=cos(theta).*cos(phi);y=cos(theta).*sin(phi);z=sin(theta);"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot3(x,y,z,'bo');grid on");		
							}

							if(problem == "SDMOP5" || problem == "SDMOP6" || problem == "SDMOP7")
							{
								engEvalString(ep,"[t1 t2]=meshgrid(linspace(0,1,20),linspace(0,1,20)) ;"); 
								engEvalString(ep,"hold on");
							    engEvalString(ep,"x=cos(0.5*pi*t1).*cos(0.5*pi*t2);  y=cos(0.5*pi*t1).*sin(0.5*pi*t2);  z=sin(0.5*pi*t1);"); 
								engEvalString(ep,"alpha(surf(x,y,z),0.8);   re=[0,1,0];  colormap(re); view(69,18);");
							}
							if(problem == "SDMOP8")
							{
								engEvalString(ep,"[t1 t2]=meshgrid(linspace(0,1,20),linspace(0,1,20)) ;"); 
								engEvalString(ep,"hold on");
							    engEvalString(ep,"x=t1;  y=t2;  z=3-0.5*x.*x.*(1+sin(3*pi*x))-0.5*y.*y.*(1+sin(3*pi*y));));"); 
								engEvalString(ep,"alpha(surf(x,y,z),0.8);   re=[0,1,0];  colormap(re); view(69,18);"); 
							}

						}
						else
						{
							engEvalString( ep ,"set(h, 'XData',TX     ,'YData',TY  ,'ZData',TZ  )");

						}
					}
				}//跑PS点图注释的最下端
				engEvalString(ep,"hold off");
				//engClose(ep);
				//system("pause");
				ipf0 += pEA->Population().Size();
				justinit = true;
				std::cout<<"T"<<gen<<"\n";
				sprintf(pfname,"Example/pf_%s_%d_%d",problem.c_str(),ccc,gen/taot);   //C://Users//ruangan//Desktop//阮干//MOEA//实验数据//rg
				std::stringstream spf;
				spf<<pfname<<".csv";
				pf.open(spf.str().c_str());
				pf<<std::scientific<<std::setprecision(5);

				for(i=0; i<pEA->Population().Size(); i++)
				{
					for(j=0; j<pEA->P().FSize(); j++) pf<<POF[i][j]<<",";
					pf<<std::endl;
				}
				
				std::stringstream ss5;
				char psfilename[100];
				sprintf(psfilename, "Example/Ps_%s_%s_%d_%d_%d_%d.csv", method.c_str(), problem.c_str(), strategy, taot, nt, gen / taot);
				//ss5 << psfilename << "last_ps_" << ir << ".csv";
				ss5 << psfilename;
				f4.open(ss5.str().c_str());
				//std::cout<<ss0.str()<<std::endl;
				f4 << std::scientific << std::setprecision(5);
				for (i = 0; i < pEA->Population().Size(); i++)
				{
					for (j = 0; j < xdim; j++)
						f4 << PS2[i][j] << ",";
					f4 << std::endl;
				}
				f4.close();


				POF.clear();
				pf.close();
			}
			cou++;
		}
		//std::cout<<ir<<" ";
		// save data
		std::stringstream ss0,ss1,ps1,ss2,ss3,ss4;
		ps1<<savefilename<<"ps_"<<ir<<".csv";
		//std::cout<<ss1.str()<<std::endl;
		f2.open(ps1.str().c_str());
		f2<<std::scientific<<std::setprecision(5);
	
		printf("%d", k);
		int loc;
		for (int e = 0; e < k; e++) {
			if (!e)
				loc = 0;
			else
				loc += count_rank_1[e - 1];
			for (i = 0; i < count_rank_1[e]; i++)
			{
				for (j = 0; j < xdim; j++) 
					f2 << PS_1[loc+i][j] << ",";
				f2 << std::endl;
			}
			f2 << std::endl;
		}
	
		f2.close();

		ss2 << savefilename << "pf_" << ir << ".csv";
		f1.open(ss2.str().c_str());
		f1 << std::scientific << std::setprecision(5);
		//std::cout<<ss1.str()<<std::endl;
	
		for (int e = 0; e < k; e++) {
			if (!e)
				loc = 0;
			else
				loc += count_rank_1[e - 1];
			for (i = 0; i < count_rank_1[e]; i++)
			{
				for (j = 0; j < pEA->P().FSize(); j++)
					f1 << PF_1[loc + i][j] << ",";
				f1 << std::endl;
			}
			f1 << std::endl;
		}

		f1.close();
		//std::stringstream ss5;
		//char psfilename[100];
		//sprintf(psfilename, "Example/Ps_%s_%s_%d_%d_%d_%d.csv", method.c_str(), problem.c_str(), strategy, taot, nt, gen / taot);
		////ss5 << psfilename << "last_ps_" << ir << ".csv";
		//ss5 << psfilename;
		//f4.open(ss5.str().c_str());
		////std::cout<<ss0.str()<<std::endl;
		//f4 << std::scientific << std::setprecision(5);
		//for (i = 0; i < ipf1; i++)
		//{
		//	for (j = 0; j < xdim; j++)
		//		f4 << PS1[i][j] << ",";
		//	f4 << std::endl;
		//}
		//f4.close();
		
		ss3 << savefilename << "last_pf_" << ir << ".csv";
		f3.open(ss3.str().c_str());
		//std::cout<<ss0.str()<<std::endl;
		f3 << std::scientific << std::setprecision(5);
		for (i = 0; i < ipf1; i++)
		{
			for (j = 0; j < pEA->P().FSize() ; j++) 
				f3 << PF1[i][j] << ",";
			f3 << std::endl;
			if(ipf1 % pEA->Population().Size() ==0)
				f3 << std::endl;
		}
		f3.close();

		
		//ss0<<savefilename<<"_"<<ir<<".pop";
		//f0.open(ss0.str().c_str());
		////std::cout<<ss0.str()<<std::endl;
		//f0<<std::scientific<<std::setprecision(5);
		//for(i=0; i<ipf0; i++)
		//{
		//	for(j=0; j<pEA->P().FSize()+xdim; j++) f0<<PF0[i][j]<<"\t";
		//	f0<<std::endl;
		//}
		//f0.close();

		//PF0.clear();
		PF_1.clear();
		PS_1.clear();
		PF1.clear();
		PS1.clear();
		//average();
		
	}
	//delete pEA;
	//system("pause");
	//std::cout<<std::endl;
	ccc++;
	}
	//runNN(cmdParser);
	return 1;
}


void average() {

	int ccc = 0;
	int nt = 10;//环境变化强度
	int run = 1;    //独立运行次数
	int predict = 120; //预测代数
	int metrics = 5; //指标数
	int testNumber = 17;//测试函数个数

	std::ofstream pf;
	std::ofstream sta;

	//char* instances[]  = {"FDA1","FDA2","FDA3","FDA4","DMOP1","DMOP2","DMOP3","DMOPA","DMOPB","DMOPC","DMOPD","DMOPE","JY1","JY2","JY3","JY4","JY5","JY6","JY7","JY8","JY9","FDA4","DMOPF"}; // names of test instances

	char* instances[] = { "FDA1","FDA2","FDA3","FDA4","FDA5","DMOP1","DMOP2","DMOP3","JY1","JY2","JY3","JY4","JY5","JY6","JY7","JY8","JY9" }; // names of test instances
	//char* instances[]= {"FDA4"};
	while (ccc < run)
	{
		evaluate eva;

		vector<double> gd;
		char filename1[1024];
		char filename2[1024];
		char    strTestInstance[256];
		for (int i = 0; i <1;i++)
		{
			if (instances[i] == "JY1" || instances[i] == "JY2" || instances[i] == "JY3" || instances[i] == "JY4" || instances[i] == "JY5"
				|| instances[i] == "JY6" || instances[i] == "JY7" || instances[i] == "JY8" || instances[i] == "JY9") {
				char pfname[1024];
				sprintf(pfname, "evaluate/data/%s_%d.dat", instances[i], ccc);

				pf.open(pfname, ios::trunc);

				double ave = 0.0, fangcha = 0.0;
				int gen = 0;
				for (gen = 1; gen <= predict; gen++)
				{
					sprintf(strTestInstance, "%s", instances[i]);
					//sprintf(filename1,"PF/pf_%s.dat",strTestInstance);
					//eva.loadpfront(filename1,eva.pf,2);	//真实的PF放到文件中	
					eva.getPOF(instances[i], nt, gen, eva.pf);  //取真实的点

					sprintf(filename2, "PF/pf_%s_%d_%d.dat", strTestInstance, ccc, gen);
					eva.loadpfront(filename2, eva.p, 2);

					//pf<<setprecision(7)<<setiosflags(ios::scientific);
					pf << eva.indicator_GD() << "               ";
					pf << eva.indicator_IGD() << "               ";
					pf << eva.indicator_SP(eva.p) << "               ";
					pf << eva.indicator_MS(eva.pf, eva.p, 2) << "               ";
					pf << eva.indicator_hvd(eva.p, gen, instances[i], nt) << "\n";
					cout << "ccc = " << ccc << "    " << instances[i] << "  gen = " << gen << "  " << endl;
					eva.p.clear();
					eva.pf.clear();
				}
				pf.close();
			}

			if (instances[i] == "FDA1" || instances[i] == "FDA2" || instances[i] == "FDA3" || instances[i] == "DMOP1" || instances[i] == "DMOP2" || instances[i] == "DMOP3" ||
				instances[i] == "DMOPA" || instances[i] == "DMOPB" || instances[i] == "DMOPC" || instances[i] == "DMOPD" || instances[i] == "DMOPE")
			{
				char pfname[1024];
				int gen = 0;
				double ave = 0.0;
				sprintf(pfname, "evaluate/data/%s_%d.dat", instances[i], ccc);

				//pf.open(pfname);
				pf.open(pfname, ios::trunc);
				for (gen = 1; gen <= predict; gen++)
				{
					sprintf(strTestInstance, "%s", instances[i]);
					eva.getPOF(instances[i], nt, gen, eva.pf);  //取真实的点

					sprintf(filename2, "PF/pf_%s_%d_%d.dat", strTestInstance, ccc, gen);
					eva.loadpfront(filename2, eva.p, 2);

					pf << eva.indicator_GD() << "               ";
					pf << eva.indicator_IGD() << "               ";
					pf << eva.indicator_SP(eva.p) << "               ";
					pf << eva.indicator_MS(eva.pf, eva.p, 2) << "               ";
					pf << eva.indicator_hvd(eva.p, gen, instances[i], nt) << "\n";
					cout << "ccc = " << ccc << "    " << instances[i] << "  gen = " << gen << "  " << endl;;

					eva.p.clear();
					eva.pf.clear();
				}
				pf.close();
			}

			if (instances[i] == "DMOPF" || instances[i] == "FDA4" || instances[i] == "FDA5")
			{
				char pfname[1024];
				double ave = 0.0;
				int gen = 0;
				sprintf(pfname, "evaluate/data/%s_%d.dat", instances[i], ccc);

				pf.open(pfname, ios::trunc);

				for (gen = 1; gen <= predict; gen++)
				{
					sprintf(strTestInstance, "%s", instances[i]);


					eva.getPOF(instances[i], nt, gen, eva.pf);  //取真实的点

					sprintf(filename2, "PF/pf_%s_%d_%d.dat", strTestInstance, ccc, gen);
					eva.loadpfront(filename2, eva.p, 3);

					pf << eva.indicator_GD() << "               ";
					pf << eva.indicator_IGD() << "               ";
					pf << eva.indicator_SP(eva.p) << "               ";
					pf << eva.indicator_MS(eva.pf, eva.p, 2) << "               ";
					pf << eva.indicator_hvd(eva.p, gen, instances[i], nt) << "\n";
					cout << "ccc = " << ccc << "    " << instances[i] << "  gen = " << gen << "  " << endl;
					eva.p.clear();
					eva.pf.clear();
				}
				pf.close();
			}

		}
		ccc++;
	}


	//求样本平均值和样本方差
	int kkk = 0, iii;
	evaluate eva;
	for (int kkk = 0; kkk < 1; kkk++)
	{
		char pfname[1024];

		double ave = 0.0;
		int t = 0;
		sprintf(pfname, "evaluate/avg/%s.dat", instances[kkk]);

		pf.open(pfname, ios::trunc);

		char    strTestInstance[256];
		char    strTestInstanceAv[256];

		cout << instances[kkk] << endl;

		pf << setprecision(5) << setiosflags(ios::scientific);
		for (ccc = 0; ccc < run; ccc++) {     //每次独立运行的求平均
			sprintf(strTestInstance, "evaluate/data/%s_%d.dat", instances[kkk], ccc);
			eva.loadpfront(strTestInstance, eva.pf, metrics);
			eva.indicator_AVG(eva.pf, eva.in_avg);
			for (iii = 0; iii < eva.in_avg.size(); iii++) {
				pf << eva.in_avg[iii] << "      ";
			}
			pf << "\n";
		}
		pf.close();

		//独立run次的平均和方差
		pf.open(pfname, ios::app);
		pf << setprecision(5) << setiosflags(ios::scientific);
		sprintf(strTestInstanceAv, "evaluate/avg/%s.dat", instances[kkk]);
		eva.loadpfront(strTestInstanceAv, eva.pf, metrics);
		eva.indicator_AVG(eva.pf, eva.in_avg);
		for (iii = 0; iii < eva.in_avg.size(); iii++) {
			pf << eva.in_avg[iii] << "      ";
		}
		pf << "\n";
		int size1 = eva.pf.size();
		int size2 = eva.pf[0].size();
		vector<double> values;
		for (iii = 0; iii < size2; iii++) {
			for (ccc = 0; ccc < size1; ccc++) {
				values.push_back(eva.pf[ccc][iii]);
			}
			pf << eva.STD(values) << "      ";
			values.clear();
		}
		pf << "\n";
		pf.close();

		//求平均IGD
		eva.GetAvgIGD(instances[kkk], run, predict);
		eva.Statistics(instances[kkk], metrics, run);//统计数据
	}
}



int runNN(cli::Parser &cmdParser)
{

	cmdParser.set_required<std::string>("d", "DataFile", "Path to training data csv file.");
	cmdParser.set_required<std::string>("o", "OutFile", "Path to training data csv file.");
	cmdParser.set_required<uint32_t>("in", "NumInputs", "Num Input neurons.");
	cmdParser.set_required<uint32_t>("hidden", "NumHidden", "Num Hidden neurons.");
	cmdParser.set_required<uint32_t>("out", "NumOutputs", "Num Output neurons.");

	if (!cmdParser.run())
	{
		std::cout << "Invalid command line arguments";
		return 1;
	}

	std::string trainingDataPath = cmdParser.get<std::string>("d").c_str();
	std::string outDataPath = cmdParser.get<std::string>("o").c_str();

	uint32_t const numInputs = cmdParser.get<uint32_t>("in");
	uint32_t const numHidden = cmdParser.get<uint32_t>("hidden");
	uint32_t const numOutputs = cmdParser.get<uint32_t>("out");

	BPN::TrainingDataReader dataReader(trainingDataPath, outDataPath, numInputs, numOutputs);
	if (!dataReader.ReadData() )
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