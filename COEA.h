#ifndef _COEA_H
#define _COEA_H

#include "EA.h"
#include "sub.h" 

class COEA {
public:
	COEA() {
		m_radius=1.0/POP2_SIZE;
		RND.seed(SEED);
		//for(int i=0;i<OBJ_NUM;i++){m_low[i]=-1.0;m_high[i]=1.0;}  //set the lower and upper bounds

		m_gens=0;
		m_acceptdist1=0;
		m_acceptdist2=0;
		m_compcount=0;
		count_enter=0;
		m_change = false;
		for(int i=0;i<GENS;i++){
			m_GDvarArray[i]=0;
			m_GDobjArray[i]=0;
			for(int j=0;j<VAR_NUM;j++){
				m_varTrace[i][j]=0;
				m_PopvarTrace[i][j]=0;
				m_Popvar_max_Trace[i][j]=0;
				m_Popvar_min_Trace[i][j]=0;
				m_Pop_Trace[i][j]=-1;
				m_repArray[j]=j;
			}
		}
		m_rand=false;
	}
	~COEA() {}

	int compare_to_archive(const IND &ind);
	int rank(IND &ind);
	double niche(IND &ind);
	void updateArchive(BinIND &ind);
	void extend(int nClone);void extend2();
	bool comp(IND &ind1,IND &ind2);
	const BinIND& tournament(int tour=2);
	void evaluate(int i);
	void evaluateC(int i,int*rnd);
	void evaluateD(int id, int *idx, int pt);
	void evaluateE(int id);
	void evaluateE4(int id,int *idx);
	//void runCCEA(char *filename, int run_count, int mode);
	//void runCOEA(char *filename, int run_count, int mode);
       // void runCOEA(char *filename, char *filename2, int run_count, int mode);
	//void runSCCEA(char *filename, int run_count, int mode);
	void runSCOEA(char *filename, int run_count, int mode);
	void writeVariables(char*filename,const POP& pop);
	void implicitCompetitionNoRand();
	void implicitCompetition();
	void reinit();
	void memory1();
	void memory2();

	void MS(POP &pop1);
	void S(POP &pop1);

	//void selectionSort(double *pA,int n,int pass);
	void extend3();
	double dynamicDistance();
	void trackingDynamicsI();
	void trackingDynamicsII();

	void runT5(char *filename);

public:
	SubPOP m_subpop[VAR_NUM];
	POP m_pop2;POP m_memory;
	double m_low[OBJ_NUM],m_high[OBJ_NUM];
	double m_radius;
	double m_acceptdist1;
	double m_acceptdist2;
	double m_GDvarArray[GENS];
	double m_GDobjArray[GENS];
	double m_S[GENS];
	double m_MS[GENS];
	double m_varTrace[GENS][VAR_NUM];
	double m_PopvarTrace[GENS][VAR_NUM];
	double m_Popvar_max_Trace[GENS][VAR_NUM];
	double m_Popvar_min_Trace[GENS][VAR_NUM];
	double m_Pop_Trace[GENS][VAR_NUM];
	int m_gens;
	int count_enter;
	int m_repArray[VAR_NUM];
	bool m_re_eval;
	bool m_change;
	bool m_rand;
	int m_compcount;
	
};

#endif