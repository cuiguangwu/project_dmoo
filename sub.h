#ifndef _SUB_h_
#define _SUB_h_

#include <math.h>
#include <iostream>
#include <stdio.h>
#include "Random.h"
#include "global.h"

///////////////////////////////////
class SubIND;
class SubBinIND;
class SubPOP;

class SubIND {
public:
	enum nRel { SMALLER, EQUAL, BIGGER, INCOMPARABLE };  // numeric comparison.

public:
	SubIND() { 
		m_var=0.0; 
		for(int i=0;i<OBJ_NUM;i++) m_obj[i]=0.0;
		m_score=0;m_gridLoc=0;m_nicheCount=0;
		m_org_pop=0;
	}
	void copy(const SubIND &other) {
		m_var=other.m_var; 
		for(int i=0;i<OBJ_NUM;i++) m_obj[i]=other.m_obj[i];
		m_score=other.m_score;
		m_gridLoc=other.m_gridLoc;m_nicheCount=other.m_nicheCount;
		m_org_pop=other.m_org_pop;
	}
	SubIND& operator= (const SubIND& ind) { copy(ind); return *this; }

	nRel compare(const SubIND& other) const {
		bool bigger = false, smaller = false; 
		bool equal = false, indiff = false;
		for(int i=0; !(indiff) && i<OBJ_NUM; i++) {
			if(m_obj[i] > other.m_obj[i]) bigger  = true;
			if(m_obj[i] < other.m_obj[i]) smaller = true;
			indiff = (bigger && smaller);
		}
		if(indiff) return  INCOMPARABLE;
		if(bigger) return  BIGGER;
		if(smaller)return  SMALLER;
		return             EQUAL;
	}

public:
	double m_var;
	double m_obj[OBJ_NUM];
	double m_score;
	int m_gridLoc;
	double m_nicheCount;
	int m_org_pop;
};

class SubBinIND : public SubIND {
public:
	SubBinIND() : SubIND() { for(int i=0;i<VAR_BIT;i++) m_chrom[i]=0;}
	void copy(const SubBinIND &other) {
		SubIND::copy(other); 
		for(int i=0;i<VAR_BIT;i++) m_chrom[i]=other.m_chrom[i];
	}
	SubBinIND& operator= (const SubBinIND& ind) { copy(ind); return *this; }

	void randomize() {
		for(int i=0;i<VAR_BIT;i++) m_chrom[i]=(RND.unit()<0.5)?0:1;
	}
	void binarycode() {
		double temp;int binary[VAR_BIT];
		temp=m_var*pow(2.0,VAR_BIT);
		for (int j=0; j<VAR_BIT; j++)				
			binary[j]=((unsigned long)(temp) >> j) & 0x1;   //%%%%%%%%%%% changed here - check %%%%%%%%%%%%
		int q=0;
		for(int j = VAR_BIT-1; j>=0; j--)
			m_chrom[j]=binary[q++];
	}
	void decode() {
		int mul = 1;
		m_var=0;
		for (int j = VAR_BIT-1; j >= 0; j--) {
			  m_var += mul * m_chrom[j];
			  mul *= 2;
		}
		m_var=m_var/pow(2.0,VAR_BIT);
	}
	void crossover(SubBinIND &dad, SubBinIND &sis, SubBinIND &bro) {
		sis=*this; bro=dad; 
		//uniform crossover
		if( RND.unit() < PC ) { // cross over occurs.
			for(int j=0; j<VAR_BIT; j++) 
				if(RND.unit() < 0.1) { 
					sis.m_chrom[j] = dad.m_chrom[j]; bro.m_chrom[j] = m_chrom[j]; 
				}
		}
	}
	void onePointCrossover(SubBinIND &dad, SubBinIND &sis, SubBinIND &bro) {
		sis=*this; bro=dad; 
		//single point crossover
		if( RND.unit() < PC ) { // cross over occurs.
			int k=(int)( RND.unit()*VAR_BIT-2)+1;
			for(int j=k; j<VAR_BIT; j++) { 
					sis.m_chrom[j] = dad.m_chrom[j]; 
					bro.m_chrom[j] = m_chrom[j]; 
				}
		}
	}
	void twoPointCrossover(SubBinIND &dad, SubBinIND &sis, SubBinIND &bro) {
		sis=*this; bro=dad; 
		//single point crossover
		if( RND.unit() < PC ) { // cross over occurs.
			int m,n,temp;
			m=(int)( RND.unit()*VAR_BIT);
			n=(int)( RND.unit()*VAR_BIT);
			if(m>n) { temp=m;m=n;n=temp;}
			for(int j=m; j<=n; j++) { 
					sis.m_chrom[j] = dad.m_chrom[j]; 
					bro.m_chrom[j] = m_chrom[j]; 
				}
		}
	}

	void mutate() {
		for(int i=0;i<VAR_BIT;i++) if( RND.unit()<PM1 ) m_chrom[i]=1-m_chrom[i];
	}

	void mutate(double pm) {
		for(int i=0;i<VAR_BIT;i++) if( RND.unit()<pm ) m_chrom[i]=1-m_chrom[i];
	}

public:
	char m_chrom[VAR_BIT];
};

class SubPOP {
public:
	SubPOP() {
		m_capacity=0;m_size=0;m_pind=NULL;
	}
	~SubPOP() { delete[] m_pind;}

	void initialize(int capacity) { 
		m_capacity=capacity;m_size=0;
		delete[] m_pind;
		m_pind=new SubBinIND[capacity];
	}
	void add(SubBinIND &ind) { if(m_size>=m_capacity) return;m_pind[m_size]=ind;m_size++;}
	bool comp(SubBinIND &ind1,SubBinIND &ind2) {
		SubIND::nRel rel;
		rel=ind1.compare(ind2);
		if(rel==SubIND::SMALLER) return true;
		if(rel==SubIND::BIGGER) return false;
		if(ind1.m_score<ind2.m_score) return true;
		if(ind1.m_score==ind2.m_score && ind1.m_nicheCount<ind2.m_nicheCount) return true;
		return false;
	}
	const SubBinIND& tournament(int tour=2) {//select smaller score
		int slct,pos;
		slct=(int)( RND.unit()*m_size);
		for(int k=1;k<tour;k++){
			pos=(int)( RND.unit()*m_size);
			if( comp(m_pind[pos],m_pind[slct]) ) slct=pos;
		}
		return m_pind[slct];
	} 
	void updateRepresentive() {
		int i,j,pos;
		int *tag=new int[m_size];
		for(i=0;i<m_size;i++) tag[i]=0;
		SubBinIND tmp;
		//
		for(i=0;i<REP_SIZE;i++) {
			for(j=0;j<m_size;j++) if(tag[j]==0) {pos=j; tmp=m_pind[j]; break;}
			for(j=0;j<m_size;j++){
				if(tag[j]!=0) continue;
				if( comp(m_pind[j],tmp) ) {
					tmp=m_pind[j]; pos=j;
				}
			}
			tag[pos]=1;
			m_rep[i]=m_pind[pos];
			m_rep_pop=m_pind[pos].m_org_pop;
		}
		//
		delete[] tag;
	}


public:
	SubBinIND* m_pind;
	SubBinIND  m_rep[REP_SIZE]; //representive of this subpopulation
	int m_size;
	int m_capacity;
	int m_rep_pop;
};


#endif

