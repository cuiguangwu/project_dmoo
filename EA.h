//EA.H
//Evolutionary algorithm defination
#ifndef _EA_h_
#define _EA_h_

#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "Random.h"
#include "global.h"

using namespace std;
///////////////////////////////////

// %%%%%%%%%%%%%%% changed here %%%%%%%%%%%%
class IND;
class BinIND;
class POP;
class EA;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class IND {
public:
	enum nRel { SMALLER, EQUAL, BIGGER, INCOMPARABLE };  // numeric comparison.
	typedef void(*Evaluator)(IND& ind); 
	typedef double (*Distance) (const IND&, const IND&);

public:
	IND() { 
		for(int i=0;i<VAR_NUM;i++) m_var[i]=0.0; 
		for(int i=0;i<OBJ_NUM;i++) m_obj[i]=0.0;
		m_gridLoc=0;m_nicheCount=0;m_genNum=0;m_GDvar=0;m_GDobj=0;
	}
	void copy(const IND &other) {
		for(int i=0;i<VAR_NUM;i++) m_var[i]=other.m_var[i]; 
		for(int i=0;i<OBJ_NUM;i++) m_obj[i]=other.m_obj[i];
		m_score=other.m_score;
		m_nicheCount=other.m_nicheCount;
		m_gridLoc=other.m_gridLoc;
		m_genNum=other.m_genNum;
		m_GDvar=other.m_GDvar;
		m_GDobj=other.m_GDobj;
	}
	IND& operator= (const IND& ind) { copy(ind); return *this; }

	nRel compare(const IND& other) const {
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
	/*+++++++++++++++++++ NOTE: != means not comparable. ++++++++++++++++++++++++++++*/
	bool operator!= (const IND &r) { return compare(r)==INCOMPARABLE? true:false; }
	bool operator> (const IND  &r) { return compare(r)==BIGGER      ? true:false; }
	bool operator< (const IND  &r) { return compare(r)==SMALLER     ? true:false; }
	bool operator== (const IND &r) { return compare(r)==EQUAL       ? true:false; }
	bool operator>= (const IND &r) { return *this > r || *this == r ; }
	bool operator<= (const IND &r) { return *this < r || *this == r ; }

	double distance(const IND &other) const {
		double sum = 0.0;
		for(int i=0; i<OBJ_NUM; i++)
		sum += (m_obj[i]-other.m_obj[i]) * (m_obj[i]-other.m_obj[i]);
		return sqrt(sum);
	}

	double distance(const IND &other, double *pLow,double *pHigh) const {
		double range,d;
		double sum = 0.0;
		for(int i=0; i<OBJ_NUM; i++){
			range=pHigh[i]-pLow[i];
			if(range<1.0e-5) range=1.0;
			d=(m_obj[i]-other.m_obj[i])/range;
			sum+=d*d;
		}
		return sqrt(sum);
	}

	void printVariables(ostream& os=cout) const{
		for(int i=0;i<VAR_NUM;i++) os<<m_var[i]<<"     ";
		os<<endl;
	};

	void printObjectives(ostream& os=cout) const{
		for(int i=0;i<OBJ_NUM;i++) os<<m_obj[i]<<"     ";
		os<<endl;
	};

public:
	double m_var[VAR_NUM];
	double m_obj[OBJ_NUM];
	static Evaluator evalFp;
	double m_score;
	double m_GDvar;
	double m_GDobj;
	double m_nicheCount;
	int m_gridLoc;
	int m_genNum;
};

class BinIND : public IND {
public:
	BinIND() : IND() { for(int i=0;i<VAR_NUM*VAR_BIT;i++) m_chrom[i]=0;}
	void copy(const BinIND &other) {
		IND::copy(other); 
		for(int i=0;i<VAR_NUM*VAR_BIT;i++) m_chrom[i]=other.m_chrom[i];
	}
	BinIND& operator= (const BinIND& ind) { copy(ind); return *this; }

	void randomize() {
		for(int i=0;i<VAR_NUM*VAR_BIT;i++) m_chrom[i]=(RND.unit()<0.5)?0:1;
	}
	void decode() {
		//m_chrom: b[0]b[1]...b[VAR_BIT-1]|b[VAR_BIT]b[VAR_BIT+1]...b[2*VAR_BIT-1]|...
		//For the segment of each variable, b[0] is the most significant bit and
		//b[VAR_BIT-1] is the least significant bit.
		int i,j,mul;
		for (i=0; i<VAR_NUM;i++) {
			m_var[i]=0;
			mul = 1;
			for (j = VAR_BIT-1; j >= 0; j--) {
				  m_var[i] += mul * m_chrom[i*VAR_BIT+j];
				  mul *= 2;
			}
			m_var[i]=m_var[i]/pow(2.0,VAR_BIT);
		}
	}
	void grayDecode() {
		//gray code to binary code
		char bin[VAR_NUM*VAR_BIT];
		int i,j;
		for (int i=0; i<VAR_NUM;i++) {
			bin[i*VAR_BIT]=m_chrom[i*VAR_BIT];
			for (j=1;j<=VAR_BIT-1; j++) {
				  bin[i*VAR_BIT+j]=bin[i*VAR_BIT+j-1]^m_chrom[i*VAR_BIT+j];
			}
		}
		//binary code to double
		int mul;
		for (i=0; i<VAR_NUM;i++) {
			m_var[i]=0;
			mul = 1;
			for (j = VAR_BIT-1; j >= 0; j--) {
				  m_var[i] += mul * bin[i*VAR_BIT+j];
				  mul *= 2;
			}
			m_var[i]=m_var[i]/pow(2.0,VAR_BIT);
		}
	}
	void evaluate() { decode();evalFp(*this);}
	void crossover(BinIND &dad, BinIND &sis, BinIND &bro) {
		sis=*this; bro=dad; 
		// uniform crossover
		if( RND.unit() < PC ) { // cross over occurs.
			for(int j=0; j<VAR_NUM*VAR_BIT; j++) 
				if(RND.unit() > 0.5) { 
					sis.m_chrom[j] = dad.m_chrom[j]; bro.m_chrom[j] = m_chrom[j]; 
				}
		}
	}
	void mutate() {
		for(int i=0;i<VAR_NUM*VAR_BIT;i++) if( RND.unit()<PM1 ) m_chrom[i]=1-m_chrom[i];
	}

public:
	char m_chrom[VAR_NUM*VAR_BIT];
};

class POP {
public:
	POP() {m_capacity=0;m_size=0;m_pind=NULL;}
	~POP() { delete[] m_pind;}

	void initialize(int capacity) { 
		m_capacity=capacity;m_size=0;
		delete[] m_pind;m_pind=new BinIND[capacity];
	}
	void add(const BinIND &ind) { if(m_size>=m_capacity) return;m_pind[m_size]=ind;m_size++;}

public:
	BinIND* m_pind;
	int m_size;
	int m_capacity;
};

class EA {
public:
	double m_dLow[OBJ_NUM], m_dHigh[OBJ_NUM];

public:
	EA(){
		for(int i=0;i<OBJ_NUM;i++){ m_dLow[i]=0;m_dHigh[i]=1;}
	}

        //../../PhD/Other-DMOO/Goh/ccea2/COEA.h:20: error: name lookup of ‘i’ changed for ISO ‘for’ scoping
public:
	void objectiveSpaceRange(POP& pop)  
	{                               
		int i,j;
		double aLow[OBJ_NUM], aHigh[OBJ_NUM];
		for (int i=0;i<OBJ_NUM;i++){ aLow[i]=1.0e+10;aHigh[i]=-1.0e+10;}
		for (int i=0;i<pop.m_size;i++){
			for(j=0;j<OBJ_NUM;j++){
				if(aLow[j]>pop.m_pind[i].m_obj[j]) aLow[j]=pop.m_pind[i].m_obj[j];
				if(aHigh[j]<pop.m_pind[i].m_obj[j]) aHigh[j]=pop.m_pind[i].m_obj[j];
			}
		}
		for (int i=0;i<OBJ_NUM;i++){ m_dLow[i]=aLow[i];m_dHigh[i]=aHigh[i];}
	}

	static void objectiveSpaceRange(POP& pop,double *pLow,double *pHigh)  
	{                               
		int i,j;
		double aLow[OBJ_NUM], aHigh[OBJ_NUM];
		for (int i=0;i<OBJ_NUM;i++){ aLow[i]=1.0e+10;aHigh[i]=-1.0e+10;}
		for (int i=0;i<pop.m_size;i++){
			for(j=0;j<OBJ_NUM;j++){
				if(aLow[j]>pop.m_pind[i].m_obj[j]) aLow[j]=pop.m_pind[i].m_obj[j];
				if(aHigh[j]<pop.m_pind[i].m_obj[j]) aHigh[j]=pop.m_pind[i].m_obj[j];
			}
		}
		for (int i=0;i<OBJ_NUM;i++){ pLow[i]=aLow[i];pHigh[i]=aHigh[i];}
	}

	//pre: n is the number of elements in pA and n>1. 
	//	   pass is the number of members that will be sorted.
	//		If pass==-1, all members will be sorted.
	//post: the array is in ascending order. idx is the index.
	static void selectionSort(double *pA,int *idx,int n,int pass=-1)
	{
		int i,j,smallIndex,itmp;
		double temp;
		//Set the index array.
		for (int i=0;i<n;i++) idx[i]=i;
		//The first (pass) members will be sorted.
		if(pass==-1) pass=n;
		//sort pA[0],...,pA[pass]
		for (int i=0;i<pass;i++)
		{
			//start the scan at index i
			smallIndex=i;
			for(j=i+1;j<n;j++) { if(pA[j]<pA[smallIndex]) smallIndex=j;}
			//place the smallest element in pA[i]
			temp=pA[i];pA[i]=pA[smallIndex];pA[smallIndex]=temp;
			itmp=idx[i];idx[i]=idx[smallIndex];idx[smallIndex]=itmp;
		}
	}

	//pre: pop is the solution set. 
	//post: front is the set of nondominated solutions in pop.
	static void paretoFront(const POP &pop,POP &front)
	{
		int i,j;
		bool nondom;
		//clear the front
		front.m_size=0;
		//add nondominated members of pop into front 
		for (i=0;i<pop.m_size;i++){
			nondom=true;
			for (j=0;j<pop.m_size; j++){
				if(pop.m_pind[i]>pop.m_pind[j]){
					nondom=false;
					break;
				}
			}
			if(nondom) front.add(pop.m_pind[i]);
		}
		//
	}

	//Output the pop to a file.
	static void write(char*filename,const POP& pop)
	{
		ofstream fout;
		fout.open(filename,ios::out);
		for (int i = 0; i < pop.m_size; i++) { 
			pop.m_pind[i].printObjectives(fout);
		}
		fout.close();
	}

        //Output the pop to a file.
	static void write2(char*filename,const POP& pop)
	{
		ofstream fout;
		fout.open(filename,ios::app);
		for (int i = 0; i < pop.m_size; i++) {
			pop.m_pind[i].printObjectives(fout);
		}
		fout.close();
	}
};

class Grid
{
public:
	Grid() { pPop=NULL;	}
	~Grid();
	
	void initialize(POP *pp,int objNum, int depth);
	void update_pop(BinIND &ind);
	int find_loc(IND &ind);

private:
	void update_grid();

public:
	int nObjNum;
	int nDepth;
	POP *pPop;
	double *gl_offset,*gl_largest,*gl_range; // the range, offset etc of the grid 
	int *grid_pop; // the array holding the population residing in each grid location
};

#endif

