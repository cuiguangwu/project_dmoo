#include "COEA.h"
#include "statistics.h"
#include <fstream>
#include <math.h>


void COEA::memory1()
{
	/*if(m_change==false){
		POP tempPop;tempPop.initialize(m_pop2.m_size);
		for(int i=0;i<m_pop2.m_size;i++){
			tempPop.add(m_pop2.m_pind[i]);
		}
		for (int i=0;i<m_pop2.m_size;i++){
			m_pop2.m_pind[i].m_genNum=m_gens;
			m_pop2.m_pind[i].evaluate();
		}
		double dist = 0;
		for (int i=0;i<m_pop2.m_size;i++){
			double dist2=tempPop.m_pind[i].distance(m_pop2.m_pind[i]);
			if(dist2>dist)
				dist=dist2;
		}
		m_acceptdist1 = dist;
		//m_change = true;
	}*/
	double *pA=new double [m_pop2.m_size];int *idx=new int [m_pop2.m_size];
	for(int i=0;i<m_pop2.m_size;i++){
		pA[i]=m_pop2.m_pind[i].m_obj[0];
		idx[i]=i;
	}
	EA::selectionSort(pA,idx,m_pop2.m_size);
	
/*	double dist = 100000;
	if(m_change==false){
		for (int i=0;i<RETENT;i++){
			for(int j=0;j<RETENT;j++){
				if(i!=j){
					int sel = (int)(m_pop2.m_size*i/(RETENT-1));
					if(sel!=0)
						sel=sel-1;
					int sel1 = idx[sel];
					sel = (int)(m_pop2.m_size*j/(RETENT-1));
					if(sel!=0)
						sel=sel-1;
					int sel2 = idx[sel];
					double dist2=m_pop2.m_pind[sel1].distance(m_pop2.m_pind[sel2]);
					if(dist2<dist && dist2!=0)
						dist=dist2;
				}
			}
		}
		m_change = true;
		m_acceptdist2 = dist;
	}*/

	if(m_memory.m_size<POP2_SIZE){
		for (int i=0;i<RETENT;i++){
			int sel1 = (int)(m_pop2.m_size*i/(RETENT-1));
			if(sel1!=0)
				sel1=sel1-1;
			int sel = idx[sel1];
			m_memory.add(m_pop2.m_pind[sel]);
		}
	}
	else{
		POP tempPop;tempPop.initialize(POP2_SIZE);
		for (int i=RETENT;i<m_memory.m_size;i++)
			tempPop.add(m_memory.m_pind[i]);
		m_memory.m_size=0;
		for (int i=0;i<tempPop.m_size;i++)
			m_memory.add(tempPop.m_pind[i]);
		for (int i=0;i<RETENT;i++){
			int sel1 = (int)(m_pop2.m_size*i/(RETENT-1));
			if(sel1!=0)
				sel1=sel1-1;
			int sel = idx[sel1];
			m_memory.add(m_pop2.m_pind[sel]);
		}
	}
	POP tempPop;tempPop.initialize(m_memory.m_size);
	for (int i=0;i<m_memory.m_size;i++){
		tempPop.add(m_memory.m_pind[i]);
		tempPop.m_pind[i].m_genNum=m_gens;
		tempPop.m_pind[i].evaluate();
	}
	/*POP tempPop2;tempPop2.initialize(m_memory.m_size);
	for (int i=0;i<tempPop.m_size;i++){
		if( compare_to_archive(tempPop.m_pind[i])==1) {
			tempPop2.add(tempPop.m_pind[i]);
		}
	}
	if(tempPop2.m_size>0){
		for (int i=0;i<tempPop.m_size;i++){
			if( compare_to_archive(tempPop.m_pind[i])!=-1) {
				bool update=true;
				double dist=10000;
				for(int j=0;j<tempPop2.m_size;j++){
					double dist2=tempPop.m_pind[i].distance(tempPop2.m_pind[j]);
					if(dist2<dist)
						dist=dist2;
				}
				if(dist<ADIST*m_acceptdist1 && dist!=0)
					tempPop2.add(tempPop.m_pind[i]);
					//updateArchive(tempPop.m_pind[i]);
			}
		}
	}*/
	m_pop2.m_size=0;
	for (int i=0;i<tempPop.m_size;i++)
		updateArchive(tempPop.m_pind[i]);

	delete[]pA;
	delete[]idx;
}

void COEA::memory2()
{
	
	POP tempPop;tempPop.initialize(m_pop2.m_size);
	for(int i=0;i<m_pop2.m_size;i++){
		tempPop.add(m_pop2.m_pind[i]);
		tempPop.m_pind[i].m_genNum=m_gens;
		tempPop.m_pind[i].evaluate();
	}

	m_pop2.m_size=0;
	for (int i=0;i<tempPop.m_size;i++)
		updateArchive(tempPop.m_pind[i]);
}

void COEA::S(POP &pop)
{
	double *arr=new double [pop.m_size];
	if(pop.m_size<3){
		m_S[m_gens]=1e7;
		return;
	}
	else {
		for(int i=0;i<pop.m_size;i++){
			double dist=1e7;
			for(int j=0;j<pop.m_size;j++){
				double tmp = 0;
				if(i!=j && !(pop.m_pind[i]==pop.m_pind[j])){
					for(int k=0;k<OBJ_NUM;k++){
						tmp = tmp + pow((pop.m_pind[i].m_obj[k]-pop.m_pind[j].m_obj[k]),2);
						if(tmp==0)
							int yy=0;
					}
					tmp = sqrt(tmp);
					if(tmp<dist) dist = tmp;
				}
			}
			arr[i]=dist;
		}
	}
	double	dt = mean(arr,pop.m_size);
	double delta = sqrt(variance(arr,pop.m_size));
	double s = delta/dt;
	m_S[m_gens] = s/pop.m_size;
	delete []arr;
}

void COEA::MS(POP &pop1)
{
	double ms=0;
	double *prange=new double[OBJ_NUM];
	double *ratio=new double[OBJ_NUM];
	double *range=new double[OBJ_NUM];
	double *max=new double[OBJ_NUM];double *min=new double[OBJ_NUM];
	for(int i=0;i<OBJ_NUM;i++){
		prange[i] = 1-0;
		max[i] = -1e7;
		min[i] = 1e7;
	}
	for(int i=0;i<pop1.m_size;i++){
		for(int j=0;j<OBJ_NUM;j++){
			if(pop1.m_pind[i].m_obj[j]<min[j])
				min[j]=pop1.m_pind[i].m_obj[j];
			if(pop1.m_pind[i].m_obj[j]>max[j])
				max[j]=pop1.m_pind[i].m_obj[j];
		}
	}

	for (int i=0;i<OBJ_NUM;i++){
		if(0>min[i])
			min[i]=0.0;
		if(1<max[i])
			max[i]=1.0;
	}
	for (int i=0;i<OBJ_NUM;i++){
		range[i] = max[i]-min[i];
		if(range[i]<0)
			range[i]=0;
		ratio[i] = range[i]/prange[i];
		ms = ms + ratio[i]*ratio[i];
	}
	ms = sqrt(ms)/sqrt((double)OBJ_NUM);
	m_MS[m_gens]=ms;

	delete []prange;
	delete []ratio;
	delete []range;
	delete []max;delete []min;
}

void COEA::updateArchive(BinIND &ind)
{
	// given a solution s, add it to the archive if
	// a) the archive is empty 
	// b) the archive is not full and s is not dominated or equal to anything currently in the archive
	// c) s dominates anything in the archive                        
	// d) the archive is full but s is nondominated and is in a no more crowded square than at least one solution
	// in addition, maintain the archive such that all solutions are nondominated. 
	
	POP *pPop=&m_pop2;
	int i,k;
	//case a)
	if (pPop->m_size==0) { 
		pPop->add(ind);
		return; 
	}
	//
	int *tag=new int[pPop->m_size];
	for (int i=0;i<pPop->m_size;i++) tag[i]=0;
	int join = 0;	//if the ind comes into the pop
	int set = 0;	//if delete members of the pop
	IND::nRel result = IND::INCOMPARABLE;//if the ind dominates the pop
	for (int i=0;i<pPop->m_size;i++)
	{
		result = ind.compare(pPop->m_pind[i]);
		if (result==IND::EQUAL || result==IND::BIGGER) break;
		if ( (result==IND::SMALLER)&&(join==0) ) { pPop->m_pind[i]=ind; join=1;}
		else if (result==IND::SMALLER) { tag[i]=1;set=1;}	    	    
	}
	//
	if (result==IND::EQUAL || result==IND::BIGGER) {
		
	}
	//case c)
	if (set==1)
	{
		int old_size = pPop->m_size;
		POP tmp;
		tmp.initialize(old_size);
		for (int i=0;i<old_size;i++) tmp.m_pind[i] = pPop->m_pind[i];
		pPop->m_size= 0;
		for (int i=0;i<old_size;i++) if (tag[i]!=1) pPop->add(tmp.m_pind[i]);
	}

	if ( (join==0)&&(result!=IND::EQUAL)&&(result!=IND::BIGGER) )  // ie solution is non-dominated by the list
	{
		//case d)
		if (pPop->m_size >= pPop->m_capacity)  
		{	  
			double nc;
			double most = -10000;
			int repl;
			for (i = 0; i < pPop->m_size; i++)
			{
				nc=niche(pPop->m_pind[i]);
				if (nc > most)
				{
					most = nc;
					repl = i;
				}
			}
			pPop->m_pind[repl] = ind;
		}
		//case b)
		else 
			pPop->add(ind);
	}
	//
	delete[] tag;
	//
	for(k=0;k<OBJ_NUM;k++) {
		m_low[k]=m_pop2.m_pind[0].m_obj[k];
		m_high[k]=m_pop2.m_pind[0].m_obj[k];
	}
	for (int i=1;i<m_pop2.m_size;i++) {
		for(k=0;k<OBJ_NUM;k++){
			if(m_low[k]>m_pop2.m_pind[i].m_obj[k]) m_low[k]=m_pop2.m_pind[i].m_obj[k];
			if(m_high[k]<m_pop2.m_pind[i].m_obj[k]) m_high[k]=m_pop2.m_pind[i].m_obj[k];
		}
	}
	//
	return;
}

// compares a solution to every member of the archive. 
//Post: Returns -1 if dominated by any member, 1 if dominates any member, and 0 otherwise
int COEA::compare_to_archive(const IND &ind)  
{                               
	int i=0;
	int result=0;
	while((i<m_pop2.m_size)&&(result!=1)&&(result!=-1))
    {
		IND::nRel rel;
		rel=ind.compare(m_pop2.m_pind[i]);
		if ( rel==IND::SMALLER) result=1;
		if ( rel==IND::BIGGER) result=-1;
		i++;
    }
	if(m_pop2.m_size==0) result=1;
	//
	return(result);
}

int COEA::rank(IND &ind)  
{                               
	int ret=1;
	for(int i=0;i<m_pop2.m_size;i++)
	{
		if(m_pop2.m_pind[i]<ind) ret++;
	}
	return ret;
}

double COEA::dynamicDistance()
{
	if(m_pop2.m_size<2) return 1.0/POP2_SIZE;
	//
	int i,j,k,ix,iy;
	double dist,dmin,dmax,d1,d2,d,radius;
	// Find the range of objective space
	for (int i=1;i<m_pop2.m_size;i++) {
		for(k=0;k<OBJ_NUM;k++){
			if(m_low[k]>m_pop2.m_pind[i].m_obj[k]) m_low[k]=m_pop2.m_pind[i].m_obj[k];
			if(m_high[k]<m_pop2.m_pind[i].m_obj[k]) m_high[k]=m_pop2.m_pind[i].m_obj[k];
		}
	}
	//dy
	dmin=0.0;
	for (int i=0;i<m_pop2.m_size;i++)
	{
		for(j=0;j<m_pop2.m_size;j++)
		{
			dist=m_pop2.m_pind[i].distance(m_pop2.m_pind[j],m_low,m_high);
			if(dist>dmin) { dmin=dist;ix=i;iy=j;}
		}
	}
	//
	IND ind;
	for(k=0;k<OBJ_NUM;k++)
	{
		ind.m_obj[k]=m_pop2.m_pind[ix].m_obj[k];
		if(ind.m_obj[k]>m_pop2.m_pind[iy].m_obj[k])
		{
			ind.m_obj[k]=m_pop2.m_pind[iy].m_obj[k];
		}
	}
	d1=ind.distance(m_pop2.m_pind[ix],m_low,m_high);
	d2=ind.distance(m_pop2.m_pind[iy],m_low,m_high);
	dmax=d1+d2;
	//
	d=(dmin+dmax)/2.0;
	//
	double psize,m;
	psize=POP2_SIZE;m=OBJ_NUM;
	radius=(d/2)/pow(psize,1/(m-1));
	m_radius=radius;
	return radius;
}

double COEA::niche(IND &ind)  
{                               
	double dist,radius,ret;
	radius=2*m_radius;
	ret=0.0;
	for(int i=0;i<m_pop2.m_size;i++)
	{
		//dist=ind.distance(m_pop2.m_pind[i]);
		dist=ind.distance(m_pop2.m_pind[i],m_low,m_high);
		if(dist<radius) ret=ret+(1-dist/radius)*(1-dist/radius);
		//if(dist<radius) ret=ret+(1-dist/radius);
	}
	return ret;
}

//Comparation operator of CCEA
bool COEA::comp(IND &ind1,IND &ind2) {
	IND::nRel rel;
	rel=ind1.compare(ind2);
	if(rel==IND::SMALLER) return true;
	if(rel==IND::BIGGER) return false;
	ind1.m_score=rank(ind1);ind2.m_score=rank(ind2);
	if(ind1.m_score<ind2.m_score) return true;
	ind1.m_nicheCount=niche(ind1);ind2.m_nicheCount=niche(ind2);
	if(ind1.m_score==ind2.m_score && ind1.m_nicheCount<ind2.m_nicheCount) return true;
	return false;
}

const BinIND& COEA::tournament(int tour) {//select smaller score
	int slct,pos;
	slct=(int)( RND.unit()*m_pop2.m_size);
	for(int k=1;k<tour;k++){
		pos=(int)( RND.unit()*m_pop2.m_size);
		if( comp(m_pop2.m_pind[pos],m_pop2.m_pind[slct]) ) slct=pos;
	}
	return m_pop2.m_pind[slct];
} 

//Pre: id the number of subpopulation
//Assign fitness to every member of subpopulation id.
void COEA::evaluate(int id)
{
	BinIND ind,best;
	int i,j,k,m,gene,rnd;
	//estimate the score of subpop[id]
	for (int i=0;i<m_subpop[id].m_size;i++) {
		//set initial values for best
		for(j=0;j<OBJ_NUM;j++) best.m_obj[j]=1.0e+10;
		//multiple tests
		for(m=0;m<1;m++){
			//assemble partial solution into complete solution
			gene=0;
			for(k=0;k<VAR_NUM;k++){
				//*
				rnd=0;
				/*if( m>0 ) rnd=(int)( RND.unit()*REP_SIZE );
				for(j=0;j<VAR_BIT;j++) 
					ind.m_chrom[gene++]=m_subpop[k].m_rep[rnd].m_chrom[j];*/
				//
				//if( m==0 ) {
					rnd=0;
					for(j=0;j<VAR_BIT;j++) 
						ind.m_chrom[gene++]=m_subpop[k].m_rep[rnd].m_chrom[j];
			//	}
				/*if( m>0 ) {
					rnd=(int)( RND.unit()*POP1_SIZE );
					for(j=0;j<VAR_BIT;j++) 
					ind.m_chrom[gene++]=m_subpop[k].m_pind[rnd].m_chrom[j];
				}*/
			}
			for(j=0;j<VAR_BIT;j++) 
				ind.m_chrom[id*VAR_BIT+j]=m_subpop[id].m_pind[i].m_chrom[j];
			ind.m_genNum=m_gens;
			ind.evaluate();
			if( compare_to_archive(ind)!=-1 ) 
				updateArchive(ind);
			//
			//if( comp(ind,best) ) best=ind;
			if( ind<best ) best=ind;
			//if( !(ind>best) ) best=ind;
		}
		//end of multiple tests
		best.m_score=rank(best);best.m_nicheCount=niche(best);
		for(j=0;j<OBJ_NUM;j++) m_subpop[id].m_pind[i].m_obj[j]=best.m_obj[j];
		m_subpop[id].m_pind[i].m_var=best.m_var[id];
		m_subpop[id].m_pind[i].m_score=best.m_score;
		m_subpop[id].m_pind[i].m_nicheCount=best.m_nicheCount;
	}
}

void COEA::evaluateC(int id,int*rnd)
{
	BinIND ind,best;
	int i,j,k,m,gene;
	
	//estimate the score of subpop[id]
	for (int i=0;i<m_subpop[id].m_size;i++) {
		//set initial values for best
		for(j=0;j<OBJ_NUM;j++) best.m_obj[j]=1.0e+10;
		//multiple tests
		for(m=0;m<1;m++){
			//assemble partial solution into complete solution
			gene=0;
			for(k=0;k<VAR_NUM;k++){
				//*
				
				for(j=0;j<VAR_BIT;j++) 
					ind.m_chrom[gene++]=m_subpop[k].m_pind[rnd[k]].m_chrom[j];
		
			}
			for(j=0;j<VAR_BIT;j++) 
				ind.m_chrom[id*VAR_BIT+j]=m_subpop[id].m_pind[i].m_chrom[j];
			ind.m_genNum=m_gens;
			ind.evaluate();
			if( compare_to_archive(ind)!=-1 ) 
				updateArchive(ind);
			//
			//if( comp(ind,best) ) best=ind;
			if( ind<best ) best=ind;
			//if( !(ind>best) ) best=ind;
		}
		//end of multiple tests
		best.m_score=rank(best);best.m_nicheCount=niche(best);
		for(j=0;j<OBJ_NUM;j++) m_subpop[id].m_pind[i].m_obj[j]=best.m_obj[j];
		m_subpop[id].m_pind[i].m_var=best.m_var[id];
		m_subpop[id].m_pind[i].m_score=best.m_score;
		m_subpop[id].m_pind[i].m_nicheCount=best.m_nicheCount;
	}
}

void COEA::evaluateE(int id)
{
	BinIND ind,best;
	int i,j,k,m,gene,rnd[VAR_NUM], type=0;
	//estimate the score of subpop[id]
	for (int i=0;i<VAR_NUM;i++){
		if(RND.unit()<0.5)
			rnd[i]=(int)( RND.unit()*POP1_SIZE);
		else
			rnd[i]=-1;
	}

	for (int i=0;i<m_subpop[id].m_size;i++) {
		//set initial values for best
		for(j=0;j<OBJ_NUM;j++) best.m_obj[j]=1.0e+10;
		//multiple tests
		for(m=0;m<1;m++){
			//assemble partial solution into complete solution
			gene=0;
			for(k=0;k<VAR_NUM;k++){
				//*
				if(rnd[k]!=-1){
					for(j=0;j<VAR_BIT;j++) 
						ind.m_chrom[gene++]=m_subpop[k].m_pind[rnd[k]].m_chrom[j];
				}
				else{
					for(j=0;j<VAR_BIT;j++) 
						ind.m_chrom[gene++]=m_subpop[k].m_rep[0].m_chrom[j];
				}
		
			}
			for(j=0;j<VAR_BIT;j++) 
				ind.m_chrom[id*VAR_BIT+j]=m_subpop[id].m_pind[i].m_chrom[j];
			ind.m_genNum=m_gens;
			ind.evaluate();
			if( compare_to_archive(ind)!=-1 ) 
				updateArchive(ind);
			if( ind<best ) best=ind;
		}
		//end of multiple tests
		best.m_score=rank(best);best.m_nicheCount=niche(best);
		for(j=0;j<OBJ_NUM;j++) m_subpop[id].m_pind[i].m_obj[j]=best.m_obj[j];
		m_subpop[id].m_pind[i].m_var=best.m_var[id];
		m_subpop[id].m_pind[i].m_score=best.m_score;
		m_subpop[id].m_pind[i].m_nicheCount=best.m_nicheCount;
	}
}

void COEA::evaluateE4(int id, int*rnd)
{
	BinIND ind,best;
	int i,j,k,m,gene, type=0;
	//estimate the score of subpop[id]

	for (int i=0;i<m_subpop[id].m_size;i++) {
		//set initial values for best
		for(j=0;j<OBJ_NUM;j++) best.m_obj[j]=1.0e+10;
		//multiple tests
		for(m=0;m<1;m++){
			//assemble partial solution into complete solution
			gene=0;
			for(k=0;k<VAR_NUM;k++){
				//*
				if(rnd[k]!=-1){
					for(j=0;j<VAR_BIT;j++) 
						ind.m_chrom[gene++]=m_subpop[k].m_pind[rnd[k]].m_chrom[j];
				}
				else{
					for(j=0;j<VAR_BIT;j++) 
						ind.m_chrom[gene++]=m_subpop[k].m_rep[0].m_chrom[j];
				}
		
			}
			for(j=0;j<VAR_BIT;j++) 
				ind.m_chrom[id*VAR_BIT+j]=m_subpop[id].m_pind[i].m_chrom[j];
			ind.m_genNum=m_gens;
			ind.evaluate();
			if( compare_to_archive(ind)!=-1 ) 
				updateArchive(ind);
			if( ind<best ) best=ind;
		}
		//end of multiple tests
		best.m_score=rank(best);best.m_nicheCount=niche(best);
		for(j=0;j<OBJ_NUM;j++) m_subpop[id].m_pind[i].m_obj[j]=best.m_obj[j];
		m_subpop[id].m_pind[i].m_var=best.m_var[id];
		m_subpop[id].m_pind[i].m_score=best.m_score;
		m_subpop[id].m_pind[i].m_nicheCount=best.m_nicheCount;
	}
}

void COEA::evaluateD(int id, int *idx, int pt)
{
	BinIND ind,best;
	int i,j,k,m,gene,rnd;
	//estimate the score of subpop[id]
	for (int i=0;i<m_subpop[id].m_size;i++) {
		//set initial values for best
		for(j=0;j<OBJ_NUM;j++) best.m_obj[j]=1.0e+10;
		//multiple tests
		for(m=0;m<1;m++){
			//assemble partial solution into complete solution
			gene=0;
			for(k=0;k<VAR_NUM;k++){
				bool support = false;
				for(int h=0;h<pt;h++){
					if(k==idx[h])
						support=true;
				}
				if(support){
					for(j=0;j<VAR_BIT;j++) 
						ind.m_chrom[gene++]=m_subpop[k].m_rep[0].m_chrom[j];
				}
				else{
					rnd=(int)( RND.unit()*POP1_SIZE);
					for(j=0;j<VAR_BIT;j++) 
						ind.m_chrom[gene++]=m_subpop[k].m_pind[rnd].m_chrom[j];
				}
			
			}
			for(j=0;j<VAR_BIT;j++) 
				ind.m_chrom[id*VAR_BIT+j]=m_subpop[id].m_pind[i].m_chrom[j];
			ind.m_genNum=m_gens;
			ind.evaluate();
			if( compare_to_archive(ind)!=-1 ) 
				updateArchive(ind);
			//
			//if( comp(ind,best) ) best=ind;
			if( ind<best ) best=ind;
			//if( !(ind>best) ) best=ind;
		}
		//end of multiple tests
		best.m_score=rank(best);best.m_nicheCount=niche(best);
		for(j=0;j<OBJ_NUM;j++) m_subpop[id].m_pind[i].m_obj[j]=best.m_obj[j];
		m_subpop[id].m_pind[i].m_var=best.m_var[id];
		m_subpop[id].m_pind[i].m_score=best.m_score;
		m_subpop[id].m_pind[i].m_nicheCount=best.m_nicheCount;
	}
}
/////////////////////////////////////////////////////////////////////
//Implicit Competition
/////////////////////////////////////////////////////////////////////
void COEA::implicitCompetition(){
	//m_re_eval=false;
	

	int com_degree = (int)(POP1_SIZE*CR1);
	SubPOP sp[VAR_NUM]; 
	for(int i=0;i<VAR_NUM;i++) {
		sp[i].initialize(POP1_SIZE); 
		for(int j=0;j<POP1_SIZE;j++)
			sp[i].add(m_subpop[i].m_pind[j]);
		
	}
	if(VAR_NUM<=com_degree){
		for(int i=0;i<VAR_NUM;i++){
			int count=POP1_SIZE-1;
			for(int j=0;j<VAR_NUM;j++){
				if(j!=i){
					int rnd=(int)( RND.unit()*REP_SIZE );
					m_subpop[i].m_pind[count]=m_subpop[j].m_rep[0];
					m_subpop[i].m_pind[count].m_org_pop = j;
					count--;
				}
			}
			int countk=1;
			for(int j=count;j>=1;j--){
				double k = count-1;
				double range = 1.0/k;
				double val = RND.unit()*(range*countk-range*(countk-1)) + range*(countk-1);
				m_subpop[i].m_pind[j].m_var=val;
				m_subpop[i].m_pind[j].binarycode();
				//m_subpop[i].m_pind[j].decode();
				m_subpop[i].m_pind[j].m_org_pop = -1;
				countk++;
			}
			for(int j=0;j<1;j++){
				m_subpop[i].m_pind[j] = m_subpop[i].m_rep[0];
				m_subpop[i].m_pind[j].m_org_pop = i;
			}
			m_subpop[i].m_rep_pop = i;
		}
	}
	else {
		for(int i=0;i<VAR_NUM;i++){
			double array[VAR_NUM];
			for(int j=0;j<VAR_NUM;j++){
				if(j!=i)
					array[j]=1;
				else
					array[j]=0;
			}
			for(int j=POP1_SIZE-1;j>=POP1_SIZE-com_degree;j--){
				int sel=i;
				while(array[sel]==0){
					sel=(int)( RND.unit()*VAR_NUM);
				}
				array[sel]=0;
				int rnd=(int)( RND.unit()*REP_SIZE );
				m_subpop[i].m_pind[j]=m_subpop[sel].m_rep[rnd];
				m_subpop[i].m_pind[j].m_org_pop = sel;
			}
			int countk=1;
			for(int j=POP1_SIZE-com_degree-1;j>=1;j--){
				double k = POP1_SIZE-com_degree-1;
				double range = 1.0/k;
				double val = RND.unit()*(range*countk-range*(countk-1)) + range*(countk-1);
				m_subpop[i].m_pind[j].m_var=val;
				m_subpop[i].m_pind[j].binarycode();
				//m_subpop[i].m_pind[j].decode();
				m_subpop[i].m_pind[j].m_org_pop = -1;
				countk++;
			}
			for(int j=0;j<1;j++){
				m_subpop[i].m_pind[j] = m_subpop[i].m_rep[0];
				m_subpop[i].m_pind[j].m_org_pop = i;
			}
			m_subpop[i].m_rep_pop = i;
		}
	}

	double array[VAR_NUM];
	int index[VAR_NUM];
		
	for (int i=0; i<VAR_NUM;i++){
		array[i] = RND.unit();
		index[i] = i;
	}
		
	/*int  smallIndex, itmp;
	double temp;
	
	for (int i=0;i<VAR_NUM;i++){
		smallIndex=i;
		for(int j=i+1;j<VAR_NUM;j++) { 
			if(array[j]<array[smallIndex]) smallIndex=j;
		}
		temp=array[i];array[i]=array[smallIndex];array[smallIndex]=temp;
		itmp=index[i];index[i]=index[smallIndex];index[smallIndex]=itmp;
	}*/
	
	int rnd[VAR_NUM];
	for (int i=0;i<VAR_NUM;i++){
		if(RND.unit()<0.5)
			rnd[i]=(int)( RND.unit()*POP1_SIZE);
		else
			rnd[i]=-1;
	}
	for (int i=0;i<VAR_NUM;i++){
		//evaluateD(index[i], index, i);
		//m_subpop[index[i]].updateRepresentive();
		/*if(m_rand){
			evaluateC(i);
			m_rand=false;
		}
		else{
			evaluate(i);
			m_rand=true;
		}*/
		//evaluate(i);
		//evaluateC(i,rnd);
		evaluateE4(i,rnd);
		m_subpop[i].updateRepresentive();
	}
	for (int i=0;i<VAR_NUM;i++){
		if(m_subpop[i].m_rep_pop!=-1){
			m_subpop[i].m_size=0;
			for(int j=0;j<10;j++){
				m_subpop[i].add(sp[m_subpop[i].m_rep_pop].m_pind[j]);
				m_subpop[i].m_pind[j].m_org_pop=i;
			}
			m_Pop_Trace[m_compcount][i]=m_subpop[i].m_rep_pop;
		}
		else{
			int countk=1;
			double k = POP1_SIZE-com_degree-1;
			double range = 1.0/k;
			double rangeII=0;
			for(int j=POP1_SIZE-com_degree-1;j>=1;j--){
				if(m_subpop[i].m_rep[0].m_var<range*countk){
					rangeII=range*countk;
					break;
				}
				countk++;
			}
			for(int j=0;j<10;j++){
				m_subpop[i].m_pind[j].m_var = RND.uniform(rangeII-range-(0.5*range), rangeII+(0.5*range));
				if(m_subpop[i].m_pind[j].m_var>1)
					m_subpop[i].m_pind[j].m_var=1;
				if(m_subpop[i].m_pind[j].m_var<0)
					m_subpop[i].m_pind[j].m_var=0;
				m_subpop[i].m_pind[j].binarycode();
				m_subpop[i].m_pind[j].m_org_pop=i;
			}
		}
		
		m_subpop[i].m_rep_pop=i;
		m_subpop[i].m_pind[0]=m_subpop[i].m_rep[0];
	}
	m_compcount++;

}

void COEA::implicitCompetitionNoRand(){
	
	int com_degree = POP1_SIZE-1;//(int)(POP1_SIZE*CR);//
	SubPOP sp[VAR_NUM]; 
	for(int i=0;i<VAR_NUM;i++) {
		sp[i].initialize(POP1_SIZE); 
		for(int j=0;j<POP1_SIZE;j++)
			sp[i].add(m_subpop[i].m_pind[j]);
		
	}
	if(VAR_NUM<=com_degree){
		for(int i=0;i<VAR_NUM;i++){
			int count=POP1_SIZE-1;
			for(int j=0;j<VAR_NUM;j++){
				if(j!=i){
					//int rnd=(int)( RND.unit()*REP_SIZE );
					m_subpop[i].m_pind[count]=m_subpop[j].m_rep[0];
					m_subpop[i].m_pind[count].m_org_pop = j;
					count--;
				}
			}
			for(int j=0;j<com_degree-VAR_NUM+1;j++){
				int rnd = i;
				while(rnd==i){
					rnd = (int)(RND.unit()*VAR_NUM);
				}
				int sel = (int)(RND.unit()*POP1_SIZE);
				m_subpop[i].m_pind[count]=m_subpop[rnd].m_pind[sel];
				m_subpop[i].m_pind[count].m_org_pop = rnd;
				count--;
			}
			for(int j=0;j<=count;j++){
				m_subpop[i].m_pind[j].m_org_pop = i;
			}
			//m_subpop[i].m_pind[j] = m_subpop[i].m_rep[0];
			m_subpop[i].m_rep_pop=i;
			m_subpop[i].m_pind[0]=m_subpop[i].m_rep[0];
		}
	}
	else {
		for(int i=0;i<VAR_NUM;i++){
			double array[VAR_NUM];
			for(int j=0;j<VAR_NUM;j++){
				if(j!=i)
					array[j]=1;
				else
					array[j]=0;
			}
			for(int j=POP1_SIZE-1;j>=POP1_SIZE-com_degree;j--){
				int sel=i;
				while(array[sel]==0){
					sel=(int)( RND.unit()*VAR_NUM);
				}
				array[sel]=0;
				int rnd=(int)( RND.unit()*REP_SIZE );
				m_subpop[i].m_pind[j]=m_subpop[sel].m_rep[0];
				m_subpop[i].m_pind[j].m_org_pop = sel;
			}

			for(int j=POP1_SIZE-com_degree-1;j>=0;j--){
				m_subpop[i].m_pind[j].m_org_pop = i;
			}
			m_subpop[i].m_pind[0] = m_subpop[i].m_rep[0];
			m_subpop[i].m_rep_pop = i;
		}
	}

	double array[VAR_NUM];
	int index[VAR_NUM];
		
	for (int i=0; i<VAR_NUM;i++){
		array[i] = RND.unit();
		index[i] = i;
	}
		
	int  smallIndex, itmp;
	double temp;
	
	for (int i=0;i<VAR_NUM;i++){
		smallIndex=i;
		for(int j=i+1;j<VAR_NUM;j++) { 
			if(array[j]<array[smallIndex]) smallIndex=j;
		}
		temp=array[i];array[i]=array[smallIndex];array[smallIndex]=temp;
		itmp=index[i];index[i]=index[smallIndex];index[smallIndex]=itmp;
	}
	int rnd[VAR_NUM];
	for (int i=0;i<VAR_NUM;i++){
		//if(RND.unit()<0.5)
			rnd[i]=(int)( RND.unit()*POP1_SIZE);
		//else
		//	rnd[i]=-1;
	}
	for (int i=0;i<VAR_NUM;i++){
		//evaluateD(index[i], index, i);
		//m_subpop[index[i]].updateRepresentive();
		/*if(m_rand){
			evaluateC(i);
			m_rand=false;
		}
		else{
			evaluate(i);
			m_rand=true;
		}*/
		//evaluate(i);
		evaluateC(i,rnd);
		//evaluateE4(i,rnd);
		m_subpop[i].updateRepresentive();
	}

	for (int i=0;i<VAR_NUM;i++){
		m_subpop[i].m_size=0;
		for(int j=0;j<10;j++){//POP1_SIZE
			m_subpop[i].add(sp[m_subpop[i].m_rep_pop].m_pind[j]);
			m_subpop[i].m_pind[j].m_org_pop=i;
		}
		/*for( j=0;j<5;j++){
			m_subpop[i].add(sp[m_subpop[i].m_rep_pop].m_pind[j]);
			m_subpop[i].m_pind[4+j].randomize();
			m_subpop[i].m_pind[4+j].m_org_pop=i;
		}*/
		m_Pop_Trace[m_compcount][i]=m_subpop[i].m_rep_pop;
		m_subpop[i].m_pind[0]=m_subpop[i].m_rep[0];
		m_subpop[i].m_rep_pop=i;
		//dynamicDistance();
	}
	m_compcount++;
}

void COEA::reinit(){
	
	for(int i=0;i<VAR_NUM;i++){
		for(int j=0;j<POP1_SIZE;j++){
			m_subpop[i].m_pind[j].randomize();
		}
	//}  // %%%%%%%%%%%%%%%%%% changed her - check
                m_subpop[i].m_pind[0]=m_subpop[i].m_rep[0];
        }
	for (int i=0;i<VAR_NUM;i++){
		evaluate(i);
		m_subpop[i].updateRepresentive();
		dynamicDistance();
	}
}

//Extending operator
//It clones the archive member with the least niche count into subpopulations.
void COEA::extend(int nClone)
{
	if(nClone<1) return;
	//if(m_pop2.m_size<m_pop2.m_capacity) return;
	if(m_pop2.m_size==0) return;
	//find the member with lowest nc
	double nc,min;
	int pose;
	min=1.0e+10;
	for (int i = 0; i < m_pop2.m_size; i++)
	{
		nc=niche(m_pop2.m_pind[i]);
		if (nc < min) { min = nc; pose=i;}
	}
	//
	double *pA=new double [m_pop2.m_size];int *idx=new int [m_pop2.m_size];
	int objsel=(int)RND.unit()*OBJ_NUM;
	for (int i=0;i<m_pop2.m_size;i++){
		pA[i]=m_pop2.m_pind[i].m_obj[objsel];
		idx[i]=i;
	}
	EA::selectionSort(pA,idx,m_pop2.m_size);
	
	int pos[2];pos[0]=idx[0];pos[1]=idx[m_pop2.m_size-1];
	for (int i=0;i<2;i++){
		int gene=0;
		int tail = POP1_SIZE-i-1;
		for(int k=0;k<VAR_NUM;k++) {
			for(int j=0;j<VAR_BIT;j++)
				m_subpop[k].m_pind[tail].m_chrom[j]=m_pop2.m_pind[pos[i]].m_chrom[gene++];
		}
	}
	delete[]pA;
	delete[]idx;
	int tail,count,gene;
	count=0;
	
	gene=0;
	for(int k=0;k<VAR_NUM;k++) 
	{
		tail=m_subpop[k].m_size-2;
		for (int i=0;i<VAR_BIT;i++)
			m_subpop[k].m_pind[tail].m_chrom[i]=m_pop2.m_pind[pose].m_chrom[gene++];
			
	}
	if(tail>count) count++;
	
}

void COEA::trackingDynamicsI()
{
	m_re_eval = false;
	if(m_pop2.m_size>0){
		for(int i=0;i<SAMPLE;i++){
			int sel = (int)(RND.unit()*m_pop2.m_size);
			BinIND ind;
			ind = m_pop2.m_pind[sel];
			ind.m_genNum=m_gens;
			ind.evaluate();
			for(int j=0;j<OBJ_NUM;j++){
				if(m_pop2.m_pind[sel].m_obj[j]!=ind.m_obj[j])
					m_re_eval=true;
			}
			if(m_re_eval)
				break;
		}
	}
}

void COEA::trackingDynamicsII()
{
	for(int i=0;i<VAR_NUM;i++){
		double maxVar =0,minVar = 999;
		for(int j=0;j<m_subpop[m_repArray[i]].m_size;j++){
			m_subpop[m_repArray[i]].m_pind[j].decode();
			m_PopvarTrace[m_gens][i]=m_PopvarTrace[m_gens][i]+ m_subpop[m_repArray[i]].m_pind[j].m_var/(double)m_subpop[m_repArray[i]].m_size;
			if(m_subpop[m_repArray[i]].m_pind[j].m_var>maxVar){
				maxVar=m_subpop[m_repArray[i]].m_pind[j].m_var;
				m_Popvar_max_Trace[m_gens][i]=maxVar;
			}
			if(m_subpop[m_repArray[i]].m_pind[j].m_var<minVar){
				minVar=m_subpop[m_repArray[i]].m_pind[j].m_var;
				m_Popvar_min_Trace[m_gens][i]=minVar;
			}
		}
	}
	
	POP temp;temp.initialize(m_pop2.m_size);
	for (int i=0;i<m_pop2.m_size;i++){
		temp.add(m_pop2.m_pind[i]);
	}
	for (int i=0;i<temp.m_size;i++){
		temp.m_pind[i].m_genNum=m_gens;
		temp.m_pind[i].evaluate();
	}
	double *arr1=new double [m_pop2.m_size];
	double *arr2=new double [m_pop2.m_size];
	if(m_re_eval){
		m_re_eval=false;
		m_pop2.m_size=0;
		for (int i=0;i<temp.m_size;i++){
			if( compare_to_archive(temp.m_pind[i])!=-1 ) 
				updateArchive(temp.m_pind[i]);
		}
		for (int i=0;i<m_pop2.m_size;i++){
			arr1[i]=m_pop2.m_pind[i].m_GDvar;
			arr2[i]=m_pop2.m_pind[i].m_GDobj;
			for(int j=0;j<VAR_NUM;j++){
				m_varTrace[m_gens][j]=m_varTrace[m_gens][j]+m_pop2.m_pind[i].m_var[j]/(double)m_pop2.m_size;
			}
		}
		m_GDvarArray[m_gens] = mean(arr1, m_pop2.m_size);
		m_GDobjArray[m_gens] = mean(arr2, m_pop2.m_size);
	}
	else {
		for (int i=0;i<temp.m_size;i++){
			arr1[i]=temp.m_pind[i].m_GDvar;
			arr2[i]=temp.m_pind[i].m_GDobj;
			for(int j=0;j<VAR_NUM;j++){
				m_varTrace[m_gens][j]=m_varTrace[m_gens][j]+temp.m_pind[i].m_var[j]/(double)temp.m_size;
			}
		}
		m_GDvarArray[m_gens] = median(arr1, temp.m_size);
		m_GDobjArray[m_gens] = median(arr2, temp.m_size);
	}
	delete []arr1;
	delete []arr2;
}

//*****************************************
//m_subpop is the array of subpopulations. 
//m_pop2 is the archive that stores the nondominated solutions found so far.
//*****************************************

//void COEA::runSCCEA(char *filename, int run_count, int mode)
//{
//	//m_subpop is the array of subpopulations. 
//	//m_pop2 is the archive that stores the nondominated solutions found so far.
//	//////////////////////////////////
//	int i,j,k;
//	
//	//Set the size of m_pop2. 
//	m_pop2.initialize(POP2_SIZE);
//	m_memory.initialize(POP2_SIZE);
//	//Set the sizes of the subpopulations. 
//	for (int i=0;i<VAR_NUM;i++) {
//		m_subpop[i].initialize(POP1_SIZE); 
//		m_subpop[i].m_size=m_subpop[i].m_capacity;
//	}
//	//Initialize the subpopulations.
//	for (int i=0;i<VAR_NUM;i++){
//		for(j=0;j<m_subpop[i].m_size;j++) {
//			m_subpop[i].m_pind[j].randomize();m_subpop[i].m_pind[j].m_org_pop = i;
//			for(k=0;k<OBJ_NUM;k++) m_subpop[i].m_pind[j].m_obj[k]=1.0e+10;
//		}
//	}
//	for (int i=0;i<VAR_NUM;i++){
//		m_subpop[i].m_rep_pop = i;
//		for(j=0;j<REP_SIZE;j++) m_subpop[i].m_rep[j]=m_subpop[i].m_pind[j];
//	}
//
//	//Begin main loop
//	SubBinIND mom,dad,sis,bro;
//	SubPOP sp; sp.initialize(POP1_SIZE);sp.m_size=sp.m_capacity;
//	m_re_eval=false;
//	bool randm=false;int gen_change=-999999999;
//	for (m_gens = 0; m_gens < GENS; m_gens++)
//    {
//		//Print out the generation number every 10 generations.
//		if (gen_change+1==m_gens) {
//			printf("Gen %d\n", m_gens);
//		}
//		
//		trackingDynamicsI();
//		if(m_re_eval){
//			gen_change=m_gens;
//			memory2();
//			//m_pop2.m_size=0;
//			/*for (int i=0;i<m_pop2.m_size;i++){
//				m_pop2.m_pind[i].m_genNum=m_gens;
//				m_pop2.m_pind[i].evaluate();
//			}*/
//
//		}
//		//Evaluate every subpopulation.
//		if(gen_change==m_gens){
//			//reinit();
//			randm = true;
//		}
//		if(gen_change==m_gens || (m_gens%mode==0 && gen_change!=m_gens && m_gens!=0)){
//			implicitCompetition();
//			//implicitCompetitionNoRand();
//			randm = true;
//		}
//		else{
//			/*if(gen_change+1==m_gens){
//				double *pA=new double [m_pop2.m_size];int *idx=new int [m_pop2.m_size];
//				for (int i=0;i<m_pop2.m_size;i++){
//					pA[i]=m_pop2.m_pind[i].m_obj[0];
//					idx[i]=i;
//				}
//				EA::selectionSort(pA,idx,m_pop2.m_size);
//				
//				int pos[2];pos[0]=idx[0];pos[1]=idx[m_pop2.m_size-1];
//				for (int i=0;i<2;i++){
//					int gene=0;
//					int tail = POP1_SIZE-i-1;
//					for(k=0;k<VAR_NUM;k++) {
//						for(j=0;j<VAR_BIT;j++)
//							m_subpop[k].m_pind[tail].m_chrom[j]=m_pop2.m_pind[pos[i]].m_chrom[gene++];
//					}
//				}
//				delete[]pA;
//				delete[]idx;
//			}*/
//			for (int i=0;i<VAR_NUM;i++){
//				evaluate(i);
//				m_subpop[i].updateRepresentive();
//				dynamicDistance();
//			}
//			randm = false;
//			extend(CLONE_NUM);
//			
//		}
//		trackingDynamicsII();
//		MS(m_pop2);
//		S(m_pop2);
//		
//		//A cycle where all the subpopulations evolve sequentially once.
//		for (int i=0;i<VAR_NUM;i++){
//			//select, crossover, mutate
//			double shuffle[POP1_SIZE];
//			int index[POP1_SIZE];
//			for(int s=0;s<POP1_SIZE;s++){
//				shuffle[s]=RND.unit();
//				index[s]=s;
//			}
//			EA::selectionSort(shuffle,index,POP1_SIZE,POP1_SIZE);
//			for(j=0;j<sp.m_size;j=j+2){
//				if(randm){
//					mom=m_subpop[i].m_pind[index[j]];
//					dad=m_subpop[i].m_pind[index[j+1]];
//				}
//				else{
//					mom=m_subpop[i].tournament();
//					dad=m_subpop[i].tournament();
//				}
//				sis=mom;bro=dad;
//				mom.crossover(dad,sis,bro);
//				sis.mutate();bro.mutate();
//				sp.m_pind[j]=sis;sp.m_pind[j+1]=bro;
//			}
//			//
//			for(j=0;j<sp.m_size;j++) m_subpop[i].m_pind[j]=sp.m_pind[j];
//			m_subpop[i].m_pind[0]=m_subpop[i].m_rep[0];
//			//
//		}
//		
//		//Output the archive in this generation
//		//char strRst[200];
//		//sprintf(strRst,"f:\\roger\\moeaprg\\tmp\\%d.txt",k);
//		//EA::write(strRst,m_pop2);
//		//printf("\n PS distance is %f, PF distance is %f \n",m_GDvarArray[k],m_GDobjArray[k] );
//    }
//	//end of main loop
//	//
//	//Output m_pop2 of the final generation
//	char out1[200];
//	sprintf(out1,"SCCEA_%s_vardistx%d.txt",PROB_NAME, run_count);
//	out.PrintoNmtxtFile(out1, m_GDvarArray, GENS);
//
//	char out2[200];
//	sprintf(out2,"SCCEA_%s_objdistx%d.txt",PROB_NAME, run_count);
//	out.PrintoNmtxtFile(out2, m_GDobjArray, GENS);
//
//	/*char out3[200];
//	sprintf(out3,"SCCEA_%s_arc_mean_tracex%d.txt",PROB_NAME,run_count);
//	out.PrintoNmtxtFile(out3,m_varTrace, GENS);
//
//	char out4[200];
//	sprintf(out4,"SCCEA_%s_pop_mean_tracex%d.txt",PROB_NAME,run_count);
//	out.PrintoNmtxtFile(out4,m_PopvarTrace, GENS);
//
//	char out5[200];
//	sprintf(out5,"SCCEA_%s_pop_max_tracex%d.txt",PROB_NAME,run_count);
//	out.PrintoNmtxtFile(out5,m_Popvar_max_Trace, GENS);
//
//	char out6[200];
//	sprintf(out6,"SCCEA_%s_pop_min_tracex%d.txt",PROB_NAME,run_count);
//	out.PrintoNmtxtFile(out6,m_Popvar_min_Trace, GENS);
//
//	char out7[200];
//	sprintf(out7,"SCCEA_%s_pop_tracex%d.txt",PROB_NAME,run_count);
//	out.PrintoNmtxtFile(out7, m_Pop_Trace, GENS);*/
//
//	char out8[200];
//	sprintf(out8,"SCCEA_%s_MSx%d.txt",PROB_NAME,run_count);
//	out.PrintoNmtxtFile(out8, m_MS, GENS);
//
//	char out9[200];
//	sprintf(out9,"SCCEA_%s_Sx%d.txt",PROB_NAME,run_count);
//	out.PrintoNmtxtFile(out9, m_S, GENS);
//
//	printf("\nThe pop is now... \n");
//	EA::write("pop.txt",m_pop2);
//	//Output the nondominated members of m_pop1 of the final generation.
//	POP tmpPop;
//	tmpPop.initialize(POP2_SIZE);
//	EA::paretoFront(m_pop2,tmpPop);
//	printf("\nThe nondominated set of pop is now... \n");
//	EA::write(filename,tmpPop);
//}

void COEA::runSCOEA(char *filename, int run_count, int mode)
{
	//m_subpop is the array of subpopulations. 
	//m_pop2 is the archive that stores the nondominated solutions found so far.
	//////////////////////////////////
	int i,j,k;
	
	//Set the size of m_pop2. 
	m_pop2.initialize(POP2_SIZE);
	//Set the sizes of the subpopulations. 
	for (int i=0;i<VAR_NUM;i++) {
		m_subpop[i].initialize(POP1_SIZE); 
		m_subpop[i].m_size=m_subpop[i].m_capacity;
	}
	//Initialize the subpopulations.
	for (int i=0;i<VAR_NUM;i++){
		for(j=0;j<m_subpop[i].m_size;j++) {
			m_subpop[i].m_pind[j].randomize();m_subpop[i].m_pind[j].m_org_pop = i;
			for(k=0;k<OBJ_NUM;k++) m_subpop[i].m_pind[j].m_obj[k]=1.0e+10;
		}
	}
	for (int i=0;i<VAR_NUM;i++){
		m_subpop[i].m_rep_pop = i;
		for(j=0;j<REP_SIZE;j++) m_subpop[i].m_rep[j]=m_subpop[i].m_pind[j];
	}

	//Begin main loop
	SubBinIND mom,dad,sis,bro;
	SubPOP sp; sp.initialize(POP1_SIZE);sp.m_size=sp.m_capacity;
	m_re_eval=false;
	bool randm=false;int gen_change=-999999999;
	for (m_gens = 0; m_gens < GENS; m_gens++)
    {
		//Print out the generation number every 10 generations.
		if (gen_change+1==m_gens) {
			printf("Gen %d\n", m_gens);
		}
		trackingDynamicsI();
		if(m_re_eval){
			gen_change=m_gens;
			//m_pop2.m_size=0;
		}
		//Evaluate every subpopulation.
		if(gen_change==m_gens){
			//reinit();
			gen_change=m_gens;
			memory1();
			//m_pop2.m_size=0;
			/*for (int i=0;i<m_pop2.m_size;i++){
				m_pop2.m_pind[i].m_genNum=m_gens;
				m_pop2.m_pind[i].evaluate();
			}*/
			//randm = true;
		}
		
		randm = false;
		for (int i=0;i<VAR_NUM;i++){
			evaluate(i);
			m_subpop[i].updateRepresentive();
			dynamicDistance();
		}
		extend(CLONE_NUM);
		trackingDynamicsII();
		MS(m_pop2);
		S(m_pop2);
		
		
		//A cycle where all the subpopulations evolve sequentially once.
		for (int i=0;i<VAR_NUM;i++){
			//select, crossover, mutate
			double shuffle[POP1_SIZE];
			int index[POP1_SIZE];
			for(int s=0;s<POP1_SIZE;s++){
				shuffle[s]=RND.unit();
				index[s]=s;
			}
			EA::selectionSort(shuffle,index,POP1_SIZE,POP1_SIZE);
			for(j=0;j<sp.m_size;j=j+2){
				if(randm){
					mom=m_subpop[i].m_pind[index[j]];
					dad=m_subpop[i].m_pind[index[j+1]];
				}
				else{
					mom=m_subpop[i].tournament();
					dad=m_subpop[i].tournament();
				}
				sis=mom;bro=dad;
				mom.crossover(dad,sis,bro);
				sis.mutate();bro.mutate();
				sp.m_pind[j]=sis;sp.m_pind[j+1]=bro;
			}
			//
			for(j=0;j<sp.m_size;j++) m_subpop[i].m_pind[j]=sp.m_pind[j];
			m_subpop[i].m_pind[0]=m_subpop[i].m_rep[0];
			//
		}
		
		//Output the archive in this generation
		//char strRst[200];
		//sprintf(strRst,"f:\\roger\\moeaprg\\tmp\\%d.txt",k);
		//EA::write(strRst,m_pop2);
		//printf("\n PS distance is %f, PF distance is %f \n",m_GDvarArray[k],m_GDobjArray[k] );
    }
	//end of main loop
	//
	//Output m_pop2 of the final generation
	//char out1[200];
	//sprintf(out1,"SCOEA_%s_vardistx%d.txt",PROB_NAME, run_count);
	//out.PrintoNmtxtFile(out1, m_GDvarArray, GENS);

	//char out2[200];
	//sprintf(out2,"SCOEA_%s_objdistx%d.txt",PROB_NAME, run_count);
	//out.PrintoNmtxtFile(out2, m_GDobjArray, GENS);

	/*char out3[200];
	sprintf(out3,"SCOEA_%s_arc_mean_tracex%d.txt",PROB_NAME,run_count);
	out.PrintoNmtxtFile(out3,m_varTrace, GENS);

	char out4[200];
	sprintf(out4,"SCOEA_%s_pop_mean_tracex%d.txt",PROB_NAME,run_count);
	out.PrintoNmtxtFile(out4,m_PopvarTrace, GENS);

	char out5[200];
	sprintf(out5,"SCOEA_%s_pop_max_tracex%d.txt",PROB_NAME,run_count);
	out.PrintoNmtxtFile(out5,m_Popvar_max_Trace, GENS);

	char out6[200];
	sprintf(out6,"SCOEA_%s_pop_min_tracex%d.txt",PROB_NAME,run_count);
	out.PrintoNmtxtFile(out6,m_Popvar_min_Trace, GENS);

	char out7[200];
	sprintf(out7,"SCOEA_%s_pop_tracex%d.txt",PROB_NAME,run_count);
	out.PrintoNmtxtFile(out7, m_Pop_Trace, GENS);*/

	/*char out8[200];
	sprintf(out8,"SCOEA_%s_MSx%d.txt",PROB_NAME,run_count);
	out.PrintoNmtxtFile(out8, m_MS, GENS);

	char out9[200];
	sprintf(out9,"SCOEA_%s_Sx%d.txt",PROB_NAME,run_count);
	out.PrintoNmtxtFile(out9, m_S, GENS);*/

	printf("\nThe pop is now... \n");
	EA::write("pop.txt",m_pop2);
	//Output the nondominated members of m_pop1 of the final generation.
	POP tmpPop;
	tmpPop.initialize(POP2_SIZE);
	EA::paretoFront(m_pop2,tmpPop);
	printf("\nThe nondominated set of pop is now... \n");
	EA::write(filename,tmpPop);
}

//void COEA::runCCEA(char *filename, int run_count, int mode)
//{
//	//m_subpop is the array of subpopulations. 
//	//m_pop2 is the archive that stores the nondominated solutions found so far.
//	//////////////////////////////////
//	int i,j,k;
//	
//	//Set the size of m_pop2. 
//	m_pop2.initialize(POP2_SIZE);
//	//Set the sizes of the subpopulations. 
//	for (int i=0;i<VAR_NUM;i++) {
//		m_subpop[i].initialize(POP1_SIZE); 
//		m_subpop[i].m_size=m_subpop[i].m_capacity;
//	}
//	//Initialize the subpopulations.
//	for (int i=0;i<VAR_NUM;i++){
//		for(j=0;j<m_subpop[i].m_size;j++) {
//			m_subpop[i].m_pind[j].randomize();m_subpop[i].m_pind[j].m_org_pop = i;
//			for(k=0;k<OBJ_NUM;k++) m_subpop[i].m_pind[j].m_obj[k]=1.0e+10;
//		}
//	}
//	for (int i=0;i<VAR_NUM;i++){
//		m_subpop[i].m_rep_pop = i;
//		for(j=0;j<REP_SIZE;j++) m_subpop[i].m_rep[j]=m_subpop[i].m_pind[j];
//	}
//
//	//Begin main loop
//	SubBinIND mom,dad,sis,bro;
//	SubPOP sp; sp.initialize(POP1_SIZE);sp.m_size=sp.m_capacity;
//	bool randm=false;
//	for (m_gens = 0; m_gens < GENS; m_gens++)
//    {
//		//Print out the generation number every 10 generations.
//		if (m_gens%10==0) {
//			printf("Gen %d\n", m_gens);
//		}
//		//trackingDynamicsI();
//		//Evaluate every subpopulation.
//		/*if(m_gens%mode==0 && m_gens!=1){
//			//implicitCompetition();
//			implicitCompetitionNoRand();
//			randm=true;
//		}
//		else{*/
//			for (int i=0;i<VAR_NUM;i++){
//				evaluate(i);
//				m_subpop[i].updateRepresentive();
//				dynamicDistance();
//			}
//			extend(CLONE_NUM);
//			randm=false;
//		//}
//		for (int i=0;i<VAR_NUM;i++){
//			m_PopvarTrace[m_gens][i]=m_subpop[i].m_rep[0].m_var;
//		}
//		//trackingDynamicsII();
//		if (m_gens%5==0) {
//			//char out1[200];
//			printf("Gen %d\n", m_gens);
//			//sprintf(out1,"%s_archivex%d_%d.txt",PROB_NAME,run_count,m_gens);
//			//EA::write(out1,m_pop2);
//		}
//		
//		//A cycle where all the subpopulations evolve sequentially once.
//		for (int i=0;i<VAR_NUM;i++){
//			//select, crossover, mutate
//			double shuffle[POP1_SIZE];
//			int index[POP1_SIZE];
//			for(int s=0;s<POP1_SIZE;s++){
//				shuffle[s]=RND.unit();	index[s]=s;
//			}
//			EA::selectionSort(shuffle,index,POP1_SIZE,POP1_SIZE);
//			for(j=0;j<sp.m_size;j=j+2){
//				if(randm){
//					mom=m_subpop[i].m_pind[index[j]];
//					dad=m_subpop[i].m_pind[index[j+1]];
//				}
//				else{
//					mom=m_subpop[i].tournament();
//					dad=m_subpop[i].tournament();
//				}
//				sis=mom;bro=dad;
//				mom.crossover(dad,sis,bro);
//				sis.mutate();bro.mutate();
//				sp.m_pind[j]=sis;sp.m_pind[j+1]=bro;
//			}
//			//
//			for(j=0;j<sp.m_size;j++) m_subpop[i].m_pind[j]=sp.m_pind[j];
//			m_subpop[i].m_pind[0]=m_subpop[i].m_rep[0];
//			//
//		}
//	}
//	//end of main loop
//	//
//	//Output m_pop2 of the final generation
//	/*char out1[200];
//	sprintf(out1,"%s_%s_var_distx%d.txt",PROB_NAME,COMPETE,run_count);
//	out.PrintoNmtxtFile(out1, m_GDvarArray, GENS);*/
//
//	/*char out2[200];
//	sprintf(out2,"COEA %s_obj_distx%d.txt",PROB_NAME,run_count);
//	out.PrintoNmtxtFile(out2, m_GDobjArray, GENS);*/
//
//	/*char out3[200];
//	sprintf(out3,"%s_%s_arc_mean_tracex%d.txt",PROB_NAME,COMPETE,run_count);
//	out.PrintoNmtxtFile(out3,m_varTrace, GENS);
//
//	char out4[200];
//	sprintf(out4,"%s_%s_pop_mean_tracex%d.txt",PROB_NAME,COMPETE,run_count);
//	out.PrintoNmtxtFile(out4,m_PopvarTrace, GENS);
//
//	char out5[200];
//	sprintf(out5,"%s_%s_pop_max_tracex%d.txt",PROB_NAME,COMPETE,run_count);
//	out.PrintoNmtxtFile(out5,m_Popvar_max_Trace, GENS);
//
//	char out6[200];
//	sprintf(out6,"%s_%s_pop_min_tracex%d.txt",PROB_NAME,COMPETE,run_count);
//	out.PrintoNmtxtFile(out6,m_Popvar_min_Trace, GENS);*/
//
//	char out7[200];
//	sprintf(out7,"%s_pop_tracex%d.txt",PROB_NAME,run_count);
//	out.PrintoNmtxtFile(out7, m_Pop_Trace, GENS);
//
//	char out8[200];
//	sprintf(out8,"%s_var_tracex%d.txt",PROB_NAME,run_count);
//	out.PrintoNmtxtFile(out8,m_PopvarTrace, GENS);
//
//
//	printf("\nThe pop is now... \n");
//	EA::write("pop.txt",m_pop2);
//	//Output the nondominated members of m_pop1 of the final generation.
//	POP tmpPop;
//	tmpPop.initialize(POP2_SIZE);
//	EA::paretoFront(m_pop2,tmpPop);
//	printf("\nThe nondominated set of pop is now... \n");
//	EA::write(filename,tmpPop);
//}
//
//void COEA::runCOEA(char *filename, char *filename2, int run_count, int mode)
//{
//	//m_subpop is the array of subpopulations. 
//	//m_pop2 is the archive that stores the nondominated solutions found so far.
//	//////////////////////////////////
//	int i,j,k;
//	
//	//Set the size of m_pop2. 
//	m_pop2.initialize(POP2_SIZE);
//	//Set the sizes of the subpopulations. 
//	for (int i=0;i<VAR_NUM;i++) {
//		m_subpop[i].initialize(POP1_SIZE); 
//		m_subpop[i].m_size=m_subpop[i].m_capacity;
//	}
//	//Initialize the subpopulations.
//	for (int i=0;i<VAR_NUM;i++){
//		for(j=0;j<m_subpop[i].m_size;j++) {
//			m_subpop[i].m_pind[j].randomize();m_subpop[i].m_pind[j].m_org_pop = i;
//			for(k=0;k<OBJ_NUM;k++) m_subpop[i].m_pind[j].m_obj[k]=1.0e+10;
//		}
//	}
//	for (int i=0;i<VAR_NUM;i++){
//		m_subpop[i].m_rep_pop = i;
//		for(j=0;j<REP_SIZE;j++) m_subpop[i].m_rep[j]=m_subpop[i].m_pind[j];
//	}
//
//	//Begin main loop
//	SubBinIND mom,dad,sis,bro;
//	SubPOP sp; sp.initialize(POP1_SIZE);sp.m_size=sp.m_capacity;
//
//	for (m_gens = 0; m_gens < GENS; m_gens++)
//    {
//		
//		//Evaluate every subpopulation.
//                //%%%%%%%%%%%%%%%%%%
//		//trackingDynamicsII();  //choose here I or II
//                //%%%%%%%%%%%%%%%%%% 
//		for (int i=0;i<VAR_NUM;i++){
//			evaluate(i);
//			m_subpop[i].updateRepresentive();
//			dynamicDistance();
//		}
//		extend(CLONE_NUM);
//		
//		//A cycle where all the subpopulations evolve sequentially once.
//		for (int i=0;i<VAR_NUM;i++){
//			for(j=0;j<sp.m_size;j=j+2){
//				mom=m_subpop[i].tournament();
//				dad=m_subpop[i].tournament();
//				sis=mom;bro=dad;
//				mom.crossover(dad,sis,bro);
//				sis.mutate();bro.mutate();
//				sp.m_pind[j]=sis;sp.m_pind[j+1]=bro;
//			}
//				//
//			for(j=0;j<sp.m_size;j++) m_subpop[i].m_pind[j]=sp.m_pind[j];
//			m_subpop[i].m_pind[0]=m_subpop[i].m_rep[0];
//		}
//		//Output the archive in this generation
//		//char strRst[200];
//		//sprintf(strRst,"f:\\roger\\moeaprg\\tmp\\%d.txt",k);
//		//EA::write(strRst,m_pop2);
//		//printf("\n PS distance is %f, PF distance is %f \n",m_GDvarArray[k],m_GDobjArray[k] );
//    }
//	//end of main loop
//	//
//	//Output m_pop2 of the final generation
//	/*char out1[200];
//	sprintf(out1,"%s_%s_var_distx%d.txt",PROB_NAME,COMPETE,run_count);
//	out.PrintoNmtxtFile(out1, m_GDvarArray, GENS);
//
//	char out2[200];
//	sprintf(out2,"%s_%s_obj_distx%d.txt",PROB_NAME,COMPETE,run_count);
//	out.PrintoNmtxtFile(out2, m_GDobjArray, GENS);
//
//	char out3[200];
//	sprintf(out3,"%s_%s_arc_mean_tracex%d.txt",PROB_NAME,COMPETE,run_count);
//	out.PrintoNmtxtFile(out3,m_varTrace, GENS);
//
//	char out4[200];
//	sprintf(out4,"%s_%s_pop_mean_tracex%d.txt",PROB_NAME,COMPETE,run_count);
//	out.PrintoNmtxtFile(out4,m_PopvarTrace, GENS);
//
//	char out5[200];
//	sprintf(out5,"%s_%s_pop_max_tracex%d.txt",PROB_NAME,COMPETE,run_count);
//	out.PrintoNmtxtFile(out5,m_Popvar_max_Trace, GENS);
//
//	char out6[200];
//	sprintf(out6,"%s_%s_pop_min_tracex%d.txt",PROB_NAME,COMPETE,run_count);
//	out.PrintoNmtxtFile(out6,m_Popvar_min_Trace, GENS);
//
//	char out7[200];
//	sprintf(out7,"%s_%s_pop_tracex%d.txt",PROB_NAME,COMPETE,run_count);
//	out.PrintoNmtxtFile(out7, m_Pop_Trace, GENS);*/
//
//	printf("\nThe pop is now... \n");
//	EA::write("pop.txt",m_pop2);
//	//Output the nondominated members of m_pop1 of the final generation.
//	POP tmpPop;
//	tmpPop.initialize(POP2_SIZE);
//	EA::paretoFront(m_pop2,tmpPop);
//	printf("\nThe nondominated set of pop is now... \n");
//	EA::write(filename,tmpPop);
//	
//}

void COEA::writeVariables(char*filename,const POP& pop)
{
	ofstream fout;
	fout.open(filename,ios::out);
	for (int i = 0; i < pop.m_size; i++) 
	{ 
		//m_pop2.m_pind[i].printObjectives();
		pop.m_pind[i].printVariables(fout);
	}
	fout.close();
}