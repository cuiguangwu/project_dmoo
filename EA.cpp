// Evolutionary algorithm class implementation

#include "EA.h"

///////////////////////////
//Initialize the static member of class IND
IND::Evaluator IND::evalFp = NULL;

////////////////////////////

void Grid::initialize(POP *pp,int objNum,int depth)
{
	pPop=pp;
	nObjNum=objNum;
	nDepth=depth;
	//
	gl_offset= new double[objNum];
	gl_largest= new double[objNum];
	gl_range= new double[objNum];
	for(int i=0;i<objNum;i++) { gl_offset[i]=0;gl_largest[i]=0;gl_range[i]=1;}
	//
	double tmp;
	tmp=pow(2.0,objNum*depth)+1;
	grid_pop= new int[(int)tmp];
	for(int i=0;i<tmp;i++) grid_pop[i]=0;
}

Grid::~Grid()
{
	if(pPop!=NULL){
		delete[] gl_offset; delete[] gl_largest; delete[] gl_range;
		delete[] grid_pop;
	}
}

//Pre:	1 the members of the pop changed
//		2 the grid_loc of new members should be set to -1
//Post: 1 the new offset,largest, range and grid_pop of the grid if necessory. 
//		2 the new grid_pop of the pop if the grid changes
//		3 the grid_loc of new members
void Grid::update_grid()
{
	// recalculate ranges for grid in the light of a new solution s
	static int change = 0;
	double *offset=new double[nObjNum];
	double *largest=new double[nObjNum];
	//
	for (int a = 0; a < nObjNum; a++) { offset[a] = 1.0e+10;largest[a] = -1.0e+10; }
	for (int b = 0; b < nObjNum; b++)
	{
		for (int a = 0; a < pPop->m_size;a++)
		{
			if (pPop->m_pind[a].m_obj[b] < offset[b])
				offset[b] = pPop->m_pind[a].m_obj[b];
			if (pPop->m_pind[a].m_obj[b] > largest[b])
				largest[b] = pPop->m_pind[a].m_obj[b];		   
		}
	}
	//Find the changes of ranges of the grid. 
	//If the changes are some large, recalculate the grid.
	double range,diff1,diff2;
	bool modify=false;
	for (int a = 0; a < nObjNum; a++)
	{
		range=largest[a]-offset[a];
		diff1=gl_offset[a] - offset[a];
		diff2=gl_largest[a]-largest[a];
		if( diff1<-0.1*range || diff1>0 ) modify=true;
		if( diff2<0 || diff2>0.1*range ) modify=true;
	}
	if (modify)	
	{                                  
		//if the summed squared error (difference) between old and new
		//minima and maxima in each of the objectives 
		//is bigger than 10 percent of the square of the size of the space
		// then renormalise the space and recalculte grid locations
		change++;                        

		for (int a = 0; a < nObjNum; a++)
		{
		  gl_largest[a] = largest[a]+0.1*(largest[a]-offset[a]);
		  gl_offset[a] = offset[a]-0.1*(largest[a]-offset[a]);
		  gl_range[a] = gl_largest[a] - gl_offset[a];
		}
		//recalculate the grid_loc and grid_pop
		int square;
		for (int a=(int)pow(2.0,(nObjNum*nDepth));a>=0;a--) grid_pop[a] = 0;
		for (int a=0; a < pPop->m_size; a++)
		{
			square = find_loc(pPop->m_pind[a]);
			pPop->m_pind[a].m_gridLoc = square;
			grid_pop[square]++;
		}
	}
	else
	{
		//recalculate the grid_pop
		int square;
		for (int a=(int)pow(2.0,(nObjNum*nDepth));a>=0;a--) grid_pop[a] = 0;
		for (int a=0; a < pPop->m_size; a++)
		{
			square=pPop->m_pind[a].m_gridLoc;
			//calculate the grid_loc of new elements
			if(square<0) {
				square = find_loc(pPop->m_pind[a]);
				pPop->m_pind[a].m_gridLoc = square;
			}
			//
			grid_pop[square]++;
		}
	}
	//
	//members outside the grid have the lowest grid_pop
	grid_pop[(int)pow(2.0,(nObjNum*nDepth))] = -5;
	//
	delete[] offset; delete[] largest;
	return;
}

//update the grid with ind and the grid location of ind
//Pre: set -1 to grid_loc of ind
//Post: new pop, new grid, grid_loc of ind
void Grid::update_pop(BinIND &ind)
{
	// given a solution s, add it to the archive if
	// a) the archive is empty 
	// b) the archive is not full and s is not dominated or equal to anything currently in the archive
	// c) s dominates anything in the archive                        
	// d) the archive is full but s is nondominated and is in a no more crowded square than at least one solution
	// in addition, maintain the archive such that all solutions are nondominated. 
 
	ind.m_gridLoc=-1;
	//case a)
	if (pPop->m_size==0) { 
		pPop->add(ind);
		update_grid();
		ind.m_gridLoc=find_loc(ind);
		return; 
	}
	//
	int *tag=new int[pPop->m_size];
	for (int i=0;i<pPop->m_size;i++) tag[i]=0;
	int join = 0;	//if the ind comes into the pop
	int set = 0;	//if delete members of the pop
	IND::nRel result = IND::INCOMPARABLE;//if the ind dominates the pop
	for(int i=0;i<pPop->m_size;i++)
	{
		result = ind.compare(pPop->m_pind[i]);
		if (result==IND::EQUAL || result==IND::BIGGER) break;
		if ( (result==IND::SMALLER)&&(join==0) ) { pPop->m_pind[i]=ind; join=1;}
		else if (result==IND::SMALLER) { tag[i]=1;set=1;}	    	    
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
			int loc;
			int most = -10000;
			int repl;
			for (int i = 0; i < pPop->m_size; i++)
			{
				loc=pPop->m_pind[i].m_gridLoc;
				if (grid_pop[loc] > most)
				{
					most = grid_pop[loc];
					repl = i;
				}
			}
			pPop->m_pind[repl] = ind;
		}
		//case b)
		else 
			pPop->add(ind);
	}
	//update grid
	if( result==IND::EQUAL||result==IND::BIGGER ) 
		ind.m_gridLoc=find_loc(ind);
	else {
		update_grid();
		ind.m_gridLoc=find_loc(ind);
	}
	//
	delete[] tag;
	return;
}

//find the grid location of a solution given a vector of its objective values
//Post: grid location of ind
int Grid::find_loc(IND &ind)
{
	// if the solution is out of range on any objective, return 1 more than the maximum possible grid location number
	for (int i = 0; i < nObjNum; i++)
	{
		if ( (ind.m_obj[i]< gl_offset[i])||(ind.m_obj[i]>gl_offset[i]+gl_range[i]) )
		return( (int)pow(2.0,(nObjNum*nDepth)) );
	}
	//set inc and grid width
	int *inc=new int[nObjNum];
	double *width= new double[nObjNum];
	double *offset= new double[nObjNum];
	int n=1;
	for (int i = 0; i < nObjNum; i++)
	{
	  inc[i] = n; n*=2;
	  width[i]=gl_range[i]; offset[i]=gl_offset[i];
	}
	//find grid_loc using binary search 	
	int loc = 0;
	for (int d = 1; d <= nDepth; d++)
	{
		for (int i = 0; i < nObjNum; i++)
		{
			if(ind.m_obj[i] < width[i]/2+offset[i])
			loc += inc[i];
			else
			offset[i] += width[i]/2;
		}
		for (int i = 0; i < nObjNum; i++)
		{
			inc[i] *= (nObjNum *2);
			width[i] /= 2;
		}
	}
	//
	delete[] inc;delete[] width;delete[] offset;
	//
	return loc;
}
