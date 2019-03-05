/* Crossover routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to cross two individuals */
void crossover (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
	if (nreal!=0)
	{
		realcross (parent1, parent2, child1, child2);
	}
	if (nbin!=0)
	{
		bincross (parent1, parent2, child1, child2);
	}
	return;
}

/* Routine for real variable SBX crossover */
void realcross (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
	int i;
	double rand;
	double y1, y2, yl, yu, a1, a2;
	double alpha, beta, betaq;
	if (randomperc() <= pcross_real)
	{
		nrealcross++;
		for (i=0; i<nreal; i++)
		{
			if (randomperc() <= 0.5 )
			{
				if (parent1->xreal[i] < parent2->xreal[i])
				{
					y1 = parent1->xreal[i];
					y2 = parent2->xreal[i];
				}
				else
				{
					y1 = parent2->xreal[i];
					y2 = parent1->xreal[i];
				}
				if (fabs(parent1->xreal[i]-parent2->xreal[i]) > EPS)
				{
					//printf("i=%d\n",i);
					//printf("min=%f\n",min_realvar[i]);
					yl = min_realvar[i];
					yu = max_realvar[i];
					rand = randomperc();
					beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
					alpha = 2.0 - pow(beta,-(eta_c+1.0));
					if (rand <= (1.0/alpha))
					{
						betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
					}
					else
					{
						betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
					}
					a1 = 0.5*((y1+y2)-betaq*(y2-y1));
					beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
					alpha = 2.0 - pow(beta,-(eta_c+1.0));
					if (rand <= (1.0/alpha))
					{
						betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
					}
					else
					{
						betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
					}
					a2 = 0.5*((y1+y2)+betaq*(y2-y1));

					if (randomperc()<0.5)
					{
						child1->xreal[i] = a1;
						child2->xreal[i] = a2;
					}
					else
					{
						child1->xreal[i] = a2;
						child2->xreal[i] = a1;
					}

					if (child1->xreal[i]<yl)
					{
						child1->xreal[i]=yl;
					}
					if (child1->xreal[i]>yu)
					{
						child1->xreal[i]=yu;
					}
					if (child2->xreal[i]<yl)
					{
						child2->xreal[i]=yl;
					}
					if (child2->xreal[i]>yu)
					{
						child2->xreal[i]=yu;
					}
				}
				else
				{
					child1->xreal[i] = parent1->xreal[i];
					child2->xreal[i] = parent2->xreal[i];
				}
			}
			else
			{
				child1->xreal[i] = parent1->xreal[i];
				child2->xreal[i] = parent2->xreal[i];
			}
		}
	}
	else
	{
		for (i=0; i<nreal; i++)
		{
			child1->xreal[i] = parent1->xreal[i];
			child2->xreal[i] = parent2->xreal[i];
		}
	}
	return;
}

/* Routine for two point binary crossover */
void bincross (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
	int i, j;
	double rand;
	int temp, site1, site2;
	for (i=0; i<nbin; i++)
	{
		rand = randomperc();
		if (rand <= pcross_bin)
		{
			nbincross++;
			site1 = rnd(0,nbits[i]-1);
			site2 = rnd(0,nbits[i]-1);
			if (site1 > site2)
			{
				temp = site1;
				site1 = site2;
				site2 = temp;
			}
			for (j=0; j<site1; j++)
			{
				child1->gene[i][j] = parent1->gene[i][j];
				child2->gene[i][j] = parent2->gene[i][j];
			}
			for (j=site1; j<site2; j++)
			{
				child1->gene[i][j] = parent2->gene[i][j];
				child2->gene[i][j] = parent1->gene[i][j];
			}
			for (j=site2; j<nbits[i]; j++)
			{
				child1->gene[i][j] = parent1->gene[i][j];
				child2->gene[i][j] = parent2->gene[i][j];
			}
		}
		else
		{
			for (j=0; j<nbits[i]; j++)
			{
				child1->gene[i][j] = parent1->gene[i][j];
				child2->gene[i][j] = parent2->gene[i][j];
			}
		}
	}
	return;
}

void crossover_arithmetic (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
	if (nreal!=0)
	{
		realcross_arithmetic (parent1, parent2, child1, child2);
	}
	return;
}



void realcross_arithmetic (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
	int i;
	double rand;
	double y1, y2, yl, yu;
	double alpha, beta, betaq;

	nrealcross++;


	for (i=0; i<nreal; i++)
	{

		if (randomperc() <= 0.5)		
		{
			rand = randomperc();
			child1->xreal[i] = rand*(parent2->xreal[i]-parent1->xreal[i])+parent1->xreal[i];
			//			rand = randomperc();
			//			child2->xreal[i] = rand*(parent1->xreal[i]-parent2->xreal[i])+parent2->xreal[i];
		}
		else
		{
			child1->xreal[i] = parent1->xreal[i];
			//			child2->xreal[i] = parent2->xreal[i];
		}
		if (randomperc() <= 0.5)		
		{
			//			rand = randomperc();
			//			child1->xreal[i] = rand*(parent2->xreal[i]-parent1->xreal[i])+parent1->xreal[i];
			rand = randomperc();
			child2->xreal[i] = rand*(parent1->xreal[i]-parent2->xreal[i])+parent2->xreal[i];
		}
		else
		{
			//			child1->xreal[i] = parent1->xreal[i];
			child2->xreal[i] = parent2->xreal[i];
		}
	}


	return;
}

void realcross_heuristic(individual *parent1, individual *parent2, individual *child1, individual *child2, double para)
{
	int i;
	double rand;
	double y1, y2, yl, yu;
	double alpha, beta, betaq;

	nrealcross++;


	for (i=0; i<nreal; i++)
	{
		if (randomperc() <= 0.5)
		{		
			rand = rndreal(0, para*2);
			child2->xreal[i] = rand*(parent1->xreal[i]-parent2->xreal[i])+parent1->xreal[i];
			if (child2->xreal[i]<min_realvar[i])
				//				child2->xreal[i]=2*min_realvar[i]-child2->xreal[i];
				child2->xreal[i]=min_realvar[i];
			if (child2->xreal[i]>max_realvar[i])
				//				child2->xreal[i]=2*max_realvar[i]-child2->xreal[i];
				child2->xreal[i]=max_realvar[i];
		}
		else
		{
			child2->xreal[i] = parent1->xreal[i];
		}


	}
	for (i=0; i<nreal; i++)
	{
		if (randomperc() <= 0.5)
		{		
			rand = rndreal(0, para*2);
			child1->xreal[i] = rand*(parent1->xreal[i]-parent2->xreal[i])+parent1->xreal[i];
			if (child1->xreal[i]<min_realvar[i])
				//			    child1->xreal[i]=2*min_realvar[i]-child1->xreal[i];
				child1->xreal[i]=min_realvar[i];
			if (child1->xreal[i]>max_realvar[i])
				//				child1->xreal[i]=2*max_realvar[i]-child1->xreal[i];
				child1->xreal[i]=max_realvar[i];
		}
		else
		{
			child1->xreal[i] = parent1->xreal[i];
		}


	}

	return;
}

