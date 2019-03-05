/* Data initializtion routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to initialize a population randomly */
void initialize_pop (population *pop)
{
	int i;
	for (i=0; i<popsize; i++)
	{
		initialize_ind (&(pop->ind[i]));
	}
	return;
}

/* Function to initialize an individual randomly */
void initialize_ind (individual *ind)
{
	int j, k;
	if (nreal!=0)
	{
		for (j=0; j<nreal; j++)
		{
			ind->xreal[j] = rndreal (min_realvar[j], max_realvar[j]);
		}
	}
	ind->boundary_ind = 0;
	if (nbin!=0)
	{
		for (j=0; j<nbin; j++)
		{
			for (k=0; k<nbits[j]; k++)
			{
				if (randomperc() <= 0.5)
				{
					ind->gene[j][k] = 0;
				}
				else
				{
					ind->gene[j][k] = 1;
				}
			}
		}
	}
	return;
}
void initialize_area (unit_area *area)
{
	area->area_exist = 0;
	area->area_indcount = 0;
	area->area_repreasent_ind = -1;
	area->area_nondominated_number = 0;
	area->area_ind->index = -1;
	area->area_ind->parent = NULL;
	area->area_ind->child = NULL;
	area->area_non_domianted->index = -1;
	area->area_non_domianted->parent = NULL;
	area->area_non_domianted->child = NULL;
}
