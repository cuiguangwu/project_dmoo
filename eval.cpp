/* Routine for evaluating population members  */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
#include "emo/IndividualMO.h"
#include "emo/PopulationMO.h"

# include "global.h"
# include "rand.h"
#include "./AlgD.h"
__declspec(dllimport) void test_problem_wfg(double *xreal,double *obj,int nreal,int nobj,int k, char *type);

/* Routine to evaluate objective function values and constraints for a population */
//void evaluate_pop (population *pop)
//{
//	int i;
//	for (i=0; i<popsize; i++)
//	{
//		evaluate_ind (&(pop->ind[i]));
//	}
//	return;
//}

/* Routine to evaluate objective function values and constraints for an individual */
void evaluate_ind (individual *ind, az::mea::CIndividualMO& eind)
{
	int j;
	//test_problem (ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
	//	test_problem_wfg(ind->xreal, ind->obj, nreal, nobj, 2, "wfg9");
	for (j=0;j<nreal;j++)
	{
		eind[j]=ind->xreal[j];
	}
	eind.Evaluate();
	for (j=0;j<nobj;j++)
	{
		ind->obj[j]=eind.F(j);
	}
	eval_time++;
	ind->boundary_ind = 0;
	if (ncon==0)
	{
		ind->constr_violation = 0.0;
	}
	else
	{
		ind->constr_violation = 0.0;
		for (j=0; j<ncon; j++)
		{
			if (ind->constr[j]<0.0)
			{
				ind->constr_violation += ind->constr[j];
			}
		}
	}
	return;
}
