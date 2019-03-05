/* Memory allocation and deallocation routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to allocate memory to a population */
/* 为种群开辟内存空间 */
void allocate_memory_pop (population *pop, int size)
{
    int i;
	//printf("size=%d\n",size);
    pop->ind = (individual *)malloc(size*sizeof(individual));
    for (i=0; i<size; i++)
    {
        allocate_memory_ind (&(pop->ind[i]));
    }
    return;
}

/* Function to allocate memory to an individual */
/* 为个体开辟内存空间 */
void allocate_memory_ind (individual *ind)
{
    int j;
//	ind->coordinate = (int *)malloc(nobj*sizeof(int));
    if (nreal != 0)
    {
        ind->xreal = (double *)malloc(nreal*sizeof(double));
    }
    if (nbin != 0)
    {
        ind->xbin = (double *)malloc(nbin*sizeof(double));
        ind->gene = (int **)malloc(nbin*sizeof(int));
        for (j=0; j<nbin; j++)
        {
            ind->gene[j] = (int *)malloc(nbits[j]*sizeof(int));
        }
    }
    ind->obj = (double *)malloc(nobj*sizeof(double));
    if (ncon != 0)
    {
        ind->constr = (double *)malloc(ncon*sizeof(double));
    }
    return;
}

/* 为单元域开辟内存空间 */
void allocate_memory_area (unit_area *area)
{
    area->area_ind = (list *)malloc(sizeof(list));
	area->area_non_domianted = (list *)malloc(sizeof(list));
    return;
}
/*
void allocate_memory_envind (env_individual *ind)
{
    int j;
    if (nreal != 0)
    {
        ind->xreal = (double *)malloc(nreal*sizeof(double));
    }
    if (nbin != 0)
    {
        ind->xbin = (double *)malloc(nbin*sizeof(double));
        ind->gene = (int **)malloc(nbin*sizeof(int));
        for (j=0; j<nbin; j++)
        {
            ind->gene[j] = (int *)malloc(nbits[j]*sizeof(int));
        }
    }
    ind->obj = (double *)malloc(nobj*sizeof(double));
    if (ncon != 0)
    {
        ind->constr = (double *)malloc(ncon*sizeof(double));
    }
    return;
}
*/
/* Function to deallocate memory to a population */
/* 释放空间 */
void deallocate_memory_pop (population *pop, int size)
{
    int i;
    for (i=0; i<size; i++)
    {
        deallocate_memory_ind (&(pop->ind[i]));
    }
    free (pop->ind);
    return;
}
/* 释放空间 */
void deallocate_memory_area (unit_area *area)
{
    list *temp;
	while (area->area_ind!=NULL)
    {
        temp = area->area_ind;
        area->area_ind = area->area_ind->child;
        free (temp);
    }
	while (area->area_non_domianted!=NULL)
    {
        temp = area->area_non_domianted;
        area->area_non_domianted = area->area_non_domianted->child;
        free (temp);
    }
    
    return;
}

/* Function to deallocate memory to an individual */
/* 释放空间 */
void deallocate_memory_ind (individual *ind)
{
    int j;
//	free(ind->coordinate);
    if (nreal != 0)
    {
        free(ind->xreal);
    }
    if (nbin != 0)
    {
        for (j=0; j<nbin; j++)
        {
            free(ind->gene[j]);
        }
        free(ind->xbin);
        free(ind->gene);
    }
    free(ind->obj);
    if (ncon != 0)
    {
        free(ind->constr);
    }
    return;
}
/*
void deallocate_memory_envind (env_individual *ind)
{
    int j;
    if (nreal != 0)
    {
        free(ind->xreal);
    }
    if (nbin != 0)
    {
        for (j=0; j<nbin; j++)
        {
            free(ind->gene[j]);
        }
        free(ind->xbin);
        free(ind->gene);
    }
    free(ind->obj);
    if (ncon != 0)
    {
        free(ind->constr);
    }
    return;
}
*/

void allocate_memory_pop_rediv (population *pop, int size, int obj_num)
{
	int i;
    pop->ind = (individual *)malloc(size*sizeof(individual));
    for (i=0; i<size; i++)
    {
         allocate_memory_ind (&(pop->ind[i]));
    }
    return;

}

void allocate_memory_ind_rediv (individual *ind, int obj_num)
{
    int j;
    if (nreal != 0)
    {
        ind->xreal = (double *)malloc(nreal*sizeof(double));
    }
    if (nbin != 0)
    {
        ind->xbin = (double *)malloc(nbin*sizeof(double));
        ind->gene = (int **)malloc(nbin*sizeof(int));
        for (j=0; j<nbin; j++)
        {
            ind->gene[j] = (int *)malloc(nbits[j]*sizeof(int));
        }
    }
    ind->obj = (double *)malloc(obj_num*sizeof(double));
    if (ncon != 0)
    {
        ind->constr = (double *)malloc(ncon*sizeof(double));
    }
    return;
}
