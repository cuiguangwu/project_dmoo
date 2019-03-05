/* Routines for storing population data into files */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "global2.h"
# include "rand.h"

/* Function to print the information of a population in a file */
void report_pop (population *pop, FILE *fpt, int size)
{
    int i, j;
    for (i=0; i<size; i++)
    {
        for (j=0; j<nobj-1; j++)
        {
            fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
		
        }
		fprintf(fpt,"%e\n",pop->ind[i].obj[nobj-1]);

    }
	
    return;
}
void report_time(double time, FILE *fpt)
{
	fprintf(fpt,"%f\n",time);
}

void report_best(population *pop, FILE *fpt)
{
	int i, j;
	for (i=0; i<popsize; i++)
	{
            for (j=0; j<nobj-1; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
            }
            fprintf(fpt,"%e\n",pop->ind[i].obj[nobj-1]);  
	}
}

void report_nondomi_set(population *pop)
{
	int i, j;
	FILE *fpt, *fpt1, *fpt2, *fpt3;
	list *temp, *temp1;
	int flag, *mark, sum, count, final_size;
	population *final_pop;
	double **front;
	int noPoints;
    fpt = fopen("nondominated_set.txt", "a");
	fpt1 = fopen("hypervolume.txt", "a");
 
	temp = beyond_ind->child;
	sum = 0;
	while (temp!=NULL)
	{
		sum ++;
		temp = temp->child;
	}
	mark = (int *)malloc(sum*sizeof(int));
	for (i=0; i<sum; i++)
	{
		mark[i] = 1;
	}
	
    i = 0;
	temp = beyond_ind->child;
	while (temp!=NULL)
	{
		temp1 = temp->child;
		j = i+1;
		while(temp1!=NULL&&mark[i]!=0)
		{
			if (mark[j]!=0)
			{
				flag = check_dominance (&pop->ind[temp->index], &pop->ind[temp1->index]);
				if (flag == -1||flag == 2)
				{
					mark[i] = 0;
				}
				if (flag == 1)
				{
					mark[j] = 0;
				}
			}
			j++;
			temp1 = temp1->child;
		}
		i++;
		temp = temp->child;		
	}
	final_size = 0;
    for (i=0; i<sum; i++)
	{
		if (mark[i]==1)
		    final_size++;	
	}

	final_pop = (population *)malloc(sizeof(population));                              // 评价并输出结果
    allocate_memory_pop (final_pop, final_size);

	

    count = 0;    
    temp = beyond_ind->child;
	for (i=0; i<sum; i++)
	{
		if (mark[i]==1)
		{
			for (j=0; j<nobj-1; j++)
			{
				fprintf(fpt,"%e\t", pop->ind[temp->index].obj[j]);
			}
//			fprintf(fpt, "%e\n", pop->ind[temp->index].obj[nobj-1]);
			if (count == final_size - 1)
				fprintf(fpt,"%e\t",pop->ind[temp->index].obj[nobj-1]);  
			else
				fprintf(fpt,"%e\t\n",pop->ind[temp->index].obj[nobj-1]);
			copy_ind(&pop->ind[temp->index], &final_pop->ind[count]);
			count++;
		}
		temp = temp->child;
	}
	fprintf(fpt,"\n\n\n");
	fflush (stdout);
	fflush (fpt);
    

	fpt3 = fopen("onebest.txt","w");
	report_best_final_hv (final_pop, fpt3, final_size);
	fflush(stdout);
	fflush(fpt3);
	fclose(fpt3);
	
	

	
	newscore_diversity(final_pop, final_size);
//	gd_igd_wfg(final_pop, final_size);
	score_pop(final_pop, final_size);
//	get_gd(final_pop, final_size);
//	get_igd(final_pop, final_size);
//	gd_igd_uf(final_pop, final_size);
//	gd_igd_cf(final_pop,final_size);
	score_extent(final_pop, final_size);
	find_extent(final_pop, final_size);
    ////////////////////////
	fpt2 = fopen("onebest.txt", "r");	
	noPoints = ReadFile(&front, fpt2, nobj);		
		
	if (nobj < 3)
		hypervolume = CalculateArea_2D(front, noPoints, nobj);
	else
		hypervolume = CalculateHypervolume(front, noPoints, nobj);
	if (noPoints == 1)
	{
		hypervolume = 0;
	}
		
	printf("\nThe hypervolume of the population is: %lf\n", hypervolume);
	fprintf(fpt1, "%lf\n", hypervolume);
    ////////////////////////

	
	fflush(fpt1);
	fflush(fpt2);
    fclose (fpt);
	fclose (fpt1);
	fclose (fpt2);
	free (mark);
	deallocate_memory_pop (final_pop, final_size);
	free (final_pop);
	return;

}
void report_best_final(population *pop, FILE *fpt, int final_size)
{
	int i, j;
	for (i=0; i<final_size; i++)
	{		
		for (j=0; j<nobj-1; j++)
		{
			fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
		}
		if (i == final_size - 1)
			fprintf(fpt,"%e\t",pop->ind[i].obj[nobj-1]);  
		else
			fprintf(fpt,"%e\t\n",pop->ind[i].obj[nobj-1]);  
   	}
}
void report_best_final_hv (population *pop, FILE *fpt, int final_size)
{
    int i, j, k, flag, count;
	int *mark;
	mark = (int *)malloc(final_size*sizeof(int));
	for (i=0; i<final_size; i++)
	{
		mark[i] = 1;
	}
	// delete repeat individuals
	for (i=0; i<final_size; i++)
	{
		for (j=i+1; j<final_size; j++)
		{
			if (check_repeat(&pop->ind[i], &pop->ind[j]) == 0)
			{
				mark[i] = 0;
				break;
			}
		}
	}
	// delete beyond reference individuals
	for (i=0; i<final_size; i++)
	{
		if (mark[i]==1)
		{
			for (j=0; j<nobj; j++)
			{
				if (point_ref[j]<pop->ind[i].obj[j])
				{
					mark[i] = 0;
					break;
				}
			}
		}
	}
	count = 0;
	for (i=0; i<final_size; i++)
	{
		if (mark[i]==1)
			count++;

	}
    k = 0;
    for (i=0; i<final_size; i++)
    {
        if (mark[i]==1)
		{
			for (j=0; j<nobj-1; j++)
			{
				fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
			}
			if (k == count - 1)
				fprintf(fpt,"%e\t",pop->ind[i].obj[nobj-1]);  
			else
				fprintf(fpt,"%e\t\n",pop->ind[i].obj[nobj-1]);
			k++;
			
		}
    }
	free (mark);
	
    return;
}


void report_nondomi_set_represent(population *pop, list *ultimate_list)
{
	int i, j;
	FILE *fpt, *fpt1, *fpt2;
	list *temp, *temp1;
	int flag, *mark, sum, count, final_size;
	population *final_pop;
	double **front;
	int noPoints;
    fpt = fopen("nondominated_set.txt", "w");
	fpt1 = fopen("hypervolume.txt", "a");
 
	temp = ultimate_list->child;
	sum = 0;
	while (temp!=NULL)
	{
		sum ++;
		temp = temp->child;
	}
	mark = (int *)malloc(sum*sizeof(int));
	for (i=0; i<sum; i++)
	{
		mark[i] = 1;
	}
	
    i = 0;
	temp = ultimate_list->child;
	while (temp!=NULL)
	{
		temp1 = temp->child;
		j = i+1;
		while(temp1!=NULL&&mark[i]!=0)
		{
			if (mark[j]!=0)
			{
				flag = check_dominance (&pop->ind[temp->index], &pop->ind[temp1->index]);
				if (flag == -1||flag == 2)
				{
					mark[i] = 0;
				}
				if (flag == 1)
				{
					mark[j] = 0;
				}
			}
			j++;
			temp1 = temp1->child;
		}
		i++;
		temp = temp->child;		
	}
	final_size = 0;
    for (i=0; i<sum; i++)
	{
		if (mark[i]==1)
		    final_size++;	
	}

	final_pop = (population *)malloc(sizeof(population));                              // 评价并输出结果
    allocate_memory_pop (final_pop, final_size);

	

    count = 0;    
    temp = ultimate_list->child;
	for (i=0; i<sum; i++)
	{
		if (mark[i]==1)
		{
			for (j=0; j<nobj-1; j++)
			{
				fprintf(fpt,"%e\t", pop->ind[temp->index].obj[j]);
			}
//			fprintf(fpt, "%e\n", pop->ind[temp->index].obj[nobj-1]);
			if (count == final_size - 1)
				fprintf(fpt,"%e\t",pop->ind[temp->index].obj[nobj-1]);  
			else
				fprintf(fpt,"%e\t\n",pop->ind[temp->index].obj[nobj-1]);
			copy_ind(&pop->ind[temp->index], &final_pop->ind[count]);
			count++;
		}
		temp = temp->child;
	}
	fflush (stdout);
	fflush (fpt);

	
	newscore_diversity(final_pop, final_size);
	score_pop(final_pop, final_size);	
	get_igd(final_pop, final_size);
	score_extent(final_pop, final_size);
	find_extent(final_pop, final_size);
    ////////////////////////
	fpt2 = fopen("D://环境程序/进化环境1.0/nondominated_set.txt", "r");	
	noPoints = ReadFile(&front, fpt2, nobj);		
		
	if (nobj < 3)
		hypervolume = CalculateArea_2D(front, noPoints, nobj);
	else
		hypervolume = CalculateHypervolume(front, noPoints, nobj);
		
	printf("\nThe hypervolume of the population is: %lf\n", hypervolume);
	fprintf(fpt1, "%lf\n", hypervolume);
    ////////////////////////

	
	fflush(fpt1);
	fflush(fpt2);
    fclose (fpt);
	fclose (fpt1);
	fclose (fpt2);
	free (mark);
	deallocate_memory_pop (final_pop, final_size);
	free (final_pop);
	return;

}

void report_best_list(population *pop, list *node, FILE *fpt)
{
	int i, j;
	list *temp, *temp1;
	int flag, *mark, sum, final_size;
	
	temp = node->child;
	sum = 0;
	while (temp!=NULL)
	{
		sum ++;
		temp = temp->child;
	}
	mark = (int *)malloc(sum*sizeof(int));
	for (i=0; i<sum; i++)
	{
		mark[i] = 1;
	}
	
    i = 0;
	temp = node->child;
	while (temp!=NULL)
	{
		temp1 = temp->child;
		j = i+1;
		while(temp1!=NULL&&mark[i]!=0)
		{
			if (mark[j]!=0)
			{
				flag = check_dominance (&pop->ind[temp->index], &pop->ind[temp1->index]);
				if (flag == -1||flag == 2)
				{
					mark[i] = 0;
				}
				if (flag == 1)
				{
					mark[j] = 0;
				}
			}
			j++;
			temp1 = temp1->child;
		}
		i++;
		temp = temp->child;		
	}
	final_size = 0;
    for (i=0; i<sum; i++)
	{
		if (mark[i]==1)
		    final_size++;	
	}

    temp = node->child;
	for (i=0; i<sum; i++)
	{
		if (mark[i]==1)
		{
			for (j=0; j<nobj-1; j++)
			{
				fprintf(fpt,"%e\t", pop->ind[temp->index].obj[j]);
			}
			fprintf(fpt,"%e\t\n",pop->ind[temp->index].obj[nobj-1]);
		}
		temp = temp->child;
	}
	free (mark);
	return;
}

void report_all(population *pop, FILE *fpt)
{
	list *temp, *temp1;
	int i, j, k;

    fprintf(fpt, "\n\n !!! generation is %d!", ngen);

	fprintf(fpt, "\nthe individuals in population are\n");
	for (i=0; i<6*popsize; i++)
	{
		for (j=0; j<nobj-1; j++)
		{
			fprintf(fpt, "%e\t", pop->ind[i].obj[j]);
		}
		fprintf(fpt, "%e\t%d\n", pop->ind[i].obj[nobj-1], pop->ind[i].boundary_ind);
	}
	fprintf(fpt, "\nthe individuals in beyond_ind are\n");
	temp = beyond_ind->child;
	while (temp!=NULL)
	{
		for (j=0; j<nobj-1; j++)
		{
			fprintf(fpt,"%e\t", pop->ind[temp->index].obj[j]);
		}
		fprintf(fpt, "%e\t%d\n", pop->ind[temp->index].obj[nobj-1], temp->index);
		temp = temp->child;	
	}

//	fprintf(fpt, "\nthe individuals in empty_ind are\n");

/*	temp = empty_ind->child;
	while (temp!=NULL)
	{
		for (j=0; j<nobj-1; j++)
		{
			fprintf(fpt,"%e\t", pop->ind[temp->index].obj[j]);
		}
		fprintf(fpt, "%e\t%d\n", pop->ind[temp->index].obj[nobj-1], temp->index);
		temp = temp->child;	
	}
*/	
	fprintf(fpt, "\nthe individuals in environment are\n");
	temp = environment_list->child;
	while (temp!=NULL)
	{
	    temp1 = area[temp->index].area_ind->child;
		while(temp1!=NULL)
		{
			for (j=0; j<nobj-1; j++)
			{
			    fprintf(fpt,"%e\t", pop->ind[temp1->index].obj[j]);
			}
		    fprintf(fpt, "%e\t%d\n", pop->ind[temp1->index].obj[nobj-1], pop->ind[temp1->index].boundary_ind);
			temp1=temp1->child;
		}
         temp = temp->child; 		
	}
    fprintf(fpt, "\nthe information environment are\n");
	fprintf(fpt, "the environment_min and the environment_max are\n");
	for (j=0; j<nobj; j++)
	{
		fprintf(fpt, "%e\t%e\n", environment_min[j], environment_max[j]);
	}

	fprintf(fpt, "\nthe individuals count in environment is %d\n", environment_indcount);
    fprintf(fpt, "\nthe environment div are\n");
	for (j=0; j<nobj; j++)
	{
		fprintf(fpt, "%e\n", environment_distance[j]);
	}
    fprintf(fpt, "\nthe boundary ind in each_objective are\n");
	for (j=0; j<nobj; j++)
	{
		fprintf(fpt, "the boundary ind in %d objective is\n", j);
		for (k=0; k<nobj; k++)
		{
			fprintf(fpt, "%e\t", pop->ind[extreme_ind[j]].obj[k]);
		}
		fprintf(fpt, "%d\n", extreme_ind[j]);
	}
    fprintf(fpt, "\nthe information of area in environment are\n");
	temp = environment_list->child;
	while (temp!=NULL)
	{
		fprintf(fpt, "\nthe number of area in environment is %d\n", temp->index);
        fprintf(fpt, "the individual count in this area  is %d\n", area[temp->index].area_indcount);
		fprintf(fpt, "the non-dominated area number for this area  is %d\n", area[temp->index].area_nondominated_number);
		fprintf(fpt, "the non-dominated area for this area are\n");
		temp1 = area[temp->index].area_non_domianted->child;
		while (temp1!=NULL)
		{
			fprintf(fpt, "->%d", temp1->index);
			temp1 = temp1->child;
		}

		fprintf(fpt, "\nthe reprensent individual  in this area  is\n");
		for (j=0; j<nobj; j++)
		{
			fprintf(fpt, "%e\t", pop->ind[area[temp->index].area_repreasent_ind].obj[j]);
		}
        fprintf(fpt, "\nthe individuals in this area are\n");        
		temp1 = area[temp->index].area_ind->child;
		
		while(temp1!=NULL)
		{
		    for (j=0; j<nobj; j++)
			{
			    fprintf(fpt, "%e\t", pop->ind[temp1->index].obj[j]);
			}
			fprintf(fpt, "\n");
			temp1 = temp1->child;
		}
		temp = temp->child;		
	}







}
/* Function to print the information of feasible and non-dominated population in a file */
void report_feasible (population *pop, FILE *fpt)
{
    int i, j, k;
    for (i=0; i<popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)
        {
            for (j=0; j<nobj; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
            }
             if (ncon!=0)
            {
                for (j=0; j<ncon; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
                }
            }
            if (nreal!=0)
            {
                for (j=0; j<nreal; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].xreal[j]);
                }
            }
            if (nbin!=0)
            {
                for (j=0; j<nbin; j++)
                {
                    for (k=0; k<nbits[j]; k++)
                    {
                        fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                    }
                }
            }
            fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
            fprintf(fpt,"%d\t",pop->ind[i].rank);
            fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
        }
    }
    return;
}

void report_best_reduce(population *pop, FILE *fpt, int obj_num)
{
	int i, j;
	for (i=0; i<popsize; i++)
	{
        for (j=0; j<obj_num-1; j++)
        {
            fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
        }
        fprintf(fpt,"%e\n",pop->ind[i].obj[obj_num-1]);  
 	}
    fprintf(fpt,"\n");
}

void report_best_extent(population *pop, FILE *fpt, int count)
{
	int i, j;
	for (i=0; i<count; i++)
	{
        for (j=0; j<nobj-1; j++)
        {
            fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
        }
        fprintf(fpt,"%e\n",pop->ind[i].obj[nobj-1]);  
 	}
    fprintf(fpt,"\n");
}


