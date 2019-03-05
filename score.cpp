# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"


//计算MS和SP
void score_pop(population *pop_final, int final_size)
{
		//SP
	double d_avg, sp, d, p, sum=0, gd, sum1, temp1, sum2, ms;
    int i, j, h, g, m, k;
	double *min, *min_ms, *max_ms;
	double s1, s2;
	int num_criterion=500;
	FILE *sp_ptr;
	FILE *ms_ptr;

	min = (double *)malloc(final_size*sizeof(double));
	min_ms = (double *)malloc(nobj*sizeof(double));
	max_ms = (double *)malloc(nobj*sizeof(double));


	

	sp_ptr = fopen("sp.txt", "a");
	ms_ptr = fopen("ms.txt", "a");
	for (m=0; m<nobj; m++)
	{
		min_ms[m] = INF;
		max_ms[m] = -100000;

	}
    
	for (h=0; h<final_size; h++)
	{
		for (m=0; m<nobj; m++)
		{
			if (pop_final->ind[h].obj[m]<min_ms[m])
				min_ms[m]= pop_final->ind[h].obj[m];
			if (pop_final->ind[h].obj[m]>max_ms[m])
				max_ms[m]= pop_final->ind[h].obj[m];

		}
	}
	sum2 = 0;
	for (m=0; m<nobj; m++)
	{
		sum2 += pow(max_ms[m]-min_ms[m], 2);
	}
	ms = sqrt(sum2);
	fprintf(ms_ptr,"%f\n",ms);



	for (h=0; h<final_size; h++)
	{  	
		min[h]=INF;//
		for (g=0; g<final_size; g++)
		 if (h!=g)
		 {
            d=0;
			for(m=0; m<nobj; m++)
			{
			  d+=(pop_final->ind[h].obj[m]-pop_final->ind[g].obj[m]>0?pop_final->ind[h].obj[m]-pop_final->ind[g].obj[m]:pop_final->ind[g].obj[m]-pop_final->ind[h].obj[m]);
			}
         //   min[h]+=d ;

			if (d<min[h])
				min[h]=d;
		 }
		//min[h]=min[h]/final_size;
		sum +=min[h];
	}
	d_avg = sum/final_size;
	sp=0; 
	for (h=0; h<final_size; h++)
		sp +=(d_avg-min[h])*(d_avg-min[h]);
    sp =sqrt(sp/(final_size-1));
    fprintf(sp_ptr,"%f\n",sp);
		                                
		                                
	                                              
	

	fclose(sp_ptr);

	fclose(ms_ptr);

	free(min);
	free(min_ms);
	free(max_ms);
}

//计算HV
void score_extent(population *pop_final, int final_size)
{
	int i, j, k, temp, m, betterInAnyObjective, count;
	double weightheart[10], ref[10], value[10], sum, weighheartvolume, sum_all, temp1, temp2, temp3, max[10], min[10];
	population* pop_rediv, *pop_nondomi;
	FILE *extent_ptr, *fpt;
	pop_rediv = (population *)malloc(sizeof(population));
    pop_nondomi = (population *)malloc(sizeof(population));
    allocate_memory_pop_rediv (pop_rediv, final_size, nobj-1);
	allocate_memory_pop_rediv (pop_nondomi, final_size, nobj-1);

	extent_ptr = fopen("extent.txt", "a");
//	fpt = fopen("reducediv.txt", "w");
	
	for (k=0; k<nobj; k++)
	{
		weighheartvolume = 1;
		for (i=0; i<final_size; i++)
		{
		    for (j=0; j<nobj-1; j++)
			{
			    temp=k+j+1;
				if (temp>=nobj)
				{
					temp = temp%nobj;
				}
    
				pop_rediv->ind[ i].obj[j]=pop_final->ind[i].obj[temp];
			}
		}
		// calculate weightheart
/*		for (j=0; j<nobj-1; j++)
		{
			sum = 0;
			for (i=0; i<final_size; i++)
			{
				sum += pop_rediv->ind[i].obj[j];
			}
			weightheart[j] = sum/final_size;

		}
*/		
		// calculate reference point
		temp1 = 1;
        for (j=1; j<nobj; j++)
		{
			temp1 *=j;
		}
		temp1 = final_size * temp1;
		temp2 = 1.0/(nobj-1.0);
		temp1 = pow(temp1, temp2);
		
        for (j=0; j<nobj-1; j++)
		{
			max[j] = 0;
			min[j] = 1000000;
			for (i=0; i<final_size; i++)
			{
				if (pop_rediv->ind[i].obj[j]>max[j])
				{
					max[j]=pop_rediv->ind[i].obj[j];
				}
				if (pop_rediv->ind[i].obj[j]<min[j])
				{
					min[j]=pop_rediv->ind[i].obj[j];
				}
			}
		    temp3 = (max[j]-min[j])/temp1;
			ref[j] = min[j] - temp3;
			
		}
		// calculate volume
/*		for (j=0; j<nobj-1; j++)
		{
			weighheartvolume *= (weightheart[j] - ref[j]);
		}
*/		

//		report_best_reduce(pop_rediv,fpt,nobj-1);

		value[k] = calhypervolume(pop_rediv, final_size, nobj-1);
		if (nobj==2)
			value[k] -= ref[0];
		if (nobj==3)
		{
			value[k] -= ((max[0]-ref[0])*ref[1]+max[1]*ref[0]);
		}
		if (nobj==4)
		{
			value[k] -= ((max[0]*max[1]*ref[2]+max[1]*max[2]*ref[0]+max[2]*max[0]*ref[1])-(max[0]*ref[1]*ref[2]+max[1]*ref[0]*ref[2]+max[2]*ref[1]*ref[0])+ref[0]*ref[1]*ref[2]);
		}

 
//		value[k] = value[k]/weighheartvolume;
//		fflush(fpt);

	}
	sum_all = 1;
	for (k=0; k<nobj; k++)
	{
		sum_all *= value[k];
	}
	sum_all = pow(sum_all, 1.0/nobj);
	fprintf(extent_ptr,"%f\n",sum_all);
//	fflush(fpt);
//	fclose(fpt);
	fclose(extent_ptr);
    deallocate_memory_pop (pop_nondomi, final_size);
	deallocate_memory_pop (pop_rediv, final_size);
	free(pop_nondomi);
	free(pop_rediv);
	return;

}

/*
void hypervolume(population *pop_final)
{
	int i,j;
	double value;
	FILE *hv_ptr;
    

	int true_max[20], true_min[20];

	hv_ptr = fopen("hv.txt", "a");
	// 真实pareto面每维最大，最小值
	// ZDT except ZDT3
	for (i=0; i<nobj; i++)
	{
		true_max[i]=0.5;
		true_min[i]=0.0;
	}
	for (i=0; i<final_size; i++)
	{
		for(j=0; j<nobj; j++)
		{
			pop_final->ind[i].obj[j]=(pop_final->ind[i].obj[j]-true_min[j])/(true_max[j]-true_min[j]);
		}
	}
    for (i=0; i<final_size; i++)
	{
		for(j=0; j<nobj; j++)
		{
			if (pop_final->ind[i].obj[j]<=1.0&&pop_final->ind[i].obj[j]>=0.0)
			{
				pop_final->ind[i].obj[j]=1.0-pop_final->ind[i].obj[j];
			}
			else
			{
				if (pop_final->ind[i].obj[j]>1.0)
					pop_final->ind[i].obj[j]=0;
				else 
					pop_final->ind[i].obj[j]=1.0;
			}

		}
	}
	value = calhypervolume(pop_final, final_size, nobj);
	fprintf(hv_ptr,"%f\n",value);
	fclose(hv_ptr);
}
*/
//计算HV
double calhypervolume(population *pop_final, int ind_num, int obj_num)
{
	int i, j;
	int n, nondomi_num;
	double volume=0, distance=0, tempvolume, tempdistance, max;
	n = ind_num;
	if (obj_num!=1)
	{
	    while (n>0)
		{
		    nondomi_num = filter_nondomi(pop_final, n, obj_num-1);
		    tempvolume = 0;
            if (obj_num < 3) 
			{
			    tempvolume = pop_final->ind[0].obj[0];
			} 
		    else
			    tempvolume = calhypervolume(pop_final, nondomi_num, obj_num - 1);

            tempdistance = surfaceUnchangedTo(pop_final, n, obj_num - 1);
            volume += tempvolume * (tempdistance - distance);
            distance = tempdistance;
            n = reduceNondominatedSet(pop_final, n, obj_num - 1, distance);
		}
	}
	else
	{
		max = 0;
		for (i=0; i<n; i++)
		{
			if (pop_final->ind[i].obj[0]>max)
			max = pop_final->ind[i].obj[0];
		}
		volume = max;
	}
	return volume;

}

//计算HV：挑选出第一层个体
int filter_nondomi(population *pop_final, int ind_num, int obj_num)
{
	int i,j;
	int n, temp;
	n = ind_num;
	i = 0;
	while (i < n)
	{
		j=i+1;
		while (j < n) 
		{
			
			temp = dominates(&(pop_final->ind[i]), &(pop_final->ind[j]), obj_num);
			if (temp == 1) 
			{
	/* remove point 'j' */
	            n--;
	            swap(pop_final, j, n);
			} 
			else
			{
				temp = dominates(&(pop_final->ind[j]), &(pop_final->ind[i]), obj_num);
			    if (temp == 1) 
				{
	/* remove point 'i'; ensure that the point copied to index 'i'
	   is considered in the next outer loop (thus, decrement i) */
	                n--;
	                swap(pop_final, i, n);
	                i--;
	                break;
				} 
				else
	                j++;
			}
		}
		i++;
	}
	return (n);

}
//计算HV：若个体1支配（>关系）个体2，则返回1
int dominates(individual *ind1, individual *ind2, int obj_num) 
{
    int i;
    int betterInAnyObjective;
	int mark;

    betterInAnyObjective = 0;
    for (i = 0; i < obj_num && ind1->obj[i] >= ind2->obj[i]; i++)  //
        if (ind1->obj[i] > ind2->obj[i])
	    betterInAnyObjective = 1;
	if ((i >= obj_num) && (betterInAnyObjective>0))//???????????????????????????????????????????????????
		mark = 1;
	else
		mark = 0;
    
    return (mark);
} //Dominates 
 
//计算HV：在种群中交换两个个体（的位置）
void  swap(population *pop_final, int  i, int  j)
{
    individual *temp;
	temp = (individual*)malloc(sizeof(individual));
	allocate_memory_ind (temp);
    
    copy_ind(&(pop_final->ind[i]), temp);
	copy_ind(&(pop_final->ind[j]), &(pop_final->ind[i]));
	copy_ind(temp, &(pop_final->ind[j]));
	deallocate_memory_ind (temp);
	free(temp);
    
  } // Swap 


//计算HV：求出obj[obj_num]最小值
double  surfaceUnchangedTo(population *pop_final, int ind_num, int obj_num)
{
    int i;
    double  minValue, value;

    minValue = pop_final->ind[0].obj[obj_num];
		
    for (i = 1; i < ind_num; i++) 
	{
        value = pop_final->ind[i].obj[obj_num];
		if (value < minValue)
            minValue = value;
    }
    return minValue;
} // SurfaceUnchangedTo 


//计算HV：保留obj[obj_num]<=threshold的个体
int  reduceNondominatedSet(population *pop_final, int ind_num, int obj_num, double threshold)
{
    int  n;
    int  i;

    n = ind_num;
    for (i = 0; i < n; i++)

        if (pop_final->ind[i].obj[obj_num]<=threshold)
		{
            n--;
            swap(pop_final, i, n);
		}
  
    return n;
} // ReduceNondominatedSet

//
void find_extent(population *pop_final, int final_size)
{
	int i, j, k, m, number, n, count = 0, betterInAnyObjective;
    population* pop_extent;
	int tag[1000];
	FILE *fpt;
	pop_extent = (population *)malloc(sizeof(population));
    allocate_memory_pop (pop_extent, final_size);
	fpt = fopen("out_extent.txt", "w");
	for (i=0; i<final_size; i++)
	{
		tag[i] = 0;
	}

	for(k=0; k<nobj; k++)
	{
	    for(i=0; i<final_size; i++)
		{
			for(j=0; j<final_size; j++)
			{
				m = k+1;
				if (m>=nobj)
				{
					m = m % nobj;
				}
				
				
				betterInAnyObjective = 0;
                for (n = 0; n < nobj-1 && pop_final->ind[j].obj[m] >= pop_final->ind[i].obj[m]; n++)
				{
                    if (pop_final->ind[j].obj[m] > pop_final->ind[i].obj[m])
	                    betterInAnyObjective = 1;
					m++;
					if (m>=nobj)
					{
					    m = m % nobj;
					}
				}
	            if ((n >= nobj-1) && (betterInAnyObjective>0))
		            break;
			}
			if (j>=final_size)
			{
// 全部边界解？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
				if (tag[i]==0)
// 部分边界解
//				if (tag[i]==0&&k==1)
				{
					copy_ind (&pop_final->ind[i], &pop_extent->ind[count]);
				    tag[i] = -1;
					count++;
				}
				
			}

		}
	}
	fflush(fpt);

	fclose(fpt);

	deallocate_memory_pop (pop_extent, final_size);
	free(pop_extent);

}



