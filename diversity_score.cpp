# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"



void newscore_diversity(population *pop_new, int final_size) //
{
	double sp, sum=0, sum1=0, countsum, value, value2;
	int i, j, k, *b, *d, count1, count2;
	double temp=0, temp1, min2, *lowcost, *e, valuesum;
	int *closest;
	FILE *sp_ptr1;

	double **c;




	lowcost = (double *)malloc(final_size*sizeof(double));
	closest = (int *)malloc(final_size*sizeof(int));	


	b = (int *)malloc(final_size*sizeof(int));
	d = (int *)malloc(final_size*sizeof(int));
	e = (double *)malloc(final_size*sizeof(double));
	sp_ptr1 = fopen("mstend_sp.txt", "a");






	c = (double **)malloc(final_size*sizeof(double));

	for (i=0; i<final_size; i++)
	{
		c[i] = (double *)malloc(final_size*sizeof(double));
	}
	for (i=0; i<final_size; i++)
	{   

		for(j=i+1; j<final_size; j++)
		{
			temp1=0;
			for(k=0; k<nobj; k++)
			{
				temp1+=(pop_new->ind[i].obj[k]-pop_new->ind[j].obj[k])*(pop_new->ind[i].obj[k]-pop_new->ind[j].obj[k]);
			}

			temp1=sqrt(temp1);


			c[i][j]=temp1;
			c[j][i]=temp1;
		}
		c[i][i]=INF;
	}


	for(i=0; i<final_size; i++)
	{
		lowcost[i]=c[0][i];
		closest[i]=0;
	}
	closest[0]=-1;


	for(i=0; i<final_size-1; i++)
	{
		min2=INF;
		k=0;
		for(j=0; j<final_size; j++)
			if(closest[j]!=-1&&lowcost[j]<min2)
			{
				min2=lowcost[j];
				k=j;
			}
			if(k)
			{			
				b[i]=closest[k];
				d[i]=k;
				e[i]=min2;


				closest[k]=-1;
				for(j=1;j<final_size;j++)
					if(closest[j]!=-1&&c[k][j]<lowcost[j])
					{
						lowcost[j]=c[k][j];
						closest[j]=k;
					}
			}
	}
	countsum=0;
	valuesum=0;
	for(i=0; i<final_size-1; i++)
	{   
		count1=0;
		count2=0;
		value=1;
		value2=0;
		for(j=0; j<final_size; j++)
		{   
			if(c[b[i]][j]<e[i])
			{
				value*=(c[b[i]][j]/e[i]);
				count1++;
			}
			if(c[d[i]][j]<e[i])
			{
				value*=(c[d[i]][j]/e[i]);
				count2++;
				value2=1;
			}
		}
		if (count1==0&&count2==0)
			countsum+=1;
		valuesum+=pow(value,1);


	}

	value2=final_size-1;
	valuesum=(valuesum-countsum)/(value2-countsum);
	valuesum=pow(valuesum,1);
	sp=valuesum;


	fprintf(sp_ptr1,"%.12f\n",sp);

	fclose(sp_ptr1);


	free(b);
	free(d);
	free(e);




	for (i=0; i<final_size; i++)
	{
		free(c[i]);
	}
	free(c);
	free(closest);
	free(lowcost);
}

