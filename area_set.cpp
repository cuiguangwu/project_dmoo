// environment setting
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"


//判断个体是否可以进入单元域：将pop[index]与area[number]中的个体进行比较（支配或者被支配）,number为单元域的index
int judge_area(population *pop, int index, int number)
{
	list *temp;
	int i;
	int flag;
	if (area[number].area_exist==0)
	{
		printf("Judge_area:this area is not exist!");
		return -1;
	}
	else
	{
		temp = area[number].area_ind->child;
		while (temp!=NULL)
		{
			flag = check_dominance (&pop->ind[index], &pop->ind[temp->index]);
		switch(flag)
		{
		case 1:
			{
				return 1;
			}
		case -1:
			{
				return -1;
			}
		case 2:					//a与b相等
			{
				return -1;
			}
		case 0:					//a与b互不支配
			{
				temp = temp->child;
				break;
			}
		}			
		}
		return 1;
	}
}


//比较：比较两个单元域，比较标准 indcount（个体数）>nondominated_number（非支配单元域数）
//1：node is better than node1
int compare_area(list *node, list *node1)
{
	if (area[node->index].area_indcount<area[node1->index].area_indcount)
	{
		return 1;
	}
	else 
	{
		if (area[node->index].area_indcount==area[node1->index].area_indcount&&area[node->index].area_nondominated_number>area[node1->index].area_nondominated_number)
			return 1;
		else 
			return 0;
	}
}


//计算个体pop[index]到单元域area[number]的距离
double cal_dist_origin(population *pop, int index, int number)
{
	int j;
	int *a_coordinate;
	int remainder;
	double distance;

		
    a_coordinate = (int *)malloc(nobj*sizeof(int));
	remainder = number;
	for (j=0; j<nobj; j++)
	{
		if (j == nobj-1)
			a_coordinate[j] = remainder;
		else
		{
			a_coordinate[j] = (int)floor(remainder/pow((double)environment_div, nobj-1-j));
		    remainder = remainder%(int)pow((double)environment_div, nobj-1-j);
		}			
	}
	distance = 1;
	for (j=0; j<nobj; j++)
	{
		distance = distance*(pop->ind[index].obj[j]-(a_coordinate[j]*environment_distance[j])); 
	}
	distance = pow(distance, 1.0/nobj);
//	printf("\n!! distance is %f", distance);
    free(a_coordinate);
	return distance;
}


//将pop[index]直接放入环境area[number]中
void enter_directly_area(int index, int number)
{
	list *temp;
	area[number].area_exist=1;
	area[number].area_indcount=1;
	temp=area[number].area_ind;
	insert(temp, index);
	area[number].area_repreasent_ind=index;
	return;
}

//删除单元域
void del_area(int number)
{
	list *temp, *temp1;
	area[number].area_exist=0;
	area[number].area_indcount=0;
	area[number].area_repreasent_ind=-1;
	
   
	
	temp = area[number].area_ind->child;	
	while (temp!=NULL)
	{
		temp = del(temp);
	    temp = temp->child;
	}
	
	temp = area[number].area_non_domianted->child;
	
	while (temp!=NULL)
	{
		area[temp->index].area_nondominated_number--;
		temp1 = findnode (area[temp->index].area_non_domianted->child, number);
		temp1 = del (temp1);
		temp = del(temp);
	    temp = temp->child;
	}
	area[number].area_nondominated_number=0;
    return; 	
}


//计算a与b之间的fitness
double cal_fitness(int a, int b)
{
	int j;
	int *a_coordinate, *b_coordinate;
	int remainder;
	double max, sum;
	if (a==b)
	{
		printf("= is not strict dominated");
		exit(1);
	}
		
    a_coordinate = (int *)malloc(nobj*sizeof(int));
	b_coordinate = (int *)malloc(nobj*sizeof(int));
	remainder = a;
	for (j=0; j<nobj; j++)
	{
		if (j == nobj-1)
			a_coordinate[j] = remainder;
		else
		{
			a_coordinate[j] = (int)floor(remainder/pow((double)environment_div, nobj-1-j));
		    remainder = remainder%(int)pow((double)environment_div, nobj-1-j);
		}			
	}
	remainder = b;
	for (j=0; j<nobj; j++)
	{
		if (j == nobj-1)
			b_coordinate[j] = remainder;
		else
		{
			b_coordinate[j] = (int)floor(remainder/pow((double)environment_div, nobj-1-j));
		    remainder = remainder%(int)pow((double)environment_div, nobj-1-j);
		}			
	}
/*	
    max = 0;
	for (j=0; j<nobj; j++)
	{
		if ((a_coordinate[j]-b_coordinate[j])>max)
			max = a_coordinate[j]-b_coordinate[j];
	}
*/
	sum = 0;
	for (j=0; j<nobj; j++)
	{
		if ((a_coordinate[j]-b_coordinate[j])>0)
			sum += a_coordinate[j]-b_coordinate[j];
		else
		{
			sum += 1.0/(2.0+b_coordinate[j]-a_coordinate[j]);
		}
	}
	
	free(a_coordinate);
	free(b_coordinate);
//	return max;
	return sum;

}
