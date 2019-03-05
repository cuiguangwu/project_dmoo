// environment setting
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"


//�жϸ����Ƿ���Խ��뵥Ԫ�򣺽�pop[index]��area[number]�еĸ�����бȽϣ�֧����߱�֧�䣩,numberΪ��Ԫ���index
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
		case 2:					//a��b���
			{
				return -1;
			}
		case 0:					//a��b����֧��
			{
				temp = temp->child;
				break;
			}
		}			
		}
		return 1;
	}
}


//�Ƚϣ��Ƚ�������Ԫ�򣬱Ƚϱ�׼ indcount����������>nondominated_number����֧�䵥Ԫ������
//1��node is better than node1
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


//�������pop[index]����Ԫ��area[number]�ľ���
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


//��pop[index]ֱ�ӷ��뻷��area[number]��
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

//ɾ����Ԫ��
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


//����a��b֮���fitness
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
