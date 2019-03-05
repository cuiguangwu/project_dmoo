/* Domination checking routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Routine for usual non-domination checking
   It will return the following values
   1 if a dominates b
   -1 if b dominates a
   0 if both a and b are non-dominated 
   2 if a = b
*/

//非网格下的支配关系判断   若返回1则表示a支配b
int check_dominance (individual *a, individual *b)
{
    int i, end;
    int flag1;
    int flag2;
    flag1 = 0;
    flag2 = 0;
	end = 0;

	for (i=0; i<nobj; i++)
	{
		if (fabs(a->obj[i]-b->obj[i]) > EPS)
			end = 1;
	}
	if (end == 0)
	{
		return 2;
	}
    if (a->constr_violation<0 && b->constr_violation<0)
    {
        if (a->constr_violation > b->constr_violation)
        {
            return (1);
        }
        else
        {
            if (a->constr_violation < b->constr_violation)
            {
                return (-1);
            }
            else
            {
                return (0);
            }
        }
    }
    else
    {
        if (a->constr_violation < 0 && b->constr_violation == 0)
        {
            return (-1);
        }
        else
        {
            if (a->constr_violation == 0 && b->constr_violation <0)
            {
                return (1);
            }
            else
            {
                for (i=0; i<nobj; i++)
                {
                    if (a->obj[i] < b->obj[i])
                    {
                        flag1 = 1;

                    }
                    else
                    {
                        if (a->obj[i] > b->obj[i])
                        {
                            flag2 = 1;
                        }
                    }
                }
                if (flag1==1 && flag2==0)
                {
                    return (1);
                }
                else
                {
                    if (flag1==0 && flag2==1)
                    {
                        return (-1);
                    }
                    else
                    {
                        return (0);
                    }
                }
            }
        }
    }
}
// a dominate b return 1; a weak dominate b return 2; b dominate a return -1; b weak dominate a return -2; a = b return 3; a non-dominate b return 0
int check_area_dominance (int a, int b)
{
	int j, flag1, flag2;
	int *a_coordinate, *b_coordinate;
	int remainder;
	if (a==b)
	{
		return 3;
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
	flag1 = 0;
	flag2 = 0;
    for (j=0; j<nobj; j++)
    {
        if (a_coordinate[j] < b_coordinate[j])
        {
            flag1 = 1;
        }
        else
        {
            if (a_coordinate[j] > b_coordinate[j])
            {
                flag2 = 1;
            }
        }
    }
    if (flag1==1 && flag2==0)
    {
        for (j=0; j<nobj; j++)
		{
			if (a_coordinate[j] == b_coordinate[j])
			{
				free(a_coordinate);
				free(b_coordinate);
				return 2;
			}
		}
		free(a_coordinate);
		free(b_coordinate);
		return 1;
    }
	if (flag1==1 && flag2==1)
	{
		free(a_coordinate);
		free(b_coordinate);
		return 0;
	}
	
    for (j=0; j<nobj; j++)
	{
		if (a_coordinate[j] == b_coordinate[j])
		{
			free(a_coordinate);
			free(b_coordinate);
			return -2;
		}	
	}
	free(a_coordinate);
	free(b_coordinate);
	return -1;
 
}

/* 检查两个个体是否相等，是则返回0 */
int check_repeat (individual *a, individual *b)
{
	int i;
    int flag;
    flag = 0;   
	for (i=0; i<nobj; i++)
	{
		if (a->obj[i]!=b->obj[i])
			flag = 1;
	}
	if (flag==0)
	{
		return (0);
	}
	return (1);
}


