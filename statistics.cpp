#include "statistics.h"
#include <math.h>

double mean(double *arr, int size)
{
	double sum=0;
	for( int i=0; i<size;i++)
		sum =sum + arr[i];
	return double(sum/size);
}

double median(double *arr, int size)
{
	if(size==1)
		return arr[0];
	else if ((size+1)%2)
		return arr[(size+1)/2];
	else
		return 0.5*(arr[size/2] + arr[1+(size/2)]);
}

double variance(double *pop_arr, int size)
{
	double ave = 0;double sum = 0;
	ave = mean(pop_arr, size);
	
	for(int i=0; i<size; i++)
		sum = sum + pow( (pop_arr[i]-ave), 2 );
	
	return sum = double(sum/(size-1));
}

double min(double A, double B)//copy
{
	if(A<B)
		return A;
	else return B;
}


double max(double A, double B)//copy
{
	if(A>B)
		return A;
	else return B;
}