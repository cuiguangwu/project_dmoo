/* this random generator is copied from genesis by Grefenstette. */
/* random number generator (copied from Grefensttete's genesis) */

#ifndef _RANDOM_H
#define _RANDOM_H

#include <time.h>
#include <stdlib.h>
#include <math.h>

#define MASK 2147483647
#define PRIME 65539
#define SCALE 0.4656612875e-9
#define DEFAULT_SEED   123456789

class Rand {
public:
	Rand() { m_seed = DEFAULT_SEED; }
	Rand(long sd) { seed(sd); }
	~Rand() {}

protected:
	unsigned long m_seed;

public:
	/**
	/* Set the seed for random generator
	*/
	void seed(long sd) {
		if(sd<0) sd=(unsigned)time( NULL );
		sd=sd & MASK;
		m_seed=sd;
	}

	/**
	/* Return the seed for random generator
	*/
	unsigned long seed() { return m_seed; }

	/**
	* Returns a psuedo-random double value in a uniform distribution in [0,1)
	* which is between 0 and 1, excluding 1.
	*/
	double unit() {
		m_seed = ( (m_seed * PRIME) & MASK);
		return m_seed * SCALE;
	}

	/**
	* Returns a random integer in a uniform distribution in [low,high]
	* @param low  Low bound
	* @param high Upper bound
	*/
	int integer(int low, int high) {
		int ret;
		ret = (int) ((low) + ((high)-(low)+1) * unit());
		return ret;
	}

	/**
	* Returns an uniform double value between a range.
	* @param _lower  The lowest value
	* @param _upper  The uppest value
	*/
	double uniform(double _lower, double _upper) {
		return (_upper-_lower)*unit()+_lower;
	}

	/**
	* Returns a random value in a normal distribution
	* @param _mean  Dsitribution mean
	* @param _sd    Standard Deviation
	*/
	double normal(double _mean=0.0, double _sd=1.0)
	{
		double rsquare, factor, var1, var2;
		do{
		var1 = 2.0 * unit() - 1.0;
		var2 = 2.0 * unit() - 1.0;
		rsquare = var1*var1 + var2*var2;
		} while(rsquare >= 1.0 || rsquare == 0.0);

		double val = -2.0 * log(rsquare) / rsquare;
		if(val > 0.0) factor = sqrt(val);
		else           factor = 0.0;	// should not happen, but might due to roundoff

		return (var1 * factor * _sd) + _mean;
	}

	/**
	* Returns a random value in a negative exponential distribution
	* @param _mean  Dsitribution mean
	*/
	double negExp(double _mean)
	{
		return -_mean*log(unit());
	}

};

///////////////////////////////////
///////////////////////////////////
/* random number generator */
/* this random generator is implemented with the help of function rand() in math library. */
/* ??the problem is that we can only have one instance of this class because of rand()*/
class Rand2{
public:
	Rand2() {srand(6764421);}
	Rand2(int sd) { seed(sd); }
	void seed(int sd) {
		if(sd<0) sd=(unsigned)time( NULL );
		srand(sd);
	}
	
	double unit() {
		int v=rand();
		return v/(RAND_MAX+1.0);
	} 

	double uniform(double _lower, double _upper) {
		return (_upper-_lower)*unit()+_lower;
	}

	double normal(double _mean=0.0, double _sd=1.0)
	{
		double rsquare, factor, var1, var2;
		do{
		var1 = 2.0 * unit() - 1.0;
		var2 = 2.0 * unit() - 1.0;
		rsquare = var1*var1 + var2*var2;
		} while(rsquare >= 1.0 || rsquare == 0.0);

		double val = -2.0 * log(rsquare) / rsquare;
		if(val > 0.0) factor = sqrt(val);
		else           factor = 0.0;	// should not happen, but might due to roundoff

		return (var1 * factor * _sd) + _mean;
	}

};

extern Rand2 RND;

#endif
