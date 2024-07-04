#include "optfem.h"

/* ���ȯ���ץ������ */
/* ¾�Υץ��������Ѥ������ seed, basic_seed */
/* �򥹥��ƥ��å�Ū���ѿ��ˤ����Ѥ��� */

/*
Caution!
Before using it, initialization of random number sequence 
is necessary. It is done by using randomize().
time function is used for the initialization.

randomize() : initialization of random number sequence
randperc() : generate random number between 0 and 1 
rand_normal(mu, sigma) : generata random number with Gaussian distribution 
with the average mu and the variance sigma
*/


double oldrand[55];   /* Array of 55 random numbers */
int jrand;				 /* current random number */

/* Create next batch of 55 random numbers */
void advance_random()
{
	int j1;
	double new_random;

	for(j1 = 0; j1 < 24; j1++)
	{
		new_random = oldrand[j1] - oldrand[j1+31];
		if (new_random < 0.0) new_random = new_random + 1.0;
		oldrand[j1] = new_random;
	}
	for(j1 = 24; j1 < 55; j1++)
	{
		new_random = oldrand [j1] - oldrand [j1-24];
		if (new_random < 0.0) new_random = new_random + 1.0;
		oldrand[j1] = new_random;
	}
}

/* Flip a biased coin - true if heads */
int flip(double prob)
{
	if (randomperc() <= prob)
		return(1);
	else
		return(0);
}

/* Get seed number for random and start it up */
void randomize()
{
	int j1;
	double seed;

	seed = (double)rand()/(double)(RAND_MAX);

	for(j1=0; j1<55; j1++) oldrand[j1] = 0.0;
	jrand=0;

	warmup_random(seed);
}

/* Fetch a single random number between 0.0 and 1.0 -  */
/* Subtractive Method . See Knuth, D. (1969), v. 2 for */
/* details.Name changed from random() to avoid library */
/* conflicts on some machines						  */
double randomperc()
{
	jrand++;
	if (jrand >= 55)
	{
		jrand = 1;
		advance_random();
	}
	return((double) oldrand[jrand]);
}

/* Pick a random integer between low and high */
int rnd(int low, int high)
{
	int i;
	if (low >= high) {
		i = low;
	}
	else
	{
		i = ((int)randomperc() * (high - low + 1)) + low;
		// *added the above code is modified to the below by muratsubaki
		// i = (int)(randomperc() * (high - low + 1)) + low; ??
		if (i > high) i = high;
	}
	return(i);
}

/* real random number between specified limits */
double rndreal(double lo, double hi)
{
	return((randomperc() * (hi - lo)) + lo);
}

/* Get random off and running */
void warmup_random(double random_seed)
{
	int j1, ii;
	double new_random, prev_random;

	oldrand[54] = random_seed;
	new_random = 0.000000001;
	prev_random = random_seed;
	for(j1 = 1 ; j1 <= 54; j1++)
	{
		ii = (21*j1)%54;
		oldrand[ii] = new_random;
		new_random = prev_random-new_random;
		if (new_random < 0.0) new_random = new_random + 1.0;
		prev_random = oldrand[ii];
	}

	advance_random();
	advance_random();
	advance_random();

	jrand = 0;
}

/* Normal distribution */
double rand_normal( double mu, double sigma ) {
	double z=sqrt( -2.0*log(randomperc()) ) * sin( 2.0*M_PI*randomperc() );
	return mu + sigma*z;
}
