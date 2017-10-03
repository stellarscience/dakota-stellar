#ifndef RANDOM_H
#define RANDOM_H


#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <sys/time.h>
#include <iostream>


class Random
	{
	public:
		static int timeSeed();
		static void initRandom(int seed);

		static double uniformDeviate(double a, double b);

		static Array<int> randomIVector(int length);
	private:
		static int globSeed[4];
		static bool globSeedReady;
	};

#endif

