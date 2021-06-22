#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

struct stat {
	int s1;
	int s2;
	int t;
	int d1;
	int d2;
	int h;
	int b;
	int h1;
	int h2;
	double L;
};

struct cstat {
	int s1;
	int s2;
	int t;
	int d1;
	int d2;
	int h;
	int b;
	int h1;
	int h2;
	double L;
};

void CalculateThresholdsNoSeq (run_params p, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc);
vector<double> CalculateThresholdsNSSpecific (run_params p, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc);
void CalculateThresholdsFullExplicit (run_params p, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc, gsl_rng *rgen);
vector<double> CalculateThresholdsSpecificDExplicit (run_params p, int d1, int d2d, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc, gsl_rng *rgen);


