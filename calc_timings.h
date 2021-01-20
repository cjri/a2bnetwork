#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void CalculateTimingStats (run_params p);
void ReadTimeData (run_params p, vector< vector<int> >& timings, vector<double>& likes);
void CalculateRelativeTimings(vector<int>& max, vector< vector<int> >& timings);
void NormaliseLikelihoods (vector<double>& likes);
void CalculateTLikelihoods (vector<int>& max, vector< vector<double> >& data, vector<double>& likes, const vector< vector<int> >& timings);


