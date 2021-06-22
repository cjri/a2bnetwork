#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void CalculateTimingStats (run_params p);
void ReadTimeData (run_params p, int& min, int& max, vector<tdata>& timing_data);
void NormaliseLikelihoods (vector<tdata>& timing_data);
void SetupLikelihoods(run_params& p, int min, int max, vector< vector<double> >& likelihoods);
void AllocateEdgeNumbers (const vector<string>& all_edge, vector<tdata>& timing_data);
void CompileLikelihoodData (run_params& p, int min, vector<tdata>& timing_data, vector< vector<double> >& likelihoods);
void OutputTimeData (int min, const vector<string>& all_edge, const vector< vector<double> >& likelihoods);



