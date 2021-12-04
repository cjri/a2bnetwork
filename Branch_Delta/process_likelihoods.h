#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void ProcessLikelihoods (run_params p);
void ReadLikelihoodFile(run_params& p, vector<double>& likes, vector< vector<string> >& chains);
void ConstructTransmissions(const vector< vector<string> >& chains, vector< vector< vector<int> > >& transmissions);
void PrintTransmissions (vector< vector< vector<int> > >& transmissions);
void ProcessLikelihoodVals (vector<double>& likes);
void MakeConnections (const run_params& p, const vector<double>& likes, vector< vector<double> >& connects, const vector< vector< vector<int> > >& transmissions);
void PrintConnections (const vector< vector<double> >& connects);
