#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void ConvertLikelihoodsAlt (run_params p);
void CollectTransmissionData (const vector<int>& in_init, const vector< vector< vector<int> > >& transmissions, vector< vector<int> >& trans_data);
void WriteLikelihoodsAlt (const vector< vector<int> >& trans_data, const vector<double>& likes);


