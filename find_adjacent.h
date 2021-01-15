#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>


void ListAdjacentK (run_params p, int k);
void RemoveDuplicatesTransIn(vector< vector<int> >& trans_in);
void AnalyseAdjacent (run_params p, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, const vector< vector<sparseseq> >& variants_sets);
void ReadBestLikelihood (run_params p, vector<double>& likelihoods, vector< vector<int> >& trans_in);
void FindDownstreamNodes (const vector< vector<int> >& trans_in, vector< vector<int> >& downstream);
void PrintDownstream (const vector< vector<int> >& downstream);
void FindInitialAdjacent (run_params p, int start, const vector< vector<int> >& downstream, vector< vector<int> >& trans_in);
void PrintAdjacentNetworks (const vector< vector<int> >& trans_in);
void RemoveIfLikelihoodKnown (run_params p, vector< vector<int> >& trans_in, vector<double>& like_before, vector< vector<int> >& c_before);
void FindTransGroup (const vector< vector<int> >& trans_in, vector< vector<int> >& trans_group);
void PrintTransGroup (const vector< vector<int> >& trans_group);
void ReadOrders (const vector< vector<int> >& trans_in, const vector< vector<int> >& trans_group, vector<tpairs>& trans_sets, vector<int>& new_subset);
void PrintTransSets (vector<tpairs>& trans_sets);
void PrintTransSetsFile (vector<tpairs>& trans_sets);
void PrintPrevious (ofstream& cp_file, const vector<double>& like_before, const vector< vector<int> >& c_before);
void CalculateLikelihoods (run_params p, ofstream& cp_file, vector<tpairs>& trans_sets, const vector<varcount>& vc_multi, const vector<varcount>& vc_single, const vector<pat>& new_pdat, vector< vector<ijlike> >& new_like_trans, const vector< vector<int> >& orders);
void CalculateLikelihoodsFromList (run_params p, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, const vector< vector<sparseseq> >& variants_sets);

