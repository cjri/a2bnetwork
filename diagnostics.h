#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void PrintPdat (const vector<pat>& pdat);
void PrintVariantLoci (const vector<allele>& allvar);
void PrintVariants (const vector<sparseseq>& variants, const vector<pat>& pdat);
void PrintSeqDistances (const vector< vector<int> >& seqdists);
void PrintSharedSingleVariants (vector<varcount>& vc_multi, vector<varcount>& vc_single);
void ChainOutput (int c, const vector<tpairs>& trans_sets, vector<pat>& new_pdat);
void OrderOutput(int c, int o, const vector<tpairs>& trans_sets);
void PrintDescendentsSourcesBelow(const vector< vector<int> >& descendents, const vector< vector< vector<int> > >& below, const vector< vector<int> >& sources);
void PrintTimings (int c, const vector<tpairs>& trans_sets, const vector< vector<ijlike> >& new_like_trans);
void PrintLikelihoods (int c, const vector<tpairs>& trans_sets, const vector< vector<ijlike> >& new_like_trans);
void PrintTimes (const vector<int>& times);
void PrintTransSampleTimes (int c, const vector<int>& times, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets);
void PrintTreespace (const vector<nbranch>& treespace);
