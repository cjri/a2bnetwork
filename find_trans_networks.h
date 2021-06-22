#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void FindPlausibleTransmissionNetworks (run_params& p, const vector< vector<int> >& new_subsets, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, vector< vector<sparseseq> >& variants_sets);
void FPTN_Set (run_params& p, int set, const vector< vector<int> >& new_subsets, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, vector< vector<sparseseq> >& variants_sets, ofstream& check_file);
void ConstructIndivSets (const int remove, const vector<int>& indivs, vector< vector<int> >& indiv_sets, vector< vector<int> >& remainders);
void SetupParamsSubsetAnalysis (int set, const vector<int>& new_subset, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, const vector< vector<sparseseq> >& variants_sets, vector<pat>& new_pdat, vector< vector<ijlike> >& new_like_trans, vector<sparseseq>& new_variants);
void PrintVariantInformation (const vector<pat>& new_pdat, const vector<sparseseq>& new_variants);
void PrintNLT (const vector< vector<ijlike> >& new_like_trans);
void FindSharedVariantsNew(vector<varcount>& vc_multi, vector<varcount>& vc_single, const vector<sparseseq>& new_variants);
void FindTransSets(run_params& p, const vector<int>& indivs, const vector< vector<ijlike> >& new_like_trans, const vector<varcount>& vc_multi, vector<tpairs>& trans_sets);
void FindPairSets (run_params& p, const vector<int>& indivs, const vector< vector<ijlike> >& like_trans, const vector<varcount>& vc_multi, vector<tpairs>& trans_sets);
void CheckLikelihood(int& check_unlike, run_params& p, const vector<int>& index, const vector< vector<tpair> >& allsets, const  vector< vector<ijlike> >& like_trans);
void CheckSpan (int& check_span, const vector<int>& index, const vector< vector<tpair> >& allsets);
void CheckCycle (int& check_cycle, const vector<int>& index, const vector< vector<tpair> >& allsets);
int TransFindCycle (vector<tpair> tpset);
void CheckSharedVariants(int& check_seq, const vector<int>& index, const vector< vector<tpair> >& allsets, const vector<varcount>& vc_multi);
void OrderBySharedVariants(const vector<int>& index, const vector< vector<tpair> >& allsets, const vector<varcount>& vc_multi, tpairs& ts);
void ImportTimings(const vector< vector<ijlike> >& like_trans, vector<tpairs>& trans_sets);
void FindOrderingConstraints(vector<tpairs>& trans_sets);
void CheckOrderingContradictions (vector<tpairs>& trans_sets);
void ProcessRemainder(int set, int iset, run_params& p, ofstream& check_file, const vector< vector<int> >& remainders, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, const vector< vector<sparseseq> >& variants_sets);
