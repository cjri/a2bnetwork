#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>


void ProcessVariantInformation (run_params p, const vector<string>& seqs, vector<sparseseq>& variants, vector<pat>& pdat);
void FindConsensus (string& consensus, const vector<string>& seqs);
void FindVariants (vector<sparseseq>& variants, string& consensus, vector<pat>& pdat);
bool compare_allele(allele v1, allele v2);
void ListAllVariantPositions (const vector<sparseseq>& variants, vector<allele>& allvar);
void FindAmbiguousVarPositions (const vector<allele>& allvar, vector<pat>& pdat, vector<int>& nloci, vector<int>& nloc_count);
void FindAmbiguousVarFreqs (const string all_consensus, const vector<allele>& allvar, const vector<sparseseq>& variants, vector<pat>& pdat, const vector<int>& nloci, vector<ipair>& freqs);
void ProcessAmbiguousFixedVariants(vector<ipair>& freqs, vector<int>& nloc_count, vector<sparseseq>& variants, vector<pat>& pdat);
void RemoveIndividualsMultipleN (const run_params p, vector<int>& nloc_count, vector<pat>& pdat, vector<sparseseq>& variants);
void FindPairwiseDistances (run_params p, vector< vector<int> >& seqdists, vector<sparseseq>& variants, vector<pat>& pdat);
void FromConsensusDistances (const vector<sparseseq>& variants, vector< vector<tpair> >& seqdists_c);





void ModifyVariantsByBin(int b, const vector< vector<int> >& bin, double& lLvar, const vector<pat>& pdat, vector<sparseseq>& variants);
void FindSharedVariants(int set, vector<varcount>& vc_multi, const vector< vector<sparseseq> >& variants_sets);
void FindSharedVariantsNew(vector<varcount>& vc_multi, vector<varcount>& vc_single, const vector<sparseseq>& new_variants);
void EditSharedVariants(const vector<tpairs>& trans_sets, vector<varcount>& vc_multi);
void GetSharedInOut (int c, const vector<tpairs>& trans_sets, const vector<varcount>& vc_multi, vector< vector<int> >& inout);
void GetNonSharedInOut (int c, const vector<tpairs>& trans_sets, const vector<varcount>& vc_single, vector< vector<int> >& inout_s);


