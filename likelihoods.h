#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

/*struct placetime {
    int date;
    string ward;
};*/


void CalculateTDLikelihoods (run_params p, const vector<pat>& pdat, const vector< vector<int> >& seqdists, const vector< vector<tpair> >& seqdists_c, vector< vector<ijlike> >& like_trans);
bool comparetprob(tprob v1, tprob v2);
void RemoveDuplicateTimes (vector<tprob>& contact_times_probs);
void FillTimes (vector<tprob>& contact_times_probs);
double LikelihoodFromItoJTimeK (int i, int j, int k, const run_params& p, const vector<tprob>& contact_times_probs, const vector< vector<int> >& seqdists, const vector< vector<tpair> >& seqdists_c, const vector<pat>& pdat);
double NoSeqLikelihoodFromItoJTimeK (int i, int j, int k, const run_params& p, const vector<tprob>& contact_times_probs, const vector<pat>& pdat);
void TrimLikelihoods (ijlike& lij);
void GetSubsetsIJ (run_params p, const vector< vector<ijlike> >& like_trans, const vector<pat>& pdat, vector< vector<int> >& subsets);
void FindLikelihoodSubsetsIJ(run_params p, const vector< vector<ijlike> >& like_trans, vector< vector<int> >& subsets);


void FindNetworkLikelihood (run_params p, vector< vector<int> >& new_subsets, vector< vector<pat> >& pdat_sets, vector< vector< vector<ijlike> > >& like_trans_sets, vector< vector<sparseseq> >& variants_sets);
void SetupTreeLikelihoodCalculation(run_params p, int set, const vector<int>& new_subset, vector<tpairs>& trans_sets, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, const vector< vector<sparseseq> >& variants_sets, vector< vector<int> >& orders, vector<pat>& new_pdat, vector< vector<ijlike> >& new_like_trans, vector<sparseseq>& new_variants, vector<varcount>& vc_multi, vector<varcount>& vc_single);
void EditSharedVariants(const vector<tpairs>& trans_sets, vector<varcount>& vc_multi);
void MakeOrderFromIndex (int c, const vector< vector<int> >& orders, vector<tpairs>& trans_sets);
void CalculateOrderLike (run_params& p, int c, int o, double& chain_L, vector<tpairs>& trans_sets, const vector<varcount>& vc_multi, const vector<varcount>& vc_single, const vector<pat>& new_pdat, vector< vector<ijlike> >& new_like_trans);
void MakeOrderedPairs (int c, int o, vector<tpairs>& trans_sets);
void MakeDescendents(int c, const vector<tpairs>& trans_sets, vector< vector<int> >& descendents);
void MakeBelow (const vector< vector<int> >& descendents, vector< vector< vector<int> > >& below);
void MakeSources (int c, const vector<tpairs>& trans_sets, vector< vector<int> >& sources);
void SetUpTimeVectors(int c, int o, const vector<tpairs>& trans_sets, const vector< vector<ijlike> >& new_like_trans, vector<int>& times, vector<int>& min_times, vector<int>& max_times, vector<int>& equiv_times);
void CheckTimeOrder(int c, int o, const vector<tpairs>& trans_sets, vector<int>& times);
void CheckMaxTimes (int c, int o, const vector<tpairs>& trans_sets, vector<int>& max_times);
int CheckFeasible (const vector<int>& times, const vector<int>& max_times);
void GenerateLikelihoodAllTimings (int c, int o, run_params& p, double& order_L, vector<int>& times, const vector<int>& min_times, const vector<int>& max_times, vector<int>& equiv_times, const vector<varcount>& vc_single, const vector<varcount>& vc_multi, const vector< vector<int> >& sources, const vector< vector<int> >& descendents, const vector< vector< vector<int> > >& below, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets, const vector< vector<ijlike> >& new_like_trans);
void CreateTimeTree(int c, const vector<int>& times, const vector< vector<int> >& sources, const vector< vector<int> >& descendents, const vector< vector< vector<int> > >& below, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets, vector<nbranch>& treespace);
void BranchTwoIndividuals(int c, int index, const vector<int>& times, const vector< vector<int> >& sources, const vector< vector<int> >& descendents, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets, vector<nbranch>& treespace);
void BranchIndividualNode(int c, int index, const vector<int>& times, const vector< vector<int> >& sources, const vector< vector<int> >& descendents, const vector< vector< vector<int> > >& below, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets, vector<nbranch>& treespace);
void BranchNodeIndividual (int c, int index, const vector<int>& times, const vector< vector<int> >& sources, const vector< vector<int> >& descendents, const vector< vector< vector<int> > >& below, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets, vector<nbranch>& treespace);
void BranchTwoNodes (int c, int index, const vector<int>& times, const vector< vector<int> >& sources, const vector< vector<int> >& descendents, const vector< vector< vector<int> > >& below, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets, vector<nbranch>& treespace);
int FindPrevious (int c, int source, int index, const vector<int>& times, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets);
void AssignMutationsToBranches (const vector<varcount>& vc_multi, const vector<varcount>& vc_single, vector<nbranch>& treespace);
void FindEquivalentTimes(const vector<int>& times, const vector<int>& min_times, vector<int>& equiv_times);
void CalculateTreeLikelihood (int c, run_params& p, double& order_L, const vector<int>& times, const vector<int>& equiv_times, const vector<nbranch>& treespace, const vector<tpairs>& trans_sets, const vector< vector<ijlike> >& new_like_trans);
void FindNextTimeSet (int c, int o, const run_params& p, int& done, int& up, int& fin, const vector<int>& min_times, const vector<int>& max_times, const vector<int>& orig_times, const vector<tpairs>& trans_sets, const vector< vector<ijlike> >& new_like_trans, vector<int>& times);








void Thresholds (run_params p, double L);

