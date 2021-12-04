#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

//void MatchSequencetoPatient (vector<pat>& pdat);

void CollectSubsetData (const vector< vector<int> >& subsets, const vector<pat>& pdat, const vector< vector<ijlike> >& like_trans, const vector<sparseseq>& variants, vector< vector<pat> >& pdat_sets, vector< vector< vector<ijlike> > >& like_trans_sets, vector< vector<sparseseq> >& variants_sets, vector< vector<int> >& new_subsets);
void MakeNewSubsets (const vector< vector<int> >& subsets, vector< vector<int> >& new_subsets);
void SplitPDatSubsets (const vector< vector<int> >& subsets, const vector<pat>& pdat, vector< vector<pat> >& pdat_sets);
void SplitLikeTransSubsets (const vector< vector<int> >& subsets, const vector< vector<ijlike> >& like_trans, vector< vector< vector<ijlike> > >& like_trans_sets);
void SplitVariantsSubsets (const vector< vector<int> >& subsets, const vector<sparseseq>& variants, vector< vector<sparseseq> >& variants_sets);

