#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "basicmodel.h"
#include "utilities.h"
#include "distributions.h"
#include "likelihoods.h"
#include "io.h"

void CollectSubsetData (const vector< vector<int> >& subsets, const vector<pat>& pdat, const vector< vector<ijlike> >& like_trans, const vector<sparseseq>& variants, vector< vector<pat> >& pdat_sets, vector< vector< vector<ijlike> > >& like_trans_sets, vector< vector<sparseseq> >& variants_sets, vector< vector<int> >& new_subsets) {
    MakeNewSubsets(subsets,new_subsets);
	SplitPDatSubsets(subsets,pdat,pdat_sets);
	SplitLikeTransSubsets (subsets,like_trans,like_trans_sets);
	SplitVariantsSubsets (subsets,variants,variants_sets);
	cout << "Number of subsets to analyse " << pdat_sets.size() << "\n";
}

void MakeNewSubsets (const vector< vector<int> >& subsets, vector< vector<int> >& new_subsets) {
    for (int i=0;i<subsets.size();i++) {
        vector<int> ss;
        if (subsets[i].size()>1) {
            ss=subsets[i];
            new_subsets.push_back(ss);
        }
    }
}

void SplitPDatSubsets (const vector< vector<int> >& subsets, const vector<pat>& pdat, vector< vector<pat> >& pdat_sets) {
	for (int i=0;i<subsets.size();i++) {
		vector<pat> pd;
		if (subsets[i].size()>1) {
			for (int j=0;j<subsets[i].size();j++) {
				pd.push_back(pdat[subsets[i][j]]);
			}
			pdat_sets.push_back(pd);
		}
	}
}

void SplitLikeTransSubsets (const vector< vector<int> >& subsets, const vector< vector<ijlike> >& like_trans, vector< vector< vector<ijlike> > >& like_trans_sets) {
	for (int i=0;i<subsets.size();i++) {
		vector< vector<ijlike> > likes;
		if (subsets[i].size()>1) {
			for (int j=0;j<subsets[i].size();j++) {
				vector<ijlike> l;
				for (int k=0;k<subsets[i].size();k++) {
					l.push_back(like_trans[subsets[i][j]][subsets[i][k]]);
				}
				likes.push_back(l);
			}
			like_trans_sets.push_back(likes);
		}
	}
}

void SplitVariantsSubsets (const vector< vector<int> >& subsets, const vector<sparseseq>& variants, vector< vector<sparseseq> >& variants_sets) {
	for (int i=0;i<subsets.size();i++) {
		vector<sparseseq> vars;
		if (subsets[i].size()>1) {
			for (int j=0;j<subsets[i].size();j++) {
				vars.push_back(variants[subsets[i][j]]);
			}
            variants_sets.push_back(vars);
		}
	}
}

