using namespace std;
#include "basicmodel.h"
#include <sstream>

void PrintPdat (const vector<pat>& pdat) {
    for (int i=0;i<pdat.size();i++) {
        cout << pdat[i].code << " " << pdat[i].time_s << " ";
        for (int j=0;j<pdat[i].locat.size();j++) {
            cout << pdat[i].locat[j].ward << " " << pdat[i].locat[j].date << " " << pdat[i].locat[j].prob << " ";
        }
        cout << "\n";
    }
}

void PrintVariantLoci (const vector<allele>& allvar) {
    cout << "Variant loci\n";
    for (int i=0;i<allvar.size();i++) {
        cout << allvar[i].loc << allvar[i].nuc << " ";
    }
    cout << "\n";
}

void PrintVariants (const vector<sparseseq>& variants, const vector<pat>& pdat) {
    cout << "Variants\n";
    for (int i=0;i<variants.size();i++) {
        cout << i << " " << pdat[i].code << " " << pdat[i].time_s << " ";
        for (int j=0;j<variants[i].locus.size();j++) {
            cout << variants[i].locus[j] << variants[i].allele[j] << " ";
        }
        cout << "\n";
    }
}

void PrintSeqDistances (const vector< vector<int> >& seqdists) {
	cout << "Sequence distances\n";
	for (int i=0;i<seqdists.size();i++) {
		for (int j=0;j<seqdists[i].size();j++) {
			cout << seqdists[i][j] << " ";
		}
		cout << "\n";
	}
}

void PrintSharedSingleVariants (vector<varcount>& vc_multi, vector<varcount>& vc_single) {
    cout << "Print vc_multi\n";
    for (int j=0;j<vc_multi.size();j++) {
        cout << vc_multi[j].locus << vc_multi[j].allele << " ";
        for (int k=0;k<vc_multi[j].individual.size();k++) {
            cout << vc_multi[j].individual[k] << " ";
        }
        cout << "\n";
    }
    
    cout << "Print vc_single\n";
    for (int j=0;j<vc_single.size();j++) {
        cout << vc_single[j].locus << vc_single[j].allele << " ";
            cout << vc_single[j].individual[0] << " ";
        cout << "\n";
    }
}

void ChainOutput (int c, const vector<tpairs>& trans_sets, vector<pat>& new_pdat) {
    cout << "Chain " << c << "\n";
    //Print chain
    for (int j=0;j<trans_sets[c].pair.size();j++) {
        cout << trans_sets[c].pair[j].from << "-" << trans_sets[c].pair[j].to << " ";
    }
    cout << "\n";

    //Times of data collection for these individuals
    for (int j=0;j<new_pdat.size();j++) {
        cout << "Time " << j << " " << new_pdat[j].time_seq << "\n";
    }
}

void OrderOutput(int c, int o, const vector<tpairs>& trans_sets) {
    cout << "Ordering " << o << "\n";
    //Print orders
    for (int k=0;k<trans_sets[c].orders[o].size();k++) {
        cout << trans_sets[c].orders[o][k] << " ";
    }
    cout << "\n";
    cout << "Size " << trans_sets[c].ordered_pair.size() << "\n";
    //Print ordered pairs
    for (int j=0;j<trans_sets[c].ordered_pair.size();j++) {
        cout << trans_sets[c].ordered_pair[j].from << "-" << trans_sets[c].ordered_pair[j].to << " ";
    }
    cout << "\n";

}

void PrintDescendentsSourcesBelow(const vector< vector<int> >& descendents, const vector< vector< vector<int> > >& below, const vector< vector<int> >& sources) {
    cout << "Descendants\n";
    for (int k=0;k<descendents.size();k++) {
        cout << k << " " << descendents[k][0] << " " << descendents[k][1] << "\n";
    }
      
    cout << "Sources\n";
    for (int k=0;k<sources.size();k++) {
          cout << k << " " << sources[k][0] << " " << sources[k][1] << "\n";
    }
    cout << "Below\n";
    for (int k=0;k<below.size();k++) {
        cout << k << " ";
        for (int l=0;l<below[k].size();l++) {
            for (int m=0;m<below[k][l].size();m++) {
                cout << below[k][l][m] << " ";
            }
            cout << "    ";
        }
        cout << "\n";
    }
}

void PrintTimings (int c, const vector<tpairs>& trans_sets, const vector< vector<ijlike> >& new_like_trans) {
    cout << "Timings\n";
    for (int k=0;k<trans_sets[c].ordered_pair.size();k++) {
        cout << new_like_trans[trans_sets[c].ordered_pair[k].from][trans_sets[c].ordered_pair[k].to].min << " " << new_like_trans[trans_sets[c].ordered_pair[k].from][trans_sets[c].ordered_pair[k].to].max << "\n";
    }
}

void PrintLikelihoods (int c, const vector<tpairs>& trans_sets, const vector< vector<ijlike> >& new_like_trans) {
    cout << "Likelihoods\n";
    for (int k=0;k<trans_sets[c].ordered_pair.size();k++) {
        for (int t=0;t<new_like_trans[trans_sets[c].ordered_pair[k].from][trans_sets[c].ordered_pair[k].to].noseq_likes.size();t++) {
            cout <<new_like_trans[trans_sets[c].ordered_pair[k].from][trans_sets[c].ordered_pair[k].to].noseq_likes[t] << " ";
        }
        cout << "\n";
    }
}

void PrintTimes (const vector<int>& times) {
    for (int k=0;k<times.size();k++) {
        cout << times[k] << " ";
    }
    cout << "\n";
}

void PrintTransSampleTimes (int c, const vector<int>& times, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets) {
    cout << "Transmission times\n";
    for (int k=0;k<times.size();k++) {//Action to be done here is outside of next loop
        cout << trans_sets[c].ordered_pair[k].from << "-" << trans_sets[c].ordered_pair[k].to << " T " << times[k] << "\n";
    }
    cout << "Sample times\n";
    for (int k=0;k<new_pdat.size();k++) {
        cout << "Time " << k << " " << new_pdat[k].time_seq << "\n";
    }
}

void PrintRelativeTTimes (int c, int o, vector<int>& times, const vector<tpairs>& trans_sets) {
    for (int i=0;i<trans_sets[c].orders[o].size();i++) {
        cout << trans_sets[c].ordered_pair[i].from << "-" << trans_sets[c].ordered_pair[i].to << " ";
    }
    cout << "T ";
    for (int i=0;i<times.size();i++) {
        cout << times[i]-times[0] << " ";
    }
}

void PrintAbsoluteTTimes (int c, int o, vector<int>& times, const vector<tpairs>& trans_sets) {
    for (int i=0;i<trans_sets[c].orders[o].size();i++) {
        cout << trans_sets[c].ordered_pair[i].from << "-" << trans_sets[c].ordered_pair[i].to << " ";
    }
    cout << "T ";
    for (int i=0;i<times.size();i++) {
        cout << times[i] << " ";
    }

}

void PrintTreespace (const vector<nbranch>& treespace) {
    for (int k=0;k<treespace.size();k++) {
        cout << k << " ";
        for (int l=0;l<treespace[k].individuals.size();l++) {
            cout <<treespace[k].individuals[l] << " ";
        }
        cout << " T " << treespace[k].time;
        cout << " M " << treespace[k].mutations;
        cout << "\n";
    }
}
