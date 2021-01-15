#include "basicmodel.h"
#include "process_likelihoods.h"
#include "convert_likelihoods_alt.h"

void ConvertLikelihoodsAlt (run_params p) {
    vector< vector<string> > chains;
    vector<double> likes;
    ReadLikelihoodFile(p,likes,chains);
    vector< vector< vector<int> > > transmissions;
    ConstructTransmissions(chains,transmissions);
    vector<int> in_init;
    for (int j=0;j<p.n_edges+1;j++)  {
        in_init.push_back(-1);
    }
    vector< vector<int> > trans_data;
    CollectTransmissionData(in_init,transmissions,trans_data);
    WriteLikelihoodsAlt(trans_data,likes);
}

void CollectTransmissionData (const vector<int>& in_init, const vector< vector< vector<int> > >& transmissions, vector< vector<int> >& trans_data) {
    for (int i=0;i<transmissions.size();i++) {
        vector<int> temp=in_init;
        for (int j=0;j<transmissions[i].size();j++) {
            temp[transmissions[i][j][1]]=transmissions[i][j][0];
        }
        trans_data.push_back(temp);
    }
}

void WriteLikelihoodsAlt (const vector< vector<int> >& trans_data, const vector<double>& likes) {
    ofstream l_file;
    l_file.open("Likelihoods_alt.out");
    for (int i=0;i<trans_data.size();i++) {
        for (int j=0;j<trans_data[i].size();j++) {
            l_file << trans_data[i][j] << " ";
        }
        l_file << likes[i] << "\n";
    }
}
