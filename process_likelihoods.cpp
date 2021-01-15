using namespace std;
#include <sstream>
#include "basicmodel.h"
#include "process_likelihoods.h"

void ProcessLikelihoods (run_params p) {
    //Read in data
    vector< vector<string> > chains;
    vector<double> likes;
    ReadLikelihoodFile(p,likes,chains);
    vector< vector< vector<int> > > transmissions;
    ConstructTransmissions(chains,transmissions);
    if (p.diagnostic==1) {
        //PrintTransmissions(transmissions);
    }
    ProcessLikelihoodVals (likes);
    vector< vector<double> > connects;
    MakeConnections (p,likes,connects,transmissions);
    PrintConnections(connects);
}

void ReadLikelihoodFile(run_params& p, vector<double>& likes, vector< vector<string> >& chains) {
    ifstream cp_file;
    string name=p.likelihood_file;
    cp_file.open(name.c_str());
    string s;
    int n;
    double x;
    for (int i=0;i<1000000;i++) {
        vector<string> tss;
        if (!(cp_file >> s)) break;
        if (!(cp_file >> n)) break;
        for (int j=0;j<p.n_edges;j++) {
            if (!(cp_file >> s)) break;
            tss.push_back(s);
        }
        if (tss.size()>0) {
            chains.push_back(tss);
        }
        if (!(cp_file >> s)) break;
        if (!(cp_file >> n)) break;
        if (!(cp_file >> x)) break;
        likes.push_back(x);
    }
    //cout << "Size " << chains.size() << "\n";
}

void ConstructTransmissions(const vector< vector<string> >& chains, vector< vector< vector<int> > >& transmissions) {
    for (int j=0;j<chains.size();j++) {
        vector< vector<int> > tt;
        for (int k=0;k<chains[j].size();k++) {
            vector<int> t;
            stringstream sss(chains[j][k]);
            while (sss.good()) {
                string sr;
                getline(sss,sr,'-');
                t.push_back(atoi(sr.c_str()));
            }
            tt.push_back(t);
        }
        transmissions.push_back(tt);
    }
    //cout << transmissions.size() << "\n";
}

void PrintTransmissions (vector< vector< vector<int> > >& transmissions) {
    for (int j=0;j<transmissions.size();j++) {
        cout << "Set " << j << "\n";
        for (int k=0;k<transmissions[j].size();k++) {
            cout << transmissions[j][k][0] << " " << transmissions[j][k][1] << "    ";
        }
        cout << "\n";
    }
}

void ProcessLikelihoodVals (vector<double>& likes) {
    double tot=0;
    for (int j=0;j<likes.size();j++) {
        likes[j]=exp(likes[j]);
        tot=tot+likes[j];
    }
    for (int j=0;j<likes.size();j++) {
        likes[j]=likes[j]/tot;
    }
}

void MakeConnections (const run_params& p, const vector<double>& likes, vector< vector<double> >& connects, const vector< vector< vector<int> > >& transmissions) {
    vector<double> ed;
    for (int j=0;j<p.n_edges+1;j++) {
        ed.push_back(0);
    }
    for (int j=0;j<p.n_edges+1;j++) {
        connects.push_back(ed);
    }
    for (int j=0;j<transmissions.size();j++) {
        for (int k=0;k<transmissions[j].size();k++) {
            //cout << transmissions[j][k][0] << " " << transmissions[j][k][1] << "\n";
            connects[transmissions[j][k][0]][transmissions[j][k][1]]=connects[transmissions[j][k][0]][transmissions[j][k][1]]+likes[j];
        }
    }
    //cout << connects.size() << "\n";
}

void PrintConnections (const vector< vector<double> >& connects) {
    ofstream edge_file;
    edge_file.open("Edge_occupancies.out");
    for (int j=0;j<connects.size();j++) {
        for (int k=0;k<connects[j].size();k++) {
            edge_file << connects[j][k] << "  ";
        }
        edge_file << "\n";
    }
}
