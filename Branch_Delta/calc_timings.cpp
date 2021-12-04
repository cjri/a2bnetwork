#include "basicmodel.h"
#include "calc_timings.h"

void CalculateTimingStats (run_params p) {
    vector<tdata> timing_data;
    int min=1000;
    int max=0;
    ReadTimeData (p,min,max,timing_data);
    NormaliseLikelihoods(timing_data);

    //Setup likelihoods;
    vector< vector<double> > likelihoods;
    SetupLikelihoods(p,min,max,likelihoods);

    vector<string> all_edge=timing_data[timing_data.size()-1].edges;
    //Allocate edge numbers
    cout << "Allocating edge numbers\n";
    AllocateEdgeNumbers (all_edge,timing_data);
    
    //Fill likelihood table
    CompileLikelihoodData (p,min,timing_data,likelihoods);

    OutputTimeData (min,all_edge,likelihoods);
    
}

void CompileLikelihoodData (run_params& p, int min, vector<tdata>& timing_data, vector< vector<double> >& likelihoods) {
    cout << "Get likelihoods\n";
    if (p.absolute==0) {//Have relative timings
        for (int i=0;i<timing_data.size();i++) {
            for (int j=0;j<timing_data[i].times.size();j++) {
                likelihoods[timing_data[i].edge_number[j]][timing_data[i].times[j]]=likelihoods[timing_data[i].edge_number[j]][timing_data[i].times[j]]+timing_data[i].L;
            }
        }
    } else {
        for (int i=0;i<timing_data.size();i++) {
            for (int j=0;j<timing_data[i].times.size();j++) {
                likelihoods[timing_data[i].edge_number[j]][timing_data[i].times[j]-min]=likelihoods[timing_data[i].edge_number[j]][timing_data[i].times[j]-min]+timing_data[i].L;
            }
        }
    }

}

void OutputTimeData (int min, const vector<string>& all_edge, const vector< vector<double> >& likelihoods) {
    cout << "X ";
    for (int i=0;i<all_edge.size();i++) {
        cout << all_edge[i] << " ";
    }
    cout << "\n";
    
    for (int i=0;i<likelihoods[0].size();i++) {
        cout << min+i << " ";
        for (int j=0;j<likelihoods.size();j++) {
            cout << likelihoods[j][i] << " ";
        }
        cout << "\n";
    }

}

void ReadTimeData (run_params p, int& min, int& max, vector<tdata>& timing_data) {
    cout << "Reading data " << p.time_file << "\n";
    ifstream time_file;
    time_file.open(p.time_file.c_str());
    int n;
    double x;
    string s;
    for (int i=0;i<500000000;i++) {
        tdata t;
        for (int j=0;j<p.n_edges;j++) {
            if (!(time_file >> s)) break;
            t.edges.push_back(s);
          //  cout << s << " ";
        }
        if (!(time_file >> s)) break;
        for (int j=0;j<p.n_edges;j++) {
            if (!(time_file >> n)) break;
            t.times.push_back(n);
           // cout << n << " ";
            if (n<min) {
                min=n;
            }
            if (n>max) {
                max=n;
            }
        }
       // cout << "\n";
        if (!(time_file >> x)) break;
        t.L=x;
        timing_data.push_back(t);
    }
}

void NormaliseLikelihoods (vector<tdata>& timing_data) {
    cout << "Calculate normalised likelihoods " << timing_data.size() << "\n";
    double tot=0;
    for (int i=0;i<timing_data.size();i++) {
        timing_data[i].L=exp(timing_data[i].L);
        tot=tot+timing_data[i].L;
    }
    for (int i=0;i<timing_data.size();i++) {
        timing_data[i].L=timing_data[i].L/tot;
    }
}

void SetupLikelihoods(run_params& p, int min, int max, vector< vector<double> >& likelihoods) {
    cout << "Setup likelihoods " << min << " " << max << "\n";
    vector<double> l_store;
    for (int i=min;i<=max;i++) {
        l_store.push_back(0);
    }
    for (int i=0;i<p.n_edges;i++) {
        likelihoods.push_back(l_store);
    }

}

void AllocateEdgeNumbers (const vector<string>& all_edge, vector<tdata>& timing_data) {
    for (int i=0;i<timing_data.size();i++) {
        for (int j=0;j<timing_data[i].edges.size();j++) {
            for (int k=0;k<all_edge.size();k++) {
                if (timing_data[i].edges[j]==all_edge[k]) {
                    timing_data[i].edge_number.push_back(k);
                    break;
                }
            }
        }
    }
}

