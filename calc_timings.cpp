#include "basicmodel.h"
#include "calc_timings.h"

void CalculateTimingStats (run_params p) {
    vector< vector<int> > timings;
    vector<double> likes;
    ReadTimeData (p,timings,likes);
    vector<int> max=timings[0];
    CalculateRelativeTimings(max,timings);
    NormaliseLikelihoods(likes);
    vector< vector<double> > data;
    CalculateTLikelihoods (max,data,likes,timings);

    //Output data
    ofstream rel_file;
    rel_file.open("Relative_times.out");
    for (int i=0;i<data.size();i++) {
        for (int j=0;j<data[i].size();j++) {
            rel_file << data[i][j] << " ";
        }
        rel_file << "\n";
    }
}

void ReadTimeData (run_params p, vector< vector<int> >& timings, vector<double>& likes) {
    cout << "Reading data\n";
    ifstream time_file;
    time_file.open(p.time_file.c_str());
    int n;
    double x;
    vector<int> init;
    for (int i=0;i<=p.n_edges;i++) {
        init.push_back(-1);
    }
    vector<int> tmp;
    int index=-1;
    for (int i=0;i<500000000;i++) {
        tmp=init;
        vector<int> tm;
        for (int j=0;j<p.n_edges;j++) {
            if (!(time_file >> n)) break;
            index=n;
            if (!(time_file >> n)) break;
            tmp[index]=n;
        }
        if (!(time_file >> x)) break;
        timings.push_back(tmp);
        likes.push_back(x);
    }
}

void CalculateRelativeTimings(vector<int>& max, vector< vector<int> >& timings) {
    cout << "Calculate relative times\n";
    for (int i=0;i<max.size();i++) {
        max[i]=0;
    }
    for (int i=0;i<timings.size();i++) {
        int zero=1000;
        for (int j=0;j<timings[i].size();j++) {
            if (timings[i][j]>=0&&timings[i][j]<zero) {
                zero=timings[i][j];
            }
        }
        for (int j=0;j<timings[i].size();j++) {
            if (timings[i][j]>=zero) {
                timings[i][j]=timings[i][j]-zero;
                if (timings[i][j]>max[j]) {
                    max[j]=timings[i][j];
                }
            }
        }
    }
}

void NormaliseLikelihoods (vector<double>& likes) {
    cout << "Calculate real likelihoods\n";
    double tot=0;
    for (int i=0;i<likes.size();i++) {
        likes[i]=exp(likes[i]);
        tot=tot+likes[i];
    }
    for (int i=0;i<likes.size();i++) {
        likes[i]=likes[i]/tot;
    }
}

void CalculateTLikelihoods (vector<int>& max, vector< vector<double> >& data, vector<double>& likes, const vector< vector<int> >& timings) {
    for (int i=0;i<max.size();i++) {
        vector<double> d;
        for (int j=0;j<max[i];j++) {
            d.push_back(0);
        }
        d.push_back(0);
        data.push_back(d);
    }
    for (int i=0;i<timings.size();i++) {
        for (int j=0;j<timings[i].size();j++) {
            if (timings[i][j]>=0) {
                data[j][timings[i][j]]=data[j][timings[i][j]]+likes[i];
            }
        }
    }
}
