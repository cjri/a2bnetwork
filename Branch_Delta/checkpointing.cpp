using namespace std;
#include "basicmodel.h"
#include "checkpointing.h"
#include <sstream>

void WriteCheckpoint1 (int set, int iset, vector<int> new_subset, vector<tpairs>& trans_sets) {
    ofstream cp_file;
    ostringstream convert;
    ostringstream convert2;
    convert << set;
    string temp=convert.str();
    convert2 << iset;
    string temp2=convert2.str();
    string name="Check1_"+temp+"_"+temp2+".out";
    cp_file.open(name.c_str());
    cp_file << new_subset.size() << " ";
    for (int i=0;i<new_subset.size();i++) {
        cp_file << new_subset[i] << " ";
    }
    cp_file << "\n";
    for (int i=0;i<trans_sets.size();i++) {
        cp_file << trans_sets[i].root << " " << trans_sets[i].pair.size() << " ";
        for (int j=0;j<trans_sets[i].pair.size();j++) {
            cp_file <<trans_sets[i].pair[j].from << " " << trans_sets[i].pair[j].to << " " << trans_sets[i].min_time[j] << " " << trans_sets[i].max_time[j] << " ";
        }
        cp_file << trans_sets[i].constraints.size() << " ";
        for (int j=0;j<trans_sets[i].constraints.size();j++) {
            cp_file << trans_sets[i].constraints[j][0] << " " << trans_sets[i].constraints[j][1] << " ";
        }
        cp_file << trans_sets[i].const_contradiction << "\n";
    }
}

void WriteCheckpoint_Rem1 (int set, int iset, vector<int> new_subset, vector<tpairs>& trans_sets) {
    ofstream cp_file;
    ostringstream convert;
    ostringstream convert2;
    convert << set;
    string temp=convert.str();
    convert2 << iset;
    string temp2=convert2.str();
    string name="Check1_"+temp+"_"+temp2+"_Rem.out";
    cp_file.open(name.c_str());
    cp_file << new_subset.size() << " ";
    for (int i=0;i<new_subset.size();i++) {
        cp_file << new_subset[i] << " ";
    }
    cp_file << "\n";
    for (int i=0;i<trans_sets.size();i++) {
        cp_file << trans_sets[i].root << " " << trans_sets[i].pair.size() << " ";
        for (int j=0;j<trans_sets[i].pair.size();j++) {
            cp_file <<trans_sets[i].pair[j].from << " " << trans_sets[i].pair[j].to << " " << trans_sets[i].min_time[j] << " " << trans_sets[i].max_time[j] << " ";
        }
        cp_file << trans_sets[i].constraints.size() << " ";
        for (int j=0;j<trans_sets[i].constraints.size();j++) {
            cp_file << trans_sets[i].constraints[j][0] << " " << trans_sets[i].constraints[j][1] << " ";
        }
        cp_file << trans_sets[i].const_contradiction << "\n";
    }
}

void ReadCheckpointInfo (vector< vector<int> >& check_info) {
    ifstream check_file;
    check_file.open("Checkpoint1.out");
    int n;
    for (int i=0;i<1000000;i++) {
        vector<int> c;
        if (!(check_file >> n)) break;
        c.push_back(n);
        if (!(check_file >> n)) break;
        c.push_back(n);
        if (!(check_file >> n)) break;
        c.push_back(n);
        check_info.push_back(c);
    }
}

void ReadCheckpoint1 (int set, int iset, vector<int>& new_subset, vector<tpairs>& trans_sets) {
    ifstream cp_file;
    ostringstream convert;
    ostringstream convert2;
    convert << set;
    string temp=convert.str();
    convert2 << iset;
    string temp2=convert2.str();
    string name="Check1_"+temp+"_"+temp2+".out";
    cp_file.open(name.c_str());
    int n;
    int n1;
    int n2;
    int n3;
    int n4;
    tpairs t;
    for (int i=0;i<10000000;i++) {
        if (i==0) {
            //Read in indiviudals in the subset
            if (!(cp_file >> n)) break;
            for (int j=0;j<n;j++) {
                if (!(cp_file >> n1)) break;
                new_subset.push_back(n1);
            }
        }
        if (!(cp_file >> n)) break;
        t.root=n;
        if (!(cp_file >> n)) break;
        t.pair.clear();
        t.min_time.clear();
        t.max_time.clear();
        for (int j=0;j<n;j++) {
            if (!(cp_file >> n1)) break;
            if (!(cp_file >> n2)) break;
            if (!(cp_file >> n3)) break;
            if (!(cp_file >> n4)) break;
            tpair tp;
            tp.from=n1;
            tp.to=n2;
            t.pair.push_back(tp);
            t.min_time.push_back(n3);
            t.max_time.push_back(n4);
        }
        if (!(cp_file >> n)) break;
        t.constraints.clear();
        for (int j=0;j<n;j++) {
            vector<int> c;
            if (!(cp_file >> n1)) break;
            if (!(cp_file >> n2)) break;
            c.push_back(n1);
            c.push_back(n2);
            t.constraints.push_back(c);
        }
        if (!(cp_file >> n)) break;
        t.const_contradiction=n;
        trans_sets.push_back(t);
    }
}


void WriteCheckpoint2 (int set, int iset, vector<int> new_subset, vector<tpairs>& trans_sets) {
    ofstream cp_file;
    ostringstream convert;
    ostringstream convert2;
    convert << set;
    string temp=convert.str();
    convert2 << iset;
    string temp2=convert2.str();
    string name="Check2_"+temp+"_"+temp2+".out";
    cp_file.open(name.c_str());
    cp_file << new_subset.size() << " ";
    for (int i=0;i<new_subset.size();i++) {
        cp_file << new_subset[i] << " ";
    }
    cp_file << "\n";
    for (int i=0;i<trans_sets.size();i++) {
        if (trans_sets[i].order_index.size()>0) {
            cp_file << trans_sets[i].root << " " << trans_sets[i].pair.size() << " ";
            for (int j=0;j<trans_sets[i].pair.size();j++) {
                cp_file <<trans_sets[i].pair[j].from << " " << trans_sets[i].pair[j].to << " " << trans_sets[i].min_time[j] << " " << trans_sets[i].max_time[j] << " ";
            }
            cp_file << trans_sets[i].order_index.size() << " ";
            for (int j=0;j<trans_sets[i].order_index.size();j++) {
                cp_file << trans_sets[i].order_index[j] << " ";
            }
            cp_file << "\n";
        }
    }
}

void ReadCheckpoint2 (int set, int iset, vector<int>& new_subset, vector<tpairs>& trans_sets) {
    ifstream cp_file;
    ostringstream convert;
    ostringstream convert2;
    convert << set;
    string temp=convert.str();
    convert2 << iset;
    string temp2=convert2.str();
    string name="Check2_"+temp+"_"+temp2+".out";
    cp_file.open(name.c_str());
    int n;
    int n1;
    int n2;
    int n3;
    int n4;
    tpairs t;
    for (int i=0;i<10000000;i++) {
        if (i==0) {
            //Read in indiviudals in the subset
            if (!(cp_file >> n)) break;
            for (int j=0;j<n;j++) {
                if (!(cp_file >> n1)) break;
                new_subset.push_back(n1);
            }
        }
        if (!(cp_file >> n)) break;
        t.root=n;
        if (!(cp_file >> n)) break;
        t.pair.clear();
        t.min_time.clear();
        t.max_time.clear();
        for (int j=0;j<n;j++) {
            if (!(cp_file >> n1)) break;
            if (!(cp_file >> n2)) break;
            if (!(cp_file >> n3)) break;
            if (!(cp_file >> n4)) break;
            tpair tp;
            tp.from=n1;
            tp.to=n2;
            t.pair.push_back(tp);
            t.min_time.push_back(n3);
            t.max_time.push_back(n4);
        }
        if (!(cp_file >> n)) break;
        t.order_index.clear();
        for (int j=0;j<n;j++) {
            if (!(cp_file >> n1)) break;
            t.order_index.push_back(n1);
        }
        trans_sets.push_back(t);
    }
}


void WriteCheckpoint2Alt (int set, int iset, vector<int> new_subset, vector<tpairs>& trans_sets) {
    ofstream cp_file;
    ostringstream convert;
    ostringstream convert2;
    convert << set;
    string temp=convert.str();
    convert2 << iset;
    string temp2=convert2.str();
    string name="Check2_"+temp+"_"+temp2+"_alt.out";
    cp_file.open(name.c_str());
    cp_file << new_subset.size() << " ";
    for (int i=0;i<new_subset.size();i++) {
        cp_file << new_subset[i] << " ";
    }
    cp_file << "\n";
    for (int i=0;i<trans_sets.size();i++) {
        if (trans_sets[i].order_index.size()>0) {
            cp_file << trans_sets[i].trans_in.size() << " ";
            for (int j=0;j<trans_sets[i].trans_in.size();j++) {
                cp_file << trans_sets[i].trans_in[j] << " ";
            }
            cp_file << trans_sets[i].root << " " << trans_sets[i].pair.size() << " ";
            for (int j=0;j<trans_sets[i].pair.size();j++) {
                cp_file <<trans_sets[i].pair[j].from << " " << trans_sets[i].pair[j].to << " " << trans_sets[i].min_time[j] << " " << trans_sets[i].max_time[j] << " ";
            }
            cp_file << trans_sets[i].order_index.size() << " ";
            for (int j=0;j<trans_sets[i].order_index.size();j++) {
                cp_file << trans_sets[i].order_index[j] << " ";
            }
            cp_file << "\n";
        }
    }
}
