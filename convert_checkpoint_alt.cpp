#include "basicmodel.h"
#include "checkpointing.h"

void ConvertCheckpointAlt() {
    vector< vector<int> > check_info;
    ReadCheckpointInfo(check_info);
    for (int i=0;i<check_info.size();i++) {
        int set=check_info[i][0];
        vector<int> new_subset;
        vector<tpairs> trans_sets;
        cout << "Reading checkpoint\n";
        ReadCheckpoint2 (set,check_info[i][1],new_subset,trans_sets);
        int n_nodes=trans_sets[0].pair.size()+1;
        vector<int> in_init;
        for (int j=0;j<n_nodes;j++)  {
            in_init.push_back(-1);
        }
        for (int j=0;j<trans_sets.size();j++) {
            trans_sets[j].trans_in=in_init;
            for (int k=0;k<trans_sets[j].pair.size();k++) {
                trans_sets[j].trans_in[trans_sets[j].pair[k].to]=trans_sets[j].pair[k].from;
            }
        }
        WriteCheckpoint2Alt (set,check_info[i][1],new_subset,trans_sets);
    }
}


