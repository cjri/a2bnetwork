#include <stdio.h>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void WriteCheckpoint1 (int set, int iset, vector<int> new_subset, vector<tpairs>& trans_sets);
void WriteCheckpoint_Rem1 (int set, int iset, vector<int> new_subset, vector<tpairs>& trans_sets);
void ReadCheckpointInfo (vector< vector<int> >& check_info);
void ReadCheckpoint1 (int set, int iset, vector<int>& new_subset, vector<tpairs>& trans_sets);
void WriteCheckpoint2 (int set, int iset, vector<int> new_subset, vector<tpairs>& trans_sets);
void ReadCheckpoint2 (int set, int iset, vector<int>& new_subset, vector<tpairs>& trans_sets);
void WriteCheckpoint2Alt (int set, int iset, vector<int> new_subset, vector<tpairs>& trans_sets);
