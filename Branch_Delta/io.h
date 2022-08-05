#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>


void GetOptions (run_params& p, int argc, const char **argv);
void GetThresholds (vector< vector<double> >& thresholds95, vector< vector<double> >& thresholds99, double& t95NS, double& t99NS, int& error);
void GetThresholdsInternal (run_params& p, vector< vector<double> >& thresholds95, vector< vector<double> >& thresholds99, double& t95NS, double& t99NS, int& error);

void SetThreshold (run_params& p);
//Basic data input
void ReadPatFromCSV(run_params p, vector<pat>& pdata);
void RemovePunc(string& str);
void SplitCommas(const string str, vector<string>& subs);
void MakeDMY (const int j, const vector<string>& subs, char delim, vector<int>& dmy);
void MakeDMY (run_params& p, string file, const int j, const vector<string>& subs, char delim, vector<int>& dmy);
void DateError (string file, string subsj, char delim);
void YearOOR (string file, string subsj);
void RemoveSpaces(string file, int i, int& pr, vector<string>& subs);
int DatetoDay (vector<int>& dmy);
//Genome sequence data input
void ReadFastaAli (run_params p, vector<string>& names, vector<string>& seqs);
void CheckBaseCase (vector<string>& seqs);
void MatchSequencesToIndividuals (const run_params p, vector<pat>& pdat, vector<string>& names, vector<string>& seqs);
void CorrectNames (vector<string>& names);
void IncorporateSequenceData (vector<pat>& pdat, vector<string>& names, vector<string>& seqs);
void RemoveIndividualsNoSequence(const run_params p, vector<pat>& pdat);
void IfNoSeqUseDefaultSequence (const run_params p, vector<pat>& pdat);
void CountNsInEachSequence (vector<pat>& pdat);
void ResolveFilterRepeatPatients (vector<pat>& pdat);
void RepairSequence(int i, int j, vector<pat>& pdat);
//Location data input
void ReadLocationData (run_params p, vector<pat>& pdat);
void EnforceUtopia (vector<pat>& pdat);
void EnforceModerateUtopia (vector<pat>& pdat);
void ReadHCWMovFromCSV(run_params& p, vector<pat>& pdat);
void ReadExtraHCWMovFromCSV(run_params& p, vector<pat>& pdat);
void ReadWardMovFromCSV(run_params p, vector<pat>& pdata);
void ReadWardMovBaysFromCSV(run_params& p, vector<pat>& pdata);
void ReadWardMovFromCSVNew(run_params& p, vector<pat>& pdata);
void RemoveIndividualsNoLocation (const run_params p, vector<pat>& pdat);
void ReadSubsets (run_params p, vector< vector<int> >& subsets);





void EditHCWMovData (run_params p, vector<pat>& pdat);

void ReadSymptomData (vector<pat>& pdat);
void ReadSeqTimes(vector<pat>& pdat);
void ReadLocationData (vector<pat>& pdat);
void ReadCommunityDistances (vector<int>& comm_dist);
//void WritePairwiseDistances (vector< vector<int> >& seqdists);
void WriteBestRoot (int index, const vector<treestore>& ts_roots);
void WriteMLDetails (run_params p, const int index, const vector< vector<int> >& bin, const vector<treestore_plus>& ts_bin, vector<pat>& pdat);
