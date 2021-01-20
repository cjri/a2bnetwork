using namespace std;
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>


#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>

struct run_params {
	int dim;	//Model dimension
	int verb;	//Verbosity
	double ia;	//Parameters for the distribution between reported symptoms
	double ib;
	double io;
	double pa;	//Parameters for the infectious potential relative to time of symptoms
	double pb;
	double po;
	double smu;	//Parameters for the time of symptoms relative to date of infection
	double ssigma;
	double ucta; //Parameters for the distribution between becoming symptomatic and being tested
	double uctb; //For use when the time of becoming symptomatic is unknown
	double ucto;
	double uct_mean;  //Mean of the above distribution
	int opt_uct;	//Flag to optimise the reconstruction over the unknown times of being symptomatic
	double rate;	//Substitutions per site per year
	double seq_noise; //Mean number of errors in sequencing
	string ali_file;	//Name of alignment file
	string pat_file;	//Name of file with patient data i.e. code, date of symptoms, HCW status, etc.
	string mov_file;	//Name of file with HCW movement data
	string ward_file;	//Name of file with patient movement data
	string sub_file;	//Name of file with subset data
    string time_file;    //Name of file with transmission timings data
	int turner;	//Don't estimate uncertainty in tree
	int space; //Include spatial data for individuals
	int read_subs; //Read in subsets - can be used to split the inferred graph if needed
	int seqvar; //Flag to account for uncertainty over sequence composition
	int max_n; //Maximum number of ambiguous sites to allow per sequence.  Delete any sequence with more than this
	double threshold; //Likelihood threshold at which to cut the network into subsets
	double thresholdns; //Likelihood threshold for case of no sequence data.  Used in tree evaluation
    double threshold_singleLL; //Likelhood threshold for a single point in the time-dependent transmission likelihood.  Value of -20 reflects maximum around -4
	int rem_no_seq; //Flag to remove individuals for whom there is no sequence data
	char pat_delim; //Delimiter for date information in the patient data file
	char mov_delim; //Delimiter for date information in the HCW movement data file
	int diagnostic; //Flag to run various diagnostics - prints out various data along the way
	int noseq; //Flag to ignore genome seuqence data
	int utopia; //Flag to ignore
	int alt_like; //Flag for alternative likelihood model.  This increases the likelihood of transmission for each additional day of contact
	int calc_thresholds; //Calculate thresholds of the likelihood function
    int calculation; //Controls what the code does.  1: Go to checkpoint 1. 2: Read from checkpoint 1.
    
    //1. Find transmission networks
    int specify_set;
    int specify_remove;
    int consistency;

    //3. Calculate likelihoods
    int c_start;
    int c_step;
    
    
    //4. Process likelihood values
    int n_edges; //Number of edges to read in likelihoods.
    string likelihood_file; //Name of the likelihood output file to process in stage 4
};


struct loc {
	string ward;
	int date;
	double prob;
};

struct allele {
	int loc;
	char nuc;
	int count;
	double freq;
};

struct ipair {
	int loc;
	double q;
};

struct idat {
	int loc;
	double q;
	char allele;
};

struct tprob {
	int time;
	double weight;
};

struct ijlike {
	double lL_tot;
	double ns_lL_tot;
	int min; //Minimum time
	int max; //Maximum time
	vector<int> contact_times;
	vector<double> contact_likes;
	vector<double> noseq_likes;
};

struct pat {
	string code;
	string seq;
	string code_match;
	int hcw; //Flag to indicate health care worker (as opposed to patient) 1 = HCW; 0 = Patient
	int type; //Flag for type of infection.  1 = Community.  2 = Patient.  3 = HCW
	int time_s;  //Time of becoming symptomatic or time of testing positive
	int time_seq; //Time of collecting sample used in sequencing
	int dtime_s;  //If the time of becoming symptomatic is optimised, this gives the change in that time from the default
	int time_s_cert; //Certainty in time_s.  Flag 1 = Treat time_s as accurate.  0 = Asymptomatic or unknown.  time_s is the date of positive test
	int seq_n; //Number of N nucleotides in sequence
	vector<int> location; //Dates on which the individual is able to transmit or be infected due to their location
	vector<loc> locat; //Dates on which the individual is able to transmit or be infected due to their location
	vector<idat> seq_uncertainty; //Sites at which the sequence is not known
};

struct sparseseq {
    vector<int> locus;
    vector<char> allele;
};

struct varcount {
	int locus;
	char allele;
	vector<int> individual;
};

struct nbranch {
    int time;
    int mutations; //#Mutations observed
    vector<int> individuals; //Individuals with these mutations
};

struct branch {
	int from;
	int nbelow; //Number of individuals in clade
	vector<allele> vars;
	vector<int> indivs;
};

struct shared_branch {
	int index;
	int nbelow;
	vector<int> dependencies;
};

struct treestore {
	vector<sparseseq> variants;
	vector<string> individuals;
	vector<int> origins;
	vector<int> dests;
	vector<double> weights;
	double logL;
};

struct treestore_plus {
	vector<treestore> ts;
	int total_inc;
};

struct ed { //Directed from i to j
	int dest; //destination
	double weight;
	int origin; //Where this came from
};


struct tpair {
	int from;
	int to;
};

struct tpairs {  //Possible transmission structure
	int root;
	vector<tpair> pair; //List of transmissions
    vector<tpair> ordered_pair; //List of transmissions
    vector<int> trans_in; //To each node, from which individual?  -1 indicates root
	vector<int> min_time;
	vector<int> max_time;
	vector< vector<int> > constraints;
    int const_contradiction;
	vector< vector<int> > orders; //Possible orders of these transmissions
    vector<int> order_index;  //Index of orders from list
};

