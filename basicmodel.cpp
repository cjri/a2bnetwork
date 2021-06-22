#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#include "basicmodel.h"
#include "calc_timings.h"
#include "checkpointing.h"
#include "convert_checkpoint_alt.h"
#include "convert_likelihoods_alt.h"
#include "diagnostics.h"
#include "distributions.h"
#include "find_adjacent.h"
#include "find_orders.h"
#include "find_trans_networks.h"
#include "io.h"
#include "likelihoods.h"
#include "process_likelihoods.h"
#include "thresholds.h"
#include "threshold_data.h"
#include "utilities.h"
#include "variants.h"

int main(int argc, const char **argv) {

    //Random number initialisation
	int seed=(int) time(NULL);
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);
	
    //Object sets basic parameters for running the code
	run_params p;
	GetOptions(p,argc,argv);

    //Next step: Some calculations don't require read-in of the original data
    
    if (p.calculation==10) {
        //Convert second checkpoint file
        ConvertCheckpointAlt();
        return 0;
    }

    //Next calculation: Convert the likelihood file output to alt format
    if (p.calculation==11) {
        ConvertLikelihoodsAlt(p);
        return 0;
    }

    //Generate a list of adjacent networks, currently within two of the maximum
    //Exclude from the list 1.  Sites with likelihoods already counted.  2.  Sites with no orders
    if (p.calculation==14) {
        ListAdjacentK (p,p.distance);
        return 0;
    }

    if (p.calculation==16) {
		//Generate statistics detailing the timing of each transmission
        CalculateTimingStats(p);
        return 0;
    }
    
    //Main routine starts here
	
	//Flag to calculate thresholds for pairwise calculation
	if (p.calc_thresholds==1) {
		//Pre-calculate likelihoods: Speed up calculation
		vector<double> LNPreCalc; //LogNormal.  Assumes p.smu, p.ssigma parameters
		vector<double> OGPreCalcP; //Offset gamma
		PreCalculateLikelihoods(p,LNPreCalc,OGPreCalcP);
		CalculateThresholdsNoSeq(p,OGPreCalcP,LNPreCalc);
		CalculateThresholdsFullExplicit(p,OGPreCalcP,LNPreCalc,rgen);
		return 0;
	}

	int error=0;
	//External files required to set thresholds for pairwise calculation
	//If these are lost, they can be re-generated using the --calc_thresholds 1 flag
	GetThresholdsInternal (p.threshold95,p.threshold99,p.t95NS,p.t99NS,error);
	if (error==1) {
		return 0;
	}
	
	//Read in basic information from pat_file
	vector<pat> pdat;
	ReadPatFromCSV(p,pdat);
	
    //Read in sequence alignment
    vector<string> names;
    vector<string> seqs;
    ReadFastaAli(p,names,seqs);
    //Match sequence data to records
    MatchSequencesToIndividuals(p,pdat,names,seqs);
	cout << "Have complete data for " << pdat.size() << " individuals\n";

	//Code to read in location data
    ReadLocationData (p,pdat);
	cout << "Number of individuals now " << pdat.size() << "\n";

	//Convert sequences to variant data
	vector<sparseseq> variants;
    ProcessVariantInformation(p,seqs,variants,pdat);
	
	//Find distances between sequences
	vector< vector<int> > seqdists;
    vector< vector<tpair> > seqdists_c; //Distances to pairwise consensus
    FindPairwiseDistances (p,seqdists,variants,pdat);
    FromConsensusDistances (variants,seqdists_c);

	//Calculate pairwise likelihoods
	vector< vector<ijlike> > like_trans;
	CalculateTDLikelihoods (p,pdat,seqdists,seqdists_c,like_trans);
	
	//Find clusters of individuals
	cout << "Get subsets\n";
	vector< vector<int> > subsets;
	GetSubsetsIJ (p,like_trans,pdat,subsets);
	cout << "Here\n";
	
    //Gather data for calculation
    vector< vector<int> > new_subsets;
	vector< vector<pat> > pdat_sets;
	vector< vector< vector<ijlike> > > like_trans_sets;
	vector< vector<sparseseq> > variants_sets;
	CollectSubsetData (subsets,pdat,like_trans,variants,pdat_sets,like_trans_sets,variants_sets,new_subsets);

    if (p.calculation==1) {
        //Find plausible transmission networks and output to checkpoint
        FindPlausibleTransmissionNetworks (p,new_subsets,pdat_sets,like_trans_sets,variants_sets);
        return 0;
        
    }
    
    if (p.calculation==2) {
        FindPlausibleOrders();
        return 0;
    }
    
    if (p.calculation==3) {
        cout << "Run calculation 3\n";
        FindNetworkLikelihood (p,new_subsets,pdat_sets,like_trans_sets,variants_sets);
        return 0;
    }
    
    if (p.calculation==4) {
        ProcessLikelihoods(p);
        return 0;
    }
    
    if (p.calculation==12) {
        AnalyseAdjacent (p,pdat_sets,like_trans_sets,variants_sets);
        return 0;
    }
        
    
    if (p.calculation==15) {
        CalculateLikelihoodsFromList (p,pdat_sets,like_trans_sets,variants_sets);
        return 0;
    }

    
    
    return 0;
}
