using namespace std;
#include "basicmodel.h"
#include "checkpointing.h"
#include "diagnostics.h"
#include "distributions.h"
#include "find_orders.h"
#include "find_trans_networks.h"
#include "io.h"
#include "likelihoods.h"


void CalculateTDLikelihoods (run_params p, const vector<pat>& pdat, const vector< vector<int> >& seqdists, const vector< vector<tpair> >& seqdists_c, vector< vector<ijlike> >& like_trans) {
	for (int i=0;i<pdat.size();i++) {
		vector<ijlike> lt;
		for (int j=0;j<pdat.size();j++) {
			if (p.diagnostic==1) {
				cout << "i " << i << " j " << j << " ";
			}
			ijlike lij;
			lij.lL_tot=0;
			lij.ns_lL_tot=0;
			vector<tprob> contact_times_probs;
			if (i==j) {
				lij.lL_tot=-1e10;
				lij.ns_lL_tot=-1e10;
				lij.min=0;
				lij.max=-1;
			} else {
				for (int k=0;k<pdat[i].locat.size();k++) {
					for (int l=0;l<pdat[j].locat.size();l++) {
						if (pdat[i].locat[k].ward.compare(pdat[j].locat[l].ward)==0&&pdat[i].locat[k].date==pdat[j].locat[l].date) {
							//Same time and place
							tprob t;
							t.time=pdat[i].locat[k].date;
							t.weight=pdat[i].locat[k].prob*pdat[j].locat[l].prob;
							contact_times_probs.push_back(t);
						}
					}
				}
				//cout << "Number of positive contacts " << contact_times_probs.size() << "\n";
				if (contact_times_probs.size()==0) {
					lij.lL_tot=-1e10;
					lij.ns_lL_tot=-1e10;
					lij.min=0;
					lij.max=-1;
				} else {
					RemoveDuplicateTimes(contact_times_probs);
					lij.min=contact_times_probs[0].time;
					lij.max=contact_times_probs[contact_times_probs.size()-1].time;
					//cout << lij.min << " " << lij.max << "\n";
				
					FillTimes(contact_times_probs);
					for (int k=0;k<contact_times_probs.size();k++) {
						double lL=LikelihoodFromItoJTimeK (i,j,k,p,contact_times_probs,seqdists,seqdists_c,pdat);
                        //Add in filter here - all likelihoods must be at least a thrshold.  Remove everything else
                        if (lL>p.threshold_singleLL) {
                            lij.lL_tot=lij.lL_tot+exp(lL);
                            lij.contact_likes.push_back(lL);
                        } else {
                            lij.contact_likes.push_back(-1e10);
                        }
						lij.contact_times.push_back(contact_times_probs[k].time);
						lij.contact_likes.push_back(lL);
						double lLnoseq=NoSeqLikelihoodFromItoJTimeK (i,j,k,p,contact_times_probs,pdat);
                        if (lL>p.threshold_singleLL){
                            lij.ns_lL_tot=lij.ns_lL_tot+exp(lLnoseq);
                            lij.noseq_likes.push_back(lLnoseq);
                        } else {
                            lij.noseq_likes.push_back(-1e10);
                        }
					}
						
					//Trim lij values.  Subtract large and negative values from the end of the time.
					//cout << "Do trim\n";
					TrimLikelihoods(lij);
					//cout << "Done\n";
					if (lij.lL_tot>0) {
						lij.lL_tot=log(lij.lL_tot);
					} else {
						lij.lL_tot=-1e10;
					}
					if (lij.ns_lL_tot>0) {
						lij.ns_lL_tot=log(lij.ns_lL_tot);
					} else {
						lij.ns_lL_tot=-1e10;
					}

					//cout << "Here\n";
				}
			}
			//Standardise impossible tranmissions - use this in optimising times to skip possibilities
			//cout << "Standardise\n";
			if (lij.lL_tot<-1e9) {
				lij.lL_tot=-1e10;
			}
			//cout << "Pij\n";
			lt.push_back(lij);
			//cout << "Done p\n";
			if (p.diagnostic==1) {
                cout << lij.min << " " << lij.max << " ";
				cout << lij.lL_tot << " " << lij.ns_lL_tot << "\n";
			}
			//cout << "Now j " << j << "\n";
		}
		like_trans.push_back(lt);
	}
}

bool comparetprob (tprob v1, tprob v2) {
    return (v1.time < v2.time);
}

void RemoveDuplicateTimes (vector<tprob>& contact_times_probs) {
	sort(contact_times_probs.begin(),contact_times_probs.end(),comparetprob);
	vector<int> to_rem;
	for (int k=0;k<contact_times_probs.size()-1;k++) {
		if (contact_times_probs[k].time==contact_times_probs[k+1].time) {
			if (contact_times_probs[k].weight<contact_times_probs[k+1].weight) {
				to_rem.push_back(k);
			} else {
				to_rem.push_back(k+1);
			}
		}
	}
	sort(to_rem.begin(),to_rem.end());
	reverse(to_rem.begin(),to_rem.end());
	for (int k=0;k<to_rem.size();k++) {
		contact_times_probs.erase(contact_times_probs.begin()+to_rem[k]);
	}
}

void FillTimes (vector<tprob>& contact_times_probs) {
	int size=contact_times_probs.size();
	int pos=0;
	while (pos<size-1) {
		if (contact_times_probs[pos+1].time!=contact_times_probs[pos].time+1) {
			tprob t;
			t.time=contact_times_probs[pos].time+1;
			t.weight=0;
			contact_times_probs.insert(contact_times_probs.begin()+pos+1,t);
			size++;
		}
		pos++;
	}
}

double LikelihoodFromItoJTimeK (int i, int j, int k, const run_params& p, const vector<tprob>& contact_times_probs, const vector< vector<int> >& seqdists, const vector< vector<tpair> >& seqdists_c, const vector<pat>& pdat) {
    //Transmission from i to j at time contact_times[k].  Integrate sequence evolution date with key likelihood.
    double L=0;
    //Sequence component of the likelihoods
    if (pdat[j].time_seq<contact_times_probs[k].time) { //Data must be collected from j after the time of transmission
        L=L-1e10; //Can't transmit after collecting data from recipient
    } else {
        if (p.noseq==0) {
            //Two sequence likelihoods -
            if (pdat[i].time_seq<contact_times_probs[k].time) {//First, case where transmission is after collecting the first sequence sample
                //cout << "First\n";
                int evo_time=(pdat[j].time_seq-pdat[i].time_seq);//Amount of time over which evolution has occurred
                double evo=evo_time*p.rate;
                if (seqdists[i][j]==-1) {
                    L=L+Poisson(9,evo+p.seq_noise);
                } else {
                    L=L+Poisson(seqdists_c[i][j].from,(p.seq_noise/2)); //D1 data - any novel variants must arise from noise
                    L=L+Poisson(seqdists_c[i][j].to,evo+(p.seq_noise/2)); //Evolution between D1 and D2
                }
            } else { //Second, case where transmission is before collecting the first sequence sample
                //cout << "Second\n";
                //Calculate times
                int evo_time1=pdat[i].time_seq-contact_times_probs[k].time;
                int evo_time2=pdat[j].time_seq-contact_times_probs[k].time;
                double evo1=evo_time1*p.rate;
                double evo2=evo_time2*p.rate;
                if (seqdists[i][j]==-1) {
                    L=L+Poisson(9,evo1+(p.seq_noise/2));
                    L=L+Poisson(9,evo2+(p.seq_noise/2));
                } else {
                    L=L+Poisson(seqdists_c[i][j].from,evo1+(p.seq_noise/2)); //Evolution from transmission to D1
                    L=L+Poisson(seqdists_c[i][j].to,evo2+(p.seq_noise/2)); //Evolution from transmission to D2
                }
            }
        }
        double OGL=OffsetGammaCDFFlex(contact_times_probs[k].time-pdat[i].time_s,p.pa,p.pb,p.po); //Time from symptom onset of i to transmission
        L=L+OGL;
        if (contact_times_probs[k].weight==0) {
            L=L-1e10; //Can't transmit when there is no contact between individuals
        } else {
            L=L+log(contact_times_probs[k].weight); //Probability of contact between individuals
        }
        double LN=LogNormal(pdat[j].time_s-contact_times_probs[k].time,p.smu,p.ssigma);//Time from transmission to j becoming symptomatic
        L=L+LN;
    }
    //cout << "Log " << L << "\n";
    return L;
}


double NoSeqLikelihoodFromItoJTimeK (int i, int j, int k, const run_params& p, const vector<tprob>& contact_times_probs, const vector<pat>& pdat) {
    //Transmission from i to j at time contact_times[k].  Integrate sequence evolution date with key likelihood.
    double L=0;
    double OGL=OffsetGammaCDFFlex(contact_times_probs[k].time-pdat[i].time_s,p.pa,p.pb,p.po); //Time from symptom onset of i to transmission
    L=L+OGL;
    if (contact_times_probs[k].weight==0) {
        L=L-1e10; //Can't transmit when there is no contact between individuals
    } else {
        L=L+log(contact_times_probs[k].weight); //Probability of contact between individuals
    }
    double LN=LogNormal(pdat[j].time_s-contact_times_probs[k].time,p.smu,p.ssigma);//Time from transmission to j becoming symptomatic
    L=L+LN;
    //cout << "Log " << L << "\n";
    return L;
}

void TrimLikelihoods (ijlike& lij) {
    //Remove impossible events from the ends : Be more efficient later
    //cout << lij.contact_times.size() << " " << lij.contact_likes.size() << " " << lij.noseq_likes.size() << "\n";
    vector<int> to_rem;
    for (int k=0;k<lij.contact_times.size();k++) {
        if (lij.noseq_likes[k]<-1e5) {
            to_rem.push_back(k);
            //cout << k << "\n";
        } else {
            break;
        }
    }
    sort(to_rem.begin(),to_rem.end());
    reverse(to_rem.begin(),to_rem.end());
    for (int k=0;k<to_rem.size();k++) {
        //cout << k << " " << to_rem[k] << " " << lij.contact_times.size() << " " << lij.contact_likes.size() << " " << lij.noseq_likes.size() << "\n";
        lij.contact_times.erase(lij.contact_times.begin()+to_rem[k]);
        lij.contact_likes.erase(lij.contact_likes.begin()+to_rem[k]);
        lij.noseq_likes.erase(lij.noseq_likes.begin()+to_rem[k]);
        lij.min++;
        //cout << lij.max << "\n";
    }
    to_rem.clear();
    for (int k=lij.contact_times.size()-1;k>=0;k--) {
        if (lij.noseq_likes[k]<-1e5) {
            to_rem.push_back(k);
            //cout << k << "\n";
        } else {
            break;
        }
    }
    sort(to_rem.begin(),to_rem.end());
    reverse(to_rem.begin(),to_rem.end());
    for (int k=0;k<to_rem.size();k++) {
        //cout << k << " " << to_rem[k] << " " << lij.contact_times.size() << " " << lij.contact_likes.size() << " " << lij.noseq_likes.size() << "\n";
        lij.contact_times.erase(lij.contact_times.begin()+to_rem[k]);
        lij.contact_likes.erase(lij.contact_likes.begin()+to_rem[k]);
        lij.noseq_likes.erase(lij.noseq_likes.begin()+to_rem[k]);
        lij.max--;
        //cout << lij.max << "\n";
    }
    //cout << lij.contact_times.size() << " " << lij.contact_likes.size() << " " << lij.noseq_likes.size() << "\n";
}

void GetSubsetsIJ (run_params p, const vector< vector<ijlike> >& like_trans, const vector<pat>& pdat, vector< vector<int> >& subsets) {
    if (p.read_subs==1) {
        cout << "Read in subset data\n";
        ReadSubsets(p,subsets);
    } else {
        FindLikelihoodSubsetsIJ(p,like_trans,subsets);
    }
    cout << "Subsets\n";
    for (int i=0;i<subsets.size();i++) {
        for (int j=0;j<subsets[i].size();j++) {
            cout << subsets[i][j] << " " << pdat[subsets[i][j]].code << " ";
        }
        cout << "\n";
    }
    cout << "Number of subsets " << subsets.size() << "\n";
}

void FindLikelihoodSubsetsIJ(run_params p, const vector< vector<ijlike> >& like_trans, vector< vector<int> >& subsets) {
    vector <int> range;
    for (int i=0;i<like_trans.size();i++) {
        range.push_back(i);
    }
    while (range.size()>0) {
        vector<int> sset;
        sset.clear();
        sset.push_back(range[0]);
        int add=1;
        while (add==1) {
            int ss=sset.size();
            add=0;
            for (int i=0;i<ss;i++) {
                for (int j=0;j<like_trans.size();j++) {
//                    if (like_trans[sset[i]][j].lL_tot>p.threshold&&like_trans[sset[i]][j].ns_lL_tot>-13.7123) {
                    if (like_trans[sset[i]][j].lL_tot>p.threshold) {
                        //cout << "Add " << sset[i] << " " << j << "\n";
                        sset.push_back(j);
                    }
    //                if (like_trans[j][sset[i]].lL_tot>p.threshold&&like_trans[j][sset[i]].ns_lL_tot>-13.7123) {
                    if (like_trans[j][sset[i]].lL_tot>p.threshold) {
                        //cout << "Add " << sset[i] << " " << j << "\n";
                        sset.push_back(j);
                    }
                }
            }
            sort(sset.begin(),sset.end());
            sset.erase(unique(sset.begin(),sset.end()),sset.end());
            //cout << sset.size() << "\n";
            if (sset.size()>ss) {
                add=1;
            }
        }
        //Add sset to subsets
        //cout << "Push back to subsets\n";
        subsets.push_back(sset);
        for (int i=0;i<sset.size();i++) {
            for (int j=0;j<range.size();j++) {
                if (range[j]==sset[i]) {
                    range.erase(range.begin()+j);
                    break;
                }
            }
        }
    }
}

void FindNetworkLikelihood (run_params p, vector< vector<int> >& new_subsets, vector< vector<pat> >& pdat_sets, vector< vector< vector<ijlike> > >& like_trans_sets, vector< vector<sparseseq> >& variants_sets) {
    //Calculation to get likelihoods for networks.  Calculation 3.
    vector< vector<int> > check_info;
    ReadCheckpointInfo(check_info);
    for (int i=0;i<check_info.size();i++) {
        int set=check_info[i][0];
        vector<int> new_subset;
        vector<tpairs> trans_sets;
        cout << "Reading checkpoint\n";
        ReadCheckpoint2 (set,check_info[i][1],new_subset,trans_sets);

        vector< vector<int> > orders;
        vector<pat> new_pdat;
        vector< vector<ijlike> > new_like_trans;
        vector<sparseseq> new_variants;
        vector<varcount> vc_multi;
        vector<varcount> vc_single;

        SetupTreeLikelihoodCalculation(p,set,new_subset,trans_sets,pdat_sets,like_trans_sets,variants_sets,orders,new_pdat,new_like_trans,new_variants,vc_multi,vc_single);
        
        //Set up output file or files
        ofstream cp_file;
        ostringstream convert;
        ostringstream convert2;
        convert << check_info[i][0];
        string temp=convert.str();
        convert2 << check_info[i][1];
        string temp2=convert2.str();
        string name="Likelihoods_"+temp+"_"+temp2+".out";
        cp_file.open(name.c_str());

        //Loop over transmission chains
        for (int c=p.c_start;c<trans_sets.size();c=c+p.c_step) {
            //Make order c from order index
            MakeOrderFromIndex (c,orders,trans_sets);
            cout << "Set " << c << "\n";
            double chain_L=0;
            if (p.diagnostic==1) {
                //ChainOutput(c,trans_sets,new_pdat);
            }
            //Loop over orderings
            for (int o=0;o<trans_sets[c].orders.size();o++) {
                CalculateOrderLike (p,c,o,chain_L,trans_sets,vc_multi,vc_single,new_pdat,new_like_trans);
            }
            cp_file << "Network " << c << " ";
            for (int k=0;k<trans_sets[c].ordered_pair.size();k++) {//Action to be done here is outside of next loop
                cp_file << trans_sets[c].ordered_pair[k].from << "-" << trans_sets[c].ordered_pair[k].to << " ";
            }
            cp_file << "Likelihood " << c << " " << log(chain_L) << "\n";
            //End of looking at chain c
            //Clear orderings
            trans_sets[c].orders.clear();
        }
        cp_file.close();
    }
}

void SetupTreeLikelihoodCalculation(run_params p, int set, const vector<int>& new_subset, vector<tpairs>& trans_sets, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, const vector< vector<sparseseq> >& variants_sets, vector< vector<int> >& orders, vector<pat>& new_pdat, vector< vector<ijlike> >& new_like_trans, vector<sparseseq>& new_variants, vector<varcount>& vc_multi, vector<varcount>& vc_single) {
    //Make orderings list.  All possible orderings
    cout << "Make orderings\n";
    ConstructOrderingsNew(trans_sets,orders);

    //Compile data
    SetupParamsSubsetAnalysis (set,new_subset,pdat_sets,like_trans_sets,variants_sets,new_pdat,new_like_trans,new_variants);
    
    //Mutations are tied to individuals - constant for all chains.
    //Identify mutations here
    FindSharedVariantsNew(vc_multi,vc_single,new_variants);
    
    //Edit multi mutations - remove mutations shared by all individuals
    EditSharedVariants(trans_sets,vc_multi);
    
    if (p.diagnostic==1){
        PrintSharedSingleVariants(vc_multi,vc_single);
    }
}

void EditSharedVariants(const vector<tpairs>& trans_sets, vector<varcount>& vc_multi) {
    vector<int> to_rem;
    for (int i=0;i<vc_multi.size();i++) {
        if (vc_multi[i].individual.size()==trans_sets[0].pair.size()+1) {
            to_rem.push_back(i);
        }
    }
    sort(to_rem.begin(),to_rem.end());
    reverse(to_rem.begin(),to_rem.end());
    for (int i=0;i<to_rem.size();i++) {
        vc_multi.erase(vc_multi.begin()+to_rem[i]);
    }
}

void MakeOrderFromIndex (int c, const vector< vector<int> >& orders, vector<tpairs>& trans_sets) {
    for (int k=0;k<trans_sets[c].order_index.size();k++){
        trans_sets[c].orders.push_back(orders[trans_sets[c].order_index[k]]);
    }
}

//Loop for single order
void CalculateOrderLike (run_params p, int c, int o, double& chain_L, vector<tpairs>& trans_sets, const vector<varcount>& vc_multi, const vector<varcount>& vc_single, const vector<pat>& new_pdat, vector< vector<ijlike> >& new_like_trans) {
   // cout << "CalculateOrderLike\n";
    double order_L=0;
    if (p.diagnostic==1) {
        //OrderOutput(c,o,trans_sets);
    }
    MakeOrderedPairs(c,o,trans_sets);
    vector< vector<int> > descendents; //Descendents of each transmission. Individual i is 100+i, node is j
    MakeDescendents(c,trans_sets,descendents);
    vector< vector< vector<int> > > below; //Which individuals are under each node
    MakeBelow (descendents,below);
    vector< vector<int> > sources; //Origin in each transmission
    MakeSources(c,trans_sets,sources);
    if (p.diagnostic==1) {
        //PrintDescendentsSourcesBelow(descendents,below,sources);
    }
    //Set up time vectors
    vector<int> times;
    vector<int> max_times;
    vector<int> min_times;
    vector<int> equiv_times;
    SetUpTimeVectors(c,o,trans_sets,new_like_trans,times,min_times,max_times,equiv_times);
    //Alter times so that they follow rules.  Must be in order.  One different for A-B then B-C
    CheckTimeOrder(c,o,trans_sets,times);

    //Print out time boundaries
    if (p.diagnostic==1) {
        //PrintTimings(c,trans_sets,new_like_trans);
        //PrintLikelihoods(c,trans_sets,new_like_trans);
        //PrintTimes(times);
    }

    //Check whether the initial time is feasible
    int feasible=CheckFeasible(times,max_times);
    if (feasible==1) {
        if (p.diagnostic==1) {
            //PrintTimes(equiv_times);
        }
        GenerateLikelihoodAllTimings (c,o,p,order_L,times,min_times,max_times,equiv_times,vc_single,vc_multi,sources,descendents,below,new_pdat,trans_sets,new_like_trans);
    }
    chain_L=chain_L+order_L;
    //cout << "Chain " << c << " likelihood is now " << chain_L << "\n";
}

void MakeOrderedPairs (int c, int o, vector<tpairs>& trans_sets) {
    trans_sets[c].ordered_pair.clear();
    for (int k=0;k<trans_sets[c].orders[o].size();k++) {
        trans_sets[c].ordered_pair.push_back(trans_sets[c].pair[trans_sets[c].orders[o][k]]);
    }
}

void MakeDescendents(int c, const vector<tpairs>& trans_sets, vector< vector<int> >& descendents) {
    for (int k=0;k<trans_sets[c].ordered_pair.size();k++) {
        vector<int> d;
        d.push_back(-1);
        d.push_back(-1);
        descendents.push_back(d);
    }
    for (int k=0;k<trans_sets[c].ordered_pair.size();k++) {
        for (int l=0;l<k;l++) {
            for (int m=0;m<=1;m++) {
                if (descendents[l][m]==trans_sets[c].ordered_pair[k].from+100) {
                    descendents[l][m]=k;
                }
            }
        }
        descendents[k][0]=trans_sets[c].ordered_pair[k].from+100;
        descendents[k][1]=trans_sets[c].ordered_pair[k].to+100;
    }
}

void MakeBelow (const vector< vector<int> >& descendents, vector< vector< vector<int> > >& below) {
    //Initialise object
    vector< vector<int> > be;
    vector<int> b;
    be.push_back(b);
    be.push_back(b);
    for (int k=0;k<descendents.size();k++) {
        below.push_back(be);
    }
    //Make object.  Individuals beneath each side of each node
    for (int k=descendents.size()-1;k>=0;k--) {
        if (descendents[k][0]>99) {
            below[k][0].push_back(descendents[k][0]-100);
        } else {
            for (int l=0;l<below[descendents[k][0]][0].size();l++) {
                below[k][0].push_back(below[descendents[k][0]][0][l]);
            }
            for (int l=0;l<below[descendents[k][0]][1].size();l++) {
                below[k][0].push_back(below[descendents[k][0]][1][l]);
            }
        }
            
        if (descendents[k][1]>99) {
            below[k][1].push_back(descendents[k][1]-100);
        } else {
            for (int l=0;l<below[descendents[k][1]][0].size();l++) {
                below[k][1].push_back(below[descendents[k][1]][0][l]);
            }
            for (int l=0;l<below[descendents[k][1]][1].size();l++) {
                below[k][1].push_back(below[descendents[k][1]][1][l]);
            }
        }
    }
    for (int i=0;i<below.size();i++) {
        for (int j=0;j<below[i].size();j++) {
            sort(below[i][j].begin(),below[i][j].end());
        }
    }
    
}

void MakeSources (int c, const vector<tpairs>& trans_sets, vector< vector<int> >& sources) {
    for (int k=0;k<trans_sets[c].ordered_pair.size();k++) {
        vector<int> s;
        s.push_back(trans_sets[c].ordered_pair[k].from);
        s.push_back(trans_sets[c].ordered_pair[k].to);
        sources.push_back(s);
    }
}

void SetUpTimeVectors(int c, int o, const vector<tpairs>& trans_sets, const vector< vector<ijlike> >& new_like_trans, vector<int>& times, vector<int>& min_times, vector<int>& max_times, vector<int>& equiv_times) {
    for (int k=0;k<trans_sets[c].ordered_pair.size();k++) {
        times.push_back(new_like_trans[trans_sets[c].ordered_pair[k].from][trans_sets[c].ordered_pair[k].to].min);
        max_times.push_back(new_like_trans[trans_sets[c].ordered_pair[k].from][trans_sets[c].ordered_pair[k].to].max);
    }
    min_times=times;

    //Initialise times vector by ordering
    CheckTimeOrder(c,o,trans_sets,times);

    //Correct the max times here to preserve ordering
    CheckMaxTimes (c,o,trans_sets,max_times);
    
    for (int k=0;k<trans_sets[c].ordered_pair.size();k++) {
        equiv_times.push_back(times[k]-min_times[k]);
    }
}

void CheckTimeOrder(int c, int o, const vector<tpairs>& trans_sets, vector<int>& times) {
    //Initial times must follow ordering rules
    for (int k=1;k<times.size();k++) {
        while (times[k]<times[k-1]) {
            times[k]++;
        }
        int onediff=0;
        if (trans_sets[c].ordered_pair[k].from==trans_sets[c].ordered_pair[k-1].to) {
            onediff=1;
        }
        if (trans_sets[c].orders[o][k]<trans_sets[c].orders[o][k-1]) {
            onediff=1;
        }
        if (onediff==1&&times[k]==times[k-1]) {
            times[k]++;
        }
    }

}

void CheckMaxTimes (int c, int o, const vector<tpairs>& trans_sets, vector<int>& max_times) {
    //Correct the max_times vector - reduce so as not to break the ordering rules
    for (int k=max_times.size()-2;k>=0;k--) {
        while (max_times[k]>max_times[k+1]) {
            max_times[k]--;
        }
        int onediff=0;
        if (trans_sets[c].ordered_pair[k].to==trans_sets[c].ordered_pair[k+1].from) {
            onediff=1;
        }
        if (trans_sets[c].orders[o][k]>trans_sets[c].orders[o][k+1]) {
            onediff=1;
        }
        if (onediff==1&&max_times[k]==max_times[k+1]) {
            max_times[k]--;
        }
    }
}

int CheckFeasible (const vector<int>& times, const vector<int>& max_times) {
    int feasible=1;
    for (int k=0;k<times.size();k++) {
        if (times[k]>max_times[k]) {
            feasible=0;
            break;
        }
    }
    return feasible;
}

void GenerateLikelihoodAllTimings (int c, int o, run_params& p, double& order_L, vector<int>& times, const vector<int>& min_times, const vector<int>& max_times, vector<int>& equiv_times, const vector<varcount>& vc_single, const vector<varcount>& vc_multi, const vector< vector<int> >& sources, const vector< vector<int> >& descendents, const vector< vector< vector<int> > >& below, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets, const vector< vector<ijlike> >& new_like_trans) {
    vector<int> orig_times=times;
    int fin=0;
    int change=-1;
    while (fin==0) {
        int done=0;
        int up=times.size()-1;
        if (p.diagnostic==1) {
            //PrintTransSampleTimes (c,times,new_pdat,trans_sets);
        }
        
        //Construct tree: Find window of time in which mutations affect each set of individuals
        vector<nbranch> treespace;
        CreateTimeTree(c,times,sources,descendents,below,new_pdat,trans_sets,treespace);
        //Assign mutations to time tree
        AssignMutationsToBranches (vc_multi,vc_single,treespace);
        if (p.diagnostic==1) {
            //PrintTreespace(treespace);
        }
            
        FindEquivalentTimes(times,min_times,equiv_times);
        //Calcualte likelihood
        CalculateTreeLikelihood (c,p,order_L,times,equiv_times,treespace,trans_sets,new_like_trans);
        //Go to the next set of times.
        FindNextTimeSet (c,o,p,done,up,fin,min_times,max_times,orig_times,trans_sets,new_like_trans,times);
        
    }
}

void CreateTimeTree(int c, const vector<int>& times, const vector< vector<int> >& sources, const vector< vector<int> >& descendents, const vector< vector< vector<int> > >& below, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets, vector<nbranch>& treespace) {
    int index=descendents.size()-1;
    while (index>=0) {
        if (descendents[index][0]>99) {
            if (descendents[index][1]>99) {
                //Branch to two individuals
                BranchTwoIndividuals(c,index,times,sources,descendents,new_pdat,trans_sets,treespace);
            } else {
                BranchIndividualNode(c,index,times,sources,descendents,below,new_pdat,trans_sets,treespace);
            }
        } else {
            if (descendents[index][1]>99) {
                BranchNodeIndividual(c,index,times,sources,descendents,below,new_pdat,trans_sets,treespace);
            } else {
                BranchTwoNodes(c,index,times,sources,descendents,below,new_pdat,trans_sets,treespace);
            }
        }
        index--;
    }
}

void BranchTwoIndividuals(int c, int index, const vector<int>& times, const vector< vector<int> >& sources, const vector< vector<int> >& descendents, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets, vector<nbranch>& treespace) {
    int source=sources[index][0];
    int dest=sources[index][1];
    nbranch n_source;
    nbranch n_dest;
    if (new_pdat[source].time_seq>=times[index]) {
        //Both measured after transmission
        n_source.time=new_pdat[source].time_seq-times[index];
        n_source.individuals.push_back(source);
        treespace.push_back(n_source);

        n_dest.time=new_pdat[dest].time_seq-times[index];
        n_dest.individuals.push_back(dest);
        treespace.push_back(n_dest);
    } else {
        //Source not measured after transmission
        n_source.time=0;
        n_source.individuals.push_back(source);
        treespace.push_back(n_source);

        //Find out where to go back to
        int e=FindPrevious (c,source,index,times,new_pdat,trans_sets);
        n_dest.time=new_pdat[dest].time_seq-e;
        n_dest.individuals.push_back(dest);
        treespace.push_back(n_dest);
    }

  
    /*cout << "BranchII\n";

    cout << source << " " << n_source.time << " {";
    for (int k=0;k<n_source.individuals.size();k++) {
        cout << n_source.individuals[k] << ", ";
    }
    cout << "}\n";
    cout << dest << " " << n_dest.time << " {";
    for (int k=0;k<n_dest.individuals.size();k++) {
        cout << n_dest.individuals[k] << ", ";
    }
    cout << "}\n";*/
}

void BranchIndividualNode(int c, int index, const vector<int>& times, const vector< vector<int> >& sources, const vector< vector<int> >& descendents, const vector< vector< vector<int> > >& below, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets, vector<nbranch>& treespace) {
    //Destination is a node, source is not
    int source=sources[index][0];
    int dest=sources[index][1];
    /*cout << "Source " << sources[index][0] << " dest " << sources[index][1] << "\n";
    cout << "Descendents " << descendents[index][0] << " " << descendents[index][1]  << "\n";
    cout << "Source measured at " << new_pdat[source].time_seq << " node time is " << times[index] << "\n";*/
    nbranch n_source;
    nbranch n_dest;
    if (new_pdat[source].time_seq>=times[index]) {
        //Source measured after transmission
        n_source.time=new_pdat[source].time_seq-times[index];
        n_source.individuals.push_back(source);
        treespace.push_back(n_source);
        
        n_dest.individuals=below[index][1];
        //cout << "Branch time is " << times[descendents[index][1]] << " node is " << times[index] << "\n";
        if (new_pdat[source].time_seq>=times[descendents[index][1]]) {
            n_dest.time=times[descendents[index][1]]-times[index];
        } else {
            n_dest.time=new_pdat[dest].time_seq-times[index];
        }
        treespace.push_back(n_dest);

    } else {
        //Source measured before transmission
        n_source.time=0;
        n_source.individuals.push_back(source);
        treespace.push_back(n_source);

        //Find where to go back to
        int e=FindPrevious (c,source,index,times,new_pdat,trans_sets);
        n_dest.individuals=below[index][1];
        //cout << "Previous time is " << e << "\n";
        //cout << "That branch " << times[descendents[index][1]] << " compare time_seq " <<new_pdat[source].time_seq << "\n";
        if (new_pdat[source].time_seq>=/*Time of that branch*/times[descendents[index][1]]) {
            n_dest.time=times[descendents[index][1]]-e;
        } else {
            n_dest.time=new_pdat[dest].time_seq-e;
        }
        treespace.push_back(n_dest);
    }
    /*cout << "BranchIB\n";

    cout << n_source.time << " {";
    for (int k=0;k<n_source.individuals.size();k++) {
        cout << n_source.individuals[k] << ", ";
    }
    cout << "}\n";
    cout << n_dest.time << " {";
    for (int k=0;k<n_dest.individuals.size();k++) {
        cout << n_dest.individuals[k] << ", ";
    }
    cout << "}\n";*/
}

void BranchNodeIndividual (int c, int index, const vector<int>& times, const vector< vector<int> >& sources, const vector< vector<int> >& descendents, const vector< vector< vector<int> > >& below, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets, vector<nbranch>& treespace) {
    //Source is a branch, destination is not
    int source=sources[index][0];
    int dest=sources[index][1];
    
    nbranch n_source;
    if (new_pdat[source].time_seq>=times[descendents[index][0]]) {
        n_source.time=times[descendents[index][0]]-times[index];
        n_source.individuals=below[index][0];
    } else if (new_pdat[source].time_seq<times[descendents[index][0]]&&new_pdat[source].time_seq>=times[index]) {
        n_source.time=new_pdat[source].time_seq-times[index];
        n_source.individuals=below[index][0];
    } else {
        int e=FindPrevious (c,source,index,times,new_pdat,trans_sets);
        n_source.time=times[index]-e;
        n_source.individuals.push_back(dest);
        for (int k=0;k<below[index][0].size();k++) {
            if (below[index][0][k]!=source) {
                n_source.individuals.push_back(below[index][0][k]);
            }
        }
    }
    treespace.push_back(n_source);
    
    nbranch n_dest;
    n_dest.time=new_pdat[dest].time_seq-times[index];
    n_dest.individuals.push_back(dest);
    treespace.push_back(n_dest);
    
    /*cout << "BranchBI\n";

    cout << n_source.time << " {";
    for (int k=0;k<n_source.individuals.size();k++) {
        cout << n_source.individuals[k] << ", ";
    }
    cout << "}\n";
    cout << n_dest.time << " {";
    for (int k=0;k<n_dest.individuals.size();k++) {
        cout << n_dest.individuals[k] << ", ";
    }
    cout << "}\n";*/

}

void BranchTwoNodes (int c, int index, const vector<int>& times, const vector< vector<int> >& sources, const vector< vector<int> >& descendents, const vector< vector< vector<int> > >& below, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets, vector<nbranch>& treespace) {
    int source=sources[index][0];
    int dest=sources[index][1];
    nbranch n_source;
    if (new_pdat[source].time_seq>=times[descendents[index][0]]) {
        n_source.time=times[descendents[index][0]]-times[index];
        n_source.individuals=below[index][0];
    } else if (new_pdat[source].time_seq<times[descendents[index][0]]&&new_pdat[source].time_seq>=times[index]) {
        n_source.time=new_pdat[source].time_seq-times[index];
        n_source.individuals=below[index][0];
    } else {
        int e=FindPrevious (c,source,index,times,new_pdat,trans_sets);
        n_source.time=times[index]-e;
        for (int k=0;k<below[index][1].size();k++) {
            n_source.individuals.push_back(below[index][1][k]);
        }
        for (int k=0;k<below[index][0].size();k++) {
            if (below[index][0][k]!=source) {
                n_source.individuals.push_back(below[index][0][k]);
            }
        }
    }
    treespace.push_back(n_source);

    nbranch n_dest;
    if (new_pdat[dest].time_seq>=times[descendents[index][1]]) {
        n_dest.time=times[descendents[index][1]]-times[index];
    } else {
        n_dest.time=new_pdat[dest].time_seq-times[index];
    }
    n_dest.individuals=below[index][1];
    treespace.push_back(n_dest);

    /*cout << "BranchBB\n";

    cout << n_source.time << " {";
    for (int k=0;k<n_source.individuals.size();k++) {
        cout << n_source.individuals[k] << ", ";
    }
    cout << "}\n";
    cout << n_dest.time << " {";
    for (int k=0;k<n_dest.individuals.size();k++) {
        cout << n_dest.individuals[k] << ", ";
    }
    cout << "}\n";*/

}


int FindPrevious (int c, int source, int index, const vector<int>& times, const vector<pat>& new_pdat, const vector<tpairs>& trans_sets) {
    int measure=new_pdat[source].time_seq;
    //cout << "Measure " << new_pdat[source].time_seq << "\n";
    int prev=-10000;
    for (int k=index-1;k>=0;k--) {
        if (trans_sets[c].ordered_pair[k].from==source) {
            prev=times[k];
        }
    }
    int e=max(measure,prev);
    return e;
}


void AssignMutationsToBranches (const vector<varcount>& vc_multi, const vector<varcount>& vc_single, vector<nbranch>& treespace) {
    for (int k=0;k<treespace.size();k++) {
        treespace[k].mutations=0;
        if (treespace[k].individuals.size()>1) {
            for (int l=0;l<vc_multi.size();l++) {
                if (vc_multi[l].individual==treespace[k].individuals) {
                    treespace[k].mutations++;
                }
            }
        } else {
            for (int l=0;l<vc_single.size();l++) {
                if (vc_single[l].individual[0]==treespace[k].individuals[0]) {
                    treespace[k].mutations++;
                }
            }
        }
    }
}

void FindEquivalentTimes(const vector<int>& times, const vector<int>& min_times, vector<int>& equiv_times) {
    for (int k=0;k<times.size();k++) {
        equiv_times[k]=times[k]-min_times[k];
    }
}

void CalculateTreeLikelihood (int c, run_params p, double& order_L, const vector<int>& times, const vector<int>& equiv_times, const vector<nbranch>& treespace, const vector<tpairs>& trans_sets, const vector< vector<ijlike> >& new_like_trans) {
    double L=0;
    //cout << "Time likelihood\n";
    //Likelihoods of transmission times
    for (int k=0;k<times.size();k++) {
        L=L+new_like_trans[trans_sets[c].ordered_pair[k].from][trans_sets[c].ordered_pair[k].to].noseq_likes[equiv_times[k]];
        //cout << new_like_trans[trans_sets[c].ordered_pair[k].from][trans_sets[c].ordered_pair[k].to].noseq_likes[equiv_times[k]] << "\n";
    }
    //Calculate sequence-based part of the likelihood
    //cout << "Vals " << p.rate << " " << p.seq_noise << "\n";
    for (int k=0;k<treespace.size();k++) {
        if (treespace[k].individuals.size()==1) {
            L=L+Poisson(treespace[k].mutations,(treespace[k].time*p.rate)+(p.seq_noise/2));
            //cout << k << " N " <<Poisson(treespace[k].mutations,(treespace[k].time*p.rate)+(p.seq_noise/2)) << "\n";
        } else {
            L=L+Poisson(treespace[k].mutations,(treespace[k].time*p.rate));
            //cout << k << " X " << Poisson(treespace[k].mutations,(treespace[k].time*p.rate)) << "\n";
        }
    }
    //Increment likelihood for this order
    order_L=order_L+exp(L);
    //cout << "Likelihood " << L << "\n";
    //cout << "Likelihood is now : " << order_L << "\n";
}

void FindNextTimeSet (int c, int o, const run_params& p, int& done, int& up, int& fin, const vector<int>& min_times, const vector<int>& max_times, const vector<int>& orig_times, const vector<tpairs>& trans_sets, const vector< vector<ijlike> >& new_like_trans, vector<int>& times) {
    while (done==0) {
        if (up<0) {
            done=1;
            fin=1;
        } else {
            if (times[up]>max_times[up]) {
                cout << "Error in FindNextTimeSet\n";
            }
            if (times[up]==max_times[up]) {
                times[up]=orig_times[up];
                up--;
            } else {
                times[up]++;
                //Ignore cases with small likelihoods
                while (new_like_trans[trans_sets[c].ordered_pair[up].from][trans_sets[c].ordered_pair[up].to].noseq_likes[times[up]-min_times[up]]<p.threshold_singleLL&&times[up]<max_times[up]) {
                    //cout << "Add to " << up << " " << times[up] << "\n";
                    times[up]++;
                }
            
                up=times.size()-1;
                CheckTimeOrder(c,o,trans_sets,times);
                done=1;
            }
        }
    }
}

					








void Thresholds (run_params p, double L) {
	if (p.noseq==1) {
		//99.999%:-13.7123
		if (L>-4.62534) {
			cout << "Consistent "; //p>0.05
		} else if (L>-6.17563) {
			cout << "Borderline "; //0.05>p>0.01
		} else {
			cout << "Unlikely "; //0.01>p
		}
	} else {
		//99.9%:-12.8883
		//99.999%:-17.6422
		if (L>-8.15176) {
			cout << "Consistent "; //p>0.05
		} else if (L>-10.1578) {
			cout << "Borderline "; //0.05>p>0.01
		} else {
			cout << "Unlikely "; //0.01>p
		}
	}
	if (L<p.threshold) {
		cout << "* ";
	}
}

