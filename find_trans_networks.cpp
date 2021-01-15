using namespace std;
#include "basicmodel.h"
#include "checkpointing.h"
#include "find_trans_networks.h"

void FindPlausibleTransmissionNetworks (run_params p, const vector< vector<int> >& new_subsets, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, vector< vector<sparseseq> >& variants_sets) {
    ofstream check_file;
    check_file.open("Checkpoint1.out");
    if (p.specify_set>-1) {
        //Run for a specific set
        FPTN_Set (p,p.specify_set,new_subsets,pdat_sets,like_trans_sets,variants_sets,check_file);
    } else {
        //Carry out loop over sets
        for (int set=0;set<pdat_sets.size();set++) {
            FPTN_Set (p,set,new_subsets,pdat_sets,like_trans_sets,variants_sets,check_file);
        }
    }
}

void FPTN_Set (run_params p, int set, const vector< vector<int> >& new_subsets, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, vector< vector<sparseseq> >& variants_sets, ofstream& check_file) {
    cout << "Set " << set << " of " << pdat_sets.size() << "\n";
    int found_complete=0;
    int remove=0;
    if (p.specify_remove>-1) {
        remove=p.specify_remove;
    }
    vector<int> subset;  //Contains an index of which individuals in the subset are being considered.
    for (int j=0;j<new_subsets[set].size();j++) {
        subset.push_back(j);
    }

    //The code proceeds as follows:
    //Look for a complete plausible reconstruction of transmission events.  If this is found it is reported.
    //If no reconstruction is found, systematically remove individuals from the set one by one, and reconstruct again eachtime.
    //If reconstructions are found for multiple sets of individuals of the same size, all of these are reported.
    //The number of individuals removed is incremented up by one each time no reconstruction is found.
    //Where more than one individual is removed, the removed individuals are evaluated for whether they form a plausible network in their own right.
    while (found_complete==0) {
        vector< vector<int> > indiv_sets;
        //Find subsets - containing all but removed individuals
        vector< vector<int> > remainders; //What is left over
        ConstructIndivSets(remove,subset,indiv_sets,remainders);
        int max=indiv_sets.size();
        for (int iset=0;iset<max;iset++) { //Loop over possible sets of individuals
            cout << "Reduced set " << iset << " of " << max << "\n";
            vector<int> new_subset=indiv_sets[iset];
            vector<pat> new_pdat;//Don't actually need this?  Only for printing...
            vector< vector<ijlike> > new_like_trans;
            vector<sparseseq> new_variants;
            SetupParamsSubsetAnalysis (set,new_subset,pdat_sets,like_trans_sets,variants_sets,new_pdat,new_like_trans,new_variants);
            
            //Which individuals are contained?
            cout << "Individuals\n";
            vector<int> indivs;
            for (int j=0;j<new_subset.size();j++) {
                indivs.push_back(new_subset[j]);
                cout << new_subset[j] << " ";
            }
            cout << "\n";

            //Variant information for these individuals...
            if (p.diagnostic==1) {
                PrintVariantInformation (new_pdat,new_variants);
                PrintNLT(new_like_trans);
            }
            
            //Find variants that occur more than once
            vector<varcount> vc_multi;
            vector<varcount> vc_single;
            FindSharedVariantsNew(vc_multi,vc_single,new_variants);
            if (p.diagnostic==1) {
                cout << "FindSharedVariants\n";
                for (int j=0;j<vc_multi.size();j++) {
                    cout << vc_multi[j].locus << vc_multi[j].allele << " ";
                    for (int k=0;k<vc_multi[j].individual.size();k++) {
                        cout << vc_multi[j].individual[k] << " ";
                    }
                    cout << "\n";
                }
                cout << "NonSharedVariants\n";
                for (int j=0;j<vc_single.size();j++) {
                    cout << vc_single[j].locus << vc_single[j].allele << " ";
                    for (int k=0;k<vc_single[j].individual.size();k++) {
                        cout << vc_single[j].individual[k] << " ";
                    }
                    cout << "\n";
                }

            }
           // cin.ignore();
            
            //Identify all plausible sets of transmissions
            //cout << "FindTransSets\n";
            vector<tpairs> trans_sets;
            FindTransSets(p,indivs,new_like_trans,vc_multi,trans_sets);
                
            cout << "TS size " << trans_sets.size() << "\n";
            if (trans_sets.size()>0) {
                found_complete=1;
                WriteCheckpoint1(set,iset,new_subset,trans_sets);
                //Write to record of checkpoint files
                check_file << set << " " << iset << " 0\n";
                //cout << "Here " << iset << "\n";
                if (remainders.size()>1) {
                    //cout << "Remainder " << remainders[iset].size() << "\n";
                    if (remainders[iset].size()>1) {
                        //Deal with remainder set
                        ProcessRemainder(set,iset,p,check_file,remainders,pdat_sets,like_trans_sets,variants_sets);
                    }
                }
            }
        }
        if (p.specify_remove>-1) {
            break;
        } else {
            remove++;
        }
    }
}

void ConstructIndivSets (const int remove, const vector<int>& indivs, vector< vector<int> >& indiv_sets, vector< vector<int> >& remainders) {
    if (remove==0) {
        indiv_sets.push_back(indivs);
    }
    if (remove==1) {
        for (int a1=0;a1<indivs.size();a1++) {
            vector<int> ind_set;
            vector<int> rem;
            rem.push_back(a1);
            for (int i=0;i<indivs.size();i++) {
                if (i!=a1) {
                    ind_set.push_back(indivs[i]);
                }
            }
            indiv_sets.push_back(ind_set);
            remainders.push_back(rem);
        }
    }
    if (remove==2) {
        for (int a1=0;a1<indivs.size()-1;a1++) {
            for (int a2=a1+1;a2<indivs.size();a2++) {
                vector<int> ind_set;
                vector<int> rem;
                rem.push_back(a1);
                rem.push_back(a2);
                for (int i=0;i<indivs.size();i++) {
                    if (i!=a1&&i!=a2) {
                        ind_set.push_back(indivs[i]);
                    }
                }
                indiv_sets.push_back(ind_set);
                remainders.push_back(rem);
            }
        }
    }
    if (remove==3) {
        for (int a1=0;a1<indivs.size()-2;a1++) {
            for (int a2=a1+1;a2<indivs.size()-1;a2++) {
                for (int a3=a2+1;a3<indivs.size();a3++) {
                    vector<int> ind_set;
                    vector<int> rem;
                    rem.push_back(a1);
                    rem.push_back(a2);
                    rem.push_back(a3);
                    for (int i=0;i<indivs.size();i++) {
                        if (i!=a1&&i!=a2&&i!=a3) {
                            ind_set.push_back(indivs[i]);
                        }
                    }
                    indiv_sets.push_back(ind_set);
                    remainders.push_back(rem);
                }
            }
        }
    }
    if (remove==4) {
        for (int a1=0;a1<indivs.size()-3;a1++) {
            for (int a2=a1+1;a2<indivs.size()-2;a2++) {
                for (int a3=a2+1;a3<indivs.size()-1;a3++) {
                    for (int a4=a3+1;a4<indivs.size();a4++) {
                        vector<int> ind_set;
                        vector<int> rem;
                        rem.push_back(a1);
                        rem.push_back(a2);
                        rem.push_back(a3);
                        rem.push_back(a4);
                        for (int i=0;i<indivs.size();i++) {
                            if (i!=a1&&i!=a2&&i!=a3&&i!=a4) {
                                ind_set.push_back(indivs[i]);
                            }
                        }
                        indiv_sets.push_back(ind_set);
                        remainders.push_back(rem);
}
                }
            }
        }
    }
    if (remove==5) {
        for (int a1=0;a1<indivs.size()-4;a1++) {
            for (int a2=a1+1;a2<indivs.size()-3;a2++) {
                for (int a3=a2+1;a3<indivs.size()-2;a3++) {
                    for (int a4=a3+1;a4<indivs.size()-1;a4++) {
                        for (int a5=a4+1;a5<indivs.size();a5++) {
                            vector<int> ind_set;
                            vector<int> rem;
                            rem.push_back(a1);
                            rem.push_back(a2);
                            rem.push_back(a3);
                            rem.push_back(a4);
                            rem.push_back(a5);
                            for (int i=0;i<indivs.size();i++) {
                                if (i!=a1&&i!=a2&&i!=a3&&i!=a4&&i!=a5) {
                                    ind_set.push_back(indivs[i]);
                                }
                            }
                            indiv_sets.push_back(ind_set);
                            remainders.push_back(rem);
                        }
                    }
                }
            }
        }
    }
    if (remove==6) {
        for (int a1=0;a1<indivs.size()-4;a1++) {
            for (int a2=a1+1;a2<indivs.size()-3;a2++) {
                for (int a3=a2+1;a3<indivs.size()-2;a3++) {
                    for (int a4=a3+1;a4<indivs.size()-1;a4++) {
                        for (int a5=a4+1;a5<indivs.size();a5++) {
                            for (int a6=a5+1;a6<indivs.size();a6++) {
                                vector<int> ind_set;
                                vector<int> rem;
                                rem.push_back(a1);
                                rem.push_back(a2);
                                rem.push_back(a3);
                                rem.push_back(a4);
                                rem.push_back(a5);
                                rem.push_back(a6);
                                for (int i=0;i<indivs.size();i++) {
                                    if (i!=a1&&i!=a2&&i!=a3&&i!=a4&&i!=a5&&i!=a6) {
                                        ind_set.push_back(indivs[i]);
                                    }
                                }
                                indiv_sets.push_back(ind_set);
                                remainders.push_back(rem);
                            }
                        }
                    }
                }
            }
        }
    }
}

void SetupParamsSubsetAnalysis (int set, const vector<int>& new_subset, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, const vector< vector<sparseseq> >& variants_sets, vector<pat>& new_pdat, vector< vector<ijlike> >& new_like_trans, vector<sparseseq>& new_variants) {
    for (int j=0;j<new_subset.size();j++) {
        //cout << "Set " << set << " " << new_subset[j] << "\n";
        new_pdat.push_back(pdat_sets[set][new_subset[j]]);
        //Problem with likelihoods - need to do double reconstruction
        vector<ijlike> nlt;
        //cout << like_trans_sets.size() << " " << like_trans_sets[0].size() << "\n";
        for (int k=0;k<new_subset.size();k++) {
            //cout << "k " << k << " " << new_subset.size() << "\n";
            nlt.push_back(like_trans_sets[set][new_subset[j]][new_subset[k]]);
        }
        new_like_trans.push_back(nlt);
        new_variants.push_back(variants_sets[set][new_subset[j]]);
    }
}

void PrintVariantInformation (const vector<pat>& new_pdat, const vector<sparseseq>& new_variants) {
    cout << "Variants\n";
    for (int i=0;i<new_variants.size();i++) {
        cout << i << " " << new_pdat[i].code << " " << new_pdat[i].time_s << " ";
        for (int j=0;j<new_variants[i].locus.size();j++) {
            cout << new_variants[i].locus[j] << new_variants[i].allele[j] << " ";
        }
        cout << "\n";
    }
}

void PrintNLT (const vector< vector<ijlike> >& new_like_trans) {
    cout << "Likelihoods\n";
    for (int j=0;j<new_like_trans.size();j++) {
        for (int k=0;k<new_like_trans[j].size();k++) {
            if (j!=k) {
                cout << j << " " << k << " " << new_like_trans[j][k].lL_tot << " " << new_like_trans[j][k].ns_lL_tot << " " <<  new_like_trans[j][k].min << " " << new_like_trans[j][k].max << "\n";
            }
        }
    }
}

void FindSharedVariantsNew(vector<varcount>& vc_multi, vector<varcount>& vc_single, const vector<sparseseq>& new_variants) {
    //Finds shared variants.  Clusters variants by sets of shared individuals
    vector<varcount> vc;
    for (int i=0;i<new_variants.size();i++) {
        //cout << "i " << i << " " << new_variants.size() << "\n";
        for (int j=0;j<new_variants[i].locus.size();j++) {
            int found=0;
            for (int k=0;k<vc.size();k++) {
                if (new_variants[i].locus[j]==vc[k].locus&&new_variants[i].allele[j]==vc[k].allele) {
                    vc[k].individual.push_back(i);
                    found=1;
                }
            }
            if (found==0) {
                //cout << "Add " << new_variants[i].locus[j] << " " << new_variants[i].allele[j] << "\n";
                varcount v;
                v.locus=new_variants[i].locus[j];
                v.allele=new_variants[i].allele[j];
                v.individual.push_back(i);
                vc.push_back(v);
            }
        }
    }

    for (int i=0;i<vc.size();i++) {
        //cout << "i " << i << " " << vc.size() << "\n";
        if (vc[i].individual.size()>1) {
            vc_multi.push_back(vc[i]);
        }
        if (vc[i].individual.size()==1) {
            vc_single.push_back(vc[i]);
        }
    }
    //Delete duplicate sets by which indiviudals they contain
    vector<int> to_rem;
    if (vc_multi.size()>0) {
        for (int i=0;i<vc_multi.size()-1;i++) {
            //Remove variants found in all individuals?
            if (vc_multi[i].individual.size()==new_variants.size()) {
                to_rem.push_back(i);
            }
            for (int j=i+1;j<vc_multi.size();j++) {
                //Remove variants with a common set of individuals
                if (vc_multi[i].individual==vc_multi[j].individual) {
                    to_rem.push_back(j);
                }
            }
        }
    }
    sort(to_rem.begin(),to_rem.end());
    to_rem.erase(unique(to_rem.begin(),to_rem.end()),to_rem.end());
    reverse(to_rem.begin(),to_rem.end());
    for (int i=0;i<to_rem.size();i++){
        vc_multi.erase(vc_multi.begin()+to_rem[i]);
    }
    /*cout << "Final variant calls\n";
    for (int i=0;i<vc.size();i++) {
        cout << vc[i].locus << vc[i].allele <<" ";
        for (int j=0;j<vc[i].individual.size();j++) {
            cout << vc[i].individual[j] << " ";
        }
        cout << "\n";
    }*/
}

void FindTransSets(run_params p, const vector<int>& indivs, const vector< vector<ijlike> >& new_like_trans, const vector<varcount>& vc_multi, vector<tpairs>& trans_sets) {
    //Identify all plausible sets of transmissions
    FindPairSets(p,indivs,new_like_trans,vc_multi,trans_sets);
    //Import timings
    ImportTimings(new_like_trans,trans_sets);
    //Identify constraints on the orderings of transmission events in each case from intrinsic events and timings
    FindOrderingConstraints(trans_sets);
    //Check for contradictions in this ordering
    CheckOrderingContradictions (trans_sets);
}

void FindPairSets (run_params p, const vector<int>& indivs, const vector< vector<ijlike> >& like_trans, const vector<varcount>& vc_multi, vector<tpairs>& trans_sets) {
    //cout << "FindPairSets\n";
    for (int i=0;i<indivs.size();i++) {
        //cout << "I= " << i << "\n";
        vector< vector<tpair> > allsets;
        for (int j=0;j<indivs.size()-1;j++) {
            vector<tpair> tp;
            allsets.push_back(tp);
        }
        int root=indivs[i];
        //cout << "Root is " << indivs[i] << "\n";
        
        //Make list of all possible transmissions conditional on the root chosen.  Reduces the number of sets
        for (int j=0;j<indivs.size();j++) {
            for (int k=0;k<indivs.size();k++) {
                if (j!=k&&indivs[k]!=root) {
                    tpair t;
                    t.from=j;
                    t.to=k;
                    if (k<i) {
                        allsets[k].push_back(t);
                    } else {
                        allsets[k-1].push_back(t);
                    }
                }
            }
        }
        vector<int> index (allsets.size(),0);
        int fin=0;
        while (fin==0) {
            int check_root=1;
            int check_cycle=1;
            int check_unlike=1;
            int check_span=0;
            int check_seq=1;
            //Must contain transmission from root to something else.  Can't contain transmission to the root from something else
            if (check_root==1) {
                //Must not contain unlikely events according to the non-sequence threshold
                CheckLikelihood(check_unlike,p,index,allsets,like_trans);
            }

            if (check_root==1&&check_unlike==1) {
                //Check that edges span all individuals
                CheckSpan (check_span,index,allsets);
            }
            
            if (check_root==1&&check_span==1&&check_unlike==1) {
                //Check there are no cycles
                CheckCycle (check_cycle,index,allsets);
            }
            
            
            if (check_root==1&&check_cycle==1&&check_span==1&&check_unlike==1) {
                //Check for phylogenetic consistency using the shared sequence variants
                CheckSharedVariants(check_seq,index,allsets,vc_multi);
            }
            
            if (check_root==1&&check_cycle==1&&check_span==1&&check_unlike==1&&check_seq==1) {  //Keep this set of transmissions
                tpairs ts;
                ts.root=i;
                //cout << "K ";
                cout << i << " " << index.size() << " ";
                for (int j=0;j<index.size();j++) {
                    cout << allsets[j][index[j]].from << "-" << allsets[j][index[j]].to << ", ";
                    ts.pair.push_back(allsets[j][index[j]]);
                }
                //Begin ordering calculation: Look here at orderings imposed by shared mutations.  If the first individual to report a variant transmits to set A of individuals without the variant and to set B of individuals with the variant
                OrderBySharedVariants(index,allsets,vc_multi,ts);
                trans_sets.push_back(ts);
                cout << "\n";
            }
            int done=0;
            int up=index.size()-1;
            while (done==0) {
                if (index[up]==allsets.size()-1) {
                    index[up]=0;
                    up--;
                    if (up<0) {
                        fin=1;
                        up=0;
                        done=1;
                    }
                } else {
                    index[up]++;
                    up=index.size()-1;
                    done=1;
                }
            }
        }
    }
}

void CheckLikelihood(int& check_unlike, run_params p, const vector<int>& index, const vector< vector<tpair> >& allsets, const  vector< vector<ijlike> >& like_trans) {
    for (int j=0;j<index.size();j++) {
        if (p.consistency==1) {
            if (like_trans[allsets[j][index[j]].from][allsets[j][index[j]].to].lL_tot<-8.15176) {
                check_unlike=0;
                break;
            }
        }
        if (p.consistency==2) {
            if (like_trans[allsets[j][index[j]].from][allsets[j][index[j]].to].lL_tot<-10.1578) {
                check_unlike=0;
                break;
            }
        }
    }
}

void CheckSpan (int& check_span, const vector<int>& index, const vector< vector<tpair> >& allsets) {
    vector<int> contains;
    contains.push_back(allsets[0][index[0]].to);
    contains.push_back(allsets[0][index[0]].from);
    int add=1;
    while (add==1) {
        add=0;
        int s=contains.size();
        for (int j=1;j<index.size();j++) {
            for (int k=0;k<s;k++) {
                if (allsets[j][index[j]].to==contains[k]) {
                    contains.push_back(allsets[j][index[j]].from);
                }
                if (allsets[j][index[j]].from==contains[k]) {
                    contains.push_back(allsets[j][index[j]].to);
                }
            }
        }
        sort(contains.begin(),contains.end());
        contains.erase(unique(contains.begin(),contains.end()),contains.end());
        if (contains.size()>s) {
            add=1;
        }
    }
    if (contains.size()==index.size()+1) {
        check_span=1;
    }
}

void CheckCycle (int& check_cycle, const vector<int>& index, const vector< vector<tpair> >& allsets) {
    vector<tpair> tpset;
    for (int j=0;j<index.size();j++) {
        tpair t;
        t.from=allsets[j][index[j]].from;
        t.to=allsets[j][index[j]].to;
        tpset.push_back(t);
    }
    int findc=TransFindCycle(tpset);
    if (findc==1) {
        check_cycle=0;
    }
}

int TransFindCycle (vector<tpair> tpset) {
    vector<int> path;
    int loop=0;
    for (int i=0;i<tpset.size();i++) {
        if (loop==0) {
            //cout << "i " << tpset[i].from << "-" << tpset[i].to << "\n";
            path.clear();
            path.push_back(tpset[i].from);
            int next=tpset[i].to;
            loop=0;
            path.push_back(next);
            int add=1;
            while (add==1) {
                add=0;
                int size=path.size();
                //Did this just complete a cycle?
                if (path[path.size()-1]==path[0]) {
                    loop=1;
                    //cout << "Found\n";
                    break;
                }
                //Did this just reach an end?
                int atend=1;
                int iadd=-1;
                for (int j=0;j<tpset.size();j++) {
                    if (tpset[j].from==path[path.size()-1]) {
                        atend=0;
                        iadd=j;
                    }
                }
                if (atend==1) {
                    loop=0;
                    break;
                } else {
                    path.push_back(tpset[iadd].to);
                    //cout << tpset[iadd].from << "-" << tpset[iadd].to << "\n";
                }
                //Did we just add something?
                if (path.size()>size) {
                    add=1;
                }
            }
        }
    }

    if (loop==0) {
        //cout << "No loop\n";
        return 0;
    } else {
        //cout << "Loop\n";
        return 1;
    }
}


void CheckSharedVariants(int& check_seq, const vector<int>& index, const vector< vector<tpair> >& allsets, const vector<varcount>& vc_multi) {
    for (int v=0;v<vc_multi.size();v++) {
        //Not more than one transmission into the set
        int into_set=0;
        int from_ind=-1;
        for (int j=0;j<index.size();j++) {
            int to_in=0;
            int from_in=0;
            for (int k=0;k<vc_multi[v].individual.size();k++) {
                if (allsets[j][index[j]].to==vc_multi[v].individual[k]) {
                    to_in=1;
                }
                if (allsets[j][index[j]].from==vc_multi[v].individual[k]) {
                    from_in=1;
                }
            }
            if (to_in==1&&from_in==0) {
                //cout << "Into set " << allsets[j][index[j]].from << " " << allsets[j][index[j]].to << "\n";
                into_set++;
                from_ind=allsets[j][index[j]].to;
            }
        }
        //cout << "Number in " << into_set << " ";
        if (into_set>1) {
            check_seq=0;
            //cout << "\n";
            break;
        } else {
            //cout << from_ind << "\n";
            //Apart from the first individual to get the virus, no transmission to anything outside the set
            for (int j=0;j<index.size();j++) {
                int to_in=0;
                int from_in=0;
                for (int k=0;k<vc_multi[v].individual.size();k++) {
                    if (allsets[j][index[j]].to==vc_multi[v].individual[k]) {
                        to_in=1;
                    }
                    if (allsets[j][index[j]].from==vc_multi[v].individual[k]) {
                        from_in=1;
                    }
                }
            
                if (from_in==1&&to_in==0&&allsets[j][index[j]].from!=from_ind) {
                    //Edge goes from something inside to something outside and is not from the original individual.
                    //cout << "Illegal step out " << allsets[j][index[j]].from << " " << allsets[j][index[j]].to << "\n";
                    check_seq=0;
                    break;
                }
                //Note here.  If it is from the original individual we want to know so as to add in constraints
                //Would imply from_in==1, to_in==0, all_sets[j][index[j]].from==from_ind.
                //Add this in at a later point
                
            }
        }
    }
}

void OrderBySharedVariants(const vector<int>& index, const vector< vector<tpair> >& allsets, const vector<varcount>& vc_multi, tpairs& ts) {
    vector<int> two;
    two.push_back(0);
    two.push_back(0);
    for (int v=0;v<vc_multi.size();v++) {
        //Count transmissions into the set
        int into_set=0;
        int from_ind=-1;
        for (int j=0;j<index.size();j++) {
            int to_in=0;
            int from_in=0;
            for (int k=0;k<vc_multi[v].individual.size();k++) {
                if (allsets[j][index[j]].to==vc_multi[v].individual[k]) {
                    to_in=1;
                }
                if (allsets[j][index[j]].from==vc_multi[v].individual[k]) {
                    from_in=1;
                }
            }
            if (to_in==1&&from_in==0) {
                into_set++;
                from_ind=allsets[j][index[j]].to;
            }
        }
        if (into_set==1) {
            //Find edges from the first individual to get the virus in the set
            vector<int> before;
            vector<int> after;
            //Every transmission to something with the mutation must be after everything before
            for (int j=0;j<index.size();j++) {
                if (allsets[j][index[j]].from!=from_ind) {
                    int to_in=0;
                    for (int k=0;k<vc_multi[v].individual.size();k++) {
                        if (allsets[j][index[j]].to==vc_multi[v].individual[k]) {
                            to_in=1;
                        }
                    }
                    if (to_in==1) {
                        after.push_back(j);
                    } else {
                        before.push_back(j);
                    }
                }
            }
            for (int j=0;j<before.size();j++) {
                for (int k=0;k<after.size();k++) {
                    two[0]=before[j];
                    two[1]=after[k];
                    ts.constraints.push_back(two);
                }
            }
        }
    }

}

void ImportTimings(const vector< vector<ijlike> >& like_trans, vector<tpairs>& trans_sets) {
    for (int i=0;i<trans_sets.size();i++) {
        for (int j=0;j<trans_sets[i].pair.size();j++) {
            trans_sets[i].min_time.push_back(like_trans[trans_sets[i].pair[j].from][trans_sets[i].pair[j].to].min);
            trans_sets[i].max_time.push_back(like_trans[trans_sets[i].pair[j].from][trans_sets[i].pair[j].to].max);
        }
    }
}

void FindOrderingConstraints(vector<tpairs>& trans_sets) {
    for (int i=0;i<trans_sets.size();i++) {
        /*cout << "Orders for " << i << "\n";
        for (int j=0;j<trans_sets[i].pair.size();j++) {
            cout << j << " " << trans_sets[i].pair[j].from << "-" << trans_sets[i].pair[j].to << " " << trans_sets[i].min_time[j] << " " << trans_sets[i].max_time[j] << "\n";
        }*/
        vector<int> two;
        two.push_back(0);
        two.push_back(0);
        for (int j=0;j<trans_sets[i].pair.size();j++) {
            for (int k=0;k<trans_sets[i].pair.size();k++) {
                if (j!=k) {
                    //Ordering related to individuals e.g. 2-3 after 1-2
                    if (trans_sets[i].pair[j].to==trans_sets[i].pair[k].from) {
                        //cout << j << " " << trans_sets[i].pair[j].from << "-" << trans_sets[i].pair[j].to << " before " << k << " " << trans_sets[i].pair[k].from << "-" << trans_sets[i].pair[k].to << "\n";
                        two[0]=j;
                        two[1]=k;
                        trans_sets[i].constraints.push_back(two);
                        if (trans_sets[i].min_time[j]>=trans_sets[i].min_time[k]) {
                            trans_sets[i].min_time[k]=trans_sets[i].min_time[j]+1;
                            //cout << "Change min time of " << trans_sets[i].pair[k].from << "-" << trans_sets[i].pair[k].to << " to " << trans_sets[i].min_time[j]+1 << "\n";
                        }
                    }
                }
            }
        }
        int change=1; //Cycle through successive events - keep going until convergence
        while (change==1) {
            change=0;
            for (int j=0;j<trans_sets[i].pair.size();j++) {
                for (int k=0;k<trans_sets[i].pair.size();k++) {
                    if (j!=k) {
                        if (trans_sets[i].pair[j].to==trans_sets[i].pair[k].from) {
                            if (trans_sets[i].min_time[j]>=trans_sets[i].min_time[k]) {
                                trans_sets[i].min_time[k]=trans_sets[i].min_time[j]+1;
                                change=1;
                            }
                        }
                    }
                }
            }
        }
        
        /*cout << "Things now stand as:\n";
        for (int j=0;j<trans_sets[i].pair.size();j++) {
            cout << j << " " << trans_sets[i].pair[j].from << "-" << trans_sets[i].pair[j].to << " " << trans_sets[i].min_time[j] << " " << trans_sets[i].max_time[j] << "\n";
        }*/
        for (int j=0;j<trans_sets[i].pair.size();j++) {
            for (int k=0;k<trans_sets[i].pair.size();k++) {
                if (j!=k) {
                    //Orderings related to time
                    if (trans_sets[i].min_time[j]>trans_sets[i].max_time[k]) {
                        //cout << k << " " << trans_sets[i].pair[k].from << "-" << trans_sets[i].pair[k].to << " " <<  trans_sets[i].min_time[k] << " " << trans_sets[i].max_time[k] << " before " << j << " " <<  trans_sets[i].pair[j].from << "-" << trans_sets[i].pair[j].to << " " << trans_sets[i].min_time[j] << " " << trans_sets[i].max_time[j] << "\n";
                        two[0]=k;
                        two[1]=j;
                        trans_sets[i].constraints.push_back(two);
                    }
                }
            }
        }
    }
}

void CheckOrderingContradictions (vector<tpairs>& trans_sets) {
    vector<tpairs> trans_sets_new;
    for (int i=0;i<trans_sets.size();i++) {
        /*for (int j=0;j<trans_sets[i].pair.size();j++) {
            cout << trans_sets[i].pair[j].from << "-" << trans_sets[i].pair[j].to << " " << trans_sets[i].min_time[j] << " " << trans_sets[i].max_time[j] << "\n";
        }
        for (int j=0;j<trans_sets[i].constraints.size();j++) {
            cout << trans_sets[i].constraints[j][0] << " " << trans_sets[i].constraints[j][1] << "\n";
        }*/

        trans_sets[i].const_contradiction=0;
        for (int j=0;j<trans_sets[i].min_time.size();j++) {
            if (trans_sets[i].min_time[j]>trans_sets[i].max_time[j]) {
                trans_sets[i].const_contradiction=1;
                break;
            }
        }
        if (trans_sets[i].const_contradiction==0) {
            for (int j=0;j<trans_sets[i].constraints.size();j++) {
                if (trans_sets[i].const_contradiction==0) {
                    for (int k=0;k<trans_sets[i].constraints.size();k++) {
                        if (j!=k) {
                            if (trans_sets[i].constraints[j][0]==trans_sets[i].constraints[k][1]&&trans_sets[i].constraints[j][1]==trans_sets[i].constraints[k][0]) {
                                trans_sets[i].const_contradiction=1;
                                //cout << "Set " << i << " found contradiction " << trans_sets[i].constraints[j][0] << " " << trans_sets[i].constraints[j][1] << "\n";
                                break;
                            }
                        }
                    }
                }
            }
        }
        if (trans_sets[i].const_contradiction==0) {
            trans_sets_new.push_back(trans_sets[i]);
        }
    }
    trans_sets=trans_sets_new;
}

void ProcessRemainder(int set, int iset, run_params p, ofstream& check_file, const vector< vector<int> >& remainders, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, const vector< vector<sparseseq> >& variants_sets) {
    //Set up parameters
    vector<int> new_subset_rem=remainders[iset];
    vector<pat> new_pdat_rem;//Don't actually need this?
    vector< vector<ijlike> > new_like_trans_rem;
    vector<sparseseq> new_variants_rem;
    SetupParamsSubsetAnalysis (set,new_subset_rem,pdat_sets,like_trans_sets,variants_sets,new_pdat_rem,new_like_trans_rem,new_variants_rem);
    
    vector<int> indivs_rem;
    for (int j=0;j<new_subset_rem.size();j++) {
        indivs_rem.push_back(new_subset_rem[j]);
    }

    //Generate transmissions
    //Find variants that occur more than once
    vector<varcount> vc_multi_rem;
    vector<varcount> vc_single_rem;
    FindSharedVariantsNew(vc_multi_rem,vc_single_rem,new_variants_rem);
    
    //Identify all plausible sets of transmissions
    vector<tpairs> trans_sets_rem;
    FindTransSets(p,indivs_rem,new_like_trans_rem,vc_multi_rem,trans_sets_rem);
    
    if (trans_sets_rem.size()>0) {
        WriteCheckpoint_Rem1(set,iset,new_subset_rem,trans_sets_rem);
        check_file << set << " " << iset << " 1\n";
    }
}

