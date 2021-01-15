using namespace std;
#include "basicmodel.h"
#include "diagnostics.h"
#include "variants.h"

void ProcessVariantInformation (run_params p, const vector<string>& seqs, vector<sparseseq>& variants, vector<pat>& pdat) {
    string all_consensus;
    FindConsensus(all_consensus,seqs);
    //Variants with respect to general consensus.
    FindVariants (variants,all_consensus,pdat);
    //Find sequences with an 'N' at sites with called variants
    vector<allele> allvar;
    ListAllVariantPositions(variants,allvar);
    if (p.diagnostic==1) {
        PrintVariantLoci(allvar);
    }
    //Next:
    //a) Find internal allele frequencies at loci for which there is at least one N
    //b) Remove sequences with too many uncertain loci in variant positions
    //c) Discard unneeded frequencies
    vector<int> nloci;
    vector<int> nloc_count;
    FindAmbiguousVarPositions (allvar,pdat,nloci,nloc_count);
    vector<ipair> freqs;
    FindAmbiguousVarFreqs (all_consensus,allvar,variants,pdat,nloci,freqs);
    ProcessAmbiguousFixedVariants(freqs,nloc_count,variants,pdat);
    //Delete sequences with too many Ns
    RemoveIndividualsMultipleN(p,nloc_count,pdat,variants);
    if (p.noseq==0&&p.diagnostic==1) {
        PrintVariants(variants,pdat);
    }
}

void FindConsensus (string& consensus, const vector<string>& seqs) {
    consensus=seqs[0];
    cout << "Find consensus of all input sequences\n";
    int nA=0;
    int nC=0;
    int nG=0;
    int nT=0;
    for (int pos=0;pos<seqs[0].size();pos++) {
        nA=0;
        nC=0;
        nG=0;
        nT=0;
        for (int seq=0;seq<seqs.size();seq++) {
            if (seqs[seq][pos]=='A') {
                nA++;
            }
            if (seqs[seq][pos]=='C') {
                nC++;
            }
            if (seqs[seq][pos]=='G') {
                nG++;
            }
            if (seqs[seq][pos]=='T') {
                nT++;
            }
        }
        int max=nA;
        consensus[pos]='A';
        if (nC>max) {
            max=nC;
            consensus[pos]='C';
        }
        if (nG>max) {
            max=nG;
            consensus[pos]='G';
        }
        if (nT>max) {
            consensus[pos]='T';
            max=nT;
        }
        if (max==0) {
            consensus[pos]='-';
        }
    }
}

void FindVariants (vector<sparseseq>& variants, string& consensus, vector<pat>& pdat) {
	for (int i=0;i<pdat.size();i++) {
        sparseseq s;
        for (int pos=0;pos<pdat[i].seq.size();pos++) {
			if (pdat[i].seq.compare(pos,1,consensus,pos,1)!=0) {
				if (pdat[i].seq.compare(pos,1,"A")==0||pdat[i].seq.compare(pos,1,"C")==0||pdat[i].seq.compare(pos,1,"G")==0||pdat[i].seq.compare(pos,1,"T")==0) {
					//cout << "Found variant " << pdat[i].code_match << " " << pos << " " << consensus[pos] << " " << pdat[i].seq[pos] << "\n";
					if (pos!=28827&&pos!=28828) {
						s.locus.push_back(pos);
						s.allele.push_back(pdat[i].seq[pos]);
					}
                }
            }
        }
        variants.push_back(s);
    }
}

bool compare_allele(allele v1, allele v2) {
    return (v1.loc < v2.loc);
}

void ListAllVariantPositions (const vector<sparseseq>& variants, vector<allele>& allvar) {
    //Put into the vector allvar
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            allele a;
            a.loc=variants[i].locus[j];
            a.nuc=variants[i].allele[j];
            allvar.push_back(a);
        }
    }
    sort(allvar.begin(),allvar.end(),compare_allele);
    vector<int> to_rem;
    for (int i=allvar.size()-1;i>0;i--) {
        if (allvar[i].loc==allvar[i-1].loc&&allvar[i].nuc==allvar[i-1].nuc) {
            to_rem.push_back(i);
            //cout << "Remove " << i << "\n";
        }
    }
    for (int i=0;i<to_rem.size();i++) {
        allvar.erase(allvar.begin()+to_rem[i]);
    }
}

void FindAmbiguousVarPositions (const vector<allele>& allvar, vector<pat>& pdat, vector<int>& nloci, vector<int>& nloc_count) {
    //Variant positions for which at least one sequence has an N
    for (int i=0;i<pdat.size();i++) {
        int nc=0;
        for (int j=0;j<allvar.size();j++) {
            if (pdat[i].seq.compare(allvar[j].loc,1,"N")==0) {
                idat id;
                id.loc=allvar[j].loc;
                id.q=0;
                pdat[i].seq_uncertainty.push_back(id);
                nc++;
                nloci.push_back(allvar[j].loc);
            }
        }
        nloc_count.push_back(nc);
    }
    sort(nloci.begin(),nloci.end());
    nloci.erase(unique(nloci.begin(),nloci.end()),nloci.end());
}

void FindAmbiguousVarFreqs (const string all_consensus, const vector<allele>& allvar, const vector<sparseseq>& variants, vector<pat>& pdat, const vector<int>& nloci, vector<ipair>& freqs) {
    //Variant allele frequencies from the population of the variants for which Ns exist
    
    //Find consensus and variant loci at each position
    vector<char> allnuc;
    vector<char> cons;
    for (int i=0;i<nloci.size();i++) {
        int done=0;
        cons.push_back(all_consensus[nloci[i]]);
        for (int j=0;j<pdat.size();j++) {
            if (done==0) {
                for (int k=0;k<variants[j].locus.size();k++) {
                    if (variants[j].locus[k]==nloci[i]) {
                        allnuc.push_back(variants[j].allele[k]);
                        done=1;
                        break;
                    }
                }
            }
        }
    }

    for (int i=0;i<nloci.size();i++) {
        ipair ip;
        ip.loc=nloci[i];
        double countc=0;
        double countv=0;
        for (int j=0;j<pdat.size();j++) {
            string c(1,cons[i]);
            string v(1,allnuc[i]);
            if (pdat[j].seq.compare(nloci[i],1,c,0,1)==0) {
                countc++;
            }
            if (pdat[j].seq.compare(nloci[i],1,v,0,1)==0) {
                countv++;
            }
        }
        ip.q=countv/(countc+countv);
        freqs.push_back(ip);
    }
    
    //Add variant frequencies and alleles into pdat
    for (int i=0;i<pdat.size();i++) {
        for (int j=0;j<pdat[i].seq_uncertainty.size();j++) {
            for (int k=0;k<freqs.size();k++) {
                if (freqs[k].loc==pdat[i].seq_uncertainty[j].loc) {
                    pdat[i].seq_uncertainty[j].q=freqs[k].q;
                    pdat[i].seq_uncertainty[j].allele=allnuc[k];
                }
            }
            
        }
    }
}

void ProcessAmbiguousFixedVariants(vector<ipair>& freqs, vector<int>& nloc_count, vector<sparseseq>& variants, vector<pat>& pdat) {
    //Repair uncertain variants that are seen in all other sequences i.e. frequency 1
    for (int i=freqs.size()-1;i>=0;i--) {
        if (freqs[i].q==1) {
            //Remove this as an uncertain site from each sequence
            for (int j=0;j<pdat.size();j++) {
                for (int k=0;k<pdat[j].seq_uncertainty.size();k++) {
                    if (pdat[j].seq_uncertainty[k].loc==freqs[i].loc) {
                        pdat[j].seq_uncertainty.erase(pdat[j].seq_uncertainty.begin()+k);
                        nloc_count[j]--;
                        break;
                    }
                }
            }
            //Add this to the list of variants
            char v;
            for (int j=0;j<variants.size();j++) {
                for (int k=0;k<variants[j].locus.size();k++) {
                    if (variants[j].locus[k]==freqs[i].loc) {
                        v=variants[j].allele[k];
                        break;
                    }
                }
            }
            //cout << "Here " << v << "\n";
            for (int j=0;j<variants.size();j++) {
                for (int k=0;k<=variants[j].locus.size();k++) {
                    if (k==0&&variants[j].locus[k]>freqs[i].loc) {
                        //cout << "Site 1\n";
                        variants[j].locus.insert(variants[j].locus.begin(),freqs[i].loc);
                        variants[j].allele.insert(variants[j].allele.begin(),v);
                        break;
                    }
                    if (k>=0&&k<variants[j].locus.size()) {
                        if (variants[j].locus[k]<freqs[i].loc&&variants[j].locus[k+1]>freqs[i].loc) {
                            //cout << "Site 2\n";
                            variants[j].locus.insert(variants[j].locus.begin()+k+1,freqs[i].loc);
                            variants[j].allele.insert(variants[j].allele.begin()+k+1,v);
                            break;

                        }
                    }
                    if (k==variants[j].locus.size()&&variants[j].locus[k-1]<freqs[i].loc) {
                        //cout << "Site 3 " << k << " " << variants[j].locus.size() << "\n";
                        variants[j].locus.insert(variants[j].locus.begin()+k,freqs[i].loc);
                        variants[j].allele.insert(variants[j].allele.begin()+k,v);
                        break;
                    }
                }
            }
            //Remove this from freqs as now dealt with
            freqs.erase(freqs.begin()+i);
        }
    }
}

void RemoveIndividualsMultipleN (const run_params p, vector<int>& nloc_count, vector<pat>& pdat, vector<sparseseq>& variants) {
    vector<int> rem;
    for (int i=pdat.size()-1;i>=0;i--) {
        if (nloc_count[i]>p.max_n) {
            rem.push_back(i);
        }
    }
    for (int i=0;i<rem.size();i++) {
        pdat.erase(pdat.begin()+rem[i]);
        variants.erase(variants.begin()+rem[i]);
        nloc_count.erase(nloc_count.begin()+rem[i]);
    }
}

void FindPairwiseDistances (run_params p, vector< vector<int> >& seqdists, vector<sparseseq>& variants, vector<pat>& pdat) {
    vector<int> zeros(pdat.size(),0);
    for (int i=0;i<pdat.size();i++) {
        seqdists.push_back(zeros);
    }
    for (int i=0;i<pdat.size();i++) {
        for (int j=i+1;j<pdat.size();j++) {
            int dist=0;
            //Find unique difference positions;
            vector<int> uniq;
            for (int k=0;k<variants[i].locus.size();k++) {
                uniq.push_back(variants[i].locus[k]);
            }
            for (int k=0;k<variants[j].locus.size();k++) {
                uniq.push_back(variants[j].locus[k]);
            }
            sort(uniq.begin(),uniq.end());
            uniq.erase(unique(uniq.begin(),uniq.end()),uniq.end());
            for (int k=0;k<uniq.size();k++) {
                if (pdat[i].seq[uniq[k]]=='A'||pdat[i].seq[uniq[k]]=='C'||pdat[i].seq[uniq[k]]=='G'||pdat[i].seq[uniq[k]]=='T') {
                    if (pdat[j].seq[uniq[k]]=='A'||pdat[j].seq[uniq[k]]=='C'||pdat[j].seq[uniq[k]]=='G'||pdat[j].seq[uniq[k]]=='T') {
                        if (pdat[i].seq[uniq[k]]!=pdat[j].seq[uniq[k]]) {
                            dist++;
                        }
                    }
                }
            }
            seqdists[i][j]=dist;
            seqdists[j][i]=dist;
        }
    }
    for (int i=0;i<pdat.size();i++) {
        if (pdat[i].seq.size()==0) {
            for (int j=0;j<pdat.size();j++){
                seqdists[i][j]=-1;
                seqdists[j][i]=-1;
            }
            seqdists[i][i]=0;
        }
    }
    if (p.diagnostic==1) {
        PrintSeqDistances(seqdists);
    }
}

void FromConsensusDistances (const vector<sparseseq>& variants, vector< vector<tpair> >& seqdists_c) {
    //Construct a pairwise consensus between i and j as having a list of variants present in both sequences.
    //Find the number of variants each sequence would need to gain in order to be constructed from that consensus.
    //Used in the calculation of the amount of evolution occurring in each sequence via noise or error after the transmission event
    for (int i=0;i<variants.size();i++) {
        vector<tpair> vc;
        for (int j=0;j<variants.size();j++) {
            tpair p;
            p.from=0;
            p.to=0;
            //We want the number of variants in i not in j and vice versa.
            for (int k=0;k<variants[i].locus.size();k++) {
                int mycount = count(variants[j].locus.begin(), variants[j].locus.end(), variants[i].locus[k]);
                if (mycount==0) {
                    p.from++;
                }
            }
            for (int k=0;k<variants[j].locus.size();k++) {
                int mycount = count(variants[i].locus.begin(), variants[i].locus.end(), variants[j].locus[k]);
                if (mycount==0) {
                    p.to++;
                }
            }
            vc.push_back(p);
        }
        seqdists_c.push_back(vc);
    }
}













void ModifyVariantsByBin(int b, const vector< vector<int> >& bin, double& lLvar, const vector<pat>& pdat, vector<sparseseq>& variants) {
	int index=0;
	for (int i=0;i<pdat.size();i++) {
		for (int j=0;j<pdat[i].seq_uncertainty.size();j++) {
			if (bin[b][index]==1) {
				//Fix variantby adding to variants[i]
				variants[i].locus.push_back(pdat[i].seq_uncertainty[j].loc);
				variants[i].allele.push_back(pdat[i].seq_uncertainty[j].allele);
				lLvar=lLvar+log(pdat[i].seq_uncertainty[j].q);
			} else {
				//Don't fix - just change probability
				lLvar=lLvar+log(1-pdat[i].seq_uncertainty[j].q);
			}
			index++;
		}
	}
}

void FindSharedVariants(int set, vector<varcount>& vc_multi, const vector< vector<sparseseq> >& variants_sets) {
    //cout << "FindSharedVariants\n";
	vector<varcount> vc;
	for (int i=0;i<variants_sets[set].size();i++) {
		for (int j=0;j<variants_sets[set][i].locus.size();j++) {
			int found=0;
			for (int k=0;k<vc.size();k++) {
				if (variants_sets[set][i].locus[j]==vc[k].locus&&variants_sets[set][i].allele[j]==vc[k].allele) {
					vc[k].individual.push_back(i);
					found=1;
				}
			}
			if (found==0) {
				//cout << "Add " << variants_sets[set][i].locus[j] << " " << variants_sets[set][i].allele[j] << "\n";
				varcount v;
				v.locus=variants_sets[set][i].locus[j];
				v.allele=variants_sets[set][i].allele[j];
				v.individual.push_back(i);
				vc.push_back(v);
			}
		}
	}
    //cout << "Now here\n";
	for (int i=0;i<vc.size();i++) {
		if (vc[i].individual.size()>1) {
			vc_multi.push_back(vc[i]);
			/*cout << vc[i].locus << vc[i].allele <<" ";
			for (int j=0;j<vc[i].individual.size();j++) {
				cout << vc[i].individual[j] << " ";
			}
			cout << "\n";*/
		}
	}
	//Delete duplicate sets
	vector<int> to_rem;
	for (int i=0;i<vc_multi.size()-1;i++) {
		for (int j=i+1;j<vc_multi.size();j++) {
			if (vc_multi[i].individual==vc_multi[j].individual) {
				to_rem.push_back(j);
			}
		}
	}
	sort(to_rem.begin(),to_rem.end());
	to_rem.erase(unique(to_rem.begin(),to_rem.end()),to_rem.end());
	reverse(to_rem.begin(),to_rem.end());
	for (int i=0;i<to_rem.size();i++){
		vc_multi.erase(vc_multi.begin()+to_rem[i]);
	}
}




/*void GetSharedInOut (int c, const vector<tpairs>& trans_sets, const vector<varcount>& vc_multi, vector< vector<int> >& inout) {
    for (int v=0;v<vc_multi.size();v++) {
        vector<int> ino;
        for (int j=0;j<trans_sets[c].pair.size();j++) {
            int to_in=0;
            int from_in=0;
            for (int k=0;k<vc_multi[v].individual.size();k++) {
                if (trans_sets[c].pair[j].to==vc_multi[v].individual[k]) {
                    to_in=1;
                }
                if (trans_sets[c].pair[j].from==vc_multi[v].individual[k]) {
                    from_in=1;
                }
            }
            if (to_in==1&&from_in==0) { //This goes into the set.  Can only be one...
                ino.push_back(trans_sets[c].pair[j].from);
                ino.push_back(trans_sets[c].pair[j].to);

            }
        }
        inout.push_back(ino);
    }
}*/

/*void GetNonSharedInOut (int c, const vector<tpairs>& trans_sets, const vector<varcount>& vc_single, vector< vector<int> >& inout_s) {
    for (int v=0;v<vc_single.size();v++) {
        vector<int> ino;
        int found=0;
        for (int j=0;j<trans_sets[c].pair.size();j++) {
            if (trans_sets[c].pair[j].to==vc_single[v].individual[0]) {
                found=1;
                ino.push_back(trans_sets[c].pair[j].from);
                ino.push_back(trans_sets[c].pair[j].to);
            }
        }
        if (found==1) {
            inout_s.push_back(ino);
        } else {
            //This variant is in the root individual
            ino.push_back(-1);
            ino.push_back(-1);
            inout_s.push_back(ino);
        }
    }
}*/


