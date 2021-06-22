#include "basicmodel.h"
#include "find_adjacent.h"
#include "likelihoods.h"

void ListAdjacentK (run_params p, int k) {
    //Outputs list of adjacent networks
    
    vector<double> likelihoods;
    vector< vector<int> > trans_in;
    ReadBestLikelihood (p,likelihoods,trans_in);

    //Find downstream nodes
    vector< vector<int> > downstream;
    FindDownstreamNodes (trans_in,downstream);

    //Find all adjacent graphs - within 1.  Adds to trans_in.
    FindInitialAdjacent (p,0,downstream,trans_in);
    if (p.diagnostic==1) {
        PrintAdjacentNetworks(trans_in);
    }

    //Find all adjacent graphs to these - within k.  Add to trans_in.
    for (int add=2;add<=k;add++) {
        int end=trans_in.size();
        for (int i=1;i<end;i++) {
            FindInitialAdjacent (p,i,downstream,trans_in);
        }
        cout << "Size is now " << trans_in.size() << "\n";
    }
    
    //Remove duplicates
    cout << "Remove duplicates\n";
    RemoveDuplicatesTransIn(trans_in);
    cout << "Size is now " << trans_in.size() << "\n";

    cout << "Remove if known\n";
    vector< vector<int> > c_before;
    vector<double> like_before;
    RemoveIfLikelihoodKnown(p,trans_in,like_before,c_before);
    cout << "Size is now " << trans_in.size() << "\n";

    //Check that there is at least one plausible order
    vector< vector<int> > trans_group;
    FindTransGroup (trans_in,trans_group);

    cout << "Read orders\n";
    vector<tpairs> trans_sets;
    vector<int> new_subset;
    ReadOrders(trans_in,trans_group,trans_sets,new_subset);
    
    cout << "Size is now " << trans_sets.size() << "\n";

    PrintTransSetsFile(trans_sets);

}

void RemoveDuplicatesTransIn(vector< vector<int> >& trans_in) {
    vector<int> trans_in_keep;
    for (int i=0;i<trans_in.size();i++) {
        trans_in_keep.push_back(1);
    }
    for (int i=0;i<trans_in.size()-1;i++) {
        if (i%1000==0) {
            cout << i << " " << trans_in.size() << "\n";
        }
        if (trans_in_keep[i]==1) {
            for (int j=i+1;j<trans_in.size();j++) {
                if (trans_in[i]==trans_in[j]) {
                    trans_in_keep[j]=0;
                }
            }
        }
    }
    vector< vector<int> > trans_in_new;
    for (int i=0;i<trans_in.size();i++) {
        if (trans_in_keep[i]==1) {
            trans_in_new.push_back(trans_in[i]);
        }
    }
    trans_in=trans_in_new;
/*
    sort(to_rem.begin(),to_rem.end());
    to_rem.erase(unique(to_rem.begin(),to_rem.end()),to_rem.end());
    reverse(to_rem.begin(),to_rem.end());
    cout << "Size " << to_rem.size() << "\n";
    for (int i=0;i<to_rem.size();i++) {
        trans_in.erase(trans_in.begin()+to_rem[i]);
    }*/
}


void AnalyseAdjacent (run_params p, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, const vector< vector<sparseseq> >& variants_sets) {
    //Calculates likelihoods of adjacent networks
    vector<double> likelihoods;
    vector< vector<int> > trans_in;
    ReadBestLikelihood (p,likelihoods,trans_in);

    //Find downstream nodes
    vector< vector<int> > downstream;
    FindDownstreamNodes (trans_in,downstream);
    if (p.diagnostic==1) {
        //PrintDownstream(downstream);
    }

    //Find all adjacent graphs
    FindInitialAdjacent (p,0,downstream,trans_in);
    if (p.diagnostic==1) {
        //PrintAdjacentNetworks(trans_in);
    }

      //Check we don't already have a likelihood
      vector< vector<int> > c_before;
      vector<double> like_before;
      RemoveIfLikelihoodKnown(p,trans_in,like_before,c_before);
    
      //Check that there is at least one plausible order
      vector< vector<int> > trans_group;
      FindTransGroup (trans_in,trans_group);

      //Print Trans_group
      if (p.diagnostic==1) {
          //PrintTransGroup(trans_group);
      }

    cout << "Size " << trans_in.size() << "\n";
    vector<tpairs> trans_sets;
    vector<int> new_subset;
    ReadOrders(trans_in,trans_group,trans_sets,new_subset);
    if (p.diagnostic==1) {
        PrintTransSets(trans_sets);
    }
    cout << "Size now " << trans_sets.size() << "\n";

    //Other prep for running calculation
    vector< vector<int> > orders;
    vector<pat> new_pdat;
    vector< vector<ijlike> > new_like_trans;
    vector<sparseseq> new_variants;
    vector<varcount> vc_multi;
    vector<varcount> vc_single;
    SetupTreeLikelihoodCalculation(p,p.specify_set,new_subset,trans_sets,pdat_sets,like_trans_sets,variants_sets,orders,new_pdat,new_like_trans,new_variants,vc_multi,vc_single);

    ofstream cp_file;
    string name="Likelihoods_search.out";
    cp_file.open(name.c_str());
    //Print to file likelihoods which have previously been calculated
    //PrintPrevious (cp_file,like_before,c_before);

    //Calculate Likelihoods
    CalculateLikelihoods (p,cp_file,trans_sets,vc_multi,vc_single,new_pdat,new_like_trans,orders);
    
    cp_file.close();

}
    


void ReadBestLikelihood (run_params p, vector<double>& likelihoods, vector< vector<int> >& trans_in) {
    ifstream b_file;
    b_file.open("Best_likelihood.out");
    int n;
    double x;
    vector<int> temp;
    for (int i=0;i<p.n_edges+1;i++) {
        if (!(b_file >> n)) break;
        temp.push_back(n);
        //cout << n << " ";
    }
    b_file >> x;
    trans_in.push_back(temp);
    likelihoods.push_back(x);
    cout << x << "\n";
}

void FindDownstreamNodes (const vector< vector<int> >& trans_in, vector< vector<int> >& downstream) {
    //Represent impossible network changes
    for (int i=0;i<trans_in[0].size();i++) {
        vector<int> down;
        for (int j=0;j<trans_in[0].size();j++) {
            if (trans_in[0][j]==i) {
                down.push_back(j);
            }
        }
        int check=0;
        do {
            check=0;
            int count=down.size();
            for (int j=0;j<trans_in[0].size();j++) {
                for (int k=0;k<down.size();k++) {
                    if (trans_in[0][j]==down[k]) {
                        down.push_back(j);
                    }
                }
            }
            sort(down.begin(),down.end());
            down.erase(unique(down.begin(),down.end()),down.end());
            if (count==down.size()) {
                check=1;
            }
        }
        while (check==0);
        downstream.push_back(down);
    }
}

void PrintDownstream (const vector< vector<int> >& downstream) {
    for (int i=0;i<downstream.size();i++) {
        cout << i << " ";
        for (int j=0;j<downstream[i].size();j++) {
            cout << downstream[i][j] << " ";
        }
        cout << "\n";
    }
}

void FindInitialAdjacent (run_params p, int start, const vector< vector<int> >& downstream, vector< vector<int> >& trans_in) {
    for (int i=0;i<trans_in[start].size();i++) {
        if (trans_in[start][i]!=-1) {
            for (int j=0;j<p.n_edges+1;j++) {
                if (j!=trans_in[start][i]) {
                    int found=0;
                    for (int k=0;k<downstream[j].size();k++) {
                        if (downstream[j][k]==j) {
                            found=1;
                            break;
                        }
                    }
                    if (found==0) {
                        //Make new
                        vector<int> temp=trans_in[start];
                        temp[i]=j;
                        trans_in.push_back(temp);
                    }
                }
            }
        }
    }
}

void PrintAdjacentNetworks (const vector< vector<int> >& trans_in) {
    cout << "Adjacent networks\n";
    for (int i=1;i<trans_in.size();i++) {
        for (int j=0;j<trans_in[i].size();j++) {
            cout << trans_in[i][j] << " ";
        }
        cout << "\n";
    }
}

void RemoveIfLikelihoodKnown (run_params& p, vector< vector<int> >& trans_in, vector<double>& like_before, vector< vector<int> >& c_before) {
    vector<int> to_rem;
    for (int i=1;i<trans_in.size();i++) {
        int found=0;
        ifstream like_file;
        ostringstream convert;
        ostringstream convert2;
        convert << trans_in[i][0];
        convert2 << trans_in[i][2];
        string temp=convert.str();
        string temp2=convert2.str();
        string name="../Likelihood_data/Likelihoods"+temp+"_"+temp2+"_alt.out";
        like_file.open(name.c_str());
        int n;
        double x;
        for (int j=0;j<1000000;j++) {
            vector<int> origins;
            for (int i=0;i<p.n_edges+1;i++) {
                if (!(like_file >> n)) break;
                origins.push_back(n);
            }
            if (!(like_file >> x)) break;
            if (origins==trans_in[i]) {
                found=1;
                c_before.push_back(origins);
                like_before.push_back(x);
                break;
            }
        }
        if (found==1) {
            to_rem.push_back(i);
        }
        like_file.close();
    }
    //cout << to_rem.size() << "\n";
    if (to_rem.size()>1) {
        sort(to_rem.begin(),to_rem.end());
        reverse(to_rem.begin(),to_rem.end());
    }
    for (int i=0;i<to_rem.size();i++) {
        trans_in.erase(trans_in.begin()+to_rem[i]);
    }
}

void FindTransGroup (const vector< vector<int> >& trans_in, vector< vector<int> >& trans_group) {
    //Groups which come from the same ordering file.  Efficient for file read
    for (int i=0;i<trans_in.size();i++) {
        vector<int> tg;
        int found=0;
        for (int j=0;j<trans_group.size();j++) {
            if (trans_group[j][0]==trans_in[i][0]&&trans_group[j][1]==trans_in[i][2]) {
                trans_group[j].push_back(i);
                found=1;
            }
        }
        if (found==0) {
            int count=0;
            int j=0;
            while (count<2) {
                if (trans_in[i][j]!=-1) {
                    tg.push_back(trans_in[i][j]);
                    count++;
                }
                j++;
            }
          //  tg.push_back(trans_in[i][0]);
          //  tg.push_back(trans_in[i][2]);
            tg.push_back(i);
            trans_group.push_back(tg);
        }
    }
}

void PrintTransGroup (const vector< vector<int> >& trans_group) {
    for (int i=0;i<trans_group.size();i++) {
        for (int j=0;j<trans_group[i].size();j++) {
            cout << trans_group[i][j] << " ";
        }
        cout << "\n";
    }
}

void ReadOrders (const vector< vector<int> >& trans_in, const vector< vector<int> >& trans_group, vector<tpairs>& trans_sets, vector<int>& new_subset) {
    int read_new_subset=0;
    for (int i=0;i<trans_group.size();i++) {
        int found=0;
        ifstream order_file;
        ostringstream convert;
        ostringstream convert2;
        convert << trans_group[i][0];
        convert2 << trans_group[i][1];
        string temp=convert.str();
        string temp2=convert2.str();
        string name="../Order_data/Check2_"+temp+"_"+temp2+"_alt.out";
        cout << "Reading " << name << "\n";
        order_file.open(name.c_str());
        int n;
        int n1;
        int n2;
        int n3;
        int n4;
        for (int l=0;l<10000000;l++) {
            tpairs t;
            if (l==0) {
                //Read in indiviudals in the subset
                if (!(order_file >> n)) break;
                for (int j=0;j<n;j++) {
                    if (!(order_file >> n1)) break;
                    if (read_new_subset==0) {
                        new_subset.push_back(n1);
                    }
                }
                if (read_new_subset==0) {
                    read_new_subset=1;
                }
            }
            if (!(order_file >> n)) break;
            for (int j=0;j<n;j++) {
                if (!(order_file >> n1)) break;
                t.trans_in.push_back(n1);
            }
            if (!(order_file >> n)) break;
            t.root=n;
            if (!(order_file >> n)) break;
            for (int j=0;j<n;j++) {
                if (!(order_file >> n1)) break;
                if (!(order_file >> n2)) break;
                if (!(order_file >> n3)) break;
                if (!(order_file >> n4)) break;
                tpair tp;
                tp.from=n1;
                tp.to=n2;
                t.pair.push_back(tp);
                t.min_time.push_back(n3);
                t.max_time.push_back(n4);
            }
            if (!(order_file >> n)) break;
            t.order_index.clear();
            for (int j=0;j<n;j++) {
                if (!(order_file >> n1)) break;
                t.order_index.push_back(n1);
            }
            //Match to trans_group data
            for (int j=2;j<trans_group[i].size();j++) {
                if (t.trans_in==trans_in[trans_group[i][j]]) {
                    trans_sets.push_back(t);
                    break;
                }
            }
        }
    }
}

void PrintTransSets (vector<tpairs>& trans_sets) {
    for (int i=0;i<trans_sets.size();i++) {
        for (int j=0;j<trans_sets[i].trans_in.size();j++) {
            cout << trans_sets[i].trans_in[j] << " ";
        }
        cout << "\n";
    }
}

void PrintTransSetsFile (vector<tpairs>& trans_sets) {
    ofstream t_file;
    t_file.open("TransmissionSets.in");
    for (int i=0;i<trans_sets.size();i++) {
        for (int j=0;j<trans_sets[i].trans_in.size();j++) {
            t_file << trans_sets[i].trans_in[j] << " ";
        }
        t_file << "\n";
    }
}


void PrintPrevious (ofstream& cp_file, const vector<double>& like_before, const vector< vector<int> >& c_before) {
    for (int i=0;i<like_before.size();i++) {
        cp_file << "Chain B ";
        for (int j=0;j<c_before[i].size();j++) {
            cp_file << c_before[i][j] << "-" << j << " ";
        }
        cout << like_before[i] << "\n";
        
    }

}



void CalculateLikelihoods (run_params& p, ofstream& cp_file, vector<tpairs>& trans_sets, const vector<varcount>& vc_multi, const vector<varcount>& vc_single, const vector<pat>& new_pdat, vector< vector<ijlike> >& new_like_trans, const vector< vector<int> >& orders) {
    for (int c=0;c<trans_sets.size();c++) {
        //Make order c from order index
        MakeOrderFromIndex (c,orders,trans_sets);
        cout << "Set " << c << "\n";
        double chain_L=0;

        //Loop over orderings
        for (int o=0;o<trans_sets[c].orders.size();o++) {
            CalculateOrderLike (p,c,o,chain_L,trans_sets,vc_multi,vc_single,new_pdat,new_like_trans);
        }

        cp_file << "Chain " << c << " ";
        for (int k=0;k<trans_sets[c].ordered_pair.size();k++) {//Action to be done here is outside of next loop
            cp_file << trans_sets[c].ordered_pair[k].from << "-" << trans_sets[c].ordered_pair[k].to << " ";
        }
        cp_file << "Likelihood " << c << " " << log(chain_L) << "\n";
        //End of looking at chain c
        //Clear orderings
        trans_sets[c].orders.clear();
        
    }
}


void CalculateLikelihoodsFromList (run_params& p, const vector< vector<pat> >& pdat_sets, const vector< vector< vector<ijlike> > >& like_trans_sets, const vector< vector<sparseseq> >& variants_sets) {
    //Calculate likelihood of a list of imported networks.
    vector< vector<int> > trans_in;
    ifstream sets_file;
    sets_file.open("TransmissionSets.in");
    int n;
    double x;
    for (int i=0;i<10000000;i++) {
        vector<int> temp;
        for (int j=0;j<p.n_edges+1;j++) {
            if (!(sets_file >> n)) break;
            temp.push_back(n);
        }
        if (temp.size()==0) {
            break;
        }
        trans_in.push_back(temp);
    }
    cout << "Transmission events\n";
    cout << trans_in.size() << "\n";
    
    for (int i=0;i<trans_in.size();i++) {
        for (int j=0;j<trans_in[i].size();j++) {
            cout << trans_in[i][j] << " ";
        }
        cout << "\n";
    }
    vector< vector<int> > trans_group;
    FindTransGroup (trans_in,trans_group);
    vector<tpairs> trans_sets;
    vector<int> new_subset;
    ReadOrders(trans_in,trans_group,trans_sets,new_subset);
    
    //Other prep for running calculation
    vector< vector<int> > orders;
    vector<pat> new_pdat;
    vector< vector<ijlike> > new_like_trans;
    vector<sparseseq> new_variants;
    vector<varcount> vc_multi;
    vector<varcount> vc_single;
    SetupTreeLikelihoodCalculation(p,p.specify_set,new_subset,trans_sets,pdat_sets,like_trans_sets,variants_sets,orders,new_pdat,new_like_trans,new_variants,vc_multi,vc_single);
    
    ofstream cp_file;
    string name="Likelihoods_search.out";
    cp_file.open(name.c_str());
    //Print to file likelihoods which have previously been calculated
    //PrintPrevious (cp_file,like_before,c_before);

    //Calculate Likelihoods
    CalculateLikelihoods (p,cp_file,trans_sets,vc_multi,vc_single,new_pdat,new_like_trans,orders);

}
