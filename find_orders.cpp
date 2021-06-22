#include "basicmodel.h"
#include "checkpointing.h"
#include "find_orders.h"

void FindPlausibleOrders () {
    vector< vector<int> > check_info;
    ReadCheckpointInfo(check_info);
    for (int i=0;i<check_info.size();i++) {
        int set=check_info[i][0];
        vector<int> new_subset;
        vector<tpairs> trans_sets;
        cout << "Read data\n";
        ReadCheckpoint1 (set,check_info[i][1],new_subset,trans_sets);
        
        cout << "Set " << i << "\n";
        cout << "Size " << trans_sets.size() << "\n";
        
        vector<int> factorial;
        MakeFactorial(trans_sets[0].pair.size(),factorial);
        for (int i=0;i<factorial.size();i++) {
            cout << i << " " << factorial[i] << "\n";
        }
        
        //Make orderings list.  All possible orderings
        vector< vector<int> > orders;
        cout << "Make orderings\n";
        ConstructOrderingsNew(trans_sets,orders);
        cout << "Number of orderings " << orders.size() << "\n";
        
        //For each possible constraint find which orders are compatible.  Store all of these
        cout << "Find order filter " << trans_sets[0].pair.size() << "\n";
        vector< vector< vector<int> > > order_filter;
        FindOrderFilter (trans_sets[0].pair.size(),factorial,orders,order_filter);
        
        //Find orderings which match all of the relevant constraint filters
        FilterPossibleOrderings(orders.size(),order_filter,trans_sets);
        
      //  CleanTransSets(trans_sets);
      //  cout << "Number of sets with orders " << trans_sets.size() << "\n";
        
        //Make new checkpoint
        WriteCheckpoint2(check_info[i][0],check_info[i][1],new_subset,trans_sets);
    }
}

void MakeFactorial(int top, vector<int>& factorial) {
    int n=1;
    factorial.push_back(n);
    for (int i=1;i<=top;i++) {
        n=n*i;
        factorial.push_back(n);
    }
}

void ConstructOrderingsNew(vector<tpairs>& trans_sets, vector< vector<int> >& orders) {
    int s=trans_sets[0].pair.size();
    int order [trans_sets[0].pair.size()];
    vector<int> temp;
    for (int k=0;k<trans_sets[0].pair.size();k++) {
        temp.push_back(k);
        order[k]=k;
    }
    do {
        for (int k=0;k<trans_sets[0].pair.size();k++) {
            temp[k]=order[k];
        }
        orders.push_back(temp);
    } while (next_permutation(order,order+s));
    cout << "Done\n";
}

void FindOrderFilter (int size, const vector<int>& factorial, const vector< vector<int> >& orders, vector< vector< vector<int> > >& order_filter) {
    //For each ordered pair, does the ordering fit.  i.e. order_filter[i][j][k]=1 if i is before j in order k
    vector< vector<int> > orf;
    for (int j=0;j<size;j++) {
        vector<int> of;
        for (int k=0;k<size;k++) {
            orf.push_back(of);
        }
        order_filter.push_back(orf);
    }
    for (int before=0;before<size-1;before++) {
        for (int after=before+1;after<size;after++) {
            for (int j=0;j<orders.size();j++) {
                int c=CheckOrder(before,after,j,size,factorial);//1 if before and after in order otherwise 0
                int d=1-c;
                order_filter[before][after].push_back(c);
                order_filter[after][before].push_back(d);
            }
        }
    }
}

int CheckOrder(int i, int j, int k, int n, vector<int> factorial) {
    //Qn: Is j>i in the ordering k of n-digit orders?
    int t=n;
    while (1==1) {
        t--;
        int r=k/factorial[t];
        if (r==i) {
            return 1;
        }
        if (r==j) {
            return 0;
        }
        if (r<i) {
            i--;
        }
        if (r<j) {
            j--;
        }
        k=k%factorial[t];
    }
}

void FilterPossibleOrderings (int n_orders, const vector< vector< vector<int> > >& order_filter, vector<tpairs>& trans_sets) {
    for (int i=0;i<trans_sets.size();i++) {
        cout << "Set " << i << "\n";
        for (int j=0;j<n_orders;j++) {
            int include=1;
            for (int k=0;k<trans_sets[i].constraints.size();k++) {
                if (order_filter[trans_sets[i].constraints[k][0]][trans_sets[i].constraints[k][1]][j]==0) {
                    include=0;
                    break;
                }
            }
            if (include==1) {
                trans_sets[i].order_index.push_back(j);
            }
        }
    }
}

void CleanTransSets(vector<tpairs>& trans_sets) {
    vector<int> to_rem;
    for (int i=0;i<trans_sets.size();i++) {
        if (trans_sets[i].order_index.size()==0) {
            cout << "Nothing in " << i << "\n";
            to_rem.push_back(i);
        }
    }
    reverse(to_rem.begin(),to_rem.end());
    for (int i=0;i<to_rem.size();i++) {
        trans_sets.erase(trans_sets.begin()+to_rem[i]);
    }
}
