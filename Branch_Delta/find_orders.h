#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void FindPlausibleOrders ();
void MakeFactorial(int top, vector<int>& factorial);
void ConstructOrderingsNew(vector<tpairs>& trans_sets, vector< vector<int> >& orders);
void FindOrderFilter (int size, const vector<int>& factorial, const vector< vector<int> >& orders, vector< vector< vector<int> > >& order_filter);
int CheckOrder(int i, int j, int k, int n, vector<int> factorial);
void FilterPossibleOrderings (int n_orders, const vector< vector< vector<int> > >& order_filter, vector<tpairs>& trans_sets);
void CleanTransSets(vector<tpairs>& trans_sets);
