using namespace std;
#include "basicmodel.h"
#include "thresholds.h"
#include "distributions.h"

void CalculateThresholdsNoSeq (run_params p, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc) {
	cout << "Calculating NoSeq thresholds\n";
	ofstream f95;
	ofstream f99;
	f95.open("../Data/Threshold95NS.out");
	f99.open("../Data/Threshold99NS.out");
	vector<double> y=CalculateThresholdsNSSpecific(p,OGPreCalcP,LNPreCalc);
	f95 << y[0] << "\n";
	f99 << y[1] << "\n";
}

vector<double> CalculateThresholdsNSSpecific (run_params p, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc) {
	vector<double> allstats;
	//Assume that S1 is at time zero
	for (int s2=-11;s2<=87;s2++) {
		double L=0;
		for (int t=-11;t<=16;t++) {
			if (t<=s2) {
				double Ln=OffsetGammaPreCalcP(t,OGPreCalcP);
				//double L=OffsetGammaCDFFlex(t,p.pa,p.pb,p.po);
				Ln=Ln+LogNormalP(s2-t,LNPreCalc);
				//L=L+LogNormal(s2-t,p.smu,p.ssigma);
				Ln=exp(Ln);
				L=L+Ln;
			}
		}
		L=log(L);
		L=L+log(p.chat);
		allstats.push_back(L);
	}
	//Calculation
	double tot=0;
	for (int i=0;i<allstats.size();i++) {
		allstats[i]=exp(allstats[i]);
		tot=tot+allstats[i];
	}
	sort(allstats.begin(),allstats.end());
	reverse(allstats.begin(),allstats.end());
	vector<int> found;
	vector<double> prob;
	prob.push_back(0.95);
	prob.push_back(0.99);
	found.push_back(0);
	found.push_back(0);
	int min=0;
	vector<double> y;
	for (int i=1;i<allstats.size();i++) {
		double unlog=log(allstats[i]);
		allstats[i]=allstats[i]+allstats[i-1];
		for (int j=min;j<prob.size();j++) {
			if (allstats[i]>tot*prob[j]&&found[j]==0) {
				found[j]=1;
				min=j;
				//cout << "Threshold " << prob[j] << " " << unlog << " " << allstats[i]/tot << "\n";
				if (abs(prob[j]-0.95)<0.001||abs(prob[j]-0.99)<0.001) {
					y.push_back(unlog);
				}
			}
		}
	}
	return y;
}

void CalculateThresholdsFullExplicit (run_params p, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc, gsl_rng *rgen) {
	ofstream f95;
	ofstream f99;
	f95.open("../Data/Thresholds95.out");
	f99.open("../Data/Thresholds99.out");
	vector< vector<double> > thresholds95;
	vector< vector<double> > thresholds99;
	cout << "Calculating thresholds: (May be a few minutes)\n";
	cout << "---------------------------------------------------\n";
	for (int d1=-10;d1<=40;d1++) {
		vector<double> t95;
		vector<double> t99;
		for (int d2d=-10;d2d<=40;d2d++) {
			vector<double> y=CalculateThresholdsSpecificDExplicit(p,d1,d2d,OGPreCalcP,LNPreCalc,rgen);
			t95.push_back(y[0]);
			t99.push_back(y[1]);
		}
		thresholds95.push_back(t95);
		thresholds99.push_back(t99);
		cout << "|" << std::flush;
	}
	cout << "\n";
	for (int i=0;i<thresholds95.size();i++) {
		for (int j=0;j<thresholds95[i].size();j++) {
			f95 << i-10 << " " << j-10 << " " << thresholds95[i][j] << "\n";
		}
	}
	for (int i=0;i<thresholds95.size();i++) {
		for (int j=0;j<thresholds95[i].size();j++) {
			f99 << i-10 << " " << j-10 << " " << thresholds99[i][j] << "\n";
		}
	}
}

vector<double> CalculateThresholdsSpecificDExplicit (run_params p, int d1, int d2d, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc, gsl_rng *rgen) {
	vector<double> allstats;
	//Assume that S1 is at time zero
	double n=0.41369;
	double r=(0.0008*29782)/365.25;
	for (int s2=-11;s2<=87;s2++) {
		int d2=s2-d2d;
		for (int h=0;h<=10;h++) {
			for (int b=0;b<=h;b++) {
				int h1=b;
				int h2=h-b;
				vector<double> tlog;
				for (int t=-11;t<=16;t++) {
					double tl=0;
					if (t<=s2&&d2>=t) {
						if (d1<t) {
							int evotime2=d2-d1;
							tl=OffsetGammaPreCalcP(t,OGPreCalcP);
							tl=tl+LogNormalP(s2-t,LNPreCalc);
							double lambda=n/2;
							tl=tl+Poisson(h1,lambda);
							lambda=(n/2)+(evotime2*r);
							tl=tl+Poisson(h2,lambda);
							tl=exp(tl);
						} else {
							int evotime1=d1-t;
							int evotime2=d2-t;
							tl=OffsetGammaPreCalcP(t,OGPreCalcP);
							tl=tl+LogNormalP(s2-t,LNPreCalc);
							double lambda=(n/2)+(evotime1*r);
							tl=tl+Poisson(h1,lambda);
							lambda=(n/2)+(evotime2*r);
							tl=tl+Poisson(h2,lambda);
							tl=exp(tl);
						}
					}
					tlog.push_back(tl);
				}
				for (int rep=0;rep<100;rep++) {
					double logL=0;
					for (int r=0;r<tlog.size();r++) {
						if (gsl_ran_bernoulli(rgen,p.chat)==1) {
							logL=logL+tlog[r];
						}
					}
					logL=log(logL);
					allstats.push_back(logL);
				}
			}
		}
	}
	//Calculation
	double tot=0;
	for (int i=0;i<allstats.size();i++) {
		allstats[i]=exp(allstats[i]);
		tot=tot+allstats[i];
	}
	sort(allstats.begin(),allstats.end());
	reverse(allstats.begin(),allstats.end());
	vector<int> found;
	vector<double> prob;
	prob.push_back(0.95);
	prob.push_back(0.99);
	found.push_back(0);
	found.push_back(0);
	vector<double> y;
	int min=0;
	for (int i=1;i<allstats.size();i++) {
		double unlog=log(allstats[i]);
		allstats[i]=allstats[i]+allstats[i-1];
		for (int j=min;j<prob.size();j++) {
			if (allstats[i]>tot*prob[j]&&found[j]==0) {
				found[j]=1;
				min=j;
				if (abs(prob[j]-0.95)<0.001||abs(prob[j]-0.99)<0.001) {
					y.push_back(unlog);
				}
			}
		}
	}
	return y;
}
