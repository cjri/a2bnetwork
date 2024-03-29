using namespace std;
#include "basicmodel.h"
#include "diagnostics.h"
#include "likelihoods.h"
#include "threshold_data.h"
#include "io.h"

void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.dim=1;
	p.verb=0;
	//Interval between symptoms parameters
	//Gamma distribution from our generated compound distribution.
	p.ia=8.37419;
	p.ib=1.54413;
	p.io=7.19477;
	
	//Infectious potential relative to time of symptoms
	//Derived from Ashcroft with Bonhoeffer
	p.pa=97.18750;
	p.pb=0.268908;
	p.po=25.625;
	
	//Time of symptoms relative to date of infection
	//Lognormal distribution cited by He et al.
	p.smu=1.434065;
	p.ssigma=0.6612;
	
	//The following are learnt from CUH data.
	//They characterise the distribution of time from positive test to date of reporting symptoms
	p.ucta=2.5932152095707406;
	p.uctb=3.7760060663975437;
	p.ucto=3.112080041460921;
	p.uct_mean=6.67992;
	p.opt_uct=0;
	//Can create an option to optimise over the unknown times of becoming symptomatic
	p.rate=0.0008;
	p.turner=1;
	p.rem_no_seq=1;
	p.seq_noise=0.41369;
	p.chat=0.5;
    p.tight=0;
	p.threshold=-17.3;
    p.threshold_singleLL=-20;
	p.thresholdns=0;
	p.space=1; //By default read in spatial information
	p.read_subs=0;
	p.max_n=100;  //N.B. Have disabled max number of Ns at variant sites with this default
	p.seqvar=0;
    //Core files for analysis
    p.ali_file="NULL";
    p.pat_file="NULL";
    p.mov_file="NULL";
	p.ward_file="NULL";
    p.ward_bay_file="NULL";  //Patient movements at bay-level resolution i.e. ward data with more info
    p.extra_mov_file="NULL"; //Time non-specific HCW ward assignments
    p.ward_format_old=0;
    p.error=0;
    //May want to read in subset and likelihood data
	p.sub_file="Subsets1.in";
    p.likelihood_file="Likelihoods_0_1.out";
    p.time_file="Timings.out";
    //Delimiters for date information
	p.pat_delim='/';
	p.mov_delim='.';
	p.diagnostic=0;
    p.hcw_location_default=0.5714286;
	p.noseq=0;
	p.utopia=0;
	p.calc_thresholds=0;
	p.hcw_gap=0.5;
    p.specify_set=-1;
    p.specify_remove=-1;
    p.c_start=0;
    p.c_step=1;
    p.calculation=1;
    p.n_edges=5;
    p.consistency=1;
    p.distance=1;
    p.absolute=0;
    p.strict_locations=0;
    p.include_HCW_U=0;
    p.delta=0;
	int x=1;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--dim")==0) {
			x++;
			p.dim=atoi(argv[x]);
		} else if (p_switch.compare("--verb")==0) {
			x++;
			p.verb=atoi(argv[x]);
		} else if (p_switch.compare("--ia")==0) {
			x++;
			p.verb=atof(argv[x]);
		} else if (p_switch.compare("--ib")==0) {
			x++;
			p.verb=atof(argv[x]);
		} else if (p_switch.compare("--io")==0) {
			x++;
			p.verb=atof(argv[x]);
		} else if (p_switch.compare("--opt_uct")==0) {
			x++;
			p.opt_uct=atoi(argv[x]);
		} else if (p_switch.compare("--evo_rate")==0) {
			x++;
			p.rate=atof(argv[x]);
		} else if (p_switch.compare("--seq_noise")==0) {
			x++;
			p.seq_noise=atof(argv[x]);
		} else if (p_switch.compare("--calc_stats")==0) {
			x++;
			p.turner=1-atoi(argv[x]);
		} else if (p_switch.compare("--diag")==0) {
			x++;
			p.diagnostic=atoi(argv[x]);
        } else if (p_switch.compare("--delta")==0) {
            x++;
            p.delta=atoi(argv[x]);
		} else if (p_switch.compare("--noseq")==0) {
			x++;
			p.noseq=atoi(argv[x]);
		} else if (p_switch.compare("--noplace")==0) {
			x++;
			p.utopia=atoi(argv[x]);
		} else if (p_switch.compare("--ali_file")==0) {
			x++;
			p.ali_file=argv[x];
		} else if (p_switch.compare("--pat_file")==0) {
			x++;
			p.pat_file=argv[x];
		} else if (p_switch.compare("--mov_file")==0) {
			x++;
			p.mov_file=argv[x];
        } else if (p_switch.compare("--extra_mov_file")==0) {
            x++;
            p.extra_mov_file=argv[x];
        } else if (p_switch.compare("--ward_file")==0) {
            x++;
            p.ward_file=argv[x];
        } else if (p_switch.compare("--ward_bay_file")==0) {
            x++;
            p.ward_bay_file=argv[x];
        } else if (p_switch.compare("--ward_format_old")==0) {
            x++;
            p.ward_format_old=atoi(argv[x]);
        } else if (p_switch.compare("--strict_locations")==0) {
            x++;
            p.strict_locations=atoi(argv[x]);
        } else if (p_switch.compare("--include_U")==0) {
            x++;
            p.include_HCW_U=atoi(argv[x]);
		} else if (p_switch.compare("--sub_file")==0) {
			x++;
			p.sub_file=argv[x];
        } else if (p_switch.compare("--time_file")==0) {
            x++;
            p.time_file=argv[x];
		} else if (p_switch.compare("--spatial_data")==0) {
			x++;
			p.space=atoi(argv[x]);
		} else if (p_switch.compare("--rem_no_seq")==0) {
			x++;
			p.rem_no_seq=atoi(argv[x]);
		} else if (p_switch.compare("--threshold")==0) {
			x++;
			p.threshold=atoi(argv[x]);
        } else if (p_switch.compare("--tight")==0) {
            x++;
            p.tight=atoi(argv[x]);
		} else if (p_switch.compare("--maxn")==0) {
			x++;
			p.max_n=atoi(argv[x]);
		} else if (p_switch.compare("--read_subs")==0) {
			x++;
			p.read_subs=atoi(argv[x]);
		} else if (p_switch.compare("--seq_var")==0) {
			x++;
			p.seqvar=atoi(argv[x]);
		} else if (p_switch.compare("--calc_thresholds")==0) {
			x++;
			p.calc_thresholds=atoi(argv[x]);
        } else if (p_switch.compare("--calculation")==0) {
            x++;
            p.calculation=atoi(argv[x]);
		} else if (p_switch.compare("--hcw_gap")==0) {
			x++;
			p.hcw_gap=atof(argv[x]);
        } else if (p_switch.compare("--specify_set")==0) {
            x++;
            p.specify_set=atoi(argv[x]);
        } else if (p_switch.compare("--specify_remove")==0) {
            x++;
            p.specify_remove=atoi(argv[x]);
        } else if (p_switch.compare("--c_start")==0) {
            x++;
            p.c_start=atoi(argv[x]);
        } else if (p_switch.compare("--c_step")==0) {
            x++;
            p.c_step=atoi(argv[x]);
        } else if (p_switch.compare("--likelihood_file")==0) {
            x++;
            p.likelihood_file=argv[x];
        } else if (p_switch.compare("--n_edges")==0) {
            x++;
            p.n_edges=atoi(argv[x]);
        } else if (p_switch.compare("--distance")==0) {
            x++;
            p.distance=atoi(argv[x]);
        } else if (p_switch.compare("--absolute")==0) {
            x++;
            p.absolute=atoi(argv[x]);
       } else {
			cout << "Incorrect usage " << argv[x] << "\n";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
	p.rate=(p.rate*29900)/365.25;
	if (p.noseq==1) {
		p.seqvar=0;
	}
    /*DELTA PARAMETERS*/
    if (p.delta==1) {
        p.pa=38.4805;
        p.pb=0.468049;
        p.po=20;
        p.smu=1.39599;
        p.ssigma=0.41354;
    }
	SetThreshold(p);
}


/*
 The following routine is deprecated in favour of reading in data directly through code, but it allows threshold data to be read in from external files.  Example files are provided in the Data directory on the GitHub page
 
void GetThresholds (vector< vector<double> >& thresholds95, vector< vector<double> >& thresholds99, double& t95NS, double& t99NS, int& error) {
	t95NS=0;
	t99NS=0;
	ifstream t95;
	ifstream t99;
	ifstream t95n;
	ifstream t99n;
	t95.open("../Data/Thresholds95.dat");
	t99.open("../Data/Thresholds99.dat");
	t95n.open("../Data/Thresholds95NS.dat");
	t99n.open("../Data/Thresholds99NS.dat");
	int n;
	double x;
	t95n >> x;
	t95NS=x;
	t99n >> x;
	t99NS=x;
	if (t95NS>-0.01) {
		cout << "Error reading NS threshold: File not found ../Data/Thresholds95NS.dat\n";
		cout << "Note: The missing file can be generated using the --calc_thresholds 1 flag\n";
		error=1;
	}
	if (t99NS>-0.01) {
		cout << "Error reading NS threshold: File not found ../Data/Thresholds99NS.dat\n";
		cout << "Note: The missing file can be generated using the --calc_thresholds 1 flag\n";
		error=1;
	}
	int index1=-10;
	vector<double> t;
	for (int i=0;i<1000000;i++) {
		if (!(t95 >> n)) break;
		if (!(t95 >> n)) break;
		if (!(t95 >> x)) break;
		t.push_back(x);
		index1++;
		if (index1>40) {
			index1=-10;
			thresholds95.push_back(t);
			t.clear();
		}
	}
	index1=-10;
	for (int i=0;i<1000000;i++) {
		if (!(t99 >> n)) break;
		if (!(t99 >> n)) break;
		if (!(t99 >> x)) break;
		t.push_back(x);
		index1++;
		if (index1>40) {
			index1=-10;
			thresholds99.push_back(t);
			t.clear();
		}
	}
	if (thresholds95.size()<1) {
		cout << "Error reading thresholds.  Can't find file ../Data/Thresholds95.dat: Check this\n";
		cout << "Note: The missing file can be generated using the --calc_thresholds 1 flag\n";
		error=1;
	}
	if (thresholds99.size()<1) {
		cout << "Error reading thresholds.  Can't find file ../Data/Thresholds99.dat: Check this\n";
		cout << "Note: The missing file can be generated using the --calc_thresholds 1 flag\n";
		error=1;
	}
}
*/


void GetThresholdsInternal (run_params& p, vector< vector<double> >& thresholds95, vector< vector<double> >& thresholds99, double& t95NS, double& t99NS, int& error) {
    if (p.delta==1) {
        if (p.seq_noise==0) {
            ReadThresholdData_Delta_NoNoise (thresholds95,thresholds99,t95NS,t99NS);
        } else {
        ReadThresholdData_Delta (thresholds95,thresholds99,t95NS,t99NS);
        }
    } else {
        ReadThresholdData (thresholds95,thresholds99,t95NS,t99NS);
    }
    if (thresholds95.size()==0|| thresholds99.size()==0) {
        cout << "Error importing threshold data.  Check code...\n";
    }
}


void SetThreshold (run_params& p) {
    //Likelihood thresholds used in subsetting
    p.threshold=-10.1578;
    p.thresholdns=-6.17563;
}


void ReadPatFromCSV(run_params p, vector<pat>& pdat) {
    //Basic data on individuals
	ifstream csv_file;
	csv_file.open(p.pat_file.c_str());
	string str;
	for (int i=0;i<100000;i++) {
		if (!(csv_file >> str)) break;
		if (i>0) {
			pat pt;
			//Edit string to remove "
			RemovePunc(str);

			//Split by commas
			vector<string> subs;
			SplitCommas(str,subs);

            pt.code=subs[0];
			pt.code_match=subs[4];
			pt.type=atoi(subs[3].c_str());
			if (pt.type==3) {
				pt.hcw=1;
			} else {
				pt.hcw=0;
			}

			//Sort out the dates
			int temp[] = {1,5};
			vector<int> d(temp,temp+sizeof(temp) / sizeof(int));
			for (int j=0;j<d.size();j++) {
				vector<int> dmy;
				MakeDMY(d[j],subs,p.pat_delim,dmy);
				int day=DatetoDay(dmy);
				cout << dmy[0] << " " << dmy[1] << " " << dmy[2] << " " << day << "\n";
				if (d[j]==1) {
					pt.time_s=day;
				}
				if (d[j]==5) {
					pt.time_seq=day;
				}
			}
			if (subs[2]=="1") {
				pt.time_s_cert=1;
			} else {
                pt.time_s_cert=0;
                pt.time_s=pt.time_s-floor(p.uct_mean+0.5);

			}
			
			cout << pt.code << " " << pt.code_match << " " << pt.hcw << " " << pt.type << " " << pt.time_s << " " << pt.time_seq << " " << pt.time_s_cert << "\n";
			pdat.push_back(pt);
		}
	}
}

void RemovePunc(string& str) {
    //Edit string to remove "
    vector<int> rem;
    for (int j=0;j<str.size();j++) {
        if (str[j]=='"') {
            rem.push_back(j);
        }
    }
    reverse(rem.begin(),rem.end());
    for (int j=0;j<rem.size();j++) {
        str.erase(str.begin()+rem[j]);
    }
}

void SplitCommas(const string str, vector<string>& subs) {
    stringstream ss(str);
    while (ss.good() ) {
        string sr;
        getline(ss,sr,',');
        subs.push_back(sr);
    }
}

void MakeDMY (const int j, const vector<string>& subs, char delim, vector<int>& dmy) {
//cout << j << " " << subs.size() << " " << subs[j] << " ";
    stringstream sss(subs[j]);
    while (sss.good()) {
        string sr;
        getline(sss,sr,delim);
        dmy.push_back(atoi(sr.c_str()));
    }
}

void MakeDMY (run_params& p, string file, const int j, const vector<string>& subs, char delim, vector<int>& dmy) {
//cout << j << " " << subs.size() << " " << subs[j] << " ";
    stringstream sss(subs[j]);
    while (sss.good()&&p.error==0) {
        string sr;
        getline(sss,sr,delim);
        vector<string> months;
        months.push_back("Jan");
        months.push_back("Feb");
        months.push_back("Mar");
        months.push_back("Apr");
        months.push_back("May");
        months.push_back("Jun");
        months.push_back("Jul");
        months.push_back("Aug");
        months.push_back("Sep");
        months.push_back("Oct");
        months.push_back("Nov");
        months.push_back("Dec");
        int found=0;
        int add=-1;
        for (int i=0;i<months.size();i++) {
            if (sr.compare(months[i])==0) {
                add=i+1;
                found=1;
                break;
            }
        }
        if (found==1) {
            dmy.push_back(add);
        } else {
            dmy.push_back(atoi(sr.c_str()));
        }
    }
    if (p.error==0) {
        if (dmy.size()!=3) {
            DateError(file,subs[j],delim);
            p.error=1;
        }
        if (dmy[2]>50&&dmy[2]<100) {
            YearOOR(file,subs[j]);
            p.error=1;
        }
        if (dmy[2]>99&&dmy[2]<2000) {
            YearOOR(file,subs[j]);
            p.error=1;
        }
        if (dmy[2]>2050) {
            YearOOR(file,subs[j]);
            p.error=1;
        }
    }
}

void DateError (string file, string subsj, char delim) {
    cout << "Error reading in date\n" << subsj << " in " << file << "\n";
    cout << "Format must be DD" << delim << "MM" << delim << "YY or DD" << delim << "MM" << delim << "YYYY\n";
}

void YearOOR (string file, string subsj) {
    cout << "Error: Year out of range in entry\n" << subsj << " in " << file << "\n";
    cout << "Possible fix: check date is in DD/MM/YY or DD/MM/YYYY format \n";
}

void RemoveSpaces(string file, int i, int& pr, vector<string>& subs) {
    for (int j=0;j<subs.size();j++) {
        //Remove spaces
        vector<int> to_rem;
        for (int k=0;k<subs[j].size();k++){
            if (subs[j].compare(k,1," ")==0) {
                to_rem.push_back(k);
            }
        }
        sort(to_rem.begin(),to_rem.end());
        reverse(to_rem.begin(),to_rem.end());
        for (int k=0;k<to_rem.size();k++) {
            if (pr==0) {
                cout << "Warning: Removing spaces from line " << i << " of " << file << "\n";
                pr=1;
            }
            subs[j].erase(subs[j].begin()+to_rem[k]);
        }
    }
}


int DatetoDay (vector<int>& dmy) {
    //cout << dmy[0] << " " << dmy[1] << " " << dmy[2] << "\n";
    int day=dmy[0];
    int lp=1;
    //Convert to YY format from YYYY
    if (dmy[2]>1999) {
        dmy[2]=dmy[2]-2000;
    }
    //Note: Here 0 indicates the year 2020.
    
    //Year calculation
    int diff=dmy[2]-20;
    day=day+(diff*365);
    if (diff>0) {
        int leap=diff/4;
        day=day+leap+1;
    }
    if (diff%4!=3) {
        lp=0;
    }

    //Month calculation
    if (dmy[1]>1) {
        day=day+31;
    }
    if (dmy[1]>2) {
        day=day+28+lp;
    }
    if (dmy[1]>3) {
        day=day+31;
    }
    if (dmy[1]>4) {
        day=day+30;
    }
    if (dmy[1]>5) {
        day=day+31;
    }
    if (dmy[1]>6) {
        day=day+30;
    }
    if (dmy[1]>7) {
        day=day+31;
    }
    if (dmy[1]>8) {
        day=day+31;
    }
    if (dmy[1]>9) {
        day=day+30;
    }
    if (dmy[1]>10) {
        day=day+31;
    }
    if (dmy[1]>11) {
        day=day+30;
    }
    if (dmy[1]>12) {
        day=day+31;
    }
    //cout << "Day " << day << "\n";
    return day;
}


void ReadFastaAli (run_params p, vector<string>& names, vector<string>& seqs) {
    ifstream ali_file;
    ali_file.open(p.ali_file.c_str());
    //cout << p.ali_file.c_str() << "\n";
    int index=0;
    for (int i=0;i<1000000;i++) {
        if (index==0) {
            string name;
            if (!(ali_file >> name)) break;
            names.push_back(name);
            index=1-index;
        } else {
            string seq;
            if (!(ali_file >> seq)) break;
            seqs.push_back(seq);
            index=1-index;
        }
    }
}

void CheckBaseCase (vector<string>& seqs) {
    cout << "Check base case\n";
    for (int i=0;i<seqs.size();i++) {
        for (int j=0;j<seqs[i].size();j++) {
            if (seqs[i].compare(j,1,"a")==0) {
                seqs[i][j]='A';
            } else if (seqs[i].compare(j,1,"c")==0) {
                seqs[i][j]='C';
            } else if (seqs[i].compare(j,1,"g")==0) {
                seqs[i][j]='G';
            } else if (seqs[i].compare(j,1,"t")==0) {
                seqs[i][j]='T';
            } else if (seqs[i].compare(j,1,"n")==0) {
                seqs[i][j]='N';
            } else if (seqs[i].compare(j,1,"-")==0) {
                seqs[i][j]='N';
            }
        }
    }
}

void MatchSequencesToIndividuals (const run_params p, vector<pat>& pdat, vector<string>& names, vector<string>& seqs) {
    //Alter names of sequences to match input format>
    CorrectNames(names);

    //Add sequence data to records
    IncorporateSequenceData(pdat,names,seqs);
    RemoveIndividualsNoSequence(p,pdat);
    IfNoSeqUseDefaultSequence (p,pdat);
    
    //Calculate number of Ns in each sequence.
    CountNsInEachSequence(pdat);
    
    //Remove duplicate records for individuals with more than one sequence.  Take the earliest, best quality
    ResolveFilterRepeatPatients(pdat); //QC for sequence data
}

void CorrectNames (vector<string>& names) {
    //Process names in FASTA file to match input names
    for (int i=0;i<names.size();i++) {
        names[i].erase(0,1);
        stringstream a(names[i]);
        string segment;
        vector<string> seglist;
        while(getline(a, segment, '/')) {
           seglist.push_back(segment);
        }
        names[i]=seglist[0];
    }
}

void IncorporateSequenceData (vector<pat>& pdat, vector<string>& names, vector<string>& seqs) {
    //Add sequence data into the main data object
    for (int i=0;i<names.size();i++) {
        //cout << "Name " << i << " of " << names.size() << " " << names[i] << "\n";
        for (int j=0;j<pdat.size();j++) {
            if (pdat[j].code_match.compare(names[i])==0) {
                cout << "Match " << j << " " << pdat[j].code << " " << pdat[j].code_match << "\n";
                pdat[j].seq=seqs[i];
            }
        }
    }
}

void RemoveIndividualsNoSequence(const run_params p, vector<pat>& pdat) {
    //Given flags, remove individuals with no sequence data
    if (p.noseq==0||p.rem_no_seq==1) { //Can flag this to not remove data for these cases
        vector<int> noseq;
        for (int i=0;i<pdat.size();i++) {
            if (pdat[i].seq.size()==0) {
                cout << "No sequence data for individual " << pdat[i].code << " " << pdat[i].code_match << ": Excluding from input\n";
                noseq.push_back(i);
            }
        }
        sort(noseq.begin(),noseq.end());
        reverse(noseq.begin(),noseq.end());
        for (int i=0;i<noseq.size();i++) {
            pdat.erase(pdat.begin()+noseq[i]);
        }
    }
}

void IfNoSeqUseDefaultSequence (const run_params p, vector<pat>& pdat) {
    //Where genome sequence data is neglected, copy the first sequence to everything else
    if (p.noseq==1) {
        for (int i=1;i<pdat.size();i++) {
            pdat[i].seq=pdat[0].seq;
        }
    }
}

void CountNsInEachSequence (vector<pat>& pdat) {
    cout << "Count N \n";
    //Number of Ns in each sequence.
    //Filter for sequence quality
    //Tiebreaker for samples when there are more than one from the same day.
    for (int i=0;i<pdat.size();i++) {
        int nn=0;
        for (int j=0;j<pdat[i].seq.size();j++) {
            if (pdat[i].seq.compare(j,1,"N")==0) {
                nn++;
            }
        }
        cout << i << " " << pdat[i].seq.size() << " " << nn << "\n";
        pdat[i].seq_n=nn;
    }
}

void ResolveFilterRepeatPatients (vector<pat>& pdat) {
    cout << "Resolve Repeat Patients: QC for sequence data\n";
    //Removes individuals with more than one listing i.e. sequence.  Takes the first sequence ech time
    //If sequences have identical dates, takes the one with the fewest N nucleotides
    vector<int> rem;
    for (int i=0;i<pdat.size();i++) {
        if (pdat[i].seq_n<0.2*pdat[i].seq.size()) {
            //Remove anything later
            for (int j=0;j<pdat.size();j++) {
                if (pdat[i].code==pdat[j].code) {
                    if (pdat[i].time_seq<pdat[j].time_seq) {
                        rem.push_back(j);
                        //cout << "Remove " << pdat[j].code << " " << pdat[j].code_match << " " << j << " " << i << "\n";
                        RepairSequence(i,j,pdat);
                    }
                    //Remove anything with the same time of worse quality
                    if (pdat[i].time_seq==pdat[j].time_seq) {
                        if (pdat[i].seq_n<pdat[j].seq_n) {
                            rem.push_back(j);
                            //cout << "Remove " << pdat[j].code << " " << pdat[j].code_match << " " << j << " " << i << "\n";
                            RepairSequence(i,j,pdat);
                        }
                    }
                    //Remove anything of poor quality
                    if (pdat[j].seq_n>0.2*pdat[i].seq.size()) {
                        rem.push_back(j);
                        //cout << "Remove " << pdat[j].code << " " << pdat[j].code_match << " " << j << " " << i << "\n";
                        RepairSequence(i,j,pdat);
                    }
                }
            }
        } else {
            //Low quality.  Remove anything of worse quality
            for (int j=0;j<pdat.size();j++) {
                if (pdat[i].code==pdat[j].code) {
                    if (pdat[i].seq_n<pdat[j].seq_n) {
                        rem.push_back(j);
                        cout << "Remove " << pdat[j].code << " " << pdat[j].code_match << " " << j << " " << i << "\n";
                        RepairSequence(i,j,pdat);
                    }
                }
            }
        }
    }
    sort(rem.begin(),rem.end());
    rem.erase(unique(rem.begin(),rem.end()),rem.end());
    reverse(rem.begin(),rem.end());
    for (int i=0;i<rem.size();i++) {
        pdat.erase(pdat.begin()+rem[i]);
    }
    cout << "Data quality\n";
    for (int i=0;i<pdat.size();i++) {
        cout << pdat[i].code << " " << pdat[i].code_match << " " << pdat[i].seq_n << "\n";
    }
}

void RepairSequence(int i, int j, vector<pat>& pdat) {
    for (int k=0;k<pdat[i].seq.size();k++) {
        if (pdat[i].seq.compare(k,1,"N")==0) {
            if (pdat[j].seq.compare(k,1,"N")!=0) {
                pdat[i].seq[k]=pdat[j].seq[k];
            }
        }
    }
}

void ReadLocationData (run_params p, vector<pat>& pdat) {
    if (p.space==1) {
        if (p.utopia==2) {
            //Removed location data
            EnforceUtopia(pdat);
        } else if (p.utopia==1) {
            //People there between -15 and 25? Check a@B            EnforceUtopia(pdat);
            EnforceModerateUtopia(pdat);
        } else {
            int pd=0;
            int hd=0;
            if (p.mov_file.compare("NULL")!=0) {
                ReadHCWMovFromCSV(p,pdat);
                EditHCWMovData(p,pdat); //12 hour window of uncertainty - days with probability p.hcw_gap by default = 0.5
                hd=1;
            }
            //HCW data where location is known but no dates
            if (p.extra_mov_file.compare("NULL")!=0) {
                ReadExtraHCWMovFromCSV(p,pdat);
                hd=1;
            }
            if (hd==0) {
                cout << "Warning: No HCW location data\n";
            }
            if (p.ward_file.compare("NULL")!=0) {
                if (p.ward_format_old==1) {
                    ReadWardMovFromCSV(p,pdat);
                } else {
                    ReadWardMovFromCSVNew(p,pdat);
                }
                pd=1;
            }
            if (p.ward_bay_file.compare("NULL")!=0) {
                if (pd==1) {
                    cout << "Warning: Have specified two types of patient location data\n";
                }
                ReadWardMovBaysFromCSV(p,pdat);
                pd=1;
            }
            if (pd==0) {
                cout << "Warning: No patient location data\n";
            }
        }
        if (p.diagnostic==1) {
            PrintPdat(pdat);
        }
        p.opt_uct=0; //Use the mean time for cases of unknown infection date.
    }
    //Remove individuals with no sequence data or no location data
    if (p.utopia==0) {
        RemoveIndividualsNoLocation (p,pdat);
    }
}

void EnforceUtopia (vector<pat>& pdat) {
    //Everyone is everywhere; extreme solution
    cout << "No place information\n";
    for (int i=0;i<pdat.size();i++) {
        loc l;
        l.ward="WARD_01";
        if (pdat[i].hcw==1) {
            l.prob=4./7.;
        } else {
            l.prob=1;
        }
        for (int j=0;j<200;j++) {
            l.date=j;
            pdat[i].locat.push_back(l);
        }
    }
}

void EnforceModerateUtopia (vector<pat>& pdat) {
    //Use patient symptom onset time to decide who was where and when
    cout << "No place information: Use symptom onset times to construct this.\n";
    for (int i=0;i<pdat.size();i++) {
        loc l;
        l.ward="WARD_01";
        if (pdat[i].hcw==1) {
            l.prob=4./7.;
        } else {
            l.prob=1;
        }
        for (int j=pdat[i].time_s-15;j<=pdat[i].time_s+25;j++) {
            l.date=j;
            pdat[i].locat.push_back(l);
        }
    }
}


/*void ReadHCWMovFromCSV(run_params p, vector<pat>& pdat) {
	ifstream csv_file;
	csv_file.open(p.mov_file.c_str());
	string str;
	vector<int> dates;
	for (int i=0;i<100000;i++) {
		if (!(csv_file >> str)) break;
		//Edit string to remove "
		RemovePunc(str);
		//Split at commas
		vector<string> subs;
		SplitCommas(str,subs);
		if (i==0) {
			//Date information
			for (int j=2;j<subs.size();j++) {
				vector<int> dmy;
				MakeDMY(j,subs,p.mov_delim,dmy);
				int day=DatetoDay(dmy);
				//cout << subs[j] << " " << day << "\n";
				dates.push_back(day);
			}
		} else {
			//Put information into pdat location information
			string pat=subs[0];
			string w;
            //Convert letter code to anonymised ward numbers
			if (subs[1].compare("A")==0) {
				w="WARD_19";
			}
            if (subs[1].compare("B")==0) {
                w="WARD_47";
            }
            if (subs[1].compare("C")==0) {
                w="WARD_10";
            }
			if (subs[1].compare("D")==0) {
				w="WARD_36";
			}
			if (subs[1].compare("E")==0) {
				w="WARD_52";
			}

			int index=-1;
			for (int j=0;j<pdat.size();j++) {
				if (pdat[j].code==pat) {
					index=j;
					break;
				}
			}
			if (index!=-1) {
				for (int j=2;j<subs.size();j++) {
					if (subs[j].compare("Y")==0) {
						loc l;
						l.ward=w;
						l.date=dates[j-2];
						l.prob=1;
						pdat[index].locat.push_back(l);
					}
				}
			}
		}
	}
}*/

void ReadHCWMovFromCSV(run_params& p, vector<pat>& pdat) {
    ifstream csv_file;
    csv_file.open(p.mov_file.c_str());
    string str;
    vector<int> dates;
    int i=-1;
    int spwarn=0;
    while (getline(csv_file,str)) {
        i++;
        //Edit string to remove "
        RemovePunc(str);
        //Split at commas
        vector<string> subs;
        SplitCommas(str,subs);
        RemoveSpaces(p.mov_file.c_str(),i,spwarn,subs);
        if (i==0) {
            //Date information
            for (int j=2;j<subs.size();j++) {
                vector<int> dmy;
                MakeDMY(p,p.mov_file.c_str(),j,subs,p.mov_delim,dmy);
                int day=DatetoDay(dmy);
                //cout << subs[j] << " " << day << "\n";
                dates.push_back(day);
            }
            cout << "Number of dates " << dates.size() << "\n";
        } else {
            //Put information into pdat location information
            string pat=subs[0];
            string w=subs[1];
            int index=-1;
            for (int j=0;j<pdat.size();j++) {
                if (pdat[j].code==pat) {
                    index=j;
                    break;
                }
            }
            if (index!=-1) {
                for (int j=2;j<subs.size();j++) {
                    if (subs[j].compare(0,1,"Y")==0||subs[j].compare(0,1,"A")==0) { //A is ambiguous ward location
                        loc l;
                        l.ward=w;
                        l.date=dates[j-2];
                        l.prob=1;
                        pdat[index].locat.push_back(l);
                    }
                    if (p.include_HCW_U==1&&subs[j].compare(0,1,"U")==0) { //Include unknown locations under flag condition
                        loc l;
                        l.ward=w;
                        l.date=dates[j-2];
                        l.prob=1;
                        pdat[index].locat.push_back(l);
                    }
                }
            }
        }
    }
}


void ReadExtraHCWMovFromCSV(run_params& p, vector<pat>& pdat) {
    //Wards in which HCWs are always located
    //First, find minimum and maximum times
    int min=100000;
    int max=-100000;
    for (int i=0;i<pdat.size();i++) {
        for (int j=0;j<pdat[i].locat.size();j++) {
            if (pdat[i].locat[j].date>max) {
                max=pdat[i].locat[j].date;
            }
            if (pdat[i].locat[j].date<min) {
                min=pdat[i].locat[j].date;
            }
        }
        if (pdat[i].time_s+25>max) {
            max=pdat[i].time_s+25;
        }
        if (pdat[i].time_s-15<min) {
            min=pdat[i].time_s-15;
        }
    }
    ifstream csv_file;
    csv_file.open(p.extra_mov_file.c_str());
    string str;
    int i=-1;
    while (getline(csv_file,str)) {
        i++;
        if (i>0) {
            RemovePunc(str);
            //Split at commas
            vector<string> subs;
            SplitCommas(str,subs);
            //Look for matching HCW
            string pat=subs[0];
            int index=-1;
            for (int j=0;j<pdat.size();j++) {
                if (pdat[j].code==pat) {
                    index=j;
                    //cout << pat << "\n";
                    break;
                }
            }
            //Assign the HCW as being on the specified ward every day of the study
            loc l;
            l.ward=subs[1];
            l.prob=p.hcw_location_default;
            if (p.diagnostic==0) {
                cout << "Here " << pdat[index].code << " ward " << subs[1] << "    prob " << l.prob << "\n";
            }
            for (int t=min;t<=max;t++) {
                l.date=t;
                pdat[index].locat.push_back(l);
            }
        }
    }
    csv_file.close();
}

void ReadWardMovFromCSV(run_params p, vector<pat>& pdat) {
	ifstream csv_file;
	csv_file.open(p.ward_file.c_str());
	string str;
	vector<int> dates;
	for (int i=0;i<100000;i++) {
		if (!(csv_file >> str)) break;
		if (i>0) {
			//Edit string to remove "
			RemovePunc(str);
			//Split at commas
			vector<string> subs;
			SplitCommas(str,subs);
			//Look for matching patient
			string pat=subs[0];
			int index=-1;
			for (int j=0;j<pdat.size();j++) {
				if (pdat[j].code==pat) {
					index=j;
					//cout << pat << "\n";
					break;
				}
			}
			if (index!=-1) {
				for (int j=4;j<subs.size();j=j+3) {
					loc l;
					l.ward=subs[j];
					if (l.ward.compare("WARD_48")==0) {
						l.ward="WARD_47";
					}
					if (l.ward.compare("WARD_49")==0) {
						l.ward="WARD_47";
					}
					if (l.ward.compare("WARD_50")==0) {
						l.ward="WARD_47";
					}
					if (l.ward.size()>0) {
						vector<int> dmy;
						MakeDMY(j+1,subs,p.pat_delim,dmy);
						int day1=DatetoDay(dmy);
						dmy.clear();
						MakeDMY(j+2,subs,p.pat_delim,dmy);
						int day2=DatetoDay(dmy);
						//cout << j << " " << l.ward << " " << day1 << " " << day2 << "\n";
						for (int k=day1;k<=day2;k++) {
							l.date=k;
							l.prob=1;
							pdat[index].locat.push_back(l);
						}
					}
				}
			}
		}
	}
}

void ReadWardMovBaysFromCSV(run_params& p, vector<pat>& pdata) {
    ifstream csv_file;
    csv_file.open(p.ward_bay_file.c_str());
    string str;
    vector<int> dates;
    int i=-1;
    int spwarn=0;
    while (getline(csv_file,str)) {
        i++;
        if (i>0&&p.error==0) {
            //cout << str << "\n";
            /*if (p.diagnostic==1) {
                cout << "Ward_file string " << str << "\n";
            }*/
            //Edit string to remove "
            RemovePunc(str);
            //Split at commas
            vector<string> subs;
            SplitCommas(str,subs);
            RemoveSpaces(p.ward_file.c_str(),i,spwarn,subs);

            //Look for matching patient
            string pat=subs[0];
            int index=-1;
            for (int j=0;j<pdata.size();j++) {
                if (pdata[j].code==pat) {
                    index=j;
                    break;
                }
            }
            //Found matching patient
            if (index!=-1) {
                int j=1;
                loc l;
                l.ward=subs[j];
                //Check the following readout...
                if (subs[j+1].compare("")!=0) {
                    l.bay=stoi(subs[j+1]);
                    l.bay_size=stoi(subs[j+2]);
                } else {
                    l.bay=-1;
                    l.bay_size=-1;
                }
                if (l.ward.size()>0) {
                    vector<int> dmy;
                    MakeDMY(p,p.ward_file.c_str(),j+3,subs,p.pat_delim,dmy);
                    int day1=DatetoDay(dmy);
                    dmy.clear();
                    if (p.error==0) {
                        MakeDMY(p,p.ward_file.c_str(),j+7,subs,p.pat_delim,dmy);
                        int day2=DatetoDay(dmy);
                        for (int k=day1;k<=day2;k++) {
                            l.date=k;
                            l.prob=1;
                            pdata[index].locat.push_back(l);
                        }
                    }
                }
            }
        } else if (i>0) {
            pat pt;
            //Edit string to remove "
            RemovePunc(str);
            //Split by commas
            vector<string> subs;
            SplitCommas(str,subs);
            RemoveSpaces(p.pat_file.c_str(),i,spwarn,subs);
            cout << "Exclude individual " << subs[0] << " due to error in data\n";
        }
    }
}

void ReadWardMovFromCSVNew(run_params& p, vector<pat>& pdata) {
    ifstream csv_file;
    csv_file.open(p.ward_file.c_str());
    string str;
    vector<int> dates;
    int i=-1;
    int spwarn=0;
    while (getline(csv_file,str)) {
        i++;
        if (i>0&&p.error==0) {
            //cout << str << "\n";
            /*if (p.diagnostic==1) {
                cout << "Ward_file string " << str << "\n";
            }*/
            //Edit string to remove "
            RemovePunc(str);
            //Split at commas
            vector<string> subs;
            SplitCommas(str,subs);
            RemoveSpaces(p.ward_file.c_str(),i,spwarn,subs);

            //Look for matching patient
            string pat=subs[0];
            int index=-1;
            for (int j=0;j<pdata.size();j++) {
                if (pdata[j].code==pat) {
                    index=j;
                    //cout << pat << "\n";
                    break;
                }
            }
            /*if (p.diagnostic==1) {
                cout << subs.size() << " " << index << "\n";
            }*/

            if (index!=-1) {
                //cout << "Patient " << pat << "\n";
                for (int j=1;j<subs.size();j=j+4) {
                    loc l;
                    l.ward=subs[j];
                    /*if (l.ward.compare("WARD_48")==0) {
                        l.ward="WARD_47";
                    }
                    if (l.ward.compare("WARD_49")==0) {
                        l.ward="WARD_47";
                    }
                    if (l.ward.compare("WARD_50")==0) {
                        l.ward="WARD_47";
                    }*/
                    //cout << "Ward " << subs[j] << " date " << subs[j+1] << "\n";
                    if (l.ward.size()>0) {
                        vector<int> dmy;
                        MakeDMY(p,p.ward_file.c_str(),j+1,subs,p.pat_delim,dmy);
                        /*cout << subs[j+1] << "\n";
                        cout << dmy[0] << " " << dmy[1] << " " << dmy[2] << "\n";*/
                        int day1=DatetoDay(dmy);
                    //    cout << day1 << "\n";
                    //    cout << "Error " << p.error << "\n";
                        dmy.clear();
                        if (p.error==0) {
                            MakeDMY(p,p.ward_file.c_str(),j+3,subs,p.pat_delim,dmy);
                            int day2=DatetoDay(dmy);
                        //    cout << subs[j+3] << "\n";
                        //    cout << dmy[0] << " " << dmy[1] << " " << dmy[2] << "\n";
                        //    cout << day2 << "\n";
                            /*if (p.diagnostic==1) {
                                cout << day1 << " " << day2 << "\n";
                            }*/
                            for (int k=day1;k<=day2;k++) {
                                l.date=k;
                                l.prob=1;
                                pdata[index].locat.push_back(l);
                            }
                        }
                    }
                }
            }
        } else if (i>0) {
            pat pt;
            //Edit string to remove "
            RemovePunc(str);
            //Split by commas
            vector<string> subs;
            SplitCommas(str,subs);
            RemoveSpaces(p.pat_file.c_str(),i,spwarn,subs);
            cout << "Exclude individual " << subs[0] << " due to error in data\n";
        }
    }
}



void EditHCWMovData (run_params p, vector<pat>& pdat) {
	//Adds a 12-hour window of uncertainty in the location of HCWs (not patients)
	//Represents e.g. transmission via touching equipment, night shift timings, etc.
	//Assign a probability of 0.5 to the days before and after a known presence
	for (int i=0;i<pdat.size();i++) {
        if (pdat[i].hcw==1) {
            vector<loc> new_loc;
            for (int j=0;j<pdat[i].locat.size();j++) {
                loc l=pdat[i].locat[j];
                int minus=1;
                int plus=1;
                if (j>0&&pdat[i].locat[j-1].date==l.date-1&&pdat[i].locat[j-1].ward==l.ward&&pdat[i].locat[j-1].prob==1) {
                    minus=0;
                }
                if (j<pdat[i].locat.size()-1&&pdat[i].locat[j+1].date==l.date+1&&pdat[i].locat[j+1].ward==l.ward&&pdat[i].locat[j+1].prob==1) {
                    plus=0;
                }
                if (minus==1) {
                    loc lm=l;
                    lm.date--;
                    lm.prob=p.hcw_gap;
                    new_loc.push_back(lm);
                }
                if (plus==1) {
                    loc lp=l;
                    lp.date++;
                    lp.prob=p.hcw_gap;
                    new_loc.push_back(lp);
                }
            }
            for (int j=0;j<new_loc.size();j++) {
                pdat[i].locat.push_back(new_loc[j]);
            }
        }
	}
}

void RemoveIndividualsNoLocation (const run_params p, vector<pat>& pdat) {
    if (p.space==1) { //A consequence of there being no location data is that transmission is considered impossible
        vector<int> noseq;
        for (int i=0;i<pdat.size();i++) {
            if (pdat[i].locat.size()==0) {
                cout << "No location data for individual " << pdat[i].code << ": Excluding from input\n";
                noseq.push_back(i);
            }
        }
        sort(noseq.begin(),noseq.end());
        reverse(noseq.begin(),noseq.end());
        for (int i=0;i<noseq.size();i++) {
            pdat.erase(pdat.begin()+noseq[i]);
        }
    }
}

void ReadSubsets (run_params p, vector< vector<int> >& subsets) {
    ifstream str_file;
    cout << p.sub_file << "\n";
    str_file.open(p.sub_file.c_str());
    string str;
    for (int i=0;i<100000;i++) {
        if (!(str_file >> str)) break;
        //cout << str << "\n";
        vector<string> subs;
        SplitCommas(str,subs);
        vector<int> vals;
        for (int j=0;j<subs.size();j++) {
            vals.push_back(atoi(subs[j].c_str()));
            //cout << atoi(subs[j].c_str()) << "\n";
        }
        subsets.push_back(vals);
    }
}













void ReadSymptomData (vector<pat>& pdat) {
	ifstream sym_file;
	sym_file.open("../Model_data/HCW_symptoms.dat");
	for (int i=0;i<1000000;i++) {
		pat p;
		string code;
		int x;
		if (!(sym_file >> code)) break;
		p.code=code;
		if (!(sym_file >> x)) break;
		p.time_s=x;
		if (!(sym_file >> x)) break;
		p.time_s_cert=x;
		p.hcw=1;
		pdat.push_back(p);
	}
	sym_file.close();
	sym_file.open("../Model_data/Patient_symptoms.dat");
	for (int i=0;i<1000000;i++) {
		pat p;
		string code;
		int x;
		if (!(sym_file >> code)) break;
		p.code=code;
		if (!(sym_file >> x)) break;
		p.time_s=x;
		if (!(sym_file >> x)) break;
		p.time_s_cert=x;
		p.hcw=0;
		pdat.push_back(p);
	}
	sym_file.close();
}

void ReadSeqTimes(vector<pat>& pdat) {
	ifstream st_file;
	st_file.open("../Model_data/Sequence_dates.dat");
	for (int i=0;i<1000000;i++) {
		string code;
		int x;
		if (!(st_file >> code)) break;
		if (!(st_file >> x)) break;
		for (int j=0;j<pdat.size();j++) {
			if (pdat[j].code_match==code) {
				pdat[j].time_seq=x;
			}
		}
	}

}
 
void ReadLocationData (vector<pat>& pdat) {
	ifstream loc_file;
	loc_file.open("../Model_data/HCW_locations.dat");
	for (int i=0;i<1000000;i++) {
		string code;
		int x;
		if (!(loc_file >> code)) break;
		for (int j=0;j<pdat.size();j++) {
			if (pdat[j].code==code) {
				if (!(loc_file >> x)) break;
				pdat[j].location.push_back(x);
				break;
			}
		}
	}
	loc_file.close();
	loc_file.open("../Model_data/Patient_locations.dat");
	for (int i=0;i<1000000;i++) {
		string code;
		int x;
		if (!(loc_file >> code)) break;
		for (int j=0;j<pdat.size();j++) {
			if (pdat[j].code==code) {
				if (!(loc_file >> x)) break;
				pdat[j].location.push_back(x);
				break;
			}
		}
	}
	loc_file.close();
}

void ReadCommunityDistances (vector<int>& comm_dist) {
	ifstream com_file;
	com_file.open("../Model_data/CUH_sequences_uniq_dist.out");
	for (int i=0;i<1000000;i++) {
		int x;
		if (!(com_file >> x)) break;
		if (!(com_file >> x)) break;
		if (!(com_file >> x)) break;
		comm_dist.push_back(x);
	}
}

void WriteBestRoot (int index, const vector<treestore>& ts_roots) {
	cout << "Maximum output from roots:\n";
	for (int i=0;i<ts_roots[index].origins.size();i++) {
		cout << ts_roots[index].origins[i] << " " << ts_roots[index].dests[i] << " " << ts_roots[index].weights[i] << "\n";
	}
	cout << "Likelihood " << ts_roots[index].logL << "\n";
}

/*
void WriteMLDetails (run_params p, const int index, const vector< vector<int> >& bin, const vector<treestore_plus>& ts_bin, vector<pat>& pdat) {
	cout << "Maximum likelihood reconstruction\n";
	for (int i=0;i<bin[index].size();i++) {
		cout << bin[index][i];
	}
	cout << "\n";
	cout << "Variants\n";
	for (int i=0;i<ts_bin[index].ts[0].variants.size();i++) {
		cout << i << " " << pdat[i].code << " " << pdat[i].time_s << " " << pdat[i].hcw << " ";
		for (int j=0;j<ts_bin[index].ts[0].variants[i].locus.size();j++) {
			cout << ts_bin[index].ts[0].variants[i].locus[j] << ts_bin[index].ts[0].variants[i].allele[j] << " ";
		}
		cout << "\n";
	}
	for (int i=0;i<ts_bin[index].ts.size();i++) {
		cout << "Cluster " << i << "\n";
		for (int j=0;j<ts_bin[index].ts[i].individuals.size();j++) {
			cout << j << " " << ts_bin[index].ts[i].individuals[j] << "\n";
		}
		for (int j=0;j<ts_bin[index].ts[i].origins.size();j++) {
			if (j<10) {
				cout << " ";
			}
			cout << j << " ";
			if (ts_bin[index].ts[i].origins[j]<10) {
				cout << " ";
			}
			cout << ts_bin[index].ts[i].origins[j] << " ";
			if (ts_bin[index].ts[i].dests[j]<10) {
				cout << " ";
			}
			cout << ts_bin[index].ts[i].dests[j] << "  ";
			cout << ts_bin[index].ts[i].individuals[ts_bin[index].ts[i].origins[j]] << " " << ts_bin[index].ts[i].individuals[ts_bin[index].ts[i].dests[j]] << " " << ts_bin[index].ts[i].weights[j] << " ";
			Thresholds(p,ts_bin[index].ts[i].weights[j]);
			cout << "\n";
		}
		cout << "Log L " << ts_bin[index].ts[i].logL << "\n";
	}
}*/
